"""Raw FASTQ to RPFs in one call (spec §7): infer the structure on a bounded
sample, check it against a reference when one is given, then apply it to
every read.

This is what ``getRPF extract --infer-structure`` runs, once per sample in a
pipeline. Every output shares the prefix of the requested output file, as the
legacy ``extract`` names them, so pipeline globs keep working:
``<prefix>.structure.{json,txt}``, ``<prefix>.seqspec.yaml``,
``<prefix>.extraction_report.json`` and the reads (the requested FASTQ, or
``<prefix>.collapsed.fa``). A withheld transform writes the reports and no
reads.
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from ..processors.collapsed import TwoStageCollapser
from .assemble import StructureInference, build_preparation_architecture, infer_structure
from .anchors import locate
from .config import InferenceConfig
from .observe import read_fastq, reservoir_fastq_sample
from .report import _plain, write_report
from .transform import (
    ExtractionSummary,
    PreparationSummary,
    extract_reads,
    prepare_reads,
    validate_extraction,
)
from .model import Alternative, Answer, Architecture, Block, Evidence, JunctionCall, Status

REPORT_SCHEMA = "getrpf.structure.extraction/1"


@dataclass(frozen=True)
class StructureExtraction:
    result: StructureInference
    summary: Optional[Any]
    """None when the transform was withheld."""
    report_path: Path
    reads_path: Optional[Path]
    """None when withheld, audited, or no insert was accepted."""

    @property
    def emitted(self) -> bool:
        return self.summary is not None

    def describe(self) -> str:
        architecture = self.result.architecture
        shape = architecture.describe() if architecture else "no architecture"
        if self.summary is None:
            reasons = "; ".join(self.result.transform.reasons)
            return f"Transform withheld ({shape}): {reasons}. Report {self.report_path}"
        if isinstance(self.summary, PreparationSummary):
            return (
                f"Alignment preparation: {self.summary.emitted_reads} of "
                f"{self.summary.input_reads} reads emitted -> "
                f"{self.reads_path or 'no reads written'}; report {self.report_path}"
            )
        return (
            f"{shape}: {self.summary.accepted} of {self.summary.input_reads} reads "
            f"accepted -> {self.reads_path or 'no reads written'}; "
            f"{len(self.result.transform.flags)} flag(s); report {self.report_path}"
        )


def output_prefix(output_file: Path) -> Path:
    """``x_trimmed.fastq.gz`` -> ``x_trimmed``: the prefix every output shares."""
    name = output_file.name
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for suffix in (".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return output_file.with_name(name)


def infer_and_extract(
    input_file: Path,
    output_file: Path,
    config: InferenceConfig,
    star_index: Optional[Path] = None,
    collapsed: bool = False,
    audit: bool = False,
    max_reads: Optional[int] = None,
    prepare_for_alignment: bool = False,
    trim_recipe: Optional[Dict[str, Any]] = None,
) -> StructureExtraction:
    headers, sequences, _ = read_fastq(input_file, config.sample_reads)
    if not sequences:
        raise ValueError(f"no reads in {input_file}")
    result = infer_structure(sequences, headers, config)
    sampling = {
        "initial": {"size": len(sequences), "method": "prefix", "seed": None},
        "retry": None,
    }
    if _pending_decision(result):
        retry_seed = 1729
        retry_headers, retry_reads = reservoir_fastq_sample(
            input_file, config.sample_reads, retry_seed
        )
        sampling["retry"] = {
            "size": len(retry_reads),
            "method": "streaming_reservoir_across_fastq",
            "seed": retry_seed,
        }
        if retry_reads:
            retry = infer_structure(retry_reads, retry_headers, config)
            result = _merge_retry(result, retry, config)
    if not prepare_for_alignment and star_index is not None and result.architecture is not None:
        result = _alignment_checked(result, sequences, star_index, config)

    prefix = output_prefix(output_file)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_report(result, prefix.parent, prefix.name)

    summary: Optional[Any] = None
    reads_path: Optional[Path] = None
    architecture = result.architecture
    preparation: Optional[PreparationSummary] = None
    if prepare_for_alignment:
        if max_reads is not None:
            raise ValueError("alignment preparation requires a complete FASTQ parse")
        recipe_decisions: Dict[str, Any] = {}
        adapter_override = None
        recipe_adapter = (trim_recipe or {}).get("adapter")
        if recipe_adapter and str(recipe_adapter).lower() != "passed":
            recipe_adapter = str(recipe_adapter).upper().replace("U", "T")
            hits = [locate(read, recipe_adapter, config) for read in sequences]
            recipe_support = sum(hit is not None for hit in hits) / len(hits)
            inferred_adapter = getattr(getattr(result.q1.value, "adapter", None), "sequence", None)
            inferred_mode = getattr(result.q1.value, "start_mode", None)
            recipe_starts = [hit.start for hit in hits if hit is not None]
            recipe_mode = None
            if recipe_starts:
                from collections import Counter as _Counter

                recipe_mode = _Counter(recipe_starts).most_common(1)[0][0]
            same_boundary = inferred_mode is None or recipe_mode is None or abs(inferred_mode - recipe_mode) <= 2
            compatible_adapter = (
                inferred_adapter is None
                or inferred_adapter.startswith(recipe_adapter)
                or recipe_adapter.startswith(inferred_adapter)
            )
            longer_inferred = bool(inferred_adapter and inferred_adapter.startswith(recipe_adapter)
                                   and len(inferred_adapter) > len(recipe_adapter))
            if recipe_support >= config.min_anchor_support and same_boundary and compatible_adapter and not longer_inferred:
                adapter_override = (
                    recipe_adapter
                    if inferred_adapter is None or len(recipe_adapter) > len(inferred_adapter)
                    else None
                )
                recipe_decisions["adapter"] = {
                    "status": "verified", "support": round(recipe_support, 4),
                    "action": "supported longer supplied adapter used" if adapter_override and inferred_adapter
                              else ("supplied adapter used" if adapter_override else "inferred adapter retained"),
                }
            else:
                recipe_decisions["adapter"] = {
                    "status": "declined", "support": round(recipe_support, 4),
                    "reason": ("observed layout support below threshold" if recipe_support < config.min_anchor_support
                               else "reported adapter boundary or sequence conflicts with observed layout"),
                }
            if longer_inferred and recipe_support >= config.min_anchor_support:
                recipe_decisions["adapter"] = {
                    "status": "verified", "support": round(recipe_support, 4),
                    "action": "longer inferred adapter cut retained",
                    "reported_adapter_length": len(recipe_adapter),
                    "inferred_adapter_length": len(inferred_adapter),
                }
        prep_architecture = build_preparation_architecture(result, config, adapter_override)
        for field, answer in (("5p_barcode_length", result.q2), ("3p_barcode_length", result.q3)):
            requested = (trim_recipe or {}).get(field)
            supported = getattr(answer.value, "technical_length", None)
            if requested is not None:
                verified = (supported == 0 if requested == 0 else
                            supported is not None and requested <= supported)
                recipe_decisions[field] = {
                    "requested": requested,
                    "supported_length": supported,
                    "status": "verified" if verified else "declined",
                    "reason": None if verified else "the full reported span is not supported by read evidence",
                    "action": "reported span used" if verified and requested else "evidence-supported inferred cuts used",
                }
                if verified and requested:
                    prep_architecture = _limit_recipe_span(prep_architecture, answer.question, requested)
        counts = Counter()
        preparation = prepare_reads(
            input_file,
            prep_architecture,
            config,
            collect=counts,
            pretrimmed=(trim_recipe or {}).get("input_state") == "pretrimmed",
            five_prime_resolved=(
                result.q2.status.value in ("resolved", "interval")
                and result.q2.value is not None
            ),
            three_prime_resolved=(
                result.q3.status.value in ("resolved", "interval")
                and result.q3.value is not None
            ) or (trim_recipe or {}).get("input_state") == "pretrimmed",
        )
        summary = preparation
        target = prefix.with_name(f"{prefix.name}.collapsed.fa")
        if counts and TwoStageCollapser().write_collapsed_fasta(counts, target)["output_path"]:
            reads_path = target
    elif result.transform.emit and architecture is not None:
        counts: Optional[Counter[str]] = Counter() if collapsed and not audit else None
        fastq = None if (collapsed or audit) else output_file
        summary = extract_reads(
            input_file, fastq, architecture, config, limit=max_reads, collect=counts
        )
        if fastq is not None:
            reads_path = fastq
        elif counts is not None:
            target = prefix.with_name(f"{prefix.name}.collapsed.fa")
            if TwoStageCollapser().write_collapsed_fasta(counts, target)["output_path"]:
                reads_path = target

    report_path = prefix.with_name(f"{prefix.name}.extraction_report.json")
    report = _report(result, summary, reads_path, input_file, len(sequences), config)
    if preparation is not None:
        report["preparation"] = {
            **preparation.to_dict(),
            "recipe": trim_recipe,
            "recipe_provenance": (trim_recipe or {}).get("source"),
            "recipe_decisions": recipe_decisions,
        }
    report["sampling"] = sampling
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return StructureExtraction(result, summary, report_path, reads_path)


def _pending_decision(result: StructureInference) -> bool:
    acceptable = {Status.RESOLVED, Status.INTERVAL}
    anchor_pending = result.q1.status not in acceptable and not (
        result.q1.status is Status.NOT_OBSERVABLE and result.observation.input_state == "trimmed"
    )
    return anchor_pending or any(answer.status not in acceptable for answer in (result.q2, result.q3))


def _limit_recipe_span(architecture: Architecture, question: str, length: int) -> Architecture:
    """Honor a reported span only inside the independently supported block run."""
    blocks = list(architecture.blocks)
    insert_at = next(i for i, block in enumerate(blocks) if block.type == "insert")
    if question == "Q2":
        indices = list(range(insert_at))
        reverse = False
    else:
        indices = [i for i in range(insert_at + 1, len(blocks))
                   if blocks[i].type not in ("tail", "adapter")]
        reverse = True
    remaining = length
    chosen = {}
    walk = list(reversed(indices)) if reverse else indices
    for index in walk:
        block = blocks[index]
        take = min(block.length[1], remaining)
        remaining -= take
        if take:
            chosen[index] = replace(block, length=(take, take),
                                    sequence=(block.sequence[:take] if block.sequence else None))
        if remaining <= 0:
            break
    for index in indices:
        blocks[index] = chosen.get(index)
    blocks = [block for block in blocks if block is not None]
    return replace(architecture, blocks=tuple(blocks))


def _merge_retry(
    first: StructureInference, retry: StructureInference, config: InferenceConfig
) -> StructureInference:
    """Use retry evidence for pending calls and retain any contested trim span."""
    q1 = _merge_answer(first.q1, retry.q1)
    q2 = _merge_answer(first.q2, retry.q2)
    q3 = _merge_answer(first.q3, retry.q3)
    junctions = replace(retry.junctions, q2=q2, q3=q3)
    from .assemble import answer_umi, build_architecture, decide_transform

    q4 = answer_umi(retry.observation, q2, q3)
    transform = decide_transform(retry.observation, q1, q2, q3, config)
    architecture = build_architecture(q1, junctions, transform, config)
    return replace(retry, q1=q1, junctions=junctions, q4=q4, transform=transform, architecture=architecture)


def _merge_answer(first: Answer, retry: Answer) -> Answer:
    acceptable = {Status.RESOLVED, Status.INTERVAL}
    first_ok, retry_ok = first.status in acceptable, retry.status in acceptable
    if not first_ok:
        return retry if retry_ok else first
    if not retry_ok:
        return first
    first_length = getattr(first.value, "technical_length", None)
    retry_length = getattr(retry.value, "technical_length", None)
    if first.question == "Q1":
        first_adapter = getattr(getattr(first.value, "adapter", None), "sequence", None)
        retry_adapter = getattr(getattr(retry.value, "adapter", None), "sequence", None)
        if first_adapter and retry_adapter and first_adapter != retry_adapter:
            return Answer(
                "Q1", None, Status.AMBIGUOUS,
                first.evidence + retry.evidence + (Evidence("retry", "contested_adapters", (first_adapter, retry_adapter), 2),),
                (Alternative(first_adapter, "untested", "supported only in the prefix sample"),
                 Alternative(retry_adapter, "untested", "supported only in the reservoir sample")),
                "The two supported samples call different adapter boundaries; the adapter cut is withheld.",
            )
        return first
    if first_length is None or retry_length is None or first_length == retry_length:
        return first
    retained = min(first_length, retry_length)
    call = first.value
    blocks = []
    remaining = retained
    for kind, length in call.blocks:
        if remaining <= 0:
            break
        take = min(length, remaining)
        blocks.append((kind, take))
        remaining -= take
    clipped = replace(call, technical_length=retained, blocks=tuple(blocks))
    return Answer(
        first.question, clipped, Status.AMBIGUOUS,
        first.evidence + retry.evidence + (Evidence(
            "retry", "contested_trim_length", {"prefix": first_length, "reservoir": retry_length, "retained": retained}, 2,
            note="only the common supported span is removed",
        ),),
        first.alternatives + retry.alternatives,
        f"Supported samples disagree on the technical span ({first_length} vs {retry_length} nt); retaining the contested span after {retained} nt.",
    )


def _alignment_checked(
    result: StructureInference,
    sequences: Sequence[str],
    star_index: Path,
    config: InferenceConfig,
) -> StructureInference:
    """Fold in the step 4 check. A check that cannot run is a flag, not a
    failure, so one broken aligner or reference does not stop a cohort."""
    from .align import apply_alignment_check

    try:
        checked, _ = apply_alignment_check(result, sequences, star_index, config)
        return checked
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        detail = str(error)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            detail = str(error.stderr).strip()[-300:]
        flag = f"alignment check not run: {detail}"
        transform = replace(result.transform, flags=result.transform.flags + (flag,))
        return replace(result, transform=transform)


def _report(
    result: StructureInference,
    summary: Optional[Any],
    reads_path: Optional[Path],
    input_file: Path,
    inference_reads: int,
    config: InferenceConfig,
) -> Dict[str, Any]:
    return {
        "schema": REPORT_SCHEMA,
        "input": input_file.name,
        "status": (
            "prepared" if isinstance(summary, PreparationSummary)
            else ("emitted" if summary is not None else "withheld")
        ),
        "inference_reads": inference_reads,
        "architecture": (
            result.architecture.describe() if result.architecture else None
        ),
        "transform": _plain(result.transform),
        "umi": result.q4.explanation,
        "extraction": summary.to_dict() if summary is not None else None,
        "validation": (
            _plain(validate_extraction(summary, config))
            if isinstance(summary, ExtractionSummary)
            else None
        ),
        "reads": reads_path.name if reads_path is not None else None,
    }
