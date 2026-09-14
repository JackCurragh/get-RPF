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
from .assemble import StructureInference, infer_structure
from .config import InferenceConfig
from .observe import read_fastq
from .report import _plain, write_report
from .transform import ExtractionSummary, extract_reads, validate_extraction

REPORT_SCHEMA = "getrpf.structure.extraction/1"


@dataclass(frozen=True)
class StructureExtraction:
    result: StructureInference
    summary: Optional[ExtractionSummary]
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
) -> StructureExtraction:
    headers, sequences, _ = read_fastq(input_file, config.sample_reads)
    if not sequences:
        raise ValueError(f"no reads in {input_file}")
    result = infer_structure(sequences, headers, config)
    if star_index is not None and result.architecture is not None:
        result = _alignment_checked(result, sequences, star_index, config)

    prefix = output_prefix(output_file)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_report(result, prefix.parent, prefix.name)

    summary: Optional[ExtractionSummary] = None
    reads_path: Optional[Path] = None
    architecture = result.architecture
    if result.transform.emit and architecture is not None:
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
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return StructureExtraction(result, summary, report_path, reads_path)


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
    summary: Optional[ExtractionSummary],
    reads_path: Optional[Path],
    input_file: Path,
    inference_reads: int,
    config: InferenceConfig,
) -> Dict[str, Any]:
    return {
        "schema": REPORT_SCHEMA,
        "input": input_file.name,
        "status": "emitted" if summary is not None else "withheld",
        "inference_reads": inference_reads,
        "architecture": (
            result.architecture.describe() if result.architecture else None
        ),
        "transform": _plain(result.transform),
        "umi": result.q4.explanation,
        "extraction": summary.to_dict() if summary is not None else None,
        "validation": (
            _plain(validate_extraction(summary, config))
            if summary is not None
            else None
        ),
        "reads": reads_path.name if reads_path is not None else None,
    }
