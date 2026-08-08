"""End-to-end per-sample pipeline: sketch -> boundary plan -> apply ->
identity screen -> evidence/release class.

Single source of truth shared by the `extract` CLI command and the M6
samplesheet driver, so audit and production modes agree by construction
(they call the same function with `apply_trims` toggled). See
docs/release_qc_and_terminal_trimming_plan.md sections 4, 10, and M6.
"""

import importlib.metadata
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .apply import override_rule_record, plan_from_file
from .processors.collapsed import parse_collapsed_fasta
from .processors.identity_screen import identity_screen
from .processors.rpf_extractor import RPFExtractor
from .processors.sketch import Sketch, classify_terminal_signal
from .processors.types import RPFExtractionResult
from .release import build_evidence, check_self_consistency, write_evidence

logger = logging.getLogger(__name__)

MAX_IDENTITY_SCREEN_READS = 50_000


def _getrpf_version() -> str:
    try:
        return importlib.metadata.version("getRPF")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _known_adapters(extractor: RPFExtractor) -> List[str]:
    """Flatten adapter sequences across every architecture in the
    extractor's loaded database into a candidate list for the b_adapter
    boundary estimator. This is deliberately the *whole* database, not
    just the architecture RPFExtractor's own detection settles on --
    b_adapter is "does a known adapter match", independent of which
    architecture, if any, wins the separate detection process."""
    adapters: List[str] = []
    seen = set()
    for arch in extractor.architecture_db.architectures:
        for adapter in getattr(arch, "adapter_sequences", None) or []:
            if adapter not in seen:
                seen.add(adapter)
                adapters.append(adapter)
    return adapters


def _weighted_extracted_population(output_file: Path) -> Tuple[List[str], Dict[int, int]]:
    """Read the collapsed FASTA written by extract_rpfs and rebuild the
    count-weighted read population: a length_distribution weighted by true
    read counts (not unique-sequence counts), and a bounded weighted read
    sample for the identity screen's sequence-level checks (duplication,
    contamination). Reading the collapsed FASTA's sequence lines directly
    would silently treat every unique sequence as exactly one read, which
    both zeroes out the duplication signal and skews the length shape.
    """
    collapsed_path = output_file.with_suffix(".collapsed.fa")
    if not collapsed_path.exists():
        return [], {}

    sequences, counts = parse_collapsed_fasta(collapsed_path, "seq{id}_x{count}")

    length_distribution: Dict[int, int] = {}
    reads: List[str] = []
    for header, seq in sequences.items():
        count = counts.get(header, 1)
        length_distribution[len(seq)] = length_distribution.get(len(seq), 0) + count
        if len(reads) < MAX_IDENTITY_SCREEN_READS:
            reads.extend([seq] * min(count, MAX_IDENTITY_SCREEN_READS - len(reads)))

    return reads, length_distribution


def _quality_evidence(sketch: Sketch) -> Dict:
    """Run the M1 poly-G/adapter discriminator against the sketch's pooled
    3' signal, so basecaller_artifact_hold can actually fire (section 7,
    section 8)."""
    entropy_3p = sketch.pooled.entropy_3p
    composition_3p = sketch.pooled.composition_3p
    if not entropy_3p or not composition_3p or not composition_3p[0]:
        return {"three_prime_q_collapse": False}

    dominant_base, _ = max(composition_3p[0].items(), key=lambda kv: kv[1])
    quality_mean = sketch.quality_3p[0] if sketch.quality_3p else None
    classification = classify_terminal_signal(entropy_3p[0], quality_mean, dominant_base)
    return {"three_prime_q_collapse": classification == "basecaller_artifact"}


def run_sample(
    input_file: Path,
    output_file: Path,
    format: str = "fastq",
    sample_id: Optional[str] = None,
    infer_reads: int = 500_000,
    max_reads: Optional[int] = None,
    apply_trims: bool = True,
    override_trims: Optional[Dict[int, Dict[str, int]]] = None,
    architecture_db: Optional[Path] = None,
    seqspec_dir: Optional[Path] = None,
    collapse_output: bool = True,
    collapsed_only: bool = False,
    generate_seqspec: bool = False,
    protocol: Optional[Dict] = None,
) -> Tuple[Dict, RPFExtractionResult]:
    """Run the full reference-free pipeline for one sample and return its
    evidence object.

    apply_trims=True (production): the boundary plan's safe_to_trim rules
        are applied via RPFExtractor.extract_rpfs.
    apply_trims=False (audit): the same plan is computed, but override_trims
        is not passed to extract_rpfs, so no reads are modified. The
        `trim_rules_proposed_not_applied`/`trim_rules_applied` distinction
        in the returned evidence reflects what *would* have been applied.

    override_trims, if given, bypasses the *inferred* plan (explicit frozen
    rules replace it), though boundary estimation still runs so the sketch
    and boundary_estimates_3p in the evidence object stay informative.
    """
    sample_id = sample_id or input_file.stem

    extractor = RPFExtractor(architecture_db_path=architecture_db, seqspec_dir=seqspec_dir)

    (
        inferred_trims, inferred_applied, inferred_proposed, decisions, sketch,
    ) = plan_from_file(
        input_file, format=format, infer_reads=infer_reads,
        known_adapters=_known_adapters(extractor),
    )

    analysis_reads = [
        read
        for reads_at_length in sketch.per_length_reads.values()
        for read in reads_at_length
    ][:50_000]

    if override_trims is not None:
        applied_rules = [
            override_rule_record(length, rule.get("trim_3p", 0))
            for length, rule in override_trims.items()
        ]
        proposed_rules: List[Dict] = []
        plan_trims = override_trims
    else:
        applied_rules = inferred_applied
        proposed_rules = inferred_proposed
        plan_trims = inferred_trims

    trims_to_apply = plan_trims if apply_trims else None

    result = extractor.extract_rpfs(
        input_file, output_file, format=format, max_reads=max_reads,
        generate_seqspec=generate_seqspec,
        collapse_output=collapse_output, collapsed_only=collapsed_only,
        override_trims=trims_to_apply,
        sample_reads=analysis_reads,
    )

    self_consistency = check_self_consistency(sketch, trims_to_apply or {})
    quality_evidence = _quality_evidence(sketch)

    extracted_reads, extracted_length_distribution = _weighted_extracted_population(output_file)
    biological_screen = identity_screen(extracted_reads, extracted_length_distribution)

    if not apply_trims:
        # Audit mode: nothing was actually applied; report the plan as
        # entirely proposed, since no reads were modified.
        applied_for_evidence: List[Dict] = []
        proposed_for_evidence = applied_rules + proposed_rules
    else:
        applied_for_evidence = applied_rules
        proposed_for_evidence = proposed_rules

    boundary_estimates_3p = {str(length): d.to_dict() for length, d in decisions.items()}

    evidence = build_evidence(
        sample_id=sample_id,
        getrpf_version=_getrpf_version(),
        read_count_input=result.input_reads,
        read_count_output=result.extracted_rpfs,
        boundary_estimates_3p=boundary_estimates_3p,
        applied_rules=applied_for_evidence,
        proposed_rules=proposed_for_evidence,
        self_consistency=self_consistency,
        biological_screen=biological_screen,
        quality_evidence=quality_evidence,
        protocol=protocol,
    )

    evidence_path = output_file.with_suffix(".evidence.json")
    write_evidence(evidence, evidence_path)

    release_qc_path = output_file.with_suffix(".release_qc.json")
    with open(release_qc_path, "w") as f:
        json.dump(
            {
                "sample_id": evidence["sample_id"],
                "release_class": evidence["release_class"],
                "biological_confirmation": evidence["biological_confirmation"],
                "warnings": evidence["warnings"],
            },
            f, indent=2,
        )

    applied_path = output_file.with_suffix(".trim_rules.applied.json")
    proposed_path = output_file.with_suffix(".trim_rules.proposed.json")
    with open(applied_path, "w") as f:
        json.dump(applied_for_evidence, f, indent=2)
    with open(proposed_path, "w") as f:
        json.dump(proposed_for_evidence, f, indent=2)

    return evidence, result
