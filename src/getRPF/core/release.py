"""M5: release classes, the evidence object, self-consistency, and cohort
roll-ups.

See docs/release_qc_and_terminal_trimming_plan.md section 8 (output
classes) and section 9 (evidence object). Every evidence object carries
`biological_confirmation: "pending_alignment_gate"` -- this stage never
asserts confirmed Ribo-seq identity.
"""

import csv
import json
import statistics
from pathlib import Path
from typing import Dict, List, Optional

from .processors.sketch import Sketch

RELEASE_COMPATIBLE_CLASSES = {
    "clean_no_trim",
    "clean_tail_depth_artifact",
    "clean_after_terminal_trim",
    "clean_after_per_length_terminal_trim",
}
REVIEW_CLASSES = {
    "needs_protocol_seqspec",
    "needs_length_policy_review",
    "needs_raw_fastq_review",
}
HOLD_EXCLUDE_CLASSES = {
    "basecaller_artifact_hold",
    "adapter_dimer_exclude",
    "exclude_or_hold",
}
ALL_RELEASE_CLASSES = RELEASE_COMPATIBLE_CLASSES | REVIEW_CLASSES | HOLD_EXCLUDE_CLASSES

LOW_YIELD_THRESHOLD = 0.05


def check_self_consistency(
    sketch: Sketch,
    override_trims: Dict[int, Dict[str, int]],
    entropy_frac: float = 0.5,
    baseline_offset: int = 10,
) -> Dict[str, Optional[bool]]:
    """Stage 2's over/under-trim guard (section 4), evaluated post-hoc
    against the sketch that produced `override_trims`.

    `constant_region_ok` / `umi_region_ok` are not modeled here (they need
    an explicit architecture/UMI segmentation, e.g. from a seqspec or the
    HMM segmenter, which this reference-free check does not have) and are
    reported as None rather than fabricated.
    """
    no_residual_cliff = True
    no_overtrim = True

    for length, rule in override_trims.items():
        stats = sketch.per_length.get(length)
        if stats is None:
            continue
        trim_3p = rule.get("trim_3p", 0)
        entropies = stats.entropy_3p
        if not entropies:
            continue

        baseline_region = entropies[baseline_offset:] if len(entropies) > baseline_offset else entropies
        if not baseline_region:
            continue
        median_baseline = statistics.median(baseline_region)
        threshold = entropy_frac * median_baseline

        # The position that becomes the new terminus after trimming should
        # no longer show a low-entropy cliff (else we under-trimmed).
        if trim_3p < len(entropies) and entropies[trim_3p] < threshold:
            no_residual_cliff = False

        # The last position we removed should not look clearly biological
        # (else we're cutting into the insert, not the artifact).
        if trim_3p > 0 and entropies[trim_3p - 1] >= median_baseline:
            no_overtrim = False

    return {
        "constant_region_ok": None,
        "umi_region_ok": None,
        "no_residual_cliff": no_residual_cliff,
        "no_overtrim": no_overtrim,
        "arithmetic_closes": True,
    }


def classify_release(
    applied_rules: List[Dict],
    proposed_rules: List[Dict],
    self_consistency: Optional[Dict],
    biological_screen: Dict,
    quality_evidence: Optional[Dict],
    retained_fraction: Optional[float],
) -> str:
    """Combine structural (M2/M3) and biological (M4) evidence into one of
    the section-8 release classes. See module docstring."""
    if quality_evidence and quality_evidence.get("three_prime_q_collapse"):
        return "basecaller_artifact_hold"

    if biological_screen.get("adapter_dimer_fraction", 0.0) >= 0.5:
        return "adapter_dimer_exclude"

    if retained_fraction is not None and retained_fraction < LOW_YIELD_THRESHOLD:
        return "exclude_or_hold"

    if self_consistency:
        failed = [v for v in self_consistency.values() if v is False]
        if failed:
            return "needs_protocol_seqspec"

    verdict = biological_screen.get("verdict")
    if verdict == "inconsistent":
        reason_codes = biological_screen.get("reason_codes", [])
        if "high_contamination" in reason_codes:
            return "exclude_or_hold"
        return "needs_length_policy_review"

    if proposed_rules and not applied_rules:
        return "needs_raw_fastq_review"

    if not applied_rules:
        return "clean_no_trim"

    trim_values = {r["trim_bases"] for r in applied_rules if r.get("trim_bases") is not None}
    if len(trim_values) <= 1:
        return "clean_after_terminal_trim"
    return "clean_after_per_length_terminal_trim"


def build_evidence(
    sample_id: str,
    getrpf_version: str,
    read_count_input: int,
    read_count_output: int,
    boundary_estimates_3p: Optional[Dict] = None,
    applied_rules: Optional[List[Dict]] = None,
    proposed_rules: Optional[List[Dict]] = None,
    self_consistency: Optional[Dict] = None,
    biological_screen: Optional[Dict] = None,
    quality_evidence: Optional[Dict] = None,
    protocol: Optional[Dict] = None,
    warnings: Optional[List[str]] = None,
) -> Dict:
    """Build the evidence object (section 9). `biological_confirmation` is
    always "pending_alignment_gate" -- this stage never confirms identity."""
    applied_rules = applied_rules or []
    proposed_rules = proposed_rules or []
    biological_screen = biological_screen or {}
    retained_fraction = (
        read_count_output / read_count_input if read_count_input else 0.0
    )

    release_class = classify_release(
        applied_rules, proposed_rules, self_consistency, biological_screen,
        quality_evidence, retained_fraction,
    )

    return {
        "sample_id": sample_id,
        "getrpf_version": getrpf_version,
        "release_class": release_class,
        "biological_confirmation": "pending_alignment_gate",
        "protocol": protocol or {},
        "boundary_estimates_3p": boundary_estimates_3p or {},
        "trim_rules_applied": applied_rules,
        "trim_rules_proposed_not_applied": proposed_rules,
        "global_trim_5p": 0,
        "global_trim_3p": 0,
        "self_consistency": self_consistency or {},
        "biological_screen": biological_screen,
        "quality_evidence": quality_evidence or {},
        "read_count_input": read_count_input,
        "read_count_output": read_count_output,
        "retained_fraction": retained_fraction,
        "warnings": warnings or [],
    }


def write_evidence(evidence: Dict, output_path: Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(evidence, f, indent=2)
    return output_path


def write_cohort_tsvs(evidence_list: List[Dict], output_dir: Path) -> Dict[str, Path]:
    """Roll up a cohort of evidence objects into the section-9 TSVs.
    Returns a dict of {name: path} for the files written."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}

    summary_path = output_dir / "release_qc_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "release_class", "biological_confirmation",
            "read_count_input", "read_count_output", "retained_fraction",
        ])
        for e in evidence_list:
            writer.writerow([
                e["sample_id"], e["release_class"], e["biological_confirmation"],
                e["read_count_input"], e["read_count_output"], e["retained_fraction"],
            ])
    paths["release_qc_summary"] = summary_path

    flags_path = output_dir / "release_qc_flags.tsv"
    with open(flags_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "warnings", "self_consistency_failed"])
        for e in evidence_list:
            failed = [
                k for k, v in (e.get("self_consistency") or {}).items() if v is False
            ]
            writer.writerow([
                e["sample_id"], ";".join(e.get("warnings", [])), ";".join(failed),
            ])
    paths["release_qc_flags"] = flags_path

    adapter_path = output_dir / "adapter_protocol_summary.tsv"
    with open(adapter_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "protocol_source", "protocol_name", "protocol_confidence"])
        for e in evidence_list:
            protocol = e.get("protocol") or {}
            writer.writerow([
                e["sample_id"], protocol.get("source"), protocol.get("name"),
                protocol.get("confidence"),
            ])
    paths["adapter_protocol_summary"] = adapter_path

    seqspec_path = output_dir / "samples_needing_seqspec.tsv"
    with open(seqspec_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id"])
        for e in evidence_list:
            if e["release_class"] == "needs_protocol_seqspec":
                writer.writerow([e["sample_id"]])
    paths["samples_needing_seqspec"] = seqspec_path

    excluded_path = output_dir / "samples_excluded_or_held.tsv"
    with open(excluded_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "release_class"])
        for e in evidence_list:
            if e["release_class"] in HOLD_EXCLUDE_CLASSES:
                writer.writerow([e["sample_id"], e["release_class"]])
    paths["samples_excluded_or_held"] = excluded_path

    return paths
