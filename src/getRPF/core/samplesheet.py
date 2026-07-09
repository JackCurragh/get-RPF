"""M6: local-FASTQ samplesheet cohort driver, and the audit/override apply
modes shared with the single-sample `extract` command.

See docs/release_qc_and_terminal_trimming_plan.md sections 10-11 and M6.
"""

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from .pipeline import run_sample
from .release import ALL_RELEASE_CLASSES, HOLD_EXCLUDE_CLASSES, REVIEW_CLASSES, write_cohort_tsvs

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = {"sample_id", "fastq_1"}


@dataclass
class SampleSheetRow:
    sample_id: str
    fastq_1: Path
    # organism/study_id are carried into the evidence object's protocol
    # field (see run_cohort) for cohort tracking. fastq_2/layout/fastqc_dir
    # are accepted for compatibility with the full samplesheet spec
    # (section 10) but not yet consumed: paired-end extraction and
    # FastQC-derived protocol evidence aren't implemented.
    fastq_2: Optional[Path] = None
    layout: Optional[str] = None
    organism: Optional[str] = None
    study_id: Optional[str] = None
    protocol_hint: Optional[str] = None
    fastqc_dir: Optional[Path] = None


@dataclass
class CohortRunResult:
    evidence: List[Dict] = field(default_factory=list)
    cohort_tsv_paths: Dict[str, Path] = field(default_factory=dict)

    @property
    def release_classes(self) -> List[str]:
        return [e["release_class"] for e in self.evidence]

    def exit_code(self, fail_on: str = "none") -> int:
        return compute_exit_code(self.release_classes, fail_on=fail_on)


def parse_samplesheet(path: Path) -> List[SampleSheetRow]:
    """Parse a local-FASTQ samplesheet. Minimum columns: sample_id,fastq_1."""
    path = Path(path)
    rows: List[SampleSheetRow] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Samplesheet missing required column(s): {sorted(missing)}")

        for record in reader:
            rows.append(
                SampleSheetRow(
                    sample_id=record["sample_id"],
                    fastq_1=Path(record["fastq_1"]),
                    fastq_2=Path(record["fastq_2"]) if record.get("fastq_2") else None,
                    layout=record.get("layout") or None,
                    organism=record.get("organism") or None,
                    study_id=record.get("study_id") or None,
                    protocol_hint=record.get("protocol_hint") or None,
                    fastqc_dir=Path(record["fastqc_dir"]) if record.get("fastqc_dir") else None,
                )
            )
    return rows


def _infer_format(fastq_path: Path) -> str:
    suffixes = "".join(fastq_path.suffixes).lower()
    if "fasta" in suffixes or suffixes.endswith((".fa", ".fa.gz")):
        return "fasta"
    return "fastq"


def compute_exit_code(release_classes: List[str], fail_on: str = "none") -> int:
    """Section 11's process exit-code mirror.

    `fail_on` sets the severity threshold that causes a nonzero exit:
    "none" (default) never fails the process -- review/hold are correct,
    expected refusal outcomes, not tool errors. "review" fails on any
    review-or-worse sample; "hold" fails only on hold/exclude samples
    (review is treated as acceptable, per the worked example in section 11).
    """
    worst = 0
    for c in release_classes:
        if c not in ALL_RELEASE_CLASSES:
            raise ValueError(f"Unknown release class: {c!r}")
        if c in HOLD_EXCLUDE_CLASSES:
            worst = max(worst, 20)
        elif c in REVIEW_CLASSES:
            worst = max(worst, 10)

    if fail_on == "none":
        return 0
    if fail_on == "review":
        return worst
    if fail_on == "hold":
        return worst if worst >= 20 else 0
    return 0


def run_cohort(
    samplesheet: Path,
    output_dir: Path,
    audit_only: bool = False,
    rules_path: Optional[Path] = None,
    infer_reads: int = 500_000,
    max_reads: Optional[int] = None,
    collapse: bool = True,
) -> CohortRunResult:
    """Run the M1-M5 pipeline over every row of a local-FASTQ samplesheet.

    audit_only=True: infer + screen, apply nothing (rule development /
        new-family validation).
    rules_path: override mode -- apply explicit frozen rules (same shape
        as a `.trim_rules.applied.json`-derived override_trims dict,
        keyed by sample_id) instead of inferring them.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = parse_samplesheet(samplesheet)

    all_rules: Dict[str, Dict[int, Dict[str, int]]] = {}
    if rules_path is not None:
        all_rules = json.loads(Path(rules_path).read_text())

    result = CohortRunResult()

    for row in rows:
        # A single-suffix base path: downstream artifacts are derived via
        # Path.with_suffix (e.g. .collapsed.fa, .evidence.json), which only
        # replaces the *last* suffix -- an embedded dot (e.g. "id.rpfs.fastq")
        # would silently misname them.
        sample_output = output_dir / f"{row.sample_id}.fastq"
        format_ = _infer_format(row.fastq_1)

        override_trims = None
        if row.sample_id in all_rules:
            override_trims = {
                int(length): rule for length, rule in all_rules[row.sample_id].items()
            }

        protocol: Optional[Dict] = None
        if row.protocol_hint or row.organism or row.study_id:
            protocol = {}
            if row.protocol_hint:
                protocol["source"] = "protocol_hint"
                protocol["name"] = row.protocol_hint
            if row.organism:
                protocol["organism"] = row.organism
            if row.study_id:
                protocol["study_id"] = row.study_id

        logger.info(f"Running sample {row.sample_id} ({row.fastq_1})...")
        evidence, _extraction_result = run_sample(
            input_file=row.fastq_1,
            output_file=sample_output,
            format=format_,
            sample_id=row.sample_id,
            infer_reads=infer_reads,
            max_reads=max_reads,
            apply_trims=not audit_only,
            override_trims=override_trims,
            collapse_output=collapse,
            protocol=protocol,
        )
        result.evidence.append(evidence)

    result.cohort_tsv_paths = write_cohort_tsvs(result.evidence, output_dir)
    return result
