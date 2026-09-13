#!/usr/bin/env python3
"""Run the observable structure panel through infer-structure's transform.

This deliberately does not call the legacy ``extract-rpf`` workflow.  Each
run first produces a resolved architecture, then audit and (only when
explicitly requested) production apply that exact architecture.  A withheld
transform or zero accepted audit is a recorded hold, never a fallback trim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from getRPF.core.structure.assemble import infer_structure
from getRPF.core.structure.benchmark import format_markdown, load_truth, score_reports
from getRPF.core.structure.config import InferenceConfig
from getRPF.core.structure.observe import read_fastq
from getRPF.core.structure.report import to_dict, write_report
from getRPF.core.structure.transform import extract_reads


def _write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def _input(data_dir: Path, accession: str) -> Path | None:
    for suffix in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        path = data_dir / f"{accession}{suffix}"
        if path.exists():
            return path
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--truth", type=Path, default=Path("validation/panel_v1/truth.yaml")
    )
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=300_000)
    parser.add_argument("--max-reads", type=int, default=300_000)
    parser.add_argument(
        "--approve-production",
        action="store_true",
        help="write production FASTQs after each non-empty audit agrees",
    )
    args = parser.parse_args()

    truth = load_truth(args.truth)
    reports_dir = args.output_dir / "reports"
    summaries_dir = args.output_dir / "summaries"
    manifests_dir = args.output_dir / "runs"
    reports: dict[str, dict[str, Any]] = {}
    summaries: dict[str, dict[str, Any]] = {}

    for row in truth["runs"]:
        accession = str(row["run"])
        manifest: dict[str, Any] = {
            "run": accession,
            "approved_production": args.approve_production,
        }
        input_path = _input(args.data_dir, accession)
        if input_path is None:
            manifest.update(
                {"status": "missing_input", "hold_reason": "panel FASTQ not found"}
            )
            _write(manifests_dir / f"{accession}.json", manifest)
            continue

        headers, reads, _ = read_fastq(input_path, args.reads)
        if not reads:
            manifest.update({"status": "hold", "hold_reason": "empty FASTQ"})
            _write(manifests_dir / f"{accession}.json", manifest)
            continue
        result = infer_structure(
            reads, headers, InferenceConfig(sample_reads=args.reads)
        )
        write_report(result, reports_dir, accession)
        report = to_dict(result, accession)
        reports[accession] = report
        manifest["transform"] = report["transform"]
        manifest["architecture"] = report["architecture"]
        if not result.transform.emit or result.architecture is None:
            manifest.update({"status": "hold", "hold_reason": "transform withheld"})
            _write(manifests_dir / f"{accession}.json", manifest)
            continue

        audit = extract_reads(
            input_path, None, result.architecture, limit=args.max_reads
        )
        audit_payload = audit.to_dict()
        _write(summaries_dir / f"{accession}.rpf.fastq.gz.summary.json", audit_payload)
        summaries[accession] = audit_payload
        manifest["audit"] = audit_payload
        if audit.accepted == 0:
            manifest.update(
                {"status": "hold", "hold_reason": "zero accepted reads in audit"}
            )
            _write(manifests_dir / f"{accession}.json", manifest)
            continue
        if not args.approve_production:
            manifest["status"] = "audit_complete"
            _write(manifests_dir / f"{accession}.json", manifest)
            continue

        output = args.output_dir / "production" / f"{accession}.rpf.fastq.gz"
        output.parent.mkdir(parents=True, exist_ok=True)
        production = extract_reads(
            input_path, output, result.architecture, limit=args.max_reads
        )
        if production.to_dict() != audit_payload:
            raise RuntimeError(f"{accession}: audit and production transforms disagree")
        manifest.update(
            {"status": "production_complete", "production": production.to_dict()}
        )
        _write(manifests_dir / f"{accession}.json", manifest)

    score = score_reports(truth, reports, summaries)
    minimum_validation_depth = InferenceConfig().sample_reads
    score["inference_reads"] = args.reads
    score["production_requested"] = args.approve_production
    if args.reads < minimum_validation_depth:
        # Small bounded runs are valuable for diagnosing stability, but the
        # real-panel expectations were calibrated at the configured depth.
        # Do not let an exploratory sample masquerade as a validation pass.
        score["validation_status"] = "diagnostic_exploratory_depth"
        score["validation_note"] = (
            f"{args.reads} reads is below the required {minimum_validation_depth}"
        )
    elif score["passed"]:
        score["validation_status"] = "audit_passed"
        score["validation_note"] = (
            "An audit pass is not a production promotion; production remains "
            "an explicit separate action."
        )
    else:
        score["validation_status"] = "diagnostic_not_passed"
    _write(args.output_dir / "benchmark.json", score)
    (args.output_dir / "benchmark.md").write_text(format_markdown(score))
    return 0 if score["passed"] and args.reads >= minimum_validation_depth else 3


if __name__ == "__main__":
    raise SystemExit(main())
