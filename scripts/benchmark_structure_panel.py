#!/usr/bin/env python3
"""Score inferred read structures against the observable panel contract."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from getRPF.core.structure.benchmark import format_markdown, load_truth, score_reports


def _reports(directory: Path, suffix: str) -> dict[str, dict]:
    result = {}
    for path in directory.glob(f"*{suffix}"):
        accession = path.name[: -len(suffix)]
        result[accession] = json.loads(path.read_text())
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--summaries-dir", type=Path)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--markdown-out", type=Path, required=True)
    args = parser.parse_args()

    reports = _reports(args.reports_dir, ".structure.json")
    summaries = (
        _reports(args.summaries_dir, ".rpf.fastq.gz.summary.json")
        if args.summaries_dir
        else {}
    )
    score = score_reports(load_truth(args.truth), reports, summaries)
    args.json_out.write_text(json.dumps(score, indent=2) + "\n")
    args.markdown_out.write_text(format_markdown(score))
    return 0 if score["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
