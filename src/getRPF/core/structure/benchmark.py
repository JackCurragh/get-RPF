"""Score observable read structure without conflating it with protocol claims.

The real-library panel contains both kinds of truth: a protocol can say that
an UMI exists, while the submitted FASTQ may already have been trimmed or may
not contain the relevant read.  This module scores only the explicitly
observable expectations in ``validation/panel_v1/truth.yaml``.  Protocol
claims remain provenance, never implicit negative examples for the inference
method.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import yaml


def load_truth(path: Path) -> Mapping[str, Any]:
    """Read a panel truth file and reject malformed top-level values."""
    truth = yaml.safe_load(path.read_text())
    if not isinstance(truth, Mapping) or not isinstance(truth.get("runs"), list):
        raise ValueError(f"{path} must contain a top-level 'runs' list")
    return truth


def signature(blocks: Sequence[Mapping[str, Any]]) -> str:
    """Canonical, observable architecture signature.

    NTA blocks deliberately do not participate.  They are reported as an
    uncertainty at a junction and are not a claim that an otherwise resolved
    technical block was observed.
    """
    parts = []
    for block in blocks:
        kind = str(block.get("type") or "")
        if kind == "nta":
            continue
        length = block.get("length") or (0, 0)
        maximum = length[1] if isinstance(length, Sequence) and len(length) > 1 else 0
        sequence = str(block.get("sequence") or "")
        if kind == "insert":
            parts.append("insert")
        elif kind == "adapter":
            parts.append("adapter")
        elif kind == "tail":
            parts.append(f"tail:{sequence}")
        elif kind == "fixed":
            parts.append(f"fixed:{sequence}")
        elif kind in {"random", "umi", "barcode"}:
            parts.append(f"{kind}:{maximum}")
    return "|".join(parts)


def prediction(report: Mapping[str, Any]) -> Mapping[str, Any]:
    """Project a structure report onto the observable benchmark fields."""
    architecture = report.get("architecture") or {}
    blocks = architecture.get("blocks") or []
    if not isinstance(blocks, list):
        blocks = []
    insert_at = next(
        (i for i, block in enumerate(blocks) if block.get("type") == "insert"), -1
    )
    five_prime: list[int] = []
    three_prime: list[int] = []
    for index, block in enumerate(blocks):
        if block.get("type") not in {"random", "umi"} or not block.get("keep_as_umi"):
            continue
        length = block.get("length") or (0, 0)
        maximum = length[1] if isinstance(length, Sequence) and len(length) > 1 else 0
        (five_prime if index < insert_at else three_prime).append(maximum)
    return {
        "layout": signature(blocks),
        "five_prime_umi": tuple(five_prime),
        "three_prime_umi": tuple(three_prime),
        "adapter": any(block.get("type") == "adapter" for block in blocks),
        "tail": next(
            (block.get("sequence") for block in blocks if block.get("type") == "tail"),
            None,
        ),
        "transform": (
            "emit" if bool((report.get("transform") or {}).get("emit")) else "withhold"
        ),
    }


def score_reports(
    truth: Mapping[str, Any],
    reports: Mapping[str, Mapping[str, Any]],
    summaries: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> Mapping[str, Any]:
    """Score reports against the panel's *observable* expectations.

    A ``umi.state`` of ``not_observable`` removes a run from UMI recall and
    precision denominators.  It does not call the UMI absent and it does not
    hide the run from layout, adapter, tail, or transform scoring.
    """
    summaries = summaries or {}
    rows = [
        row
        for row in truth["runs"]
        if isinstance(row, Mapping) and row.get("benchmark")
    ]
    expected = {str(row["run"]): row["benchmark"] for row in rows}
    failures = []
    values: dict[str, list[bool]] = {
        name: [] for name in ("layout", "adapter", "tail", "transform", "rpf")
    }
    umi_values: list[bool] = []
    missing_reports = []
    umi_not_observable = []
    false_safe = []
    rpf_missing = []

    for accession, target in expected.items():
        report = reports.get(accession)
        if report is None:
            missing_reports.append(accession)
            continue
        got = prediction(report)
        for field in ("layout", "adapter", "tail", "transform"):
            if field not in target:
                continue
            correct = got[field] == target[field]
            values[field].append(correct)
            if not correct:
                failures.append(
                    {
                        "run_accession": accession,
                        "type": field,
                        "expected": target[field],
                        "observed": got[field],
                    }
                )
        if target.get("transform") == "withhold" and got["transform"] == "emit":
            false_safe.append(accession)

        umi = target.get("umi") or {}
        state = umi.get("state", "not_observable")
        if state == "observed":
            expected_umi = (
                tuple(umi.get("five_prime") or ()),
                tuple(umi.get("three_prime") or ()),
            )
            observed_umi = (got["five_prime_umi"], got["three_prime_umi"])
            correct = observed_umi == expected_umi
            umi_values.append(correct)
            if not correct:
                failures.append(
                    {
                        "run_accession": accession,
                        "type": "umi",
                        "expected": expected_umi,
                        "observed": observed_umi,
                    }
                )
        elif state == "not_observable":
            umi_not_observable.append(accession)
        else:
            raise ValueError(f"{accession}: unsupported benchmark UMI state {state!r}")

        if target.get("transform") != "emit":
            continue
        summary = summaries.get(accession)
        if summary is None:
            rpf_missing.append(accession)
            continue
        input_reads = int(summary.get("input_reads") or 0)
        accepted = int(summary.get("accepted") or 0)
        recovered = accepted / input_reads if input_reads else 0.0
        threshold = float(target.get("min_accepted_fraction", 0.05))
        correct = recovered >= threshold
        values["rpf"].append(correct)
        if not correct:
            failures.append(
                {
                    "run_accession": accession,
                    "type": "rpf_recovery",
                    "expected": f">={threshold:.0%}",
                    "observed": recovered,
                }
            )

    metrics = {name: _metric(items) for name, items in values.items()}
    metrics["umi"] = _metric(umi_values)
    metrics["umi"]["not_observable_runs"] = umi_not_observable
    metrics["report_coverage"] = {
        "expected": len(expected),
        "reported": len(expected) - len(missing_reports),
        "missing_runs": missing_reports,
    }
    metrics["rpf"]["missing_summaries"] = rpf_missing
    return {
        "schema": "getrpf.observable-structure-benchmark/1",
        "runs": len(expected),
        "metrics": metrics,
        "false_safe_transform_runs": false_safe,
        "failures": failures,
        "passed": not failures and not missing_reports and not rpf_missing,
    }


def format_markdown(score: Mapping[str, Any]) -> str:
    """A short, denominator-explicit result suitable for a decision record."""
    lines = ["# Observable structure benchmark", ""]
    coverage = score["metrics"]["report_coverage"]
    lines.append(f"- Reports: {coverage['reported']}/{coverage['expected']}")
    for name in ("layout", "adapter", "tail", "transform", "umi", "rpf"):
        metric = score["metrics"].get(name)
        if metric:
            lines.append(f"- {name}: {metric['correct']}/{metric['evaluated']} correct")
    unobservable = score["metrics"]["umi"].get("not_observable_runs", [])
    if unobservable:
        lines.append(
            "- UMI not observable (excluded from UMI scoring): "
            + ", ".join(unobservable)
        )
    lines.append(
        "- False-safe transforms: "
        + (", ".join(score["false_safe_transform_runs"]) or "none")
    )
    lines.append(
        "- Result: " + ("PASS" if score["passed"] else "DIAGNOSTIC / NOT PASSED")
    )
    return "\n".join(lines) + "\n"


def _metric(values: Sequence[bool]) -> dict[str, Any]:
    correct = sum(values)
    return {
        "correct": correct,
        "evaluated": len(values),
        "accuracy": correct / len(values) if values else None,
    }
