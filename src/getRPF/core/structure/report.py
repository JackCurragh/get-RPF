"""Structure reports: machine-readable JSON and a readable text summary.

The JSON keeps every answer's evidence and rejected alternatives (spec §4.1)
and the per-position agreement profiles behind them. Per-read arrays are left
out; they belong to extraction, not to the report.
"""

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .assemble import StructureInference
from .pileup import FrameProfile, Pileup

SCHEMA = "getrpf.structure/1"
_PER_READ = {"adapter_starts", "insert_ends"}


def to_dict(result: StructureInference, sample: str) -> Dict[str, Any]:
    architecture = result.architecture
    return {
        "schema": SCHEMA,
        "sample": sample,
        "observation": _plain(result.observation),
        "window": result.window,
        "answers": {
            answer.question: _plain(answer)
            for answer in (result.q1, result.q2, result.q3, result.q4)
        },
        "architecture": (
            None
            if architecture is None
            else {"describe": architecture.describe(), **_plain(architecture)}
        ),
        "transform": _plain(result.transform),
        "profiles": {
            "first_pass": _profiles(result.junctions.first_pass),
            "final": _profiles(result.junctions.pileup),
        },
    }


def format_text(result: StructureInference, sample: str) -> str:
    observation = result.observation
    transform = result.transform
    lines = [
        f"get-RPF structure inference: {sample}",
        f"Input: {observation.input_state}, modal length {observation.modal_length} nt "
        f"({observation.modal_fraction:.0%} of {observation.reads} reads); "
        f"read-name UMI: {observation.header_umi or 'none'}",
        "Architecture: "
        + (result.architecture.describe() if result.architecture else "not built"),
        f"Transform: {'emitted' if transform.emit else 'withheld'} "
        f"(bounded by {transform.bound.value})",
    ]
    lines += [f"  withheld because {reason}" for reason in transform.reasons]
    lines += [f"  convention: {convention}" for convention in transform.conventions]
    for answer in (result.q1, result.q2, result.q3, result.q4):
        lines.append("")
        lines.append(f"{answer.question} [{answer.status.value}] {answer.explanation}")
        for alternative in answer.alternatives:
            lines.append(
                f"    {alternative.value}: {alternative.verdict} - {alternative.reason}"
            )
    return "\n".join(lines) + "\n"


def write_report(
    result: StructureInference, output_dir: Path, sample: str
) -> Tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{sample}.structure.json"
    text_path = output_dir / f"{sample}.structure.txt"
    json_path.write_text(json.dumps(to_dict(result, sample), indent=2) + "\n")
    text_path.write_text(format_text(result, sample))
    return json_path, text_path


def _profiles(pileup: Pileup) -> Dict[str, Any]:
    return {
        "groups": pileup.groups,
        "reads_grouped": pileup.reads_grouped,
        "reads_used": pileup.reads_used,
        "consensus_trim": list(pileup.consensus_trim),
        "read_start": _frame(pileup.read_start),
        "anchor": _frame(pileup.anchor),
    }


def _frame(profile: FrameProfile) -> List[Dict[str, Any]]:
    return [
        {
            "position": stats.position,
            "agreement": None if stats.agreement is None else round(stats.agreement, 4),
            "scored": stats.scored,
            "dominant_base": max(stats.composition, key=lambda b: stats.composition[b]),
            "dominant_fraction": round(stats.dominant_fraction, 4),
            "entropy_bits": round(stats.entropy, 4),
        }
        for stats in profile.positions
    ]


def _plain(value: Any) -> Any:
    """Dataclasses, enums and tuples as JSON-ready values."""
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _plain(getattr(value, field.name))
            for field in fields(value)
            if field.name not in _PER_READ
        }
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, float):
        return round(value, 4)
    return value
