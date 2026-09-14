"""Data model shared by every inference step (spec §8)."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Literal, Optional, Tuple

Frame = Literal["read_start", "anchor"]
"""Coordinate frame: positions counted inward from the read start, or back
from the anchor (the adapter start, or the read end for trimmed reads)."""


class Status(str, Enum):
    """How far a question is answered (spec §4)."""

    RESOLVED = "resolved"
    INTERVAL = "interval"
    AMBIGUOUS = "ambiguous"
    CONFLICTING = "conflicting"
    UNDERPOWERED = "underpowered"
    NOT_OBSERVABLE = "not_observable"


@dataclass(frozen=True)
class Evidence:
    """One measurement behind an answer. ``value`` is the raw metric."""

    source: str
    metric: str
    value: Any
    n: int
    note: str = ""


@dataclass(frozen=True)
class Alternative:
    """An answer that was considered, and why it was not chosen (spec §4.1)."""

    value: Any
    verdict: Literal["contradicted", "weakly_compatible", "untested"]
    reason: str


@dataclass(frozen=True)
class Answer:
    """The answer to one question, with its provenance."""

    question: Literal["Q1", "Q2", "Q3", "Q4", "Q5"]
    value: Any
    status: Status
    evidence: Tuple[Evidence, ...]
    alternatives: Tuple[Alternative, ...]
    explanation: str


@dataclass(frozen=True)
class JunctionCall:
    """What sits between a frame's edge and the insert: the value of Q2 or Q3.

    ``technical_length`` counts bases that are technical and would be removed,
    split into ``blocks`` of (kind, length) ordered from the frame edge inward.
    ``nta_length`` counts the junction bases after them whose agreement lies
    between the random and biology levels, which is what a non-templated
    addition in some reads looks like. v1 reports them and does not remove them.
    """

    frame: Frame
    technical_length: int
    blocks: Tuple[Tuple[str, int], ...]
    nta_length: int
    nta_rates: Tuple[Optional[float], ...]
    nta_bases: Dict[str, float]


@dataclass(frozen=True)
class Block:
    """One element of the read, in read order (spec §8).

    ``remove`` says whether the transform cuts the block out of the emitted
    insert; ``keep_as_umi`` whether the removed bases are carried as a UMI.
    """

    type: Literal[
        "umi", "random", "fixed", "barcode", "nta", "tail", "insert", "adapter"
    ]
    frame: Frame
    length: Tuple[int, int]
    sequence: Optional[str]
    keep_as_umi: bool
    remove: bool
    status: Status
    note: str = ""
    evidence: Tuple[Evidence, ...] = ()


@dataclass(frozen=True)
class Architecture:
    """The ordered blocks of a read: the contract with extraction (spec §7.1)."""

    blocks: Tuple[Block, ...]
    source: Literal["inferred", "template", "template_completed"]
    status: Status
    fragment_policy: str

    def describe(self) -> str:
        """One line, e.g. ``[random 5 UMI][insert][adapter AGATCGGAAGAG...]``."""
        parts = []
        for block in self.blocks:
            low, high = block.length
            sequence = block.sequence or ""
            if block.type == "insert":
                parts.append("[insert]")
            elif block.type == "adapter":
                shown = sequence[:12] + ("..." if len(sequence) > 12 else "")
                parts.append(f"[adapter {shown}]")
            elif block.type == "tail":
                parts.append(f"[poly({sequence})]")
            elif block.type == "fixed":
                parts.append(f"[fixed {sequence}]")
            elif block.type == "nta":
                parts.append(
                    f"[nta {low}-{high} {'removed' if block.remove else 'kept'}]"
                )
            else:
                size = f"{low}" if low == high else f"{low}-{high}"
                parts.append(
                    f"[{block.type} {size}{' UMI' if block.keep_as_umi else ''}]"
                )
        return "".join(parts)


@dataclass(frozen=True)
class TransformDecision:
    """Whether a per-read transform may be emitted (spec §7.2)."""

    emit: bool
    bound: Status
    """The least-resolved status among the questions that change emitted bases."""
    reasons: Tuple[str, ...]
    """Why the transform is withheld; empty when it is emitted."""
    conventions: Tuple[str, ...]
    """Named conventions the emitted transform relies on."""
    flags: Tuple[str, ...] = ()
    """Findings a reviewer should see that never withhold the transform: junction
    bases kept in the insert although they look non-templated in most reads, an
    alignment check that disagrees about one junction base, or a check that could
    not run (spec §7.2)."""
