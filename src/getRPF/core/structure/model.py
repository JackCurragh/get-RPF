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
