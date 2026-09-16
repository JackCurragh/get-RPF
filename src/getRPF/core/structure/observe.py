"""Steps 0 and 1: what was sequenced, and what kind of input it is (spec §5.0-5.1).

Nothing here trims or decides structure. It records which information
channels exist and sets what later steps may assume.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import IO, List, Optional, Sequence, Tuple

from .config import InferenceConfig
from .model import Evidence

_COLON_UMI = re.compile(r"^[ACGTN+]{4,}$")
_SUFFIX_UMI = re.compile(r"_[ACGTN]{4,}$")


@dataclass(frozen=True)
class Observation:
    reads: int
    modal_length: int
    modal_fraction: float
    length_range: Tuple[int, int]
    """1st and 99th percentile read length."""
    input_state: str
    """fixed_length_footprint_candidate | raw_fixed_length | trimmed | mixed"""
    header_umi: Optional[str]
    """None, or the header convention carrying a UMI: colon_field | underscore_suffix."""
    terminal_poly_g: float
    """Fraction of reads ending in ten or more G (no-signal cycles on
    two-colour instruments)."""
    evidence: Tuple[Evidence, ...]


def read_fastq(path: Path, limit: int) -> Tuple[List[str], List[str], List[str]]:
    """The first ``limit`` records as (headers, sequences, qualities)."""
    headers: List[str] = []
    reads: List[str] = []
    quals: List[str] = []
    with _open_text(path) as handle:
        while len(reads) < limit:
            header = handle.readline()
            if not header:
                break
            reads.append(handle.readline().rstrip("\n"))
            handle.readline()
            quals.append(handle.readline().rstrip("\n"))
            headers.append(header.rstrip("\n"))
    return headers, reads, quals


def observe(
    reads: Sequence[str],
    headers: Sequence[str] = (),
    config: Optional[InferenceConfig] = None,
) -> Observation:
    config = config or InferenceConfig()
    lengths = Counter(len(read) for read in reads)
    total = len(reads)
    modal_length, modal_count = lengths.most_common(1)[0]
    modal_fraction = modal_count / total
    length_range = (_percentile(lengths, 0.01), _percentile(lengths, 0.99))
    if modal_fraction >= config.raw_modal_fraction:
        # A 35-cycle library is short enough that a normal 20--40 nt
        # footprint plus the conservative 10 nt inference overlap cannot be
        # observed in most reads.  Keep this as an uncertainty state: length
        # alone must never select zero trim or a UMI.
        if modal_length < config.typical_footprint + config.min_partial_overlap:
            input_state = "fixed_length_footprint_candidate"
        else:
            input_state = "raw_fixed_length"
    elif modal_fraction <= config.trimmed_modal_fraction:
        input_state = "trimmed"
    else:
        input_state = "mixed"
    header_umi = _header_umi(headers[:1000])
    poly_g = sum(1 for read in reads if read.endswith("G" * 10)) / total
    evidence = (
        Evidence(
            "profile",
            "modal_length_fraction",
            round(modal_fraction, 4),
            total,
            note=f"modal length {modal_length} nt; 1-99% range "
            f"{length_range[0]}-{length_range[1]} nt",
        ),
        Evidence("profile", "input_state", input_state, total),
        Evidence(
            "profile",
            "header_umi",
            header_umi,
            min(len(headers), 1000),
            note="read-name convention carrying a UMI, if any",
        ),
        Evidence("profile", "terminal_poly_g", round(poly_g, 4), total),
    )
    return Observation(
        reads=total,
        modal_length=modal_length,
        modal_fraction=modal_fraction,
        length_range=length_range,
        input_state=input_state,
        header_umi=header_umi,
        terminal_poly_g=poly_g,
        evidence=evidence,
    )


def _open_text(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _percentile(lengths: Counter[int], q: float) -> int:
    total = sum(lengths.values())
    running = 0
    for length in sorted(lengths):
        running += lengths[length]
        if running >= q * total:
            return length
    return max(lengths)


def _header_umi(headers: Sequence[str]) -> Optional[str]:
    """Detect a UMI carried in read names (bcl-convert or umi_tools style).

    The Illumina comment field ``1:N:0:ACGTACGT`` holds the sample index, not
    a UMI, and has too few colon fields to match here.
    """
    if not headers:
        return None
    kinds: Counter[str] = Counter()
    for header in headers:
        for token in header.lstrip("@").split():
            if token.count(":") >= 7 and _COLON_UMI.match(token.rsplit(":", 1)[1]):
                kinds["colon_field"] += 1
                break
            if _SUFFIX_UMI.search(token):
                kinds["underscore_suffix"] += 1
                break
    if not kinds:
        return None
    kind, count = kinds.most_common(1)[0]
    return kind if count >= 0.9 * len(headers) else None
