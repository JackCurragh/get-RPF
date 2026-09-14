"""The per-read transform (spec §7), the only code that removes bases, and the
length-based part of Q5 validation (spec §5.7).

Audit and production call the same ``Transform.apply``; audit simply writes no
FASTQ. Read names and the qualities of the retained bases are preserved, and
a kept UMI is appended to the read name as ``_<UMI>`` (the umi_tools
convention).
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any, Dict, Iterator, List, Optional, Tuple, Union

from .anchors import _run_before, locate
from .config import InferenceConfig
from .model import Answer, Architecture, Block, Evidence, Status
from .observe import _open_text

TRANSFORM_SCHEMA = "getrpf.structure.transform/3"
"""Bump whenever ``Transform.apply`` semantics change.

Architecture and configuration alone are not a sufficient transform identity:
the tail rule is executable semantics.  Including this marker prevents an
audit made by an older implementation being presented as identical to a later
production run with the same architecture file.
"""


@dataclass(frozen=True)
class Accepted:
    insert: str
    quality: str
    umi: str
    boundary: str


@dataclass(frozen=True)
class Rejected:
    reason: str
    insert_length: Optional[int] = None
    boundary: Optional[str] = None


Outcome = Union[Accepted, Rejected]


def transform_hash(architecture: Architecture, config: InferenceConfig) -> str:
    """Identity of the per-read transform: architecture, policy and anchoring."""
    payload = {
        "schema": TRANSFORM_SCHEMA,
        "blocks": [
            [b.type, b.frame, list(b.length), b.sequence, b.keep_as_umi, b.remove]
            for b in architecture.blocks
        ],
        "fragment_policy": architecture.fragment_policy,
        "anchoring": [
            config.anchor_seed,
            config.min_insert,
            config.min_partial_overlap,
            config.max_mismatch_per_10nt,
            config.transform_min_overlap,
            config.fixed_max_mismatches,
        ],
    }
    encoded = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()[:16]


class Transform:
    """Apply one architecture to reads."""

    def __init__(
        self, architecture: Architecture, config: Optional[InferenceConfig] = None
    ) -> None:
        self.architecture = architecture
        self.config = config or InferenceConfig()
        blocks = architecture.blocks
        inserts = [i for i, block in enumerate(blocks) if block.type == "insert"]
        if len(inserts) != 1:
            raise ValueError("an architecture needs exactly one insert block")
        at = inserts[0]
        after = blocks[at + 1 :]
        self.five: Tuple[Block, ...] = blocks[:at]
        self.three: Tuple[Block, ...] = tuple(
            b for b in after if b.type not in ("adapter", "tail")
        )
        self.adapter = next(
            (b.sequence for b in after if b.type == "adapter" and b.sequence), None
        )
        self.tail = next(
            (b.sequence for b in after if b.type == "tail" and b.sequence), None
        )
        self.min_length, self.max_length = blocks[at].length
        self.hash = transform_hash(architecture, self.config)

    def apply(self, read: str, quality: str) -> Outcome:
        config = self.config
        end = len(read)
        boundary = "read_end"
        if self.adapter is not None:
            anchor = self._anchor(read)
            if anchor is None:
                return Rejected("anchor_not_found")
            if anchor < config.min_insert:
                return Rejected("adapter_dimer")
            end = anchor
            boundary = "adapter"
            if self.tail:
                # Spec §6: once a tail architecture is established, the insert
                # ends at the first base of the run before the adapter -- the
                # same rule inference used for the insert-side anchor. Inserts
                # that genuinely end in the tail base lose those bases, which
                # is why the boundary is an interval by convention.
                tail_length = _run_before(read, anchor, self.tail)
                if tail_length:
                    end -= tail_length
                    boundary = f"tail_{self.tail}"

        start = 0
        umi5: List[str] = []
        for block in self.five:
            if not block.remove:
                continue  # kept in the insert (non-templated-like, v1)
            segment = read[start : start + block.length[1]]
            if not self._fixed_ok(block, segment):
                return Rejected("fixed_mismatch")
            if block.keep_as_umi:
                umi5.append(segment)
            start += block.length[1]

        umi3: List[str] = []
        for block in reversed(self.three):
            if not block.remove:
                continue
            segment = read[max(0, end - block.length[1]) : end]
            if not self._fixed_ok(block, segment):
                return Rejected("fixed_mismatch")
            if block.keep_as_umi:
                umi3.insert(0, segment)
            end -= block.length[1]

        size = end - start
        if size <= 0:
            return Rejected("no_insert", size, boundary)
        if size < self.min_length:
            return Rejected("shorter_than_policy", size, boundary)
        if size > self.max_length:
            return Rejected("longer_than_policy", size, boundary)
        return Accepted(
            read[start:end], quality[start:end], "".join(umi5 + umi3), boundary
        )

    def _anchor(self, read: str) -> Optional[int]:
        """Adapter start: a full or long partial match, else a short exact
        adapter prefix at the read end (at least ``transform_min_overlap``)."""
        assert self.adapter is not None
        hit = locate(read, self.adapter, self.config)
        if hit is not None:
            return hit.start
        longest = min(len(self.adapter), self.config.min_partial_overlap - 1)
        for k in range(longest, self.config.transform_min_overlap - 1, -1):
            if read.endswith(self.adapter[:k]):
                return len(read) - k
        return None

    def _fixed_ok(self, block: Block, segment: str) -> bool:
        if block.type != "fixed" or not block.sequence:
            return True
        mismatches = sum(1 for a, b in zip(segment, block.sequence) if a != b)
        mismatches += abs(len(segment) - len(block.sequence))
        return mismatches <= self.config.fixed_max_mismatches


@dataclass
class ExtractionSummary:
    transform_hash: str
    architecture: str
    fragment_policy: str
    input_reads: int = 0
    accepted: int = 0
    rejected: Counter[str] = field(default_factory=Counter)
    lengths: Counter[int] = field(default_factory=Counter)
    candidate_lengths: Counter[int] = field(default_factory=Counter)
    boundaries: Counter[str] = field(default_factory=Counter)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "transform_schema": TRANSFORM_SCHEMA,
            "transform_hash": self.transform_hash,
            "architecture": self.architecture,
            "fragment_policy": self.fragment_policy,
            "input_reads": self.input_reads,
            "accepted": self.accepted,
            "accepted_fraction": (
                round(self.accepted / self.input_reads, 4) if self.input_reads else 0.0
            ),
            "rejected": dict(sorted(self.rejected.items())),
            "insert_lengths": {str(k): v for k, v in sorted(self.lengths.items())},
            # Includes otherwise valid, anchored inserts outside the fragment
            # policy.  This distinguishes a bad boundary from a genuinely
            # long-fragment library without weakening the closed-fail gate.
            "candidate_insert_lengths": {
                str(k): v for k, v in sorted(self.candidate_lengths.items())
            },
            "boundaries": dict(sorted(self.boundaries.items())),
        }


def iter_fastq(path: Path) -> Iterator[Tuple[str, str, str]]:
    """Stream (header, sequence, quality) records."""
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip("\n")
            handle.readline()
            quality = handle.readline().rstrip("\n")
            yield header.rstrip("\n"), sequence, quality


def extract_reads(
    input_path: Path,
    output_path: Optional[Path],
    architecture: Architecture,
    config: Optional[InferenceConfig] = None,
    skip: int = 0,
    limit: Optional[int] = None,
    collect: Optional[Counter[str]] = None,
) -> ExtractionSummary:
    """Apply the transform to a FASTQ. ``output_path=None`` is an audit;
    ``collect`` counts the accepted inserts (for collapsed output)."""
    transform = Transform(architecture, config)
    summary = ExtractionSummary(
        transform.hash, architecture.describe(), architecture.fragment_policy
    )
    handle: Optional[IO[str]] = None
    if output_path is not None:
        if str(output_path).endswith(".gz"):
            handle = gzip.open(output_path, "wt")
        else:
            handle = open(output_path, "w")
    try:
        for index, (header, sequence, quality) in enumerate(iter_fastq(input_path)):
            if index < skip:
                continue
            if limit is not None and summary.input_reads >= limit:
                break
            summary.input_reads += 1
            outcome = transform.apply(sequence, quality)
            if isinstance(outcome, Rejected):
                summary.rejected[outcome.reason] += 1
                if outcome.insert_length is not None and outcome.insert_length > 0:
                    summary.candidate_lengths[outcome.insert_length] += 1
                if outcome.boundary is not None:
                    summary.boundaries[outcome.boundary] += 1
                continue
            summary.accepted += 1
            summary.lengths[len(outcome.insert)] += 1
            summary.candidate_lengths[len(outcome.insert)] += 1
            summary.boundaries[outcome.boundary] += 1
            if collect is not None:
                collect[outcome.insert] += 1
            if handle is not None:
                handle.write(
                    f"{_named(header, outcome.umi)}\n{outcome.insert}\n+\n"
                    f"{outcome.quality}\n"
                )
    finally:
        if handle is not None:
            handle.close()
    return summary


def validate_extraction(
    summary: ExtractionSummary, config: Optional[InferenceConfig] = None
) -> Answer:
    """Q5 from insert lengths alone (spec §5.7). Periodicity needs alignment."""
    config = config or InferenceConfig()
    total, accepted = summary.input_reads, summary.accepted
    fraction = accepted / total if total else 0.0
    core = sum(
        count
        for length, count in summary.lengths.items()
        if config.core_min <= length <= config.core_max
    )
    core_fraction = core / accepted if accepted else 0.0
    mode = summary.lengths.most_common(1)[0][0] if accepted else None
    too_long = summary.rejected.get("longer_than_policy", 0) / total if total else 0.0
    evidence = (
        Evidence("validation", "accepted_fraction", round(fraction, 4), total),
        Evidence("validation", "insert_mode", mode, accepted),
        Evidence(
            "validation",
            "core_fraction",
            round(core_fraction, 4),
            accepted,
            note=f"inserts of {config.core_min}-{config.core_max} nt",
        ),
        Evidence(
            "validation", "longer_than_policy_fraction", round(too_long, 4), total
        ),
    )
    if total < config.min_validation_reads:
        return Answer(
            "Q5",
            "insufficient_evidence",
            Status.UNDERPOWERED,
            evidence,
            (),
            f"Only {total} reads were transformed "
            f"(need {config.min_validation_reads}).",
        )
    if fraction < config.min_accepted_fraction:
        identity = "rna_like" if too_long >= 0.5 else "technical_failure"
    elif (
        core_fraction >= config.core_fraction_likely
        and mode is not None
        and config.core_min <= mode <= config.core_max
    ):
        identity = "riboseq_likely"
    else:
        identity = "riboseq_possible"
    return Answer(
        "Q5",
        identity,
        Status.RESOLVED,
        evidence,
        (),
        f"{identity}: {fraction:.0%} of reads give an insert within "
        f"{summary.fragment_policy}; mode {mode} nt; {core_fraction:.0%} of "
        f"inserts are {config.core_min}-{config.core_max} nt. Length evidence "
        "only; frame periodicity needs alignment.",
    )


def _named(header: str, umi: str) -> str:
    if not umi:
        return header
    name, _, rest = header.partition(" ")
    return f"{name}_{umi}" + (f" {rest}" if rest else "")
