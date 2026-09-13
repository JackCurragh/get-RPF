"""Library-as-reference pileup and the Q2/Q3 junction calls (spec §5.3).

Reads that share an abundant k-mer come from the same biological fragment.
Lined up on that k-mer, their biological bases agree, random technical bases
(UMIs, randomised linker bases) agree only by chance, and non-templated
additions agree partially. Scoring each read position against the
leave-one-out consensus of the other distinct sequences in its group
therefore separates biology from technical sequence without a genome.

That within-fragment agreement is one axis of the spec §2 table. The other is
the base mix across the whole library at the same position: a position that
is constant across the library is fixed technical sequence, whatever its
agreement.

Positions are counted inward from two edges: the read start (Q2) and the
anchor (Q3) -- the adapter start, or the read end for reads that were already
trimmed. Reads whose anchor was not located are not used.

Inference runs in up to two passes. The first builds each group's consensus
from every base. If it finds technical blocks, a second pass rebuilds the
consensus without them, so that one read's UMI cannot blur another read's
insert boundary. Scoring is identical in both passes.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Sequence, Tuple

from .config import InferenceConfig
from .model import Alternative, Answer, Evidence, Frame, JunctionCall, Status

BASES = "ACGT"
CHANCE_AGREEMENT = 0.25
"""Agreement expected for a base unrelated to the consensus."""

_WHERE: Dict[str, str] = {
    "read_start": "from the read start",
    "anchor": "before the anchor",
}


@dataclass(frozen=True)
class PositionStats:
    """One position in one frame, 1-based, counted inward from the frame edge."""

    position: int
    scored: int
    agreement: Optional[float]
    composition: Dict[str, float]
    entropy: float
    dominant_fraction: float
    disagreeing: Dict[str, float]


@dataclass(frozen=True)
class FrameProfile:
    frame: Frame
    positions: Tuple[PositionStats, ...]


@dataclass(frozen=True)
class Pileup:
    groups: int
    reads_grouped: int
    reads_used: int
    consensus_trim: Tuple[int, int]
    read_start: FrameProfile
    anchor: FrameProfile


@dataclass(frozen=True)
class JunctionInference:
    q2: Answer
    q3: Answer
    pileup: Pileup
    first_pass: Pileup


@dataclass
class _Tally:
    scored: int = 0
    agree: int = 0
    disagreeing: Counter[str] = field(default_factory=Counter)


Group = List[Tuple[str, int]]
"""Distinct insert regions in one group, each with the position of its seed."""


def infer_junctions(
    reads: Sequence[str],
    anchors: Optional[Sequence[Optional[int]]] = None,
    config: Optional[InferenceConfig] = None,
) -> JunctionInference:
    """Answer Q2 and Q3 from the reads alone.

    ``anchors`` gives each read's adapter start (``None`` where no adapter was
    located; those reads are not used). Pass ``anchors=None`` for reads that
    contain no adapter at all, such as already-trimmed input; the read end is
    then the anchor.
    """
    config = config or InferenceConfig()
    regions = _regions(reads, anchors, config)
    groups, grouped = _group(regions, config)
    sample = regions[: config.few_values_sample]

    first = _pileup(regions, groups, grouped, config, 0, 0)
    q2 = call_junction(first.read_start, first, sample, config)
    q3 = call_junction(first.anchor, first, sample, config)
    trim5, trim3 = _technical_length(q2), _technical_length(q3)
    if not (trim5 or trim3):
        return JunctionInference(q2, q3, first, first)

    final = _pileup(regions, groups, grouped, config, trim5, trim3)
    return JunctionInference(
        call_junction(final.read_start, final, sample, config),
        call_junction(final.anchor, final, sample, config),
        final,
        first,
    )


def build_pileup(
    reads: Sequence[str],
    anchors: Optional[Sequence[Optional[int]]] = None,
    config: Optional[InferenceConfig] = None,
    consensus_trim: Tuple[int, int] = (0, 0),
) -> Pileup:
    """Measure within-fragment agreement and library composition per position."""
    config = config or InferenceConfig()
    regions = _regions(reads, anchors, config)
    groups, grouped = _group(regions, config)
    return _pileup(regions, groups, grouped, config, *consensus_trim)


def call_junction(
    profile: FrameProfile,
    pileup: Pileup,
    sample: Sequence[str],
    config: InferenceConfig,
) -> Answer:
    """Turn one frame's profile into a Q2 (read start) or Q3 (anchor) answer."""
    frame = profile.frame
    question: Literal["Q2", "Q3"] = "Q2" if frame == "read_start" else "Q3"
    where = _WHERE[frame]
    evidence = [
        Evidence(
            "pileup",
            "informative_groups",
            pileup.groups,
            pileup.reads_grouped,
            note=f"{pileup.reads_grouped} of {pileup.reads_used} reads grouped",
        )
    ]
    if pileup.consensus_trim != (0, 0):
        evidence.append(
            Evidence(
                "pileup",
                "consensus_excludes",
                pileup.consensus_trim,
                pileup.groups,
                note="second pass: consensus built without the 5' and 3' "
                "technical bases found in the first pass",
            )
        )

    if pileup.reads_used == 0:
        return Answer(
            question,
            None,
            Status.NOT_OBSERVABLE,
            tuple(evidence),
            (Alternative(None, "untested", "no read has a located anchor"),),
            f"No read has a located anchor, so the bases {where} cannot be examined.",
        )

    classes = [_classify(stats, config) for stats in profile.positions]
    if pileup.groups < config.min_groups or classes[0] == "unscored":
        return Answer(
            question,
            None,
            Status.UNDERPOWERED,
            tuple(evidence),
            (
                Alternative(
                    None,
                    "untested",
                    f"only {pileup.groups} informative groups "
                    f"(need {config.min_groups})",
                ),
            ),
            f"Too few repeated-fragment groups ({pileup.groups}, need "
            f"{config.min_groups}) to tell technical from biological bases "
            f"{where}; not called.",
        )

    technical = [i for i, cls in enumerate(classes) if cls in ("random", "fixed")]
    tech_len = technical[-1] + 1 if technical else 0
    first_bio = next(
        (i for i in range(tech_len, len(classes)) if classes[i] == "biology"), None
    )
    shown = len(classes) if first_bio is None else min(len(classes), first_bio + 2)
    evidence.extend(
        _position_evidence(profile.positions[i], classes[i]) for i in range(shown)
    )
    if first_bio is None:
        return Answer(
            question,
            None,
            Status.AMBIGUOUS,
            tuple(evidence),
            (
                Alternative(
                    tech_len,
                    "untested",
                    f"no biology-level agreement within {config.window} positions",
                ),
            ),
            f"No position within {config.window} nt {where} reaches "
            f"biology-level agreement (>= {config.agree_bio:.2f}); the technical "
            "block may be longer than the window.",
        )

    blocks = _technical_blocks(classes, tech_len, sample, frame, config)
    if blocks:
        evidence.append(Evidence("blocks", "technical_blocks", blocks, len(sample)))

    bio_level = statistics.median(
        stats.agreement
        for stats, cls in zip(profile.positions[first_bio:], classes[first_bio:])
        if cls == "biology" and stats.agreement is not None
    )
    zone = profile.positions[tech_len:first_bio]
    nta_rates = tuple(
        (
            None
            if stats.agreement is None
            else round(
                min(
                    1.0,
                    max(
                        0.0,
                        (bio_level - stats.agreement) / (bio_level - CHANCE_AGREEMENT),
                    ),
                ),
                3,
            )
        )
        for stats in zone
    )
    nta_bases = dict(zone[0].disagreeing) if zone else {}
    if zone:
        evidence.append(
            Evidence(
                "pileup",
                "nta_rate_estimate",
                nta_rates,
                zone[0].scored,
                note="(biology - agreement) / (biology - 0.25) per junction position",
            )
        )

    inside_bio = [i + 1 for i in range(tech_len) if classes[i] == "biology"]
    if inside_bio:
        status = Status.AMBIGUOUS
    elif zone:
        status = Status.INTERVAL
    else:
        status = Status.RESOLVED
    call = JunctionCall(frame, tech_len, blocks, len(zone), nta_rates, nta_bases)
    return Answer(
        question,
        call,
        status,
        tuple(evidence),
        _alternatives(profile, classes, tech_len, first_bio, config),
        _explain(call, profile, first_bio, bio_level, pileup.groups, inside_bio, where),
    )


def _regions(
    reads: Sequence[str],
    anchors: Optional[Sequence[Optional[int]]],
    config: InferenceConfig,
) -> List[str]:
    """The part of each read upstream of its anchor, within the bounded sample."""
    if anchors is None:
        return list(reads[: config.sample_reads])
    regions = [
        read[:anchor] for read, anchor in zip(reads, anchors) if anchor is not None
    ]
    return regions[: config.sample_reads]


def _is_low_complexity(kmer: str, config: InferenceConfig) -> bool:
    if max(Counter(kmer).values()) / len(kmer) > config.max_seed_base_fraction:
        return True
    run = 1
    for previous, current in zip(kmer, kmer[1:]):
        run = run + 1 if current == previous else 1
        if run > config.max_seed_homopolymer:
            return True
    return False


def _group(regions: Sequence[str], config: InferenceConfig) -> Tuple[List[Group], int]:
    """Assign each read to its highest-ranked seed; keep informative groups."""
    k = config.seed_k
    seed_counts: Counter[str] = Counter()
    for region in regions:
        seed_counts.update({region[p : p + k] for p in range(len(region) - k + 1)})
    rank: Dict[str, int] = {}
    for kmer, _count in seed_counts.most_common(config.top_seeds * 3):
        if len(rank) == config.top_seeds:
            break
        if not _is_low_complexity(kmer, config):
            rank[kmer] = len(rank)

    by_seed: Dict[str, Dict[str, int]] = defaultdict(dict)
    grouped = 0
    for region in regions:
        best: Optional[Tuple[int, int, str]] = None
        for p in range(len(region) - k + 1):
            kmer = region[p : p + k]
            seed_rank = rank.get(kmer)
            if seed_rank is not None and (best is None or seed_rank < best[0]):
                best = (seed_rank, p, kmer)
        if best is not None:
            # Identical regions collapse: PCR duplicates carry no information.
            by_seed[best[2]].setdefault(region, best[1])
            grouped += 1
    groups = [
        list(members.items())
        for members in by_seed.values()
        if len(members) >= config.min_group_distinct
    ]
    return groups, grouped


def _score(
    column: Counter[str],
    base: str,
    own_counted: bool,
    tally: _Tally,
    min_others: int,
) -> None:
    """Score one base against the leave-one-out consensus of its column."""
    own = 1 if own_counted else 0
    if sum(column.values()) - own < min_others:
        return
    best_other = max(
        count - (own if other == base else 0) for other, count in column.items()
    )
    tally.scored += 1
    if column[base] - own == best_other:
        tally.agree += 1
    else:
        tally.disagreeing[base] += 1


def _pileup(
    regions: Sequence[str],
    groups: Sequence[Group],
    grouped: int,
    config: InferenceConfig,
    trim5: int,
    trim3: int,
) -> Pileup:
    window = config.window
    start = [_Tally() for _ in range(window)]
    end = [_Tally() for _ in range(window)]
    for members in groups:
        columns: Dict[int, Counter[str]] = defaultdict(Counter)
        for region, seed_at in members:
            for j in range(trim5, len(region) - trim3):
                columns[j - seed_at][region[j]] += 1
        for region, seed_at in members:
            length = len(region)
            for i in range(min(window, length)):
                _score(
                    columns[i - seed_at],
                    region[i],
                    trim5 <= i < length - trim3,
                    start[i],
                    config.min_others,
                )
                j = length - 1 - i
                _score(
                    columns[j - seed_at],
                    region[j],
                    trim5 <= j < length - trim3,
                    end[i],
                    config.min_others,
                )

    start_counts: List[Counter[str]] = [Counter() for _ in range(window)]
    end_counts: List[Counter[str]] = [Counter() for _ in range(window)]
    for region in regions[: config.composition_sample]:
        for i in range(min(window, len(region))):
            start_counts[i][region[i]] += 1
            end_counts[i][region[-1 - i]] += 1

    return Pileup(
        groups=len(groups),
        reads_grouped=grouped,
        reads_used=len(regions),
        consensus_trim=(trim5, trim3),
        read_start=_frame_profile("read_start", start, start_counts, config),
        anchor=_frame_profile("anchor", end, end_counts, config),
    )


def _frame_profile(
    frame: Frame,
    tallies: Sequence[_Tally],
    counts: Sequence[Counter[str]],
    config: InferenceConfig,
) -> FrameProfile:
    positions = []
    for i, (tally, count) in enumerate(zip(tallies, counts)):
        total = sum(count[base] for base in BASES)
        composition = {base: (count[base] / total if total else 0.0) for base in BASES}
        disagreeing_total = sum(tally.disagreeing.values())
        positions.append(
            PositionStats(
                position=i + 1,
                scored=tally.scored,
                agreement=(
                    tally.agree / tally.scored
                    if tally.scored >= config.min_scored_per_position
                    else None
                ),
                composition=composition,
                entropy=-sum(p * math.log2(p) for p in composition.values() if p > 0),
                dominant_fraction=max(composition.values()),
                disagreeing=(
                    {
                        base: round(tally.disagreeing[base] / disagreeing_total, 3)
                        for base in BASES
                    }
                    if disagreeing_total
                    else {}
                ),
            )
        )
    return FrameProfile(frame, tuple(positions))


def _classify(stats: PositionStats, config: InferenceConfig) -> str:
    """Place a position in the spec §2 table from its two axes."""
    if stats.agreement is None:
        return "unscored"
    if stats.dominant_fraction >= config.fixed_base_fraction:
        return "fixed"
    if stats.agreement >= config.agree_bio:
        return "biology"
    if stats.agreement <= config.agree_random:
        return "random" if stats.entropy >= config.random_min_entropy else "skewed"
    return "intermediate"


def _technical_blocks(
    classes: Sequence[str],
    tech_len: int,
    sample: Sequence[str],
    frame: Frame,
    config: InferenceConfig,
) -> Tuple[Tuple[str, int], ...]:
    """Split the technical run into fixed, barcode and random blocks."""
    blocks: List[Tuple[str, int]] = []
    i = 0
    while i < tech_len:
        j = i
        if classes[i] == "fixed":
            while j < tech_len and classes[j] == "fixed":
                j += 1
            blocks.append(("fixed", j - i))
        else:
            while j < tech_len and classes[j] != "fixed":
                j += 1
            blocks.extend(_split_variable(i, j, sample, frame, config))
        i = j
    return tuple(blocks)


def _split_variable(
    start: int,
    end: int,
    sample: Sequence[str],
    frame: Frame,
    config: InferenceConfig,
) -> List[Tuple[str, int]]:
    """Find barcode-like segments (few distinct values) inside a variable run.

    Single positions cannot show this: four barcodes can look perfectly
    random base by base. Whole segments can, because a barcode segment takes
    a handful of values while a random one takes 4**length.
    """
    blocks: List[Tuple[str, int]] = []
    random_from = pos = start
    width = config.min_value_segment
    while pos + width <= end:
        needed = _values_needed(sample, frame, pos, pos + width, config)
        if needed > config.few_values_max:
            pos += 1
            continue
        stop = pos + width
        while (
            stop < end
            and _values_needed(sample, frame, pos, stop + 1, config)
            <= needed * config.few_values_growth
        ):
            stop += 1
        if random_from < pos:
            blocks.append(("random", pos - random_from))
        blocks.append(("barcode", stop - pos))
        random_from = pos = stop
    if random_from < end:
        blocks.append(("random", end - random_from))
    return blocks


def _values_needed(
    sample: Sequence[str],
    frame: Frame,
    start: int,
    end: int,
    config: InferenceConfig,
) -> int:
    """How many distinct segment values it takes to cover most reads."""
    values: Counter[str] = Counter()
    for region in sample:
        length = len(region)
        if length < end:
            continue
        if frame == "read_start":
            values[region[start:end]] += 1
        else:
            values[region[length - end : length - start]] += 1
    total = sum(values.values())
    if not total:
        return 1 << 30
    covered = 0
    for needed, (_value, count) in enumerate(values.most_common(), 1):
        covered += count
        if covered >= config.few_values_coverage * total:
            return needed
    return len(values)


def _technical_length(answer: Answer) -> int:
    if isinstance(answer.value, JunctionCall):
        return answer.value.technical_length
    return 0


def _fmt(value: Optional[float]) -> str:
    return "n/a" if value is None else f"{value:.2f}"


def _position_evidence(stats: PositionStats, cls: str) -> Evidence:
    return Evidence(
        "pileup",
        f"agreement_at_{stats.position}",
        None if stats.agreement is None else round(stats.agreement, 3),
        stats.scored,
        note=f"{cls}; library entropy {stats.entropy:.2f} bits; "
        f"dominant base {stats.dominant_fraction:.2f}",
    )


def _alternatives(
    profile: FrameProfile,
    classes: Sequence[str],
    tech_len: int,
    first_bio: int,
    config: InferenceConfig,
) -> Tuple[Alternative, ...]:
    """Every nearby boundary, with the reason it was not chosen (spec §4.1)."""
    candidates = sorted(
        {0, tech_len - 1, tech_len + 1, first_bio, first_bio + 1} - {tech_len}
    )
    alternatives = []
    for k in candidates:
        if k < 0 or k > len(classes):
            continue
        if k < tech_len:
            stats = profile.positions[k]
            alternatives.append(
                Alternative(
                    k,
                    "contradicted",
                    f"position {k + 1} would be kept as insert, but it is "
                    f"{classes[k]} (agreement {_fmt(stats.agreement)}, library "
                    f"entropy {stats.entropy:.2f} bits)",
                )
            )
        elif k <= first_bio:
            values = ", ".join(
                _fmt(profile.positions[i].agreement) for i in range(tech_len, k)
            )
            alternatives.append(
                Alternative(
                    k,
                    "weakly_compatible",
                    f"positions {tech_len + 1}-{k} have intermediate agreement "
                    f"({values}); consistent with non-templated bases in some "
                    "reads, which v1 reports but does not remove",
                )
            )
        else:
            stats = profile.positions[first_bio]
            alternatives.append(
                Alternative(
                    k,
                    "contradicted",
                    f"position {first_bio + 1} would be removed, but its agreement "
                    f"{_fmt(stats.agreement)} is at biology level "
                    f"(>= {config.agree_bio:.2f})",
                )
            )
    return tuple(alternatives)


def _explain(
    call: JunctionCall,
    profile: FrameProfile,
    first_bio: int,
    bio_level: float,
    groups: int,
    inside_bio: Sequence[int],
    where: str,
) -> str:
    """One readable sentence built only from the evidence (spec §7.4)."""
    parts = []
    if inside_bio:
        parts.append(
            f"Ambiguous: position(s) {', '.join(map(str, inside_bio))} inside the "
            "technical block agree at biology level"
        )
    if call.technical_length:
        described = " + ".join(f"{length} nt {kind}" for kind, length in call.blocks)
        values = [
            stats.agreement
            for stats in profile.positions[: call.technical_length]
            if stats.agreement is not None
        ]
        span = f"; within-fragment agreement {min(values):.2f}-{max(values):.2f}"
        parts.append(
            f"{call.technical_length} technical bases {where} ({described})"
            + (span if values else "")
        )
    else:
        parts.append(f"No technical bases {where}")
    if call.nta_length:
        zone = profile.positions[call.technical_length : first_bio]
        rate = max((r for r in call.nta_rates if r is not None), default=0.0)
        top = max(call.nta_bases.items(), key=lambda item: item[1], default=None)
        parts.append(
            f"then {call.nta_length} junction base(s) with intermediate agreement "
            f"({', '.join(_fmt(stats.agreement) for stats in zone)}; biology "
            f"{bio_level:.2f}), consistent with a non-templated base in about "
            f"{rate:.0%} of reads"
            + (f" (disagreeing bases {top[1]:.0%} {top[0]})" if top else "")
            + "; reported, not removed (v1)"
        )
    parts.append(
        f"agreement {_fmt(profile.positions[first_bio].agreement)} at position "
        f"{first_bio + 1} ({groups} groups)"
    )
    return "; ".join(parts) + "."
