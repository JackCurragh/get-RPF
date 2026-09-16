"""Step 2: locate the 3' anchor in each read (spec §5.2, Q1).

An adapter is accepted wherever it lies in the read, provided the read agrees
with it over the whole overlap: the full adapter when the read extends past
it, or the adapter's prefix up to the read end. A short internal resemblance
therefore cannot cut an insert, while an adapter in the middle of a long read
is still found. Platform sequence that follows the adapter (P7) is never an
anchor.

A homopolymer run directly before the adapter in most reads is a tail
(tailing-based libraries). The insert-side anchor is then the start of that
run, under the spec §6 convention.
"""

from __future__ import annotations

import statistics
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, Sequence, Tuple

from .config import InferenceConfig
from .model import Alternative, Answer, Evidence, Status
from .observe import Observation


@dataclass(frozen=True)
class Adapter:
    name: str
    sequence: str
    role: Literal["insert_3p", "downstream"] = "insert_3p"


CATALOGUE: Tuple[Adapter, ...] = (
    Adapter("truseq_3p", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
    Adapter("truseq_small_rna_3p", "TGGAATTCTCGGGTGCCAAGG"),
    Adapter("ingolia_linker_3p", "CTGTAGGCACCATCAAT"),
    Adapter("nextera_3p", "CTGTCTCTTATACACATCT"),
)
"""Insert-proximal 3' adapters: the first constant sequence after the insert."""

DOWNSTREAM: Tuple[Adapter, ...] = (
    Adapter("illumina_p7", "ATCTCGTATGCCGTCTTCTGCTTG", "downstream"),
)
"""Platform sequence that follows the adapter and index. Never an anchor."""


@dataclass(frozen=True)
class Hit:
    start: int
    overlap: int
    mismatches: int
    full: bool


@dataclass(frozen=True)
class ShortAdapterProbe:
    """Audit-only evidence for an exact adapter prefix at the read end."""

    adapter: Adapter
    overlap: int
    support: float
    mismatches: int
    competing: Tuple[Tuple[str, int, float], ...] = ()


def probe_short_end(
    reads: Sequence[str],
    config: Optional[InferenceConfig] = None,
    catalogue: Optional[Sequence[Adapter]] = None,
) -> Tuple[ShortAdapterProbe, ...]:
    """Record 5--9 nt end-anchored adapter candidates without resolving Q1.

    The normal ``locate`` path intentionally keeps its 10 nt threshold.  This
    probe is used by recovery reports only and requires an exact prefix at the
    3' read end, so an internal biological motif cannot become a boundary.
    """
    config = config or InferenceConfig()
    catalogue = catalogue or CATALOGUE
    if not reads:
        return ()
    candidates = []
    for adapter in catalogue:
        best = None
        for overlap in range(config.short_adapter_max_overlap, config.short_adapter_min_overlap - 1, -1):
            if sum(read.endswith(adapter.sequence[:overlap]) for read in reads) / len(reads) >= config.short_adapter_min_support:
                support = sum(read.endswith(adapter.sequence[:overlap]) for read in reads) / len(reads)
                best = ShortAdapterProbe(adapter, overlap, support, 0)
                break
        if best is not None:
            candidates.append(best)
    return tuple(sorted(candidates, key=lambda item: (-item.support, -item.overlap, item.adapter.name)))


@dataclass(frozen=True)
class AnchorCall:
    """The value of Q1."""

    adapter: Adapter
    source: Literal["catalogue", "de_novo"]
    support: float
    """Fraction of sampled reads with the adapter located after ``min_insert``."""
    dimer_fraction: float
    full_fraction: float
    """Of located reads, the fraction containing the whole adapter."""
    start_mode: int
    start_range: Tuple[int, int]
    """5th-95th percentile adapter start among located reads."""
    tail_base: Optional[str]
    tail_fraction: float
    downstream_consistency: Optional[float]
    """Of reads long enough to show it, the fraction where P7 follows the
    adapter within an index length. None when too few reads extend that far."""
    adapter_starts: Tuple[Optional[int], ...]
    insert_ends: Tuple[Optional[int], ...]
    """Per read: the insert-side anchor (tail start when there is a tail,
    otherwise the adapter start); None for unlocated reads and dimers."""


def locate(read: str, adapter: str, config: InferenceConfig) -> Optional[Hit]:
    """The leftmost acceptable match of ``adapter`` in ``read``."""
    seed = config.anchor_seed
    n, length = len(read), len(adapter)
    starts = set()
    for offset in (0, seed):
        probe = adapter[offset : offset + seed]
        if len(probe) < seed:
            continue
        i = read.find(probe)
        while i != -1:
            if i >= offset:
                starts.add(i - offset)
            i = read.find(probe, i + 1)
    for start in sorted(starts):
        overlap = min(length, n - start)
        if overlap < config.min_partial_overlap:
            continue
        mismatches = sum(
            1
            for base, expected in zip(read[start : start + overlap], adapter)
            if base != expected
        )
        if mismatches * 10 <= config.max_mismatch_per_10nt * overlap:
            return Hit(start, overlap, mismatches, overlap == length)
    return None


def find_anchor(
    reads: Sequence[str],
    observation: Observation,
    config: Optional[InferenceConfig] = None,
    catalogue: Sequence[Adapter] = CATALOGUE,
) -> Answer:
    """Answer Q1: which adapter anchors the reads, and where, per read."""
    config = config or InferenceConfig()
    total = len(reads)
    scans: List[Tuple[Adapter, Literal["catalogue", "de_novo"], List[Optional[Hit]]]]
    scans = [(a, "catalogue", _scan(reads, a.sequence, config)) for a in catalogue]
    if max((_support(hits, config) for _, _, hits in scans), default=0.0) < (
        config.min_anchor_support
    ):
        found = discover(reads, config)
        if found is not None:
            scans.append((found, "de_novo", _scan(reads, found.sequence, config)))
    scans.sort(key=lambda scan: _support(scan[2], config), reverse=True)

    downstream_alternatives = tuple(
        Alternative(
            element.name,
            "contradicted",
            f"{element.name} ({element.sequence[:12]}...) is platform sequence "
            "after the adapter and index; never an anchor (located in "
            f"{_support(_scan(reads, element.sequence, config), config):.1%} "
            "of reads)",
        )
        for element in DOWNSTREAM
    )
    observation_evidence = Evidence(
        "profile",
        "input_state",
        observation.input_state,
        observation.reads,
        note=f"modal length {observation.modal_length} nt "
        f"({observation.modal_fraction:.0%} of reads)",
    )
    short_candidates = probe_short_end(reads, config)

    best_adapter, best_source, best_hits = scans[0] if scans else (None, None, [])
    best_support = _support(best_hits, config)
    other_alternatives = tuple(
        Alternative(
            adapter.name,
            "contradicted",
            f"located in only {_support(hits, config):.1%} of reads "
            f"(need {config.min_anchor_support:.0%})",
        )
        for adapter, _, hits in scans[1:]
        if _support(hits, config) < config.min_anchor_support
    )

    if best_adapter is None or best_support < config.min_anchor_support:
        reads_too_short = (
            observation.modal_length
            < config.typical_footprint + config.min_partial_overlap
        )
        if observation.input_state == "trimmed":
            status = Status.NOT_OBSERVABLE
            reason = "the reads are already trimmed, so no adapter remains"
        elif reads_too_short:
            status = Status.NOT_OBSERVABLE
            reason = (
                f"reads ({observation.modal_length} nt) end before the adapter "
                "would start in most molecules"
            )
        else:
            status = Status.UNDERPOWERED
            reason = (
                "no known or de novo adapter reaches "
                f"{config.min_anchor_support:.0%} support"
            )
        support_text = (
            f"best candidate {best_adapter.name} in {best_support:.1%} of reads"
            if best_adapter is not None
            else "no candidate"
        )
        return Answer(
            "Q1",
            None,
            status,
            (
                observation_evidence,
                Evidence("anchor", "best_support", round(best_support, 4), total),
            ),
            other_alternatives
            + tuple(
                Alternative(
                    f"{item.adapter.name}:{item.overlap}nt",
                    "untested",
                    f"short exact 3' prefix in {item.support:.1%} of reads; audit evidence only",
                )
                for item in short_candidates
            )
            + downstream_alternatives,
            f"No 3' anchor: {reason} ({support_text}).",
        )

    competing = [
        (adapter, hits)
        for adapter, _, hits in scans[1:]
        if _support(hits, config) >= config.min_anchor_support
        and not _related(adapter.sequence, best_adapter.sequence)
    ]
    assert best_source is not None
    call = _anchor_call(reads, best_adapter, best_source, best_hits, config)
    evidence = [observation_evidence, *_anchor_evidence(call, best_hits, config)]
    evidence.extend(
        Evidence(
            "anchor_audit",
            "short_end_adapter_candidate",
            (item.adapter.name, item.overlap, round(item.support, 4)),
            len(reads),
            note="candidate only; not used to resolve ordinary Q1",
        )
        for item in short_candidates
    )
    alternatives = list(other_alternatives + downstream_alternatives)
    for adapter, hits in competing:
        alternatives.append(
            Alternative(
                adapter.name,
                "weakly_compatible",
                f"also located in {_support(hits, config):.1%} of reads: "
                "possibly a mixed library",
            )
        )
    if competing:
        status = Status.AMBIGUOUS
    elif call.tail_base is not None:
        status = Status.INTERVAL
    else:
        status = Status.RESOLVED
    return Answer(
        "Q1", call, status, tuple(evidence), tuple(alternatives), _explain(call)
    )


def search_window(
    call: Optional[AnchorCall], observation: Observation, config: InferenceConfig
) -> int:
    """Spec §2 rule 3: how far to look for technical sequence.

    The insert-side anchor mode minus a typical footprint estimates the total
    technical length. This sizes the window only; it never places a boundary.
    """
    if call is None:
        reach = observation.modal_length
    else:
        ends = [end for end in call.insert_ends if end is not None]
        reach = Counter(ends).most_common(1)[0][0] if ends else observation.modal_length
    technical = max(0, reach - config.typical_footprint)
    return min(config.max_window, max(config.window, technical + 4))


def discover(reads: Sequence[str], config: InferenceConfig) -> Optional[Adapter]:
    """De novo adapter discovery.

    An adapter is a constant element at a variable position, preceded by
    diverse sequence (many different inserts) and followed by consistent
    sequence. A biological fragment fails the first test: its upstream bases
    are the same molecule, so extension runs back to the read start.
    """
    sample = list(reads[: config.denovo_sample])
    k = config.denovo_k
    counts: Counter[str] = Counter()
    for read in sample:
        counts.update({read[p : p + k] for p in range(len(read) - k + 1)})
    minimum = config.denovo_min_support * len(sample)
    candidates = [
        kmer
        for kmer, count in counts.most_common(50)
        if count >= minimum
        and not _low_complexity(kmer)
        and not any(_related(kmer, element.sequence) for element in DOWNSTREAM)
    ]
    for kmer in candidates:
        occurrences = [(read, read.find(kmer)) for read in sample if kmer in read]
        positions = [p for _, p in occurrences]
        if len(positions) < 2 or statistics.pstdev(positions) < (
            config.denovo_min_position_sd
        ):
            continue
        element = _extend(kmer, occurrences, config)
        if element is not None:
            return Adapter("de_novo", element)
    return None


def _extend(
    kmer: str, occurrences: Sequence[Tuple[str, int]], config: InferenceConfig
) -> Optional[str]:
    sequence = kmer
    left = 0
    while len(sequence) < config.denovo_max_length:
        bases = Counter(
            read[p - left - 1] for read, p in occurrences if p - left - 1 >= 0
        )
        covered = sum(bases.values())
        if covered < 0.5 * len(occurrences):
            return None  # runs back to the read start: a 5' element or biology
        base, count = bases.most_common(1)[0]
        if count / covered < config.denovo_extend_dominance:
            break  # diverse upstream: this is where the adapter starts
        if sequence[:2] == base * 2:
            break  # a growing homopolymer is a tail, not the adapter
        sequence = base + sequence
        left += 1
    right = 0
    while len(sequence) < config.denovo_max_length:
        bases = Counter(
            read[p + len(kmer) + right]
            for read, p in occurrences
            if p + len(kmer) + right < len(read)
        )
        covered = sum(bases.values())
        if covered < 0.5 * len(occurrences):
            break
        base, count = bases.most_common(1)[0]
        if count / covered < config.denovo_extend_dominance:
            break
        sequence = sequence + base
        right += 1
    return sequence


def _scan(
    reads: Sequence[str], adapter: str, config: InferenceConfig
) -> List[Optional[Hit]]:
    return [locate(read, adapter, config) for read in reads]


def _support(hits: Sequence[Optional[Hit]], config: InferenceConfig) -> float:
    if not hits:
        return 0.0
    located = sum(
        1 for hit in hits if hit is not None and hit.start >= config.min_insert
    )
    return located / len(hits)


def _related(a: str, b: str, k: int = 10) -> bool:
    """Do two sequences share a k-mer (nested catalogue variants, or P7)?"""
    kmers = {a[i : i + k] for i in range(len(a) - k + 1)}
    return any(b[i : i + k] in kmers for i in range(len(b) - k + 1))


def _low_complexity(kmer: str) -> bool:
    if max(Counter(kmer).values()) / len(kmer) > 0.5:
        return True
    run = 1
    for previous, current in zip(kmer, kmer[1:]):
        run = run + 1 if current == previous else 1
        if run >= 5:
            return True
    return False


def _run_before(read: str, end: int, base: str) -> int:
    """Length of the ``base`` run ending at ``end``, read backwards.

    An isolated mismatch is stepped over when the run resumes for at least
    three bases beyond it (a sequencing error inside a tail). The first one is
    always allowed; further ones at most one per ten bases. The run always
    starts on a ``base``.
    """
    run = 0
    mismatches = 0
    i = end - 1
    while i >= 0:
        if read[i] == base:
            run = end - i
        elif (
            i >= 3
            and read[i - 3 : i] == base * 3
            and (mismatches == 0 or (mismatches + 1) * 10 <= end - i)
        ):
            mismatches += 1
        else:
            break
        i -= 1
    return run


def _anchor_call(
    reads: Sequence[str],
    adapter: Adapter,
    source: Literal["catalogue", "de_novo"],
    hits: Sequence[Optional[Hit]],
    config: InferenceConfig,
) -> AnchorCall:
    total = len(reads)
    starts: List[Optional[int]] = [
        hit.start if hit is not None and hit.start >= config.min_insert else None
        for hit in hits
    ]
    located = [s for s in starts if s is not None]
    dimers = sum(1 for hit in hits if hit is not None and hit.start < config.min_insert)
    full = sum(
        1
        for hit in hits
        if hit is not None and hit.start >= config.min_insert and hit.full
    )
    ordered = sorted(located)

    tail_base: Optional[str] = None
    tail_fraction = 0.0
    run_fractions: Dict[str, float] = {}
    for base in "ACGT":
        long_runs = sum(
            1
            for read, start in zip(reads, starts)
            if start is not None
            and _run_before(read, start, base) >= config.tail_min_run
        )
        run_fractions[base] = long_runs / len(located) if located else 0.0
    best_base = max(run_fractions, key=lambda base: run_fractions[base])
    if run_fractions[best_base] >= config.tail_min_fraction:
        tail_base, tail_fraction = best_base, run_fractions[best_base]
    insert_ends: List[Optional[int]] = [
        (
            None
            if start is None
            else start - (_run_before(read, start, tail_base) if tail_base else 0)
        )
        for read, start in zip(reads, starts)
    ]

    return AnchorCall(
        adapter=adapter,
        source=source,
        support=len(located) / total,
        dimer_fraction=dimers / total,
        full_fraction=full / len(located) if located else 0.0,
        start_mode=Counter(located).most_common(1)[0][0],
        start_range=(
            ordered[int(0.05 * (len(ordered) - 1))],
            ordered[int(0.95 * (len(ordered) - 1))],
        ),
        tail_base=tail_base,
        tail_fraction=tail_fraction,
        downstream_consistency=_downstream_consistency(reads, hits, adapter, config),
        adapter_starts=tuple(starts),
        insert_ends=tuple(insert_ends),
    )


def _downstream_consistency(
    reads: Sequence[str],
    hits: Sequence[Optional[Hit]],
    adapter: Adapter,
    config: InferenceConfig,
) -> Optional[float]:
    probe = DOWNSTREAM[0].sequence[: config.anchor_seed]
    reach = config.max_index_length + config.anchor_seed
    eligible = consistent = 0
    for read, hit in zip(reads, hits):
        if hit is None or not hit.full or hit.start < config.min_insert:
            continue
        after = read[hit.start + len(adapter.sequence) :]
        if len(after) < reach:
            continue
        eligible += 1
        if probe in after[:reach]:
            consistent += 1
    if eligible < config.min_scored_per_position:
        return None
    return consistent / eligible


def _anchor_evidence(
    call: AnchorCall, hits: Sequence[Optional[Hit]], config: InferenceConfig
) -> List[Evidence]:
    located = [h for h in hits if h is not None and h.start >= config.min_insert]
    mismatch_rate = (
        sum(h.mismatches for h in located) / sum(h.overlap for h in located)
        if located
        else 0.0
    )
    evidence = [
        Evidence(
            "anchor",
            "support",
            round(call.support, 4),
            len(hits),
            note=f"{call.adapter.name} ({call.source}): {call.adapter.sequence}",
        ),
        Evidence(
            "anchor",
            "start_distribution",
            (call.start_mode, call.start_range),
            len(located),
            note="mode and 5-95% range of the adapter start",
        ),
        Evidence(
            "anchor",
            "full_adapter_fraction",
            round(call.full_fraction, 4),
            len(located),
        ),
        Evidence("anchor", "dimer_fraction", round(call.dimer_fraction, 4), len(hits)),
        Evidence(
            "anchor",
            "mismatch_rate",
            round(mismatch_rate, 4),
            len(located),
            note="mismatches per adapter base over the matched overlap",
        ),
    ]
    if call.downstream_consistency is not None:
        evidence.append(
            Evidence(
                "anchor",
                "downstream_p7_consistency",
                round(call.downstream_consistency, 4),
                len(located),
                note="reads where P7 follows the adapter within an index length",
            )
        )
    if call.tail_base is not None:
        runs = [
            start - end
            for start, end in zip(call.adapter_starts, call.insert_ends)
            if start is not None and end is not None
        ]
        evidence.append(
            Evidence(
                "anchor",
                "tail",
                (call.tail_base, round(call.tail_fraction, 4), statistics.median(runs)),
                len(runs),
                note="tail base, fraction of anchored reads with a run of "
                f">= {config.tail_min_run}, median run length; insert end is "
                "placed at the first base of the run (spec §6)",
            )
        )
    return evidence


def _explain(call: AnchorCall) -> str:
    parts = [
        f"{call.adapter.name} ({call.source}) located in {call.support:.0%} of "
        f"reads, starting at position {call.start_mode + 1} "
        f"(5-95%: {call.start_range[0] + 1}-{call.start_range[1] + 1})",
        f"{call.full_fraction:.0%} of those contain the whole adapter",
    ]
    if call.downstream_consistency is not None:
        parts.append(
            f"P7 follows it in {call.downstream_consistency:.0%} of reads long "
            "enough to show it"
        )
    if call.dimer_fraction:
        parts.append(f"{call.dimer_fraction:.1%} of reads are adapter dimers")
    if call.tail_base is not None:
        parts.append(
            f"a poly({call.tail_base}) tail precedes the adapter in "
            f"{call.tail_fraction:.0%} of anchored reads; the insert end is "
            "placed at the first base of the run (interval by convention)"
        )
    return "; ".join(parts) + "."
