"""Boundary estimators and the concordance decision function.

Implements docs/release_qc_and_terminal_trimming_plan.md section 6: up to
four independent per-length-class estimators of how many bases to trim from
an end, reconciled via a concordance/support-weighted confidence score into
a three-way decision (safe_to_trim / safe_to_infer_but_not_trim / hold).

All estimators and decisions are expressed in a single unit: trim_bases,
the number of bases to remove from the given end of that length class. This
keeps spread/concordance directly comparable across estimators without a
separate position<->trim_bases conversion.
"""

import statistics
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .signals import SignalStats

# Below this per-length-class read count, no rule may fire regardless of
# how concordant the estimators are -- the evidence is statistical noise.
DEFAULT_SUPPORT_FLOOR = 200

# Support at/above which support_factor saturates to 1.0.
DEFAULT_S_MIN = 10_000.0

# Confidence thresholds. 5' trims are riskier (they corrupt frame/P-site
# interpretation downstream) so the bar is strictly higher than 3'.
DEFAULT_TAU_HIGH_3P = 0.7
DEFAULT_TAU_HIGH_5P = 0.85

DEFAULT_KMER_MIN_FREQ = 0.5
DEFAULT_ENTROPY_FRAC = 0.5
DEFAULT_MIN_ADAPTER_OVERLAP = 4


@dataclass
class LengthClassDecision:
    """The concordance decision for one (end, read-length) class."""

    end: str
    length: int
    estimates: Dict[str, int]
    support: int
    spread: Optional[int]
    n_indep: int
    concordance: float
    confidence: float
    decision: str  # "safe_to_trim" | "safe_to_infer_but_not_trim" | "hold"
    trim_bases: Optional[int] = None
    dominant_base: Optional[str] = None
    dominant_fraction: Optional[float] = None

    def to_dict(self) -> Dict:
        return {
            "end": self.end,
            "length": self.length,
            "estimates": self.estimates,
            "support": self.support,
            "spread": self.spread,
            "n_indep": self.n_indep,
            "concordance": self.concordance,
            "confidence": self.confidence,
            "decision": self.decision,
            "trim_bases": self.trim_bases,
            "dominant_base": self.dominant_base,
            "dominant_fraction": self.dominant_fraction,
        }


def concordance_and_confidence(
    estimates: Dict[str, int],
    support: int,
    s_min: float = DEFAULT_S_MIN,
    support_floor: int = DEFAULT_SUPPORT_FLOOR,
) -> Tuple[Optional[int], int, float, float]:
    """Section 6's concordance/confidence formula.

    Returns (spread, n_indep, concordance, confidence).
    """
    n_indep = len(estimates)
    if n_indep == 0 or support < support_floor:
        return None, n_indep, 0.0, 0.0

    values = list(estimates.values())
    spread = max(values) - min(values)

    if spread == 0:
        concordance = 1.0
    elif spread == 1:
        concordance = 0.7
    elif spread == 2:
        concordance = 0.3
    else:
        concordance = 0.0

    support_factor = min(1.0, support / s_min)
    confidence = concordance * support_factor
    return spread, n_indep, concordance, confidence


def _consensus_trim(estimates: Dict[str, int]) -> int:
    """Pick the representative trim_bases when estimators agree closely:
    the mode, tie-broken by the smallest value (the more conservative
    edit)."""
    counts = Counter(estimates.values())
    max_count = max(counts.values())
    candidates = [v for v, c in counts.items() if c == max_count]
    return min(candidates)


def _leading_run_length(values: List[float], continues: "callable") -> int:
    """Count leading values for which `continues(value)` holds, stopping at
    the first value that breaks it."""
    onset = 0
    for value in values:
        if not continues(value):
            break
        onset += 1
    return onset


def _b_kmer(compositions: List[Dict[str, float]], min_freq: float) -> Optional[int]:
    """De novo terminal k-mer/consensus-base enrichment: how many leading
    positions have a single base dominating at >= min_freq."""
    if not compositions:
        return None
    freqs = [max(c.values()) if c else 0.0 for c in compositions]
    onset = _leading_run_length(freqs, lambda f: f >= min_freq)
    return onset if onset > 0 else None


def _b_entropy(entropies: List[float], frac: float, baseline_offset: int = 10) -> Optional[int]:
    """Entropy-cliff onset: how many leading positions have entropy below
    `frac` of the median "insert region" entropy (positions beyond
    `baseline_offset`, used as a proxy for biological baseline entropy)."""
    if not entropies:
        return None
    baseline_region = entropies[baseline_offset:] if len(entropies) > baseline_offset else entropies
    if not baseline_region:
        return None
    median_baseline = statistics.median(baseline_region)
    if median_baseline <= 0:
        return None
    threshold = frac * median_baseline
    onset = _leading_run_length(entropies, lambda e: e < threshold)
    return onset if onset > 0 else None


def _b_adapter(
    reads: Optional[List[str]],
    known_elements: Optional[List[str]],
    anchor: str,
    min_overlap: int = DEFAULT_MIN_ADAPTER_OVERLAP,
    min_hit_fraction: float = 0.5,
    sample_cap: int = 500,
) -> Optional[int]:
    """Exact/fuzzy match against known adapter/element sequences.

    anchor="3p": match element prefixes against read suffixes (3' adapter).
    anchor="5p": match element suffixes against read prefixes (5' element,
    e.g. a UMI/barcode/linker that precedes the insert).
    """
    if not reads or not known_elements:
        return None

    sample = reads[:sample_cap]
    overlap_hits: Counter = Counter()
    for read in sample:
        best_overlap = 0
        for element in known_elements:
            max_overlap = min(len(element), len(read))
            for overlap in range(max_overlap, min_overlap - 1, -1):
                if anchor == "3p":
                    matched = read[-overlap:] == element[:overlap]
                else:
                    matched = read[:overlap] == element[-overlap:]
                if matched:
                    best_overlap = max(best_overlap, overlap)
                    break
        if best_overlap:
            overlap_hits[best_overlap] += 1

    if not overlap_hits:
        return None

    modal_overlap, hit_count = overlap_hits.most_common(1)[0]
    if hit_count / len(sample) < min_hit_fraction:
        return None
    return modal_overlap


def _b_lenmode(length: int, other_trim_amounts: Dict[int, int], min_agreeing: int = 2) -> Optional[int]:
    """Cross-length consistency: if several *other* length classes agree on
    a trim amount, that amount is a candidate for this class too (a
    footprint library's adapter onset is a fixed offset from read length,
    not a fixed absolute position)."""
    if not other_trim_amounts:
        return None
    counts = Counter(other_trim_amounts.values())
    mode_trim, mode_count = counts.most_common(1)[0]
    if mode_count < min_agreeing:
        return None
    if mode_trim < 0 or mode_trim >= length:
        return None
    return mode_trim


def _decide(
    estimates: Dict[str, int],
    support: int,
    end: str,
    named_element: bool,
    tau_high_3p: float,
    tau_high_5p: float,
    s_min: float,
    support_floor: int,
) -> Tuple[str, Optional[int], int, float, float]:
    spread, n_indep, concordance, confidence = concordance_and_confidence(
        estimates, support, s_min, support_floor
    )

    if n_indep == 0 or support < support_floor:
        return "hold", spread, n_indep, concordance, confidence

    tau_high = tau_high_5p if end == "5p" else tau_high_3p
    strong = confidence >= tau_high and n_indep >= 2

    if strong:
        if end == "5p" and not named_element:
            # 5' trimming is never performed on de novo terminal bias
            # alone -- it needs a named element from a seqspec or a
            # high-posterior HMM segment.
            return "safe_to_infer_but_not_trim", spread, n_indep, concordance, confidence
        return "safe_to_trim", spread, n_indep, concordance, confidence

    return "safe_to_infer_but_not_trim", spread, n_indep, concordance, confidence


def estimate_boundaries(
    per_length_stats: Dict[int, SignalStats],
    per_length_support: Dict[int, int],
    end: str = "3p",
    per_length_reads: Optional[Dict[int, List[str]]] = None,
    known_elements: Optional[List[str]] = None,
    named_element_lengths: Optional[Set[int]] = None,
    s_min: float = DEFAULT_S_MIN,
    support_floor: int = DEFAULT_SUPPORT_FLOOR,
    tau_high_3p: float = DEFAULT_TAU_HIGH_3P,
    tau_high_5p: float = DEFAULT_TAU_HIGH_5P,
    kmer_min_freq: float = DEFAULT_KMER_MIN_FREQ,
    entropy_frac: float = DEFAULT_ENTROPY_FRAC,
) -> Dict[int, LengthClassDecision]:
    """Run the boundary estimators and concordance decision for every
    length class, for one end ("3p" or "5p"). See module docstring and
    docs/release_qc_and_terminal_trimming_plan.md section 6.
    """
    assert end in ("3p", "5p")
    named_element_lengths = named_element_lengths or set()

    raw_estimates: Dict[int, Dict[str, int]] = {}
    for length, stats in per_length_stats.items():
        compositions = stats.composition_3p if end == "3p" else stats.composition_5p
        entropies = stats.entropy_3p if end == "3p" else stats.entropy_5p

        estimates: Dict[str, int] = {}

        reads = per_length_reads.get(length) if per_length_reads else None
        adapter_trim = _b_adapter(reads, known_elements, anchor=end)
        if adapter_trim is not None:
            estimates["b_adapter"] = adapter_trim

        kmer_trim = _b_kmer(compositions, kmer_min_freq)
        if kmer_trim is not None:
            estimates["b_kmer"] = kmer_trim

        entropy_trim = _b_entropy(entropies, entropy_frac)
        if entropy_trim is not None:
            estimates["b_entropy"] = entropy_trim

        raw_estimates[length] = estimates

    # Cross-length consensus (b_lenmode) needs the other classes' strongly
    # concordant estimates first, so it's a second pass.
    provisional_trims: Dict[int, int] = {}
    for length, estimates in raw_estimates.items():
        values = list(estimates.values())
        if len(values) >= 2 and max(values) - min(values) == 0:
            provisional_trims[length] = values[0]

    decisions: Dict[int, LengthClassDecision] = {}
    for length, stats in per_length_stats.items():
        estimates = dict(raw_estimates[length])
        other_trims = {l: t for l, t in provisional_trims.items() if l != length}
        lenmode_trim = _b_lenmode(length, other_trims)
        if lenmode_trim is not None:
            estimates["b_lenmode"] = lenmode_trim

        support = per_length_support.get(length, 0)
        named_element = length in named_element_lengths

        decision, spread, n_indep, concordance, confidence = _decide(
            estimates, support, end, named_element,
            tau_high_3p, tau_high_5p, s_min, support_floor,
        )

        trim_bases = None
        dominant_base = None
        dominant_fraction = None
        if decision == "safe_to_trim":
            trim_bases = _consensus_trim(estimates)
            compositions = stats.composition_3p if end == "3p" else stats.composition_5p
            if trim_bases > 0 and compositions and compositions[0]:
                dominant_base, dominant_fraction = max(
                    compositions[0].items(), key=lambda kv: kv[1]
                )

        decisions[length] = LengthClassDecision(
            end=end,
            length=length,
            estimates=estimates,
            support=support,
            spread=spread,
            n_indep=n_indep,
            concordance=concordance,
            confidence=confidence,
            decision=decision,
            trim_bases=trim_bases,
            dominant_base=dominant_base,
            dominant_fraction=dominant_fraction,
        )

    return decisions
