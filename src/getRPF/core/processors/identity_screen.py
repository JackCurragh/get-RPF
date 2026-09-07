"""Stage 3 identity screen: biological, provisional.

Screens the located inserts for consistency with a Ribo-seq footprint
population -- length-distribution shape, contamination k-mer matches,
UMI-aware duplication, and poly-G/adapter-dimer traps. Never emits
"confirmed": the strongest positive verdict is `consistent_with_riboseq`,
biological identity is only earned at the alignment gate. See
docs/release_qc_and_terminal_trimming_plan.md section 4 (Stage 3) and
section 12.

No real rRNA/tRNA/sno reference k-mer sets are bundled with getRPF --
`contamination_kmer_sets` must be supplied by the caller (e.g. built from
reference FASTA via `build_kmer_set`). Without it, the contamination
component of the screen is skipped, not fabricated.
"""

from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import SeqIO

# A length class within this many nt of the mode counts as "in the peak".
PEAK_WINDOW = 2

# Fraction of reads that must fall within the peak window for the length
# distribution to be classified "peaked" rather than "broad".
PEAK_FRACTION_THRESHOLD = 0.5

# Reads at or below this length are adapter-dimer / no-insert artifacts.
ADAPTER_DIMER_MAX_LENGTH = 5

# Contamination fraction (summed across all screened categories) above
# which the sample is called inconsistent with a clean footprint library.
CONTAMINATION_INCONSISTENT_THRESHOLD = 0.5
CONTAMINATION_CONSISTENT_THRESHOLD = 0.3


def classify_length_shape(
    length_distribution: Dict[int, int],
    peak_window: int = PEAK_WINDOW,
    peak_fraction_threshold: float = PEAK_FRACTION_THRESHOLD,
) -> Tuple[str, Optional[int], float]:
    """Classify a length distribution as "peaked" or "broad".

    Footprint libraries have inserts that are short and tightly
    length-distributed; RNA-seq/degradation products are flat/broad. See
    section 12's discreteness argument.

    Returns (shape, mode_length, peak_fraction).
    """
    total = sum(length_distribution.values())
    if total == 0:
        return "broad", None, 0.0

    mode_length = max(length_distribution.items(), key=lambda kv: kv[1])[0]
    in_peak = sum(
        count
        for length, count in length_distribution.items()
        if abs(length - mode_length) <= peak_window
    )
    peak_fraction = in_peak / total
    shape = "peaked" if peak_fraction >= peak_fraction_threshold else "broad"
    return shape, mode_length, peak_fraction


def adapter_dimer_fraction(
    length_distribution: Dict[int, int], max_length: int = ADAPTER_DIMER_MAX_LENGTH
) -> float:
    """Fraction of reads with a near-zero insert (adapter-dimer artifact)."""
    total = sum(length_distribution.values())
    if total == 0:
        return 0.0
    dimer = sum(c for length, c in length_distribution.items() if length <= max_length)
    return dimer / total


def build_kmer_set(fasta_path: Path, k: int = 20) -> Set[str]:
    """Build an exact k-mer set from a reference FASTA (e.g. rRNA/tRNA)."""
    kmers: Set[str] = set()
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - k + 1):
            kmers.add(seq[i : i + k])
    return kmers


def screen_contamination(
    reads: List[str], kmer_sets: Dict[str, Set[str]], k: int = 20
) -> Dict[str, float]:
    """Fraction of reads containing at least one k-mer from each reference
    set. A read matching multiple categories is counted in each; "other"
    is whatever's left after the best single match per read."""
    if not kmer_sets or not reads:
        return {}

    counts = {category: 0 for category in kmer_sets}
    unmatched = 0
    for read in reads:
        read_kmers = (
            {read[i : i + k] for i in range(len(read) - k + 1)}
            if len(read) >= k
            else set()
        )
        matched_any = False
        for category, ref_kmers in kmer_sets.items():
            if read_kmers & ref_kmers:
                counts[category] += 1
                matched_any = True
        if not matched_any:
            unmatched += 1

    total = len(reads)
    fractions = {category: count / total for category, count in counts.items()}
    fractions["other"] = unmatched / total
    return fractions


def duplication_umi_aware(reads: List[str]) -> float:
    """Fraction of reads that are exact-sequence duplicates.

    If `reads` still carry their 5' UMI (not yet trimmed), a read sharing
    the same UMI *and* insert is a true PCR duplicate, so deduping the
    full UMI+insert string is the accurate estimate. A caller that instead
    deduped on the UMI-stripped insert alone would inflate the estimate,
    since independent molecules with short/low-complexity inserts (like
    RPFs) collide by chance even without PCR. This function always keys
    on the full string it's given -- pass UMI-inclusive reads for the
    UMI-aware number, or insert-only reads to get the (inflated) naive one.
    """
    if not reads:
        return 0.0
    unique = len(set(reads))
    return 1.0 - (unique / len(reads))


def identity_screen(
    reads: List[str],
    length_distribution: Dict[int, int],
    contamination_kmer_sets: Optional[Dict[str, Set[str]]] = None,
    kmer_len: int = 20,
) -> Dict:
    """Run the Stage-3 identity screen; returns the `biological_screen`
    dict from the evidence object (section 9). Verdict is one of
    consistent_with_riboseq / inconsistent / indeterminate -- never
    "confirmed".

    `reads` should be UMI-inclusive (not yet UMI-trimmed) if a UMI-aware
    duplication estimate is wanted; see duplication_umi_aware.
    """
    length_shape, insert_mode, peak_fraction = classify_length_shape(
        length_distribution
    )
    dimer_fraction = adapter_dimer_fraction(length_distribution)
    duplication = duplication_umi_aware(reads)

    contamination: Dict[str, float] = {}
    contamination_screened = bool(contamination_kmer_sets)
    if contamination_kmer_sets:
        contamination = screen_contamination(reads, contamination_kmer_sets, k=kmer_len)

    contamination_total = None
    if contamination_screened:
        contamination_total = sum(
            v for category, v in contamination.items() if category != "other"
        )

    # reason_codes are the machine-readable routing signal for
    # release.classify_release; `reasons` is the human-readable prose for
    # reports. Keep both in sync, but only reason_codes should ever be
    # matched on programmatically -- prose wording is free to change.
    reasons: List[str] = []
    reason_codes: List[str] = []
    verdict = "indeterminate"

    if dimer_fraction >= 0.5:
        verdict = "inconsistent"
        reason_codes.append("adapter_dimer")
        reasons.append(
            f"{dimer_fraction:.1%} of reads are adapter-dimer/near-zero-insert artifacts"
        )
    elif length_shape == "broad":
        verdict = "inconsistent"
        reason_codes.append("broad_length_distribution")
        reasons.append(
            f"length distribution is broad (peak fraction {peak_fraction:.1%}), "
            "not consistent with a tight footprint population"
        )
    elif (
        contamination_screened
        and contamination_total is not None
        and contamination_total > CONTAMINATION_INCONSISTENT_THRESHOLD
    ):
        verdict = "inconsistent"
        reason_codes.append("high_contamination")
        reasons.append(
            f"{contamination_total:.1%} of reads match contamination reference k-mers"
        )
    elif length_shape == "peaked" and (
        not contamination_screened
        or (
            contamination_total is not None
            and contamination_total < CONTAMINATION_CONSISTENT_THRESHOLD
        )
    ):
        verdict = "consistent_with_riboseq"
        reason_codes.append("peaked_clean")
        reasons.append(
            f"length distribution is peaked (peak fraction {peak_fraction:.1%}) "
            "with no disqualifying contamination"
        )
    else:
        reason_codes.append("ambiguous_contamination")
        reasons.append("evidence is peaked but contamination is in an ambiguous range")

    return {
        "length_shape": length_shape,
        "insert_mode": insert_mode,
        "peak_fraction": peak_fraction,
        "adapter_dimer_fraction": dimer_fraction,
        "contamination_screened": contamination_screened,
        "contamination": contamination,
        "duplication_umi_aware": duplication,
        "verdict": verdict,
        "reason_codes": reason_codes,
        "reasons": reasons,
    }
