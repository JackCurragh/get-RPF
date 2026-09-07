"""Signal processing for read architecture detection.

This module provides the Metric Engine that calculates per-base statistics
from a sample of reads to drive architecture detection.
"""

import math
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List


@dataclass
class SignalStats:
    """Container for per-base signal statistics."""

    # 5' aligned metrics (index 0 is read start)
    entropy_5p: List[float]
    composition_5p: List[Dict[str, float]]
    dinucleotide_5p: List[Dict[str, float]]

    # 3' aligned metrics (index 0 is read end)
    entropy_3p: List[float]
    composition_3p: List[Dict[str, float]]
    dinucleotide_3p: List[Dict[str, float]]

    sample_size: int


class SignalProcessor:
    """Calculates informative signals from read populations."""

    def process_reads(
        self, reads: List[str], compute_dinucleotide: bool = True
    ) -> SignalStats:
        """Calculate signal statistics for a set of reads.

        Args:
            reads: List of DNA sequences
            compute_dinucleotide: Whether to compute per-position dinucleotide
                frequencies. No current caller reads dinucleotide_5p/3p, so
                callers that don't need them (e.g. sketch/boundary
                estimation) should pass False to skip the O(reads * length)
                work. Defaults to True to preserve prior behavior.

        Returns:
            SignalStats object with calculated metrics
        """
        if not reads:
            return self._empty_stats()

        # Filter for reasonable length reads to avoid skewing alignment
        # We generally expect ribosome profiling reads to be > 20nt
        valid_reads = [r for r in reads if len(r) >= 20]
        if not valid_reads:
            return self._empty_stats()

        # Limit processing depth for speed (though we aim for high N)
        # 50k reads is robust enough for high confidence
        sample_reads = valid_reads[:50000]

        entropy_5p, composition_5p = self._calculate_metrics(sample_reads, "5p")
        entropy_3p, composition_3p = self._calculate_metrics(sample_reads, "3p")

        return SignalStats(
            entropy_5p=entropy_5p,
            composition_5p=composition_5p,
            dinucleotide_5p=(
                self._calculate_dinucleotides(sample_reads, align="5p")
                if compute_dinucleotide
                else []
            ),
            entropy_3p=entropy_3p,
            composition_3p=composition_3p,
            dinucleotide_3p=(
                self._calculate_dinucleotides(sample_reads, align="3p")
                if compute_dinucleotide
                else []
            ),
            sample_size=len(sample_reads),
        )

    def _calculate_metrics(
        self, reads: List[str], align: str
    ) -> tuple[List[float], List[Dict[str, float]]]:
        """Calculate entropy and composition in one positional pass.

        Entropy and composition use the same per-position base counts. Keeping
        them together avoids traversing every read twice for each end.
        """
        length = self._get_max_len(reads)
        entropies: List[float] = []
        compositions: List[Dict[str, float]] = []

        for position in range(length):
            counts: Counter = Counter()
            for read in reads:
                index = position if align == "5p" else -(position + 1)
                if (align == "5p" and index < len(read)) or (
                    align == "3p" and abs(index) <= len(read)
                ):
                    counts[read[index]] += 1

            total = sum(counts.values())
            if total == 0:
                entropies.append(0.0)
                compositions.append({})
                continue

            entropy = 0.0
            frequencies = {}
            for base, count in counts.items():
                frequency = count / total
                frequencies[base] = frequency
                entropy -= frequency * math.log2(frequency)

            entropies.append(entropy)
            compositions.append(frequencies)

        return entropies, compositions

    def _calculate_entropy(self, reads: List[str], align: str = "5p") -> List[float]:
        """Calculate per-position Shannon entropy."""
        length = self._get_max_len(reads)
        entropies = []

        for i in range(length):
            bases = self._get_bases_at_pos(reads, i, align)
            if not bases:
                entropies.append(0.0)
                continue

            counts = Counter(bases)
            total = len(bases)
            entropy = 0.0

            for count in counts.values():
                p = count / total
                entropy -= p * math.log2(p)

            entropies.append(entropy)

        return entropies

    def _calculate_composition(
        self, reads: List[str], align: str = "5p"
    ) -> List[Dict[str, float]]:
        """Calculate per-position nucleotide frequencies."""
        length = self._get_max_len(reads)
        compositions: List[Dict[str, float]] = []

        for i in range(length):
            bases = self._get_bases_at_pos(reads, i, align)
            if not bases:
                compositions.append({})
                continue

            counts = Counter(bases)
            total = len(bases)
            freqs = {base: count / total for base, count in counts.items()}
            compositions.append(freqs)

        return compositions

    def _calculate_dinucleotides(
        self, reads: List[str], align: str = "5p"
    ) -> List[Dict[str, float]]:
        """Calculate per-position dinucleotide frequencies."""
        length = self._get_max_len(reads)
        di_freqs: List[Dict[str, float]] = []

        for i in range(length - 1):
            # For dinucleotides, we need position i and i+1
            start_pos = i

            dinucs = []
            for read in reads:
                idx = start_pos if align == "5p" else (len(read) - 1 - start_pos)
                next_idx = idx + 1 if align == "5p" else idx - 1

                # Check bounds
                if 0 <= idx < len(read) and 0 <= next_idx < len(read):
                    if align == "5p":
                        d = read[idx : idx + 2]
                    else:
                        # For 3', we iterate backwards but read sequence forwards
                        # If index is 0 (last base), we want base at -2,-1
                        # idx is distance from end.
                        # Real indices:
                        # 5' aligned: 0,1 -> 1,2 ...
                        # 3' aligned: (len-1, len) -> (len-2, len-1)
                        # Actually simpler: just grab the slice relative to end
                        # idx is 0-based index from end. 0 is last base.
                        # so pos i means (len - 1 - i).
                        # We want dinuc starting at that pos?
                        # Let's define 3' dinuc at pos i as: base at -(i+2) and -(i+1)
                        # e.g. pos 0 is last 2 bases.
                        p1 = -(i + 2)
                        p2 = -(i)
                        # slice notation [low:high]
                        if abs(p1) <= len(read):
                            d = read[p1:p2] if p2 != 0 else read[p1:]
                        else:
                            continue

                    if len(d) == 2:
                        dinucs.append(d)

            if not dinucs:
                di_freqs.append({})
                continue

            counts = Counter(dinucs)
            total = len(dinucs)
            df = {d: count / total for d, count in counts.items()}
            di_freqs.append(df)

        return di_freqs

    def _get_bases_at_pos(self, reads: List[str], pos: int, align: str) -> List[str]:
        """Get all bases at a specific position across reads."""
        bases = []
        for read in reads:
            if align == "5p":
                if pos < len(read):
                    bases.append(read[pos])
            else:  # 3p
                # pos 0 is the last base (index -1)
                # pos 1 is second to last (index -2)
                idx = -(pos + 1)
                if abs(idx) <= len(read):
                    bases.append(read[idx])
        return bases

    def _get_max_len(self, reads: List[str]) -> int:
        """Get reasonable max length for analysis (95th percentile)."""
        if not reads:
            return 0
        lens = sorted([len(r) for r in reads])
        return lens[int(len(lens) * 0.95)]

    def _empty_stats(self) -> SignalStats:
        """Return empty statistics object."""
        return SignalStats([], [], [], [], [], [], 0)


def process_reads(reads: List[str], compute_dinucleotide: bool = True) -> SignalStats:
    """Calculate signal statistics for a read population.

    Signal calculation is stateless. This function is the preferred API for
    new code; ``SignalProcessor`` remains as a compatibility wrapper for
    callers that used the original class-based interface.
    """
    return SignalProcessor().process_reads(
        reads, compute_dinucleotide=compute_dinucleotide
    )


def empty_stats() -> SignalStats:
    """Return an empty signal result for an empty read population."""
    return SignalProcessor()._empty_stats()
