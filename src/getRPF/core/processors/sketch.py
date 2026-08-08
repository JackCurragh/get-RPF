"""Stage-0 sketch: the reference-free substrate for all downstream stages.

Computes, once, on a bounded subsample: length distribution, coverage-normalized
per-position entropy/composition (pooled and per read-length class), a
per-position quality profile, and top 3' terminal k-mers. See
docs/release_qc_and_terminal_trimming_plan.md sections 4-5.
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from ...utils.file_utils import get_file_opener
from .collapsed import parse_collapsed_fasta
from .signals import SignalStats, empty_stats, process_reads

# Base quality below this (mean Phred) at a low-entropy terminal position is
# characteristic of 2-colour-chemistry dark cycles (poly-G), not a real
# adapter/constant region.
POLYG_QUALITY_COLLAPSE_THRESHOLD = 10.0


@dataclass
class Sketch:
    """Stage-0 sketch: substrate for boundary estimation and identity screening."""

    length_distribution: Dict[int, int]
    sample_size: int

    # Pooled (all lengths together), 5'/3' anchored.
    pooled: SignalStats

    # Per read-length class, 5'/3' anchored. Only length classes with
    # sufficient reads to be meaningful are included (see filter in build()).
    per_length: Dict[int, SignalStats] = field(default_factory=dict)
    per_length_support: Dict[int, int] = field(default_factory=dict)

    # Raw sequences per length class (in-memory only, not part of to_dict --
    # would bloat the JSON sketch output). Needed by the b_adapter boundary
    # estimator, which matches literal adapter sequences against reads and
    # can't work from aggregated per-position stats alone.
    per_length_reads: Dict[int, List[str]] = field(default_factory=dict, repr=False)

    # Mean Phred quality per position, 5'- and 3'-anchored. Empty if the
    # input format carries no quality (FASTA/collapsed).
    quality_5p: List[float] = field(default_factory=list)
    quality_3p: List[float] = field(default_factory=list)

    # Top 3' terminal k-mers (k=4 by default) with their read fraction,
    # most common first. Feeds the M2 de novo k-mer boundary estimator.
    terminal_kmers_3p: List[Tuple[str, float]] = field(default_factory=list)

    def to_dict(self) -> Dict:
        """JSON-serializable summary, for `getRPF sketch` output and reuse
        by downstream boundary estimation / evidence building."""

        def stats_dict(stats: SignalStats) -> Dict:
            return {
                "entropy_5p": stats.entropy_5p,
                "composition_5p": stats.composition_5p,
                "entropy_3p": stats.entropy_3p,
                "composition_3p": stats.composition_3p,
                "sample_size": stats.sample_size,
            }

        return {
            "sample_size": self.sample_size,
            "length_distribution": self.length_distribution,
            "pooled": stats_dict(self.pooled),
            "per_length": {
                str(length): stats_dict(stats)
                for length, stats in self.per_length.items()
            },
            "per_length_support": self.per_length_support,
            "quality_5p": self.quality_5p,
            "quality_3p": self.quality_3p,
            "terminal_kmers_3p": self.terminal_kmers_3p,
        }


class SketchBuilder:
    """Builds a Stage-0 Sketch from a FASTQ/FASTA/collapsed file."""

    def __init__(self, max_reads: int = 500_000, terminal_kmer_len: int = 4):
        self.max_reads = max_reads
        self.terminal_kmer_len = terminal_kmer_len

    def build_from_file(
        self, input_path: Path, format: str, count_pattern: Optional[str] = None
    ) -> Sketch:
        """Load reads (+ qualities, if available) and build a Sketch."""
        reads, qualities = self._load(input_path, format, count_pattern)
        return self.build_from_reads(reads, qualities)

    def build_from_reads(
        self, reads: List[str], qualities: Optional[List[List[int]]] = None
    ) -> Sketch:
        if not reads:
            return Sketch(
                length_distribution={},
                sample_size=0,
                pooled=empty_stats(),
            )

        length_distribution: Dict[int, int] = {}
        for r in reads:
            length_distribution[len(r)] = length_distribution.get(len(r), 0) + 1

        pooled = process_reads(reads, compute_dinucleotide=False)

        by_length: Dict[int, List[str]] = {}
        for r in reads:
            by_length.setdefault(len(r), []).append(r)

        per_length: Dict[int, SignalStats] = {}
        per_length_support: Dict[int, int] = {}
        for length, length_reads in by_length.items():
            per_length_support[length] = len(length_reads)
            per_length[length] = process_reads(
                length_reads, compute_dinucleotide=False
            )

        quality_5p, quality_3p = self._quality_profiles(reads, qualities)
        terminal_kmers_3p = self._terminal_kmers(reads)

        return Sketch(
            length_distribution=length_distribution,
            sample_size=len(reads),
            pooled=pooled,
            per_length=per_length,
            per_length_support=per_length_support,
            per_length_reads=by_length,
            quality_5p=quality_5p,
            quality_3p=quality_3p,
            terminal_kmers_3p=terminal_kmers_3p,
        )

    def _load(
        self, input_path: Path, format: str, count_pattern: Optional[str]
    ) -> Tuple[List[str], Optional[List[List[int]]]]:
        reads: List[str] = []
        qualities: Optional[List[List[int]]] = [] if format == "fastq" else None

        if format == "collapsed":
            if count_pattern is None:
                count_pattern = "seq{id}_x{count}"
            sequences, counts = parse_collapsed_fasta(
                input_path, count_pattern, self.max_reads
            )
            for header, seq in sequences.items():
                count = counts.get(header, 1)
                seq_upper = seq.upper()
                n = min(count, self.max_reads - len(reads))
                if n <= 0:
                    break
                reads.extend([seq_upper] * n)
            return reads, None

        opener = get_file_opener(input_path)
        with opener(str(input_path), "rt") as handle:
            for record in SeqIO.parse(handle, format):
                if len(reads) >= self.max_reads:
                    break
                reads.append(str(record.seq).upper())
                if qualities is not None:
                    qualities.append(
                        record.letter_annotations.get("phred_quality", [])
                    )
        return reads, qualities

    def _quality_profiles(
        self, reads: List[str], qualities: Optional[List[List[int]]]
    ) -> Tuple[List[float], List[float]]:
        if not qualities:
            return [], []

        max_len = max(len(r) for r in reads)
        sums_5p = [0.0] * max_len
        counts_5p = [0] * max_len
        sums_3p = [0.0] * max_len
        counts_3p = [0] * max_len

        for q in qualities:
            n = len(q)
            for pos in range(n):
                sums_5p[pos] += q[pos]
                counts_5p[pos] += 1
                rpos = n - 1 - pos
                sums_3p[rpos] += q[pos]
                counts_3p[rpos] += 1

        quality_5p = [
            sums_5p[i] / counts_5p[i] if counts_5p[i] else 0.0 for i in range(max_len)
        ]
        quality_3p = [
            sums_3p[i] / counts_3p[i] if counts_3p[i] else 0.0 for i in range(max_len)
        ]
        return quality_5p, quality_3p

    def _terminal_kmers(self, reads: List[str]) -> List[Tuple[str, float]]:
        k = self.terminal_kmer_len
        counter: Counter = Counter()
        total = 0
        for r in reads:
            if len(r) < k:
                continue
            counter[r[-k:]] += 1
            total += 1

        if total == 0:
            return []

        return [(kmer, count / total) for kmer, count in counter.most_common(10)]


def build_sketch(
    input_path: Path,
    format: str,
    max_reads: int = 500_000,
    count_pattern: Optional[str] = None,
) -> Sketch:
    """Build a bounded input sketch from a sequence file."""
    return SketchBuilder(max_reads=max_reads).build_from_file(
        input_path, format=format, count_pattern=count_pattern
    )


def classify_terminal_signal(
    entropy: float,
    quality_mean: Optional[float],
    dominant_base: Optional[str],
    entropy_threshold: float = 0.5,
) -> str:
    """Discriminate a constant terminal region from a base-caller artifact.

    A poly-G dark-cycle run (2-colour chemistry no-calls) mimics a constant
    3' "adapter": near-zero entropy, dominant base often G. The discriminator
    is quality: a real adapter/constant region is sequenced at normal
    quality; a dark-cycle run co-occurs with a quality collapse. See
    docs/release_qc_and_terminal_trimming_plan.md section 7.

    Returns one of: "biological", "adapter_or_constant", "basecaller_artifact".
    """
    if entropy >= entropy_threshold:
        return "biological"

    if (
        quality_mean is not None
        and dominant_base == "G"
        and quality_mean < POLYG_QUALITY_COLLAPSE_THRESHOLD
    ):
        return "basecaller_artifact"

    return "adapter_or_constant"
