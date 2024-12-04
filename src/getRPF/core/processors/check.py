"""Cleanliness and cleanliness checking for sequence data.

This module provides functionality for analyzing sequence quality, composition,
and potential contaminants in high-throughput sequencing data.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from Bio import SeqIO

from ...utils.file_utils import get_file_opener
from ..processors.collapsed import parse_collapsed_fasta


@dataclass
class CleanlinessResults:
    """Container for quality analysis results.

    Attributes:
        length_distribution: Distribution of read lengths
        nucleotide_frequencies: Per-position nucleotide frequencies
        quality_scores: Per-position quality scores (FASTQ only)
        gc_content: Overall GC content percentage
        complexity_scores: Sequence complexity measures
    """

    length_distribution: Dict[int, int]
    nucleotide_frequencies: Dict[str, List[float]]
    quality_scores: Optional[Dict[int, List[float]]] = None
    gc_content: Optional[float] = None
    complexity_scores: Optional[Dict[str, float]] = None

    def write_report(self, output_path: Path) -> None:
        """Write analysis results to a file.

        Args:
            output_path: Path where to write the report
        """
        with open(output_path, "w") as f:
            # Write read length distribution
            f.write("=== Read Length Distribution ===\n")
            for length, count in sorted(self.length_distribution.items()):
                f.write(f"{length}\t{count}\n")

            # Write nucleotide frequencies
            f.write("\n=== Nucleotide Frequencies ===\n")
            positions = len(next(iter(self.nucleotide_frequencies.values())))
            f.write("Position\tA\tC\tG\tT\n")
            for pos in range(positions):
                f.write(f"{pos}\t")
                f.write(
                    "\t".join(
                        f"{self.nucleotide_frequencies[nt][pos]:.3f}" for nt in "ACGT"
                    )
                )
                f.write("\n")

            # Write quality scores if available
            if self.quality_scores:
                f.write("\n=== Quality Scores ===\n")
                f.write("Position\tMean\tStd\n")
                for pos, scores in sorted(self.quality_scores.items()):
                    mean = np.mean(scores)
                    std = np.std(scores)
                    f.write(f"{pos}\t{mean:.2f}\t{std:.2f}\n")

            # Write GC content if available
            if self.gc_content is not None:
                f.write("\n=== GC Content ===\n")
                f.write(f"Overall GC%: {self.gc_content:.2f}\n")

            # Write complexity scores if available
            if self.complexity_scores:
                f.write("\n=== Sequence Complexity ===\n")
                for metric, score in self.complexity_scores.items():
                    f.write(f"{metric}: {score:.3f}\n")


class CleanlinessChecker:
    """Analyzes sequence quality and composition."""

    def __init__(
        self,
        format: str,
        min_quality: int = 20,
        threads: int = 1,
        count_pattern: Optional[str] = None,
        max_reads: Optional[int] = None,
    ):
        """Initialize the CleanlinessChecker.

        Args:
            format: Input file format (fastq/fasta/collapsed)
            min_quality: Minimum quality score threshold
            threads: Number of processing threads
        """
        self.format = format
        self.min_quality = min_quality
        self.threads = threads
        self.count_pattern = count_pattern
        self.max_reads = max_reads

    def analyze_file(
        self, input_path: Path, count_pattern: Optional[str] = None
    ) -> CleanlinessResults:
        """Analyze sequence file for quality metrics.

        Args:
            input_path: Path to input sequence file
            count_pattern: For collapsed format, pattern to extract read counts

        Returns:
            CleanlinessResults object containing analysis results
        """
        # Initialize counters
        length_dist = {}
        nuc_freqs = {nt: [] for nt in "ACGT"}
        quality_scores = {} if self.format == "fastq" else None
        gc_count = 0
        total_bases = 0

        if self.format == "collapsed":
            if count_pattern is None:
                count_pattern = "read_{prefix}_{count}"

            sequences, counts = parse_collapsed_fasta(
                input_path, count_pattern, self.max_reads
            )

            for header, seq in sequences.items():
                count = counts.get(header, 1)
                length = len(seq)
                # Update length distribution considering count
                length_dist[length] = length_dist.get(length, 0) + count

                # Update nucleotide frequencies
                while max(len(nuc_freqs[nt]) for nt in "ACGT") < length:
                    for nt in "ACGT":
                        nuc_freqs[nt].append(0)

                for pos, nt in enumerate(seq.upper()):
                    if nt in nuc_freqs:
                        nuc_freqs[nt][pos] += count
                        if nt in "GC":
                            gc_count += count
                total_bases += length * count

        else:
            # Existing code for FASTQ/FASTA processing
            file_opener = get_file_opener(input_path)
            with file_opener(
                str(input_path), "rt"
            ) as handle:  # 'rt' mode for text reading
                record_count = 0
                for record in SeqIO.parse(handle, self.format):
                    # Check if we've hit the max_reads limit
                    if self.max_reads is not None and record_count >= self.max_reads:
                        break

                    length = len(record.seq)
                    length_dist[length] = length_dist.get(length, 0) + 1

                    seq_str = str(record.seq).upper()
                    while max(len(nuc_freqs[nt]) for nt in "ACGT") < length:
                        for nt in "ACGT":
                            nuc_freqs[nt].append(0)

                    for pos, nt in enumerate(seq_str):
                        if nt in nuc_freqs:
                            nuc_freqs[nt][pos] += 1
                            if nt in "GC":
                                gc_count += 1
                    total_bases += length

                    if self.format == "fastq" and hasattr(record, "letter_annotations"):
                        phred_scores = record.letter_annotations.get(
                            "phred_quality", []
                        )
                        for pos, score in enumerate(phred_scores):
                            if pos not in quality_scores:
                                quality_scores[pos] = []
                            quality_scores[pos].append(score)

                    record_count += 1

        # Normalize nucleotide frequencies
        total_reads = sum(length_dist.values())
        for nt in "ACGT":
            nuc_freqs[nt] = [count / total_reads for count in nuc_freqs[nt]]

        # Calculate GC content
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else None

        return CleanlinessResults(
            length_distribution=length_dist,
            nucleotide_frequencies=nuc_freqs,
            quality_scores=quality_scores,
            gc_content=gc_content,
        )
