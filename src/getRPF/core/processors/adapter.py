"""Adapter detection and analysis for sequence data.

This module provides functionality for identifying and characterizing
adapter sequences in high-throughput sequencing data, including support
for collapsed read formats.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

from Bio import SeqIO

from ..processors.collapsed import parse_collapsed_fasta


@dataclass
class AdapterResults:
    """Container for adapter detection results.

    Attributes:
        adapter_positions: Position-wise adapter occurrence counts
        partial_matches: Counts of partial adapter matches
        variant_counts: Counts of common adapter variants
        contamination_rate: Overall adapter contamination rate
    """

    adapter_positions: Dict[int, int]
    partial_matches: Dict[int, int]  # key is match length
    variant_counts: Dict[str, int]
    contamination_rate: float

    def write_report(self, output_path: Path) -> None:
        """Write adapter analysis results to a file.

        Args:
            output_path: Path where to write the report
        """
        with open(output_path, "w") as f:
            # Write overall statistics
            f.write("=== Adapter Contamination Summary ===\n")
            f.write(f"Overall contamination rate: {self.contamination_rate:.2f}%\n\n")

            # Write position-wise occurrence
            f.write("=== Adapter Positions ===\n")
            f.write("Position\tCount\n")
            for pos, count in sorted(self.adapter_positions.items()):
                f.write(f"{pos}\t{count}\n")

            # Write partial match distribution
            f.write("\n=== Partial Matches ===\n")
            f.write("Match Length\tCount\n")
            for length, count in sorted(self.partial_matches.items()):
                f.write(f"{length}\t{count}\n")

            # Write variant frequencies
            f.write("\n=== Common Variants ===\n")
            f.write("Variant\tCount\n")
            for variant, count in sorted(
                self.variant_counts.items(), key=lambda x: x[1], reverse=True
            ):
                f.write(f"{variant}\t{count}\n")


class AdapterDetector:
    """Detects and analyzes adapter sequences."""

    def __init__(
        self,
        adapter: str,
        format: str,
        min_overlap: int = 10,
        max_mismatches: int = 1,
        threads: int = 1,
        count_pattern: Optional[str] = None,
        max_reads: Optional[int] = None,
    ):
        """Initialize the AdapterDetector.

        Args:
            adapter: Adapter sequence to search for
            format: Input file format (fastq/fasta/collapsed)
            min_overlap: Minimum overlap for adapter matching
            max_mismatches: Maximum allowed mismatches
            threads: Number of processing threads
            count_pattern: For collapsed format, pattern to extract read counts
        """
        self.adapter = adapter.upper()
        self.format = format
        self.min_overlap = min_overlap
        self.max_mismatches = max_mismatches
        self.threads = threads
        self.count_pattern = count_pattern
        self.max_reads = max_reads

    def _count_mismatches(self, seq1: str, seq2: str) -> int:
        """Count mismatches between two sequences.

        Args:
            seq1: First sequence
            seq2: Second sequence

        Returns:
            Number of mismatches
        """
        return sum(1 for a, b in zip(seq1, seq2) if a != b)

    def _find_best_match(
        self, sequence: str
    ) -> Tuple[Optional[int], Optional[int], Optional[str]]:
        """Find the best adapter match in a sequence.

        Args:
            sequence: Sequence to search in

        Returns:
            Tuple of (position, match_length, variant)
        """
        best_pos = None
        best_length = None
        best_variant = None
        min_mismatches = float("inf")

        seq_len = len(sequence)
        adapter_len = len(self.adapter)

        # Check each possible position
        for pos in range(seq_len - self.min_overlap + 1):
            # Check different match lengths
            for length in range(self.min_overlap, min(adapter_len, seq_len - pos) + 1):
                seq_part = sequence[pos : pos + length]
                adapter_part = self.adapter[:length]

                mismatches = self._count_mismatches(seq_part, adapter_part)

                if mismatches <= self.max_mismatches and mismatches < min_mismatches:
                    min_mismatches = mismatches
                    best_pos = pos
                    best_length = length
                    best_variant = seq_part

        return best_pos, best_length, best_variant

    def analyze_file(self, input_path: Path) -> AdapterResults:
        """Analyze sequence file for adapter content.

        Args:
            input_path: Path to input sequence file

        Returns:
            AdapterResults object containing analysis results
        """
        # Initialize counters
        positions: Dict[int, int] = {}
        partial_matches: Dict[int, int] = {}
        variants: Dict[str, int] = {}
        total_reads = 0
        contaminated_reads = 0

        if self.format == "collapsed":
            # Handle collapsed format
            if self.count_pattern is None:
                self.count_pattern = "read_{id}_{count}"

            sequences, counts = parse_collapsed_fasta(
                input_path, self.count_pattern, self.max_reads
            )

            for header, sequence in sequences.items():
                count = counts.get(header, 1)
                total_reads += count
                sequence = sequence.upper()

                # Find best adapter match
                pos, length, variant = self._find_best_match(sequence)

                if pos is not None:
                    contaminated_reads += count
                    positions[pos] = positions.get(pos, 0) + count
                    partial_matches[length] = partial_matches.get(length, 0) + count
                    if variant:
                        variants[variant] = variants.get(variant, 0) + count

        else:
            # Handle standard FASTQ/FASTA formats
            for record in SeqIO.parse(str(input_path), self.format):
                total_reads += 1
                sequence = str(record.seq).upper()

                # Find best adapter match
                pos, length, variant = self._find_best_match(sequence)

                if pos is not None:
                    contaminated_reads += 1
                    positions[pos] = positions.get(pos, 0) + 1
                    partial_matches[length] = partial_matches.get(length, 0) + 1
                    if variant:
                        variants[variant] = variants.get(variant, 0) + 1

        # Calculate contamination rate
        contamination_rate = (
            contaminated_reads / total_reads * 100 if total_reads > 0 else 0
        )

        return AdapterResults(
            adapter_positions=positions,
            partial_matches=partial_matches,
            variant_counts=variants,
            contamination_rate=contamination_rate,
        )
