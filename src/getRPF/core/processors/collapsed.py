"""Handler for collapsed FASTA format processing."""

import bz2
import gzip
import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Optional, TextIO, Tuple, Union

from getRPF.utils.file_utils import check_file_readability

logger = logging.getLogger(__name__)


@contextmanager
def fasta_opener(file_path: Path) -> TextIO:
    """Context manager for reading FASTA files (compressed or uncompressed).

    Args:
        file_path: Path to FASTA file (can be .gz or .bz2)

    Yields:
        File handle for reading

    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file cannot be read
    """
    check_file_readability(file_path)

    suffix = file_path.suffix.lower()
    if suffix == ".gz":
        handle = gzip.open(file_path, "rt", encoding="utf-8")
    elif suffix == ".bz2":
        handle = bz2.open(file_path, "rt", encoding="utf-8")
    else:
        handle = open(file_path, "r", encoding="utf-8")

    try:
        yield handle
    finally:
        handle.close()


class CollapsedHeaderParser:
    """Parser for extracting counts from collapsed FASTA headers based on custom formatting."""

    def __init__(self, format: str = "read_{id}_{count}"):
        """
        Initialize the parser with a custom header format.

        Args:
            format: A string indicating the layout of the header.
                    Use '{id}' for non-count parts and '{count}' for the numeric count.
        """
        if "{count}" not in format:
            raise ValueError("Format must contain '{count}' placeholder")

        # Split on {count} to get left and right delimiters
        self.left_part, *right_parts = format.split("{count}")
        self.right_delim = right_parts[0] if right_parts else ""

        # Get the delimiter before the count
        if "{id}" in self.left_part:
            # If there's an {id}, split on that and take what's after it
            *_, self.left_delim = self.left_part.split("{id}")
        else:
            # Otherwise use the whole left part
            self.left_delim = self.left_part

    def extract_count(self, header: str) -> Optional[int]:
        """
        Extract the read count from the given header string based on the initialized format.

        Args:
            header: FASTA header string (without '>').

        Returns:
            The extracted read count or None if the format doesn't match.
        """
        try:
            # Split based on left delimiter if it exists
            if self.left_delim:
                _, count_part = header.rsplit(self.left_delim, 1)
            else:
                count_part = header

            # Split based on right delimiter if it exists
            if self.right_delim:
                count_str, _ = count_part.split(self.right_delim, 1)
            else:
                count_str = count_part

            return int(count_str.strip("_"))
        except (ValueError, IndexError):
            return None


def parse_collapsed_fasta(
    file_path: Union[str, Path],
    count_pattern: str = "read_{count}",
    max_reads: Optional[int] = None,
) -> Tuple[dict, dict]:
    """Parse collapsed FASTA file with flexible header format.

    Supports both compressed (.gz, .bz2) and uncompressed FASTA files.

    Args:
        file_path: Path to collapsed FASTA file
        count_pattern: Pattern for extracting count from headers
        max_reads: Maximum number of FASTA entries to process. None means process all entries.
                  Note: Each entry may represent multiple reads in collapsed format.

    Returns:
        Tuple of (sequences dict, counts dict)

    Example:
        >>> seqs, counts = parse_collapsed_fasta('reads.fa.gz')
        >>> all(isinstance(count, int) for count in counts.values())
        True
    """
    if isinstance(file_path, str):
        file_path = Path(file_path)

    parser = CollapsedHeaderParser(count_pattern)
    sequences = {}
    counts = {}
    entries_processed = 0

    with fasta_opener(file_path) as f:
        current_header = None
        current_seq = []
        partial_line = ""

        for line in f:
            # Check if we've hit the max entries
            if max_reads is not None and entries_processed >= max_reads:
                break

            # Handle any partial line from previous iteration
            if partial_line:
                line = partial_line + line
                partial_line = ""

            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Process previous sequence if exists
                if current_header is not None:
                    seq = "".join(current_seq)
                    sequences[current_header] = seq

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
                entries_processed += 1

                # Extract count
                count = parser.extract_count(current_header)
                if count is not None:
                    counts[current_header] = count
                else:
                    counts[current_header] = 1  # Default to 1 if no count found
            else:
                current_seq.append(line)

        # Process last sequence if we haven't hit max entries
        if current_header is not None and (
            max_reads is None or entries_processed < max_reads
        ):
            seq = "".join(current_seq)
            sequences[current_header] = seq

    return sequences, counts


class CollapsedFASTAProcessor:
    """Processor for collapsed FASTA format files."""
    
    def __init__(self, count_pattern: Optional[str] = None):
        """Initialize processor with count pattern.
        
        Args:
            count_pattern: Pattern for extracting counts from headers
        """
        self.count_pattern = count_pattern or "read_{count}"
        self.parser = CollapsedHeaderParser(self.count_pattern)
    
    def expand_to_fastq(
        self, 
        input_file: Path, 
        output_file: Path, 
        max_reads: Optional[int] = None
    ) -> None:
        """Expand collapsed FASTA to individual FASTQ entries.
        
        Args:
            input_file: Input collapsed FASTA file
            output_file: Output FASTQ file
            max_reads: Maximum reads to process
        """
        sequences, counts = parse_collapsed_fasta(
            input_file, 
            count_pattern=self.count_pattern, 
            max_reads=max_reads
        )
        
        total_written = 0
        with open(output_file, 'w') as fout:
            for header, sequence in sequences.items():
                count = counts.get(header, 1)
                
                # Write each sequence 'count' times
                for i in range(count):
                    if max_reads and total_written >= max_reads:
                        break
                        
                    # Create unique FASTQ header
                    fastq_header = f"@{header}_copy_{i+1}"
                    quality = 'I' * len(sequence)  # High quality scores
                    
                    fout.write(f"{fastq_header}\n")
                    fout.write(f"{sequence}\n")
                    fout.write("+\n")
                    fout.write(f"{quality}\n")
                    
                    total_written += 1
                
                if max_reads and total_written >= max_reads:
                    break
