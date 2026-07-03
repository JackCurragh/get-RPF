"""Handler for collapsed FASTA format processing."""

import bz2
import gzip
import logging
from contextlib import contextmanager
from pathlib import Path
from collections import Counter
from typing import Dict, Optional, TextIO, Tuple, Union

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

    def __init__(self, format: str = "seq{id}_x{count}"):
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
    count_pattern: str = "seq{id}_x{count}",
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
    """Legacy wrapper for backward compatibility."""
    def __init__(self, count_pattern: Optional[str] = None):
        self.collapser = TwoStageCollapser()
        self.count_pattern = count_pattern or "seq{id}_x{count}"
    
    def expand_to_fastq(self, input_file: Path, output_file: Path, max_reads: Optional[int] = None) -> None:
        """Expand collapsed FASTA to FASTQ."""
        raw_counts = self.collapser.collapse_raw(input_file, format="collapsed", max_reads=max_reads)
        with open(output_file, 'w') as fout:
            for idx, (seq, count) in enumerate(raw_counts.items(), 1):
                for i in range(count):
                    fout.write(f"@seq{idx}_c{i+1}\n{seq}\n+\n{'I' * len(seq)}\n")


class TwoStageCollapser:
    """High-performance collapser that implements Two-Stage Collapsing.
    
    Stage 1: Raw collapse (count unique raw reads)
    Stage 2: Unique trimming (trim each unique sequence once)
    Stage 3: Final merge (aggregate counts of identical trimmed sequences)
    """

    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)

    def collapse_raw(
        self, 
        input_file: Path, 
        format: str = "fastq", 
        max_reads: Optional[int] = None
    ) -> Counter:
        """Stage 1: Extract and count raw sequences from file."""
        counts: Counter = Counter()
        processed = 0
        
        self.logger.info(f"Stage 1: Collapsing raw reads from {input_file}...")
        
        file_opener = gzip.open if str(input_file).endswith(".gz") else open
        with file_opener(str(input_file), "rt") as fh:
            if format == "fastq":
                iterator = _iter_fastq(fh)
            elif format in ["fasta", "collapsed"]:
                # For collapsed, we need to respect existing counts
                iterator = self._iter_fasta_or_collapsed(fh, format == "collapsed")
            else:
                raise ValueError(f"Unsupported format: {format}")

            for seq in iterator:
                if isinstance(seq, tuple): # (sequence, count) from collapsed
                    counts[seq[0]] += seq[1]
                else:
                    counts[seq] += 1
                
                processed += 1
                if max_reads and processed >= max_reads:
                    break
                    
        self.logger.info(f"  Processed {processed} reads -> {len(counts)} unique raw sequences")
        return counts

    def _iter_fasta_or_collapsed(self, fh, is_collapsed: bool):
        """Internal iterator for FASTA/Collapsed FASTA."""
        parser = CollapsedHeaderParser() if is_collapsed else None
        current_seq = []
        current_count = 1
        
        for line in fh:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                if current_seq:
                    yield ("".join(current_seq).upper(), current_count)
                
                if is_collapsed:
                    cnt = parser.extract_count(line[1:])
                    current_count = cnt if cnt is not None else 1
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq:
            yield ("".join(current_seq).upper(), current_count)

    def apply_trimming(
        self, 
        raw_counts: Counter, 
        trim_func: callable,
        min_length: int = 20,
        max_length: Optional[int] = None,
    ) -> Counter:
        """Stage 2 & 3: Trim unique sequences and merge results."""
        self.logger.info("Stage 2: Trimming unique sequences and merging...")
        final_counts: Counter = Counter()
        
        for raw_seq, count in raw_counts.items():
            trimmed_seq = trim_func(raw_seq)
            if trimmed_seq and len(trimmed_seq) >= min_length and (
                max_length is None or len(trimmed_seq) <= max_length
            ):
                final_counts[trimmed_seq] += count
                
        self.logger.info(f"  Post-trimming: {len(final_counts)} unique sequences")
        return final_counts

    def write_collapsed_fasta(
        self, 
        counts: Counter, 
        output_file: Path,
        prefix: str = "seq"
    ) -> Dict[str, Union[int, str]]:
        """Stage 4: Write final collapsed results to FASTA."""
        total_reads = sum(counts.values())
        unique = len(counts)

        # Avoid producing empty placeholder files that confuse downstream steps
        if unique == 0 or total_reads == 0:
            self.logger.warning(
                "No RPFs to write (0 unique / 0 reads). Skipping creation of %s",
                output_file,
            )
            return {
                "unique_sequences": 0,
                "total_reads": 0,
                "output_path": None,
            }

        with open(output_file, "w") as fout:
            # Sort by count descending for better readability
            for idx, (seq, count) in enumerate(counts.most_common(), 1):
                fout.write(f">{prefix}{idx}_x{count}\n")
                fout.write(f"{seq}\n")
        
        stats = {
            "unique_sequences": unique,
            "total_reads": total_reads,
            "output_path": str(output_file)
        }
        self.logger.info(f"  Wrote {stats['unique_sequences']} sequences to {output_file}")
        return stats


def collapse_fastq_to_fasta(
    input_file: Path,
    output_file: Path,
) -> Dict[str, Union[int, str]]:
    """Legacy wrapper for backward compatibility using TwoStageCollapser."""
    collapser = TwoStageCollapser()
    raw_counts = collapser.collapse_raw(input_file, format="fastq")
    # No trimming in this basic collapse function
    return collapser.write_collapsed_fasta(raw_counts, output_file)


def _iter_fastq(handle) -> str:
    """Yield uppercase sequences from a FASTQ file handle."""
    line_num = 0
    for line in handle:
        line_num += 1
        if line_num % 4 == 2:  # sequence line
            yield line.strip().upper()
