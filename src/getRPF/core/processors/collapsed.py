"""Handler for collapsed FASTA format processing."""

import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


class CollapsedHeaderParser:
    """Parser for extracting counts from collapsed FASTA headers based on custom formatting."""

    def __init__(self, format: str = "read_{prefix}_{count}"):
        """
        Initialize the parser with a custom header format.

        Args:
            format: A string indicating the layout of the header.
                    Use '{prefix}' for non-count parts and '{count}' for the numeric count.
        """
        if "{count}" not in format:
            raise ValueError("Format must contain '{count}' placeholder.")

        self.format = format
        self.prefix_placeholder = "{prefix}"
        self.count_placeholder = "{count}"

    def extract_count(self, header: str) -> Optional[int]:
        """
        Extract the read count from the given header string based on the initialized format.

        Args:
            header: FASTA header string (without '>').

        Returns:
            The extracted read count or None if the format doesn't match.

        Raises:
            ValueError: If the count cannot be parsed as an integer.
        """
        try:
            # Split the format into parts
            if self.prefix_placeholder in self.format:
                prefix_part, count_part = self.format.split(self.prefix_placeholder)
            else:
                prefix_part, count_part = self.format.split(self.count_placeholder)

            # Remove placeholders and match parts
            prefix = prefix_part.strip("_")
            count_delimiter = count_part.strip("_")

            # Extract the count from the header
            if prefix:
                header = header.replace(prefix, "", 1).strip("_")
            if count_delimiter:
                header = header.rsplit(count_delimiter, 1)[-1]

            return int(header)
        except (IndexError, ValueError) as e:
            logger.warning(
                f"Failed to parse header '{header}' with format '{self.format}': {e}"
            )
            return None


def parse_collapsed_fasta(
    file_path: Path, count_pattern: str = "read_{count}"
) -> Tuple[dict, dict]:
    """Parse collapsed FASTA file with flexible header format.

    Args:
        file_path: Path to collapsed FASTA file
        count_pattern: Pattern for extracting count from headers

    Returns:
        Tuple of (sequences dict, counts dict)
    """
    parser = CollapsedHeaderParser(count_pattern)
    sequences = {}
    counts = {}

    with open(file_path) as f:
        current_header = None
        current_seq = []

        for line in f:
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

                # Extract count
                count = parser.extract_count(current_header)
                if count is not None:
                    counts[current_header] = count
                else:
                    counts[current_header] = 1  # Default to 1 if no count found
            else:
                current_seq.append(line)

        # Process last sequence
        if current_header is not None:
            seq = "".join(current_seq)
            sequences[current_header] = seq

    return sequences, counts
