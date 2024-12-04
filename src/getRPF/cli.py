"""Command Line Interface for getRPF.

This module implements the command-line interface for getRPF
(Get Ribosome Protected Fragment features),
a comprehensive tool for analyzing Ribosome Protected Fragments (RPFs)
from Ribo-seq experiments.

The CLI is built using Click and provides a hierarchical command structure:

Main Commands:
    check: Analyze read quality and nucleotide composition
    detect-adapter: Identify and characterize adapter sequences

Key Features:
    - Supports multiple input formats (FASTQ, FASTA, Collapsed FASTA)
    - Handles both gzipped and uncompressed files
    - Provides detailed quality metrics and visualizations
    - Implements efficient adapter detection algorithms

Examples:
    Basic quality analysis:
        $ getRPF check input.fastq --format fastq --output quality_report.txt

    Adapter detection with custom sequence:
        $ getRPF detect-adapter input.fastq \
            --format fastq --adapter AGATCGGAAGAG \
            --output adapters.txt

Notes:
    - All file paths can be either relative or absolute
    - For gzipped files, compression is automatically detected
    - Memory usage scales with read length, not file size
    - Temporary files are cleaned up automatically

See Also:
    - Documentation: https://getRPF.readthedocs.io
    - Source Code: https://github.com/yourusername/getRPF
    - Bug Reports: https://github.com/yourusername/getRPF/issues
"""

from enum import Enum
from pathlib import Path
from typing import Optional

import click

from .core.handlers import (
    handle_adapter_detection,
    handle_cleanliness_check,
)


class InputFormat(str, Enum):
    """Supported input format types for sequence data.

    This enum defines the valid input formats that getRPF can process.
    Each format has specific characteristics and requirements.

    Attributes:
        FASTQ: Standard FASTQ format
            - Contains sequence and quality scores
            - Four lines per record
            - Quality scores in Phred+33 format

        FASTA: Standard FASTA format
            - Contains sequence only
            - Two lines per record
            - No quality information

        COLLAPSED: Collapsed FASTA format
            - Modified FASTA where headers contain read counts
            - Format: >sequence_count_N
            - Used for deduplicated data

    Example:
        >>> format = InputFormat.FASTQ
        >>> format == "fastq"
        True
        >>> format in InputFormat
        True
    """

    FASTQ = "fastq"
    FASTA = "fasta"
    COLLAPSED = "collapsed"


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """getRPF - Comprehensive Ribosome Protected Fragment Analysis.

    This is the main entry point for the getRPF command-line interface.
    It provides access to various analysis tools for Ribo-seq data processing.

    The tool focuses on:
        - Quality assessment of RPF reads
        - Adapter sequence detection and analysis
        - Read length distribution analysis
        - Nucleotide composition profiling

    For detailed documentation, visit: https://getRPF.readthedocs.io
    """
    pass


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fastq", "fasta", "collapsed"]),
    help="Input file format",
    required=True,
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output file path",
    required=True,
)
@click.option(
    "--count-pattern",
    "-p",
    help="Pattern for extracting read count from collapsed FASTA headers. "
    "Use {count} to mark where the count appears. "
    'Examples: "read_{count}", "read\\d+_x{count}", "{count}_seq"',
    default="read_{count}",
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Maximum number of reads to process. Default is all reads.",
    default=None,
)
def check(
    input_file: Path,
    format: str,
    output: Path,
    count_pattern: Optional[str] = None,
    max_reads: Optional[int] = None,
):
    """Check read quality and composition.

    For collapsed FASTA format, specify how to extract read counts from headers
    using the --count-pattern option. The pattern should include {count} where
    the number appears.

    Examples:
        # Header format: >read_123_500 (count is 500)
        getRPF check input.fasta --format collapsed\
              --count-pattern "read_{prefix}_{count}"

        # Header format: >read1_x100 (count is 100)
        getRPF check input.fasta --format collapsed \
            --count-pattern "read_{prefix}_{count}"

        # Process only first 1000 reads
        getRPF check input.fastq --format fastq \
            --output report.txt --max-reads 1000
    """
    handle_cleanliness_check(
        input_file=input_file,
        format=format,
        output=output,
        count_pattern=count_pattern if format == "collapsed" else None,
        max_reads=max_reads,
    )


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fastq", "fasta", "collapsed"]),
    help="Input file format",
    required=True,
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output file path",
    required=True,
)
@click.option("--adapter", "-a", help="Adapter sequence", required=True)
@click.option(
    "--min-overlap",
    "-m",
    help="Minimum overlap for adapter matching",
    default=10,
    type=int,
)
@click.option(
    "--max-mismatches", "-M", help="Maximum allowed mismatches", default=1, type=int
)
@click.option(
    "--count-pattern",
    "-p",
    help="Pattern for extracting read count from collapsed FASTA headers. "
    "Use {count} to mark where the count appears. "
    'Examples: "read_{count}", "read\\d+_x{count}", "{count}_seq"',
    default="read_{count}",
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Maximum number of reads to process. Default is all reads.",
    default=None,
)
def detect_adapter(
    input_file: Path,
    format: str,
    output: Path,
    adapter: str,
    min_overlap: int = 10,
    max_mismatches: int = 1,
    count_pattern: Optional[str] = None,
    max_reads: Optional[int] = None,
):
    """Detect adapter sequences in reads.

    For collapsed FASTA format, specify how to extract read counts from headers
    using the --count-pattern option. The pattern should include {count} where
    the number appears.

    Examples:
        # Standard FASTQ
        getRPF detect-adapter input.fastq -f fastq -a AGATCGGAAGAG \
            -o report.txt

        # Collapsed FASTA with format >read_500
        getRPF detect-adapter input.fasta -f collapsed -a AGATCGGAAGAG \
            -o report.txt -p "read_{count}"

        # Collapsed FASTA with format >read1_x100
        getRPF detect-adapter input.fasta -f collapsed -a AGATCGGAAGAG \
            -o report.txt -p "read\\d+_x{count}"

        # Process only first 1000 reads
        getRPF detect-adapter input.fastq -f fastq -a AGATCGGAAGAG \
            -o report.txt --max-reads 1000
    """
    handle_adapter_detection(
        input_file=input_file,
        format=format,
        output=output,
        adapter=adapter,
        min_overlap=min_overlap,
        max_mismatches=max_mismatches,
        count_pattern=count_pattern if format == "collapsed" else None,
        max_reads=max_reads,
    )


if __name__ == "__main__":
    cli()
