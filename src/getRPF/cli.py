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
    handle_align_detect,
    handle_extract_rpf,
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
    help="Output directory for categorized reports",
    required=True,
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Maximum number of reads to process for testing",
    default=1000,
)
def check_cleanliness(
    input_file: Path,
    format: str,
    output: Path,
    max_reads: int = 1000,
):
    """Enhanced cleanliness checking with failure categorization.
    
    This command runs comprehensive RPF cleanliness checks and categorizes
    failures by type for batch seqspec generation. Essential for scaling
    to thousands of samples.
    
    The system checks:
    - Read length distribution (RPF size expectations)
    - Information content uniformity (no repetitive sequences)  
    - End bias detection (5'/3' nucleotide bias)
    - Base composition uniformity (no positional bias)
    - GC content within normal range
    
    Results:
    - CLEAN samples: Pass all checks, ready for analysis
    - NEEDS_SEQSPEC samples: Categorized by failure type for batch processing
    
    Examples:
        # Check sample cleanliness with categorization
        getRPF check-cleanliness input.fastq -f fastq -o reports/
        
        # Check collapsed format with limited reads
        getRPF check-cleanliness input.fasta -f collapsed -o reports/ --max-reads 5000
    """
    from .core.processors.check import CleanlinessChecker
    from .core.checkers import (
        run_all_cleanliness_checks, 
        categorize_failures, 
        write_check_report,
        LengthDistributionCheck,
        BaseCompositionCheck,
        GCContentCheck
    )
    
    # Create output directory
    output.mkdir(exist_ok=True)
    
    # Run sequence analysis
    checker = CleanlinessChecker(format=format, max_reads=max_reads)
    sequence_results = checker.analyze_file(input_file)
    
    # Run all cleanliness checks (enhanced version)
    check_results = run_all_cleanliness_checks(sequence_results)
    
    # Also run basic RPF checks for compatibility (.rpf_checks.txt format)
    basic_checks = {
        "Length Distribution": LengthDistributionCheck(),
        "Base Composition": BaseCompositionCheck(),
        "GC Content": GCContentCheck(),
    }
    basic_check_results = {
        name: check.check(sequence_results) for name, check in basic_checks.items()
    }
    
    # Categorize failures
    categories = categorize_failures(check_results)
    
    # Write detailed report (enhanced format)
    report_path = output / f"{input_file.stem}_cleanliness_report.txt"
    write_check_report(check_results, report_path)
    
    # Write basic RPF check report (.rpf_checks.txt format for compatibility)
    rpf_report_path = output / f"{input_file.stem}_cleanliness_report.rpf_checks.txt"
    write_check_report(basic_check_results, rpf_report_path)
    
    # Print summary
    if categories['is_clean']:
        click.echo(f"✅ CLEAN: {input_file.name} passed all cleanliness checks")
        click.echo(f"📄 Report: {report_path}")
    else:
        click.echo(f"❌ NEEDS_SEQSPEC: {input_file.name}")
        click.echo(f"🔍 Primary failure: {categories['primary_failure']}")
        click.echo(f"📋 All failures: {', '.join(categories['failure_categories'])}")
        click.echo(f"📄 Report: {report_path}")
        click.echo(f"💡 Batch process with similar failures for seqspec generation")


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
              --count-pattern "read_{id}_{count}"

        # Header format: >read1_x100 (count is 100)
        getRPF check input.fasta --format collapsed \
            --count-pattern "read_{id}_{count}"

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


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--star-index",
    "-s",
    type=click.Path(exists=True, path_type=Path),
    help="Path to STAR index directory",
    required=True,
)
@click.option(
    "--star-threads",
    "-t",
    type=int,
    help="Number of threads for STAR alignment",
    default=1,
)
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
    "--output-format",
    "-of",
    type=click.Choice(["json", "csv"]),
    help="Output file format",
    default="json",
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
    "--save-bam",
    "-b",
    type=click.Path(path_type=Path),
    help="Path to save the alignment BAM file",
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Maximum number of reads to process. Default is 100k.",
    default=100000,
)
def align_detect(
    input_file: Path,
    star_index: Path,
    star_threads: int,
    format: str,
    output: Path,
    output_format: str = "json",
    count_pattern: Optional[str] = None,
    save_bam: Optional[Path] = None,
    max_reads: Optional[int] = 100000,
):
    """Align reads with STAR and detect features.

    This command performs STAR alignment followed by feature detection
    on ribosome profiling reads. Essential for RPF extraction pipeline.

    Examples:
        # Align FASTQ reads and output JSON
        getRPF align-detect input.fastq -s /path/to/star/index -f fastq -o output.json

        # Align Collapsed FASTA reads with custom count pattern
        getRPF align-detect input.fasta -s /path/to/star/index -f collapsed -o output.csv \
            -p "read_{count}"

        # Process only first 5000 reads
        getRPF align-detect input.fastq -s /path/to/star/index -f fastq -o output.json --max-reads 5000
    """
    handle_align_detect(
        input_file=input_file,
        star_index=star_index,
        star_threads=star_threads,
        format=format,
        output=output,
        output_format=output_format,
        count_pattern=count_pattern if format == "collapsed" else None,
        save_bam_path=save_bam,
        max_reads=max_reads,
    )


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fastq", "fasta", "collapsed"]),
    help="Input file format",
    required=True,
)
@click.option(
    "--architecture-db",
    "-a",
    type=click.Path(exists=True, path_type=Path),
    help="Path to custom architecture database (JSON file)",
)
@click.option(
    "--seqspec-dir",
    "-s",
    type=click.Path(exists=True, path_type=Path),
    help="Directory containing seqspec files for novel protocols",
)
@click.option(
    "--generate-seqspec",
    "-g",
    is_flag=True,
    help="Generate seqspec file for detected architecture",
)
@click.option(
    "--output-format",
    "-of",
    type=click.Choice(["json", "csv"]),
    help="Format for detection report",
    default="json",
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Maximum number of reads to process",
    default=None,
)
def detect_architecture(
    input_file: Path,
    output_file: Path,
    format: str,
    architecture_db: Optional[Path] = None,
    seqspec_dir: Optional[Path] = None,
    generate_seqspec: bool = False,
    output_format: str = "json",
    max_reads: Optional[int] = None,
):
    """Detect read architecture and optionally extract clean RPFs.

    This command analyzes ribosome profiling reads to identify their structure
    using pattern matching against known architectures or de novo detection.
    The primary output is a seqspec file that describes the detected architecture.

    The system automatically:
    - Detects read architecture (UMI, barcodes, adapters) using known patterns
    - Falls back to de novo detection for unknown protocols  
    - Generates seqspec files for discovered architectures
    - Optionally extracts clean RPF sequences
    - Provides detailed detection reports

    Examples:
        # Detect architecture and generate seqspec
        getRPF detect-architecture input.fastq output_rpfs.fastq -f fastq --generate-seqspec

        # Use custom architecture database  
        getRPF detect-architecture input.fastq output_rpfs.fastq -f fastq \\
            -a custom_architectures.json --generate-seqspec

        # Load novel protocols from seqspec directory
        getRPF detect-architecture input.fastq output_rpfs.fastq -f fastq \\
            --seqspec-dir my_protocols/ --generate-seqspec

        # Architecture detection only (no RPF extraction)
        getRPF detect-architecture input.fastq temp_output.fastq -f fastq \\
            --generate-seqspec --max-reads 5000

        # Process collapsed FASTA with processing limit
        getRPF detect-architecture input.fasta output_rpfs.fasta -f collapsed \\
            --generate-seqspec --max-reads 50000
    """
    handle_extract_rpf(
        input_file=input_file,
        output_file=output_file,
        format=format,
        architecture_db=architecture_db,
        seqspec_dir=seqspec_dir,
        generate_seqspec=generate_seqspec,
        output_format=output_format,
        max_reads=max_reads,
    )



@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option("--star-index", "-s", type=click.Path(exists=True, path_type=Path), required=True)
@click.option("--format", "-f", type=click.Choice(["fastq", "fasta", "collapsed"]), required=True)
@click.option("--output", "-o", type=click.Path(path_type=Path), required=True)
@click.option("--max-reads", "-n", type=int, default=10000)
def decide_trim(input_file, star_index, format, output, max_reads):
    """Auto-configure trimming by combining methods.
    
    Runs both architecture detection and alignment verification to determine
    the optimal trimming parameters.
    """
    from .core.handlers import handle_decide_trim
    
    handle_decide_trim(
        input_file=input_file,
        star_index=star_index,
        format=format,
        output=output,
        max_reads=max_reads
    )

if __name__ == "__main__":
    cli()
