"""Core handlers for getRPF functionality.

This module implements the main processing logic for getRPF commands, serving as an
intermediary layer between the CLI and the core processing modules.

The handlers in this module follow a consistent pattern:
1. Validate inputs and options
2. Initialize appropriate processor objects
3. Process data in a memory-efficient manner
4. Aggregate and format results
5. Clean up temporary resources
"""

import logging
from pathlib import Path
from typing import Optional

from ..core.checkers import (
    BaseCompositionCheck,
    GCContentCheck,
    LengthDistributionCheck,
    write_check_report,
)
from ..core.processors.adapter import AdapterDetector
from ..core.processors.check import CleanlinessChecker
from ..core.processors.alignment import STARAligner
from ..core.processors.rpf_extractor import RPFExtractor
from ..utils.validation import validate_adapter_sequence

# Configure logging
logger = logging.getLogger(__name__)


def handle_cleanliness_check(
    input_file: Path,
    format: str,
    output: Path,
    count_pattern: Optional[str] = None,
    min_quality: Optional[int] = 20,
    threads: Optional[int] = 1,
    rpf_checks: bool = True,
    max_reads: Optional[int] = None,
) -> None:
    """Handle the quality check command workflow.

    Args:
        input_file: Path to input sequence file
        format: Format of input file (fastq/fasta/collapsed)
        output: Path where report will be written
        count_pattern: For collapsed format, pattern to extract read counts
        min_quality: Minimum quality score threshold
        threads: Number of processing threads
        rpf_checks: Whether to run RPF-specific checks
    """
    logger.info(f"Starting cleanliness check for {input_file}")
    try:
        # Initialize checker with parameters
        checker = CleanlinessChecker(
            format=format,
            min_quality=min_quality,
            threads=threads,
            max_reads=max_reads,
        )

        # Process file and get basic results
        results = checker.analyze_file(
            input_file,
            count_pattern=count_pattern if format == "collapsed" else None,
        )

        # Write basic report
        results.write_report(output)

        # Run RPF-specific checks if requested
        if rpf_checks:
            checks = {
                "Length Distribution": LengthDistributionCheck(),
                "Base Composition": BaseCompositionCheck(),
                "GC Content": GCContentCheck(),
            }

            # Run checks
            check_results = {
                name: check.check(results) for name, check in checks.items()
            }

            # Write RPF check report
            rpf_report = output.with_suffix(".rpf_checks.txt")
            write_check_report(check_results, rpf_report)
            logger.info(f"RPF check report written to {rpf_report}")

        logger.info("Cleanliness check completed successfully")

    except Exception as e:
        logger.error(f"Cleanliness check failed: {str(e)}")
        raise RuntimeError(f"Cleanliness check failed: {str(e)}") from e


def handle_adapter_detection(
    input_file: Path,
    format: str,
    output: Path,
    adapter: str,
    min_overlap: Optional[int] = 10,
    max_mismatches: Optional[int] = 1,
    threads: Optional[int] = 1,
    count_pattern: Optional[str] = None,
    max_reads: Optional[int] = None,
) -> None:
    """Handle the adapter detection command workflow.

    Args:
        input_file: Path to input sequence file
        format: Format of input file (fastq/fasta/collapsed)
        output: Path where report will be written
        adapter: Sequence of adapter to detect
        min_overlap: Minimum overlap for adapter detection
        max_mismatches: Maximum allowed mismatches
        threads: Number of processing threads
        count_pattern: Pattern to count in header
    """
    logger.info(f"Starting adapter detection for {input_file}")

    try:
        # Validate adapter sequence
        validate_adapter_sequence(adapter)

        # Initialize detector with parameters
        detector = AdapterDetector(
            adapter=adapter,
            format=format,
            min_overlap=min_overlap,
            max_mismatches=max_mismatches,
            threads=threads,
            count_pattern=count_pattern,
            max_reads=max_reads,
        )

        # Process file and get results
        results = detector.analyze_file(input_file)

        # Write report
        results.write_report(output)

        logger.info("Adapter detection completed successfully")

    except Exception as e:
        logger.error(f"Adapter detection failed: {str(e)}")
        raise RuntimeError(f"Adapter detection failed: {str(e)}") from e


def handle_align_detect(
    input_file: Path,
    star_index: Path,
    star_threads: int,
    format: str,
    output: Path,
    output_format: str = "json",
    count_pattern: Optional[str] = None,
    max_reads: Optional[int] = None,
) -> None:
    """Handle the align and detect command workflow.

    Performs STAR alignment followed by feature detection on aligned reads.
    This is the foundation for RPF extraction and analysis.

    Args:
        input_file: Path to input sequence file
        star_index: Path to STAR index directory
        star_threads: Number of threads for STAR alignment
        format: Format of input file (fastq/fasta/collapsed)
        output: Path where results will be written
        output_format: Output format (json/csv)
        count_pattern: Pattern to extract read counts from collapsed headers
        max_reads: Maximum number of reads to process
    """
    logger.info(f"Starting STAR alignment and feature detection for {input_file}")

    try:
        # Initialize STAR aligner with parameters
        aligner = STARAligner(
            star_index=star_index,
            threads=star_threads,
            max_reads=max_reads,
        )

        # Process input file and perform alignment
        alignment_results = aligner.align_reads(
            input_file,
            format=format,
            count_pattern=count_pattern if format == "collapsed" else None,
        )

        # Write results in requested format
        alignment_results.write_report(output, format=output_format)

        logger.info(f"Alignment and feature detection completed successfully")
        logger.info(f"Results written to {output}")

    except Exception as e:
        logger.error(f"Alignment and detection failed: {str(e)}")
        raise RuntimeError(f"Alignment and detection failed: {str(e)}") from e


def handle_extract_rpf(
    input_file: Path,
    output_file: Path,
    format: str,
    architecture_db: Optional[Path] = None,
    seqspec_dir: Optional[Path] = None,
    generate_seqspec: bool = False,
    output_format: str = "json",
    max_reads: Optional[int] = None,
) -> None:
    """Handle the RPF extraction command workflow.

    Automatically extracts ribosome protected fragments from raw sequencing
    reads using pattern matching or de novo detection.

    Args:
        input_file: Path to input sequence file
        output_file: Path for extracted RPF sequences
        format: Format of input file (fastq/fasta/collapsed)
        architecture_db: Path to custom architecture database
        seqspec_dir: Directory containing seqspec files for novel protocols
        generate_seqspec: Whether to generate seqspec file for detected architecture
        output_format: Format for extraction report (json/csv)
        max_reads: Maximum number of reads to process
    """
    logger.info(f"Starting RPF extraction from {input_file}")

    try:
        # Initialize RPF extractor with seqspec directory support
        if seqspec_dir and architecture_db:
            raise ValueError("Cannot specify both --architecture-db and --seqspec-dir")
        
        extractor = RPFExtractor(architecture_db_path=architecture_db, seqspec_dir=seqspec_dir)

        # Extract RPFs
        results = extractor.extract_rpfs(
            input_file=input_file,
            output_file=output_file,
            format=format,
            max_reads=max_reads,
            generate_seqspec=generate_seqspec,
        )

        # Seqspec generation is handled internally by the extractor
        # based on the generate_seqspec parameter passed to extract_rpfs()
        if generate_seqspec and results.seqspec_data:
            seqspec_path = output_file.with_suffix(".seqspec.yaml")
            logger.info(f"Generated seqspec file: {seqspec_path}")

        # Write extraction report
        report_path = output_file.with_suffix(f".extraction_report.{output_format}")
        results.write_report(report_path, format=output_format)

        logger.info(f"RPF extraction completed successfully")
        logger.info(f"Extracted {results.extracted_rpfs} RPFs from {results.input_reads} reads")
        logger.info(f"Extraction method: {results.extraction_method}")
        if results.architecture_match:
            logger.info(f"Matched architecture: {results.architecture_match}")
        logger.info(f"Results written to {output_file}")
        logger.info(f"Report written to {report_path}")

    except Exception as e:
        logger.error(f"RPF extraction failed: {str(e)}")
        raise RuntimeError(f"RPF extraction failed: {str(e)}") from e
