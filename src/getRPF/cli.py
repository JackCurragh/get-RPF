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
    handle_align_detect,
    handle_cleanliness_check,
)

# NOTE: Plotting modules pull in heavy optional deps (matplotlib, seaborn).
# Import them lazily inside the specific subcommands so that core commands
# like `extract`, `check`, and `align-detect` do not require those packages.
# Avoid importing DuckDB client at module import time; only needed by ingest-duckdb


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
@click.version_option(package_name="getRPF")
def cli():
    """getRPF - Comprehensive Ribosome Protected Fragment Analysis.

    This is the main entry point for the getRPF command-line interface.
    It provides access to various analysis tools for Ribo-seq data processing.

    The tool focuses on:
        - RPF extraction with single-nucleotide precision
        - Quality assessment of RPF reads
        - Adapter sequence detection and analysis
        - Read length distribution analysis
        - Nucleotide composition profiling

    For detailed documentation, visit: https://getRPF.readthedocs.io
    """
    pass


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fastq", "fasta", "collapsed"]),
    help="Input file format",
    default="fastq",
    show_default=True,
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
@click.option("--generate-seqspec", "-g", is_flag=True)
@click.option(
    "--output-format",
    "-of",
    type=click.Choice(["json", "csv"]),
    default="json",
    show_default=True,
)
@click.option("--max-reads", "-n", type=int, default=None)
@click.option("--star-index", type=click.Path(exists=True, path_type=Path))
@click.option("--star-threads", type=int, default=None, help="Threads for STAR verification")
@click.option("--threads", type=int, default=1, help="Alias for --star-threads")
@click.option("--output-report", type=click.Path(path_type=Path), help="Legacy report path alias")
@click.option("--collapse/--no-collapse", default=True)
@click.option("--collapsed-only", is_flag=True)
@click.option(
    "--infer-reads",
    type=int,
    default=500_000,
    show_default=True,
    help="Subsample depth for boundary-estimation inference (Stage 0/1)",
)
@click.option(
    "--apply-boundary-trims/--no-apply-boundary-trims",
    default=True,
    help="Apply the M2/M3 per-length boundary-estimated trims (default: on)",
)
@click.option(
    "--audit-only",
    is_flag=True,
    help="Infer + screen, apply nothing (rule development / new-family validation)",
)
@click.option(
    "--rules",
    type=click.Path(exists=True, path_type=Path),
    help="Override mode: JSON file of explicit frozen per-length trim rules",
)
@click.option(
    "--fail-on",
    type=click.Choice(["none", "review", "hold"]),
    default="none",
    show_default=True,
    help="Process exit-code threshold",
)
def extract(
    input_file: Path,
    output_file: Path,
    format: str = "fastq",
    architecture_db: Optional[Path] = None,
    seqspec_dir: Optional[Path] = None,
    generate_seqspec: bool = False,
    output_format: str = "json",
    max_reads: Optional[int] = None,
    star_index: Optional[Path] = None,
    star_threads: Optional[int] = None,
    threads: int = 1,
    output_report: Optional[Path] = None,
    collapse: bool = True,
    collapsed_only: bool = False,
    infer_reads: int = 500_000,
    apply_boundary_trims: bool = True,
    audit_only: bool = False,
    rules: Optional[Path] = None,
    fail_on: str = "none",
):
    """Extract trimmed reads without applying the final RPF length gate."""
    from shutil import move

    from .core.handlers import handle_extract_rpf
    from .core.samplesheet import compute_exit_code

    evidence = handle_extract_rpf(
        input_file=input_file,
        output_file=output_file,
        format=format,
        architecture_db=architecture_db,
        seqspec_dir=seqspec_dir,
        generate_seqspec=generate_seqspec,
        output_format=output_format,
        max_reads=max_reads,
        star_index=star_index,
        star_threads=star_threads if star_threads is not None else threads,
        collapse_output=collapse,
        collapsed_only=collapsed_only,
        infer_reads=infer_reads,
        apply_boundary_trims=apply_boundary_trims and not audit_only,
        rules=rules,
    )

    if output_report:
        default_report = output_file.with_suffix(f".extraction_report.{output_format}")
        if default_report != output_report and default_report.exists():
            move(str(default_report), str(output_report))

    exit_code = compute_exit_code(
        [evidence["release_class"]] if evidence else [], fail_on=fail_on
    )
    if exit_code:
        raise SystemExit(exit_code)


cli.add_command(extract, "extract-rpf")


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fastq", "fasta", "collapsed"]),
    help="Input file format",
    default="fastq",
    show_default=True,
)
@click.option(
    "--count-pattern",
    "-p",
    help="Pattern for extracting read count from collapsed FASTA headers.",
    default="seq{id}_x{count}",
)
@click.option(
    "--max-reads",
    "-n",
    type=int,
    help="Subsample depth for Stage-0 sketch inference",
    default=500_000,
    show_default=True,
)
def sketch(
    input_file: Path,
    output_file: Path,
    format: str,
    count_pattern: str,
    max_reads: int,
):
    """Stage-0 sketch: reference-free per-position/per-length diagnostics.

    Computes length distribution, coverage-normalized 5'/3' entropy and
    composition (pooled and per read-length class), a per-position quality
    profile, and top 3' terminal k-mers, on a bounded subsample. This is
    the substrate every later stage (boundary estimation, identity screen,
    release classification) reads from; see
    docs/release_qc_and_terminal_trimming_plan.md sections 4-5.

    Examples:
        getRPF sketch input.fastq sketch.json --format fastq
    """
    import json

    from .core.processors.sketch import build_sketch

    result = build_sketch(
        input_file,
        format=format,
        max_reads=max_reads,
        count_pattern=count_pattern,
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(result.to_dict(), f, indent=2)

    click.echo(f"✅ Wrote sketch ({result.sample_size} reads): {output_file}")


@cli.command()
@click.argument("samplesheet", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output", "-o", "output_dir",
    type=click.Path(path_type=Path), required=True,
    help="Output directory for per-sample outputs and cohort roll-ups",
)
@click.option(
    "--infer-reads", type=int, default=500_000, show_default=True,
    help="Subsample depth for boundary-estimation inference (Stage 0/1)",
)
@click.option("--max-reads", type=int, default=None, help="Cap on reads streamed per sample")
@click.option(
    "--audit-only", is_flag=True,
    help="Infer + screen, apply nothing (rule development / new-family validation)",
)
@click.option(
    "--rules", type=click.Path(exists=True, path_type=Path),
    help="Override mode: JSON file of {sample_id: {length: {trim_5p, trim_3p}}}",
)
@click.option(
    "--fail-on", type=click.Choice(["none", "review", "hold"]), default="none",
    show_default=True, help="Process exit-code threshold",
)
@click.option("--collapse/--no-collapse", default=True)
def run(
    samplesheet: Path,
    output_dir: Path,
    infer_reads: int,
    max_reads: Optional[int],
    audit_only: bool,
    rules: Optional[Path],
    fail_on: str,
    collapse: bool,
):
    """Cohort driver: run the full pipeline over a local-FASTQ samplesheet.

    Minimum samplesheet columns: sample_id,fastq_1. Reuses the same
    per-sample pipeline as `extract`, so audit and production modes agree
    by construction. Writes per-sample outputs plus cohort roll-up TSVs
    (release_qc_summary.tsv, release_qc_flags.tsv, adapter_protocol_summary.tsv,
    samples_needing_seqspec.tsv, samples_excluded_or_held.tsv).
    """
    from .core.samplesheet import run_cohort

    result = run_cohort(
        samplesheet, output_dir,
        audit_only=audit_only, rules_path=rules,
        infer_reads=infer_reads, max_reads=max_reads, collapse=collapse,
    )

    for e in result.evidence:
        click.echo(f"{e['sample_id']}: {e['release_class']}")
    click.echo(f"✅ Processed {len(result.evidence)} sample(s): {output_dir}")
    for name, path in result.cohort_tsv_paths.items():
        click.echo(f"  {name}: {path}")

    exit_code = result.exit_code(fail_on=fail_on)
    if exit_code:
        raise SystemExit(exit_code)


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
    "--count-pattern",
    "-p",
    help="Pattern for extracting read count from collapsed FASTA headers. "
    "Use {count} to mark where the count appears. "
    'Examples: "seq{id}_x{count}", "read_{id}_{count}"',
    default="seq{id}_x{count}",
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
    count_pattern: Optional[str] = None,
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
    from .core.checkers import (
        BaseCompositionCheck,
        GCContentCheck,
        LengthDistributionCheck,
        categorize_failures,
        run_all_cleanliness_checks,
        write_check_report,
    )
    from .core.processors.check import analyze_file

    # Create output directory
    output.mkdir(exist_ok=True)

    # Run sequence analysis
    sequence_results = analyze_file(
        input_file,
        format=format,
        max_reads=max_reads,
        count_pattern=count_pattern if format == "collapsed" else None,
    )

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
        click.echo("💡 Batch process with similar failures for seqspec generation")


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
    default="seq{id}_x{count}",
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
    default="seq{id}_x{count}",
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
    default="seq{id}_x{count}",
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
    help="Output PNG path",
    required=True,
)
@click.option("--max-reads", "-n", type=int, default=20000, help="Max reads to sample")
@click.option("--title", type=str, default=None, help="Optional plot title")
@click.option("--show-segments/--no-show-segments", default=False, help="Overlay HMM segments")
@click.option("--show-freq/--no-show-freq", default=True, help="Overlay A/C/G/T frequencies")
@click.option("--show-posteriors/--no-show-posteriors", default=False, help="Overlay per-state posterior ribbons")
def plot_hmm(input_file: Path, format: str, output: Path, max_reads: int, title: str, show_segments: bool, show_freq: bool, show_posteriors: bool):
    """Plot per-position entropy with HMM segment overlays.

    Example:
        getRPF plot-hmm input.fastq -f fastq -o hmm.png -n 20000
    """
    try:
        try:
            from .viz.hmm_plot import plot_hmm_entropy
        except ModuleNotFoundError as e:
            raise click.ClickException(
                "plot-hmm requires optional visualization dependencies.\n"
                "Install with: pip install getRPF[viz]  (or)  conda install matplotlib seaborn"
            ) from e
        out = plot_hmm_entropy(
            input_file, format, output,
            max_reads=max_reads, title=title,
            show_segments=show_segments, show_freq=show_freq, show_posteriors=show_posteriors
        )
        click.echo(f"✅ Wrote HMM entropy plot: {out}")
    except Exception as e:
        click.echo(f"❌ Plot failed: {e}", err=True)
        raise click.Abort()


@cli.command()
@click.argument("sample_id")
@click.argument("db_path", type=click.Path(path_type=Path))
@click.option("--source", type=click.Path(path_type=Path), required=True, help="Source FASTQ/FASTA path")
@click.option("--clean-report", type=click.Path(exists=True, path_type=Path), help="check-cleanliness report .txt")
@click.option("--align-json", type=click.Path(exists=True, path_type=Path), help="align-detect JSON report")
@click.option("--extract-json", type=click.Path(exists=True, path_type=Path), help="extract/extract_rpf JSON")
def ingest_duckdb(sample_id: str, db_path: Path, source: Path, clean_report: Path, align_json: Path, extract_json: Path):
    """Ingest getRPF/STAR outputs into a DuckDB for QC review.

    Example:
        getRPF ingest-duckdb SAMPLE1 qc.db \
            --source SRR.fastq.gz \
            --clean-report reports/SRR_cleanliness_report.txt \
            --align-json SRR_align.json \
            --extract-json SRR.extraction_report.json
    """
    try:
        try:
            from .duckdb.ingest import ingest_all
        except ModuleNotFoundError as e:
            raise click.ClickException(
                "ingest-duckdb requires the 'duckdb' package.\n"
                "Install with: pip install duckdb  (or) conda install duckdb"
            ) from e
        out_db = ingest_all(
            db_path=db_path,
            sample_id=sample_id,
            source_path=source,
            cleanliness_report=clean_report,
            alignment_json=align_json,
            extraction_json=extract_json,
        )
        click.echo(f"✅ Ingested into DuckDB: {out_db}")
        click.echo("   Tables: samples, getrpf_checks, alignment_stats, extraction_summary; view: qc_overview")
    except Exception as e:
        click.echo(f"❌ Ingest failed: {e}", err=True)
        raise click.Abort()


@cli.command(name="plot-softclips")
@click.option("--align-json", type=click.Path(exists=True, path_type=Path), required=True, help="align-detect JSON report")
@click.option("--bam", type=click.Path(exists=True, path_type=Path), help="Aligned BAM (optional for heatmap)")
@click.option("--no-heatmap", is_flag=True, help="Disable heatmap even if BAM provided")
@click.option("--output", "-o", type=click.Path(path_type=Path), required=True, help="Output PNG path")
@click.option("--title", type=str, default=None, help="Optional figure title")
def plot_softclips(align_json: Path, bam: Path, no_heatmap: bool, output: Path, title: str):
    """Plot alignment soft-clipping summary and optional heatmap.

    Examples:
      getRPF plot-softclips --align-json SRR_align.json -o softclips.png
      getRPF plot-softclips --align-json SRR_align.json --bam subset.bam -o softclips_heat.png
    """
    try:
        try:
            from .viz.softclip_plot import plot_softclips as render_softclips
        except ModuleNotFoundError as e:
            raise click.ClickException(
                "plot-softclips requires optional visualization dependencies.\n"
                "Install with: pip install getRPF[viz]  (or)  conda install matplotlib seaborn"
            ) from e
        out = render_softclips(align_json, output, bam, not no_heatmap, title)
        click.echo(f"✅ Wrote soft-clipping plot: {out}")
    except Exception as e:
        click.echo(f"❌ Soft-clipping plot failed: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    cli()
