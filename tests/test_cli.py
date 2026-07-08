"""Test suite for getRPF CLI functionality."""

import gzip

import pytest
from click.testing import CliRunner
from getRPF.cli import cli


@pytest.fixture
def test_data(tmp_path):
    """Create test FASTQ data in a temporary directory."""
    test_file = tmp_path / "test_data.fastq"
    reads = [
        # Normal read
        "@read1\n" "ATCGATCGATCGATCG\n" "+\n" "FFFFFFFFFFFFFFFA\n",
        # Read with adapter
        "@read2\n" "ATCGATCGAGATCGGAAGAG\n" "+\n" "FFFFFFFFFFFFFFFFF###\n",
        # Short read
        "@read3\n" "ATCGATCG\n" "+\n" "FFFFFFFF\n",
        # Read with low quality
        "@read4\n" "ATCGATCGATCGATCG\n" "+\n" "FF##FFFFFFFFFF##\n",
    ]

    test_file.write_text("".join(reads))
    return test_file


@pytest.fixture
def runner():
    """Create a CLI test runner."""
    return CliRunner()


def test_version(runner):
    """Test CLI version command."""
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "0.2.4" in result.output


def test_help(runner):
    """Test CLI help command."""
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "getRPF" in result.output
    assert "check" in result.output
    assert "detect-adapter" in result.output


def test_extract_help_accepts_pipeline_flags(runner):
    """The pipeline calls getRPF extract with these options."""
    result = runner.invoke(cli, ["extract", "--help"])

    assert result.exit_code == 0
    assert "--format" in result.output
    assert "--generate-seqspec" in result.output
    assert "--output-format" in result.output
    assert "--star-index" in result.output
    assert "--star-threads" in result.output
    assert "--collapsed-only" in result.output


def test_extract_rpf_is_extract_alias(runner):
    """extract-rpf is a compatibility alias, not a second implementation."""
    assert cli.commands["extract"].callback is cli.commands["extract-rpf"].callback

    result = runner.invoke(cli, ["extract-rpf", "--help"])

    assert result.exit_code == 0
    assert "--format" in result.output
    assert "--generate-seqspec" in result.output
    assert "--output-format" in result.output
    assert "--star-index" in result.output
    assert "--star-threads" in result.output
    assert "--collapsed-only" in result.output


class TestQualityCheck:
    """Test suite for quality check command."""

    def test_basic_check(self, runner, test_data, tmp_path):
        """Test basic quality check functionality."""
        output_file = tmp_path / "quality_report.txt"
        result = runner.invoke(
            cli,
            [
                "check",
                str(test_data),
                "--format",
                "fastq",
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code == 0
        assert output_file.exists()

        # Check report content
        content = output_file.read_text()
        assert "Read Length Distribution" in content
        assert "Nucleotide Frequencies" in content
        assert "Quality Scores" in content

    def test_invalid_format(self, runner, test_data, tmp_path):
        """Test quality check with invalid format."""
        output_file = tmp_path / "quality_report.txt"
        result = runner.invoke(
            cli,
            [
                "check",
                str(test_data),
                "--format",
                "invalid",
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code != 0
        assert "Invalid value for '--format'" in result.output

    def test_missing_file(self, runner, tmp_path):
        """Test quality check with non-existent file."""
        output_file = tmp_path / "quality_report.txt"
        result = runner.invoke(
            cli,
            [
                "check",
                "nonexistent.fastq",
                "--format",
                "fastq",
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code != 0
        assert "does not exist" in result.output


class TestAdapterDetection:
    """Test suite for adapter detection command."""

    def test_basic_detection(self, runner, test_data, tmp_path):
        """Test basic adapter detection functionality."""
        output_file = tmp_path / "adapter_report.txt"
        result = runner.invoke(
            cli,
            [
                "detect-adapter",
                str(test_data),
                "--format",
                "fastq",
                "--adapter",
                "AGATCGGAAGAG",
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code == 0
        assert output_file.exists()

        # Check report content
        content = output_file.read_text()
        assert "Adapter Contamination Summary" in content
        assert "contamination rate" in content.lower()

    def test_gzipped_fastq_detection(self, runner, test_data, tmp_path):
        """Test adapter detection on gzipped FASTQ input."""
        gzipped_input = tmp_path / "test_data.fastq.gz"
        gzipped_input.write_bytes(gzip.compress(test_data.read_bytes()))

        output_file = tmp_path / "adapter_report.txt"
        result = runner.invoke(
            cli,
            [
                "detect-adapter",
                str(gzipped_input),
                "--format",
                "fastq",
                "--adapter",
                "AGATCGGAAGAG",
                "--output",
                str(output_file),
                "--max-reads",
                "2",
            ],
        )

        assert result.exit_code == 0
        content = output_file.read_text()
        assert "Overall contamination rate: 50.00%" in content

    def test_invalid_adapter(self, runner, test_data, tmp_path):
        """Test adapter detection with invalid adapter sequence."""
        output_file = tmp_path / "adapter_report.txt"
        result = runner.invoke(
            cli,
            [
                "detect-adapter",
                str(test_data),
                "--format",
                "fastq",
                "--adapter",
                "X",  # Invalid base
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code != 0
        assert "Invalid adapter sequence" in str(result.exception)

    @pytest.mark.parametrize(
        "adapter",
        [
            "AGATCGGAAGAG",  # Common Illumina adapter
            "CTGTCTCTTATACACATCT",  # Another common adapter
            "ACGT",  # Short adapter
            "A" * 30,  # Long adapter
        ],
    )
    def test_various_adapters(self, runner, test_data, tmp_path, adapter):
        """Test adapter detection with various adapter sequences."""
        output_file = tmp_path / "adapter_report.txt"
        result = runner.invoke(
            cli,
            [
                "detect-adapter",
                str(test_data),
                "--format",
                "fastq",
                "--adapter",
                adapter,
                "--output",
                str(output_file),
            ],
        )
        assert result.exit_code == 0
        assert output_file.exists()


def test_integration(runner, test_data, tmp_path):
    """Test running multiple commands in sequence."""
    # Run quality check
    quality_output = tmp_path / "quality_report.txt"
    result1 = runner.invoke(
        cli,
        ["check", str(test_data), "--format", "fastq", "--output", str(quality_output)],
    )
    assert result1.exit_code == 0

    # Run adapter detection
    adapter_output = tmp_path / "adapter_report.txt"
    result2 = runner.invoke(
        cli,
        [
            "detect-adapter",
            str(test_data),
            "--format",
            "fastq",
            "--adapter",
            "AGATCGGAAGAG",
            "--output",
            str(adapter_output),
        ],
    )
    assert result2.exit_code == 0

    # Verify both outputs exist and have content
    assert quality_output.exists()
    assert adapter_output.exists()
    assert quality_output.stat().st_size > 0
    assert adapter_output.stat().st_size > 0
