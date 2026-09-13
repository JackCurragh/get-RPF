"""`getRPF infer-structure` end to end on small synthetic FASTQs."""

from click.testing import CliRunner

from getRPF.cli import cli

from .sim import random_bases, simulate


def write_fastq(path, reads):
    with open(path, "w") as handle:
        for i, read in enumerate(reads):
            handle.write(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n")


def test_infer_structure_writes_reports_and_exits_zero(tmp_path):
    fastq = tmp_path / "sim.fastq"
    write_fastq(fastq, simulate(five_prime=[random_bases(5)]).reads)
    result = CliRunner().invoke(
        cli, ["infer-structure", str(fastq), "-o", str(tmp_path / "out")]
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "sim.structure.json").exists()
    assert (tmp_path / "out" / "sim.structure.txt").exists()
    assert "[random 5 UMI][insert]" in result.output


def test_withheld_transform_exits_three(tmp_path):
    fastq = tmp_path / "shallow.fastq"
    library = simulate(n_reads=20_000, n_transcripts=5_000, zipf_s=0.0)
    write_fastq(fastq, library.reads)
    result = CliRunner().invoke(
        cli, ["infer-structure", str(fastq), "-o", str(tmp_path / "out")]
    )
    assert result.exit_code == 3, result.output
    assert "Transform: withheld" in result.output
