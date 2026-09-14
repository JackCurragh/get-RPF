"""`getRPF extract --infer-structure`: one command from FASTQ to RPFs."""

import json

from click.testing import CliRunner

from getRPF.cli import cli
from getRPF.core.structure.workflow import output_prefix

from .sim import nta, random_bases, simulate
from .test_cli_structure import write_fastq


def run(directory, reads, *args):
    directory.mkdir(parents=True, exist_ok=True)
    fastq = directory / "in.fastq"
    write_fastq(fastq, reads)
    out = directory / "out" / "s_trimmed.fastq"
    result = CliRunner().invoke(
        cli, ["extract", str(fastq), str(out), "--infer-structure", *args]
    )
    return result, out.parent


def report(out):
    return json.loads((out / "s_trimmed.extraction_report.json").read_text())


def test_output_prefix_matches_the_legacy_naming(tmp_path):
    assert output_prefix(tmp_path / "x_trimmed.fastq").name == "x_trimmed"
    assert output_prefix(tmp_path / "x.rpf.fastq.gz").name == "x.rpf"
    assert output_prefix(tmp_path / "x.fq.gz").name == "x"


def test_collapsed_output_counts_every_accepted_insert(tmp_path):
    library = simulate(five_prime=[random_bases(5)], n_reads=5_000)
    result, out = run(tmp_path, library.reads, "--collapsed-only")
    assert result.exit_code == 0, result.output
    extraction = report(out)
    assert extraction["status"] == "emitted"
    assert extraction["architecture"].startswith("[random 5 UMI][insert]")
    headers = [
        line
        for line in (out / "s_trimmed.collapsed.fa").read_text().splitlines()
        if line.startswith(">")
    ]
    counts = [int(header.rsplit("_x", 1)[1]) for header in headers]
    assert sum(counts) == extraction["extraction"]["accepted"] > 0
    for name in ("structure.json", "structure.txt", "seqspec.yaml"):
        assert (out / f"s_trimmed.{name}").exists()
    assert not (out / "s_trimmed.fastq").exists()


def test_fastq_output_keeps_names_qualities_and_umis(tmp_path):
    library = simulate(five_prime=[random_bases(5)], n_reads=5_000)
    result, out = run(tmp_path, library.reads)
    assert result.exit_code == 0, result.output
    name, sequence, _, quality = (out / "s_trimmed.fastq").read_text().splitlines()[:4]
    assert name.startswith("@r") and len(name.split("_")[-1]) == 5
    assert len(sequence) == len(quality)


def test_common_junction_bases_are_kept_and_flagged(tmp_path):
    skewed = {"T": 0.7, "C": 0.1, "A": 0.1, "G": 0.1}
    library = simulate(five_prime=[nta({2: 1.0}, skewed)])
    result, out = run(tmp_path, library.reads, "--collapsed-only")
    assert result.exit_code == 0, result.output
    extraction = report(out)
    assert extraction["status"] == "emitted"
    assert any(flag.startswith("Q2:") for flag in extraction["transform"]["flags"])
    assert "flag: Q2:" in (out / "s_trimmed.structure.txt").read_text()
    review, _ = run(tmp_path / "review", library.reads, "--fail-on", "review")
    assert review.exit_code == 3
    hold, _ = run(tmp_path / "hold", library.reads, "--fail-on", "hold")
    assert hold.exit_code == 0


def test_withheld_structure_writes_reports_and_no_reads(tmp_path):
    library = simulate(n_reads=20_000, n_transcripts=5_000, zipf_s=0.0)
    result, out = run(tmp_path, library.reads, "--collapsed-only")
    assert result.exit_code == 0, result.output
    extraction = report(out)
    assert extraction["status"] == "withheld"
    assert extraction["transform"]["reasons"]
    assert (out / "s_trimmed.structure.txt").exists()
    assert not (out / "s_trimmed.collapsed.fa").exists()
    held, _ = run(
        tmp_path / "hold", library.reads, "--collapsed-only", "--fail-on", "hold"
    )
    assert held.exit_code == 3
