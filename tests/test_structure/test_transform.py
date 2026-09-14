"""The per-read transform, extraction, seqspec round trip and templates."""

import gzip

import yaml
from click.testing import CliRunner

from getRPF.cli import cli
from getRPF.core.structure.assemble import infer_structure
from getRPF.core.structure.model import Architecture, Block, Status
from getRPF.core.structure.seqspec_io import (
    compare_catalogue,
    from_seqspec,
    load_catalogue,
    to_seqspec,
)
from getRPF.core.structure.transform import (
    TRANSFORM_SCHEMA,
    Accepted,
    Transform,
    extract_reads,
)

from .sim import fixed, nta, random_bases, simulate
from .test_cli_structure import write_fastq


def test_transform_recovers_true_inserts_and_umis():
    library = simulate(five_prime=[random_bases(5)], three_prime=[random_bases(4)])
    transform = Transform(infer_structure(library.reads).architecture)
    accepted = correct = 0
    for read, (start, end) in zip(library.reads, library.inserts):
        outcome = transform.apply(read, "I" * len(read))
        if isinstance(outcome, Accepted):
            accepted += 1
            umi = read[:start] + read[end : end + 4]
            correct += outcome.insert == read[start:end] and outcome.umi == umi
    assert accepted / len(library.reads) >= 0.95
    assert correct / accepted >= 0.99


def test_tail_library_cuts_at_the_run_start():
    tail = nta({n: 1.0 for n in range(8, 25)}, {"A": 1.0})
    library = simulate(
        five_prime=[random_bases(13), fixed("GGG")],
        three_prime=[tail],
        read_length=100,
        n_reads=20_000,
    )
    result = infer_structure(library.reads)
    described = result.architecture.describe()
    assert described.startswith("[random 13 UMI][fixed GGG][insert]")
    assert described.endswith("[poly(A)][adapter AGATCGGAAGAG...]")
    transform = Transform(result.architecture)
    exact = accepted = 0
    for read, (start, end) in zip(library.reads, library.inserts):
        outcome = transform.apply(read, "I" * len(read))
        if isinstance(outcome, Accepted):
            accepted += 1
            exact += outcome.insert == read[start:end]
    assert accepted / len(library.reads) >= 0.9
    assert exact / accepted >= 0.7  # inserts ending in A lose those bases


def test_tail_architecture_cuts_at_the_run_start_even_for_one_base():
    """Spec §6: the transform uses the same tail rule as inference."""
    architecture = Architecture(
        (
            Block(
                "insert", "read_start", (20, 40), None, False, False, Status.INTERVAL
            ),
            Block("tail", "anchor", (0, 30), "A", False, True, Status.INTERVAL),
            Block(
                "adapter",
                "anchor",
                (34, 34),
                "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                False,
                True,
                Status.RESOLVED,
            ),
        ),
        "inferred",
        Status.INTERVAL,
        "monosome_20_40",
    )
    read = "C" * 29 + "A" + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

    outcome = Transform(architecture).apply(read, "I" * len(read))

    assert isinstance(outcome, Accepted)
    assert outcome.insert == "C" * 29
    assert outcome.boundary == "tail_A"


def test_audit_retains_candidate_lengths_outside_the_fragment_policy(tmp_path):
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    architecture = Architecture(
        (
            Block(
                "insert", "read_start", (20, 40), None, False, False, Status.RESOLVED
            ),
            Block("adapter", "anchor", (34, 34), adapter, False, True, Status.RESOLVED),
        ),
        "inferred",
        Status.RESOLVED,
        "monosome_20_40",
    )
    fastq = tmp_path / "lengths.fastq"
    write_fastq(fastq, ["C" * 19 + adapter, "C" * 30 + adapter, "C" * 41 + adapter])

    summary = extract_reads(fastq, None, architecture)

    assert summary.accepted == 1
    assert summary.rejected == {"longer_than_policy": 1, "shorter_than_policy": 1}
    assert summary.to_dict()["candidate_insert_lengths"] == {"19": 1, "30": 1, "41": 1}
    assert summary.to_dict()["boundaries"] == {"adapter": 3}
    assert summary.to_dict()["transform_schema"] == TRANSFORM_SCHEMA


def test_audit_and_production_agree_and_keep_names(tmp_path):
    library = simulate(five_prime=[random_bases(5)], n_reads=5_000)
    fastq = tmp_path / "in.fastq"
    write_fastq(fastq, library.reads)
    architecture = infer_structure(library.reads).architecture
    audit = extract_reads(fastq, None, architecture)
    output = tmp_path / "out.fastq.gz"
    production = extract_reads(fastq, output, architecture)
    assert audit.to_dict() == production.to_dict()
    with gzip.open(output, "rt") as handle:
        header, sequence, _, quality = (handle.readline().strip() for _ in range(4))
    assert header.startswith("@r") and "_" in header
    assert len(sequence) == len(quality)


def test_seqspec_round_trip_is_exact():
    architecture = infer_structure(
        simulate(five_prime=[random_bases(5)], three_prime=[random_bases(4)]).reads
    ).architecture
    loaded = from_seqspec(yaml.safe_load(yaml.safe_dump(to_seqspec(architecture))))
    assert loaded == architecture


def test_catalogue_templates_load_and_dplex_template_is_supported():
    names = {name for name, _ in load_catalogue()}
    assert {"riboflow_template_switch", "mcglincy_ingolia_2017"} <= names
    tail = nta({n: 1.0 for n in range(8, 25)}, {"A": 1.0})
    library = simulate(
        five_prime=[random_bases(13), fixed("GGG")],
        three_prime=[tail],
        read_length=100,
        n_reads=20_000,
    )
    comparisons = {
        c.template: c.verdict
        for c in compare_catalogue(infer_structure(library.reads).architecture)
    }
    assert comparisons["riboflow_template_switch"] == "supported"
    assert comparisons["mcglincy_ingolia_2017"] == "rejected"


def test_extract_refuses_a_withheld_architecture(tmp_path):
    # Too few informative groups: Q2/Q3 underpowered, so the transform is
    # withheld (junction bases alone no longer withhold, spec §7.2).
    library = simulate(n_reads=20_000, n_transcripts=5_000, zipf_s=0.0)
    fastq = tmp_path / "shallow.fastq"
    write_fastq(fastq, library.reads)
    runner = CliRunner()
    inferred = runner.invoke(cli, ["infer-structure", str(fastq), "-o", str(tmp_path)])
    assert inferred.exit_code == 3
    seqspec = tmp_path / "shallow.seqspec.yaml"
    if not seqspec.exists():
        return  # no architecture at all: there is nothing extract could apply
    extracted = runner.invoke(
        cli,
        [
            "extract",
            str(fastq),
            str(tmp_path / "out.fastq.gz"),
            "--architecture",
            str(seqspec),
        ],
    )
    assert extracted.exit_code == 3


def test_infer_then_extract_end_to_end(tmp_path):
    fastq = tmp_path / "umi.fastq"
    write_fastq(fastq, simulate(five_prime=[random_bases(5)], n_reads=5_000).reads)
    runner = CliRunner()
    assert (
        runner.invoke(
            cli, ["infer-structure", str(fastq), "-o", str(tmp_path)]
        ).exit_code
        == 0
    )
    output = tmp_path / "umi.rpf.fastq.gz"
    result = runner.invoke(
        cli,
        [
            "extract",
            str(fastq),
            str(output),
            "--architecture",
            str(tmp_path / "umi.seqspec.yaml"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert output.exists()
    assert (tmp_path / "umi.rpf.fastq.gz.summary.json").exists()
