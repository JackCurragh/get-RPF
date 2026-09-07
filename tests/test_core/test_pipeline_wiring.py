"""Regression tests for the pipeline wiring fixes found in code review:

1. identity_screen must run on the count-weighted read population, not the
   deduplicated collapsed-FASTA sequence set (duplication was structurally
   always 0.0; length_shape/insert_mode were computed on unique sequences).
2. quality_evidence must actually be computed and passed to build_evidence,
   so basecaller_artifact_hold is reachable.
3. b_adapter must actually be able to fire: known_adapters (from the
   architecture DB) and per_length_reads (from the sketch) must be threaded
   into boundary estimation.
4. Rule records must have a consistent schema regardless of whether they
   came from boundary estimation or an explicit override.
5. boundary_estimates_3p must be populated in the evidence object.
"""

import random

from getRPF.core.apply import RULE_RECORD_KEYS
from getRPF.core.pipeline import run_sample

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def _write_fastq(path, sequences, qualities=None):
    lines = []
    for i, seq in enumerate(sequences):
        quals = (
            "".join(chr(33 + q) for q in qualities[i]) if qualities else "I" * len(seq)
        )
        lines.append(f"@r{i}\n{seq}\n+\n{quals}\n")
    path.write_text("".join(lines))


def test_duplication_and_length_shape_use_weighted_counts_not_unique_sequences(
    tmp_path,
):
    # 900 raw reads that are all the *same* 30nt sequence (heavy true
    # duplication), plus 100 raw reads that are each a *distinct* 40nt
    # sequence (no duplication). Unique-sequence counting would see the
    # length distribution as 1x30 vs 100x40 (mode=40); read-count-weighted
    # counting correctly sees 900x30 vs 100x40 (mode=30). Uses a random
    # (not tandem-periodic) 30-mer so the pre-existing architecture
    # detector doesn't mistake it for a homopolymer/adapter artifact --
    # that's a separate, real detection path this test isn't about.
    rng = random.Random(42)
    duplicated_seq = "".join(rng.choice(BASES) for _ in range(30))
    sequences = [duplicated_seq] * 900
    sequences += [_cycled_seq(i, 40) for i in range(100)]

    input_file = tmp_path / "dup.fastq"
    _write_fastq(input_file, sequences)

    # run_sample builds its own RPFExtractor internally; an empty
    # architecture-db JSON keeps its own detection from interfering with
    # either class, isolating the identity-screen wiring under test.
    empty_db = tmp_path / "empty_architectures.json"
    empty_db.write_text('{"architectures": []}')

    evidence, _result = run_sample(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        format="fastq",
        infer_reads=10_000,
        collapsed_only=True,
        architecture_db=empty_db,
    )

    screen = evidence["biological_screen"]
    assert screen["insert_mode"] == 30  # would be 40 under the old bug
    assert screen["duplication_umi_aware"] > 0.5  # would be 0.0 under the old bug


def test_basecaller_artifact_hold_is_reachable(tmp_path):
    sequences = []
    qualities = []
    for i in range(500):
        body = _cycled_seq(i, 27)
        sequences.append(body + "G")  # constant terminal G, adapter-cliff-like
        qualities.append([38] * 27 + [3])  # quality collapses at the terminal G

    input_file = tmp_path / "polyg.fastq"
    _write_fastq(input_file, sequences, qualities)

    empty_db = tmp_path / "empty_architectures.json"
    empty_db.write_text('{"architectures": []}')

    evidence, _result = run_sample(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        format="fastq",
        infer_reads=10_000,
        collapsed_only=True,
        architecture_db=empty_db,
    )

    assert evidence["quality_evidence"]["three_prime_q_collapse"] is True
    assert evidence["release_class"] == "basecaller_artifact_hold"


def test_b_adapter_fires_using_the_architecture_database(tmp_path):
    # CTGTAGGCACCATCAAT is a real adapter from the built-in
    # ingolia_2009.yaml architecture. Embed its first 10 bases as a 3'
    # terminal artifact so b_adapter (which matches known adapter
    # prefixes against read suffixes) can find it.
    adapter_prefix = "CTGTAGGCAC"
    sequences = [_cycled_seq(i, 30) + adapter_prefix for i in range(300)]

    input_file = tmp_path / "adapter.fastq"
    _write_fastq(input_file, sequences)

    evidence, _result = run_sample(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        format="fastq",
        infer_reads=10_000,
        collapsed_only=True,
        # architecture_db=None -> built-in architectures load, including
        # ingolia_2009, giving b_adapter a real known-adapter list.
    )

    length = 30 + len(adapter_prefix)
    decision = evidence["boundary_estimates_3p"].get(str(length))
    assert decision is not None
    assert "b_adapter" in decision["estimates"]


def test_rule_record_schema_is_consistent_for_overrides(tmp_path):
    sequences = [_cycled_seq(i, 29) for i in range(50)]
    input_file = tmp_path / "sample.fastq"
    _write_fastq(input_file, sequences)

    empty_db = tmp_path / "empty_architectures.json"
    empty_db.write_text('{"architectures": []}')

    evidence, _result = run_sample(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        format="fastq",
        infer_reads=10_000,
        collapsed_only=True,
        architecture_db=empty_db,
        override_trims={29: {"trim_5p": 0, "trim_3p": 1}},
    )

    assert evidence["trim_rules_applied"]
    for record in evidence["trim_rules_applied"]:
        assert set(record.keys()) == set(RULE_RECORD_KEYS)


def test_boundary_estimates_3p_is_populated(tmp_path):
    sequences = []
    for i in range(20_000):
        body = _cycled_seq(i, 29)
        tail = "A" if i < 19_800 else BASES[i % 4]
        sequences.append(body + tail)

    input_file = tmp_path / "artifact.fastq"
    _write_fastq(input_file, sequences)

    empty_db = tmp_path / "empty_architectures.json"
    empty_db.write_text('{"architectures": []}')

    evidence, _result = run_sample(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        format="fastq",
        infer_reads=100_000,
        collapsed_only=True,
        architecture_db=empty_db,
    )

    assert evidence["boundary_estimates_3p"]
    decision = evidence["boundary_estimates_3p"]["30"]
    assert decision["decision"] == "safe_to_trim"
    assert decision["trim_bases"] == 1
