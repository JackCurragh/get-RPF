"""Real-library validation of steps 0-3 against validation/panel_v1/truth.yaml.

These expectations were written from the gold set before the runs were
analysed. Where the gold set's expected structure is marked "derived" or
"unconfirmed", a failure here is a finding about the claim, not a licence to
change the test: update truth.yaml with the measurement and say why.

Skipped unless GETRPF_PANEL_DATA points at a directory of <run>.fastq.gz read
subsets (scripts/fetch_panel_subsets.sh).
"""

import os
from functools import lru_cache
from pathlib import Path

import pytest

from getRPF.core.structure import Status
from getRPF.core.structure.assemble import infer_structure
from getRPF.core.structure.observe import read_fastq

DATA = Path(os.environ.get("GETRPF_PANEL_DATA", "/nonexistent"))


@lru_cache(maxsize=None)
def run(accession):
    path = DATA / f"{accession}.fastq.gz"
    if not path.exists():
        pytest.skip(f"{path} not present; set GETRPF_PANEL_DATA")
    headers, reads, _ = read_fastq(path, 300_000)
    return infer_structure(reads, headers)


def test_srr1944950_trimmed_insert_with_candidate_5p_nta():
    result = run("SRR1944950")
    assert result.observation.input_state == "trimmed"
    assert result.q1.status is Status.NOT_OBSERVABLE
    assert result.q2.status is Status.INTERVAL
    assert result.q2.value.technical_length == 0
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.technical_length == 0


def test_srr3945920_truseq_adapter_directly_after_the_footprint():
    result = run("SRR3945920")
    assert result.observation.input_state == "raw_fixed_length"
    assert result.q1.value.adapter.name == "truseq_3p"
    assert result.q1.value.support >= 0.9
    assert 27 <= result.q1.value.start_mode <= 32
    assert result.q2.status is not Status.UNDERPOWERED
    assert result.q3.status is not Status.UNDERPOWERED


def test_srr3945930_four_plus_three_randomised_bases():
    result = run("SRR3945930")
    assert result.q1.status in (Status.RESOLVED, Status.INTERVAL)
    assert result.q2.value.blocks == (("random", 3),)
    assert result.q3.value.blocks == (("random", 4),)


def test_srr12693498_dplex_umi_spacer_and_tail():
    result = run("SRR12693498")
    assert result.q1.value.adapter.role == "insert_3p"
    assert result.q1.value.tail_base == "A"
    assert result.q2.value.technical_length == 16
    kind, length = result.q2.value.blocks[0]
    assert kind == "random" and length >= 12
    assert result.q3.value.technical_length == 0


def test_srr23242345_adapter_mid_read_and_never_p7():
    result = run("SRR23242345")
    assert result.q1.value.adapter.name == "truseq_3p"
    assert result.q1.value.downstream_consistency >= 0.8


@pytest.mark.xfail(
    strict=True,
    reason="gold claim contradicted by the data (truth.yaml, 2026-09-13): read "
    "positions 6-7 are skewed (60% T, 78% C), not random, and no technical "
    "base precedes the adapter",
)
def test_srr23242345_claimed_seven_nt_umi_and_last_base():
    result = run("SRR23242345")
    assert result.q2.value.blocks == (("random", 7),)
    q3 = result.q3.value
    assert q3.technical_length + q3.nta_length == 1


@pytest.mark.parametrize(
    "accession, emit, architecture",
    [
        ("SRR1944950", True, "[nta 0-1 kept][insert]"),
        (
            "SRR3945920",
            True,
            "[insert][nta 0-1 kept][adapter AGATCGGAAGAG...]",
        ),
        (
            "SRR3945930",
            True,
            "[random 3 UMI][nta 0-1 kept][insert][random 4 UMI][adapter CTGTAGGCACCA...]",
        ),
        (
            "SRR12693498",
            True,
            "[random 13 UMI][fixed GGG][insert][nta 0-1 kept][poly(A)]"
            "[adapter AGATCGGAAGAG...]",
        ),
        (
            "SRR23242345",
            False,
            "[random 5 UMI][nta 0-2 kept][insert][adapter AGATCGGAAGAG...]",
        ),
    ],
)
def test_measured_architecture_and_transform_decision(accession, emit, architecture):
    """Regression guard on the measured architectures (2026-09-13), not a
    validation claim. SRR23242345 is withheld because read positions 6-7 look
    non-templated in most reads (spec §7.2)."""
    result = run(accession)
    assert result.architecture.describe() == architecture
    assert result.transform.emit is emit


def test_srr23242345_measured_structure():
    """Regression guard on the measured structure, not a validation claim."""
    result = run("SRR23242345")
    q2 = result.q2.value
    assert result.q2.status is Status.INTERVAL
    assert q2.blocks == (("random", 5),)
    assert q2.technical_length + q2.nta_length == 7
    assert result.q3.value.technical_length == 0
