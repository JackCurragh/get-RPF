"""P3 gate on real data: SRR1944950 must reproduce the spec §5.3 table.

Skipped unless GETRPF_SRR1944950_FASTQ points at the run's FASTQ (not in the
repo; ENA SRR1944950, md5 recorded in the validation run manifests).
"""

import gzip
import os
import random
from pathlib import Path

import pytest

from getRPF.core.structure import Status, infer_junctions

FASTQ = os.environ.get("GETRPF_SRR1944950_FASTQ", "")

pytestmark = pytest.mark.skipif(
    not FASTQ or not Path(FASTQ).exists(),
    reason="set GETRPF_SRR1944950_FASTQ to run the real-data gate",
)


@pytest.fixture(scope="module")
def reads():
    out = []
    with gzip.open(FASTQ, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                out.append(line.strip())
                if len(out) == 300_000:
                    break
    return out


def agreements(profile, n=10):
    return [stats.agreement for stats in profile.positions[:n]]


def top_base(call):
    return max(call.nta_bases, key=call.nta_bases.get)


def test_real_reads(reads):
    result = infer_junctions(reads)
    start = agreements(result.first_pass.read_start)
    end = agreements(result.first_pass.anchor)
    assert 0.58 <= start[0] <= 0.70
    assert all(a >= 0.90 for a in start[1:])
    assert all(a >= 0.90 for a in end)
    assert result.q2.status is Status.INTERVAL
    assert (result.q2.value.technical_length, result.q2.value.nta_length) == (0, 1)
    assert top_base(result.q2.value) == "T"
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.technical_length == 0


def test_spiked_five_prime_umi(reads):
    rng = random.Random(1)
    result = infer_junctions(["".join(rng.choices("ACGT", k=5)) + r for r in reads])
    assert all(a <= 0.32 for a in agreements(result.first_pass.read_start, 5))
    assert result.q2.value.technical_length == 5
    assert result.q2.value.blocks == (("random", 5),)
    # The real non-templated first base now sits just after the UMI.
    assert result.q2.value.nta_length <= 1


def test_spiked_three_prime_random(reads):
    rng = random.Random(2)
    result = infer_junctions([r + "".join(rng.choices("ACGT", k=4)) for r in reads])
    assert all(a <= 0.32 for a in agreements(result.first_pass.anchor, 4))
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.blocks == (("random", 4),)


def test_spiked_untemplated_g(reads):
    rng = random.Random(3)
    result = infer_junctions(["G" + r if rng.random() < 0.7 else r for r in reads])
    assert result.q2.status is Status.INTERVAL
    assert result.q2.value.technical_length == 0
    assert 1 <= result.q2.value.nta_length <= 2
    assert top_base(result.q2.value) == "G"
