"""P7: the yeast panel runs checked against the reference genome.

Written before any alignment was run. Skipped unless GETRPF_YEAST_STAR_INDEX
points at a STAR index of S. cerevisiae R64-1-1 (scripts/fetch_yeast_reference.sh)
and the panel read subsets are present (GETRPF_PANEL_DATA).

Expectations:
- for all three yeast runs, the emitted inserts show no missed technical block
  (clips of >=2 nt at either end in under 10% of aligned inserts);
- SRR1944950's 5' junction base, found by the pileup, is confirmed by the
  genome (spec §10.2: "confirm the SRR1944950 NTA by alignment").
"""

import os
from functools import lru_cache
from pathlib import Path

import pytest

from getRPF.core.structure.align import apply_alignment_check

from .test_real_panel import DATA, run

INDEX = Path(os.environ.get("GETRPF_YEAST_STAR_INDEX", "/nonexistent"))
YEAST = ["SRR1944950", "SRR3945920", "SRR3945930"]


@lru_cache(maxsize=None)
def checked(accession):
    if not (INDEX / "SA").exists():
        pytest.skip("set GETRPF_YEAST_STAR_INDEX to a yeast STAR index")
    from getRPF.core.structure.observe import read_fastq

    path = DATA / f"{accession}.fastq.gz"
    if not path.exists():
        pytest.skip(f"{path} not present")
    _, reads, _ = read_fastq(path, 300_000)
    return apply_alignment_check(run(accession), reads, INDEX)


@pytest.mark.parametrize("accession", YEAST)
def test_no_missed_technical_block(accession):
    _, (q2, q3) = checked(accession)
    # Clip fractions are 0 when nothing aligns: require a real verdict first.
    assert q2.verdict != "underpowered", f"only {q2.aligned} inserts aligned"
    assert q2.longer_clip < 0.10
    assert q3.longer_clip < 0.10


def test_srr1944950_five_prime_junction_base_is_confirmed():
    _, (q2, _) = checked("SRR1944950")
    assert q2.verdict == "confirmed"
