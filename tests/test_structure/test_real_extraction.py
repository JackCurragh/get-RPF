"""Real-library checks of the transform (P5) on held-out reads.

Written before these checks were run. Held-out sets come from the same
1M-read subsets as test_real_panel.py: reads 300k-600k for architecture
stability, and reads 600k-1M for extraction and Q5. None of them overlaps the
inference sample (the first 300k).
"""

import itertools
from functools import lru_cache

import pytest

from getRPF.core.structure.assemble import infer_structure
from getRPF.core.structure.transform import (
    extract_reads,
    iter_fastq,
    validate_extraction,
)

from .test_real_panel import DATA, run

ALL = ["SRR1944950", "SRR3945920", "SRR3945930", "SRR12693498", "SRR23242345"]
EMITTED = ["SRR1944950", "SRR3945920", "SRR3945930", "SRR12693498"]


def _path(accession):
    path = DATA / f"{accession}.fastq.gz"
    if not path.exists():
        pytest.skip(f"{path} not present; set GETRPF_PANEL_DATA")
    return path


@lru_cache(maxsize=None)
def held_out_inference(accession):
    records = itertools.islice(iter_fastq(_path(accession)), 300_000, 600_000)
    headers, reads = [], []
    for header, sequence, _ in records:
        headers.append(header)
        reads.append(sequence)
    return infer_structure(reads, headers)


@lru_cache(maxsize=None)
def held_out_extraction(accession):
    architecture = run(accession).architecture
    return extract_reads(
        _path(accession), None, architecture, skip=600_000, limit=400_000
    )


@pytest.mark.parametrize("accession", ALL)
def test_architecture_is_stable_on_held_out_reads(accession):
    first, second = run(accession), held_out_inference(accession)
    assert second.architecture.describe() == first.architecture.describe()
    assert second.transform.emit is first.transform.emit


@pytest.mark.parametrize("accession", EMITTED)
def test_emitted_inserts_look_like_riboseq(accession):
    summary = held_out_extraction(accession)
    assert summary.input_reads >= 100_000
    assert validate_extraction(summary).value == "riboseq_likely"


def test_audit_and_production_agree_on_real_reads(tmp_path):
    path = _path("SRR3945920")
    architecture = run("SRR3945920").architecture
    audit = extract_reads(path, None, architecture, limit=50_000)
    production = extract_reads(
        path, tmp_path / "out.fastq.gz", architecture, limit=50_000
    )
    assert audit.to_dict() == production.to_dict()
