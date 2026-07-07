"""Tests for packaged observed protocol architectures."""

import csv
import importlib.resources
from pathlib import Path

from getRPF.core.processors.rpf_extractor import ArchitectureDatabase


def test_observed_protocols_load_as_builtin_architectures():
    db = ArchitectureDatabase()

    observed = [
        arch
        for arch in db.architectures
        if arch.protocol_name.startswith("observed_")
    ]

    assert len(observed) == 18
    assert all(arch.expected_rpf_length == (15, 40) for arch in observed)
    assert all(arch.adapter_sequences for arch in observed)
    assert all("U" not in adapter for arch in observed for adapter in arch.adapter_sequences)
    assert "observed_ctgtaggc_3p_adapter" in {
        arch.protocol_name for arch in observed
    }


def test_observed_dual_ligation_loader_keeps_5p_adapter_out_of_trim_targets():
    db = ArchitectureDatabase()
    dual = next(
        arch
        for arch in db.architectures
        if arch.protocol_name == "observed_dual_ligation_adapter_pair"
    )

    assert dual.rpf_start == len("GTTCAGAGTTCTACAGTCCGACGATC")
    assert dual.rpf_end == -1
    assert dual.adapter_sequences == [
        "GTTCAGAGTTCTACAGTCCGACGATC",
        "TCGTATGCCGTCTTCTGCTTG",
    ]
    assert dual.trim_adapter_sequences == ["TCGTATGCCGTCTTCTGCTTG"]


def test_external_validation_manifest_is_not_packaged_with_getrpf():
    packaged_manifest = (
        importlib.resources.files("getRPF")
        .joinpath("architectures")
        .joinpath("ribobase_prior_manifest.csv")
    )
    assert not packaged_manifest.is_file()

    external_manifest = Path(
        "/Users/jackt/projects/all-RiboSeq/6k_runs/analysis/"
        "ribobase_external_validation_manifest.csv"
    )
    if not external_manifest.exists():
        return

    with external_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 904
    assert sum(1 for row in rows if row["seqspec_id"]) == 898
    assert sum(1 for row in rows if not row["seqspec_id"]) == 6

    srr3306589 = next(row for row in rows if row["sample"] == "SRR3306589")
    assert srr3306589["seqspec_id"] == "observed_ctgtaggc_3p_adapter"
    assert srr3306589["ribobase_threep_adapter_dna"] == "CTGTAGGCACCATCAAT"
    assert srr3306589["ribobase_boundary_min"] == "26"
    assert srr3306589["ribobase_boundary_max"] == "29"
