"""Assembly: Q4, the transform decision, the architecture and the reports."""

import json

from getRPF.core.structure import Status
from getRPF.core.structure.assemble import infer_structure
from getRPF.core.structure.report import format_text, write_report

from .sim import nta, random_bases, simulate


def test_umi_library_is_emitted_with_the_umi_kept():
    result = infer_structure(simulate(five_prime=[random_bases(5)]).reads)
    assert result.transform.emit
    assert (
        result.architecture.describe()
        == "[random 5 UMI][insert][adapter AGATCGGAAGAG...]"
    )
    assert result.q4.status is Status.RESOLVED
    assert result.q4.value["inline"] == (("Q2", 5),)


def test_common_junction_bases_withhold_the_transform():
    skewed = {"T": 0.7, "C": 0.1, "A": 0.1, "G": 0.1}
    result = infer_structure(simulate(five_prime=[nta({2: 1.0}, skewed)]).reads)
    assert not result.transform.emit
    assert any(reason.startswith("Q2:") for reason in result.transform.reasons)


def test_trimmed_library_uses_the_read_end_as_anchor():
    result = infer_structure(simulate(read_length=None).reads)
    assert result.transform.emit
    assert any("read end is the anchor" in c for c in result.transform.conventions)
    assert result.q4.status is Status.NOT_OBSERVABLE


def test_underpowered_library_is_withheld_and_claims_no_absence():
    library = simulate(n_reads=20_000, n_transcripts=5_000, zipf_s=0.0)
    result = infer_structure(library.reads)
    assert not result.transform.emit
    assert result.q4.status is Status.UNDERPOWERED


def test_reports_are_written_and_json_ready(tmp_path):
    result = infer_structure(simulate(five_prime=[random_bases(5)]).reads)
    json_path, text_path = write_report(result, tmp_path, "sim")
    report = json.loads(json_path.read_text())
    assert report["architecture"]["describe"].startswith("[random 5 UMI]")
    assert report["answers"]["Q2"]["alternatives"]
    assert "insert_ends" not in report["answers"]["Q1"]["value"]
    assert "Architecture:" in text_path.read_text()
    assert "Transform: emitted" in format_text(result, "sim")
