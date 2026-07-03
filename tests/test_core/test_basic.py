from collections import Counter

from getRPF.core.processors.collapsed import TwoStageCollapser
from getRPF.core.processors.rpf_extractor import RPFExtractor
from getRPF.core.processors.types import ExtractionEmptyError, ReadArchitecture


def test_initial_setup():
    assert True


def test_two_stage_collapser_can_apply_upper_length_bound():
    collapser = TwoStageCollapser()
    raw_counts = Counter(
        {
            "A" * 19: 1,
            "C" * 20: 2,
            "G" * 40: 3,
            "T" * 41: 4,
        }
    )

    final_counts = collapser.apply_trimming(
        raw_counts,
        lambda seq: seq,
        min_length=20,
        max_length=40,
    )

    assert final_counts == Counter({"C" * 20: 2, "G" * 40: 3})


def test_adapter_prefix_match_finds_partial_adapter():
    sequence = "A" * 32 + "TGGAATTCTC" + "C" * 18
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) == 32


def test_adapter_prefix_match_requires_minimum_overlap():
    sequence = "A" * 32 + "TGGAATTC" + "C" * 20
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) is None


def test_adapter_evidence_prefers_observed_protocol_over_generic_catalog():
    extractor = RPFExtractor()
    adapter = "AGATCGGAAGAGCACACGTCT"
    generic = ReadArchitecture(
        protocol_name="comprehensive_adapter_check",
        lab_source="Generic",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[adapter],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    observed = ReadArchitecture(
        protocol_name="observed_truseq_21nt_3p_adapter",
        lab_source="Observed protocol catalog",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[adapter],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    extractor.architecture_db.architectures = [generic, observed]

    reads = [
        "A" * 28 + adapter,
        "C" * 31 + adapter[:14],
        "G" * 32,
    ]

    evidence = extractor._score_adapter_evidence(reads)

    assert evidence[0]["protocol_name"] == "observed_truseq_21nt_3p_adapter"
    assert evidence[0]["hit_count"] == 2
    assert evidence[0]["hit_fraction"] == 2 / 3


def test_adapter_reporting_flags_unsupported_selected_architecture_conflict():
    extractor = RPFExtractor()
    illumina_adapter = "AGATCGGAAGAGCACACGTCT"
    selected_adapter = "CTGTAGGCACCATCAAT"
    illumina = ReadArchitecture(
        protocol_name="comprehensive_adapter_check",
        lab_source="Generic",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[illumina_adapter],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    selected = ReadArchitecture(
        protocol_name="observed_ctgtaggc_3p_adapter",
        lab_source="Observed protocol catalog",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[selected_adapter],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    extractor.architecture_db.architectures = [illumina, selected]
    reads = [
        "A" * 28 + illumina_adapter,
        "C" * 31 + illumina_adapter[:14],
        "G" * 32 + illumina_adapter[:12],
    ]

    evidence = extractor._score_adapter_evidence(reads)
    report = extractor._adapter_reporting(selected, "strict_pattern_match", evidence)

    assert report["adapter_source"] == "architecture_match"
    assert report["adapter_conflict"]["has_conflict"] is True
    assert report["adapter_conflict"]["type"] == "selected_architecture_unsupported"
    assert report["adapter_conflict"]["best_supported_protocol"] == (
        "comprehensive_adapter_check"
    )
    assert report["adapter_evidence_candidates"][0]["adapter_source"] == (
        "generic_adapter_catalog"
    )


def test_adapter_reporting_marks_de_novo_without_dominant_adapter():
    extractor = RPFExtractor()
    arch = ReadArchitecture(
        protocol_name="de_novo_inferred",
        lab_source="Probabilistic Segmenter",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )

    report = extractor._adapter_reporting(arch, "probabilistic_hmm", [])

    assert report["adapter_source"] == "no_dominant_adapter"
    assert report["adapter_conflict"]["has_conflict"] is False
    assert report["adapter_evidence_candidates"] == []


def test_pretrimmed_length_filter_policy_requires_footprint_like_reads():
    extractor = RPFExtractor()

    footprint_profile = extractor._length_profile(
        ["A" * 28] * 60 + ["C" * 32] * 25 + ["G" * 60] * 15
    )
    long_profile = extractor._length_profile(["A" * 76] * 90 + ["C" * 30] * 10)

    assert extractor._looks_pretrimmed(footprint_profile) is True
    assert footprint_profile["frac_20_40"] == 0.85
    assert extractor._looks_pretrimmed(long_profile) is False


def test_extract_rpfs_reports_pretrimmed_length_filter_without_seqspec(tmp_path):
    input_file = tmp_path / "pretrimmed.fastq"
    reads = []
    for i in range(85):
        seq = ("ACGT" * 10)[: 28 + (i % 5)]
        reads.append(f"@rpf{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    for i in range(15):
        seq = "T" * 60
        reads.append(f"@long{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    input_file.write_text("".join(reads))

    extractor = RPFExtractor()
    extractor.architecture_db.architectures = []
    output_file = tmp_path / "out.fastq"

    result = extractor.extract_rpfs(
        input_file,
        output_file,
        format="fastq",
        collapsed_only=True,
    )

    assert result.architecture_match == "pretrimmed_no_dominant_adapter"
    assert result.extraction_method == "pretrimmed_length_filter"
    assert result.adapter_source == "no_dominant_adapter"
    assert result.input_reads == 100
    assert result.extracted_rpfs == 85
    assert result.quality_metrics["retained_fraction"] == 0.85
    assert result.quality_metrics["extraction_class"] == "pretrimmed_rpf_length_filter"
    assert result.quality_metrics["raw_length_profile"]["frac_gt40"] == 0.15
    assert result.quality_metrics["extracted_length_profile"]["frac_20_40"] == 1.0


def test_pretrimmed_length_filter_has_no_adapter_conflict():
    extractor = RPFExtractor()
    pretrimmed = ReadArchitecture(
        protocol_name="pretrimmed_no_dominant_adapter",
        lab_source="Length-filtered raw reads",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )

    report = extractor._adapter_reporting(
        pretrimmed,
        "pretrimmed_length_filter",
        [],
    )

    assert report["adapter_source"] == "no_dominant_adapter"
    assert report["adapter_conflict"]["has_conflict"] is False
    assert (
        extractor._extraction_class("pretrimmed_length_filter", report, [])
        == "pretrimmed_rpf_length_filter"
    )


def test_extraction_class_flags_low_yield_before_success_class():
    extractor = RPFExtractor()
    report = {
        "adapter_conflict": {"has_conflict": False},
        "adapter_source": "raw_adapter_evidence",
    }

    assert (
        extractor._extraction_class(
            "adapter_evidence_match",
            report,
            ["low_retained_fraction"],
        )
        == "adapter_supported_low_yield"
    )
    assert (
        extractor._extraction_class(
            "adapter_evidence_match",
            report,
            ["low_26_34_fraction"],
        )
        == "extracted_rpf_length_suspicious"
    )


def test_extract_rpfs_raises_when_no_rpfs_survive(tmp_path):
    input_file = tmp_path / "no_rpfs.fastq"
    seq = "A" * 60
    input_file.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i in range(20)))

    extractor = RPFExtractor()
    extractor.architecture_db.architectures = []

    try:
        extractor._extract_rpfs_from_reads(
            input_file=input_file,
            output_file=tmp_path / "out.fastq",
            segments=[],
            format="fastq",
            max_reads=None,
            adapters=[],
            collapsed_only=True,
        )
    except ExtractionEmptyError as exc:
        assert "zero RPF reads" in str(exc)
    else:
        raise AssertionError("Expected ExtractionEmptyError")


def test_adapter_evidence_selection_has_no_external_prior_conflict():
    extractor = RPFExtractor()
    illumina_adapter = "AGATCGGAAGAGCACACGTCT"
    raw_supported = ReadArchitecture(
        protocol_name="mcglincy_ingolia_2017",
        lab_source="Known protocol",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=[illumina_adapter],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    extractor.architecture_db.architectures = [raw_supported]
    evidence = extractor._score_adapter_evidence(
        [
            "A" * 28 + illumina_adapter,
            "C" * 31 + illumina_adapter[:14],
            "G" * 32 + illumina_adapter[:12],
        ]
    )

    report = extractor._adapter_reporting(
        raw_supported,
        "adapter_evidence_match",
        evidence,
    )

    assert report["adapter_source"] == "raw_adapter_evidence"
    assert report["adapter_conflict"]["has_conflict"] is False
    assert set(report["adapter_conflict"]) == {
        "has_conflict",
        "type",
        "selected_protocol",
        "selected_hit_fraction",
        "best_supported_protocol",
        "best_supported_hit_fraction",
        "message",
    }
    assert report["adapter_conflict"]["best_supported_protocol"] == "mcglincy_ingolia_2017"
