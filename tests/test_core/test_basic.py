from collections import Counter

from getRPF.core.processors.collapsed import TwoStageCollapser
from getRPF.core.processors.matcher import MatchResult
from getRPF.core.processors.rpf_extractor import (
    RPFExtractor,
    resolve_architecture_choice,
)
from getRPF.core.processors.types import ReadArchitecture


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
    sequence = "A" * 32 + "TGGAATTCTC"
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) == 32


def test_adapter_prefix_match_requires_minimum_overlap():
    sequence = "A" * 32 + "TGGAATTC"
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) is None


def test_adapter_prefix_match_rejects_internal_match():
    sequence = "A" * 32 + "TGGAATTCTC" + "C" * 18
    adapter = "TGGAATTCTCGGGTGCCAAGG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) is None


def test_adapter_prefix_match_supports_short_observed_adapter():
    sequence = "A" * 29 + "AGATCGGAG"
    adapter = "AGATCGGAG"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) == 29


def test_adapter_prefix_match_treats_adapter_n_as_wildcard():
    sequence = "A" * 30 + "GATTACCACTCGGGCACCAAGGA"
    adapter = "NNNNNNCACTCGGGCACCAAGGA"

    assert RPFExtractor._find_adapter_prefix(sequence, adapter) == 30


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


def test_specific_adapter_evidence_overrides_generic_strict_match():
    generic = ReadArchitecture(
        protocol_name="comprehensive_adapter_check",
        lab_source="Generic",
        umi_positions=[],
        barcode_positions=[],
        adapter_sequences=["AAAAAAAAAA"],
        rpf_start=0,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    riboflow = ReadArchitecture(
        protocol_name="riboflow_template_switch",
        lab_source="RiboFlow",
        umi_positions=[(0, 12)],
        barcode_positions=[],
        adapter_sequences=["AAAAAAAAAACAAAAAAAAAA"],
        rpf_start=16,
        rpf_end=-1,
        expected_rpf_length=(20, 40),
        quality_markers={},
    )
    match = MatchResult(generic, True, 1.0, ["generic match"])
    evidence = [
        {
            "architecture": riboflow,
            "protocol_name": riboflow.protocol_name,
            "adapter": riboflow.adapter_sequences[0],
            "hit_count": 80,
            "hit_fraction": 0.8,
        },
        {
            "architecture": generic,
            "protocol_name": generic.protocol_name,
            "adapter": generic.adapter_sequences[0],
            "hit_count": 80,
            "hit_fraction": 0.8,
        },
    ]

    architecture, method, trace = resolve_architecture_choice(match, evidence)

    assert architecture is riboflow
    assert method == "adapter_evidence_match"
    assert "riboflow_template_switch" in trace[0]


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


def test_length_profile_surfaces_disome_candidate_lengths_without_filtering():
    extractor = RPFExtractor()

    profile = extractor._length_profile(["A" * 60] * 80 + ["C" * 30] * 20)

    assert profile["frac_50_80"] == 0.8
    warnings = extractor._extraction_warnings(
        {"retained_fraction": 1.0, "extracted_length_profile": profile}
    )
    assert "disome_like_length_distribution" in warnings


def test_disome_candidate_lengths_require_explicit_policy_review():
    extractor = RPFExtractor()
    report = {"adapter_conflict": {"has_conflict": False}}

    assert (
        extractor._extraction_class(
            "adapter_evidence_match",
            report,
            ["disome_like_length_distribution"],
        )
        == "needs_length_policy_review"
    )


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
    assert result.extracted_rpfs == 100
    assert result.quality_metrics["retained_fraction"] == 1.0
    assert result.quality_metrics["rpf_length_gated_reads"] == 85
    assert result.quality_metrics["rpf_length_gated_fraction"] == 0.85
    assert result.quality_metrics["extraction_class"] == "pretrimmed_rpf_length_filter"
    assert result.quality_metrics["raw_length_profile"]["frac_gt40"] == 0.15
    assert result.quality_metrics["extracted_length_profile"]["frac_20_40"] == 0.85
    assert result.quality_metrics["rpf_length_profile"]["frac_20_40"] == 1.0


def test_extract_rpfs_trims_only_post_rpf_adapter_for_dual_ligation(tmp_path):
    input_file = tmp_path / "dual.fastq"
    five_prime_adapter = "GTTCAGAGTTCTACAGTCCGACGATC"
    three_prime_adapter = "TCGTATGCCGTCTTCTGCTTG"
    rpf = "ACGTACGTACGTACGTACGTACGTACGT"
    reads = []
    for i in range(20):
        seq = f"{five_prime_adapter}{rpf}{three_prime_adapter}"
        reads.append(f"@dual{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    input_file.write_text("".join(reads))

    extractor = RPFExtractor()
    dual = next(
        arch
        for arch in extractor.architecture_db.architectures
        if arch.protocol_name == "observed_dual_ligation_adapter_pair"
    )
    extractor.architecture_db.architectures = [dual]
    output_file = tmp_path / "out.fastq"

    result = extractor.extract_rpfs(
        input_file,
        output_file,
        format="fastq",
        collapsed_only=True,
    )

    assert result.architecture_match == "observed_dual_ligation_adapter_pair"
    assert result.extraction_method in {
        "strict_pattern_match",
        "adapter_evidence_match",
    }
    assert result.extracted_rpfs == 20
    assert result.quality_metrics["unique_extracted_sequences"] == 1
    assert (tmp_path / "out.collapsed.fa").read_text().splitlines()[1] == rpf


def test_extract_rpfs_removes_post_rpf_technical_bases(tmp_path):
    input_file = tmp_path / "post_rpf_technical.fastq"
    rpf = "ACGTACGTACGTACGTACGTACGTACGT"
    technical_tail = "NNNNN".replace("N", "A") + "CCCCC"
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    reads = [rpf + technical_tail + adapter] * 12
    input_file.write_text(
        "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads))
    )

    extractor = RPFExtractor()
    stats = extractor._extract_rpfs_from_reads(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        segments=[],
        format="fastq",
        max_reads=None,
        adapters=[adapter],
        post_rpf_trim_bases=10,
        collapsed_only=True,
    )

    assert stats["extracted_rpfs"] == 12
    assert stats["extracted_length_profile"]["mode_length"] == len(rpf)


def test_terminal_override_does_not_cut_rpf_when_adapter_is_known(tmp_path):
    input_file = tmp_path / "known_adapter_with_override.fastq"
    rpf = "ACGTACGTACGTACGTACGTACGTACGT"
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    reads = [rpf + adapter] * 12
    input_file.write_text(
        "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads))
    )

    stats = RPFExtractor()._extract_rpfs_from_reads(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        segments=[],
        format="fastq",
        max_reads=None,
        adapters=[adapter],
        override_trims={len(reads[0]): {"trim_5p": 0, "trim_3p": 1}},
        collapsed_only=True,
    )

    assert stats["extracted_length_profile"]["mode_length"] == len(rpf)


def test_terminal_override_trims_tail_when_no_adapter_is_found(tmp_path):
    input_file = tmp_path / "unknown_adapter_with_override.fastq"
    rpf = "ACGTACGTACGTACGTACGTACGTACGT"
    reads = [rpf + "AAA"] * 12
    input_file.write_text(
        "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads))
    )

    stats = RPFExtractor()._extract_rpfs_from_reads(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        segments=[],
        format="fastq",
        max_reads=None,
        adapters=["AGATCGGAAGAGCACACGTCT"],
        override_trims={len(reads[0]): {"trim_5p": 0, "trim_3p": 3}},
        collapsed_only=True,
    )

    assert stats["extracted_length_profile"]["mode_length"] == len(rpf)


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


def test_extract_rpfs_preserves_trimmed_reads_outside_rpf_window(tmp_path):
    input_file = tmp_path / "long_trimmed.fastq"
    seq = "A" * 60
    input_file.write_text(
        "".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i in range(20))
    )

    extractor = RPFExtractor()
    extractor.architecture_db.architectures = []

    stats = extractor._extract_rpfs_from_reads(
        input_file=input_file,
        output_file=tmp_path / "out.fastq",
        segments=[],
        format="fastq",
        max_reads=None,
        adapters=[],
        collapsed_only=True,
    )

    assert stats["extracted_rpfs"] == 20
    assert stats["rpf_length_gated_reads"] == 0
    assert (tmp_path / "out.collapsed.fa").exists()


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
    assert (
        report["adapter_conflict"]["best_supported_protocol"] == "mcglincy_ingolia_2017"
    )
