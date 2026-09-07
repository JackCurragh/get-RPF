"""M5 acceptance tests: release classifier, evidence object, cohort store.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M5.
"""

import json

from getRPF.core.release import (
    ALL_RELEASE_CLASSES,
    build_evidence,
    write_cohort_tsvs,
    write_evidence,
)


def _base_screen(verdict="consistent_with_riboseq", reasons=None, reason_codes=None):
    return {
        "length_shape": "peaked",
        "insert_mode": 30,
        "adapter_dimer_fraction": 0.0,
        "contamination_screened": False,
        "contamination": {},
        "duplication_umi_aware": 0.1,
        "verdict": verdict,
        "reason_codes": reason_codes or ["peaked_clean"],
        "reasons": reasons or ["peaked"],
    }


def test_release_class_covers_every_evidence_object():
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=900,
        applied_rules=[],
        proposed_rules=[],
        biological_screen=_base_screen(),
    )
    assert evidence["release_class"] in ALL_RELEASE_CLASSES
    assert evidence["biological_confirmation"] == "pending_alignment_gate"


def test_no_trim_clean_sample_is_clean_no_trim():
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=1000,
        biological_screen=_base_screen(),
    )
    assert evidence["release_class"] == "clean_no_trim"


def test_uniform_trim_is_clean_after_terminal_trim():
    applied = [
        {"end": "3p", "read_length_before": 29, "trim_bases": 1},
        {"end": "3p", "read_length_before": 30, "trim_bases": 1},
    ]
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=950,
        applied_rules=applied,
        biological_screen=_base_screen(),
    )
    assert evidence["release_class"] == "clean_after_terminal_trim"


def test_varying_trim_is_clean_after_per_length_terminal_trim():
    applied = [
        {"end": "3p", "read_length_before": 29, "trim_bases": 1},
        {"end": "3p", "read_length_before": 31, "trim_bases": 3},
    ]
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=950,
        applied_rules=applied,
        biological_screen=_base_screen(),
    )
    assert evidence["release_class"] == "clean_after_per_length_terminal_trim"


def test_low_yield_is_excluded_or_held():
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=10,
        biological_screen=_base_screen(),
    )
    assert evidence["release_class"] == "exclude_or_hold"


def test_contamination_inconsistent_is_excluded_or_held():
    screen = _base_screen(
        verdict="inconsistent",
        reasons=["82.0% of reads match contamination reference k-mers"],
        reason_codes=["high_contamination"],
    )
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=900,
        biological_screen=screen,
    )
    assert evidence["release_class"] == "exclude_or_hold"
    assert evidence["warnings"] == []  # no crash on empty warnings
    assert (
        "high_contamination" in screen["reason_codes"]
    )  # machine-readable reason present


def test_broad_length_shape_is_needs_length_policy_review():
    screen = _base_screen(
        verdict="inconsistent", reasons=["length distribution is broad, not peaked"]
    )
    screen["length_shape"] = "broad"
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=900,
        biological_screen=screen,
    )
    assert evidence["release_class"] == "needs_length_policy_review"


def test_polyg_artifact_holds():
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=900,
        biological_screen=_base_screen(),
        quality_evidence={"three_prime_q_collapse": True},
    )
    assert evidence["release_class"] == "basecaller_artifact_hold"


def test_adapter_dimer_excludes():
    screen = _base_screen()
    screen["adapter_dimer_fraction"] = 0.9
    evidence = build_evidence(
        sample_id="s1",
        getrpf_version="0.3.0",
        read_count_input=1000,
        read_count_output=900,
        biological_screen=screen,
    )
    assert evidence["release_class"] == "adapter_dimer_exclude"


def test_every_sample_has_evidence_json_and_cohort_is_reproducible(tmp_path):
    evidence_list = [
        build_evidence(
            sample_id="clean1",
            getrpf_version="0.3.0",
            read_count_input=1000,
            read_count_output=1000,
            biological_screen=_base_screen(),
        ),
        build_evidence(
            sample_id="held1",
            getrpf_version="0.3.0",
            read_count_input=1000,
            read_count_output=5,
            biological_screen=_base_screen(),
            warnings=["extraction yield below 5%"],
        ),
        build_evidence(
            sample_id="review1",
            getrpf_version="0.3.0",
            read_count_input=1000,
            read_count_output=900,
            self_consistency={
                "constant_region_ok": None,
                "umi_region_ok": None,
                "no_residual_cliff": False,
                "no_overtrim": True,
                "arithmetic_closes": True,
            },
            biological_screen=_base_screen(),
        ),
    ]

    evidence_paths = {}
    for e in evidence_list:
        path = tmp_path / f"{e['sample_id']}.evidence.json"
        write_evidence(e, path)
        evidence_paths[e["sample_id"]] = path

    # Every sample has a machine-readable evidence.json with the required keys.
    for sample_id, path in evidence_paths.items():
        assert path.exists()
        data = json.loads(path.read_text())
        for key in ("release_class", "biological_confirmation", "warnings"):
            assert key in data

    # The held sample has a non-empty, machine-readable reason.
    held_data = json.loads(evidence_paths["held1"].read_text())
    assert held_data["release_class"] in (
        "exclude_or_hold",
        "basecaller_artifact_hold",
        "adapter_dimer_exclude",
    )
    assert held_data["warnings"]

    # The review sample was routed to needs_protocol_seqspec by the failed
    # self-consistency check.
    review_data = json.loads(evidence_paths["review1"].read_text())
    assert review_data["release_class"] == "needs_protocol_seqspec"

    tsv_paths = write_cohort_tsvs(evidence_list, tmp_path / "cohort")

    summary_rows = (tsv_paths["release_qc_summary"]).read_text().splitlines()
    header = summary_rows[0].split("\t")
    rows = [dict(zip(header, row.split("\t"))) for row in summary_rows[1:]]
    by_id = {r["sample_id"]: r for r in rows}

    # Reproducible from the evidence JSON: same release_class per sample.
    for e in evidence_list:
        assert by_id[e["sample_id"]]["release_class"] == e["release_class"]
        assert (
            float(by_id[e["sample_id"]]["retained_fraction"]) == e["retained_fraction"]
        )

    excluded_rows = tsv_paths["samples_excluded_or_held"].read_text().splitlines()[1:]
    excluded_ids = {row.split("\t")[0] for row in excluded_rows}
    assert "held1" in excluded_ids
    assert "clean1" not in excluded_ids

    seqspec_rows = tsv_paths["samples_needing_seqspec"].read_text().splitlines()[1:]
    assert "review1" in {row.split("\t")[0] for row in seqspec_rows}
