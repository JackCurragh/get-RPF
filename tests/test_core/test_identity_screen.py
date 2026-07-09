"""M4 acceptance tests: identity screen (biological, provisional).

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M4. Never
asserts "confirmed" -- only consistent_with_riboseq / inconsistent /
indeterminate.
"""

from getRPF.core.processors.identity_screen import build_kmer_set, identity_screen

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def test_clean_footprint_library_is_consistent_with_riboseq():
    length_distribution = {}
    reads = []
    for length in (28, 29, 30, 31, 32):
        count = 180 if length in (29, 30) else 20
        for i in range(count):
            reads.append(_cycled_seq(i, length))
        length_distribution[length] = count

    result = identity_screen(reads, length_distribution)

    assert result["length_shape"] == "peaked"
    assert result["verdict"] == "consistent_with_riboseq"
    assert result["insert_mode"] in (29, 30)


def test_rnaseq_like_broad_library_is_inconsistent():
    length_distribution = {length: 20 for length in range(20, 61)}  # flat 20-60nt
    reads = [_cycled_seq(i, length) for length in length_distribution for i in range(5)]

    result = identity_screen(reads, length_distribution)

    assert result["length_shape"] == "broad"
    assert result["verdict"] == "inconsistent"
    assert any("broad" in reason for reason in result["reasons"])


def test_degradome_like_contaminated_library_is_inconsistent():
    # Peaked length distribution (footprint-like) but heavily rRNA-derived.
    length_distribution = {29: 200, 30: 700, 31: 100}

    rrna_kmer = "GGGCTACATTTTCACAGGCTA"  # 21nt synthetic rRNA reference motif
    reads = []
    for i in range(700):
        if i < 600:  # 600/700 ~ 86% carry the rRNA motif
            reads.append(rrna_kmer[:20] + "AAAAAAAAAA"[: 30 - 20])
        else:
            reads.append(_cycled_seq(i, 30))
    for i in range(200):
        reads.append(_cycled_seq(i, 29))
    for i in range(100):
        reads.append(_cycled_seq(i, 31))

    kmer_sets = {"rRNA": {rrna_kmer[i : i + 20] for i in range(len(rrna_kmer) - 19)}}

    result = identity_screen(reads, length_distribution, contamination_kmer_sets=kmer_sets)

    assert result["length_shape"] == "peaked"
    assert result["contamination_screened"] is True
    assert result["contamination"]["rRNA"] > 0.5
    assert result["verdict"] == "inconsistent"
    assert any("contamination" in reason for reason in result["reasons"])
    assert "high_contamination" in result["reason_codes"]


def test_build_kmer_set_from_reference_fasta(tmp_path):
    fasta_path = tmp_path / "rrna_reference.fasta"
    fasta_path.write_text(">rRNA_fragment\nGGGCTACATTTTCACAGGCTAAAAAA\n")

    kmers = build_kmer_set(fasta_path, k=20)

    assert len(kmers) == len("GGGCTACATTTTCACAGGCTAAAAAA") - 20 + 1
    assert "GGGCTACATTTTCACAGGCT" in kmers
    assert "AAAAAA" not in kmers  # shorter than k, can't be a k-mer

    # A read carrying one of these k-mers should be flagged by the screen.
    reads = ["GGGCTACATTTTCACAGGCTA" + "TTTTTTTT"]
    result = identity_screen(
        reads, {len(reads[0]): 1}, contamination_kmer_sets={"rRNA": kmers}
    )
    assert result["contamination"]["rRNA"] == 1.0


def test_verdict_is_never_confirmed():
    length_distribution = {30: 1000}
    reads = [_cycled_seq(i, 30) for i in range(1000)]
    result = identity_screen(reads, length_distribution)
    assert result["verdict"] != "confirmed"
    assert "confirmed" not in result["verdict"]
