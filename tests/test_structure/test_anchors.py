"""Step 2 behaviours that real data cannot isolate (spec §10.1 cases 11, 12,
15, 16 and tails). Real-library validation is in test_real_panel.py."""

from getRPF.core.structure import Status
from getRPF.core.structure.anchors import find_anchor, probe_short_end
from getRPF.core.structure.observe import observe

from .sim import nta, simulate


def anchor(reads):
    return find_anchor(reads, observe(reads))


def truth_agreement(call, truth):
    pairs = [
        (got, want) for got, want in zip(call.insert_ends, truth) if want is not None
    ]
    return sum(got == want for got, want in pairs) / len(pairs)


def test_adapter_mid_read_is_found_and_p7_is_never_the_anchor():
    library = simulate(read_length=100, n_reads=10_000)
    answer = anchor(library.reads)
    assert answer.status is Status.RESOLVED
    assert answer.value.adapter.name == "truseq_3p"
    assert answer.value.downstream_consistency >= 0.95
    assert truth_agreement(answer.value, library.anchors) >= 0.98
    assert any(
        alt.value == "illumina_p7" and alt.verdict == "contradicted"
        for alt in answer.alternatives
    )


def test_unknown_adapter_is_discovered_de_novo():
    library = simulate(adapter="TCGTACGGAGTTCAGCATGCACTGA", n_reads=10_000)
    answer = anchor(library.reads)
    assert answer.status is Status.RESOLVED
    assert answer.value.source == "de_novo"
    assert answer.value.adapter.sequence.startswith("TCGTACGGAG")
    assert truth_agreement(answer.value, library.anchors) >= 0.95


def test_two_adapters_make_q1_ambiguous():
    first = simulate(n_reads=6_000, seed=1)
    second = simulate(n_reads=4_000, seed=2, adapter="TGGAATTCTCGGGTGCCAAGG")
    answer = anchor(first.reads + second.reads)
    assert answer.status is Status.AMBIGUOUS


def test_reads_shorter_than_the_insert_have_no_anchor():
    answer = anchor(simulate(read_length=26, n_reads=5_000).reads)
    assert answer.status is Status.NOT_OBSERVABLE
    assert answer.value is None


def test_trimmed_reads_have_no_anchor():
    answer = anchor(simulate(read_length=None, n_reads=5_000).reads)
    assert answer.status is Status.NOT_OBSERVABLE


def test_fixed_length_35nt_is_a_recovery_candidate_not_raw_structure():
    reads = ["C" * 28 + "AGATCGG"] * 100
    observed = observe(reads)
    assert observed.input_state == "fixed_length_footprint_candidate"
    # Five-to-nine bases are evidence for recovery only; Q1 remains closed.
    assert find_anchor(reads, observed).status is Status.NOT_OBSERVABLE
    candidates = probe_short_end(reads)
    assert candidates[0].adapter.name == "truseq_3p"
    assert candidates[0].overlap == 7
    assert candidates[0].support == 1.0


def test_poly_a_tail_places_the_insert_end_at_the_run_start():
    tail = nta({n: 1.0 for n in range(8, 25)}, {"A": 1.0})
    library = simulate(three_prime=[tail], read_length=100, n_reads=10_000)
    answer = anchor(library.reads)
    assert answer.status is Status.INTERVAL
    assert answer.value.tail_base == "A"
    ends = [end for _, end in library.inserts]
    pairs = [
        (got, want)
        for got, want, located in zip(answer.value.insert_ends, ends, library.anchors)
        if got is not None and located is not None
    ]
    # The convention cuts early where an insert ends in A. It cuts late only
    # when a sequencing error hits the run's first base, which is rare.
    assert sum(got > want for got, want in pairs) / len(pairs) <= 0.005
    assert sum(got == want for got, want in pairs) / len(pairs) >= 0.7
