"""M1 acceptance tests: Stage-0 sketch with per-length and quality diagnostics.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M1.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from getRPF.core.processors.sketch import SketchBuilder, classify_terminal_signal

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def _write_fastq(tmp_path, name, sequences, qualities=None):
    records = []
    for i, seq in enumerate(sequences):
        quals = qualities[i] if qualities else [40] * len(seq)
        records.append(
            SeqRecord(
                Seq(seq),
                id=f"read{i}",
                description="",
                letter_annotations={"phred_quality": quals},
            )
        )
    path = tmp_path / name
    SeqIO.write(records, path, "fastq")
    return path


def test_per_length_artifact_invisible_in_pooled_view(tmp_path):
    """A 3' terminal artifact confined to one length class should be diluted
    away in the pooled view but clearly visible per-length."""
    sequences = []
    for i in range(900):
        sequences.append(_cycled_seq(i, 28))
    for i in range(900, 1000):
        body = _cycled_seq(i, 31)
        tail = "A" if i < 999 else "C"  # 99/100 = 99% terminal A
        sequences.append(body + tail)

    path = _write_fastq(tmp_path, "artifact.fastq", sequences)
    sketch = SketchBuilder(max_reads=2000).build_from_file(path, format="fastq")

    pooled_pos0 = sketch.pooled.composition_3p[0]
    pooled_max = max(pooled_pos0.values())
    assert pooled_max < 0.5  # diluted by the 900 diverse 28-mers

    assert 32 in sketch.per_length
    per_length_pos0 = sketch.per_length[32].composition_3p[0]
    assert per_length_pos0["A"] >= 0.95
    assert sketch.per_length_support[32] == 100


def test_classify_terminal_signal_distinguishes_polyg_from_adapter():
    # Low entropy + G-dominant + collapsed quality -> base-caller artifact.
    assert (
        classify_terminal_signal(entropy=0.05, quality_mean=4.0, dominant_base="G")
        == "basecaller_artifact"
    )
    # Low entropy + normal quality -> a real constant/adapter region.
    assert (
        classify_terminal_signal(entropy=0.05, quality_mean=36.0, dominant_base="G")
        == "adapter_or_constant"
    )
    # High entropy -> biological, regardless of quality.
    assert (
        classify_terminal_signal(entropy=1.9, quality_mean=36.0, dominant_base="A")
        == "biological"
    )


def test_sketch_quality_profile_flags_polyg_dark_cycle(tmp_path):
    """A poly-G run co-occurring with a quality collapse should be
    identifiable purely from the sketch's composition + quality profile."""
    sequences = []
    qualities = []
    for i in range(500):
        body = _cycled_seq(i, 27)
        seq = body + "G"  # constant terminal G, mimics an adapter cliff
        sequences.append(seq)
        # Quality collapses specifically at the terminal G (dark cycle).
        qualities.append([38] * 27 + [3])

    path = _write_fastq(tmp_path, "polyg.fastq", sequences, qualities)
    sketch = SketchBuilder(max_reads=1000).build_from_file(path, format="fastq")

    entropy_tail = sketch.pooled.entropy_3p[0]
    dominant_base = max(sketch.pooled.composition_3p[0].items(), key=lambda kv: kv[1])[
        0
    ]
    quality_tail = sketch.quality_3p[0]

    assert entropy_tail < 0.5
    assert dominant_base == "G"
    assert quality_tail < 10.0
    assert (
        classify_terminal_signal(entropy_tail, quality_tail, dominant_base)
        == "basecaller_artifact"
    )
