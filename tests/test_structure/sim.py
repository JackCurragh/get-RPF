"""Synthetic Ribo-seq reads with known structure (spec §10.1).

Footprints are drawn from a random transcriptome with Zipf-skewed abundance,
so fragments recur as they do in real libraries. Each read is built as
5' technical + insert + 3' technical + adapter + downstream construct, then
cut to the read length, so every read carries its true insert coordinates.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Element = Callable[[random.Random], str]

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
# i7 index, P7, then no-signal cycles (poly-G on two-colour instruments).
DOWNSTREAM = "GATCAGTC" + "ATCTCGTATGCCGTCTTCTGCTTG" + "G" * 100


def random_bases(n: int) -> Element:
    return lambda rng: "".join(rng.choices("ACGT", k=n))


def fixed(sequence: str) -> Element:
    return lambda rng: sequence


def barcode(values: Sequence[str]) -> Element:
    choices = list(values)
    return lambda rng: rng.choice(choices)


def nta(lengths: Dict[int, float], bases: Dict[str, float]) -> Element:
    """Non-templated addition: variable length, skewed base composition."""
    length_values, length_weights = list(lengths), list(lengths.values())
    base_values, base_weights = list(bases), list(bases.values())

    def make(rng: random.Random) -> str:
        n = rng.choices(length_values, length_weights)[0]
        return "".join(rng.choices(base_values, base_weights, k=n))

    return make


@dataclass(frozen=True)
class Library:
    reads: List[str]
    anchors: List[Optional[int]]
    """Adapter start per read; None where the read ends before the adapter."""
    inserts: List[Tuple[int, int]]
    """True insert [start, end) per read."""


def simulate(
    five_prime: Sequence[Element] = (),
    three_prime: Sequence[Element] = (),
    *,
    n_reads: int = 30_000,
    read_length: Optional[int] = 75,
    n_transcripts: int = 300,
    transcript_length: int = 1200,
    zipf_s: float = 1.1,
    footprint: Tuple[int, int] = (28, 32),
    error_rate: float = 0.001,
    seed: int = 0,
) -> Library:
    """Build a library; ``read_length=None`` gives already-trimmed reads."""
    rng = random.Random(seed)
    transcripts = [
        "".join(rng.choices("ACGT", k=transcript_length)) for _ in range(n_transcripts)
    ]
    weights = [1 / (i + 1) ** zipf_s for i in range(n_transcripts)]
    lengths = list(range(footprint[0], footprint[1] + 1))

    reads: List[str] = []
    anchors: List[Optional[int]] = []
    inserts: List[Tuple[int, int]] = []
    for t in rng.choices(range(n_transcripts), weights, k=n_reads):
        size = rng.choice(lengths)
        at = rng.randrange(transcript_length - size + 1)
        five = "".join(element(rng) for element in five_prime)
        three = "".join(element(rng) for element in three_prime)
        body = five + transcripts[t][at : at + size] + three
        if read_length is None:
            read, anchor = body, len(body)
        else:
            read = (body + ADAPTER + DOWNSTREAM)[:read_length]
            anchor = len(body) if len(body) < read_length else None
        if rng.random() < 1 - (1 - error_rate) ** len(read):
            i = rng.randrange(len(read))
            read = (
                read[:i]
                + rng.choice([b for b in "ACGT" if b != read[i]])
                + read[i + 1 :]
            )
        reads.append(read)
        anchors.append(anchor)
        inserts.append((len(five), len(five) + size))
    return Library(reads, anchors, inserts)
