"""Create test datasets for getRPF testing.

This script generates various test datasets in different formats
(FASTQ, FASTA, collapsed)
with different characteristics to test getRPF functionality.
"""

import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def ensure_output_dir(base_dir: str = "test_data") -> Path:
    """Create and return output directory."""
    out_dir = Path(base_dir)
    out_dir.mkdir(exist_ok=True)
    return out_dir


def create_rpf_read(length: int = 30, quality: int = 40) -> SeqRecord:
    """Create a single RPF-like read."""
    bases = ["A", "T", "G", "C"]
    seq = "".join(random.choices(bases, k=length))
    return SeqRecord(
        Seq(seq),
        id=f"read_{length}_{random.randint(1,1000)}",
        description="",
        letter_annotations={"phred_quality": [quality] * length},
    )


def create_good_rpf_dataset(
    output_dir: Path, n_reads: int = 1000, length_range: tuple = (28, 32)
) -> None:
    """Create dataset with good RPF characteristics."""
    reads = []
    for _ in range(n_reads):
        length = random.randint(*length_range)
        reads.append(create_rpf_read(length=length))

    # Write in different formats
    SeqIO.write(reads, output_dir / "good_rpf.fastq", "fastq")
    SeqIO.write(reads, output_dir / "good_rpf.fasta", "fasta")

    # Create collapsed FASTA (with read counts in headers)
    seq_counts = {}
    for read in reads:
        seq = str(read.seq)
        seq_counts[seq] = seq_counts.get(seq, 0) + 1

    with open(output_dir / "good_rpf.collapsed", "w") as f:
        for seq, count in seq_counts.items():
            f.write(f">read_{count}\n{seq}\n")


def create_adapter_contaminated_dataset(
    output_dir: Path,
    adapter: str = "AGATCGGAAGAG",
    contamination_rate: float = 0.3,
) -> None:
    """Create dataset with adapter contamination."""
    reads = []
    n_reads = 1000

    for i in range(n_reads):
        if random.random() < contamination_rate:
            # Create read with adapter
            base_len = random.randint(25, 35)
            base_seq = "".join(random.choices(["A", "T", "G", "C"], k=base_len))
            full_seq = base_seq + adapter[: random.randint(6, len(adapter))]
        else:
            # Create normal read
            full_seq = "".join(
                random.choices(["A", "T", "G", "C"], k=random.randint(28, 32))
            )

        reads.append(
            SeqRecord(
                Seq(full_seq),
                id=f"read_{i}",
                description="",
                letter_annotations={"phred_quality": [40] * len(full_seq)},
            )
        )

    SeqIO.write(reads, output_dir / "adapter_contaminated.fastq", "fastq")


def create_length_outliers_dataset(output_dir: Path) -> None:
    """Create dataset with abnormal length distribution."""
    reads = []

    # Add reads with various lengths
    length_distributions = {15: 100, 20: 200, 40: 300, 50: 100}

    for length, count in length_distributions.items():
        for _ in range(count):
            reads.append(create_rpf_read(length=length))

    SeqIO.write(reads, output_dir / "length_outliers.fastq", "fastq")


def create_base_biased_dataset(output_dir: Path) -> None:
    """Create dataset with position-specific base bias."""
    reads = []
    n_reads = 1000
    read_length = 30

    for i in range(n_reads):
        # Create read with strong A bias at the start
        biased_part = "A" * 10
        random_part = "".join(random.choices(["A", "T", "G", "C"], k=read_length - 10))

        reads.append(
            SeqRecord(
                Seq(biased_part + random_part),
                id=f"read_{i}",
                description="",
                letter_annotations={"phred_quality": [40] * read_length},
            )
        )

    SeqIO.write(reads, output_dir / "base_biased.fastq", "fastq")


def create_mixed_quality_dataset(output_dir: Path) -> None:
    """Create dataset with varying quality scores."""
    reads = []
    n_reads = 1000

    for i in range(n_reads):
        length = random.randint(28, 32)
        qualities = []

        # Create varying quality pattern
        for _ in range(length):
            if random.random() < 0.2:  # 20% chance of low quality
                qualities.append(random.randint(10, 20))
            else:
                qualities.append(random.randint(30, 40))

        reads.append(
            SeqRecord(
                Seq("".join(random.choices(["A", "T", "G", "C"], k=length))),
                id=f"read_{i}",
                description="",
                letter_annotations={"phred_quality": qualities},
            )
        )

    SeqIO.write(reads, output_dir / "mixed_quality.fastq", "fastq")


def main():
    """Generate all test datasets."""
    out_dir = ensure_output_dir()

    # Create different test datasets
    create_good_rpf_dataset(out_dir)
    create_adapter_contaminated_dataset(out_dir)
    create_length_outliers_dataset(out_dir)
    create_base_biased_dataset(out_dir)
    create_mixed_quality_dataset(out_dir)

    print(f"Created test datasets in {out_dir}")
    print("\nTest commands:")
    print("\n1. Check good RPF data:")
    print(
        "getRPF check test_data/good_rpf.fastq \
            --format fastq \
            --output good_rpf_report.txt"
    )

    print("\n2. Check adapter contamination:")
    print(
        "getRPF detect-adapter test_data/adapter_contaminated.fastq \
            --format fastq --adapter AGATCGGAAGAG --output adapter_report.txt"
    )

    print("\n3. Check length distribution issues:")
    print(
        "getRPF check test_data/length_outliers.fastq --format fastq \
            --output length_report.txt"
    )

    print("\n4. Check base composition bias:")
    print(
        "getRPF check test_data/base_biased.fastq --format fastq \
            --output bias_report.txt"
    )

    print("\n5. Check different input formats:")
    print(
        "getRPF check test_data/good_rpf.fasta --format fasta \
            --output fasta_report.txt"
    )
    print(
        "getRPF check test_data/good_rpf.collapsed --format collapsed \
            --output collapsed_report.txt"
    )


if __name__ == "__main__":
    main()
