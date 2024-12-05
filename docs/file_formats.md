# Supported File Formats in **getRPF**

**getRPF** supports a variety of input file formats to accommodate different sequencing workflows. This flexibility ensures compatibility with raw, processed, and deduplicated sequencing data.

---

## File Format Details

### 1. **FASTQ Format**
The **FASTQ** format is a standard for raw sequencing data. Each record consists of four lines:
1. **Header**: Starts with `@` and contains the read ID.
2. **Sequence**: The DNA or RNA sequence.
3. **Separator**: A `+` character, optionally followed by the same read ID.
4. **Quality Scores**: Encoded scores for each nucleotide in the sequence.

#### Example:
```text
@read_id
SEQUENCE
+
QUALITY_SCORES
```

Example FASTQ record:
```text
@read1
ATCGATCGATCGATCG
+
FFFFFFFFFFFFFFFFA
```

#### Key Features:
- Most commonly used for raw sequencing data.
- Includes **quality scores** for each nucleotide.
- Occupies **four lines per record**.

---

### 2. **FASTA Format**
The **FASTA** format is a simpler file format that contains only sequence data. Each record consists of two lines:
1. **Header**: Starts with `>` followed by the read ID.
2. **Sequence**: The DNA or RNA sequence.

#### Example:
```text
>read_id
SEQUENCE
```

Example FASTA record:
```text
>read1
ATCGATCGATCGATCG
```

#### Key Features:
- Contains **sequence-only data**.
- **No quality information** is included.
- Occupies **two lines per record**.

---

### 3. **Collapsed FASTA Format**
The **Collapsed FASTA** format is a variation of the standard FASTA format. It is typically used for deduplicated data where headers include a count of reads.

#### Example:
```text
>sequence_count_N
SEQUENCE
```

Example Collapsed FASTA record:
```text
>read_42
ATCGATCGATCGATCG
```

#### Key Features:
- Used for **processed or deduplicated data**.
- The header includes a **read count** (e.g., `read_42` indicates 42 reads).
- Facilitates downstream analysis by reducing redundancy.

---

## Handling Compressed Files
All supported file formats can be compressed with gzip (`.gz` extension). This allows for efficient storage and processing of large datasets.

#### Example Command:
```bash
getRPF check input.fastq.gz --format fastq --output report.txt
```

---

## Specifying the Input Format

Use the `--format` option to specify the file format when running **getRPF**. Supported formats include:
- `fastq` (for FASTQ files)
- `fasta` (for standard FASTA files)
- `collapsed` (for Collapsed FASTA files)

### Example Commands:

#### For FASTQ Format:
```bash
getRPF check input.fastq --format fastq --output report.txt
```

#### For FASTA Format:
```bash
getRPF check input.fasta --format fasta --output report.txt
```

#### For Collapsed FASTA Format:
```bash
getRPF check input.collapsed --format collapsed --output report.txt --count-pattern "read{id}_x{count}
```

---

## Summary Table

| **Format**          | **File Extension**   | **Features**                             | **Lines per Record** |
|----------------------|----------------------|------------------------------------------|-----------------------|
| FASTQ               | `.fastq`, `.fastq.gz`| Sequence and quality scores              | 4                     |
| FASTA               | `.fasta`, `.fa.gz`   | Sequence-only                            | 2                     |
| Collapsed FASTA     | `.collapsed`, `.fa.gz`| Sequence with read counts in the header  | 2                     |

---
