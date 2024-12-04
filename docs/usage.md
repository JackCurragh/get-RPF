# Usage Guide

## Installation

Install getRPF using pip:

```bash
pip install getRPF
```

## Basic Commands

getRPF provides two main commands:

### Check Command

The `check` command analyzes read quality and composition:

```bash
getRPF check input.fastq --format fastq --output report.txt
```

Options:
- `--format`: Input file format (fastq/fasta/collapsed)
- `--output`: Path for the output report
- `--min-quality`: Minimum quality score threshold (default: 20)
- `--threads`: Number of processing threads (default: 1)

### Adapter Detection

The `detect-adapter` command identifies adapter sequences:

```bash
getRPF detect-adapter input.fastq --format fastq --adapter AGATCGGAAGAG --output report.txt
```

Options:
- `--format`: Input file format (fastq/fasta/collapsed)
- `--output`: Path for the output report
- `--adapter`: Adapter sequence to search for
- `--min-overlap`: Minimum overlap for adapter matching (default: 10)
- `--max-mismatches`: Maximum allowed mismatches (default: 1)

## Example Usage

### Checking RPF Data Quality

```bash
# Basic quality check
getRPF check reads.fastq --format fastq --output quality_report.txt

# Check with custom quality threshold
getRPF check reads.fastq --format fastq --min-quality 25 --output report.txt

# Process collapsed FASTA
getRPF check reads.collapsed --format collapsed --output report.txt
```

### Detecting Adapters

```bash
# Basic adapter detection
getRPF detect-adapter reads.fastq --format fastq --adapter AGATCGGAAGAG --output report.txt

# Strict adapter matching
getRPF detect-adapter reads.fastq --format fastq --adapter AGATCGGAAGAG --min-overlap 15 --max-mismatches 0 --output report.txt
```
