# getRPF Documentation

`getRPF` is a tool for analyzing Ribosome Protected Fragments (RPFs) from Ribo-seq experiments. It provides comprehensive quality control and adapter detection capabilities.

## Quick Start

```bash
# Install the package
pip install getRPF

# Check read quality
getRPF check input.fastq --format fastq --output quality_report.txt

# Detect adapters
getRPF detect-adapter input.fastq --format fastq --adapter AGATCGGAAGAG --output adapter_report.txt
```

## Contents

```{toctree}
:maxdepth: 2

usage
rpf_checks
file_formats
api
```
