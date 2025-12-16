# getRPF - Automated Ribosome Protected Fragment Analysis

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

getRPF is a comprehensive tool for automated analysis and extraction of Ribosome Protected Fragments (RPFs) from Ribo-seq experiments. It combines intelligent architecture detection, quality assessment, and extensible protocol support to handle diverse ribosome profiling datasets at scale.

## 🚀 Key Features

- **🔍 Automated Architecture Detection**: Pattern matching against known protocols + de novo detection for novel structures
- **📊 Comprehensive Quality Assessment**: Multi-dimensional cleanliness checking with failure categorization
- **🧬 seqspec Integration**: Load novel protocols from standard seqspec files for unlimited extensibility  
- **⚡ Scalable Processing**: Designed to handle thousands of samples with efficient batch processing workflows
- **⚡ Scalable Processing**: Designed to handle thousands of samples with efficient batch processing workflows
- **📈 Insightful Reporting**: Interactive HTML reports with architecture flowcharts and alignment sparklines
- **🧠 Intelligent Detection**: Strict pattern matching + Probabilistic HMM segmentation for novel reads

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands Overview](#commands-overview)
- [User Workflows](#user-workflows)
- [Input/Output Formats](#inputoutput-formats)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)

## 🔧 Installation

### Prerequisites

- Python 3.10 or higher
- Optional: STAR aligner for alignment-based quality checks
- Optional: conda/mamba for environment management

### Install from Source

```bash
# Clone the repository
git clone https://github.com/yourusername/getRPF.git
cd getRPF

# Install in development mode
pip install -e .

# Or install with all dependencies
pip install -e .[dev,docs]
```

### Using conda/mamba

```bash
# Create environment with required dependencies
conda create -n getrf python=3.10
conda activate getrf

# Install from source
pip install -e .
```

## 🚀 Quick Start

### 1. Basic RPF Extraction

Extract RPFs using built-in architecture database:

```bash
# Extract from FASTQ file
getRPF extract-rpf input_reads.fastq output_rpfs.fastq -f fastq

# Extract with detailed report in CSV format
getRPF extract-rpf input_reads.fastq output_rpfs.fastq -f fastq --output-format csv
```

### 2. Quality Assessment

Check if your data is clean and ready for analysis:

```bash
# Basic cleanliness check with categorized failure reporting
getRPF check-cleanliness sample.fastq -f fastq -o reports/

# Result examples:
# ✅ CLEAN: sample.fastq passed all cleanliness checks  
# ❌ NEEDS_SEQSPEC: Primary failure: end_bias
```

### 3. Novel Protocol Support

Load custom protocols from seqspec files:

```bash
# Create directory with your seqspec files
mkdir my_protocols/
# Add your protocol.yaml files to my_protocols/

# Extract using novel protocols
getRPF extract-rpf input_reads.fastq output_rpfs.fastq -f fastq --seqspec-dir my_protocols/
```



### 4. Architecture Detection & Reporting

Designed for automated processing with "clean/fail" logic, but generates interactive reports for **post-hoc investigation**:

```bash
# Run automated detection
getRPF detect-architecture input.fastq output.fastq --generate-seqspec

# View report for detailed investigation: output.report.html
# Features for post-hoc analysis:
# - Interactive signal plots to verify UMI/Adapter boundaries
# - Decision trace to understand why a specific architecture was inferred
```

## 📋 Commands Overview

getRPF provides five main commands for different aspects of RPF analysis:


| Command | Purpose | Use Case |
|---------|---------|----------|
| `extract-rpf` | **Main extraction workflow** | Core RPF isolation from raw reads |
| `check-cleanliness` | **Quality assessment with categorization** | Scale to thousands of samples |
| `check` | Basic quality analysis | Simple quality metrics |
| `detect-adapter` | Adapter sequence detection | Contamination analysis |
| `align-detect` | STAR alignment + feature detection | Alignment-based quality checks |

## 🔄 User Workflows

### Workflow 1: Single Clean Sample

For samples that pass quality checks:

```bash
# 1. Check sample quality
getRPF check-cleanliness sample.fastq -f fastq -o reports/
# Result: ✅ CLEAN

# 2. Extract RPFs directly
getRPF extract-rpf sample.fastq rpfs.fastq -f fastq --generate-seqspec

# 3. Outputs:
#    - rpfs.fastq (clean RPF sequences)
#    - rpfs.seqspec.yaml (detected architecture)
#    - rpfs.extraction_report.json (statistics)
```

### Workflow 2: Large-Scale Quality Screening (Recommended for 4k+ samples)

Process thousands of samples efficiently:

```bash
# 1. Batch quality screening
for sample in samples/*.fastq; do
    getRPF check-cleanliness "$sample" -f fastq -o quality_reports/
done

# 2. Results automatically categorized:
#    - Clean samples: Ready for extraction
#    - length_distribution_failures/: Need protocol-specific seqspecs  
#    - end_bias_failures/: Need bias-aware seqspecs
#    - low_complexity_failures/: Need adapter contamination removal
```

### Workflow 3: Novel Protocol Integration

When you have custom/novel protocols:

```bash
# 1. Create seqspec files for your protocols
# Example: my_lab_protocol.yaml
cat > my_protocols/custom_protocol.yaml << EOF
assay_id: custom_protocol_2024
name: Custom Lab Protocol
sequence_spec:
  - region_id: umi
    region_type: umi
    sequence: NNNNNN
    min_len: 6
    max_len: 6
  - region_id: rpf
    region_type: cdna
    min_len: 26
    max_len: 34
  - region_id: adapter
    region_type: adapter
    sequence: AGATCGGAAGAG
    min_len: 12
    max_len: 12
EOF

# 2. Use custom protocols for extraction
getRPF extract-rpf input.fastq output.fastq -f fastq --seqspec-dir my_protocols/

# 3. System automatically detects and uses your custom protocol
```

### Workflow 4: Comprehensive Analysis Pipeline

Full analysis with all features:

```bash
#!/bin/bash
# Complete getRPF analysis pipeline

INPUT_FILE="sample.fastq"
OUTPUT_DIR="analysis_results"
mkdir -p "$OUTPUT_DIR"

echo "=== getRPF Analysis Pipeline ==="

# Step 1: Quality Assessment
echo "1. Checking sample cleanliness..."
getRPF check-cleanliness "$INPUT_FILE" -f fastq -o "$OUTPUT_DIR/quality/"

# Step 2: Extract RPFs with seqspec generation  
echo "2. Extracting RPFs..."
getRPF extract-rpf "$INPUT_FILE" "$OUTPUT_DIR/clean_rpfs.fastq" \
    -f fastq --generate-seqspec --output-format json

# Step 3: Adapter Analysis
echo "3. Analyzing adapter contamination..."
getRPF detect-adapter "$INPUT_FILE" -f fastq \
    -a "AGATCGGAAGAGCACACGTCT" -o "$OUTPUT_DIR/adapter_report.txt"

# Step 4: Alignment-based validation (if STAR index available)
if [ -d "/path/to/star/index" ]; then
    echo "4. Alignment validation..."
    getRPF align-detect "$OUTPUT_DIR/clean_rpfs.fastq" \
        -s "/path/to/star/index" -f fastq -o "$OUTPUT_DIR/alignment_report.json"
fi

echo "=== Analysis Complete ==="
echo "Results in: $OUTPUT_DIR/"
```

## 📁 Input/Output Formats

### Supported Input Formats

| Format | Extension | Description | Usage |
|--------|-----------|-------------|--------|
| **FASTQ** | `.fastq`, `.fq` | Standard format with quality scores | Most common, recommended |
| **FASTA** | `.fasta`, `.fa` | Sequences only, no quality | When quality not needed |
| **Collapsed** | `.fasta`, `.fa` | Deduplicated with counts in headers | Memory-efficient for large datasets |

### Input Format Examples

**FASTQ Format:**
```
@read_1
ATGCGATCGCTAGCGATCGCTAGCAGATCGGAAGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Collapsed FASTA Format:**
```
>read_1_count_150
ATGCGATCGCTAGCGATCGCTAGC
>read_2_count_200  
GCGATCGCTAGCGATCGCTAGCA
```

### Output Files

getRPF generates several output files for comprehensive analysis:

```
analysis_results/
├── clean_rpfs.fastq                    # Extracted RPF sequences
├── clean_rpfs.seqspec.yaml             # Architecture specification  
├── clean_rpfs.extraction_report.json   # Detailed extraction statistics
├── quality/
│   └── sample_cleanliness_report.txt   # Quality assessment results
└── adapter_report.txt                  # Adapter contamination analysis
```

## 🔧 Advanced Usage

### Custom Architecture Databases

Load your own architecture definitions:

```bash
# Using JSON architecture database
getRPF extract-rpf input.fastq output.fastq -f fastq --architecture-db my_architectures.json

# Using seqspec directory (recommended)
getRPF extract-rpf input.fastq output.fastq -f fastq --seqspec-dir protocols/
```

### Performance Optimization

For large datasets:

```bash
# Limit reads for testing
getRPF extract-rpf input.fastq output.fastq -f fastq --max-reads 10000

# Process collapsed format for memory efficiency  
getRPF extract-rpf collapsed.fasta output.fastq -f collapsed

# Parallel processing example
parallel -j 8 getRPF extract-rpf {} {.}_rpfs.fastq -f fastq ::: samples/*.fastq
```

### Integration with Other Tools

getRPF works seamlessly with standard bioinformatics tools:

```bash
# Combine with cutadapt for pre-processing
cutadapt -a AGATCGGAAGAG -o preprocessed.fastq input.fastq
getRPF extract-rpf preprocessed.fastq rpfs.fastq -f fastq

# Chain with downstream analysis
getRPF extract-rpf input.fastq rpfs.fastq -f fastq
bowtie2 -x genome_index -U rpfs.fastq -S aligned.sam

# Quality control with FastQC
getRPF extract-rpf input.fastq rpfs.fastq -f fastq
fastqc rpfs.fastq
```

## 🔍 Quality Assessment Details

getRPF implements comprehensive quality checks to ensure clean RPF data:

### Cleanliness Criteria (All Must Pass)

1. **Length Distribution**: Reads within expected RPF size range (26-35 nt)
2. **Information Content**: Uniform complexity across read positions (Shannon entropy)
3. **End Bias**: No excessive nucleotide bias at 5' or 3' ends  
4. **Base Composition**: No extreme positional nucleotide bias
5. **GC Content**: Within expected range for ribosomal sequences

### Failure Categorization

Failed samples are automatically categorized for efficient batch processing:

- **`length_distribution_failures/`**: Wrong size distribution → Find protocol-specific seqspecs
- **`end_bias_failures/`**: 5'/3' nucleotide bias → Protocol has adapter sequences
- **`low_complexity_failures/`**: Repetitive sequences → Contamination or degradation
- **`base_composition_failures/`**: Positional bias → Systematic sequencing artifacts

## 🧬 seqspec Integration

getRPF supports the standard seqspec format for unlimited protocol extensibility:

### Creating seqspec Files

```yaml
# example_protocol.yaml
assay_id: lab_protocol_2024
name: "Custom Ribosome Profiling Protocol"
description: "Lab-specific protocol with unique barcode structure"
seqspec_version: 0.3.0

sequence_spec:
  - region_id: umi
    region_type: umi
    sequence: NNNNNN
    min_len: 6
    max_len: 6
    
  - region_id: sample_barcode
    region_type: barcode  
    sequence: NNNNNNNN
    min_len: 8
    max_len: 8
    
  - region_id: ribosome_protected_fragment
    region_type: cdna
    min_len: 26
    max_len: 34
    
  - region_id: adapter
    region_type: adapter
    sequence: AGATCGGAAGAGCACACGTCT
    min_len: 21
    max_len: 21
```

### Using seqspec Protocols

```bash
# Load all protocols from directory
getRPF extract-rpf input.fastq output.fastq -f fastq --seqspec-dir my_protocols/

# System automatically:
# 1. Loads all .yaml files in directory
# 2. Creates detection profiles for each protocol  
# 3. Matches input reads against all protocols
# 4. Uses best-matching protocol for extraction
```

## 🐛 Troubleshooting

### Common Issues

**Issue**: `pysam not available - soft-clipping analysis will be disabled`
```bash
# Solution: Install pysam
pip install pysam>=0.19.0
```

**Issue**: `No architecture detected` 
```bash
# Solutions:
# 1. Check if your protocol is in built-in database
getRPF extract-rpf --help

# 2. Create seqspec file for your protocol
getRPF extract-rpf input.fastq output.fastq -f fastq --seqspec-dir my_protocols/

# 3. Use de novo detection (automatic fallback)
```

**Issue**: Sample fails cleanliness checks
```bash
# Check specific failure type
getRPF check-cleanliness sample.fastq -f fastq -o reports/
cat reports/sample_cleanliness_report.txt

# Address by failure category:
# - length_distribution: Need protocol-specific extraction
# - end_bias: Adapter contamination, create seqspec
# - low_complexity: Pre-filter with cutadapt
```

### Performance Issues

**Large files taking too long:**
```bash  
# Use collapsed format to reduce memory usage
# Test with limited reads first
getRPF extract-rpf input.fastq output.fastq -f fastq --max-reads 10000

# Process in parallel
parallel -j 4 getRPF extract-rpf {} {.}_rpfs.fastq -f fastq ::: *.fastq
```

**Memory usage too high:**
```bash
# Convert to collapsed format first
# Use streaming processing for large datasets
getRPF extract-rpf input.fasta output.fastq -f collapsed
```

## 📊 Output Interpretation

### Extraction Report Example

```json
{
  "extraction_statistics": {
    "input_reads": 50000,
    "extracted_rpfs": 45230,
    "extraction_rate": 0.9046,
    "mean_rpf_length": 28.4
  },
  "architecture_detection": {
    "method": "pattern_matching",
    "matched_architecture": "mcglincy_2017",
    "confidence": 0.95,
    "lab_source": "McGlincy Lab"
  },
  "quality_metrics": {
    "length_distribution_valid": true,
    "adapter_contamination_rate": 0.02,
    "mean_quality_score": 35.2
  }
}
```

### Cleanliness Report Example

```
=== RPF Data Check Results ===

Sample Status: CLEAN

Check Summary:
length_distribution            [PASS]
information_content            [PASS]  
end_bias                       [PASS]
base_composition               [PASS]
gc_content                     [PASS]
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/yourusername/getRPF.git
cd getRPF
pip install -e .[dev]

# Run tests
pytest tests/

# Code formatting
black src/ tests/
isort src/ tests/

# Type checking  
mypy src/
```

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use getRPF in your research, please cite:

```bibtex
@software{getRPF,
  title={getRPF: Automated Ribosome Protected Fragment Analysis},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/getRPF}
}
```

## 🔗 Links

- **Documentation**: [Full documentation](https://getRPF.readthedocs.io)
- **Bug Reports**: [GitHub Issues](https://github.com/yourusername/getRPF/issues)
- **seqspec Format**: [Lior Pachter Lab](https://pachterlab.github.io/seqspec/)
- **Ribo-seq Resources**: [Ribo-seq.org](https://ribo-seq.org)

---

**getRPF** - Making ribosome profiling analysis accessible, scalable, and extensible. 🧬✨
