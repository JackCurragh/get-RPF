# RPF-Specific Checks

getRPF performs several specific checks designed for Ribosome Protected Fragment (RPF) data.

## Length Distribution Check

RPFs typically have a characteristic length distribution due to the physical size of the ribosome:

- Expected range: 26-35 nucleotides
- Typical mode: ~28-30 nucleotides

The check will:
- PASS if mode is within range and distribution is clean
- WARN if mode is outside range but distribution is clean
- FAIL if no clear mode or most reads outside range

## Base Composition Check

Checks for concerning patterns in nucleotide composition:

- FAIL if any position shows extreme base bias (>85% single base)
- Reports position-specific nucleotide frequencies
- Identifies potential systematic biases

## GC Content Check

Monitors overall GC content:

- Expected range: 35-65%
- WARN if outside expected range
- Helps identify potential contamination

## Understanding Reports

The check command generates two reports:

1. Basic Statistics Report (`report.txt`):
   - Read length distribution
   - Nucleotide frequencies
   - Quality scores (FASTQ only)
   - GC content

2. RPF Checks Report (`report.rpf_checks.txt`):
   - Overall PASS/WARN/FAIL status
   - Detailed check results
   - Specific issues identified
   - Recommendations if problems found

Example Report:
```
=== RPF Data Check Results ===

Summary:
Length Distribution            [PASS]
Base Composition              [PASS]
GC Content                    [PASS]

Detailed Results:
Length Distribution:
Status: PASS
Message: Read length distribution matches RPF expectations
Details:
  mode_length: 30
  fraction_in_range: 0.95
  total_reads: 1000
```
