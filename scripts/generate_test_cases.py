#!/usr/bin/env python3
"""Generate comprehensive test cases for RPF extraction validation.

Creates diverse test datasets with various contaminating factors to validate
the robustness of the RPF extraction system.
"""

import random
import sys
from pathlib import Path
from typing import List, Dict, Tuple

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def generate_rpf_sequence(length: int) -> str:
    """Generate a realistic RPF sequence with biological nucleotide composition."""
    # Biological composition: slightly GC-rich for ribosome binding sites
    nucleotides = ['A', 'T', 'G', 'C']
    weights = [0.25, 0.25, 0.27, 0.23]  # Slightly GC-rich
    
    return ''.join(random.choices(nucleotides, weights=weights, k=length))

def add_umi(sequence: str, umi_length: int = 6) -> str:
    """Add random UMI to 5' end."""
    umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=umi_length))
    return umi + sequence

def add_barcode(sequence: str, barcode_length: int = 5) -> str:
    """Add sample barcode."""
    # Use realistic barcode with some structure
    barcodes = ['ATCGTA', 'TGCAGT', 'AGCTGA', 'CGATCG', 'GTACGT']
    barcode = random.choice(barcodes)[:barcode_length]
    return barcode + sequence

def add_adapter_contamination(sequence: str, adapter_type: str = 'illumina') -> str:
    """Add adapter contamination to 3' end."""
    adapters = {
        'illumina': 'AGATCGGAAGAGCAC',
        'nextera': 'CTGTCTCTTATACACATCT',
        'truseq': 'AGATCGGAAGAGCGTCGTG',
        'custom': 'TGGAATTCTCGGGTGCCAAGG'
    }
    
    adapter = adapters.get(adapter_type, adapters['illumina'])
    
    # Add partial or full adapter
    if random.random() < 0.3:  # 30% chance of partial adapter
        adapter_length = random.randint(5, len(adapter))
        adapter = adapter[:adapter_length]
    
    return sequence + adapter

def add_quality_scores(sequence: str, quality_pattern: str = 'high') -> str:
    """Generate quality scores for FASTQ format."""
    length = len(sequence)
    
    if quality_pattern == 'high':
        # High quality (mostly I = 40)
        qualities = ['I'] * length
    elif quality_pattern == 'mixed':
        # Mixed quality with some low regions
        qualities = []
        for i in range(length):
            if random.random() < 0.1:  # 10% low quality
                qualities.append(random.choice(['#', '$', '%', '&']))
            else:
                qualities.append(random.choice(['H', 'I', 'J', 'K']))
    elif quality_pattern == 'degraded':
        # Degraded quality toward 3' end
        qualities = []
        for i in range(length):
            base_qual = max(35 - (i * 2), 20)  # Decrease with position
            qual_char = chr(base_qual + 33)
            qualities.append(qual_char)
    else:
        qualities = ['I'] * length
    
    return ''.join(qualities)

def generate_test_case(case_name: str, num_reads: int = 1000) -> Dict[str, any]:
    """Generate a specific test case scenario."""
    
    reads = []
    metadata = {'case_name': case_name, 'num_reads': num_reads}
    
    if case_name == 'clean_rpfs':
        # Perfect RPF data - no contamination
        metadata.update({
            'description': 'Clean RPF sequences, no UMI/barcode/adapter',
            'expected_success_rate': 1.0,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            rpf_length = random.randint(26, 34)
            sequence = generate_rpf_sequence(rpf_length)
            quality = add_quality_scores(sequence, 'high')
            reads.append((f"@clean_read_{i}", sequence, quality))
    
    elif case_name == 'mcglincy_protocol':
        # McGlincy & Ingolia protocol: UMI + barcode + RPF
        metadata.update({
            'description': 'McGlincy protocol: 5nt UMI + 5nt barcode + RPF',
            'expected_success_rate': 1.0,
            'expected_method': 'pattern_match',
            'expected_architecture': 'mcglincy_ingolia_2017'
        })
        
        for i in range(num_reads):
            rpf_length = random.randint(28, 32)
            rpf_seq = generate_rpf_sequence(rpf_length)
            
            # Add UMI + barcode
            full_sequence = add_umi(rpf_seq, 5)
            full_sequence = add_barcode(full_sequence, 5)
            
            quality = add_quality_scores(full_sequence, 'high')
            reads.append((f"@mcglincy_read_{i}", full_sequence, quality))
    
    elif case_name == 'adapter_contaminated':
        # Mixed clean + adapter contaminated (like our test data)
        metadata.update({
            'description': '75% clean RPFs + 25% with Illumina adapter contamination',
            'expected_success_rate': 1.0,
            'expected_method': 'mixed_strategy'
        })
        
        for i in range(num_reads):
            rpf_length = random.randint(28, 35)
            sequence = generate_rpf_sequence(rpf_length)
            
            if random.random() < 0.25:  # 25% contaminated
                sequence = add_adapter_contamination(sequence, 'illumina')
            
            quality = add_quality_scores(sequence, 'mixed')
            reads.append((f"@mixed_read_{i}", sequence, quality))
    
    elif case_name == 'multi_adapter_types':
        # Multiple different adapter types mixed together
        metadata.update({
            'description': 'Mix of different adapter types (Illumina, Nextera, TruSeq)',
            'expected_success_rate': 0.9,  # Some adapters might not be recognized
            'expected_method': 'mixed_strategy'
        })
        
        adapter_types = ['illumina', 'nextera', 'truseq', 'custom']
        
        for i in range(num_reads):
            rpf_length = random.randint(26, 36)
            sequence = generate_rpf_sequence(rpf_length)
            
            if random.random() < 0.4:  # 40% contaminated
                adapter_type = random.choice(adapter_types)
                sequence = add_adapter_contamination(sequence, adapter_type)
            
            quality = add_quality_scores(sequence, 'mixed')
            reads.append((f"@multi_adapter_read_{i}", sequence, quality))
    
    elif case_name == 'length_variants':
        # RPFs with unusual length distribution
        metadata.update({
            'description': 'RPFs with wide length distribution (20-50nt)',
            'expected_success_rate': 0.85,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            # Wider length distribution
            rpf_length = random.choice([20, 21, 22] + list(range(25, 38)) + [42, 45, 50])
            sequence = generate_rpf_sequence(rpf_length)
            quality = add_quality_scores(sequence, 'high')
            reads.append((f"@length_variant_read_{i}", sequence, quality))
    
    elif case_name == 'low_quality':
        # Good sequences but with poor quality scores
        metadata.update({
            'description': 'Good RPF sequences with degraded quality scores',
            'expected_success_rate': 0.95,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            rpf_length = random.randint(28, 34)
            sequence = generate_rpf_sequence(rpf_length)
            quality = add_quality_scores(sequence, 'degraded')
            reads.append((f"@low_qual_read_{i}", sequence, quality))
    
    elif case_name == 'extreme_contamination':
        # Heavy contamination with multiple issues
        metadata.update({
            'description': 'Heavy contamination: UMI + adapter + quality issues',
            'expected_success_rate': 0.7,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            rpf_length = random.randint(26, 32)
            sequence = generate_rpf_sequence(rpf_length)
            
            # Add multiple contaminants
            if random.random() < 0.6:  # 60% have UMI
                sequence = add_umi(sequence, random.randint(4, 8))
            
            if random.random() < 0.5:  # 50% have adapter
                sequence = add_adapter_contamination(sequence, random.choice(['illumina', 'truseq']))
            
            quality = add_quality_scores(sequence, random.choice(['mixed', 'degraded']))
            reads.append((f"@extreme_contam_read_{i}", sequence, quality))
    
    elif case_name == 'collapsed_format':
        # Collapsed FASTA with count annotations
        metadata.update({
            'description': 'Collapsed FASTA format with _x count annotations',
            'expected_success_rate': 1.0,
            'expected_method': 'de_novo',
            'format': 'collapsed'
        })
        
        unique_sequences = {}
        # Generate fewer unique sequences with counts
        for i in range(num_reads // 4):  # 250 unique sequences
            rpf_length = random.randint(28, 34)
            sequence = generate_rpf_sequence(rpf_length)
            count = random.randint(1, 20)
            unique_sequences[f"seq{i}_x{count}"] = sequence
        
        # Convert to read format for consistency
        reads = [(f">{header}", sequence, "") for header, sequence in unique_sequences.items()]
        metadata['num_reads'] = len(reads)
    
    else:
        raise ValueError(f"Unknown test case: {case_name}")
    
    metadata['actual_num_reads'] = len(reads)
    return {'reads': reads, 'metadata': metadata}

def write_fastq(reads: List[Tuple[str, str, str]], output_file: Path):
    """Write reads to FASTQ format."""
    with open(output_file, 'w') as f:
        for header, sequence, quality in reads:
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")

def write_fasta(reads: List[Tuple[str, str, str]], output_file: Path):
    """Write reads to FASTA format."""
    with open(output_file, 'w') as f:
        for header, sequence, _ in reads:
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")

def main():
    """Generate all test cases."""
    
    test_cases = [
        'clean_rpfs',
        'mcglincy_protocol', 
        'adapter_contaminated',
        'multi_adapter_types',
        'length_variants',
        'low_quality',
        'extreme_contamination',
        'collapsed_format'
    ]
    
    output_dir = Path('comprehensive_test_data')
    output_dir.mkdir(exist_ok=True)
    
    # Generate metadata file
    all_metadata = {}
    
    print("🧪 Generating comprehensive test cases...")
    
    for case_name in test_cases:
        print(f"  Creating {case_name}...")
        
        test_data = generate_test_case(case_name, num_reads=1000)
        reads = test_data['reads']
        metadata = test_data['metadata']
        
        # Write files
        if metadata.get('format') == 'collapsed':
            output_file = output_dir / f"{case_name}.fasta"
            write_fasta(reads, output_file)
        else:
            output_file = output_dir / f"{case_name}.fastq"
            write_fastq(reads, output_file)
        
        metadata['file_path'] = str(output_file)
        all_metadata[case_name] = metadata
        
        print(f"    → {output_file} ({metadata['actual_num_reads']} reads)")
    
    # Write metadata summary
    import json
    metadata_file = output_dir / "test_cases_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(all_metadata, f, indent=2)
    
    print(f"\n✅ Generated {len(test_cases)} test cases in {output_dir}/")
    print(f"📋 Metadata saved to {metadata_file}")
    
    # Print summary
    print("\n📊 Test Case Summary:")
    for case_name, meta in all_metadata.items():
        expected_rate = meta['expected_success_rate'] * 100
        print(f"  {case_name:20s}: {meta['actual_num_reads']:4d} reads, expected {expected_rate:4.0f}% success")

if __name__ == "__main__":
    random.seed(42)  # Reproducible results
    main()