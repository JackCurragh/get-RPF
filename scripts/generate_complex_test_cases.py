#!/usr/bin/env python3
"""Generate complex ribosome profiling test cases with multiple structural elements.

Creates realistic test datasets representing modern ribosome profiling protocols
with UMIs, spacers, adapters, barcodes, untemplated additions, and variable regions.
"""

import random
import sys
from pathlib import Path
from typing import List, Dict, Tuple

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def generate_rpf_sequence(length: int) -> str:
    """Generate realistic RPF with codon bias."""
    # More realistic biological composition with codon bias
    codons = ['ATG', 'TGA', 'TAA', 'TAG', 'AAA', 'GAA', 'GGA', 'CGA', 'AGA', 'TTA']
    weights = [0.15, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2]  # Start codon bias
    
    sequence = ""
    while len(sequence) < length:
        codon = random.choices(codons, weights=weights, k=1)[0]
        sequence += codon
    
    return sequence[:length]

def add_untemplated_nucleotides(sequence: str, position: str = "3prime", length: int = None) -> str:
    """Add untemplated nucleotides (poly-A tails, random additions)."""
    if length is None:
        length = random.randint(1, 5)
    
    untemplated_types = {
        'poly_a': 'A' * length,
        'poly_t': 'T' * length,
        'random': ''.join(random.choices(['A', 'T', 'G', 'C'], k=length)),
        'gc_rich': ''.join(random.choices(['G', 'C'], k=length)),
    }
    
    addition = random.choice(list(untemplated_types.values()))
    
    if position == "5prime":
        return addition + sequence
    else:  # 3prime
        return sequence + addition

def add_spacer_sequence(sequence: str, spacer_type: str = "random") -> str:
    """Add spacer sequences between functional elements."""
    spacer_sequences = {
        'poly_t': 'TTTT',
        'poly_a': 'AAAA', 
        'linker': 'TGACT',
        'random': ''.join(random.choices(['A', 'T', 'G', 'C'], k=random.randint(3, 6))),
        'structured': random.choice(['TATA', 'GAGA', 'CTCT', 'GCGC']),
    }
    
    spacer = spacer_sequences.get(spacer_type, spacer_sequences['random'])
    return sequence + spacer

def generate_complex_test_case(case_name: str, num_reads: int = 1000) -> Dict[str, any]:
    """Generate complex ribosome profiling test cases."""
    
    reads = []
    metadata = {'case_name': case_name, 'num_reads': num_reads}
    
    if case_name == 'modern_protocol_v1':
        # Modern protocol: 5'UMI + spacer + barcode + adapter + RPF + 3'adapter + untemplated
        metadata.update({
            'description': '5\'-UMI-spacer-barcode-adapter-RPF-adapter-untemplated-3\'',
            'structure': 'UMI(8) + spacer(4) + barcode(6) + adapter(12) + RPF(28-32) + adapter(10) + polyA(1-4)',
            'expected_success_rate': 1.0,
            'expected_method': 'pattern_match'
        })
        
        for i in range(num_reads):
            # Build complex read structure
            sequence_parts = []
            
            # 1. 5' UMI (8nt random)
            umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=8))
            sequence_parts.append(('umi', umi))
            
            # 2. Spacer (4nt structured)
            spacer1 = random.choice(['TGAC', 'CTAG', 'GCAT', 'ATGC'])
            sequence_parts.append(('spacer', spacer1))
            
            # 3. Sample barcode (6nt semi-random from list)
            barcodes = ['ATCGTA', 'TGCAGT', 'AGCTGA', 'CGATCG', 'GTACGT', 'TAGCTG']
            barcode = random.choice(barcodes)
            sequence_parts.append(('barcode', barcode))
            
            # 4. 5' Adapter (12nt fixed)
            adapter5 = 'AGATCGGAAGAG'
            sequence_parts.append(('adapter_5', adapter5))
            
            # 5. RPF sequence (28-32nt biological)
            rpf_length = random.randint(28, 32)
            rpf_seq = generate_rpf_sequence(rpf_length)
            sequence_parts.append(('rpf', rpf_seq))
            
            # 6. 3' Adapter (10nt fixed)
            adapter3 = 'CACGTGCTAG'
            sequence_parts.append(('adapter_3', adapter3))
            
            # 7. Untemplated addition (1-4nt poly-A)
            untemplated_length = random.randint(1, 4)
            untemplated = 'A' * untemplated_length
            sequence_parts.append(('untemplated', untemplated))
            
            # Combine all parts
            full_sequence = ''.join(part[1] for part in sequence_parts)
            quality = 'I' * len(full_sequence)
            
            # Store structure info in header
            structure_info = '_'.join([f"{name}:{len(seq)}" for name, seq in sequence_parts])
            header = f"@complex_read_{i}_{structure_info}"
            
            reads.append((header, full_sequence, quality, sequence_parts))
    
    elif case_name == 'variable_structure':
        # Variable structure reads - different combinations per read
        metadata.update({
            'description': 'Variable structure: some reads have different element combinations',
            'structure': 'Variable combinations of UMI, spacers, barcodes, adapters',
            'expected_success_rate': 0.85,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            sequence_parts = []
            
            # Randomly include/exclude elements
            if random.random() < 0.8:  # 80% have UMI
                umi_length = random.choice([6, 8, 10])
                umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=umi_length))
                sequence_parts.append(('umi', umi))
                
                if random.random() < 0.6:  # 60% have spacer after UMI
                    spacer = random.choice(['TGAC', 'CTAG', 'TTTT'])
                    sequence_parts.append(('spacer', spacer))
            
            if random.random() < 0.7:  # 70% have barcode
                barcode_length = random.choice([4, 6, 8])
                barcode = ''.join(random.choices(['A', 'T', 'G', 'C'], k=barcode_length))
                sequence_parts.append(('barcode', barcode))
            
            if random.random() < 0.9:  # 90% have some adapter
                adapter_types = ['AGATCGGAAGAG', 'CTGTCTCTTATACACATCT', 'TGGAATTCTCG']
                adapter = random.choice(adapter_types)
                sequence_parts.append(('adapter', adapter))
            
            # Always have RPF
            rpf_length = random.randint(25, 35)
            rpf_seq = generate_rpf_sequence(rpf_length)
            sequence_parts.append(('rpf', rpf_seq))
            
            if random.random() < 0.4:  # 40% have untemplated addition
                untemplated_length = random.randint(1, 6)
                untemplated_type = random.choice(['poly_a', 'poly_t', 'random'])
                if untemplated_type == 'poly_a':
                    untemplated = 'A' * untemplated_length
                elif untemplated_type == 'poly_t':
                    untemplated = 'T' * untemplated_length
                else:
                    untemplated = ''.join(random.choices(['A', 'T', 'G', 'C'], k=untemplated_length))
                sequence_parts.append(('untemplated', untemplated))
            
            # Combine parts
            full_sequence = ''.join(part[1] for part in sequence_parts)
            quality = 'I' * len(full_sequence)
            
            structure_info = '_'.join([f"{name}:{len(seq)}" for name, seq in sequence_parts])
            header = f"@variable_read_{i}_{structure_info}"
            
            reads.append((header, full_sequence, quality, sequence_parts))
    
    elif case_name == 'paired_end_complex':
        # Paired-end with different structures in R1 vs R2
        metadata.update({
            'description': 'Paired-end: R1 has UMI+barcode, R2 has RPF+adapter',
            'structure': 'R1: UMI(10) + barcode(8) + spacer(4); R2: RPF(30) + adapter(15)',
            'expected_success_rate': 1.0,
            'expected_method': 'pattern_match',
            'paired_end': True
        })
        
        for i in range(num_reads):
            # R1: Index read with UMI + barcode
            r1_parts = []
            
            # UMI (10nt random)
            umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=10))
            r1_parts.append(('umi', umi))
            
            # Sample barcode (8nt from list)
            barcodes = ['ATCGTAAC', 'TGCAGTGA', 'AGCTGACT', 'CGATCGTT']
            barcode = random.choice(barcodes)
            r1_parts.append(('barcode', barcode))
            
            # Spacer
            spacer = 'TGAC'
            r1_parts.append(('spacer', spacer))
            
            r1_sequence = ''.join(part[1] for part in r1_parts)
            r1_quality = 'I' * len(r1_sequence)
            
            # R2: Biological read with RPF
            r2_parts = []
            
            # RPF sequence (30nt)
            rpf_seq = generate_rpf_sequence(30)
            r2_parts.append(('rpf', rpf_seq))
            
            # 3' adapter (15nt)
            adapter = 'AGATCGGAAGAGCAC'
            r2_parts.append(('adapter', adapter))
            
            r2_sequence = ''.join(part[1] for part in r2_parts)
            r2_quality = 'I' * len(r2_sequence)
            
            # Create paired reads
            r1_structure = '_'.join([f"{name}:{len(seq)}" for name, seq in r1_parts])
            r2_structure = '_'.join([f"{name}:{len(seq)}" for name, seq in r2_parts])
            
            r1_header = f"@paired_read_{i}_R1_{r1_structure}"
            r2_header = f"@paired_read_{i}_R2_{r2_structure}"
            
            # Store both reads with their structure info
            reads.append((r1_header, r1_sequence, r1_quality, r1_parts))
            reads.append((r2_header, r2_sequence, r2_quality, r2_parts))
    
    elif case_name == 'degraded_complex':
        # Complex reads with degradation, incomplete elements
        metadata.update({
            'description': 'Degraded complex reads: partial adapters, truncated elements',
            'structure': 'Degraded: partial UMI + incomplete adapter + RPF + random tail',
            'expected_success_rate': 0.7,
            'expected_method': 'de_novo'
        })
        
        for i in range(num_reads):
            sequence_parts = []
            
            # Partial UMI (should be 8nt, but degraded to 3-6nt)
            expected_umi_length = 8
            actual_umi_length = random.randint(3, 6)
            umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=actual_umi_length))
            sequence_parts.append(('umi_partial', umi))
            
            # Incomplete adapter (should be AGATCGGAAGAG, but truncated)
            full_adapter = 'AGATCGGAAGAG'
            truncate_length = random.randint(5, len(full_adapter))
            partial_adapter = full_adapter[:truncate_length]
            sequence_parts.append(('adapter_partial', partial_adapter))
            
            # RPF (may be shortened due to degradation)
            rpf_length = random.randint(20, 32)  # Shorter than ideal
            rpf_seq = generate_rpf_sequence(rpf_length)
            sequence_parts.append(('rpf', rpf_seq))
            
            # Random tail (degradation artifact)
            if random.random() < 0.6:  # 60% have random tail
                tail_length = random.randint(2, 8)
                tail = ''.join(random.choices(['A', 'T', 'G', 'C'], k=tail_length))
                sequence_parts.append(('random_tail', tail))
            
            # Combine parts
            full_sequence = ''.join(part[1] for part in sequence_parts)
            
            # Degraded quality (lower toward ends)
            quality_chars = []
            for pos in range(len(full_sequence)):
                base_qual = max(35 - (pos // 3), 20)  # Degrade with position
                qual_char = chr(base_qual + 33)
                quality_chars.append(qual_char)
            quality = ''.join(quality_chars)
            
            structure_info = '_'.join([f"{name}:{len(seq)}" for name, seq in sequence_parts])
            header = f"@degraded_read_{i}_{structure_info}"
            
            reads.append((header, full_sequence, quality, sequence_parts))
    
    elif case_name == 'ultra_complex':
        # Ultra-complex: all possible elements in single reads
        metadata.update({
            'description': 'Ultra-complex: UMI + spacer1 + barcode + spacer2 + adapter1 + RPF + adapter2 + untemplated',
            'structure': 'All elements: 8 different functional regions per read',
            'expected_success_rate': 0.9,
            'expected_method': 'pattern_match'
        })
        
        for i in range(num_reads):
            sequence_parts = []
            
            # 1. 5' UMI
            umi = ''.join(random.choices(['A', 'T', 'G', 'C'], k=8))
            sequence_parts.append(('umi', umi))
            
            # 2. First spacer
            spacer1 = 'TGAC'
            sequence_parts.append(('spacer1', spacer1))
            
            # 3. Sample barcode
            barcode = random.choice(['ATCGTA', 'TGCAGT', 'AGCTGA', 'CGATCG'])
            sequence_parts.append(('barcode', barcode))
            
            # 4. Second spacer
            spacer2 = 'CTAG'
            sequence_parts.append(('spacer2', spacer2))
            
            # 5. 5' adapter
            adapter1 = 'AGATCGGAAGAG'
            sequence_parts.append(('adapter1', adapter1))
            
            # 6. RPF sequence
            rpf_seq = generate_rpf_sequence(30)
            sequence_parts.append(('rpf', rpf_seq))
            
            # 7. 3' adapter
            adapter2 = 'CACGTGCTAG'
            sequence_parts.append(('adapter2', adapter2))
            
            # 8. Untemplated addition
            untemplated_length = random.randint(2, 5)
            untemplated = random.choice(['A' * untemplated_length, 
                                      'T' * untemplated_length,
                                      ''.join(random.choices(['A', 'T'], k=untemplated_length))])
            sequence_parts.append(('untemplated', untemplated))
            
            # Combine all parts
            full_sequence = ''.join(part[1] for part in sequence_parts)
            quality = 'I' * len(full_sequence)
            
            structure_info = '_'.join([f"{name}:{len(seq)}" for name, seq in sequence_parts])
            header = f"@ultra_complex_read_{i}_{structure_info}"
            
            reads.append((header, full_sequence, quality, sequence_parts))
    
    else:
        raise ValueError(f"Unknown complex test case: {case_name}")
    
    metadata['actual_num_reads'] = len(reads)
    return {'reads': reads, 'metadata': metadata}

def write_complex_fastq(reads: List[Tuple], output_file: Path):
    """Write complex reads to FASTQ format with structure annotations."""
    with open(output_file, 'w') as f:
        for read_info in reads:
            if len(read_info) >= 4:  # Has structure info
                header, sequence, quality, structure = read_info[:4]
            else:  # Simple format
                header, sequence, quality = read_info[:3]
            
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")

def main():
    """Generate complex ribosome profiling test cases."""
    
    complex_test_cases = [
        'modern_protocol_v1',
        'variable_structure', 
        'paired_end_complex',
        'degraded_complex',
        'ultra_complex'
    ]
    
    output_dir = Path('complex_test_data')
    output_dir.mkdir(exist_ok=True)
    
    all_metadata = {}
    
    print("🧬 Generating COMPLEX ribosome profiling test cases...")
    
    for case_name in complex_test_cases:
        print(f"  Creating {case_name}...")
        
        test_data = generate_complex_test_case(case_name, num_reads=500)  # Smaller for complexity
        reads = test_data['reads']
        metadata = test_data['metadata']
        
        # Write files
        output_file = output_dir / f"{case_name}.fastq"
        write_complex_fastq(reads, output_file)
        
        metadata['file_path'] = str(output_file)
        all_metadata[case_name] = metadata
        
        print(f"    → {output_file} ({metadata['actual_num_reads']} reads)")
        print(f"      Structure: {metadata['structure']}")
    
    # Write metadata summary
    import json
    metadata_file = output_dir / "complex_test_cases_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(all_metadata, f, indent=2)
    
    print(f"\n✅ Generated {len(complex_test_cases)} COMPLEX test cases in {output_dir}/")
    print(f"📋 Metadata saved to {metadata_file}")
    
    # Print summary
    print("\n📊 Complex Test Case Summary:")
    for case_name, meta in all_metadata.items():
        expected_rate = meta['expected_success_rate'] * 100
        structure = meta['structure'][:60] + "..." if len(meta['structure']) > 60 else meta['structure']
        print(f"  {case_name:20s}: {meta['actual_num_reads']:4d} reads, {expected_rate:3.0f}% expected")
        print(f"    {structure}")

if __name__ == "__main__":
    random.seed(42)  # Reproducible results
    main()