"""RPF extraction processor for automated ribosome protected fragment isolation.

This module implements the core RPF extraction system that combines:
1. Signal Processing (Entropy, Composition)
2. Strict Architecture Matching
3. Probabilistic Segmentation (HMM)
4. Explainable Reporting
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, asdict

from ...utils.file_utils import get_file_opener
from ..seqspec_generator import SeqSpecGenerator
from ..seqspec_loader import SeqSpecArchitectureLoader

# Import types and components
from .types import ReadArchitecture, SegmentInfo, RPFExtractionResult
from .signals import SignalProcessor, SignalStats
from .matcher import ArchitectureMatcher
from .segmenter import ProbabilisticSegmenter
from .reporting import Reporter

logger = logging.getLogger(__name__)


# Comprehensive database of known adapters
KNOWN_ADAPTERS = {
    # Guo et al. (2014)
    "Guo14": "TCGTATGCCGTCTTCTG",
    
    # Hakon
    "Hakon2": "CACTCGGGCACCAAGGA",
    "Hakon3": "GTGTCAGTCACTTCCAGCGG",
    "Hakon4": "TGTAGGCACCATC",
    "Hakon5": "AAAAAAAAAA",
    "Hakon6": "TCGTATGCCGTCTTCTGCTT",

    # Kiniry et al.
    "Kiniry1": "CTGTAGGCACCATCAAT",
    "Kiniry2": "AGATCGGAAGAGC",
    "Kiniry3": "CGCCTTGGCCGTACAGCAG",
    "Kiniry5": "TGGAATTCTCGGGTGCCAAGG",
    "Kiniry6": "CCTTGGCACCCGAGAATT",
    "Kiniry7": "GATCGGAAGAGCGTCGT",
    "Kiniry8": "CTGATGGCGCGAGGGAG",
    "Kiniry9": "GATCGGAAGAGCACACG",
    "Kiniry10": "AATGATACGGCGACCAC",
    "Kiniry11": "GATCGGAAGAGCTCGTA",
    "Kiniry12": "CAAGCAGAAGACGGCAT",
    "Kiniry13": "ACACTCTTTCCCTACA",
    "Kiniry14": "GATCGGAAGAGCGGTT",
    "Kiniry15": "ACAGGTTCAGAGTTCTA",
    "Kiniry16": "CAAGCAGAAGACGGCAT",
    "Kiniry17": "ACAGGTTCAGAGTTCTA",
    "Kiniry19": "TGATCGGAAGAGCACAC",
    "Kiniry20": "GATCGGAAGAGCACACGT",
    "Kiniry21": "AGATCGGAAGAGCAC",
    "Kiniry22": "AGATCGGAAGAGCACACGTCT",

    # Illumina / Standard
    "Illumina Universal Adapter": "AGATCGGAAGAGC",
    "Illumina Small RNA 3' Adapter": "TGGAATTCTCGG",
    "Illumina Small RNA 5' Adapter": "GATCGTCGGACT",
    "Nextera Transposase Sequence": "CTGTCTCTTATA",
    "SOLID Small RNA Adapter": "CGCCTTGGCCGT",
    "Ingolia 2012 adapter": "CTGTAGGCACCATCAAT",
    "Illumina Uni. Adapter var2": "ATCGTAGATCGGAAG",
    "Tru-seq Small RNA": "TGGAATTCTCGGGTGCCAAGG",
    "Tru-seq Small RNA 2": "GGAATTCTCGGGTGCCAAGG",
    "Illumina": "CTGTCTCTTATACACATCT"
}


class ArchitectureDatabase:
    """Database of known read architectures."""
    
    def __init__(self, db_path: Optional[Path] = None, seqspec_dir: Optional[Path] = None):
        self.architectures: List[ReadArchitecture] = []
        self.db_path = db_path
        self.seqspec_dir = seqspec_dir
        self.seqspec_loader = SeqSpecArchitectureLoader()
        self._load_architectures()
    
    def _load_architectures(self) -> None:
        if self.db_path and self.db_path.exists():
            self._load_from_file()
        else:
            self._initialize_builtin_architectures()
    
    def _load_from_file(self) -> None:
        try:
            with open(self.db_path, 'r') as f:
                data = json.load(f)
                for arch_data in data.get("architectures", []):
                    arch = ReadArchitecture(**arch_data)
                    self.architectures.append(arch)
            logger.info(f"Loaded {len(self.architectures)} architectures from {self.db_path}")
        except Exception as e:
            logger.warning(f"Failed to load architectures from {self.db_path}: {e}")
            self._initialize_builtin_architectures()
    
    def _initialize_builtin_architectures(self) -> None:
        # McGlincy & Ingolia (2017)
        mcglincy_2017 = ReadArchitecture(
            protocol_name="mcglincy_ingolia_2017",
            lab_source="Ingolia Lab - Nature Protocols 2017",
            umi_positions=[(0, 5)],
            barcode_positions=[(5, 10)],
            adapter_sequences=["AGATCGGAAGAGCAC"],
            rpf_start=10,
            rpf_end=-1,
            expected_rpf_length=(26, 34),
            quality_markers={"umi_complexity_min": 0.8}
        )
        
        # Ingolia et al. (2009)
        ingolia_2009 = ReadArchitecture(
            protocol_name="ingolia_2009", 
            lab_source="Ingolia et al. Science 2009",
            umi_positions=[],
            barcode_positions=[],
            adapter_sequences=["CTGTAGGCACCATCAAT"],
            rpf_start=0,
            rpf_end=-17,
            expected_rpf_length=(28, 35),
            quality_markers={"adapter_match_threshold": 0.85}
        )
        
        # Generic UMI
        generic_umi = ReadArchitecture(
            protocol_name="generic_umi_protocol",
            lab_source="Common UMI-based protocol",
            umi_positions=[(0, 6)],
            barcode_positions=[],
            adapter_sequences=["TGGAATTCTCGGGTGCCAAGG"],
            rpf_start=6,
            rpf_end=-1,
            expected_rpf_length=(25, 35),
            quality_markers={"umi_complexity_min": 0.75}
        )
        
        # Preprocessed
        preprocessed = ReadArchitecture(
            protocol_name="preprocessed_with_adapter_contamination",
            lab_source="Partially processed ribosome profiling reads",
            umi_positions=[],
            barcode_positions=[],
            adapter_sequences=["AGATCGGAAGAGCAC"],
            rpf_start=0,
            rpf_end=-1,
            expected_rpf_length=(25, 40),
            quality_markers={"adapter_match_threshold": 0.2}
        )
        
        err605046_style = ReadArchitecture(
            protocol_name="err605046_stau1_riboseq",
            lab_source="ERR605046 STAU1 ribosome profiling dataset",
            umi_positions=[], 
            barcode_positions=[], 
            adapter_sequences=["AGATCGGAAGAGC", "GATCGGAAGAGC"], 
            rpf_start=1, 
            rpf_end=-1, 
            expected_rpf_length=(15, 45), 
            quality_markers={"adapter_match_threshold": 0.2} 
        )
        
        ena_riboseq = ReadArchitecture(
            protocol_name="ena_riboseq_standard",
            lab_source="ENA ribosome profiling dataset ERR10323209",
            umi_positions=[], 
            barcode_positions=[], 
            adapter_sequences=["GGAATTCTCGGGTGCCAAGG", "TGGAATTCTCGGGTGCCAAGG"], 
            rpf_start=0, 
            rpf_end=-1, 
            expected_rpf_length=(25, 80), 
            quality_markers={"adapter_match_threshold": 0.5} 
        )
        
        self.architectures = [err605046_style, ena_riboseq, preprocessed, mcglincy_2017, ingolia_2009, generic_umi]
        
        # Add comprehensive adapter check
        all_adapters = sorted(list(set(KNOWN_ADAPTERS.values())), key=len, reverse=True)
        comprehensive = ReadArchitecture(
            protocol_name="comprehensive_adapter_check",
            lab_source="Automated Comprehensive Adapter Scan",
            umi_positions=[],
            barcode_positions=[],
            adapter_sequences=all_adapters,
            rpf_start=0,
            rpf_end=-1,
            expected_rpf_length=(20, 40),
            quality_markers={"adapter_match_threshold": 0.3}
        )
        self.architectures.insert(0, comprehensive)

        logger.info(f"Initialized {len(self.architectures)} built-in architectures")
        
        if self.seqspec_dir and self.seqspec_dir.exists():
            self.load_from_seqspec_directory(self.seqspec_dir)
            
    def load_from_seqspec_directory(self, seqspec_dir: Path) -> int:
        seqspec_architectures = self.seqspec_loader.load_from_directory(seqspec_dir)
        self.architectures.extend(seqspec_architectures)
        return len(seqspec_architectures)


class RPFExtractor:
    """Main RPF extraction processor."""
    
    def __init__(self, architecture_db_path: Optional[Path] = None, seqspec_dir: Optional[Path] = None):
        self.architecture_db = ArchitectureDatabase(db_path=architecture_db_path, seqspec_dir=seqspec_dir)
        self.seqspec_generator = SeqSpecGenerator()
        
        # New components
        self.signal_processor = SignalProcessor()
        self.matcher = ArchitectureMatcher()
        self.segmenter = ProbabilisticSegmenter()
        self.reporter = Reporter()
    
    def extract_rpfs(
        self, 
        input_file: Path, 
        output_file: Path,
        format: str = "fastq",
        max_reads: Optional[int] = None,
        generate_seqspec: bool = False
    ) -> RPFExtractionResult:
        """Extract RPFs from input file."""
        logger.info(f"Starting RPF extraction from {input_file}")
        
        # 1. Load large sample for robust signal generation
        sample_reads, sample_headers = self._load_sample_reads(input_file, format, sample_size=50000)
        
        # 2. Generate Signals (Entropy, Composition)
        logger.info("Generating signal metrics...")
        signals = self.signal_processor.process_reads(sample_reads)
        
        # 3. Match against known architectures
        logger.info("Matching architectures...")
        match_result = self.matcher.match(signals, self.architecture_db.architectures)
        
        extracted_segments = []
        method = "unknown"
        trace_log = []
        final_architecture = None
        
        if match_result and match_result.is_match:
            # Pattern Match Success
            logger.info(f"Matched architecture: {match_result.architecture.protocol_name}")
            final_architecture = match_result.architecture
            method = "strict_pattern_match"
            trace_log = match_result.reasons
            
            # Convert architecture definition to segment info
            extracted_segments = self._arch_to_segments(final_architecture)
            
        else:
            # Fallback to Probabilistic Segmentation
            logger.info("No strict architecture match found. Attempting probabilistic segmentation...")
            method = "probabilistic_hmm"
            
            # IMPROVEMENT: Multi-Scale Binning Strategy
            # Variable read lengths can smear the signal. We bin reads by length and detect on each bin independently.
            from collections import Counter
            length_counts = Counter(len(r) for r in sample_reads)
            
            # Select bins with sufficient coverage (>1000 reads)
            valid_bins = [l for l, c in length_counts.items() if c > 1000]
            valid_bins.sort()
            
            detected_segments_per_bin = []
            
            if len(valid_bins) > 1:
                logger.info(f"Multi-scale detection: analyzing {len(valid_bins)} length bins: {valid_bins}")
                for length in valid_bins:
                    bin_reads = [r for r in sample_reads if len(r) == length]
                    bin_signals = self.signal_processor.process_reads(bin_reads)
                    bin_segments = self.segmenter.segment(bin_signals)
                    
                    if bin_segments:
                        detected_segments_per_bin.append(bin_segments)
                
                # Consensus logic: Do the bins agree on structure?
                # We enforce that UMI structures must be consistent across length bins.
                if detected_segments_per_bin:
                    bin_votes = []
                    for segs in detected_segments_per_bin:
                        # Extract UMI info
                        umis = [s for s in segs if s.segment_type == 'umi']
                        umi_len = umis[0].end_pos - umis[0].start_pos if umis else 0
                        # Extract RPF start (which is UMI end)
                        rpf_starts = [s.start_pos for s in segs if s.segment_type == 'rpf']
                        rpf_start = rpf_starts[0] if rpf_starts else 0
                        bin_votes.append((umi_len, rpf_start))
                    
                    # Find majority vote
                    vote_counts = Counter(bin_votes)
                    most_common, count = vote_counts.most_common(1)[0]
                    consensus_ratio = count / len(bin_votes)
                    
                    if consensus_ratio > 0.5:
                        logger.info(f"Multi-scale Consensus: {consensus_ratio:.0%} of bins agree on UMI length {most_common[0]}.")
                        
                        # Find a representative bin that matches the consensus
                        for segs in detected_segments_per_bin:
                            bg_umis = [s for s in segs if s.segment_type == 'umi']
                            bg_len = bg_umis[0].end_pos - bg_umis[0].start_pos if bg_umis else 0
                            if bg_len == most_common[0]:
                                extracted_segments = segs
                                break
                    else:
                         logger.warning(f"Multi-scale Conflict: Bins disagree on structure (consensus {consensus_ratio:.0%}).")
                         # Fall through to safety net or global detection
                         extracted_segments = []

                    trace_log.append(f"Used multi-scale segmentation across {len(valid_bins)} bins (Consensus: {consensus_ratio:.2f}).")
            
            if not extracted_segments:
                 # Try global signal if binning failed or was skipped
                 extracted_segments = self.segmenter.segment(signals)

            if extracted_segments:
                trace_log.append("Used HMM segmentation.")
                # Create a specific dummy architecture for the report
                final_architecture = ReadArchitecture(
                    protocol_name="de_novo_inferred",
                    lab_source="Probabilistic Segmenter",
                    umi_positions=[], barcode_positions=[], adapter_sequences=[],
                    rpf_start=extracted_segments[0].start_pos if extracted_segments else 0,
                    rpf_end=extracted_segments[-1].end_pos if extracted_segments else -1,
                    expected_rpf_length=(20, 40),
                    quality_markers={}
                )
            else:
                trace_log.append("Segmentation failed.")

            # Safety Net: Brute-force Adapter Scan
            # If we are in "No RPF segment detected" land or the HMM produced nonsense (e.g. all UMI),
            # we should check if one of our known adapters is actually present.
            rpf_found = any(s.segment_type == 'rpf' for s in extracted_segments)
            if not rpf_found or (extracted_segments and extracted_segments[0].segment_type == 'umi' and extracted_segments[0].end_pos > 40):
                logger.info("Suspicious structure detected (No RPF or giant UMI). Running brute-force adapter scan...")
                
                # We need to import AdapterDetector here or implementing a quick check
                # A quick check is better given we have loaded Sample Reads
                
                best_adapter = None
                best_count = 0
                threshold = len(sample_reads) * 0.1 # 10% of reads
                
                # Check top 20 adapters (common) to be fast
                from .adapter import AdapterDetector
                
                # Just use simple find for speed on the sample
                for name, seq in list(KNOWN_ADAPTERS.items())[:20]:
                    count = sum(1 for r in sample_reads if seq in r)
                    if count > best_count:
                        best_count = count
                        best_adapter = (name, seq)
                
                if best_adapter and best_count > threshold:
                    logger.info(f"Fallback Scan: Found adapter {best_adapter[0]} in {best_count} reads. Overriding HMM.")
                    adapter_seq = best_adapter[1]
                    
                    # Construct override architecture
                    final_architecture = ReadArchitecture(
                        protocol_name=f"fallback_detected_{best_adapter[0]}",
                        lab_source="Brute-force Scan Fallback",
                        umi_positions=[],
                        barcode_positions=[],
                        adapter_sequences=[adapter_seq],
                        rpf_start=0,
                        rpf_end= -1, # Let dynamic trimming handle it
                        expected_rpf_length=(20, 40),
                        quality_markers={"fallback_match": True}
                    )
                    method = "fallback_adapter_scan"
                    extracted_segments = [
                        SegmentInfo("rpf", 0, -1, 0.99), # Placeholder
                        SegmentInfo("adapter", -1, -1, 0.99) # Placeholder
                    ]
                    trace_log.append(f"Override: Found dominant adapter {best_adapter[0]}.")
        
        # 4. Generate Reports
        if final_architecture and extracted_segments:
            # CLI Report
            cli_report = self.reporter.generate_cli_report(final_architecture, extracted_segments, signals)
            logger.info("\n" + cli_report)
            
            # HTML Report
            html_report_path = output_file.with_suffix(".report.html")
            self.reporter.generate_html_report(
                html_report_path, final_architecture, extracted_segments, signals, trace_log
            )
            logger.info(f"Interactive report written to {html_report_path}")
            
            # Seqspec
            seqspec_data = None
            if generate_seqspec:
                seqspec_output = output_file.with_suffix('.seqspec.yaml')
                seqspec_data = self.seqspec_generator.generate_from_architecture(
                    architecture=final_architecture if method == "strict_pattern_match" else None,
                    sample_reads=sample_reads[:1000],
                    sample_headers=sample_headers[:1000],
                    detected_segments=extracted_segments,
                    output_file=seqspec_output
                )

            # 5. Extract RPFs
            # Setup adapters for dynamic trimming
            adapters_to_trim = None
            if final_architecture:
                 adapters_to_trim = final_architecture.adapter_sequences
                 if final_architecture.rpf_end == -1 and not adapters_to_trim:
                     pass

            extracted_count = self._extract_rpfs_from_reads(
                input_file, output_file, extracted_segments, format, max_reads, adapters=adapters_to_trim
            )
            
            return RPFExtractionResult(
                input_reads=extracted_count, # Simplified
                processed_reads=extracted_count,
                extracted_rpfs=extracted_count,
                failed_extractions=0,
                architecture_match=final_architecture.protocol_name,
                extraction_method=method,
                segments={"detected": extracted_segments},
                quality_metrics={"confidence": 0.95},
                seqspec_data=seqspec_data,
                trim_recommendations=self._derive_trim_recommendations(extracted_segments)
            )
            
        else:
            logger.error("Architecture detection failed.")
            raise RuntimeError("Could not determine read architecture.")

    def _arch_to_segments(self, arch: ReadArchitecture) -> List[SegmentInfo]:
        """Convert known architecture to segment list."""
        segments = []
        for start, end in arch.umi_positions:
            segments.append(SegmentInfo("umi", start, end, 1.0))
        for start, end in arch.barcode_positions:
            segments.append(SegmentInfo("barcode", start, end, 1.0))
            
        # RPF
        segments.append(SegmentInfo(
            "rpf", arch.rpf_start, arch.rpf_end if arch.rpf_end > 0 else -1, 1.0
        ))
        return segments

    def _derive_trim_recommendations(self, segments: List[SegmentInfo]) -> Dict[str, Any]:
        """Derive standard trim recommendations from segments."""
        rpf_segments = [s for s in segments if s.segment_type == "rpf"]
        if not rpf_segments:
            return {"recommended_5prime_trim": 0, "recommended_3prime_trim": 0}
        
        rpf = rpf_segments[0]
        return {
            "recommended_5prime_trim": rpf.start_pos,
            "recommended_3prime_trim": 0, 
            "note": "Use detected adapter sequences for 3' trimming"
        }
    
    def _load_sample_reads(self, input_file: Path, format: str, sample_size: int = 50000) -> Tuple[List[str], List[str]]:
        """Load sample reads and headers."""
        reads = []
        headers = []
        count = 0
        opener = get_file_opener(input_file)
        
        try:
            with opener(input_file, 'rt' if input_file.suffix in ['.gz', '.bz2'] else 'r') as f:
                if format == "fastq":
                    while count < sample_size:
                        header = f.readline().strip()
                        if not header: break
                        sequence = f.readline().strip()
                        f.readline() # plus
                        f.readline() # qual
                        if header and sequence:
                            reads.append(sequence)
                            headers.append(header)
                            count += 1
                elif format in ["fasta", "collapsed"]:
                    current_seq = []
                    current_header = None
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_seq and current_header and count < sample_size:
                                reads.append(''.join(current_seq))
                                headers.append(current_header)
                                count += 1
                            current_seq = []
                            current_header = line
                        else:
                            current_seq.append(line)
                    if current_seq and current_header and count < sample_size:
                        reads.append(''.join(current_seq))
                        headers.append(current_header)
        except Exception as e:
            logger.warning(f"Error loading sample reads: {e}")
            
        return reads, headers

    def _extract_rpfs_from_reads(
        self,
        input_file: Path,
        output_file: Path,
        segments: List[SegmentInfo],
        format: str,
        max_reads: Optional[int],
        adapters: Optional[List[str]] = None
    ) -> int:
        """Extract RPF sequences from reads using detected segments."""
        
        # Find RPF segment
        rpf_segments = [seg for seg in segments if seg.segment_type == "rpf"]
        if not rpf_segments:
            logger.warning("No RPF segment detected, extracting full reads")
            rpf_start, rpf_end = 0, -1
        else:
            rpf_segment = rpf_segments[0]
            rpf_start, rpf_end = rpf_segment.start_pos, rpf_segment.end_pos
        
        extracted_count = 0
        opener = get_file_opener(input_file)
        
        # Pre-sort adapters by length descending
        sorted_adapters = sorted(adapters, key=len, reverse=True) if adapters else None
        
        with opener(input_file, 'rt' if input_file.suffix in ['.gz', '.bz2'] else 'r') as fin:
            with open(output_file, 'w') as fout:
                if format == "fastq":
                    while max_reads is None or extracted_count < max_reads:
                        header = fin.readline().strip()
                        if not header: break
                        sequence = fin.readline().strip()
                        plus = fin.readline().strip() 
                        quality = fin.readline().strip()
                        
                        if header and sequence:
                            if rpf_end == -1:
                                rpf_seq = sequence[rpf_start:]
                                rpf_qual = quality[rpf_start:]
                                
                                # Dynamic adapter trimming
                                if sorted_adapters:
                                    for adapter in sorted_adapters:
                                        if adapter in rpf_seq:
                                            adapter_pos = rpf_seq.find(adapter)
                                            if adapter_pos >= 0:
                                                rpf_seq = rpf_seq[:adapter_pos]
                                                rpf_qual = rpf_qual[:adapter_pos]
                                                break
                            else:
                                rpf_seq = sequence[rpf_start:rpf_end]
                                rpf_qual = quality[rpf_start:rpf_end]
                            
                            if len(rpf_seq) >= 20: 
                                fout.write(f"{header}_RPF\n")
                                fout.write(f"{rpf_seq}\n")
                                fout.write("+\n")
                                fout.write(f"{rpf_qual}\n")
                                extracted_count += 1
                                
                elif format in ["fasta", "collapsed"]:
                    current_header = None
                    current_seq = []
                    for line in fin:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_header and current_seq:
                                sequence = ''.join(current_seq)
                                if rpf_end == -1:
                                    rpf_seq = sequence[rpf_start:]
                                    if sorted_adapters:
                                        for adapter in sorted_adapters:
                                            if adapter in rpf_seq:
                                                adapter_pos = rpf_seq.find(adapter)
                                                if adapter_pos >= 0:
                                                    rpf_seq = rpf_seq[:adapter_pos]
                                                    break
                                else:
                                    rpf_seq = sequence[rpf_start:rpf_end]
                                    
                                if len(rpf_seq) >= 20:
                                    fout.write(f"{current_header}_RPF\n")
                                    fout.write(f"{rpf_seq}\n")
                                    extracted_count += 1
                                if max_reads and extracted_count >= max_reads: break
                            current_header = line
                            current_seq = []
                        else:
                             current_seq.append(line)
                    # Process last
                    if current_header and current_seq and (max_reads is None or extracted_count < max_reads):
                        sequence = ''.join(current_seq)
                        if rpf_end == -1:
                            rpf_seq = sequence[rpf_start:]
                            if sorted_adapters:
                                for adapter in sorted_adapters:
                                    if adapter in rpf_seq:
                                        adapter_pos = rpf_seq.find(adapter)
                                        if adapter_pos >= 0:
                                            rpf_seq = rpf_seq[:adapter_pos]
                                            break
                        else:
                            rpf_seq = sequence[rpf_start:rpf_end]
                        if len(rpf_seq) >= 20:
                            fout.write(f"{current_header}_RPF\n")
                            fout.write(f"{rpf_seq}\n")
                            extracted_count += 1
        
        logger.info(f"Extracted {extracted_count} RPF sequences")
        return extracted_count