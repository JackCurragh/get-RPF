"""RPF extraction processor for automated ribosome protected fragment isolation.

This module implements the core RPF extraction system that combines:
1. Pattern matching against known read architectures  
2. De novo detection using multi-signal analysis
3. Segment classification and RPF isolation
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from collections import Counter
import re

from ...utils.file_utils import get_file_opener

logger = logging.getLogger(__name__)


@dataclass
class ReadArchitecture:
    """Represents a known read architecture from ribosome profiling protocols."""
    
    protocol_name: str                    # e.g., "ingolia_2009", "mcglincy_2014"
    lab_source: str                       # originating lab/publication
    umi_positions: List[Tuple[int, int]]  # list of (start, end) pairs for UMI regions
    barcode_positions: List[Tuple[int, int]]  # sample barcode positions
    adapter_sequences: List[str]          # known adapter sequences (5' and 3')
    rpf_start: int                        # start position of RPF
    rpf_end: int                          # end position (-1 for variable length)
    expected_rpf_length: Tuple[int, int]  # (min_length, max_length)
    quality_markers: Dict[str, Any]       # protocol-specific validation criteria
    confidence: float = 1.0               # confidence in this architecture


@dataclass
class SegmentInfo:
    """Information about a detected segment in reads."""
    
    segment_type: str        # "umi", "barcode", "adapter", "rpf", "unknown"
    start_pos: int          # start position in read
    end_pos: int            # end position in read
    confidence: float       # confidence score for this classification
    consensus: Optional[str] = None  # consensus sequence if applicable


@dataclass
class RPFExtractionResult:
    """Results from RPF extraction process."""
    
    input_reads: int
    processed_reads: int
    extracted_rpfs: int
    failed_extractions: int
    architecture_match: Optional[str]
    extraction_method: str  # "pattern_match" or "de_novo"
    segments: Dict[str, List[SegmentInfo]]
    quality_metrics: Dict[str, float]
    
    def write_report(self, output_path: Path, format: str = "json") -> None:
        """Write extraction results to file."""
        report_data = asdict(self)
        
        if format.lower() == "json":
            with open(output_path, 'w') as f:
                json.dump(report_data, f, indent=2)
        elif format.lower() == "csv":
            import csv
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["metric", "value"])
                for key, value in report_data.items():
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            writer.writerow([f"{key}_{subkey}", subvalue])
                    else:
                        writer.writerow([key, value])


class ArchitectureDatabase:
    """Database of known read architectures from ribosome profiling protocols."""
    
    def __init__(self, db_path: Optional[Path] = None):
        """Initialize architecture database.
        
        Args:
            db_path: Path to JSON database file, or None for built-in
        """
        self.architectures: List[ReadArchitecture] = []
        self.db_path = db_path
        self._load_architectures()
    
    def _load_architectures(self) -> None:
        """Load architectures from database file or initialize built-in ones."""
        if self.db_path and self.db_path.exists():
            self._load_from_file()
        else:
            self._initialize_builtin_architectures()
    
    def _load_from_file(self) -> None:
        """Load architectures from JSON file."""
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
        """Initialize built-in architectures from literature."""
        
        # McGlincy & Ingolia (2017) - Standard protocol
        mcglincy_2017 = ReadArchitecture(
            protocol_name="mcglincy_ingolia_2017",
            lab_source="Ingolia Lab - Nature Protocols 2017",
            umi_positions=[(0, 5)],  # 5nt UMI at start
            barcode_positions=[(5, 10)],  # 5nt sample barcode  
            adapter_sequences=["AGATCGGAAGAGCAC"],  # 3' adapter
            rpf_start=10,
            rpf_end=-1,  # Variable length
            expected_rpf_length=(26, 34),
            quality_markers={
                "umi_complexity_min": 0.8,
                "adapter_match_threshold": 0.9
            }
        )
        
        # Ingolia et al. (2009) - Original protocol
        ingolia_2009 = ReadArchitecture(
            protocol_name="ingolia_2009", 
            lab_source="Ingolia et al. Science 2009",
            umi_positions=[],  # No UMI in original
            barcode_positions=[],
            adapter_sequences=["CTGTAGGCACCATCAAT"],  # 3' linker
            rpf_start=0,
            rpf_end=-17,  # Remove last 17nt (adapter)
            expected_rpf_length=(28, 35),
            quality_markers={
                "adapter_match_threshold": 0.85,
                "rpf_gc_content_range": (0.3, 0.7)
            }
        )
        
        # Generic UMI protocol (common variant)
        generic_umi = ReadArchitecture(
            protocol_name="generic_umi_protocol",
            lab_source="Common UMI-based protocol",
            umi_positions=[(0, 6)],  # 6nt UMI
            barcode_positions=[],
            adapter_sequences=["TGGAATTCTCGGGTGCCAAGG"],  # Common TruSeq adapter
            rpf_start=6,
            rpf_end=-1,
            expected_rpf_length=(25, 35),
            quality_markers={
                "umi_complexity_min": 0.75
            }
        )
        
        self.architectures = [mcglincy_2017, ingolia_2009, generic_umi]
        logger.info(f"Initialized {len(self.architectures)} built-in architectures")
    
    def match_architecture(self, reads: List[str]) -> Tuple[Optional[ReadArchitecture], float]:
        """Find best matching architecture for sample reads.
        
        Args:
            reads: Sample of reads to analyze (first 1000-10000)
            
        Returns:
            Tuple of (best_architecture, confidence_score)
        """
        best_arch = None
        best_score = 0.0
        
        for architecture in self.architectures:
            score = self._score_architecture(reads, architecture)
            logger.debug(f"Architecture {architecture.protocol_name} scored {score:.3f}")
            
            if score > best_score:
                best_score = score
                best_arch = architecture
        
        # Require minimum confidence for positive match
        confidence_threshold = 0.7
        if best_score < confidence_threshold:
            return None, best_score
        
        logger.info(f"Best architecture match: {best_arch.protocol_name} (score: {best_score:.3f})")
        return best_arch, best_score
    
    def _score_architecture(self, reads: List[str], architecture: ReadArchitecture) -> float:
        """Score how well reads match a given architecture.
        
        Args:
            reads: Sample reads
            architecture: Architecture to test
            
        Returns:
            Confidence score 0.0-1.0
        """
        scores = []
        
        # Score adapter presence
        adapter_score = self._score_adapter_presence(reads, architecture)
        scores.append(("adapter", adapter_score, 3.0))  # Weight 3
        
        # Score length distribution
        length_score = self._score_length_distribution(reads, architecture) 
        scores.append(("length", length_score, 2.0))  # Weight 2
        
        # Score UMI complexity if UMIs expected
        if architecture.umi_positions:
            umi_score = self._score_umi_complexity(reads, architecture)
            scores.append(("umi", umi_score, 2.0))  # Weight 2
        
        # Calculate weighted average
        total_score = sum(score * weight for _, score, weight in scores)
        total_weight = sum(weight for _, _, weight in scores)
        
        return total_score / total_weight if total_weight > 0 else 0.0
    
    def _score_adapter_presence(self, reads: List[str], architecture: ReadArchitecture) -> float:
        """Score presence of expected adapter sequences."""
        if not architecture.adapter_sequences:
            return 1.0  # No adapters expected
        
        matches = 0
        for read in reads[:1000]:  # Sample first 1000 reads
            for adapter in architecture.adapter_sequences:
                if self._fuzzy_match(read, adapter, max_mismatches=2):
                    matches += 1
                    break
        
        return matches / min(len(reads), 1000)
    
    def _score_length_distribution(self, reads: List[str], architecture: ReadArchitecture) -> float:
        """Score how well read lengths match expected RPF lengths."""
        min_len, max_len = architecture.expected_rpf_length
        
        # Estimate RPF lengths after removing known elements
        rpf_lengths = []
        for read in reads[:1000]:
            estimated_rpf_len = len(read)
            
            # Subtract known element lengths
            for start, end in architecture.umi_positions:
                estimated_rpf_len -= (end - start)
            for start, end in architecture.barcode_positions:
                estimated_rpf_len -= (end - start)
            for adapter in architecture.adapter_sequences:
                estimated_rpf_len -= len(adapter)
                
            if estimated_rpf_len > 0:
                rpf_lengths.append(estimated_rpf_len)
        
        if not rpf_lengths:
            return 0.0
        
        # Score based on how many fall in expected range
        in_range = sum(1 for length in rpf_lengths if min_len <= length <= max_len)
        return in_range / len(rpf_lengths)
    
    def _score_umi_complexity(self, reads: List[str], architecture: ReadArchitecture) -> float:
        """Score UMI region complexity (should be random)."""
        if not architecture.umi_positions:
            return 1.0
        
        umi_sequences = []
        for read in reads[:1000]:
            for start, end in architecture.umi_positions:
                if end <= len(read):
                    umi = read[start:end]
                    umi_sequences.append(umi)
        
        if not umi_sequences:
            return 0.0
        
        # Calculate complexity as unique fraction
        unique_umis = len(set(umi_sequences))
        total_umis = len(umi_sequences)
        complexity = unique_umis / total_umis
        
        # Compare to expected minimum
        min_complexity = architecture.quality_markers.get("umi_complexity_min", 0.7)
        return min(complexity / min_complexity, 1.0)
    
    def _fuzzy_match(self, text: str, pattern: str, max_mismatches: int = 2) -> bool:
        """Check if pattern matches text with up to max_mismatches differences."""
        # Try exact pattern match first
        if pattern in text:
            return True
        
        # Try partial pattern matches (at least 8 characters)
        min_match_len = min(8, len(pattern))
        for start in range(len(pattern) - min_match_len + 1):
            partial_pattern = pattern[start:start + min_match_len]
            if partial_pattern in text:
                return True
        
        # Try fuzzy matching with mismatches
        for i in range(len(text) - len(pattern) + 1):
            substring = text[i:i+len(pattern)]
            mismatches = sum(c1 != c2 for c1, c2 in zip(substring, pattern))
            if mismatches <= max_mismatches:
                return True
        return False


class DeNovoDetector:
    """De novo architecture detection using multi-signal analysis."""
    
    def __init__(self):
        """Initialize de novo detector."""
        pass
    
    def detect_architecture(self, reads: List[str]) -> Tuple[List[SegmentInfo], float]:
        """Detect read architecture using multi-signal analysis.
        
        Args:
            reads: Sample reads for analysis
            
        Returns:
            Tuple of (detected_segments, confidence)
        """
        # Combine all detection signals
        signals = self._analyze_all_signals(reads)
        
        # Find segment boundaries
        boundaries = self._find_segment_boundaries(signals, reads)
        
        # Classify segments
        segments = self._classify_segments(boundaries, reads)
        
        # Calculate overall confidence
        confidence = self._calculate_confidence(segments, signals)
        
        return segments, confidence
    
    def _analyze_all_signals(self, reads: List[str]) -> Dict[str, Any]:
        """Analyze all detection signals."""
        return {
            "kmer_complexity": self._analyze_kmer_complexity(reads),
            "composition": self._analyze_composition(reads), 
            "quality": self._analyze_quality_patterns(reads),
            "motifs": self._analyze_motif_boundaries(reads)
        }
    
    def _analyze_kmer_complexity(self, reads: List[str]) -> List[float]:
        """Analyze k-mer complexity across read positions."""
        if not reads:
            return []
        
        max_len = max(len(read) for read in reads[:100])
        complexity_scores = []
        
        for pos in range(max_len):
            kmers = []
            for read in reads[:100]:  # Sample reads
                if pos + 3 <= len(read):
                    kmer = read[pos:pos+3]
                    kmers.append(kmer)
            
            if kmers:
                unique_kmers = len(set(kmers))
                total_kmers = len(kmers)
                complexity = unique_kmers / total_kmers
            else:
                complexity = 0.0
            
            complexity_scores.append(complexity)
        
        return complexity_scores
    
    def _analyze_composition(self, reads: List[str]) -> List[Dict[str, float]]:
        """Analyze nucleotide composition across positions."""
        if not reads:
            return []
        
        max_len = max(len(read) for read in reads[:100])
        composition_scores = []
        
        for pos in range(max_len):
            nucleotides = []
            for read in reads[:100]:
                if pos < len(read):
                    nucleotides.append(read[pos])
            
            if nucleotides:
                counter = Counter(nucleotides)
                total = len(nucleotides)
                gc_content = (counter.get('G', 0) + counter.get('C', 0)) / total
                entropy = self._calculate_entropy(counter.values(), total)
            else:
                gc_content = 0.5
                entropy = 0.0
            
            composition_scores.append({
                "gc_content": gc_content,
                "entropy": entropy
            })
        
        return composition_scores
    
    def _analyze_quality_patterns(self, reads: List[str]) -> List[float]:
        """Analyze quality score patterns (placeholder - would use real FASTQ quality)."""
        # For now, return uniform scores since we're working with sequences only
        max_len = max(len(read) for read in reads[:10]) if reads else 0
        return [0.5] * max_len  # Neutral quality signal
    
    def _analyze_motif_boundaries(self, reads: List[str]) -> Dict[str, Any]:
        """Analyze sequence motifs and their boundaries."""
        # Analyze 5' and 3' ends for consensus sequences
        five_prime_motifs = self._find_consensus_motifs([read[:15] for read in reads[:100]])
        three_prime_motifs = self._find_consensus_motifs([read[-15:] for read in reads[:100]])
        
        return {
            "five_prime": five_prime_motifs,
            "three_prime": three_prime_motifs
        }
    
    def _find_consensus_motifs(self, sequences: List[str]) -> Dict[str, float]:
        """Find consensus motifs in sequence list."""
        if not sequences:
            return {}
        
        motifs = {}
        for i in range(min(len(seq) for seq in sequences if seq)):
            nucleotides = [seq[i] for seq in sequences if i < len(seq)]
            if nucleotides:
                counter = Counter(nucleotides)
                most_common = counter.most_common(1)[0]
                consensus_fraction = most_common[1] / len(nucleotides)
                motifs[f"pos_{i}"] = consensus_fraction
        
        return motifs
    
    def _calculate_entropy(self, counts, total: int) -> float:
        """Calculate Shannon entropy for nucleotide distribution."""
        import math
        entropy = 0.0
        for count in counts:
            if count > 0:
                p = count / total
                entropy -= p * math.log2(p)
        return entropy
    
    def _find_segment_boundaries(self, signals: Dict[str, Any], reads: List[str]) -> List[int]:
        """Find segment boundaries using combined signals."""
        if not reads:
            return []
        
        kmer_complexity = signals["kmer_complexity"]
        composition = signals["composition"]
        
        boundaries = []
        
        # Look for sharp transitions in k-mer complexity (use moving window)
        window_size = 3
        for i in range(window_size, len(kmer_complexity) - window_size):
            # Calculate average complexity before and after this position
            before_avg = sum(kmer_complexity[i-window_size:i]) / window_size
            after_avg = sum(kmer_complexity[i:i+window_size]) / window_size
            complexity_change = abs(after_avg - before_avg)
            
            if complexity_change > 0.4:  # Strong transition
                boundaries.append(i)
        
        # Look for composition changes (also use moving window)
        for i in range(window_size, len(composition) - window_size):
            before_gc = sum(comp["gc_content"] for comp in composition[i-window_size:i]) / window_size
            after_gc = sum(comp["gc_content"] for comp in composition[i:i+window_size]) / window_size
            gc_change = abs(after_gc - before_gc)
            
            if gc_change > 0.3:  # Significant GC change
                boundaries.append(i)
        
        # Remove boundaries that are too close together (minimum 5nt segments)
        filtered_boundaries = []
        for boundary in sorted(set(boundaries)):
            if not filtered_boundaries or boundary - filtered_boundaries[-1] >= 5:
                filtered_boundaries.append(boundary)
        
        return filtered_boundaries
    
    def _classify_segments(self, boundaries: List[int], reads: List[str]) -> List[SegmentInfo]:
        """Classify segments between boundaries."""
        if not boundaries or not reads:
            return []
        
        # Add start and end boundaries
        all_boundaries = [0] + boundaries + [len(reads[0]) if reads else 0]
        segments = []
        
        for i in range(len(all_boundaries) - 1):
            start = all_boundaries[i]
            end = all_boundaries[i + 1]
            
            # Analyze segment to classify type
            segment_seqs = [read[start:end] for read in reads[:100] if end <= len(read)]
            if segment_seqs:
                segment_type, confidence = self._classify_segment_type(segment_seqs, start, end)
                segments.append(SegmentInfo(
                    segment_type=segment_type,
                    start_pos=start,
                    end_pos=end,
                    confidence=confidence
                ))
        
        return segments
    
    def _classify_segment_type(self, sequences: List[str], start_pos: int, end_pos: int) -> Tuple[str, float]:
        """Classify what type of segment this is."""
        if not sequences:
            return "unknown", 0.0
        
        avg_length = sum(len(seq) for seq in sequences) / len(sequences)
        unique_fraction = len(set(sequences)) / len(sequences)
        segment_length = end_pos - start_pos
        
        # UMI: short (4-12nt), high complexity, at start
        if segment_length <= 12 and unique_fraction > 0.7 and start_pos < 20:
            return "umi", 0.9
        
        # Barcode: short-medium (4-12nt), medium complexity, early position
        elif segment_length <= 12 and 0.3 < unique_fraction < 0.8 and start_pos < 25:
            return "barcode", 0.7
        
        # Adapter: low complexity (consensus), any position
        elif unique_fraction < 0.4:
            return "adapter", 0.8
        
        # RPF: long segment (20+ nt), biological complexity, not at very start
        elif segment_length >= 20 and 0.4 < unique_fraction < 0.95:
            return "rpf", 0.95
        
        # Large unknown segment likely RPF if reasonable complexity
        elif segment_length >= 25 and unique_fraction > 0.3:
            return "rpf", 0.7
        
        else:
            return "unknown", 0.3
    
    def _calculate_confidence(self, segments: List[SegmentInfo], signals: Dict[str, Any]) -> float:
        """Calculate overall confidence in the detected architecture."""
        if not segments:
            return 0.0
        
        # Base confidence on segment classification confidence
        segment_confidences = [seg.confidence for seg in segments]
        avg_confidence = sum(segment_confidences) / len(segment_confidences)
        
        # Boost if we found expected segments
        has_rpf = any(seg.segment_type == "rpf" for seg in segments)
        has_structure = len(segments) > 1
        
        confidence = avg_confidence
        if has_rpf:
            confidence += 0.1
        if has_structure:
            confidence += 0.1
        
        return min(confidence, 1.0)


class RPFExtractor:
    """Main RPF extraction processor combining pattern matching and de novo detection."""
    
    def __init__(self, architecture_db_path: Optional[Path] = None):
        """Initialize RPF extractor.
        
        Args:
            architecture_db_path: Path to architecture database
        """
        self.architecture_db = ArchitectureDatabase(architecture_db_path)
        self.de_novo_detector = DeNovoDetector()
    
    def extract_rpfs(
        self, 
        input_file: Path, 
        output_file: Path,
        format: str = "fastq",
        max_reads: Optional[int] = None
    ) -> RPFExtractionResult:
        """Extract RPFs from input file.
        
        Args:
            input_file: Input sequence file
            output_file: Output file for extracted RPFs  
            format: Input file format
            max_reads: Maximum reads to process
            
        Returns:
            RPFExtractionResult with extraction statistics
        """
        logger.info(f"Starting RPF extraction from {input_file}")
        
        # Load sample reads for architecture detection
        sample_reads = self._load_sample_reads(input_file, format, sample_size=10000)
        
        # Try pattern matching first
        matched_arch, match_confidence = self.architecture_db.match_architecture(sample_reads)
        
        if matched_arch and match_confidence > 0.7:
            logger.info(f"Using pattern matching with {matched_arch.protocol_name}")
            return self._extract_with_pattern_matching(
                input_file, output_file, matched_arch, format, max_reads
            )
        else:
            logger.info("Pattern matching failed, using de novo detection")
            return self._extract_with_de_novo(
                input_file, output_file, format, max_reads
            )
    
    def _load_sample_reads(self, input_file: Path, format: str, sample_size: int = 10000) -> List[str]:
        """Load sample reads for analysis."""
        reads = []
        count = 0
        
        opener = get_file_opener(input_file)
        with opener(input_file, 'rt' if input_file.suffix in ['.gz', '.bz2'] else 'r') as f:
            if format == "fastq":
                # Read FASTQ format
                while count < sample_size:
                    header = f.readline().strip()
                    if not header:
                        break
                    sequence = f.readline().strip()
                    plus = f.readline().strip()
                    quality = f.readline().strip()
                    
                    if header and sequence:
                        reads.append(sequence)
                        count += 1
                        
            elif format in ["fasta", "collapsed"]:
                # Read FASTA format
                current_seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq and count < sample_size:
                            reads.append(''.join(current_seq))
                            count += 1
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Add last sequence
                if current_seq and count < sample_size:
                    reads.append(''.join(current_seq))
        
        logger.info(f"Loaded {len(reads)} sample reads for analysis")
        return reads
    
    def _extract_with_pattern_matching(
        self,
        input_file: Path,
        output_file: Path, 
        architecture: ReadArchitecture,
        format: str,
        max_reads: Optional[int]
    ) -> RPFExtractionResult:
        """Extract RPFs using known architecture pattern."""
        
        # Create segments from architecture
        segments = []
        for start, end in architecture.umi_positions:
            segments.append(SegmentInfo("umi", start, end, 0.95))
        for start, end in architecture.barcode_positions:
            segments.append(SegmentInfo("barcode", start, end, 0.95))
        
        # Add RPF segment
        rpf_segment = SegmentInfo(
            "rpf", 
            architecture.rpf_start,
            architecture.rpf_end if architecture.rpf_end > 0 else -1,
            0.95
        )
        segments.append(rpf_segment)
        
        # Process reads and extract RPFs
        extracted_count = self._extract_rpfs_from_reads(
            input_file, output_file, segments, format, max_reads
        )
        
        return RPFExtractionResult(
            input_reads=extracted_count,  # Would count in real implementation
            processed_reads=extracted_count,
            extracted_rpfs=extracted_count,
            failed_extractions=0,
            architecture_match=architecture.protocol_name,
            extraction_method="pattern_match",
            segments={"detected": segments},
            quality_metrics={"extraction_rate": 1.0}
        )
    
    def _extract_with_de_novo(
        self,
        input_file: Path,
        output_file: Path,
        format: str, 
        max_reads: Optional[int]
    ) -> RPFExtractionResult:
        """Extract RPFs using de novo detection."""
        
        # Load sample for de novo analysis
        sample_reads = self._load_sample_reads(input_file, format, 1000)
        
        # Detect architecture
        segments, confidence = self.de_novo_detector.detect_architecture(sample_reads)
        
        # Extract RPFs using detected segments
        extracted_count = self._extract_rpfs_from_reads(
            input_file, output_file, segments, format, max_reads
        )
        
        return RPFExtractionResult(
            input_reads=extracted_count,
            processed_reads=extracted_count, 
            extracted_rpfs=extracted_count,
            failed_extractions=0,
            architecture_match=None,
            extraction_method="de_novo",
            segments={"detected": segments},
            quality_metrics={"detection_confidence": confidence}
        )
    
    def _extract_rpfs_from_reads(
        self,
        input_file: Path,
        output_file: Path,
        segments: List[SegmentInfo],
        format: str,
        max_reads: Optional[int]
    ) -> int:
        """Extract RPF sequences from reads using detected segments."""
        
        # Find RPF segment
        rpf_segments = [seg for seg in segments if seg.segment_type == "rpf"]
        if not rpf_segments:
            logger.warning("No RPF segment detected, extracting full reads")
            # Fallback: use full reads as RPFs
            rpf_start, rpf_end = 0, -1
        else:
            rpf_segment = rpf_segments[0]  # Use first RPF segment
            rpf_start, rpf_end = rpf_segment.start_pos, rpf_segment.end_pos
        
        extracted_count = 0
        opener = get_file_opener(input_file)
        
        with opener(input_file, 'rt' if input_file.suffix in ['.gz', '.bz2'] else 'r') as fin:
            with open(output_file, 'w') as fout:
                if format == "fastq":
                    # Process FASTQ
                    while max_reads is None or extracted_count < max_reads:
                        header = fin.readline().strip()
                        if not header:
                            break
                        sequence = fin.readline().strip()
                        plus = fin.readline().strip() 
                        quality = fin.readline().strip()
                        
                        if header and sequence:
                            # Extract RPF portion
                            if rpf_end == -1:
                                rpf_seq = sequence[rpf_start:]
                                rpf_qual = quality[rpf_start:]
                            else:
                                rpf_seq = sequence[rpf_start:rpf_end]
                                rpf_qual = quality[rpf_start:rpf_end]
                            
                            # Write extracted RPF
                            if len(rpf_seq) >= 20:  # Minimum RPF length filter
                                fout.write(f"{header}_RPF\n")
                                fout.write(f"{rpf_seq}\n")
                                fout.write("+\n")
                                fout.write(f"{rpf_qual}\n")
                                extracted_count += 1
                
                elif format in ["fasta", "collapsed"]:
                    # Process FASTA
                    current_header = None
                    current_seq = []
                    
                    for line in fin:
                        line = line.strip()
                        if line.startswith('>'):
                            # Process previous sequence
                            if current_header and current_seq:
                                sequence = ''.join(current_seq)
                                
                                # Extract RPF portion
                                if rpf_end == -1:
                                    rpf_seq = sequence[rpf_start:]
                                else:
                                    rpf_seq = sequence[rpf_start:rpf_end]
                                
                                # Write extracted RPF
                                if len(rpf_seq) >= 20:
                                    fout.write(f"{current_header}_RPF\n")
                                    fout.write(f"{rpf_seq}\n")
                                    extracted_count += 1
                                
                                if max_reads and extracted_count >= max_reads:
                                    break
                            
                            current_header = line
                            current_seq = []
                        else:
                            current_seq.append(line)
                    
                    # Process last sequence
                    if current_header and current_seq and (max_reads is None or extracted_count < max_reads):
                        sequence = ''.join(current_seq)
                        
                        if rpf_end == -1:
                            rpf_seq = sequence[rpf_start:]
                        else:
                            rpf_seq = sequence[rpf_start:rpf_end]
                        
                        if len(rpf_seq) >= 20:
                            fout.write(f"{current_header}_RPF\n")
                            fout.write(f"{rpf_seq}\n")
                            extracted_count += 1
        
        logger.info(f"Extracted {extracted_count} RPF sequences")
        return extracted_count