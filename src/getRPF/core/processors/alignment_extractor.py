"""Alignment-based RPF extraction using STAR as ground truth.

This is the primary extraction method for getRPF. It uses reference genome
alignment to determine biological sequence boundaries with reliable structure
learning based on per-position entropy profiles.

Algorithm:
1. Sample subset of reads (default 10k)
2. Align subset with STAR
3. Learn read structure from soft-clipping patterns (per-position entropy)
4. Validate structure against known adapter database
5. Extract RPF from all reads using learned structure
6. Export learned structure as seqspec YAML

Author: getRPF team
"""

import logging
import math
import yaml
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

from Bio import SeqIO

from .alignment import STARAligner
from .types import (
    PositionProfile,
    RegionClassification,
    LearnedStructure,
    TrimmerConfig,
    StructureLearningError,
)
from ...utils.file_utils import create_temp_file, get_file_opener

logger = logging.getLogger(__name__)


# =============================================================================
# Data Structures (Local to this module)
# =============================================================================

@dataclass
class AdapterInfo:
    """Adapter detection results for reporting."""
    detected_adapters: List[Dict]
    total_reads_scanned: int
    no_adapter_fraction: float

    @property
    def top_adapter(self) -> Optional[str]:
        if self.detected_adapters:
            return self.detected_adapters[0]['name']
        return None

    def to_dict(self) -> Dict:
        return {
            'detected_adapters': self.detected_adapters,
            'total_reads_scanned': self.total_reads_scanned,
            'no_adapter_fraction': self.no_adapter_fraction
        }


@dataclass
class TrimBoundaries:
    """Consensus trim boundaries from alignment."""
    consensus_5p: int
    consensus_3p: int
    per_length: Dict[int, Dict[str, int]]
    confidence: str

    def to_dict(self) -> Dict:
        return {
            'consensus_5p': self.consensus_5p,
            'consensus_3p': self.consensus_3p,
            'per_length': self.per_length,
            'confidence': self.confidence
        }


@dataclass
class UMIInfo:
    """UMI detection results."""
    detected: bool
    positions: Optional[Tuple[int, int]]
    entropy: float
    confidence: str
    sample_sequences: List[str]

    def to_dict(self) -> Dict:
        return {
            'umi_detected': self.detected,
            'positions': list(self.positions) if self.positions else None,
            'entropy': self.entropy,
            'confidence': self.confidence,
            'sample_sequences': self.sample_sequences[:10]
        }


@dataclass
class ExtractionResult:
    """Complete extraction results."""
    input_reads: int
    extracted_rpfs: int
    extraction_rate: float
    sample_size: int
    trim_boundaries: TrimBoundaries
    learned_structure: LearnedStructure
    adapter_info: Optional[AdapterInfo]
    umi_info: Optional[UMIInfo]
    alignment_stats: Dict
    method: str = "structure_learning"

    def to_dict(self) -> Dict:
        return {
            'extraction_summary': {
                'input_reads': self.input_reads,
                'extracted_rpfs': self.extracted_rpfs,
                'extraction_rate': self.extraction_rate,
                'sample_size': self.sample_size,
                'method': self.method
            },
            'trim_boundaries': self.trim_boundaries.to_dict(),
            'learned_structure': self.learned_structure.to_dict(),
            'adapter_scan': self.adapter_info.to_dict() if self.adapter_info else None,
            'umi_detection': self.umi_info.to_dict() if self.umi_info else None,
            'alignment_statistics': self.alignment_stats
        }


# =============================================================================
# Adapter Database
# =============================================================================

KNOWN_ADAPTERS = [
    ("Illumina TruSeq Universal", "AGATCGGAAGAGCACACGTCT"),
    ("Illumina Small RNA 3'", "TGGAATTCTCGGGTGCCAAGG"),
    ("Illumina TruSeq (partial)", "GATCGGAAGAGCACACGT"),
    ("Illumina Universal (short)", "AGATCGGAAGAGC"),
    ("Illumina Multiplexing", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
    ("Ingolia 2009", "CTGTAGGCACCATCAAT"),
    ("McGlincy 2017 (with barcode)", "AGATCGGAAGAGCACACGTCTGAA"),
    ("Nextera Transposase", "CTGTCTCTTATACACATCT"),
    ("PolyA tail", "AAAAAAAAAA"),
    ("NEBNext Small RNA", "AGATCGGAAGAGCACACGTCT"),
    ("QIAseq miRNA", "AACTGTAGGCACCATCAAT"),
    ("Generic Illumina 1", "GATCGGAAGAGCGTCGT"),
    ("Generic Illumina 2", "GATCGGAAGAGCTCGTA"),
    ("Generic Illumina 3", "TGATCGGAAGAGCACAC"),
]


# =============================================================================
# Main Extractor Class
# =============================================================================

class AlignmentBasedExtractor:
    """
    Primary RPF extraction using alignment as ground truth.

    This extractor uses reference genome alignment to determine where
    biological sequence (RPF) ends and synthetic sequences (adapters,
    UMIs) begin. Structure is learned using per-position entropy profiles
    for reliable UMI vs adapter classification.
    """

    # Thresholds for structure classification
    ENTROPY_RANDOM_THRESHOLD = 1.7    # Above this = random (UMI-like)
    ENTROPY_CONSERVED_THRESHOLD = 0.5  # Below this = conserved (adapter-like)
    DOMINANT_FREQ_THRESHOLD = 0.8      # Dominant base frequency for conserved
    MIN_COVERAGE_PER_POSITION = 500    # Minimum reads to trust per-position stats
    ADAPTER_MATCH_THRESHOLD = 0.8      # Minimum identity for adapter database match

    def __init__(self, adapters: Optional[List[Tuple[str, str]]] = None):
        """Initialize extractor."""
        self.adapters = adapters or KNOWN_ADAPTERS
        logger.info(f"Initialized with {len(self.adapters)} known adapters")

    def extract(
        self,
        input_file: Path,
        output_file: Path,
        star_index: Path,
        preserve_umi: bool = False,
        sample_size: int = 10000,
        report_adapters: bool = True,
        star_threads: int = 1,
    ) -> ExtractionResult:
        """
        Extract RPF sequences using alignment-based structure learning.

        Args:
            input_file: Path to input FASTQ/FASTA file
            output_file: Path for output RPF sequences
            star_index: Path to STAR genome index
            preserve_umi: Whether to detect and preserve UMI sequences
            sample_size: Number of reads to analyze for boundaries
            report_adapters: Whether to scan for known adapters
            star_threads: Number of threads for STAR alignment

        Returns:
            ExtractionResult with statistics and metadata

        Raises:
            StructureLearningError: If structure cannot be reliably determined
        """
        logger.info("=" * 60)
        logger.info("Starting Alignment-Based RPF Extraction")
        logger.info("=" * 60)

        # Phase 1: Sample reads
        logger.info(f"Phase 1/4: Sampling {sample_size} reads...")
        subset_reads = self._sample_reads(input_file, sample_size)
        logger.info(f"  Sampled {len(subset_reads)} reads")

        # Phase 2: Characterization
        logger.info("Phase 2/4: Characterizing sample...")

        # 2a. Adapter scan (informational only)
        adapter_info = None
        if report_adapters:
            logger.info("  Scanning for known adapters...")
            adapter_info = self._scan_adapters(subset_reads)
            if adapter_info.top_adapter:
                top = adapter_info.detected_adapters[0]
                logger.info(f"    Detected: {top['name']} ({top['frequency']:.1%})")
            else:
                logger.info("    No known adapters detected")

        # 2b. STAR alignment
        logger.info("  Aligning subset with STAR...")
        subset_bam_path = output_file.with_suffix('.subset.bam')
        alignment_result = self._align_subset(
            subset_reads, star_index, star_threads, subset_bam_path
        )
        logger.info(f"    Aligned: {alignment_result['aligned_reads']}/{alignment_result['total_reads']} "
                   f"({alignment_result['alignment_rate']:.1%})")

        # 2c. Analyze soft-clipping for boundaries
        logger.info("  Analyzing soft-clipping patterns...")
        boundaries = self._analyze_soft_clipping(alignment_result['bam_file'])
        logger.info(f"    Consensus: 5'={boundaries.consensus_5p}nt, 3'={boundaries.consensus_3p}nt")

        # Phase 3: Structure Learning (THE KEY IMPROVEMENT)
        logger.info("Phase 3/4: Learning read structure (per-position entropy)...")
        learned_structure = self._learn_read_structure_reliable(
            alignment_result['bam_file']
        )

        # Log structure details
        self._log_learned_structure(learned_structure)

        # Validate structure - fail fast if low confidence
        if learned_structure.overall_confidence == 'low':
            error_msg = (
                f"Structure learning failed with low confidence. "
                f"Warnings: {learned_structure.validation_warnings}. "
                f"Please provide a custom seqspec file or manually inspect the data."
            )
            logger.error(error_msg)
            raise StructureLearningError(error_msg)

        # Phase 4: Extraction
        logger.info("Phase 4/4: Extracting RPFs using learned structure...")
        config = learned_structure.to_trimmer_config()
        extraction_stats = self._extract_with_config(
            input_file, output_file, config, preserve_umi
        )
        logger.info(f"  Extracted {extraction_stats['extracted_rpfs']} RPFs from "
                   f"{extraction_stats['total_reads']} reads ({extraction_stats['extraction_rate']:.1%})")

        # Export seqspec YAML
        seqspec_path = output_file.with_suffix('.seqspec.yaml')
        self._export_seqspec(learned_structure, seqspec_path)
        logger.info(f"  Exported learned structure to {seqspec_path}")

        # UMI detection (for reporting)
        umi_info = None
        if preserve_umi and learned_structure.five_prime:
            if learned_structure.five_prime.region_type == 'umi':
                umi_info = UMIInfo(
                    detected=True,
                    positions=(0, learned_structure.five_prime.length),
                    entropy=learned_structure.five_prime.evidence.get('avg_entropy', 0.0),
                    confidence=learned_structure.overall_confidence,
                    sample_sequences=learned_structure.five_prime.evidence.get('sample_sequences', [])
                )

        logger.info("=" * 60)
        logger.info("Extraction Complete!")
        logger.info("=" * 60)

        return ExtractionResult(
            input_reads=extraction_stats['total_reads'],
            extracted_rpfs=extraction_stats['extracted_rpfs'],
            extraction_rate=extraction_stats['extraction_rate'],
            sample_size=len(subset_reads),
            trim_boundaries=boundaries,
            learned_structure=learned_structure,
            adapter_info=adapter_info,
            umi_info=umi_info,
            alignment_stats=alignment_result,
            method="structure_learning"
        )

    # =========================================================================
    # Phase 1: Sampling
    # =========================================================================

    def _sample_reads(self, input_file: Path, sample_size: int) -> List[Tuple[str, str]]:
        """Sample reads from input file."""
        reads = []
        format_type = self._detect_format(input_file)

        opener = get_file_opener(input_file)
        with opener(input_file, 'rt') as f:
            for i, record in enumerate(SeqIO.parse(f, format_type)):
                if i >= sample_size:
                    break
                reads.append((record.id, str(record.seq)))

        return reads

    def _detect_format(self, input_file: Path) -> str:
        """Detect file format from extension."""
        name = input_file.name.lower()
        if name.endswith('.gz'):
            name = name[:-3]
        if name.endswith(('.fq', '.fastq')):
            return 'fastq'
        return 'fasta'

    # =========================================================================
    # Phase 2a: Adapter Scanning
    # =========================================================================

    def _scan_adapters(self, reads: List[Tuple[str, str]]) -> AdapterInfo:
        """Scan for known adapters (informational only)."""
        adapter_hits = defaultdict(list)
        reads_with_adapter = set()

        for adapter_name, adapter_seq in self.adapters:
            for read_id, sequence in reads:
                hit = self._find_adapter_in_read(sequence, adapter_seq)
                if hit:
                    adapter_hits[adapter_name].append({
                        'read_id': read_id,
                        'position': hit[0],
                        'match_length': hit[1],
                        'mismatches': hit[2]
                    })
                    reads_with_adapter.add(read_id)

        detected = []
        for adapter_name, hits in adapter_hits.items():
            frequency = len(hits) / len(reads)
            if frequency > 0.10:
                positions = [h['position'] for h in hits]
                detected.append({
                    'name': adapter_name,
                    'sequence': next(seq for name, seq in self.adapters if name == adapter_name),
                    'frequency': frequency,
                    'mean_position': sum(positions) / len(positions),
                    'position_std': self._std(positions) if len(positions) > 1 else 0.0
                })

        detected.sort(key=lambda x: x['frequency'], reverse=True)
        no_adapter_fraction = 1.0 - (len(reads_with_adapter) / len(reads))

        return AdapterInfo(
            detected_adapters=detected,
            total_reads_scanned=len(reads),
            no_adapter_fraction=no_adapter_fraction
        )

    def _find_adapter_in_read(
        self,
        sequence: str,
        adapter: str,
        min_overlap: int = 8,
        max_mismatches: int = 2
    ) -> Optional[Tuple[int, int, int]]:
        """Find adapter in read using seed-and-extend."""
        seq_len = len(sequence)
        adapter_len = len(adapter)

        for pos in range(seq_len - min_overlap + 1):
            for match_len in range(min_overlap, min(adapter_len, seq_len - pos) + 1):
                seq_part = sequence[pos:pos + match_len]
                adapter_part = adapter[:match_len]
                mismatches = sum(1 for a, b in zip(seq_part, adapter_part) if a != b)
                if mismatches <= max_mismatches:
                    return (pos, match_len, mismatches)

        return None

    # =========================================================================
    # Phase 2b: STAR Alignment
    # =========================================================================

    def _align_subset(
        self,
        reads: List[Tuple[str, str]],
        star_index: Path,
        threads: int,
        save_bam_path: Path
    ) -> Dict:
        """Align subset of reads with STAR."""
        temp_fastq = create_temp_file(suffix='.fastq')
        with open(temp_fastq, 'w') as f:
            for read_id, sequence in reads:
                f.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")

        aligner = STARAligner(
            star_index=star_index,
            threads=threads,
            max_reads=None
        )

        result = aligner.align_reads(
            input_file=temp_fastq,
            format='fastq',
            save_bam_path=save_bam_path
        )

        return {
            'total_reads': len(reads),
            'aligned_reads': result.aligned_reads,
            'alignment_rate': result.alignment_rate,
            'bam_file': str(result.output_files[0]) if result.output_files else None
        }

    # =========================================================================
    # Phase 2c: Soft-Clipping Analysis
    # =========================================================================

    def _analyze_soft_clipping(self, bam_file: Path) -> TrimBoundaries:
        """Analyze soft-clipping patterns to determine trim boundaries."""
        if not PYSAM_AVAILABLE:
            raise RuntimeError("pysam required for soft-clipping analysis")

        MIN_READS_PER_LENGTH = 100
        MIN_MODE_FREQUENCY = 0.70
        MIN_RPF_LENGTH = 20
        MAX_RPF_LENGTH = 40

        clip_5p = Counter()
        clip_3p = Counter()
        per_length_clips = defaultdict(lambda: {'5p': [], '3p': []})

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

                read_length = read.query_length
                cigar = read.cigartuples

                if not cigar:
                    continue

                # 5' soft-clipping
                if cigar[0][0] == 4:
                    clip_5p[cigar[0][1]] += 1
                    per_length_clips[read_length]['5p'].append(cigar[0][1])
                else:
                    clip_5p[0] += 1
                    per_length_clips[read_length]['5p'].append(0)

                # 3' soft-clipping
                if cigar[-1][0] == 4:
                    clip_3p[cigar[-1][1]] += 1
                    per_length_clips[read_length]['3p'].append(cigar[-1][1])
                else:
                    clip_3p[0] += 1
                    per_length_clips[read_length]['3p'].append(0)

        # Per-length consensus
        per_length_dict = {}
        high_confidence_count = 0

        for length in sorted(per_length_clips.keys()):
            clips_5p = per_length_clips[length]['5p']
            clips_3p = per_length_clips[length]['3p']
            n_reads = len(clips_5p)

            if n_reads < MIN_READS_PER_LENGTH:
                continue

            count_5p = Counter(clips_5p)
            count_3p = Counter(clips_3p)

            mode_5p, freq_5p = count_5p.most_common(1)[0]
            mode_3p, freq_3p = count_3p.most_common(1)[0]

            freq_5p_pct = freq_5p / n_reads
            freq_3p_pct = freq_3p / n_reads
            min_freq = min(freq_5p_pct, freq_3p_pct)

            rpf_length = length - mode_5p - mode_3p
            if rpf_length < MIN_RPF_LENGTH or rpf_length > MAX_RPF_LENGTH:
                continue

            if min_freq >= 0.80:
                conf_level = 'high'
                high_confidence_count += 1
            elif min_freq >= 0.70:
                conf_level = 'medium'
            else:
                conf_level = 'low'

            per_length_dict[length] = {
                'trim_5p': mode_5p,
                'trim_3p': mode_3p,
                'rpf_length': rpf_length,
                'n_reads': n_reads,
                'confidence_level': conf_level
            }

        consensus_5p = clip_5p.most_common(1)[0][0] if clip_5p else 0
        consensus_3p = clip_3p.most_common(1)[0][0] if clip_3p else 0

        if per_length_dict:
            high_conf_pct = high_confidence_count / len(per_length_dict)
            if high_conf_pct >= 0.80:
                overall_confidence = 'high'
            elif high_conf_pct >= 0.50:
                overall_confidence = 'medium'
            else:
                overall_confidence = 'low'
        else:
            overall_confidence = 'low'

        return TrimBoundaries(
            consensus_5p=consensus_5p,
            consensus_3p=consensus_3p,
            per_length=per_length_dict,
            confidence=overall_confidence
        )

    # =========================================================================
    # Phase 3: Structure Learning (THE KEY IMPROVEMENT)
    # =========================================================================

    def _learn_read_structure_reliable(self, bam_file: Path) -> LearnedStructure:
        """
        Learn read structure using per-position entropy profiles.

        This is the core improvement: we analyze entropy AT EACH POSITION
        in the soft-clipped regions, not entropy of whole sequences.

        High entropy per-position = random bases = UMI
        Low entropy per-position = conserved bases = Adapter

        CRITICAL: For 3' clips, we must align sequences from the 3' end
        (reverse them) because the adapter starts at a fixed position from
        the 3' end, but the RPF length varies. This ensures position 0 in
        our analysis is always the last base of the read.
        """
        if not PYSAM_AVAILABLE:
            raise RuntimeError("pysam required for structure learning")

        # Step 1: Extract soft-clip sequences
        # 5' clips: grouped by length, aligned from 5' end (forward)
        # 3' clips: ALL clips collected, will be aligned from 3' end (reversed)
        clips_5p_by_length = defaultdict(list)  # length -> list of sequences
        clips_3p_all = []  # All 3' clips (variable length)
        rpf_lengths = Counter()

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples:
                    continue

                seq = read.query_sequence
                cigar = read.cigartuples

                # Track aligned (RPF) length
                aligned_len = sum(length for op, length in cigar if op in (0, 7, 8))  # M, =, X
                rpf_lengths[aligned_len] += 1

                # 5' soft-clip: aligned from 5' end (position 0 = first base of read)
                if cigar[0][0] == 4:
                    clip_len = cigar[0][1]
                    clip_seq = seq[:clip_len]
                    clips_5p_by_length[clip_len].append(clip_seq)

                # 3' soft-clip: collect ALL clips (will align from 3' end)
                if cigar[-1][0] == 4:
                    clip_len = cigar[-1][1]
                    clip_seq = seq[-clip_len:]
                    clips_3p_all.append(clip_seq)

        # Step 2: Find dominant 5' clip length (need sufficient coverage)
        dominant_5p_len = self._find_dominant_length(clips_5p_by_length)

        validation_warnings = []

        # Step 3: Classify 5' region (aligned from 5' end - forward)
        five_prime = self._classify_region_5prime(
            clips_5p_by_length.get(dominant_5p_len, []),
            dominant_5p_len,
            validation_warnings
        )

        # Step 4: Classify 3' region (aligned from 3' end - REVERSED)
        three_prime = self._classify_region_3prime(
            clips_3p_all,
            validation_warnings
        )

        # Step 5: Determine overall confidence
        if five_prime.confidence > 0.8 and three_prime.confidence > 0.8:
            overall_confidence = 'high'
        elif five_prime.confidence > 0.5 or three_prime.confidence > 0.5:
            overall_confidence = 'medium'
        else:
            overall_confidence = 'low'

        # Additional validation: RPF length distribution
        if rpf_lengths:
            median_rpf = sorted(rpf_lengths.keys())[len(rpf_lengths) // 2]
            if median_rpf < 20 or median_rpf > 40:
                validation_warnings.append(f"Unusual median RPF length: {median_rpf}")
                overall_confidence = 'low'

        return LearnedStructure(
            five_prime=five_prime,
            three_prime=three_prime,
            rpf_length_distribution=dict(rpf_lengths),
            overall_confidence=overall_confidence,
            validation_warnings=validation_warnings
        )

    def _find_dominant_length(
        self,
        clips_by_length: Dict[int, List[str]]
    ) -> Optional[int]:
        """Find the most common clip length with sufficient coverage."""
        if not clips_by_length:
            return None

        # Find length with most sequences that meets minimum coverage
        best_len = None
        best_count = 0

        for length, sequences in clips_by_length.items():
            if len(sequences) >= self.MIN_COVERAGE_PER_POSITION and len(sequences) > best_count:
                best_count = len(sequences)
                best_len = length

        return best_len

    def _classify_region_5prime(
        self,
        sequences: List[str],
        length: Optional[int],
        warnings: List[str]
    ) -> RegionClassification:
        """
        Classify 5' region using per-position entropy (aligned from 5' end).

        For 5' clips, position 0 is the first base of the read.
        All sequences should have the same length (grouped by length).
        """
        region_name = "5'"

        # No clips = no region
        if not sequences or length is None or length == 0:
            return RegionClassification(
                region_type='none',
                length=0,
                confidence=1.0,
                evidence={'reason': 'no_clips_detected'}
            )

        # Insufficient coverage
        if len(sequences) < self.MIN_COVERAGE_PER_POSITION:
            warnings.append(f"{region_name} has only {len(sequences)} sequences (need {self.MIN_COVERAGE_PER_POSITION})")
            return RegionClassification(
                region_type='unknown',
                length=length,
                confidence=0.3,
                evidence={'reason': 'insufficient_coverage', 'count': len(sequences)}
            )

        # Build per-position profiles (forward direction)
        profiles = []
        for pos in range(length):
            bases = [seq[pos] for seq in sequences if len(seq) > pos]
            if len(bases) >= self.MIN_COVERAGE_PER_POSITION // 2:
                profile = PositionProfile.from_bases(pos, bases)
                profiles.append(profile)

        return self._classify_from_profiles(profiles, length, region_name, sequences, warnings)

    def _classify_region_3prime(
        self,
        sequences: List[str],
        warnings: List[str]
    ) -> RegionClassification:
        """
        Classify 3' region by trying BOTH orientations.

        For 3' soft-clips, the structure could be:
        1. ADAPTER at start of clip (most common): [RPF][ADAPTER...]
           - Align from 5' end of clip (forward)
        2. UMI at end of read: [RPF][...][UMI]
           - Align from 3' end of clip (reversed)

        We try both and pick the one with clearer signal (more conserved or
        more random positions, not ambiguous).
        """
        region_name = "3'"

        # No clips = no region
        if not sequences:
            return RegionClassification(
                region_type='none',
                length=0,
                confidence=1.0,
                evidence={'reason': 'no_clips_detected'}
            )

        # Insufficient coverage
        if len(sequences) < self.MIN_COVERAGE_PER_POSITION:
            warnings.append(f"{region_name} has only {len(sequences)} sequences (need {self.MIN_COVERAGE_PER_POSITION})")
            return RegionClassification(
                region_type='unknown',
                length=0,
                confidence=0.3,
                evidence={'reason': 'insufficient_coverage', 'count': len(sequences)}
            )

        # Find the dominant length for reporting
        length_counts = Counter(len(s) for s in sequences)
        dominant_length = length_counts.most_common(1)[0][0]

        # TRY FORWARD (5' end of clip): Assumes adapter at RPF boundary
        # This is most common for Ribo-seq: [RPF][ADAPTER]
        forward_profiles = self._build_profiles_forward(sequences)
        forward_result = self._classify_from_profiles(
            forward_profiles, dominant_length, f"{region_name} (forward)", sequences, []
        )

        # TRY REVERSED (3' end of read): Assumes UMI at 3' end
        # Less common but possible: [RPF][linker][UMI]
        reversed_seqs = [seq[::-1] for seq in sequences]
        reversed_profiles = self._build_profiles_forward(reversed_seqs)
        reversed_result = self._classify_from_profiles(
            reversed_profiles, dominant_length, f"{region_name} (reversed)", reversed_seqs, []
        )

        # Pick the better result based on confidence and clarity
        # Prefer adapter detection (most common) if both are similar
        logger.debug(f"  3' forward: {forward_result.region_type} (conf={forward_result.confidence:.2f})")
        logger.debug(f"  3' reversed: {reversed_result.region_type} (conf={reversed_result.confidence:.2f})")

        # Decision logic:
        # 1. If one is clearly adapter and other is unknown/low-conf, use the adapter
        # 2. If one is clearly UMI and other is unknown, use the UMI
        # 3. If both similar, prefer forward (adapter at boundary is more common)

        if forward_result.region_type == 'adapter' and forward_result.confidence > 0.5:
            # Adapter detected at RPF boundary - most common case
            return forward_result
        elif reversed_result.region_type == 'umi' and reversed_result.confidence > 0.7:
            # UMI at 3' end of read - need to un-reverse consensus if any
            return reversed_result
        elif forward_result.confidence >= reversed_result.confidence:
            return forward_result
        else:
            # Reversed was better - need to un-reverse any consensus
            if reversed_result.consensus_sequence:
                return RegionClassification(
                    region_type=reversed_result.region_type,
                    length=reversed_result.length,
                    confidence=reversed_result.confidence,
                    consensus_sequence=reversed_result.consensus_sequence[::-1],
                    adapter_name=reversed_result.adapter_name,
                    evidence=reversed_result.evidence
                )
            return reversed_result

    def _build_profiles_forward(self, sequences: List[str]) -> List[PositionProfile]:
        """Build per-position profiles aligned from start of sequences."""
        max_len = max(len(s) for s in sequences) if sequences else 0
        profiles = []
        for pos in range(min(max_len, 30)):  # Limit to first 30 positions
            bases = [seq[pos] for seq in sequences if len(seq) > pos]
            if len(bases) >= self.MIN_COVERAGE_PER_POSITION // 2:
                profile = PositionProfile.from_bases(pos, bases)
                profiles.append(profile)
        return profiles

    def _classify_from_profiles(
        self,
        profiles: List[PositionProfile],
        length: int,
        region_name: str,
        sequences: List[str],
        warnings: List[str]
    ) -> RegionClassification:
        """
        Classify a region based on per-position entropy profiles.

        Shared logic for both 5' and 3' classification.
        """
        if not profiles:
            warnings.append(f"{region_name} has no valid position profiles")
            return RegionClassification(
                region_type='unknown',
                length=length,
                confidence=0.2,
                evidence={'reason': 'no_valid_profiles'}
            )

        # Analyze entropy distribution
        entropies = [p.entropy for p in profiles]
        avg_entropy = sum(entropies) / len(entropies)
        random_positions = sum(1 for p in profiles if p.is_random)
        conserved_positions = sum(1 for p in profiles if p.is_conserved)
        random_fraction = random_positions / len(profiles)
        conserved_fraction = conserved_positions / len(profiles)

        # Classify based on entropy pattern
        if random_fraction > 0.7:
            # Most positions are random -> UMI
            confidence = min(0.5 + random_fraction / 2, 0.95)
            return RegionClassification(
                region_type='umi',
                length=length,
                confidence=confidence,
                evidence={
                    'avg_entropy': avg_entropy,
                    'random_fraction': random_fraction,
                    'conserved_fraction': conserved_fraction,
                    'sample_sequences': sequences[:20]
                }
            )

        elif conserved_fraction > 0.5:
            # Most positions are conserved -> Adapter
            consensus = ''.join(p.dominant_base for p in profiles)

            # Validate against known adapter database
            adapter_match = self._match_adapter_database(consensus)
            if adapter_match:
                confidence = min(0.6 + conserved_fraction / 2.5, 0.95)
                return RegionClassification(
                    region_type='adapter',
                    length=length,
                    confidence=confidence,
                    consensus_sequence=consensus,
                    adapter_name=adapter_match[0],
                    evidence={
                        'avg_entropy': avg_entropy,
                        'random_fraction': random_fraction,
                        'conserved_fraction': conserved_fraction,
                        'database_match': adapter_match[0],
                        'match_score': adapter_match[1]
                    }
                )
            else:
                # Conserved but not in database - warn but accept
                warnings.append(f"{region_name} adapter not found in database: {consensus[:20]}...")
                confidence = conserved_fraction * 0.7  # Lower confidence
                return RegionClassification(
                    region_type='adapter',
                    length=length,
                    confidence=confidence,
                    consensus_sequence=consensus,
                    evidence={
                        'avg_entropy': avg_entropy,
                        'random_fraction': random_fraction,
                        'conserved_fraction': conserved_fraction,
                        'database_match': None
                    }
                )

        else:
            # Ambiguous entropy pattern
            warnings.append(f"{region_name} has ambiguous entropy (random={random_fraction:.1%}, conserved={conserved_fraction:.1%})")
            return RegionClassification(
                region_type='unknown',
                length=length,
                confidence=0.4,
                evidence={
                    'avg_entropy': avg_entropy,
                    'random_fraction': random_fraction,
                    'conserved_fraction': conserved_fraction
                }
            )

    def _match_adapter_database(
        self,
        consensus: str
    ) -> Optional[Tuple[str, float]]:
        """Match consensus sequence against known adapter database."""
        best_match = None
        best_score = 0.0

        for adapter_name, adapter_seq in self.adapters:
            score = self._sequence_identity(consensus, adapter_seq)
            if score > best_score and score >= self.ADAPTER_MATCH_THRESHOLD:
                best_score = score
                best_match = (adapter_name, score)

        return best_match

    def _sequence_identity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence identity (allowing partial overlap)."""
        if not seq1 or not seq2:
            return 0.0

        # Try aligning seq1 at different positions of seq2
        best_identity = 0.0

        for offset in range(-len(seq1) + 5, len(seq2) - 5):
            matches = 0
            compared = 0

            for i, base in enumerate(seq1):
                j = i + offset
                if 0 <= j < len(seq2):
                    compared += 1
                    if base == seq2[j]:
                        matches += 1

            if compared >= 5:
                identity = matches / compared
                if identity > best_identity:
                    best_identity = identity

        return best_identity

    def _log_learned_structure(self, structure: LearnedStructure):
        """Log details of learned structure."""
        logger.info(f"  Learned Structure ({structure.overall_confidence} confidence):")

        # 5' region
        if structure.five_prime:
            fp = structure.five_prime
            if fp.region_type == 'umi':
                logger.info(f"    5' End: UMI (length={fp.length}bp, entropy={fp.evidence.get('avg_entropy', 0):.2f})")
            elif fp.region_type == 'adapter':
                logger.info(f"    5' End: ADAPTER ({fp.adapter_name or 'unknown'})")
                logger.info(f"             Sequence: {fp.consensus_sequence}")
            elif fp.region_type == 'none':
                logger.info(f"    5' End: None (no soft-clipping)")
            else:
                logger.info(f"    5' End: UNKNOWN (ambiguous)")

        # 3' region
        if structure.three_prime:
            tp = structure.three_prime
            if tp.region_type == 'umi':
                logger.info(f"    3' End: UMI (length={tp.length}bp)")
            elif tp.region_type == 'adapter':
                logger.info(f"    3' End: ADAPTER ({tp.adapter_name or 'unknown'})")
                logger.info(f"             Sequence: {tp.consensus_sequence}")
            elif tp.region_type == 'none':
                logger.info(f"    3' End: None (no soft-clipping)")
            else:
                logger.info(f"    3' End: UNKNOWN (ambiguous)")

        # Warnings
        for warning in structure.validation_warnings:
            logger.warning(f"    Warning: {warning}")

    # =========================================================================
    # Phase 4: Extraction with Learned Configuration
    # =========================================================================

    def _extract_with_config(
        self,
        input_file: Path,
        output_file: Path,
        config: TrimmerConfig,
        preserve_umi: bool
    ) -> Dict:
        """
        Extract RPFs by applying learned structure configuration.

        This is the clean, correct implementation with proper index tracking.
        No index mismatch bugs!
        """
        extracted_rpfs = 0
        total_reads = 0
        too_short = 0
        too_long = 0
        adapter_found = 0

        format_type = self._detect_format(input_file)
        opener = get_file_opener(input_file)

        with opener(input_file, "rt") as f_in, open(output_file, "w") as f_out:
            for record in SeqIO.parse(f_in, format_type):
                total_reads += 1

                seq = str(record.seq)
                qual = list(record.letter_annotations.get("phred_quality", []))

                # Track indices for clean slicing
                start_idx = 0
                end_idx = len(seq)
                umi_seq = None

                # Phase A: 5' Trimming (fixed length for UMI or adapter)
                if config.trim_5p_fixed and config.trim_5p_fixed > 0:
                    trim_5p = min(config.trim_5p_fixed, len(seq))
                    if preserve_umi and config.trim_5p_is_umi:
                        umi_seq = seq[:trim_5p]
                    start_idx = trim_5p

                # Phase B: 3' Trimming (search for adapter)
                if config.trim_3p_adapter:
                    # Search for adapter in the REMAINING sequence (after 5' trim)
                    remaining_seq = seq[start_idx:]
                    adapter_pos = self._find_adapter_best_match(
                        remaining_seq,
                        config.trim_3p_adapter,
                        config.trim_3p_adapter_min_overlap,
                        config.trim_3p_adapter_max_mismatches
                    )
                    if adapter_pos >= 0:
                        # Adjust to original coordinates
                        end_idx = start_idx + adapter_pos
                        adapter_found += 1

                # Phase C: Extract RPF
                rpf_seq = seq[start_idx:end_idx]

                # Length filter
                if len(rpf_seq) < config.min_rpf_length:
                    too_short += 1
                    continue
                if len(rpf_seq) > config.max_rpf_length:
                    too_long += 1
                    continue

                # Phase D: Handle quality scores (SAME INDICES - no mismatch!)
                if qual:
                    rpf_qual = qual[start_idx:end_idx]
                    qual_str = ''.join(chr(q + 33) for q in rpf_qual)
                else:
                    qual_str = 'I' * len(rpf_seq)

                # Phase E: Write output
                header = record.id
                if umi_seq:
                    header = f"{header} UMI:{umi_seq}"

                f_out.write(f"@{header}\n{rpf_seq}\n+\n{qual_str}\n")
                extracted_rpfs += 1

        logger.debug(f"  Adapter found in {adapter_found} reads, "
                    f"{too_short} too short, {too_long} too long")

        return {
            'total_reads': total_reads,
            'extracted_rpfs': extracted_rpfs,
            'extraction_rate': extracted_rpfs / total_reads if total_reads > 0 else 0,
            'adapter_found': adapter_found,
            'too_short': too_short,
            'too_long': too_long
        }

    def _find_adapter_best_match(
        self,
        sequence: str,
        adapter: str,
        min_overlap: int,
        max_mismatches: int
    ) -> int:
        """
        Find adapter position using BEST (longest) match, not first match.

        Returns position of adapter start, or -1 if not found.
        """
        best_pos = -1
        best_length = 0

        seq_len = len(sequence)
        adapter_len = len(adapter)

        if seq_len < min_overlap:
            return -1

        for pos in range(seq_len - min_overlap + 1):
            max_match_len = min(adapter_len, seq_len - pos)

            # Try longest match first (greedy from long to short)
            for match_len in range(max_match_len, min_overlap - 1, -1):
                seq_part = sequence[pos:pos + match_len]
                adapter_part = adapter[:match_len]

                mismatches = sum(1 for a, b in zip(seq_part, adapter_part) if a != b)

                if mismatches <= max_mismatches:
                    if match_len > best_length:
                        best_length = match_len
                        best_pos = pos
                    break  # Found best match at this position

        return best_pos

    # =========================================================================
    # Seqspec Export
    # =========================================================================

    def _export_seqspec(self, structure: LearnedStructure, output_path: Path):
        """Export learned structure as seqspec YAML."""
        seqspec_data = structure.to_seqspec_yaml()

        with open(output_path, 'w') as f:
            yaml.dump(seqspec_data, f, default_flow_style=False, sort_keys=False)

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def _std(self, values: List[float]) -> float:
        """Calculate standard deviation."""
        if len(values) < 2:
            return 0.0
        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
        return math.sqrt(variance)
