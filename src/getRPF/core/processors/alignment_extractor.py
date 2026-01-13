"""Alignment-based RPF extraction using STAR as ground truth.

This is the primary extraction method for getRPF. It uses reference genome
alignment to determine biological sequence boundaries, with optional adapter
detection for reporting.

Algorithm:
1. Sample subset of reads (default 10k)
2. Scan for known adapters (reporting only)
3. Align subset with STAR
4. Analyze soft-clipping patterns to determine trim boundaries
5. Extract RPF from all reads using consensus boundaries
6. Optionally detect UMI in soft-clipped regions

Author: getRPF team
"""

import logging
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

from Bio import SeqIO

from .alignment import STARAligner
from ...utils.file_utils import create_temp_file

logger = logging.getLogger(__name__)


# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class AdapterHit:
    """Single adapter detection hit."""
    adapter_sequence: str
    adapter_name: str
    read_id: str
    position: int
    match_length: int
    mismatches: int


@dataclass
class AdapterInfo:
    """Adapter detection results for reporting."""
    detected_adapters: List[Dict]  # List of {sequence, name, frequency, position}
    total_reads_scanned: int
    no_adapter_fraction: float

    @property
    def top_adapter(self) -> Optional[str]:
        """Return most common adapter name."""
        if self.detected_adapters:
            return self.detected_adapters[0]['name']
        return None

    def to_dict(self) -> Dict:
        """Convert to dictionary for reporting."""
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
    per_length: Dict[int, Dict[str, int]]  # {length: {trim_5p, trim_3p}}
    confidence: str  # 'high', 'medium', 'low'

    def to_dict(self) -> Dict:
        """Convert to dictionary for reporting."""
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
    positions: Optional[Tuple[int, int]]  # (start, end)
    entropy: float
    confidence: str  # 'high', 'medium', 'low', 'none'
    sample_sequences: List[str]

    def to_dict(self) -> Dict:
        """Convert to dictionary for reporting."""
        return {
            'umi_detected': self.detected,
            'positions': list(self.positions) if self.positions else None,
            'entropy': self.entropy,
            'confidence': self.confidence,
            'sample_sequences': self.sample_sequences[:10]  # Limit to 10 examples
        }


@dataclass
class ExtractionResult:
    """Complete extraction results."""
    input_reads: int
    extracted_rpfs: int
    extraction_rate: float
    sample_size: int
    trim_boundaries: TrimBoundaries
    adapter_info: Optional[AdapterInfo]
    umi_info: Optional[UMIInfo]
    alignment_stats: Dict
    method: str = "alignment_based"

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON reporting."""
        return {
            'extraction_summary': {
                'input_reads': self.input_reads,
                'extracted_rpfs': self.extracted_rpfs,
                'extraction_rate': self.extraction_rate,
                'sample_size': self.sample_size,
                'method': self.method
            },
            'trim_boundaries': self.trim_boundaries.to_dict(),
            'adapter_scan': self.adapter_info.to_dict() if self.adapter_info else None,
            'umi_detection': self.umi_info.to_dict() if self.umi_info else None,
            'alignment_statistics': self.alignment_stats
        }


# ============================================================================
# Adapter Database
# ============================================================================

KNOWN_ADAPTERS = [
    # Illumina adapters (most common)
    ("Illumina TruSeq Universal", "AGATCGGAAGAGCACACGTCT"),
    ("Illumina Small RNA 3'", "TGGAATTCTCGGGTGCCAAGG"),
    ("Illumina TruSeq (partial)", "GATCGGAAGAGCACACGT"),
    ("Illumina Universal (short)", "AGATCGGAAGAGC"),
    ("Illumina Multiplexing", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),

    # Ingolia lab protocols
    ("Ingolia 2009", "CTGTAGGCACCATCAAT"),
    ("McGlincy 2017 (with barcode)", "AGATCGGAAGAGCACACGTCTGAA"),

    # Nextera
    ("Nextera Transposase", "CTGTCTCTTATACACATCT"),

    # Other common
    ("PolyA tail", "AAAAAAAAAA"),
    ("NEBNext Small RNA", "AGATCGGAAGAGCACACGTCT"),
    ("QIAseq miRNA", "AACTGTAGGCACCATCAAT"),

    # Additional partial sequences for detection
    ("Generic Illumina 1", "GATCGGAAGAGCGTCGT"),
    ("Generic Illumina 2", "GATCGGAAGAGCTCGTA"),
    ("Generic Illumina 3", "TGATCGGAAGAGCACAC"),
]


# ============================================================================
# Main Extractor Class
# ============================================================================

class AlignmentBasedExtractor:
    """
    Primary RPF extraction using alignment as ground truth.

    This extractor uses reference genome alignment to determine where
    biological sequence (RPF) ends and synthetic sequences (adapters,
    UMIs) begin. It's more reliable than pattern matching or HMM
    segmentation because the reference genome is ground truth.
    """

    def __init__(self, adapters: Optional[List[Tuple[str, str]]] = None):
        """
        Initialize extractor.

        Args:
            adapters: List of (name, sequence) tuples for adapter detection.
                     If None, uses built-in KNOWN_ADAPTERS.
        """
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
        Extract RPF sequences using alignment-based approach.

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
        """
        logger.info("=" * 60)
        logger.info("Starting Alignment-Based RPF Extraction")
        logger.info("=" * 60)

        # Phase 1: Sample reads
        logger.info(f"Phase 1/3: Sampling {sample_size} reads for characterization...")
        subset_reads = self._sample_reads(input_file, sample_size)
        logger.info(f"  → Sampled {len(subset_reads)} reads")

        # Phase 2: Characterization (parallel operations)
        logger.info("Phase 2/3: Characterizing sample...")

        # 2a. Adapter scan (fast, informational only)
        adapter_info = None
        if report_adapters:
            logger.info("  → Scanning for known adapters...")
            adapter_info = self._scan_adapters(subset_reads)
            if adapter_info.top_adapter:
                top = adapter_info.detected_adapters[0]
                logger.info(f"    ✓ Detected: {top['name']} ({top['frequency']:.1%} of reads)")
            else:
                logger.info(f"    ✓ No known adapters detected")

        # 2b. STAR alignment (ground truth)
        logger.info("  → Aligning subset with STAR...")
        alignment_result = self._align_subset(subset_reads, star_index, star_threads)
        logger.info(f"    ✓ Aligned: {alignment_result['aligned_reads']}/{alignment_result['total_reads']} "
                   f"({alignment_result['alignment_rate']:.1%})")

        # 2c. Analyze soft-clipping
        logger.info("  → Analyzing soft-clipping patterns...")
        boundaries = self._analyze_soft_clipping(alignment_result['bam_file'])
        logger.info(f"    ✓ Consensus boundaries: 5'={boundaries.consensus_5p}nt, "
                   f"3'={boundaries.consensus_3p}nt (confidence: {boundaries.confidence})")

        # Phase 3: Extraction
        logger.info("Phase 3/3: Extracting RPF from all reads...")
        extraction_stats = self._extract_rpf_all(
            input_file,
            output_file,
            boundaries,
            preserve_umi
        )
        logger.info(f"  → Extracted {extraction_stats['extracted_rpfs']} RPFs from "
                   f"{extraction_stats['total_reads']} reads ({extraction_stats['extraction_rate']:.1%})")

        # Optional: UMI detection
        umi_info = None
        if preserve_umi:
            logger.info("  → Detecting UMI in soft-clipped regions...")
            umi_info = self._detect_umi(subset_reads, boundaries)
            if umi_info.detected:
                logger.info(f"    ✓ UMI detected: positions {umi_info.positions}, "
                           f"entropy={umi_info.entropy:.2f} (confidence: {umi_info.confidence})")
            else:
                logger.info(f"    ✗ No UMI detected (entropy too low)")

        logger.info("=" * 60)
        logger.info("Extraction Complete!")
        logger.info("=" * 60)

        return ExtractionResult(
            input_reads=extraction_stats['total_reads'],
            extracted_rpfs=extraction_stats['extracted_rpfs'],
            extraction_rate=extraction_stats['extraction_rate'],
            sample_size=len(subset_reads),
            trim_boundaries=boundaries,
            adapter_info=adapter_info,
            umi_info=umi_info,
            alignment_stats=alignment_result,
            method="alignment_based"
        )

    # ========================================================================
    # Phase 1: Sampling
    # ========================================================================

    def _sample_reads(self, input_file: Path, sample_size: int) -> List[Tuple[str, str]]:
        """
        Sample reads from input file.

        Args:
            input_file: Path to input file
            sample_size: Number of reads to sample

        Returns:
            List of (read_id, sequence) tuples
        """
        reads = []

        # Detect format
        if input_file.suffix in ['.fq', '.fastq']:
            format_type = 'fastq'
        else:
            format_type = 'fasta'

        with open(input_file) as f:
            for i, record in enumerate(SeqIO.parse(f, format_type)):
                if i >= sample_size:
                    break
                reads.append((record.id, str(record.seq)))

        return reads

    # ========================================================================
    # Phase 2a: Adapter Scanning
    # ========================================================================

    def _scan_adapters(self, reads: List[Tuple[str, str]]) -> AdapterInfo:
        """
        Scan for known adapters using k-mer search.

        This is INFORMATIONAL ONLY - we don't use adapter detection
        for trimming, only for reporting what protocol was used.

        Args:
            reads: List of (read_id, sequence) tuples

        Returns:
            AdapterInfo with detection results
        """
        adapter_hits = defaultdict(list)
        reads_with_adapter = set()

        # Search for each adapter
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

        # Compile results
        detected = []
        for adapter_name, hits in adapter_hits.items():
            frequency = len(hits) / len(reads)
            if frequency > 0.10:  # 10% threshold
                positions = [h['position'] for h in hits]
                detected.append({
                    'name': adapter_name,
                    'sequence': next(seq for name, seq in self.adapters if name == adapter_name),
                    'frequency': frequency,
                    'mean_position': sum(positions) / len(positions),
                    'position_std': self._std(positions) if len(positions) > 1 else 0.0
                })

        # Sort by frequency
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
        """
        Find adapter in read using seed-and-extend.

        Args:
            sequence: Read sequence
            adapter: Adapter sequence
            min_overlap: Minimum overlap length
            max_mismatches: Maximum mismatches allowed

        Returns:
            (position, match_length, mismatches) or None
        """
        seq_len = len(sequence)
        adapter_len = len(adapter)

        # Search for adapter starting at each position
        for pos in range(seq_len - min_overlap + 1):
            # Try different match lengths
            for match_len in range(min_overlap, min(adapter_len, seq_len - pos) + 1):
                seq_part = sequence[pos:pos + match_len]
                adapter_part = adapter[:match_len]

                mismatches = sum(1 for a, b in zip(seq_part, adapter_part) if a != b)

                if mismatches <= max_mismatches:
                    return (pos, match_len, mismatches)

        return None

    # ========================================================================
    # Phase 2b: STAR Alignment
    # ========================================================================

    def _align_subset(
        self,
        reads: List[Tuple[str, str]],
        star_index: Path,
        threads: int
    ) -> Dict:
        """
        Align subset of reads with STAR.

        Args:
            reads: List of (read_id, sequence) tuples
            star_index: Path to STAR index
            threads: Number of threads

        Returns:
            Dict with alignment stats and BAM file path
        """
        # Write reads to temp FASTQ
        temp_fastq = create_temp_file(suffix='.fastq')
        with open(temp_fastq, 'w') as f:
            for read_id, sequence in reads:
                f.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")

        # Run STAR alignment
        aligner = STARAligner(
            star_index=star_index,
            threads=threads,
            max_reads=None
        )

        result = aligner.align_reads(
            input_file=temp_fastq,
            format='fastq',
            save_bam_path=None  # Uses temp file
        )

        return {
            'total_reads': len(reads),
            'aligned_reads': result.aligned_reads,
            'alignment_rate': result.alignment_rate,
            'bam_file': result.output_files[0] if result.output_files else None
        }

    # ========================================================================
    # Phase 2c: Soft-Clipping Analysis
    # ========================================================================

    def _analyze_soft_clipping(self, bam_file: Path) -> TrimBoundaries:
        """
        Analyze soft-clipping patterns to determine trim boundaries.

        Uses mode (most common) clipping position per read length with
        strict confidence requirements for single-nucleotide precision.

        Soft-clipped bases = non-biological sequence (UMI, adapter, linker)
        Aligned bases = biological RPF sequence

        Args:
            bam_file: Path to BAM file

        Returns:
            TrimBoundaries with consensus positions and confidence metrics
        """
        if not PYSAM_AVAILABLE:
            raise RuntimeError("pysam required for soft-clipping analysis")

        MIN_READS_PER_LENGTH = 100  # Require 100 reads for robust mode
        MIN_MODE_FREQUENCY = 0.70   # Mode must be 70%+ for high confidence
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
                if cigar[0][0] == 4:  # 4 = soft clip
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

        # Per-length consensus with confidence metrics
        per_length_dict = {}
        low_confidence_lengths = []
        high_confidence_count = 0

        for length in sorted(per_length_clips.keys()):
            clips_5p = per_length_clips[length]['5p']
            clips_3p = per_length_clips[length]['3p']
            n_reads = len(clips_5p)

            # Skip lengths with insufficient reads
            if n_reads < MIN_READS_PER_LENGTH:
                logger.warning(f"Length {length}nt: Only {n_reads} reads "
                             f"(need ≥{MIN_READS_PER_LENGTH}) - skipping")
                low_confidence_lengths.append(length)
                continue

            # Calculate mode and frequency
            count_5p = Counter(clips_5p)
            count_3p = Counter(clips_3p)

            mode_5p, freq_5p = count_5p.most_common(1)[0]
            mode_3p, freq_3p = count_3p.most_common(1)[0]

            freq_5p_pct = freq_5p / n_reads
            freq_3p_pct = freq_3p / n_reads
            min_freq = min(freq_5p_pct, freq_3p_pct)

            # Check mode frequency confidence
            if freq_5p_pct < MIN_MODE_FREQUENCY or freq_3p_pct < MIN_MODE_FREQUENCY:
                logger.warning(f"Length {length}nt: Low mode consensus "
                             f"(5': {freq_5p_pct:.1%}, 3': {freq_3p_pct:.1%}) "
                             f"- possible ±1nt ambiguity")
                low_confidence_lengths.append(length)

            # Biological sanity check on RPF length
            rpf_length = length - mode_5p - mode_3p

            if rpf_length < MIN_RPF_LENGTH or rpf_length > MAX_RPF_LENGTH:
                logger.error(f"Length {length}nt: Invalid RPF length {rpf_length}nt "
                           f"(trim_5p={mode_5p}, trim_3p={mode_3p})")
                logger.error(f"Expected RPF: {MIN_RPF_LENGTH}-{MAX_RPF_LENGTH}nt")
                low_confidence_lengths.append(length)
                continue

            # Determine confidence level
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
                'confidence_5p': freq_5p_pct,
                'confidence_3p': freq_3p_pct,
                'consensus': min_freq,
                'confidence_level': conf_level
            }

            logger.info(f"Length {length}nt: 5'={mode_5p}nt, 3'={mode_3p}nt, "
                       f"RPF={rpf_length}nt, consensus={min_freq:.1%} ({conf_level})")

        # Global consensus (fallback for missing lengths)
        consensus_5p = clip_5p.most_common(1)[0][0] if clip_5p else 0
        consensus_3p = clip_3p.most_common(1)[0][0] if clip_3p else 0

        # Overall confidence assessment
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
            logger.error("No read lengths met minimum requirements!")

        # Report low confidence lengths
        if low_confidence_lengths:
            logger.warning(f"Low confidence lengths: {low_confidence_lengths}")
            logger.warning("These may have ±1nt precision errors")

        return TrimBoundaries(
            consensus_5p=consensus_5p,
            consensus_3p=consensus_3p,
            per_length=per_length_dict,
            confidence=overall_confidence
        )

    # ========================================================================
    # Phase 3: Extraction
    # ========================================================================

    def _extract_rpf_all(
        self,
        input_file: Path,
        output_file: Path,
        boundaries: TrimBoundaries,
        preserve_umi: bool
    ) -> Dict:
        """
        Extract RPF from all reads using consensus boundaries.

        Args:
            input_file: Input file path
            output_file: Output file path
            boundaries: Trim boundaries from alignment
            preserve_umi: Whether to add UMI to headers

        Returns:
            Dict with extraction statistics
        """
        total_reads = 0
        extracted_rpfs = 0

        # Detect format
        if input_file.suffix in ['.fq', '.fastq']:
            format_type = 'fastq'
        else:
            format_type = 'fasta'

        with open(input_file) as f_in, open(output_file, 'w') as f_out:
            for record in SeqIO.parse(f_in, format_type):
                total_reads += 1

                seq = str(record.seq)
                read_len = len(seq)

                # Use per-length boundaries if available, else consensus
                if read_len in boundaries.per_length:
                    trim_5p = boundaries.per_length[read_len]['trim_5p']
                    trim_3p = boundaries.per_length[read_len]['trim_3p']
                else:
                    trim_5p = boundaries.consensus_5p
                    trim_3p = boundaries.consensus_3p

                # Extract RPF
                rpf_end = read_len - trim_3p if trim_3p > 0 else read_len
                rpf_seq = seq[trim_5p:rpf_end]

                # Skip if too short
                if len(rpf_seq) < 20:
                    continue

                # Write output
                if format_type == 'fastq':
                    qual = record.letter_annotations['phred_quality']
                    rpf_qual = qual[trim_5p:rpf_end]

                    # Format: @read_id UMI:sequence (if preserve_umi)
                    header = record.id
                    if preserve_umi and trim_5p > 0:
                        umi_seq = seq[:trim_5p]
                        header = f"{record.id} UMI:{umi_seq}"

                    f_out.write(f"@{header}\n{rpf_seq}\n+\n")
                    f_out.write(''.join(chr(q + 33) for q in rpf_qual) + '\n')
                else:
                    f_out.write(f">{record.id}\n{rpf_seq}\n")

                extracted_rpfs += 1

        return {
            'total_reads': total_reads,
            'extracted_rpfs': extracted_rpfs,
            'extraction_rate': extracted_rpfs / total_reads if total_reads > 0 else 0
        }

    # ========================================================================
    # UMI Detection
    # ========================================================================

    def _detect_umi(
        self,
        reads: List[Tuple[str, str]],
        boundaries: TrimBoundaries
    ) -> UMIInfo:
        """
        Detect UMI in soft-clipped 5' regions.

        UMI signature: High entropy (nearly random nucleotide distribution)

        Args:
            reads: List of (read_id, sequence) tuples
            boundaries: Trim boundaries

        Returns:
            UMIInfo with detection results
        """
        trim_5p = boundaries.consensus_5p

        if trim_5p == 0:
            return UMIInfo(
                detected=False,
                positions=None,
                entropy=0.0,
                confidence='none',
                sample_sequences=[]
            )

        # Extract 5' clipped regions
        clipped_seqs = [seq[:trim_5p] for _, seq in reads]

        # Calculate per-position entropy
        entropies = []
        for pos in range(trim_5p):
            bases = [seq[pos] for seq in clipped_seqs if len(seq) > pos]
            entropy = self._calculate_shannon_entropy(bases)
            entropies.append(entropy)

        avg_entropy = sum(entropies) / len(entropies) if entropies else 0.0

        # Determine if it's a UMI based on entropy
        if avg_entropy > 1.8:
            detected = True
            confidence = 'high'
        elif avg_entropy > 1.5:
            detected = True
            confidence = 'medium'
        else:
            detected = False
            confidence = 'low'

        return UMIInfo(
            detected=detected,
            positions=(0, trim_5p) if detected else None,
            entropy=avg_entropy,
            confidence=confidence,
            sample_sequences=clipped_seqs[:20]  # First 20 examples
        )

    def _calculate_shannon_entropy(self, bases: List[str]) -> float:
        """Calculate Shannon entropy for a list of bases."""
        if not bases:
            return 0.0

        counts = Counter(bases)
        total = len(bases)
        entropy = 0.0

        for count in counts.values():
            p = count / total
            entropy -= p * math.log2(p)

        return entropy

    def _std(self, values: List[float]) -> float:
        """Calculate standard deviation."""
        if len(values) < 2:
            return 0.0
        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
        return math.sqrt(variance)
