"""Shared types and dataclasses for RPF processing."""

from dataclasses import dataclass, asdict, field
from typing import List, Tuple, Dict, Optional, Any, Literal
from pathlib import Path
import json
import math


# =============================================================================
# Structure Learning Types (New Reliable Algorithm)
# =============================================================================

@dataclass
class PositionProfile:
    """Per-position nucleotide statistics for structure learning.

    This is the key data structure for reliable UMI/adapter detection.
    We compute entropy and composition AT EACH POSITION, not per-sequence.
    """
    position: int
    base_counts: Dict[str, int]  # {'A': 500, 'C': 300, 'G': 150, 'T': 50}
    total: int
    entropy: float  # Shannon entropy in bits (0-2)
    dominant_base: str
    dominant_freq: float

    @property
    def is_random(self) -> bool:
        """High entropy indicates random sequence (UMI-like)."""
        return self.entropy > 1.7

    @property
    def is_conserved(self) -> bool:
        """Low entropy with high dominant frequency indicates conserved (adapter-like)."""
        return self.entropy < 0.5 and self.dominant_freq > 0.8

    @classmethod
    def from_bases(cls, position: int, bases: List[str]) -> 'PositionProfile':
        """Compute profile from list of bases at this position."""
        from collections import Counter
        counts = Counter(bases)
        total = len(bases)

        if total == 0:
            return cls(
                position=position,
                base_counts={},
                total=0,
                entropy=0.0,
                dominant_base='N',
                dominant_freq=0.0
            )

        # Shannon entropy
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * math.log2(p)

        # Dominant base
        dominant_base, dominant_count = counts.most_common(1)[0]
        dominant_freq = dominant_count / total

        return cls(
            position=position,
            base_counts=dict(counts),
            total=total,
            entropy=entropy,
            dominant_base=dominant_base,
            dominant_freq=dominant_freq
        )


@dataclass
class RegionClassification:
    """Classification of a contiguous region (5' or 3' soft-clipped region).

    This represents what we learned about one end of the read structure.
    """
    region_type: Literal['umi', 'adapter', 'none', 'unknown']
    length: int
    confidence: float  # 0-1, statistically calibrated
    consensus_sequence: Optional[str] = None  # For adapters
    adapter_name: Optional[str] = None  # Matched adapter from database
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'region_type': self.region_type,
            'length': self.length,
            'confidence': self.confidence,
            'consensus_sequence': self.consensus_sequence,
            'adapter_name': self.adapter_name,
            'evidence': self.evidence
        }


@dataclass
class LearnedStructure:
    """Complete learned library structure from alignment analysis.

    This is the output of the structure learning phase, capturing
    everything we know about the read architecture.
    """
    five_prime: Optional[RegionClassification]
    three_prime: Optional[RegionClassification]
    rpf_length_distribution: Dict[int, int]  # length -> count
    overall_confidence: Literal['high', 'medium', 'low']
    validation_warnings: List[str] = field(default_factory=list)

    def to_trimmer_config(self) -> 'TrimmerConfig':
        """Convert learned structure to actionable trimmer configuration."""
        # 5' trimming
        trim_5p_fixed = None
        trim_5p_is_umi = False
        if self.five_prime and self.five_prime.region_type in ('umi', 'adapter'):
            trim_5p_fixed = self.five_prime.length
            trim_5p_is_umi = (self.five_prime.region_type == 'umi')

        # 3' trimming
        trim_3p_adapter = None
        if self.three_prime and self.three_prime.region_type == 'adapter':
            trim_3p_adapter = self.three_prime.consensus_sequence

        return TrimmerConfig(
            trim_5p_fixed=trim_5p_fixed,
            trim_5p_is_umi=trim_5p_is_umi,
            trim_3p_adapter=trim_3p_adapter,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'five_prime': self.five_prime.to_dict() if self.five_prime else None,
            'three_prime': self.three_prime.to_dict() if self.three_prime else None,
            'rpf_length_distribution': self.rpf_length_distribution,
            'overall_confidence': self.overall_confidence,
            'validation_warnings': self.validation_warnings
        }

    def to_seqspec_yaml(self) -> Dict[str, Any]:
        """Generate seqspec-compatible YAML structure."""
        sequence_spec = []

        # 5' region
        if self.five_prime and self.five_prime.region_type != 'none':
            if self.five_prime.region_type == 'umi':
                sequence_spec.append({
                    'region_id': 'umi_5p',
                    'region_type': 'umi',
                    'sequence': 'N' * self.five_prime.length,
                    'min_len': self.five_prime.length,
                    'max_len': self.five_prime.length,
                })
            elif self.five_prime.region_type == 'adapter':
                sequence_spec.append({
                    'region_id': 'adapter_5p',
                    'region_type': 'adapter',
                    'sequence': self.five_prime.consensus_sequence,
                    'min_len': self.five_prime.length,
                    'max_len': self.five_prime.length,
                })

        # RPF region
        if self.rpf_length_distribution:
            lengths = list(self.rpf_length_distribution.keys())
            sequence_spec.append({
                'region_id': 'rpf',
                'region_type': 'cdna',
                'sequence': 'N' * 28,  # Placeholder
                'min_len': min(lengths),
                'max_len': max(lengths),
            })

        # 3' region
        if self.three_prime and self.three_prime.region_type != 'none':
            if self.three_prime.region_type == 'adapter':
                sequence_spec.append({
                    'region_id': 'adapter_3p',
                    'region_type': 'adapter',
                    'sequence': self.three_prime.consensus_sequence,
                    'min_len': len(self.three_prime.consensus_sequence or ''),
                    'max_len': len(self.three_prime.consensus_sequence or ''),
                })
            elif self.three_prime.region_type == 'umi':
                sequence_spec.append({
                    'region_id': 'umi_3p',
                    'region_type': 'umi',
                    'sequence': 'N' * self.three_prime.length,
                    'min_len': self.three_prime.length,
                    'max_len': self.three_prime.length,
                })

        return {
            'seqspec_version': '0.3.0',
            'assay_id': 'learned_structure',
            'name': 'Learned Read Structure',
            'sequence_spec': sequence_spec,
            'quality_markers': {
                'confidence': self.overall_confidence,
                'warnings': self.validation_warnings,
            }
        }


@dataclass
class TrimmerConfig:
    """Actionable configuration for read trimming and RPF extraction.

    This is passed to the extraction phase to process all reads.
    The key insight: learn once, apply to all reads efficiently.
    """
    trim_5p_fixed: Optional[int] = None  # Fixed bases to trim from 5'
    trim_5p_is_umi: bool = False  # If True, extract UMI to header
    trim_3p_adapter: Optional[str] = None  # Adapter sequence to search
    trim_3p_adapter_min_overlap: int = 8  # Minimum bases to match
    trim_3p_adapter_max_mismatches: int = 1  # Allow some errors
    min_rpf_length: int = 20
    max_rpf_length: int = 40

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'trim_5p_fixed': self.trim_5p_fixed,
            'trim_5p_is_umi': self.trim_5p_is_umi,
            'trim_3p_adapter': self.trim_3p_adapter,
            'trim_3p_adapter_min_overlap': self.trim_3p_adapter_min_overlap,
            'trim_3p_adapter_max_mismatches': self.trim_3p_adapter_max_mismatches,
            'min_rpf_length': self.min_rpf_length,
            'max_rpf_length': self.max_rpf_length,
        }


class StructureLearningError(Exception):
    """Raised when structure learning fails with low confidence."""
    pass


# =============================================================================
# Original Types (Kept for Compatibility)
# =============================================================================

@dataclass
class ReadArchitecture:
    """Represents a known read architecture from ribosome profiling protocols."""
    
    protocol_name: str
    lab_source: str
    umi_positions: List[Tuple[int, int]]
    barcode_positions: List[Tuple[int, int]]
    adapter_sequences: List[str]
    rpf_start: int
    rpf_end: int
    expected_rpf_length: Tuple[int, int]
    quality_markers: Dict[str, Any]
    confidence: float = 1.0


@dataclass
class SegmentInfo:
    """Information about a detected segment in reads."""
    
    segment_type: str        # "umi", "barcode", "adapter", "rpf", "unknown"
    start_pos: int          # start position in read
    end_pos: int            # end position in read
    confidence: float       # confidence score for this classification
    consensus: Optional[str] = None


@dataclass
class RPFExtractionResult:
    """Results from RPF extraction process."""
    
    input_reads: int
    processed_reads: int
    extracted_rpfs: int
    failed_extractions: int
    architecture_match: Optional[str]
    extraction_method: str  # "pattern_match" or "probabilistic_segmentation"
    segments: Dict[str, List[SegmentInfo]]
    quality_metrics: Dict[str, float]
    seqspec_data: Optional[Dict[str, Any]] = None
    trim_recommendations: Optional[Dict[str, Any]] = None
    
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
