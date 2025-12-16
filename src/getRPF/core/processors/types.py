"""Shared types and dataclasses for RPF processing."""

from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Optional, Any
from pathlib import Path
import json

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
