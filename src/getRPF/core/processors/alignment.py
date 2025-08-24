"""STAR alignment processor for ribosome profiling reads.

This module implements STAR alignment functionality for processing ribosome
protected fragments (RPFs) and detecting alignment features.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, List
import json

from ...utils.file_utils import create_temp_file, cleanup_temp_files
from ..processors.collapsed import CollapsedFASTAProcessor

logger = logging.getLogger(__name__)


class AlignmentResults:
    """Container for STAR alignment results and statistics."""
    
    def __init__(self):
        self.input_reads: int = 0
        self.aligned_reads: int = 0
        self.uniquely_aligned: int = 0
        self.multimapping_reads: int = 0
        self.unmapped_reads: int = 0
        self.alignment_rate: float = 0.0
        self.features: Dict[str, Any] = {}
        self.output_files: List[Path] = []
        
    def write_report(self, output_path: Path, format: str = "json") -> None:
        """Write alignment results to file in specified format."""
        report_data = {
            "alignment_statistics": {
                "input_reads": self.input_reads,
                "aligned_reads": self.aligned_reads,
                "uniquely_aligned": self.uniquely_aligned,
                "multimapping_reads": self.multimapping_reads,
                "unmapped_reads": self.unmapped_reads,
                "alignment_rate": self.alignment_rate,
            },
            "detected_features": self.features,
            "output_files": [str(f) for f in self.output_files],
        }
        
        if format.lower() == "json":
            with open(output_path, 'w') as f:
                json.dump(report_data, f, indent=2)
        elif format.lower() == "csv":
            import csv
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["metric", "value"])
                for key, value in report_data["alignment_statistics"].items():
                    writer.writerow([key, value])
        else:
            raise ValueError(f"Unsupported output format: {format}")


class STARAligner:
    """STAR alignment processor for ribosome profiling reads."""
    
    def __init__(
        self,
        star_index: Path,
        threads: int = 1,
        max_reads: Optional[int] = None,
        star_executable: str = "STAR",
    ):
        """Initialize STAR aligner.
        
        Args:
            star_index: Path to STAR genome index directory
            threads: Number of threads for alignment
            max_reads: Maximum reads to process (for testing/development)
            star_executable: Path to STAR executable
        """
        self.star_index = star_index
        self.threads = threads
        self.max_reads = max_reads
        self.star_executable = star_executable
        self.temp_files: List[Path] = []
        
        # Validate STAR index
        self._validate_star_index()
        
    def _validate_star_index(self) -> None:
        """Validate STAR index directory exists and has required files."""
        if not self.star_index.exists():
            raise FileNotFoundError(f"STAR index directory not found: {self.star_index}")
            
        required_files = [
            "chrName.txt", "chrLength.txt", "chrStart.txt",
            "Genome", "genomeParameters.txt", "SA", "SAindex"
        ]
        
        for required_file in required_files:
            if not (self.star_index / required_file).exists():
                logger.warning(f"Missing STAR index file: {required_file}")
                
    def _prepare_input_file(
        self, 
        input_file: Path, 
        format: str, 
        count_pattern: Optional[str] = None
    ) -> Path:
        """Prepare input file for STAR alignment.
        
        Converts various input formats to FASTQ for STAR alignment.
        
        Args:
            input_file: Original input file
            format: Input file format
            count_pattern: Pattern for collapsed format
            
        Returns:
            Path to prepared FASTQ file
        """
        if format == "fastq":
            return input_file
            
        elif format == "collapsed":
            # Convert collapsed FASTA to FASTQ
            processor = CollapsedFASTAProcessor(count_pattern=count_pattern)
            temp_fastq = create_temp_file(suffix=".fastq")
            self.temp_files.append(temp_fastq)
            
            # Process collapsed reads and expand to FASTQ
            processor.expand_to_fastq(input_file, temp_fastq, max_reads=self.max_reads)
            return temp_fastq
            
        elif format == "fasta":
            # Convert FASTA to FASTQ with dummy quality scores
            temp_fastq = create_temp_file(suffix=".fastq")
            self.temp_files.append(temp_fastq)
            
            with open(input_file, 'r') as fin, open(temp_fastq, 'w') as fout:
                lines = fin.readlines()
                for i in range(0, len(lines), 2):
                    if i + 1 < len(lines):
                        header = lines[i].strip()
                        sequence = lines[i + 1].strip()
                        
                        # Convert to FASTQ format
                        fastq_header = header.replace('>', '@', 1)
                        quality = 'I' * len(sequence)  # High quality scores
                        
                        fout.write(f"{fastq_header}\n")
                        fout.write(f"{sequence}\n")
                        fout.write("+\n")
                        fout.write(f"{quality}\n")
                        
            return temp_fastq
            
        else:
            raise ValueError(f"Unsupported input format: {format}")
    
    def _run_star_alignment(self, input_fastq: Path) -> Dict[str, Any]:
        """Run STAR alignment and return statistics.
        
        Args:
            input_fastq: Path to input FASTQ file
            
        Returns:
            Dictionary containing alignment statistics
        """
        # Create temporary directory for STAR output
        temp_dir = Path(tempfile.mkdtemp(prefix="star_align_"))
        self.temp_files.append(temp_dir)
        
        # STAR command with ribosome profiling optimized parameters
        star_cmd = [
            self.star_executable,
            "--runThreadN", str(self.threads),
            "--genomeDir", str(self.star_index),
            "--readFilesIn", str(input_fastq),
            "--outFileNamePrefix", str(temp_dir / "Aligned_"),
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--alignIntronMax", "1",  # No introns for most prokaryotes
            "--alignEndsType", "EndToEnd",  # Require end-to-end alignment
            "--outFilterMultimapNmax", "1",  # Only unique alignments
            "--outFilterMismatchNmax", "1",  # Allow 1 mismatch
            "--seedSearchStartLmax", "20",  # Shorter seed for short reads
        ]
        
        logger.info(f"Running STAR alignment: {' '.join(star_cmd)}")
        
        try:
            result = subprocess.run(
                star_cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600,  # 1 hour timeout
            )
            
            logger.info("STAR alignment completed successfully")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"STAR alignment failed: {e}")
            logger.error(f"STAR stderr: {e.stderr}")
            raise RuntimeError(f"STAR alignment failed: {e.stderr}")
        
        except subprocess.TimeoutExpired:
            logger.error("STAR alignment timed out")
            raise RuntimeError("STAR alignment timed out after 1 hour")
        
        # Parse STAR log files for statistics
        log_file = temp_dir / "Aligned_Log.final.out"
        stats = self._parse_star_log(log_file)
        
        # Store output files
        bam_file = temp_dir / "Aligned_Aligned.sortedByCoord.out.bam"
        if bam_file.exists():
            stats["bam_file"] = bam_file
        
        return stats
    
    def _parse_star_log(self, log_file: Path) -> Dict[str, Any]:
        """Parse STAR log file for alignment statistics.
        
        Args:
            log_file: Path to STAR log file
            
        Returns:
            Dictionary of parsed statistics
        """
        stats = {
            "input_reads": 0,
            "uniquely_mapped": 0,
            "multimapping": 0,
            "unmapped": 0,
        }
        
        if not log_file.exists():
            logger.warning(f"STAR log file not found: {log_file}")
            return stats
        
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if "Number of input reads" in line:
                        stats["input_reads"] = int(line.split('\t')[-1])
                    elif "Uniquely mapped reads number" in line:
                        stats["uniquely_mapped"] = int(line.split('\t')[-1])
                    elif "Number of reads mapped to multiple loci" in line:
                        stats["multimapping"] = int(line.split('\t')[-1])
                    elif "Number of reads unmapped" in line:
                        stats["unmapped"] = int(line.split('\t')[-1])
                        
        except Exception as e:
            logger.warning(f"Error parsing STAR log: {e}")
            
        return stats
    
    def align_reads(
        self, 
        input_file: Path, 
        format: str, 
        count_pattern: Optional[str] = None
    ) -> AlignmentResults:
        """Perform STAR alignment on input reads.
        
        Args:
            input_file: Path to input sequence file
            format: Input file format
            count_pattern: Pattern for collapsed format
            
        Returns:
            AlignmentResults object with statistics and features
        """
        try:
            # Prepare input file for STAR
            fastq_file = self._prepare_input_file(input_file, format, count_pattern)
            
            # Run STAR alignment
            star_stats = self._run_star_alignment(fastq_file)
            
            # Create results object
            results = AlignmentResults()
            results.input_reads = star_stats.get("input_reads", 0)
            results.uniquely_aligned = star_stats.get("uniquely_mapped", 0)
            results.multimapping_reads = star_stats.get("multimapping", 0)
            results.unmapped_reads = star_stats.get("unmapped", 0)
            results.aligned_reads = results.uniquely_aligned + results.multimapping_reads
            
            if results.input_reads > 0:
                results.alignment_rate = results.aligned_reads / results.input_reads
            
            # Add BAM file to output if available
            if "bam_file" in star_stats:
                results.output_files.append(star_stats["bam_file"])
            
            # TODO: Add feature detection on aligned reads
            results.features = {
                "feature_detection": "not_implemented",
                "note": "Feature detection will be added in next iteration"
            }
            
            return results
            
        finally:
            # Clean up temporary files
            cleanup_temp_files(self.temp_files)