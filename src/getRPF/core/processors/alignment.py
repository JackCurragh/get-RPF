"""STAR alignment processor for ribosome profiling reads.

This module implements STAR alignment functionality for processing ribosome
protected fragments (RPFs) and detecting alignment features.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, List
import shutil
import json
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.warning("pysam not available - soft-clipping analysis will be disabled")

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
        self.trim_recommendations: Dict[str, Any] = {}
        self.per_length_analysis: Dict[str, Any] = {}
        self.pysam_available: bool = PYSAM_AVAILABLE

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
            "trim_recommendations": self.trim_recommendations,
            "per_length_analysis": self.per_length_analysis,
            "detected_features": self.features,
            "detected_features": self.features,
            "output_files": [str(f) for f in self.output_files if Path(f).exists()],
            "metadata": {
                "soft_clipping_analysis_available": self.pysam_available
            }
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
                # Add trim recommendations to CSV
                if self.trim_recommendations:
                    writer.writerow([])
                    writer.writerow(["Trim Recommendations", ""])
                    for key, value in self.trim_recommendations.items():
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

        # Handle compressed input
        if str(input_fastq).endswith(".gz"):
            star_cmd.extend(["--readFilesCommand", "zcat"])
        
        logger.info(f"Running STAR alignment: {' '.join(star_cmd)}")
        
        try:
            result = subprocess.run(
                star_cmd,
                capture_output=True,
                text=True,
                errors="replace",  # Handle potential binary garbage in logs
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
            # Add soft-clipping analysis
            stats.update(self._analyze_soft_clipping(bam_file))
        
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
    
    def _analyze_soft_clipping(self, bam_file: Path, threshold: float = 0.85) -> Dict[str, Any]:
        """Analyze soft-clipping patterns in aligned reads with per-length analysis.

        Args:
            bam_file: Path to sorted BAM file
            threshold: Percentage of reads to clean (default 0.85 = 85%)

        Returns:
            Dictionary with soft-clipping statistics and trim recommendations
        """
        soft_clip_stats = {
            "total_reads_analyzed": 0,
            "reads_with_soft_clips": 0,
            "mean_5prime_clips": 0.0,
            "mean_3prime_clips": 0.0,
            "max_5prime_clips": 0,
            "max_3prime_clips": 0,
        }

        if not PYSAM_AVAILABLE:
            logger.warning("pysam not available - skipping soft-clipping analysis")
            return soft_clip_stats

        try:
            with pysam.AlignmentFile(str(bam_file), "rb") as bamfile:
                clip_5prime = []
                clip_3prime = []
                reads_with_clips = 0

                # Per-read-length tracking
                per_length_data = {}  # {length: {'5p': [clips], '3p': [clips]}}

                for read in bamfile:
                    if read.is_unmapped:
                        continue

                    soft_clip_stats["total_reads_analyzed"] += 1

                    # Get read length
                    read_length = read.query_length

                    # Parse CIGAR for soft clips
                    cigar = read.cigartuples
                    if cigar:
                        # 5' soft clip (start of read)
                        clip_5 = cigar[0][1] if cigar[0][0] == 4 else 0  # 4 = soft clip
                        # 3' soft clip (end of read)
                        clip_3 = cigar[-1][1] if cigar[-1][0] == 4 else 0

                        clip_5prime.append(clip_5)
                        clip_3prime.append(clip_3)

                        # Track per length
                        if read_length not in per_length_data:
                            per_length_data[read_length] = {'5p': [], '3p': []}
                        per_length_data[read_length]['5p'].append(clip_5)
                        per_length_data[read_length]['3p'].append(clip_3)

                        if clip_5 > 0 or clip_3 > 0:
                            reads_with_clips += 1

                if clip_5prime:
                    soft_clip_stats["reads_with_soft_clips"] = reads_with_clips
                    soft_clip_stats["mean_5prime_clips"] = sum(clip_5prime) / len(clip_5prime)
                    soft_clip_stats["mean_3prime_clips"] = sum(clip_3prime) / len(clip_3prime)
                    soft_clip_stats["max_5prime_clips"] = max(clip_5prime)
                    soft_clip_stats["max_3prime_clips"] = max(clip_3prime)

                    # Add trim recommendations
                    trim_recommendations = self._calculate_trim_recommendations(
                        per_length_data, threshold
                    )
                    soft_clip_stats.update(trim_recommendations)

        except Exception as e:
            logger.warning(f"Error analyzing soft clipping: {e}")

        return soft_clip_stats

    def _calculate_trim_recommendations(
        self,
        per_length_data: Dict[int, Dict[str, List[int]]],
        threshold: float = 0.85
    ) -> Dict[str, Any]:
        """Calculate trim recommendations from per-length soft-clipping data.

        Args:
            per_length_data: Dictionary mapping read length to clipping data
            threshold: Percentage of reads that should be cleaned

        Returns:
            Dictionary with trim recommendations and per-length analysis
        """
        recommendations = {
            "trim_recommendations": {
                "global_pattern_detected": False,
                "recommended_5prime_trim": 0,
                "recommended_3prime_trim": 0,
                "threshold_used": threshold,
                "consensus_level": 0.0,  # How consistent is the pattern across lengths
            },
            "per_length_analysis": {},
        }

        if not per_length_data:
            return recommendations

        # Analyze each read length
        length_recommendations = {}

        for length in sorted(per_length_data.keys()):
            clips_5p = per_length_data[length]['5p']
            clips_3p = per_length_data[length]['3p']

            if not clips_5p:
                continue

            n_reads = len(clips_5p)

            # Calculate percentile-based trim (what trim cleans threshold% of reads)
            trim_5p = self._calculate_percentile_trim(clips_5p, threshold)
            trim_3p = self._calculate_percentile_trim(clips_3p, threshold)

            # Calculate what percentage would be cleaned by this trim
            pct_cleaned_5p = sum(1 for c in clips_5p if c <= trim_5p) / n_reads
            pct_cleaned_3p = sum(1 for c in clips_3p if c <= trim_3p) / n_reads

            length_recommendations[length] = {
                "5prime_trim": trim_5p,
                "3prime_trim": trim_3p,
                "n_reads": n_reads,
                "mean_5prime_clips": sum(clips_5p) / n_reads,
                "mean_3prime_clips": sum(clips_3p) / n_reads,
                "pct_cleaned_by_trim": min(pct_cleaned_5p, pct_cleaned_3p),
            }

            recommendations["per_length_analysis"][str(length)] = length_recommendations[length]

        # Detect global consensus pattern
        if length_recommendations:
            global_pattern = self._detect_global_consensus(length_recommendations)
            recommendations["trim_recommendations"].update(global_pattern)

        return recommendations

    def _calculate_percentile_trim(self, clips: List[int], threshold: float) -> int:
        """Calculate minimum trim needed to clean threshold% of reads.

        Args:
            clips: List of clipping amounts
            threshold: Target percentage of reads to clean

        Returns:
            Minimum bases to trim
        """
        if not clips:
            return 0

        # Sort clips to find the value at threshold percentile
        sorted_clips = sorted(clips)
        target_idx = int(len(sorted_clips) * threshold)

        # Return the clip value at this percentile
        # This means: trim this many bases and you'll clean threshold% of reads
        return sorted_clips[min(target_idx, len(sorted_clips) - 1)]

    def _detect_global_consensus(
        self,
        length_recommendations: Dict[int, Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Detect if there's a global consensus trim pattern across read lengths.

        Args:
            length_recommendations: Per-length trim recommendations

        Returns:
            Dictionary with global pattern detection results
        """
        if not length_recommendations:
            return {
                "global_pattern_detected": False,
                "recommended_5prime_trim": 0,
                "recommended_3prime_trim": 0,
                "consensus_level": 0.0,
            }

        # Extract trim values weighted by read count
        trim_5p_values = []
        trim_3p_values = []
        total_reads = 0

        for data in length_recommendations.values():
            n_reads = data["n_reads"]
            trim_5p_values.extend([data["5prime_trim"]] * n_reads)
            trim_3p_values.extend([data["3prime_trim"]] * n_reads)
            total_reads += n_reads

        # Calculate mode (most common value) weighted by reads
        from collections import Counter
        counter_5p = Counter(trim_5p_values)
        counter_3p = Counter(trim_3p_values)

        mode_5p = counter_5p.most_common(1)[0][0] if counter_5p else 0
        mode_3p = counter_3p.most_common(1)[0][0] if counter_3p else 0

        # Calculate consensus level (what % of reads agree with mode)
        consensus_5p = counter_5p[mode_5p] / total_reads if total_reads > 0 else 0
        consensus_3p = counter_3p[mode_3p] / total_reads if total_reads > 0 else 0
        consensus_level = min(consensus_5p, consensus_3p)

        # Consider pattern "global" if >75% of reads agree
        global_pattern_detected = consensus_level >= 0.75

        # Additional check: variance across read lengths should be low
        unique_5p = len(counter_5p)
        unique_3p = len(counter_3p)

        # If we have many read lengths but only 1-2 unique trim values, that's strong consensus
        n_lengths = len(length_recommendations)
        if n_lengths >= 3:
            variance_check = (unique_5p <= 2 and unique_3p <= 2)
            global_pattern_detected = global_pattern_detected and variance_check

        result = {
            "global_pattern_detected": global_pattern_detected,
            "recommended_5prime_trim": mode_5p,
            "recommended_3prime_trim": mode_3p,
            "consensus_level": round(consensus_level, 3),
        }

        if global_pattern_detected:
            logger.info(
                f"Global trim pattern detected: {mode_5p}bp 5', {mode_3p}bp 3' "
                f"(consensus: {consensus_level:.1%})"
            )
        else:
            logger.warning(
                f"No strong global pattern detected (consensus: {consensus_level:.1%}). "
                f"Read lengths may require length-specific trimming."
            )

        return result

    def align_reads(
        self, 
        input_file: Path, 
        format: str, 
        count_pattern: Optional[str] = None,
        save_bam_path: Optional[Path] = None
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

            # Handle BAM file
            if "bam_file" in star_stats:
                bam_source = star_stats["bam_file"]
                if save_bam_path:
                    # Move BAM to requested location
                    save_bam_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(bam_source, save_bam_path)
                    results.output_files.append(save_bam_path)
                    logger.info(f"Saved BAM file to {save_bam_path}")
                else:
                    # Keep track of temp BAM (will be invalid after cleanup)
                    # We might want to warn about this or just not include it in final output
                    # The updated write_report handles this by checking existence
                    results.output_files.append(bam_source)

            # Add trim recommendations from soft-clipping analysis
            if "trim_recommendations" in star_stats:
                results.trim_recommendations = star_stats["trim_recommendations"]
            if "per_length_analysis" in star_stats:
                results.per_length_analysis = star_stats["per_length_analysis"]

            # TODO: Add feature detection on aligned reads
            results.features = {
                "feature_detection": "not_implemented",
                "note": "Feature detection will be added in next iteration"
            }

            return results
            
        finally:
            # Clean up temporary files
            cleanup_temp_files(self.temp_files)