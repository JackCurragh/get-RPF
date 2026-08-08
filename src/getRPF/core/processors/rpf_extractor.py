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
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter

from ...utils.file_utils import get_file_opener
from ..seqspec_generator import SeqSpecGenerator
from ..seqspec_loader import SeqSpecArchitectureLoader

# Import types and components
from .types import (
    ExtractionEmptyError,
    ReadArchitecture,
    SegmentInfo,
    RPFExtractionResult,
)
from .signals import process_reads
from .matcher import MatchResult, match_architecture
from .segmenter import segment_reads
from .reporting import generate_cli_report, write_html_report

logger = logging.getLogger(__name__)

MIN_ADAPTER_PREFIX_OVERLAP = 10
MIN_ADAPTER_EVIDENCE_FRACTION = 0.10
LOW_ADAPTER_EVIDENCE_FRACTION = 0.05
MAX_ADAPTER_EVIDENCE_CANDIDATES = 5
MIN_RPF_LENGTH = 20
MAX_RPF_LENGTH = 40
FOOTPRINT_CORE_MIN = 26
FOOTPRINT_CORE_MAX = 34
DISOME_CANDIDATE_MIN = 50
DISOME_CANDIDATE_MAX = 80
MIN_DISOME_CANDIDATE_FRACTION = 0.20
MIN_PRETRIMMED_RPF_FRACTION = 0.80
MIN_PRETRIMMED_CORE_FRACTION = 0.50
MIN_RETAINED_FRACTION_WARN = 0.05
MIN_EXTRACTED_CORE_FRACTION_WARN = 0.10
GENERIC_ARCHITECTURE = "comprehensive_adapter_check"
MIN_STRONG_ADAPTER_EVIDENCE_FRACTION = 0.50
MIN_ADAPTER_EVIDENCE_MARGIN = 0.10


def resolve_architecture_choice(
    match_result: Optional[MatchResult],
    adapter_evidence: List[Dict[str, Any]],
) -> Tuple[Optional[ReadArchitecture], str, List[str]]:
    """Choose an architecture without allowing generic matches to mask evidence.

    The comprehensive catalogue is useful as a fallback, but it is not a
    protocol identification result. A strong, specific adapter signal must
    therefore override a generic strict match. For named architectures, a
    competing specific adapter only wins when it is both strongly supported
    and clearly stronger than the selected architecture's own evidence.
    """
    best_evidence = adapter_evidence[0] if adapter_evidence else None
    best_specific = next(
        (
            item
            for item in adapter_evidence
            if item["protocol_name"] != GENERIC_ARCHITECTURE
        ),
        None,
    )

    if match_result and match_result.is_match:
        selected = match_result.architecture
        selected_evidence = next(
            (
                item
                for item in adapter_evidence
                if item["protocol_name"] == selected.protocol_name
            ),
            None,
        )
        selected_fraction = (
            selected_evidence["hit_fraction"] if selected_evidence else 0.0
        )
        specific_fraction = best_specific["hit_fraction"] if best_specific else 0.0
        generic_match = selected.protocol_name == GENERIC_ARCHITECTURE
        specific_wins = (
            best_specific is not None
            and specific_fraction >= MIN_STRONG_ADAPTER_EVIDENCE_FRACTION
            and (
                generic_match
                or specific_fraction >= selected_fraction + MIN_ADAPTER_EVIDENCE_MARGIN
            )
        )
        if specific_wins:
            return (
                best_specific["architecture"],
                "adapter_evidence_match",
                [
                    "Replaced generic/unsupported strict architecture with "
                    f"specific adapter evidence: {best_specific['protocol_name']} "
                    f"({specific_fraction:.1%} of sampled reads)."
                ],
            )
        return selected, "strict_pattern_match", list(match_result.reasons)

    if best_evidence and best_evidence["hit_fraction"] >= MIN_ADAPTER_EVIDENCE_FRACTION:
        return (
            best_evidence["architecture"],
            "adapter_evidence_match",
            [
                "Selected architecture by adapter evidence: "
                f"{best_evidence['protocol_name']} "
                f"({best_evidence['hit_fraction']:.1%} of sampled reads)."
            ],
        )

    return None, "unknown", []





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
        """Initialize architectures from the internal package directory."""
        try:
            import importlib.resources
            # For Python 3.9+, use files()
            if hasattr(importlib.resources, 'files'):
                arch_dir = importlib.resources.files('getRPF').joinpath('architectures')
            else:
                # Fallback maybe not needed if >=3.10 is required
                raise ImportError("Requires Python 3.9+ for importlib.resources.files")
                
            self.load_from_seqspec_directory(arch_dir)
            
            # Ensure comprehensive check is first (if loaded)
            comprehensive = next((a for a in self.architectures if a.protocol_name == "comprehensive_adapter_check"), None)
            if comprehensive:
                self.architectures.remove(comprehensive)
                self.architectures.insert(0, comprehensive)
                
            logger.info(f"Initialized {len(self.architectures)} built-in architectures from package resources")
            
        except Exception as e:
            logger.error(f"Failed to initialize built-in architectures: {e}")
            self.architectures = []

        if self.seqspec_dir and self.seqspec_dir.exists():
            self.load_from_seqspec_directory(self.seqspec_dir)
            
    def load_from_seqspec_directory(self, seqspec_dir: Any) -> int:
        seqspec_architectures = self.seqspec_loader.load_from_directory(seqspec_dir)
        self.architectures.extend(seqspec_architectures)
        return len(seqspec_architectures)


class RPFExtractor:
    """Main RPF extraction processor."""
    
    def __init__(self, architecture_db_path: Optional[Path] = None, seqspec_dir: Optional[Path] = None):
        self.architecture_db = ArchitectureDatabase(db_path=architecture_db_path, seqspec_dir=seqspec_dir)
        self.seqspec_generator = SeqSpecGenerator()
        
        # New components

    def infer_architecture(
        self,
        sample_reads: List[str],
        signals=None,
    ) -> Tuple[Optional[ReadArchitecture], str, List[str], List[Dict[str, Any]]]:
        """Infer architecture once for orchestration and extraction.

        ``signals`` may be supplied by the Stage-0 sketch so callers do not
        recompute the pooled signal solely to decide whether the expensive
        boundary planner is needed.
        """
        signals = signals or process_reads(sample_reads, compute_dinucleotide=False)
        adapter_evidence = self._score_adapter_evidence(sample_reads)
        match_result = match_architecture(signals, self.architecture_db.architectures)
        architecture, method, trace_log = resolve_architecture_choice(
            match_result, adapter_evidence
        )
        return architecture, method, trace_log, adapter_evidence
    
    def extract_rpfs(
        self,
        input_file: Path,
        output_file: Path,
        format: str = "fastq",
        max_reads: Optional[int] = None,
        generate_seqspec: bool = False,
        collapse_output: bool = True,
        collapsed_only: bool = False,
        override_trims: Optional[Dict[int, Dict[str, int]]] = None,
        sample_reads: Optional[List[str]] = None,
        sample_headers: Optional[List[str]] = None,
    ) -> RPFExtractionResult:
        """Extract RPFs from input file.

        Args:
            override_trims: Optional per-raw-read-length trim rules, e.g.
                {29: {"trim_5p": 0, "trim_3p": 1}}. When a read's pre-trim
                length has an entry here, that exact trim is applied instead
                of the single global architecture-derived trim. This is the
                apply path for boundary.py's per-length-class decisions; see
                docs/release_qc_and_terminal_trimming_plan.md M3.
        """
        logger.info(f"Starting RPF extraction from {input_file}")
        
        # 1. Reuse the bounded analysis sample when the caller already built
        # one (the pipeline does this after sketching). Direct callers retain
        # the original file-loading behavior.
        if sample_reads is None:
            sample_reads, loaded_headers = self._load_sample_reads(
                input_file, format, sample_size=50000
            )
            sample_headers = loaded_headers
        else:
            sample_reads = sample_reads[:50000]
            sample_headers = sample_headers or []
        
        # 2. Generate Signals (Entropy, Composition)
        logger.info("Generating signal metrics...")
        signals = process_reads(sample_reads, compute_dinucleotide=False)
        sample_length_counts = Counter(len(read) for read in sample_reads)
        raw_length_profile = self._length_profile_from_length_counts(
            sample_length_counts
        )
        (
            final_architecture,
            method,
            trace_log,
            adapter_evidence,
        ) = self.infer_architecture(sample_reads, signals)
        best_adapter_evidence = adapter_evidence[0] if adapter_evidence else None
        
        # 3. Match against known architectures
        logger.info("Matching architectures...")
        extracted_segments = []
        if final_architecture is not None:
            logger.info(
                "Selected architecture: "
                f"{final_architecture.protocol_name} ({method})"
            )
            extracted_segments = self._arch_to_segments(final_architecture)

        if final_architecture is None and self._looks_pretrimmed(raw_length_profile):
            logger.info(
                "No dominant adapter evidence and reads are already footprint-like. "
                "Using conservative length-filter extraction."
            )
            method = "pretrimmed_length_filter"
            final_architecture = ReadArchitecture(
                protocol_name="pretrimmed_no_dominant_adapter",
                lab_source="Length-filtered raw reads",
                umi_positions=[],
                barcode_positions=[],
                adapter_sequences=[],
                rpf_start=0,
                rpf_end=-1,
                expected_rpf_length=(MIN_RPF_LENGTH, MAX_RPF_LENGTH),
                quality_markers={
                    "pretrimmed_rpf_fraction": raw_length_profile["frac_20_40"],
                    "pretrimmed_core_fraction": raw_length_profile["frac_26_34"],
                },
            )
            extracted_segments = [SegmentInfo("rpf", 0, -1, 0.8)]
            trace_log.append(
                "Reads are mostly within the RPF length window and no adapter "
                "candidate is strongly supported; extracted by length filter only."
            )

        if final_architecture is None:
            # Fallback to Probabilistic Segmentation
            logger.info("No strict architecture match found. Attempting probabilistic segmentation...")
            method = "probabilistic_hmm"
            
            # IMPROVEMENT: Multi-Scale Binning Strategy
            # Variable read lengths can smear the signal. We bin reads by length and detect on each bin independently.
            reads_by_length: Dict[int, List[str]] = {}
            for read in sample_reads:
                reads_by_length.setdefault(len(read), []).append(read)
            length_counts = sample_length_counts
            
            # Select bins with sufficient coverage (>1000 reads)
            valid_bins = [l for l, c in length_counts.items() if c > 1000]
            valid_bins.sort()
            
            detected_segments_per_bin = []
            
            if len(valid_bins) > 1:
                logger.info(f"Multi-scale detection: analyzing {len(valid_bins)} length bins: {valid_bins}")
                for length in valid_bins:
                    bin_reads = reads_by_length[length]
                    bin_signals = process_reads(bin_reads, compute_dinucleotide=False)
                    bin_segments = segment_reads(bin_signals)
                    
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
                 extracted_segments = segment_reads(signals)

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
                
                if (
                    best_adapter_evidence
                    and best_adapter_evidence["hit_fraction"] >= MIN_ADAPTER_EVIDENCE_FRACTION
                ):
                    logger.info(
                        "Fallback Scan: Found adapter "
                        f"{best_adapter_evidence['adapter']} in "
                        f"{best_adapter_evidence['hit_count']} reads. Overriding HMM."
                    )
                    adapter_seq = best_adapter_evidence["adapter"]
                    
                    # Construct override architecture
                    final_architecture = ReadArchitecture(
                        protocol_name=f"fallback_detected_{best_adapter_evidence['protocol_name']}",
                        lab_source="Brute-force Scan Fallback",
                        umi_positions=[],
                        barcode_positions=[],
                        adapter_sequences=[adapter_seq],
                        trim_adapter_sequences=[adapter_seq],
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
                    trace_log.append(
                        "Override: Found dominant adapter "
                        f"{best_adapter_evidence['protocol_name']}."
                    )
        
        # 4. Generate Reports
        if final_architecture and extracted_segments:
            # CLI Report
            cli_report = generate_cli_report(
                final_architecture, extracted_segments, signals
            )
            logger.info("\n" + cli_report)
            
            # HTML Report
            html_report_path = output_file.with_suffix(".report.html")
            write_html_report(
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
                 trim_adapters = getattr(
                     final_architecture,
                     "trim_adapter_sequences",
                     None,
                 )
                 adapters_to_trim = (
                     trim_adapters
                     if trim_adapters is not None
                     else final_architecture.adapter_sequences
                 )
                 if final_architecture.rpf_end == -1 and not adapters_to_trim:
                     pass

            extraction_stats = self._extract_rpfs_from_reads(
                input_file, output_file, extracted_segments, format, max_reads,
                adapters=adapters_to_trim,
                post_rpf_trim_bases=final_architecture.post_rpf_trim_bases,
                collapse_output=collapse_output,
                collapsed_only=collapsed_only,
                override_trims=override_trims,
            )
            
            adapter_report = self._adapter_reporting(
                final_architecture,
                method,
                adapter_evidence,
            )
            adapter_quality = self._adapter_quality_metrics(best_adapter_evidence)
            confidence = 0.95
            if adapter_report["adapter_source"] == "no_dominant_adapter":
                confidence = 0.5
                adapter_quality["architecture_confidence_note"] = (
                    "No dominant adapter evidence in raw reads; sample may already "
                    "be footprint-like, pretrimmed, or mixed-length."
                )
            if method == "pretrimmed_length_filter":
                confidence = 0.75
                adapter_quality["architecture_confidence_note"] = (
                    "No seqspec was inferred. Reads were treated as already "
                    "footprint-like and only the broad RPF length window was applied."
                )
            extraction_warnings = self._extraction_warnings(extraction_stats)
            extraction_class = self._extraction_class(
                method,
                adapter_report,
                extraction_warnings,
            )

            return RPFExtractionResult(
                input_reads=extraction_stats["input_reads"],
                processed_reads=extraction_stats["input_reads"],
                extracted_rpfs=extraction_stats["extracted_rpfs"],
                failed_extractions=0,
                architecture_match=final_architecture.protocol_name,
                extraction_method=method,
                segments={"detected": extracted_segments},
                quality_metrics={
                    "confidence": confidence,
                    **adapter_quality,
                    "raw_length_profile": raw_length_profile,
                    "extracted_length_profile": extraction_stats["extracted_length_profile"],
                    "retained_fraction": extraction_stats["retained_fraction"],
                    "unique_extracted_sequences": extraction_stats["unique_extracted_sequences"],
                    "rpf_length_gated_reads": extraction_stats["rpf_length_gated_reads"],
                    "rpf_length_gated_unique_sequences": extraction_stats["rpf_length_gated_unique_sequences"],
                    "rpf_length_gated_fraction": extraction_stats["rpf_length_gated_fraction"],
                    "rpf_length_profile": extraction_stats["rpf_length_profile"],
                    "extraction_class": extraction_class,
                    "extraction_warnings": extraction_warnings,
                },
                seqspec_data=seqspec_data,
                trim_recommendations=self._derive_trim_recommendations(extracted_segments),
                adapter_source=adapter_report["adapter_source"],
                adapter_conflict=adapter_report["adapter_conflict"],
                adapter_evidence_candidates=adapter_report["adapter_evidence_candidates"],
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

    def _score_adapter_evidence(self, reads: List[str]) -> List[Dict[str, Any]]:
        """Score observed adapter-prefix evidence for every loaded architecture."""
        if not reads:
            return []

        # The same adapter sequence is present in multiple architecture
        # records. Match each distinct sequence once, then reuse its count for
        # the architecture-specific report entries.
        normalized_reads = [read.upper() for read in reads]
        adapter_sequences = {
            adapter.upper()
            for arch in self.architecture_db.architectures
            for adapter in arch.adapter_sequences
            if adapter
        }
        hit_counts = {
            adapter: sum(
                1
                for read in normalized_reads
                if self._find_adapter_prefix(read, adapter) is not None
            )
            for adapter in adapter_sequences
        }

        scored = []
        seen = set()
        for arch in self.architecture_db.architectures:
            for adapter in arch.adapter_sequences:
                adapter = adapter.upper()
                if not adapter:
                    continue

                key = (arch.protocol_name, adapter)
                if key in seen:
                    continue
                seen.add(key)

                hit_count = hit_counts[adapter]
                if hit_count == 0:
                    continue

                scored.append(
                    {
                        "architecture": arch,
                        "protocol_name": arch.protocol_name,
                        "adapter": adapter,
                        "hit_count": hit_count,
                        "hit_fraction": hit_count / len(reads),
                    }
                )

        scored.sort(
            key=lambda item: (
                item["hit_count"],
                item["hit_fraction"],
                item["protocol_name"] != "comprehensive_adapter_check",
                len(item["adapter"]),
                item["protocol_name"],
            ),
            reverse=True,
        )
        return scored

    def _adapter_quality_metrics(
        self, best_adapter_evidence: Optional[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Return serializable quality metrics for adapter-evidence reporting."""
        if not best_adapter_evidence:
            return {
                "adapter_evidence_protocol": None,
                "adapter_evidence_sequence": None,
                "adapter_evidence_hit_count": 0,
                "adapter_evidence_hit_fraction": 0.0,
            }

        return {
            "adapter_evidence_protocol": best_adapter_evidence["protocol_name"],
            "adapter_evidence_sequence": best_adapter_evidence["adapter"],
            "adapter_evidence_hit_count": best_adapter_evidence["hit_count"],
            "adapter_evidence_hit_fraction": best_adapter_evidence["hit_fraction"],
        }

    def _looks_pretrimmed(self, length_profile: Dict[str, Any]) -> bool:
        """Return True when sampled reads are already mostly plausible RPFs."""
        return (
            length_profile["frac_20_40"] >= MIN_PRETRIMMED_RPF_FRACTION
            and length_profile["frac_26_34"] >= MIN_PRETRIMMED_CORE_FRACTION
        )

    def _length_profile(self, reads: List[str]) -> Dict[str, Any]:
        """Summarize sampled read lengths for reporting and extraction policy."""
        counts = Counter(len(read) for read in reads)
        return self._length_profile_from_length_counts(counts)

    def _length_profile_from_sequence_counts(self, sequence_counts: Counter) -> Dict[str, Any]:
        """Summarize sequence length distribution from collapsed sequence counts."""
        counts = Counter()
        for sequence, count in sequence_counts.items():
            counts[len(sequence)] += count
        return self._length_profile_from_length_counts(counts)

    def _length_profile_from_length_counts(self, counts: Counter) -> Dict[str, Any]:
        """Summarize a weighted length-count mapping."""
        total = sum(counts.values())

        def fraction(min_len: int, max_len: Optional[int] = None) -> float:
            if total == 0:
                return 0.0
            if max_len is None:
                n = sum(count for length, count in counts.items() if length >= min_len)
            else:
                n = sum(
                    count
                    for length, count in counts.items()
                    if min_len <= length <= max_len
                )
            return n / total

        mode_length = None
        mode_count = 0
        if counts:
            mode_length, mode_count = counts.most_common(1)[0]

        return {
            "reads_scanned": total,
            "mode_length": mode_length,
            "mode_count": mode_count,
            "frac_20_40": fraction(MIN_RPF_LENGTH, MAX_RPF_LENGTH),
            "frac_26_34": fraction(FOOTPRINT_CORE_MIN, FOOTPRINT_CORE_MAX),
            "frac_35_40": fraction(35, MAX_RPF_LENGTH),
            "frac_50_80": fraction(DISOME_CANDIDATE_MIN, DISOME_CANDIDATE_MAX),
            "frac_gt40": fraction(MAX_RPF_LENGTH + 1),
            "top_lengths": [
                {"length": length, "count": count}
                for length, count in counts.most_common(10)
            ],
        }

    def _adapter_reporting(
        self,
        final_architecture: ReadArchitecture,
        method: str,
        adapter_evidence: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Summarize adapter evidence and conflicts for JSON reports."""
        candidates = [
            self._serialize_adapter_candidate(item)
            for item in adapter_evidence[:MAX_ADAPTER_EVIDENCE_CANDIDATES]
        ]
        selected_evidence = self._selected_adapter_evidence(
            final_architecture, adapter_evidence
        )
        best_evidence = adapter_evidence[0] if adapter_evidence else None
        selected_fraction = (
            selected_evidence["hit_fraction"] if selected_evidence else 0.0
        )

        if method in {"adapter_evidence_match", "fallback_adapter_scan"}:
            adapter_source = "raw_adapter_evidence"
        elif method == "pretrimmed_length_filter":
            adapter_source = "no_dominant_adapter"
        elif (
            best_evidence is None
            or best_evidence["hit_fraction"] < MIN_ADAPTER_EVIDENCE_FRACTION
        ) and final_architecture.protocol_name == "de_novo_inferred":
            adapter_source = "no_dominant_adapter"
        elif method == "strict_pattern_match":
            adapter_source = "architecture_match"
        else:
            adapter_source = "de_novo_inferred"

        conflict = {
            "has_conflict": False,
            "type": None,
            "selected_protocol": final_architecture.protocol_name,
            "selected_hit_fraction": selected_fraction,
            "best_supported_protocol": best_evidence["protocol_name"] if best_evidence else None,
            "best_supported_hit_fraction": best_evidence["hit_fraction"] if best_evidence else 0.0,
            "message": None,
        }
        if (
            adapter_source == "architecture_match"
            and best_evidence
            and best_evidence["protocol_name"] != final_architecture.protocol_name
            and best_evidence["hit_fraction"] >= MIN_ADAPTER_EVIDENCE_FRACTION
            and selected_fraction < LOW_ADAPTER_EVIDENCE_FRACTION
        ):
            conflict.update(
                {
                    "has_conflict": True,
                    "type": "selected_architecture_unsupported",
                    "message": (
                        "Selected architecture has low/no raw adapter support, "
                        "while another adapter candidate is strongly supported."
                    ),
                }
            )

        return {
            "adapter_source": adapter_source,
            "adapter_conflict": conflict,
            "adapter_evidence_candidates": candidates,
        }

    def _extraction_warnings(self, extraction_stats: Dict[str, Any]) -> List[str]:
        """Return stable quality warnings for downstream gating."""
        warnings = []
        if extraction_stats["retained_fraction"] < MIN_RETAINED_FRACTION_WARN:
            warnings.append("low_retained_fraction")
        extracted_profile = extraction_stats.get("extracted_length_profile") or {}
        if extracted_profile.get("frac_50_80", 0.0) >= MIN_DISOME_CANDIDATE_FRACTION:
            warnings.append("disome_like_length_distribution")
        if extracted_profile.get("frac_26_34", 0.0) < MIN_EXTRACTED_CORE_FRACTION_WARN:
            warnings.append("low_26_34_fraction")
        return warnings

    def _extraction_class(
        self,
        method: str,
        adapter_report: Dict[str, Any],
        extraction_warnings: Optional[List[str]] = None,
    ) -> str:
        """Return a stable high-level class for downstream QC gating."""
        extraction_warnings = extraction_warnings or []
        if "low_retained_fraction" in extraction_warnings:
            return "adapter_supported_low_yield"
        if "disome_like_length_distribution" in extraction_warnings:
            return "needs_length_policy_review"
        if "low_26_34_fraction" in extraction_warnings:
            return "extracted_rpf_length_suspicious"
        conflict = adapter_report.get("adapter_conflict") or {}
        if conflict.get("has_conflict"):
            return "adapter_conflict_raw_adapter_wins"
        if method == "pretrimmed_length_filter":
            return "pretrimmed_rpf_length_filter"
        if method in {"adapter_evidence_match", "fallback_adapter_scan"}:
            return "raw_adapter_supported_trimmed"
        if adapter_report.get("adapter_source") == "no_dominant_adapter":
            return "no_dominant_adapter_low_confidence"
        return "architecture_inferred"

    def _selected_adapter_evidence(
        self,
        final_architecture: ReadArchitecture,
        adapter_evidence: List[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        for item in adapter_evidence:
            if item["protocol_name"] == final_architecture.protocol_name:
                return item
        return None

    def _adapter_evidence_for_protocol(
        self, protocol_name: str, adapter_evidence: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        for item in adapter_evidence:
            if item["protocol_name"] == protocol_name:
                return item
        return None

    def _serialize_adapter_candidate(self, item: Dict[str, Any]) -> Dict[str, Any]:
        arch = item["architecture"]
        return {
            "protocol_name": item["protocol_name"],
            "adapter": item["adapter"],
            "hit_count": item["hit_count"],
            "hit_fraction": item["hit_fraction"],
            "adapter_source": self._architecture_source(arch),
        }

    def _architecture_source(self, arch: ReadArchitecture) -> str:
        if arch.protocol_name.startswith("observed_"):
            return "observed_protocol_catalog"
        if arch.protocol_name == "comprehensive_adapter_check":
            return "generic_adapter_catalog"
        return "known_architecture"

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
        adapters: Optional[List[str]] = None,
        post_rpf_trim_bases: int = 0,
        collapse_output: bool = True,
        collapsed_only: bool = False,
        override_trims: Optional[Dict[int, Dict[str, int]]] = None,
    ) -> Dict[str, Any]:
        """Extract RPF sequences from reads using Two-Stage Collapsing for performance."""
        from .collapsed import TwoStageCollapser

        # Find RPF segment boundaries
        rpf_segments = [seg for seg in segments if seg.segment_type == "rpf"]
        if not rpf_segments:
            rpf_start, rpf_end = 0, -1
        else:
            rpf_segment = rpf_segments[0]
            rpf_start, rpf_end = rpf_segment.start_pos, rpf_segment.end_pos

        sorted_adapters = sorted(adapters, key=len, reverse=True) if adapters else None

        def trim_logic(sequence: str) -> Optional[str]:
            """Inner function to trim a single unique sequence."""
            rule = override_trims.get(len(sequence)) if override_trims else None
            trim_5p = rule.get("trim_5p", 0) if rule is not None else 0
            trim_3p = rule.get("trim_3p", 0) if rule is not None else 0

            if rpf_end == -1:
                # Architecture-derived boundaries take precedence over a
                # terminal artifact rule. If an adapter is present, trimming
                # it already removes the terminal bases covered by trim_3p;
                # applying both would cut into the RPF.
                rpf_seq = sequence[rpf_start + trim_5p:]
                adapter_found = False
                if sorted_adapters:
                    for adapter in sorted_adapters:
                        adapter_pos = self._find_adapter_prefix(rpf_seq, adapter)
                        if adapter_pos is not None:
                            rpf_seq = rpf_seq[:adapter_pos]
                            adapter_found = True
                            break
                if not adapter_found and trim_3p:
                    rpf_seq = rpf_seq[:-trim_3p]
                if adapter_found and post_rpf_trim_bases:
                    rpf_seq = rpf_seq[:-post_rpf_trim_bases]
                return rpf_seq
            else:
                return sequence[rpf_start + trim_5p:rpf_end]

        collapser = TwoStageCollapser(logger=logger)
        
        # Stage 1: Raw collapse
        raw_counts = collapser.collapse_raw(input_file, format=format, max_reads=max_reads)
        
        # Stage 2 & 3: Trim unique and merge. Do not apply the RPF length gate
        # here; callers/workflows should decide which trimmed reads are released.
        final_counts = collapser.apply_trimming(
            raw_counts,
            trim_logic,
            min_length=1,
            max_length=None,
        )
        rpf_length_counts = Counter(
            {
                seq: count
                for seq, count in final_counts.items()
                if MIN_RPF_LENGTH <= len(seq) <= MAX_RPF_LENGTH
            }
        )
        
        # Stage 4: Write outputs (skip empty)
        if collapse_output:
            collapsed_path = output_file.with_suffix('.collapsed.fa')
            collapser.write_collapsed_fasta(final_counts, collapsed_path)
            
        total_extracted = sum(final_counts.values())
        raw_total = sum(raw_counts.values())
        
        if not collapsed_only:
            logger.info(f"Writing expanded FASTQ output to {output_file}...")
            # We need to write the expanded FASTQ. Since we already have final_counts,
            # this is just regenerating the file from the unique set.
            with open(output_file, 'w') as fout:
                for idx, (seq, count) in enumerate(final_counts.items(), 1):
                    for i in range(count):
                        fout.write(f"@seq{idx}_c{i+1}_RPF\n{seq}\n+\n{'I' * len(seq)}\n")
        
        logger.info(f"Extracted {total_extracted} RPF sequences ({len(final_counts)} unique)")
        if total_extracted == 0:
            raise ExtractionEmptyError(
                "Extraction produced zero reads after architecture trimming."
            )
        return {
            "input_reads": raw_total,
            "extracted_rpfs": total_extracted,
            "retained_fraction": total_extracted / raw_total if raw_total else 0.0,
            "unique_extracted_sequences": len(final_counts),
            "extracted_length_profile": self._length_profile_from_sequence_counts(final_counts),
            "rpf_length_gated_reads": sum(rpf_length_counts.values()),
            "rpf_length_gated_unique_sequences": len(rpf_length_counts),
            "rpf_length_gated_fraction": (
                sum(rpf_length_counts.values()) / total_extracted if total_extracted else 0.0
            ),
            "rpf_length_profile": self._length_profile_from_sequence_counts(rpf_length_counts),
        }

    @staticmethod
    def _find_adapter_prefix(sequence: str, adapter: str) -> Optional[int]:
        """Return the 3' boundary position of a suffix adapter overlap.

        Adapter evidence may be inspected in other ways, but extraction must
        only trim when the read ends with an adapter prefix. Searching for an
        arbitrary internal match can truncate genuine RPF sequence.
        """
        sequence = sequence.upper()
        adapter = adapter.upper()
        max_k = min(len(adapter), len(sequence))
        min_overlap = min(MIN_ADAPTER_PREFIX_OVERLAP, len(adapter))
        for k in range(max_k, min_overlap - 1, -1):
            adapter_prefix = adapter[:k]
            if "N" in adapter_prefix:
                pos = RPFExtractor._find_iupac_suffix(sequence, adapter_prefix)
            else:
                pos = len(sequence) - k if sequence.endswith(adapter_prefix) else -1
            if pos >= 0:
                return pos
        return None

    @staticmethod
    def _find_iupac_suffix(sequence: str, adapter_prefix: str) -> int:
        """Find an adapter prefix at the read's 3' boundary."""
        prefix_len = len(adapter_prefix)
        pos = len(sequence) - prefix_len
        if pos < 0:
            return -1
        for idx, expected in enumerate(adapter_prefix):
            if expected != "N" and sequence[pos + idx] != expected:
                return -1
        return pos
