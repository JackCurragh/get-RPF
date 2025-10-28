"""seqspec format generator for ribosome profiling read structures.

Generates machine-readable seqspec YAML files describing the detected
structure of ribosome profiling reads including UMI, barcode, adapter,
and RPF regions.
"""

import yaml
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class SeqSpecRegion:
    """Represents a seqspec Region object."""
    region_id: str
    region_type: str  # barcode, umi, adapter, cdna, etc.
    sequence_type: str  # fixed, random, onlist, joined
    name: str
    sequence: Optional[str] = None
    min_len: int = 0
    max_len: int = 0
    strand: str = "pos"
    regions: List['SeqSpecRegion'] = None
    
    def __post_init__(self):
        if self.regions is None:
            self.regions = []


@dataclass
class SeqSpecRead:
    """Represents a seqspec Read object."""
    read_id: str
    name: str
    modality: str
    primer_id: str
    strand: str
    min_len: int
    max_len: int
    files: List[Dict[str, str]]


@dataclass
class SeqSpecAssay:
    """Represents a complete seqspec Assay object."""
    seqspec_version: str = "0.3.0"
    assay_id: str = "ribosome_profiling_detected"
    name: str = "Ribosome Profiling - Detected Structure"
    doi: str = "auto-detected"
    publication_date: str = ""
    description: str = ""
    modalities: List[str] = None
    lib_struct: str = ""
    library_protocol: str = ""
    library_kit: str = ""
    sequence_protocol: str = ""
    sequence_kit: str = ""
    sequence_spec: List[SeqSpecRegion] = None
    library_spec: List[SeqSpecRegion] = None
    
    def __post_init__(self):
        if self.modalities is None:
            self.modalities = ["rna"]
        if self.sequence_spec is None:
            self.sequence_spec = []
        if self.library_spec is None:
            self.library_spec = []
        if not self.publication_date:
            self.publication_date = datetime.now().strftime("%Y-%m-%d")


class SeqSpecGenerator:
    """Generate seqspec YAML files from detected read structures."""
    
    def __init__(self):
        """Initialize seqspec generator."""
        pass
    
    def generate_from_architecture(
        self,
        architecture,
        sample_reads: List[str],
        detected_segments: List = None,
        output_file: Optional[Path] = None,
        sample_headers: List[str] = None
    ) -> Dict[str, Any]:
        """Generate seqspec from detected architecture or segments.
        
        Args:
            architecture: ReadArchitecture object or None for de novo
            sample_reads: Sample reads for analysis
            detected_segments: List of SegmentInfo objects from de novo detection
            output_file: Optional path to write YAML file
            sample_headers: Optional headers for structure parsing
            
        Returns:
            Dictionary representation of seqspec
        """
        # Check if headers contain structure information
        if sample_headers and self._has_structure_annotations(sample_headers):
            return self._generate_from_annotated_headers(
                sample_headers, sample_reads, output_file
            )
        elif architecture:
            return self._generate_from_known_architecture(
                architecture, sample_reads, output_file
            )
        else:
            return self._generate_from_detected_segments(
                detected_segments or [], sample_reads, output_file
            )
    
    def _generate_from_known_architecture(
        self,
        architecture,
        sample_reads: List[str],
        output_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Generate seqspec from known ReadArchitecture."""
        
        # Create assay metadata
        assay = SeqSpecAssay(
            assay_id=f"ribosome_profiling_{architecture.protocol_name}",
            name=f"Ribosome Profiling - {architecture.protocol_name}",
            description=f"Read structure detected as {architecture.protocol_name} protocol from {architecture.lab_source}",
            library_protocol=architecture.protocol_name,
            library_kit=architecture.lab_source
        )
        
        # Build sequence spec regions
        regions = []
        current_pos = 0
        
        # Add UMI regions
        for i, (start, end) in enumerate(architecture.umi_positions):
            umi_sequence = self._extract_consensus(sample_reads, start, end)
            regions.append(SeqSpecRegion(
                region_id=f"umi_{i+1}",
                region_type="umi",
                sequence_type="random",
                name=f"UMI {i+1}",
                sequence=umi_sequence if umi_sequence else "N" * (end - start),
                min_len=end - start,
                max_len=end - start
            ))
            current_pos = end
        
        # Add barcode regions
        for i, (start, end) in enumerate(architecture.barcode_positions):
            barcode_sequence = self._extract_consensus(sample_reads, start, end)
            regions.append(SeqSpecRegion(
                region_id=f"barcode_{i+1}",
                region_type="barcode",
                sequence_type="onlist" if barcode_sequence else "random",
                name=f"Sample Barcode {i+1}",
                sequence=barcode_sequence,
                min_len=end - start,
                max_len=end - start
            ))
            current_pos = end
        
        # Add adapters as separate regions
        for i, adapter_seq in enumerate(architecture.adapter_sequences):
            regions.append(SeqSpecRegion(
                region_id=f"adapter_{i+1}",
                region_type="illumina_p7" if "AGATCGG" in adapter_seq else "custom_primer",
                sequence_type="fixed",
                name=f"Adapter {i+1}",
                sequence=adapter_seq,
                min_len=len(adapter_seq),
                max_len=len(adapter_seq)
            ))
        
        # Add RPF region
        rpf_min, rpf_max = architecture.expected_rpf_length
        regions.append(SeqSpecRegion(
            region_id="rpf_sequence",
            region_type="cdna",
            sequence_type="joined",
            name="Ribosome Protected Fragment",
            min_len=rpf_min,
            max_len=rpf_max
        ))
        
        # Calculate total read structure
        read_lengths = [len(read) for read in sample_reads[:100]]
        min_read_len = min(read_lengths) if read_lengths else 20
        max_read_len = max(read_lengths) if read_lengths else 50
        
        # Create read specification
        read_spec = SeqSpecRead(
            read_id="R1",
            name="Read 1",
            modality="rna",
            primer_id="primer_5",
            strand="pos",
            min_len=min_read_len,
            max_len=max_read_len,
            files=[{"file_id": "R1.fastq.gz", "filename": "input.fastq"}]
        )
        
        assay.sequence_spec = regions
        
        # Generate the complete seqspec structure
        seqspec_dict = self._create_seqspec_dict(assay, [read_spec])
        
        if output_file:
            self._write_yaml(seqspec_dict, output_file)
            logger.info(f"seqspec file written to {output_file}")
        
        return seqspec_dict
    
    def _generate_from_detected_segments(
        self,
        detected_segments: List,
        sample_reads: List[str],
        output_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Generate seqspec from de novo detected segments."""
        
        assay = SeqSpecAssay(
            assay_id="ribosome_profiling_denovo",
            name="Ribosome Profiling - De Novo Detected",
            description="Read structure detected using de novo multi-signal analysis",
            library_protocol="unknown_detected",
            library_kit="auto-detected"
        )
        
        regions = []
        
        # Convert detected segments to seqspec regions
        for segment in detected_segments:
            # Extract consensus sequence for this segment
            segment_sequences = []
            for read in sample_reads[:50]:
                if segment.end_pos <= len(read):
                    seg_seq = read[segment.start_pos:segment.end_pos]
                    segment_sequences.append(seg_seq)
            
            consensus = self._get_consensus_sequence(segment_sequences)
            
            # Map segment types to seqspec region types
            region_type_map = {
                "umi": "umi",
                "barcode": "barcode", 
                "adapter": "illumina_p7",
                "rpf": "cdna",
                "unknown": "unknown"
            }
            
            sequence_type_map = {
                "umi": "random",
                "barcode": "onlist",
                "adapter": "fixed",
                "rpf": "joined",
                "unknown": "random"
            }
            
            region = SeqSpecRegion(
                region_id=f"{segment.segment_type}_{segment.start_pos}_{segment.end_pos}",
                region_type=region_type_map.get(segment.segment_type, "unknown"),
                sequence_type=sequence_type_map.get(segment.segment_type, "random"),
                name=f"{segment.segment_type.title()} region ({segment.start_pos}-{segment.end_pos})",
                sequence=consensus,
                min_len=segment.end_pos - segment.start_pos,
                max_len=segment.end_pos - segment.start_pos
            )
            regions.append(region)
        
        # If no segments detected, create a single RPF region
        if not regions:
            read_lengths = [len(read) for read in sample_reads[:100]]
            avg_len = sum(read_lengths) // len(read_lengths) if read_lengths else 30
            
            regions.append(SeqSpecRegion(
                region_id="full_rpf_sequence",
                region_type="cdna",
                sequence_type="joined",
                name="Complete Ribosome Protected Fragment",
                min_len=avg_len - 5,
                max_len=avg_len + 5
            ))
        
        # Create read specification
        read_lengths = [len(read) for read in sample_reads[:100]]
        min_read_len = min(read_lengths) if read_lengths else 20
        max_read_len = max(read_lengths) if read_lengths else 50
        
        read_spec = SeqSpecRead(
            read_id="R1",
            name="Read 1", 
            modality="rna",
            primer_id="unknown_primer",
            strand="pos",
            min_len=min_read_len,
            max_len=max_read_len,
            files=[{"file_id": "R1.fastq.gz", "filename": "input.fastq"}]
        )
        
        assay.sequence_spec = regions
        
        seqspec_dict = self._create_seqspec_dict(assay, [read_spec])
        
        if output_file:
            self._write_yaml(seqspec_dict, output_file)
            logger.info(f"seqspec file written to {output_file}")
        
        return seqspec_dict
    
    def _extract_consensus(self, sample_reads: List[str], start: int, end: int) -> Optional[str]:
        """Extract consensus sequence from a region across sample reads."""
        if start >= end or not sample_reads:
            return None
        
        sequences = []
        for read in sample_reads[:50]:  # Sample first 50 reads
            if end <= len(read):
                sequences.append(read[start:end])
        
        return self._get_consensus_sequence(sequences)
    
    def _get_consensus_sequence(self, sequences: List[str]) -> Optional[str]:
        """Get consensus sequence from list of sequences."""
        if not sequences:
            return None
        
        # Check if all sequences are identical (fixed sequence)
        if len(set(sequences)) == 1:
            return sequences[0]
        
        # Check if it looks random (high diversity)
        if len(set(sequences)) / len(sequences) > 0.8:
            return "N" * len(sequences[0]) if sequences else None
        
        # Build position-wise consensus
        if not sequences or not sequences[0]:
            return None
        
        consensus = []
        seq_len = len(sequences[0])
        
        for pos in range(seq_len):
            nucleotides = [seq[pos] for seq in sequences if pos < len(seq)]
            if nucleotides:
                # Find most common nucleotide
                from collections import Counter
                counts = Counter(nucleotides)
                most_common = counts.most_common(1)[0][0]
                
                # If >70% consensus, use it; otherwise use ambiguous code
                if counts[most_common] / len(nucleotides) > 0.7:
                    consensus.append(most_common)
                else:
                    consensus.append("N")
            else:
                consensus.append("N")
        
        return "".join(consensus)
    
    def _create_seqspec_dict(self, assay: SeqSpecAssay, reads: List[SeqSpecRead]) -> Dict[str, Any]:
        """Create complete seqspec dictionary structure."""
        
        # Convert regions to dict format
        def region_to_dict(region: SeqSpecRegion) -> Dict[str, Any]:
            region_dict = {
                "!Region": None,
                "region_id": region.region_id,
                "region_type": region.region_type,
                "sequence_type": region.sequence_type,
                "name": region.name,
                "min_len": region.min_len,
                "max_len": region.max_len,
                "strand": region.strand
            }
            
            if region.sequence:
                region_dict["sequence"] = region.sequence
            
            if region.regions:
                region_dict["regions"] = [region_to_dict(r) for r in region.regions]
            
            return region_dict
        
        # Convert reads to dict format
        def read_to_dict(read: SeqSpecRead) -> Dict[str, Any]:
            return {
                "!Read": None,
                "read_id": read.read_id,
                "name": read.name,
                "modality": read.modality,
                "primer_id": read.primer_id,
                "strand": read.strand,
                "min_len": read.min_len,
                "max_len": read.max_len,
                "files": [{"!File": None, **f} for f in read.files]
            }
        
        seqspec = {
            "!Assay": None,
            "seqspec_version": assay.seqspec_version,
            "assay_id": assay.assay_id,
            "name": assay.name,
            "doi": assay.doi,
            "publication_date": assay.publication_date,
            "description": assay.description,
            "modalities": assay.modalities,
            "lib_struct": assay.lib_struct,
            "library_protocol": assay.library_protocol,
            "library_kit": assay.library_kit,
            "sequence_protocol": assay.sequence_protocol,
            "sequence_kit": assay.sequence_kit,
            "sequence_spec": [region_to_dict(r) for r in assay.sequence_spec],
            "library_spec": [region_to_dict(r) for r in assay.library_spec]
        }
        
        if reads:
            seqspec["reads"] = [read_to_dict(r) for r in reads]
        
        return seqspec
    
    def _write_yaml(self, seqspec_dict: Dict[str, Any], output_file: Path):
        """Write seqspec dictionary to YAML file with proper tags."""

        class TaggedDict(dict):
            """Dictionary subclass that holds a YAML tag."""
            yaml_tag = None

        def dict_representer(dumper, data):
            """Represent tagged dictionaries with proper YAML tags."""
            if isinstance(data, TaggedDict) and data.yaml_tag:
                return dumper.represent_mapping(data.yaml_tag, dict(data))
            return dumper.represent_mapping('tag:yaml.org,2002:map', data)

        # Add representer for TaggedDict
        yaml.add_representer(TaggedDict, dict_representer)

        def convert_to_tagged(obj: Any) -> Any:
            """Recursively convert dicts with !Tag markers to TaggedDict."""
            if isinstance(obj, dict):
                tag = obj.get('!Assay') or obj.get('!Region') or obj.get('!Read') or obj.get('!File')
                clean_dict = {k: convert_to_tagged(v) for k, v in obj.items() if not k.startswith('!')}

                if tag is not None:
                    # Extract the tag name
                    for key in obj.keys():
                        if key.startswith('!'):
                            tagged = TaggedDict(clean_dict)
                            tagged.yaml_tag = key
                            return tagged

                return clean_dict
            elif isinstance(obj, list):
                return [convert_to_tagged(item) for item in obj]
            else:
                return obj

        # Convert the seqspec dict to use TaggedDict
        converted = convert_to_tagged(seqspec_dict)

        with open(output_file, 'w') as f:
            yaml.dump(
                converted,
                f,
                default_flow_style=False,
                indent=2,
                sort_keys=False,
                allow_unicode=True
            )

        logger.info(f"seqspec YAML written to {output_file}")
    
    def generate_mixed_seqspec(
        self,
        clean_reads: List[str],
        contaminated_reads: List[str], 
        architecture,
        output_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Generate seqspec for mixed clean/contaminated data."""
        
        assay = SeqSpecAssay(
            assay_id="ribosome_profiling_mixed",
            name="Ribosome Profiling - Mixed Clean/Contaminated",
            description=f"Mixed dataset: {len(clean_reads)} clean RPFs + {len(contaminated_reads)} contaminated reads",
            library_protocol="mixed_processing_states"
        )
        
        regions = []
        
        # Define two main structural variants
        # Variant 1: Clean RPFs (no preprocessing needed)
        clean_variant = SeqSpecRegion(
            region_id="clean_rpf_variant",
            region_type="cdna",
            sequence_type="joined",
            name="Clean RPF sequences (pre-processed)",
            min_len=min(len(r) for r in clean_reads) if clean_reads else 25,
            max_len=max(len(r) for r in clean_reads) if clean_reads else 40
        )
        
        # Variant 2: Contaminated reads needing processing
        if architecture and contaminated_reads:
            contam_regions = []
            
            # Add adapter region for contaminated variant
            for adapter in architecture.adapter_sequences:
                contam_regions.append(SeqSpecRegion(
                    region_id=f"contaminant_adapter",
                    region_type="illumina_p7",
                    sequence_type="fixed",
                    name="Adapter contamination (needs removal)",
                    sequence=adapter,
                    min_len=len(adapter),
                    max_len=len(adapter)
                ))
            
            # Add RPF region for contaminated reads
            contam_regions.append(SeqSpecRegion(
                region_id="contaminated_rpf",
                region_type="cdna", 
                sequence_type="joined",
                name="RPF sequence (adapter-contaminated)",
                min_len=25,
                max_len=35
            ))
            
            contaminated_variant = SeqSpecRegion(
                region_id="contaminated_variant",
                region_type="joined",
                sequence_type="joined", 
                name="Contaminated reads requiring processing",
                regions=contam_regions
            )
            
            regions = [clean_variant, contaminated_variant]
        else:
            regions = [clean_variant]
        
        assay.sequence_spec = regions
        
        # Create read spec for mixed data
        all_reads = clean_reads + contaminated_reads
        read_lengths = [len(read) for read in all_reads[:100]]
        
        read_spec = SeqSpecRead(
            read_id="R1_mixed",
            name="Mixed Read 1",
            modality="rna", 
            primer_id="mixed_primers",
            strand="pos",
            min_len=min(read_lengths) if read_lengths else 20,
            max_len=max(read_lengths) if read_lengths else 50,
            files=[{"file_id": "R1.fastq.gz", "filename": "mixed_input.fastq"}]
        )
        
        seqspec_dict = self._create_seqspec_dict(assay, [read_spec])
        
        if output_file:
            self._write_yaml(seqspec_dict, output_file)
        
        return seqspec_dict
    
    def _has_structure_annotations(self, headers: List[str]) -> bool:
        """Check if headers contain structure annotations."""
        if not headers:
            return False
        
        # Look for structure patterns in headers
        structure_keywords = ['umi:', 'barcode:', 'adapter:', 'spacer:', 'rpf:', 'untemplated:']
        
        for header in headers[:10]:  # Check first 10 headers
            if any(keyword in header.lower() for keyword in structure_keywords):
                return True
        
        return False
    
    def _generate_from_annotated_headers(
        self,
        sample_headers: List[str],
        sample_reads: List[str],
        output_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Generate seqspec from structure-annotated headers."""
        
        # Parse structure information from headers
        structure_info = self._parse_structure_annotations(sample_headers, sample_reads)
        
        assay = SeqSpecAssay(
            assay_id="ribosome_profiling_annotated",
            name="Ribosome Profiling - Structure Annotated",
            description=f"Complex read structure parsed from annotations. {len(structure_info['regions'])} distinct regions detected.",
            library_protocol="complex_annotated_structure"
        )
        
        # Build regions from parsed structure
        regions = []
        for region_info in structure_info['regions']:
            region = SeqSpecRegion(
                region_id=region_info['id'],
                region_type=self._map_annotation_to_seqspec_type(region_info['type']),
                sequence_type=self._infer_sequence_type(region_info),
                name=region_info['name'],
                sequence=region_info.get('consensus_sequence'),
                min_len=region_info['min_len'],
                max_len=region_info['max_len'],
                strand="pos"
            )
            regions.append(region)
        
        assay.sequence_spec = regions
        
        # Create read specification
        read_lengths = [len(read) for read in sample_reads[:100]]
        read_spec = SeqSpecRead(
            read_id="R1_annotated",
            name="Annotated Read 1",
            modality="rna",
            primer_id="complex_primer",
            strand="pos", 
            min_len=min(read_lengths) if read_lengths else 50,
            max_len=max(read_lengths) if read_lengths else 100,
            files=[{"file_id": "R1.fastq.gz", "filename": "annotated_input.fastq"}]
        )
        
        seqspec_dict = self._create_seqspec_dict(assay, [read_spec])
        
        if output_file:
            self._write_yaml(seqspec_dict, output_file)
            logger.info(f"seqspec file written to {output_file}")
        
        return seqspec_dict
    
    def _parse_structure_annotations(self, headers: List[str], reads: List[str]) -> Dict[str, Any]:
        """Parse structure annotations from FASTQ headers.
        
        Expected format: @read_id_element1:length_element2:length_...
        """
        region_data = {}  # region_type -> {positions: [], lengths: [], sequences: []}
        total_structures = []
        
        for i, header in enumerate(headers[:100]):  # Sample first 100
            if i >= len(reads):
                break
                
            # Parse structure from header
            structure = self._extract_structure_from_header(header)
            if structure:
                total_structures.append(structure)
                
                # Extract sequences for each region
                read_seq = reads[i]
                current_pos = 0
                
                for region_type, length in structure:
                    if region_type not in region_data:
                        region_data[region_type] = {
                            'positions': [],
                            'lengths': [],
                            'sequences': []
                        }
                    
                    end_pos = current_pos + length
                    if end_pos <= len(read_seq):
                        region_seq = read_seq[current_pos:end_pos]
                        region_data[region_type]['positions'].append((current_pos, end_pos))
                        region_data[region_type]['lengths'].append(length)
                        region_data[region_type]['sequences'].append(region_seq)
                    
                    current_pos = end_pos
        
        # Build summary regions
        regions = []
        for region_type, data in region_data.items():
            if data['lengths']:
                consensus_seq = self._get_consensus_sequence(data['sequences'])
                
                region_info = {
                    'id': f"region_{region_type}",
                    'type': region_type,
                    'name': f"{region_type.title()} Region",
                    'min_len': min(data['lengths']),
                    'max_len': max(data['lengths']),
                    'consensus_sequence': consensus_seq,
                    'occurrence_count': len(data['lengths']),
                    'sequences': data['sequences'][:10]  # Sample sequences
                }
                regions.append(region_info)
        
        return {
            'regions': regions,
            'total_reads_analyzed': len(total_structures),
            'unique_structures': len(set(str(s) for s in total_structures))
        }
    
    def _extract_structure_from_header(self, header: str) -> List[Tuple[str, int]]:
        """Extract structure information from FASTQ header.
        
        Args:
            header: FASTQ header like '@read_0_umi:8_spacer1:4_barcode:6_adapter1:12_rpf:30'
            
        Returns:
            List of (element_type, length) tuples
        """
        import re
        
        # Find all element:length patterns
        pattern = r'([a-zA-Z_][a-zA-Z0-9_]*):(\d+)'
        matches = re.findall(pattern, header)
        
        structure = []
        for element_type, length_str in matches:
            try:
                length = int(length_str)
                structure.append((element_type.lower(), length))
            except ValueError:
                continue
        
        return structure
    
    def _map_annotation_to_seqspec_type(self, annotation_type: str) -> str:
        """Map annotation type to seqspec region type."""
        mapping = {
            'umi': 'umi',
            'barcode': 'barcode',
            'adapter': 'illumina_p7',
            'adapter1': 'illumina_p5',
            'adapter2': 'illumina_p7',
            'adapter_5': 'illumina_p5',
            'adapter_3': 'illumina_p7',
            'spacer': 'linker',
            'spacer1': 'linker',
            'spacer2': 'linker',
            'rpf': 'cdna',
            'untemplated': 'poly_tail',
            'poly_a': 'poly_tail',
            'poly_t': 'poly_tail',
            'random_tail': 'poly_tail',
            'umi_partial': 'umi',
            'adapter_partial': 'illumina_p7'
        }
        
        return mapping.get(annotation_type.lower(), 'unknown')
    
    def _infer_sequence_type(self, region_info: Dict) -> str:
        """Infer sequence type from region characteristics."""
        region_type = region_info['type'].lower()
        sequences = region_info.get('sequences', [])
        
        if not sequences:
            return 'unknown'
        
        # Calculate sequence diversity
        unique_sequences = len(set(sequences))
        total_sequences = len(sequences)
        diversity = unique_sequences / total_sequences if total_sequences > 0 else 0
        
        # Type-specific inference
        if region_type in ['umi', 'umi_partial']:
            return 'random' if diversity > 0.7 else 'onlist'
        elif region_type in ['barcode']:
            return 'onlist' if diversity < 0.5 else 'random'
        elif region_type in ['adapter', 'adapter1', 'adapter2', 'adapter_5', 'adapter_3', 'adapter_partial']:
            return 'fixed' if diversity < 0.3 else 'joined'
        elif region_type in ['spacer', 'spacer1', 'spacer2']:
            return 'fixed' if diversity < 0.5 else 'random'
        elif region_type in ['rpf']:
            return 'joined'  # Biological sequence
        elif region_type in ['untemplated', 'poly_a', 'poly_t', 'random_tail']:
            return 'random'
        else:
            return 'joined'