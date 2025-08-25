#!/usr/bin/env python3

"""
seqspec-based Architecture Loader
Allows loading architectures from seqspec files in a directory
"""

import yaml
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class ReadArchitecture:
    """Represents a read architecture loaded from seqspec."""
    protocol_name: str
    lab_source: str
    umi_positions: List[tuple] = None
    barcode_positions: List[tuple] = None
    adapter_sequences: List[str] = None
    rpf_start: int = 0
    rpf_end: int = -1
    expected_rpf_length: tuple = (20, 40)
    quality_markers: Dict[str, Any] = None

    def __post_init__(self):
        if self.umi_positions is None:
            self.umi_positions = []
        if self.barcode_positions is None:
            self.barcode_positions = []
        if self.adapter_sequences is None:
            self.adapter_sequences = []
        if self.quality_markers is None:
            self.quality_markers = {}

class SeqSpecArchitectureLoader:
    """Loads read architectures from seqspec files."""
    
    def __init__(self):
        self.loaded_architectures = []
    
    def load_from_directory(self, seqspec_dir: Path) -> List[ReadArchitecture]:
        """Load all seqspec files from a directory."""
        
        seqspec_files = list(seqspec_dir.glob("*.yaml")) + list(seqspec_dir.glob("*.yml"))
        
        logger.info(f"Found {len(seqspec_files)} seqspec files in {seqspec_dir}")
        
        architectures = []
        for seqspec_file in seqspec_files:
            try:
                arch = self.load_from_seqspec(seqspec_file)
                if arch:
                    architectures.append(arch)
                    logger.info(f"Loaded architecture: {arch.protocol_name}")
            except Exception as e:
                logger.warning(f"Failed to load {seqspec_file}: {e}")
        
        self.loaded_architectures.extend(architectures)
        return architectures
    
    def load_from_seqspec(self, seqspec_file: Path) -> Optional[ReadArchitecture]:
        """Load a single architecture from a seqspec file."""
        
        logger.debug(f"Loading seqspec: {seqspec_file}")
        
        try:
            with open(seqspec_file, 'r') as f:
                content = f.read()
            
            # Try to parse as YAML (ignoring custom tags for now)
            # This is a simplified parser - real seqspec parsing would be more complex
            yaml_content = self._clean_seqspec_yaml(content)
            data = yaml.safe_load(yaml_content)
            
            if not data:
                logger.warning(f"Empty or invalid YAML in {seqspec_file}")
                return None
            
            # Extract architecture information
            protocol_name = data.get('assay_id', seqspec_file.stem)
            lab_source = data.get('library_kit', 'Unknown')
            description = data.get('description', '')
            
            # Parse sequence_spec for regions
            sequence_spec = data.get('sequence_spec', [])
            
            adapter_sequences = []
            umi_positions = []
            barcode_positions = []
            rpf_regions = []
            
            current_pos = 0
            
            for region in sequence_spec:
                if isinstance(region, dict):
                    region_type = region.get('region_type', '').lower()
                    region_id = region.get('region_id', '').lower()
                    sequence = region.get('sequence', '')
                    min_len = region.get('min_len', 0)
                    max_len = region.get('max_len', min_len)
                    
                    region_length = max(min_len, max_len, len(sequence)) if sequence else max_len
                    
                    # Categorize regions
                    if 'adapter' in region_type or 'adapter' in region_id:
                        if sequence:
                            adapter_sequences.append(sequence)
                    elif 'umi' in region_type or 'umi' in region_id:
                        if region_length > 0:
                            umi_positions.append((current_pos, current_pos + region_length))
                    elif 'barcode' in region_type or 'barcode' in region_id:
                        if region_length > 0:
                            barcode_positions.append((current_pos, current_pos + region_length))
                    elif 'rpf' in region_type or 'rpf' in region_id or 'cdna' in region_type:
                        rpf_regions.append((current_pos, current_pos + region_length, min_len, max_len))
                    
                    current_pos += region_length
            
            # Determine RPF length range
            if rpf_regions:
                rpf_min = min(region[2] for region in rpf_regions)
                rpf_max = max(region[3] for region in rpf_regions)
                expected_rpf_length = (rpf_min, rpf_max)
            else:
                expected_rpf_length = (20, 40)  # Default
            
            # Create architecture
            architecture = ReadArchitecture(
                protocol_name=protocol_name,
                lab_source=lab_source,
                umi_positions=umi_positions,
                barcode_positions=barcode_positions,
                adapter_sequences=adapter_sequences,
                rpf_start=0,  # Will be calculated based on positions
                rpf_end=-1,   # Will be calculated based on adapters
                expected_rpf_length=expected_rpf_length,
                quality_markers={
                    "adapter_match_threshold": 0.3,  # Conservative for user-defined
                    "description": description
                }
            )
            
            logger.info(f"Created architecture from seqspec: {protocol_name}")
            logger.debug(f"  Adapters: {adapter_sequences}")
            logger.debug(f"  UMI positions: {umi_positions}")
            logger.debug(f"  Barcode positions: {barcode_positions}")
            logger.debug(f"  RPF length: {expected_rpf_length}")
            
            return architecture
            
        except Exception as e:
            logger.error(f"Error loading seqspec {seqspec_file}: {e}")
            return None
    
    def _clean_seqspec_yaml(self, content: str) -> str:
        """Clean seqspec YAML content to make it parseable by standard YAML."""
        
        # Remove custom Python object tags
        lines = content.split('\n')
        cleaned_lines = []
        
        for line in lines:
            # Skip lines with Python object tags
            if '!!python/object/apply:' in line:
                continue
            # Convert custom tags to simple strings
            if line.strip().startswith('!!'):
                continue
            cleaned_lines.append(line)
        
        return '\n'.join(cleaned_lines)

    def create_sample_seqspec(self, output_dir: Path):
        """Create sample seqspec files for testing."""
        
        output_dir.mkdir(exist_ok=True)
        
        # Sample seqspec 1: Novel protocol
        sample1 = {
            'seqspec_version': '0.3.0',
            'assay_id': 'novel_protocol_2024',
            'name': 'Novel Ribosome Profiling Protocol',
            'description': 'Custom protocol with unique barcode structure',
            'library_kit': 'Custom Lab Protocol',
            'modalities': ['rna'],
            'sequence_spec': [
                {
                    'region_id': 'custom_prefix',
                    'region_type': 'technical',
                    'sequence': 'TAG',
                    'min_len': 3,
                    'max_len': 3
                },
                {
                    'region_id': 'sample_barcode',
                    'region_type': 'barcode',
                    'sequence': 'NNNNNNNN',
                    'min_len': 8,
                    'max_len': 8
                },
                {
                    'region_id': 'ribosome_protected_fragment',
                    'region_type': 'cdna',
                    'sequence_type': 'joined',
                    'min_len': 20,
                    'max_len': 35
                },
                {
                    'region_id': 'custom_adapter',
                    'region_type': 'adapter',
                    'sequence': 'GGCCTTAAGCCCGGAA',
                    'min_len': 16,
                    'max_len': 16
                }
            ]
        }
        
        with open(output_dir / 'novel_protocol.yaml', 'w') as f:
            yaml.dump(sample1, f, indent=2)
        
        # Sample seqspec 2: Modified McGlincy protocol
        sample2 = {
            'seqspec_version': '0.3.0',
            'assay_id': 'modified_mcglincy_2024',
            'name': 'Modified McGlincy Protocol',
            'description': 'McGlincy protocol with extended UMI',
            'library_kit': 'Modified McGlincy Lab',
            'modalities': ['rna'],
            'sequence_spec': [
                {
                    'region_id': 'extended_umi',
                    'region_type': 'umi',
                    'sequence': 'NNNNNNNN',
                    'min_len': 8,
                    'max_len': 8
                },
                {
                    'region_id': 'sample_barcode',
                    'region_type': 'barcode',
                    'sequence': 'NNNNN',
                    'min_len': 5,
                    'max_len': 5
                },
                {
                    'region_id': 'rpf_sequence',
                    'region_type': 'cdna',
                    'min_len': 24,
                    'max_len': 36
                },
                {
                    'region_id': 'illumina_adapter',
                    'region_type': 'adapter',
                    'sequence': 'AGATCGGAAGAGCAC',
                    'min_len': 15,
                    'max_len': 15
                }
            ]
        }
        
        with open(output_dir / 'modified_mcglincy.yaml', 'w') as f:
            yaml.dump(sample2, f, indent=2)
        
        logger.info(f"Created sample seqspec files in {output_dir}")
        
        return [output_dir / 'novel_protocol.yaml', output_dir / 'modified_mcglincy.yaml']