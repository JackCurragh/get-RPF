"""Read-structure inference (docs/read_structure_inference_spec.md).

Implemented so far: the library-as-reference pileup and the Q2/Q3 junction
calls it supports (spec §11, phase P3 -- the go/no-go gate).
"""

from .config import InferenceConfig
from .model import Alternative, Answer, Evidence, Frame, JunctionCall, Status
from .pileup import (
    FrameProfile,
    JunctionInference,
    Pileup,
    PositionStats,
    build_pileup,
    call_junction,
    infer_junctions,
)

__all__ = [
    "Alternative",
    "Answer",
    "Evidence",
    "Frame",
    "FrameProfile",
    "InferenceConfig",
    "JunctionCall",
    "JunctionInference",
    "Pileup",
    "PositionStats",
    "Status",
    "build_pileup",
    "call_junction",
    "infer_junctions",
]
