"""Put the answers together (spec §5.5).

So far this runs steps 0-3 in order on the bounded sample: observe, find the
anchor (Q1), size the search window from it, and answer Q2/Q3 with the
pileup relative to the read start and the insert-side anchor.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, Optional, Sequence

from .anchors import AnchorCall, find_anchor, search_window
from .config import InferenceConfig
from .model import Alternative, Answer, Evidence, Status
from .observe import Observation, observe
from .pileup import JunctionInference, infer_junctions


@dataclass(frozen=True)
class StructureInference:
    observation: Observation
    q1: Answer
    junctions: JunctionInference
    window: int

    @property
    def q2(self) -> Answer:
        return self.junctions.q2

    @property
    def q3(self) -> Answer:
        return self.junctions.q3


def infer_structure(
    reads: Sequence[str],
    headers: Sequence[str] = (),
    config: Optional[InferenceConfig] = None,
) -> StructureInference:
    config = config or InferenceConfig()
    reads = list(reads[: config.sample_reads])
    observation = observe(reads, headers, config)
    q1 = find_anchor(reads, observation, config)
    call = q1.value if isinstance(q1.value, AnchorCall) else None

    anchors: Optional[List[Optional[int]]]
    if call is not None:
        anchors = list(call.insert_ends)
    else:
        # Trimmed reads: the read end is the anchor. Raw reads without an
        # anchor: only the read-start frame is meaningful (Q3 is overridden).
        anchors = None
    window = search_window(call, observation, config)
    junctions = infer_junctions(reads, anchors, config, window=window)

    if call is None and observation.input_state != "trimmed":
        junctions = replace(
            junctions,
            q3=Answer(
                "Q3",
                None,
                Status.NOT_OBSERVABLE,
                (Evidence("anchor", "q1_status", q1.status.value, len(reads)),),
                (Alternative(None, "untested", "no 3' anchor in raw reads"),),
                "The insert's 3' side cannot be examined: raw reads with no "
                "located anchor end at an unknown point in the construct.",
            ),
        )
    return StructureInference(observation, q1, junctions, window)
