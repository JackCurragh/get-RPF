"""Put the answers together (spec §5.5).

Runs steps 0-3 on the bounded sample, answers Q4, decides whether a transform
may be emitted, and builds the architecture that extraction will apply.
There is no global score: transform confidence is bounded by the
least-resolved question that changes emitted bases (Q1-Q3).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, Optional, Sequence, Tuple

from .anchors import AnchorCall, find_anchor, search_window
from .config import InferenceConfig
from .model import (
    Alternative,
    Answer,
    Architecture,
    Block,
    Evidence,
    JunctionCall,
    Status,
    TransformDecision,
)
from .observe import Observation, observe
from .pileup import FrameProfile, JunctionInference, infer_junctions

_ACCEPTABLE = (Status.RESOLVED, Status.INTERVAL)


@dataclass(frozen=True)
class StructureInference:
    observation: Observation
    q1: Answer
    junctions: JunctionInference
    q4: Answer
    window: int
    transform: TransformDecision
    architecture: Optional[Architecture]

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

    # Trimmed reads: the read end is the anchor. Raw reads without an anchor:
    # only the read-start frame is meaningful, and Q3 is overridden below.
    anchors: Optional[List[Optional[int]]] = (
        list(call.insert_ends) if call is not None else None
    )
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

    q4 = answer_umi(observation, junctions.q2, junctions.q3)
    transform = decide_transform(observation, q1, junctions.q2, junctions.q3, config)
    architecture = build_architecture(q1, junctions, transform, config)
    return StructureInference(
        observation, q1, junctions, q4, window, transform, architecture
    )


def answer_umi(observation: Observation, q2: Answer, q3: Answer) -> Answer:
    """Q4: which bases to keep as a UMI. Absence is claimed only where every
    channel that could carry one was examined (spec §4)."""
    inline = tuple(
        (answer.question, length)
        for answer in (q2, q3)
        if isinstance(answer.value, JunctionCall)
        for kind, length in answer.value.blocks
        if kind == "random"
    )
    header = observation.header_umi
    evidence = (
        Evidence("blocks", "inline_random_blocks", inline, observation.reads),
        Evidence("profile", "header_umi", header, observation.reads),
    )
    value = {"inline": inline, "header": header}
    if inline or header:
        parts = [
            f"{length} nt random "
            f"{'at the read start' if question == 'Q2' else 'before the anchor'}"
            for question, length in inline
        ]
        if header:
            parts.append(f"a UMI in the read names ({header})")
        return Answer(
            "Q4",
            value,
            Status.RESOLVED,
            evidence,
            (),
            "Keep as UMI: " + "; ".join(parts) + ".",
        )
    if observation.input_state == "trimmed":
        return Answer(
            "Q4",
            None,
            Status.NOT_OBSERVABLE,
            evidence,
            (Alternative("absent", "untested", "reads were trimmed before deposit"),),
            "No UMI can be seen, but the reads were trimmed before deposit, "
            "which may have removed one.",
        )
    if q2.status in _ACCEPTABLE and q3.status in _ACCEPTABLE:
        return Answer(
            "Q4",
            value,
            Status.RESOLVED,
            evidence,
            (),
            "No inline UMI in R1 and none in the read names. Index reads were "
            "not supplied, so a UMI in an index read cannot be excluded.",
        )
    return Answer(
        "Q4",
        None,
        _least_resolved((q2.status, q3.status)),
        evidence,
        (),
        "Whether R1 carries a UMI is unresolved because Q2 or Q3 is unresolved.",
    )


def decide_transform(
    observation: Observation,
    q1: Answer,
    q2: Answer,
    q3: Answer,
    config: InferenceConfig,
) -> TransformDecision:
    """Spec §7.2: emit only when every question that changes emitted bases is
    resolved, or an interval under a named convention."""
    reasons: List[str] = []
    conventions: List[str] = []
    q1_status = q1.status
    if q1.status is Status.NOT_OBSERVABLE and observation.input_state == "trimmed":
        q1_status = Status.RESOLVED
        conventions.append("trimmed input: the read end is the anchor")
    if isinstance(q1.value, AnchorCall) and q1.value.tail_base is not None:
        conventions.append(
            f"poly({q1.value.tail_base}) tail: the insert ends at the first base "
            "of the run (spec §6)"
        )
    statuses = (q1_status, q2.status, q3.status)
    for question, status in zip(("Q1", "Q2", "Q3"), statuses):
        if status not in _ACCEPTABLE:
            reasons.append(f"{question} is {status.value}")
    # Junction bases never withhold the transform: they stay in the insert, so
    # at most a base or two per read is at stake. Common ones are flagged.
    flags: List[str] = []
    for answer in (q2, q3):
        call = answer.value
        if not isinstance(call, JunctionCall) or not call.nta_length:
            continue
        rate = max((r for r in call.nta_rates if r is not None), default=0.0)
        where = "read start" if call.frame == "read_start" else "anchor"
        conventions.append(
            f"{answer.question}: {call.nta_length} junction base(s) at the "
            f"{where} kept in the insert (non-templated-like in about "
            f"{rate:.0%} of reads; v1 does not remove them)"
        )
        if rate >= config.nta_flag_rate:
            flags.append(
                f"{answer.question}: the {call.nta_length} kept junction base(s) "
                f"at the {where} look non-templated in about {rate:.0%} of reads, "
                f"so most inserts carry up to {call.nta_length} extra base(s) there"
            )
    return TransformDecision(
        emit=not reasons,
        bound=_least_resolved(statuses),
        reasons=tuple(reasons),
        conventions=tuple(conventions),
        flags=tuple(flags),
    )


def build_architecture(
    q1: Answer,
    junctions: JunctionInference,
    transform: TransformDecision,
    config: InferenceConfig,
) -> Optional[Architecture]:
    """Read-order blocks from the answers; None when Q2 or Q3 has no value."""
    q2, q3 = junctions.q2, junctions.q3
    if not isinstance(q2.value, JunctionCall) or not isinstance(q3.value, JunctionCall):
        return None
    five = _frame_blocks(q2.value, junctions.pileup.read_start, q2.status)
    three = _frame_blocks(q3.value, junctions.pileup.anchor, q3.status)
    insert = Block(
        "insert",
        "read_start",
        (config.fragment_min, config.fragment_max),
        None,
        keep_as_umi=False,
        remove=False,
        status=transform.bound,
        note=f"accepted under {config.fragment_policy}",
    )
    blocks = [*five, insert, *reversed(three)]
    call = q1.value if isinstance(q1.value, AnchorCall) else None
    if call is not None:
        if call.tail_base is not None:
            runs = [
                start - end
                for start, end in zip(call.adapter_starts, call.insert_ends)
                if start is not None and end is not None
            ]
            blocks.append(
                Block(
                    "tail",
                    "anchor",
                    (0, max(runs, default=0)),
                    call.tail_base,
                    keep_as_umi=False,
                    remove=True,
                    status=Status.INTERVAL,
                    note="insert end placed at the first base of the run",
                )
            )
        length = len(call.adapter.sequence)
        blocks.append(
            Block(
                "adapter",
                "anchor",
                (length, length),
                call.adapter.sequence,
                keep_as_umi=False,
                remove=True,
                status=q1.status,
                note=f"{call.adapter.name} ({call.source})",
            )
        )
    return Architecture(
        tuple(blocks), "inferred", transform.bound, config.fragment_policy
    )


def _frame_blocks(
    call: JunctionCall, profile: FrameProfile, status: Status
) -> List[Block]:
    """One frame's blocks, ordered from the frame edge inward.

    In the anchor frame that is the reverse of read order; the caller reverses
    the list, and fixed sequences are written here in read order.
    """
    blocks: List[Block] = []
    position = 0
    for kind, length in call.blocks:
        sequence: Optional[str] = None
        if kind == "fixed":
            bases = [
                max(stats.composition, key=lambda base: stats.composition[base])
                for stats in profile.positions[position : position + length]
            ]
            sequence = "".join(bases if call.frame == "read_start" else reversed(bases))
        blocks.append(
            Block(
                (
                    "random"
                    if kind == "random"
                    else ("fixed" if kind == "fixed" else "barcode")
                ),
                call.frame,
                (length, length),
                sequence,
                keep_as_umi=kind == "random",
                remove=True,
                status=status,
            )
        )
        position += length
    if call.nta_length:
        rate = max((r for r in call.nta_rates if r is not None), default=0.0)
        blocks.append(
            Block(
                "nta",
                call.frame,
                (0, call.nta_length),
                None,
                keep_as_umi=False,
                remove=False,
                status=Status.INTERVAL,
                note=f"non-templated-like in about {rate:.0%} of reads; kept (v1)",
            )
        )
    return blocks


def _least_resolved(statuses: Tuple[Status, ...]) -> Status:
    for status in statuses:
        if status not in _ACCEPTABLE:
            return status
    return Status.INTERVAL if Status.INTERVAL in statuses else Status.RESOLVED
