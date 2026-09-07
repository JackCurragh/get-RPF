"""The M3 apply path: turn M2 boundary decisions into override_trims plus
provenance records, reference-free (no alignment required).

See docs/release_qc_and_terminal_trimming_plan.md sections 4 (Stage 1) and
M3. This is what connects `estimate_boundaries` (boundary.py) to
`RPFExtractor.extract_rpfs`'s `override_trims` parameter.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .processors.boundary import LengthClassDecision, estimate_boundaries
from .processors.sketch import Sketch, build_sketch

OverrideTrims = Dict[int, Dict[str, int]]


DecisionsByLength = Dict[int, LengthClassDecision]

# Every trim-rule record (whether boundary-derived or from an explicit
# override) shares this key set, so downstream consumers (evidence object,
# cohort TSVs) never have to guard for a field being present on one path
# and missing on another.
RULE_RECORD_KEYS = (
    "end",
    "read_length_before",
    "trim_bases",
    "dominant_sequence",
    "dominant_fraction",
    "supporting_reads",
    "reason",
    "estimators_agree",
    "concordance",
    "confidence",
    "decision",
)


def plan_from_sketch(
    sketch: Sketch,
    known_adapters: Optional[List[str]] = None,
) -> Tuple[OverrideTrims, List[Dict], List[Dict], DecisionsByLength]:
    """Compute the 3' apply plan from a Stage-0 sketch.

    5' trimming is intentionally not attempted here: it requires a named
    element (seqspec or high-posterior HMM segment), which this reference-
    free, sketch-only entry point does not have. See section 2's asymmetric
    5'/3' risk principle.

    Returns (override_trims, applied_rules, proposed_rules, decisions). The
    rule dicts match the trim_rules_applied / trim_rules_proposed_not_applied
    shape from the evidence object (section 9); `decisions` is the full
    per-length LengthClassDecision set, for populating boundary_estimates_3p.
    """
    decisions = estimate_boundaries(
        sketch.per_length,
        sketch.per_length_support,
        end="3p",
        per_length_reads=sketch.per_length_reads,
        known_elements=known_adapters,
    )

    override_trims: OverrideTrims = {}
    applied_rules: List[Dict] = []
    proposed_rules: List[Dict] = []

    for length in sorted(decisions):
        decision = decisions[length]
        record = _rule_record(decision)

        if decision.decision == "safe_to_trim" and decision.trim_bases:
            override_trims[length] = {"trim_5p": 0, "trim_3p": decision.trim_bases}
            applied_rules.append(record)
        elif decision.decision == "safe_to_infer_but_not_trim" and decision.estimates:
            proposed_rules.append(record)

    return override_trims, applied_rules, proposed_rules, decisions


def plan_from_file(
    input_file: Path,
    format: str,
    infer_reads: int = 500_000,
    count_pattern: Optional[str] = None,
    known_adapters: Optional[List[str]] = None,
) -> Tuple[OverrideTrims, List[Dict], List[Dict], DecisionsByLength, Sketch]:
    """Build a Stage-0 sketch from a file and compute its apply plan.

    Inference runs on the bounded `infer_reads` subsample; the caller is
    responsible for applying the resulting override_trims on a full stream
    pass (RPFExtractor.extract_rpfs already streams the whole file). This
    is the single source of truth for "sketch then plan" -- callers should
    use this rather than re-deriving a Sketch and calling plan_from_sketch
    separately, so there's one code path to keep correct.
    """
    sketch = build_sketch(
        input_file,
        format=format,
        max_reads=infer_reads,
        count_pattern=count_pattern,
    )
    override_trims, applied_rules, proposed_rules, decisions = plan_from_sketch(
        sketch, known_adapters=known_adapters
    )
    return override_trims, applied_rules, proposed_rules, decisions, sketch


def _build_rule_record(values: Dict) -> Dict:
    """Project `values` onto exactly RULE_RECORD_KEYS -- the single schema
    every rule record shares. Raises KeyError if a caller forgets a key,
    rather than letting the two record shapes silently drift apart."""
    return {key: values[key] for key in RULE_RECORD_KEYS}


def _rule_record(decision: LengthClassDecision) -> Dict:
    return _build_rule_record(
        {
            "end": decision.end,
            "read_length_before": decision.length,
            "trim_bases": decision.trim_bases,
            "dominant_sequence": decision.dominant_base,
            "dominant_fraction": decision.dominant_fraction,
            "supporting_reads": decision.support,
            "reason": f"{decision.end}_terminal_bias",
            "estimators_agree": sorted(decision.estimates.keys()),
            "concordance": decision.concordance,
            "confidence": decision.confidence,
            "decision": decision.decision,
        }
    )


def override_rule_record(length: int, trim_3p: int, end: str = "3p") -> Dict:
    """Build a rule record for an explicit user-supplied override (no
    boundary estimation involved), with the same key set as
    `_rule_record` so trim_rules_applied has one consistent shape
    regardless of which code path produced it."""
    return _build_rule_record(
        {
            "end": end,
            "read_length_before": length,
            "trim_bases": trim_3p,
            "dominant_sequence": None,
            "dominant_fraction": None,
            "supporting_reads": None,
            "reason": "override_rule",
            "estimators_agree": [],
            "concordance": None,
            "confidence": None,
            "decision": "override_applied",
        }
    )
