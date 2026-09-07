"""Probabilistic read segmentation using Viterbi algorithm.

This module implements an HMM-based segmenter to identify read components
(UMI, Barcode, RPF, Adapter) by decoding observed signal statistics.
"""

import logging
import math
from dataclasses import dataclass
from typing import Dict, List, Optional

from .signals import SignalStats
from .types import SegmentInfo

logger = logging.getLogger(__name__)


# HMM States
STATE_START = 0
STATE_UMI = 1
STATE_BARCODE = 2
STATE_RPF = 3
STATE_ADAPTER = 4
STATE_END = 5

STATE_NAMES = {
    STATE_START: "START",
    STATE_UMI: "UMI",
    STATE_BARCODE: "BARCODE",
    STATE_RPF: "RPF",
    STATE_ADAPTER: "ADAPTER",
    STATE_END: "END",
}


@dataclass
class SegmenterConfig:
    """Configuration for HMM transition and emission probabilities."""

    # Priors
    # Informative priors (can be overridden by adaptive start logic):
    p_umi_start: float = 0.2  # Baseline prior for 5' UMI presence
    p_rpf_start: float = 0.7  # Baseline prior for direct RPF start
    p_adapter_start: float = 0.1  # Rare (short read or pre-trimmed tail)
    p_umi_len_mu: float = 6.0  # Expected UMI length (if present)
    p_rpf_len_mu: float = 30.0  # Expected RPF length
    min_umi_len: int = 4  # Do not allow UMI->RPF before this many bases
    min_rpf_len: int = 22  # (Advisory) minimal RPF before adapter

    # Emission profiles (Entropy means)
    mu_entropy_umi: float = 1.8  # High entropy
    mu_entropy_rpf: float = 1.2  # Medium/Variable entropy
    mu_entropy_adapter: float = 0.2  # Low entropy (consensus)

    # Composition markers
    adapter_consensus_threshold: float = (
        0.8  # Min freq to be considered 'consed' adapter base
    )


class ProbabilisticSegmenter:
    """HMM-based segmenter for read architecture discovery."""

    def __init__(self, config: SegmenterConfig = SegmenterConfig()):
        self.config = config

    def segment(self, stats: SignalStats) -> List[SegmentInfo]:
        """Infer segmentation from signal stats using Viterbi decoding.

        Args:
            stats: Signal statistics object

        Returns:
            List of detected segments
        """
        # We work primarily with 5' aligned signals for the forward pass
        obs_seq = self._prepare_observations(stats)

        if not obs_seq:
            return []

        return self._path_to_segments(self._viterbi_path(obs_seq))

    def _viterbi_path(self, obs_seq: List[Dict]) -> List[int]:
        """Decode the highest-scoring state at every observed position."""
        T = len(obs_seq)
        if not T:
            return []

        N = len(STATE_NAMES)
        viterbi = [[-float("inf")] * N for _ in range(T)]
        backpointer = [[0] * N for _ in range(T)]

        # Adaptive start prior: high early entropy favours a UMI, otherwise RPF.
        early_k = min(5, T)
        early_entropy = sum(o["entropy"] for o in obs_seq[:early_k]) / early_k
        umi_prior = 0.05 + 0.90 / (1.0 + math.exp(-(early_entropy - 1.3) / 0.3))
        rpf_prior = max(1e-6, 1.0 - umi_prior)
        viterbi[0][STATE_UMI] = math.log(umi_prior) + self._log_emission(
            STATE_UMI, obs_seq[0]
        )
        viterbi[0][STATE_RPF] = math.log(rpf_prior) + self._log_emission(
            STATE_RPF, obs_seq[0]
        )

        for t in range(1, T):
            for state in range(N):
                candidates = [
                    viterbi[t - 1][previous] + self._log_transition(previous, state, t)
                    for previous in range(N)
                ]
                best_previous = max(range(N), key=candidates.__getitem__)
                viterbi[t][state] = candidates[best_previous] + self._log_emission(
                    state, obs_seq[t]
                )
                backpointer[t][state] = best_previous

        path = [max(range(N), key=viterbi[-1].__getitem__)]
        for t in range(T - 1, 0, -1):
            path.append(backpointer[t][path[-1]])
        return path[::-1]

    # --- New: Viterbi + Forward-Backward posteriors for visualization ---
    def decode_with_posteriors(self, stats: SignalStats):
        """Return Viterbi path, posterior state probabilities, and segments.

        Produces per-position posteriors via forward-backward in log space,
        useful for plotting what the HMM is doing beyond the hard path.
        """
        obs_seq = self._prepare_observations(stats)
        if not obs_seq:
            return [], [], []

        T = len(obs_seq)
        N = 6

        # Emission log-probabilities
        emit = [[self._log_emission(s, obs_seq[t]) for s in range(N)] for t in range(T)]

        def logsumexp(vals):
            m = max(vals)
            if m == -float("inf"):
                return m
            return m + math.log(sum(math.exp(v - m) for v in vals))

        # Forward (time-inhomogeneous transitions allowed via t)
        fwd = [[-float("inf")] * N for _ in range(T)]
        # init: START not explicit; allow UMI/RPF starts
        fwd[0][STATE_UMI] = math.log(self.config.p_umi_start) + emit[0][STATE_UMI]
        fwd[0][STATE_RPF] = math.log(self.config.p_rpf_start) + emit[0][STATE_RPF]

        for t in range(1, T):
            for s in range(N):
                prevs = [
                    fwd[t - 1][ps] + self._log_transition(ps, s, t) for ps in range(N)
                ]
                fwd[t][s] = emit[t][s] + logsumexp(prevs)

        logZ = logsumexp(fwd[-1])

        # Backward
        bwd = [[-float("inf")] * N for _ in range(T)]
        for s in range(N):
            bwd[T - 1][s] = 0.0  # log(1)
        for t in range(T - 2, -1, -1):
            for s in range(N):
                nexts = [
                    self._log_transition(s, ns, t + 1)
                    + emit[t + 1][ns]
                    + bwd[t + 1][ns]
                    for ns in range(N)
                ]
                bwd[t][s] = logsumexp(nexts)

        # Posteriors gamma[t][s]
        post = []
        for t in range(T):
            gammas = [math.exp(fwd[t][s] + bwd[t][s] - logZ) for s in range(N)]
            # Normalize for numerical stability
            ssum = sum(gammas) or 1.0
            gammas = [g / ssum for g in gammas]
            post.append(gammas)

        # Keep the diagnostic path aligned with the posterior positions. The
        # previous implementation returned only the best terminal state.
        viterbi_path = self._viterbi_path(obs_seq)
        return viterbi_path, post, self._path_to_segments(viterbi_path)

    def _prepare_observations(self, stats: SignalStats) -> List[Dict]:
        """Convert stats to observation sequence."""
        obs = []
        for i in range(len(stats.entropy_5p)):
            o = {
                "entropy": stats.entropy_5p[i],
                # Add composition info?
                # e.g. "max_freq" to detect low entropy conservation
                "max_freq": (
                    max(stats.composition_5p[i].values())
                    if stats.composition_5p[i]
                    else 0.0
                ),
            }
            obs.append(o)
        return obs

    def _log_emission(self, state: int, obs: Dict) -> float:
        """Calculate log emission probability P(observation | state)."""
        entropy = obs["entropy"]
        max_freq = obs["max_freq"]

        # Gaussian approx for entropy
        if state == STATE_UMI:
            # Expect high entropy
            mu = self.config.mu_entropy_umi
            sigma = 0.5
        elif state == STATE_RPF:
            # Expect medium entropy
            mu = self.config.mu_entropy_rpf
            sigma = 0.6
        elif state == STATE_ADAPTER:
            # Expect low entropy (consensus)
            mu = self.config.mu_entropy_adapter
            sigma = 0.4
        else:
            return -10.0  # Unlikely observation for placeholder states

        # Log Gaussian PDF: -0.5 * ((x-mu)/sigma)^2 - log(sigma * sqrt(2pi))
        # Ignore constant terms for comparison
        log_prob = -0.5 * ((entropy - mu) / sigma) ** 2

        # Boost adapter probability if high single-nucleotide frequency (composition consensus)
        if (
            state == STATE_ADAPTER
            and max_freq > self.config.adapter_consensus_threshold
        ):
            log_prob += 2.0

        return log_prob

    def _log_transition(self, prev_s: int, curr_s: int, t: int) -> float:
        """Calculate log transition probability P(curr_s | prev_s)."""
        # Simple State Machine topology:
        # UMI -> UMI (stay)
        # UMI -> RPF
        # RPF -> RPF (stay)
        # RPF -> ADAPTER
        # ADAPTER -> ADAPTER (stay)

        if prev_s == STATE_UMI:
            # Enforce minimal UMI span before moving to RPF
            if curr_s == STATE_UMI:
                # Stronger self-loop before min_umi_len, weaker after
                stay = 0.95 if t < self.config.min_umi_len else 0.6
                return math.log(stay)
            elif curr_s == STATE_RPF:
                if t < self.config.min_umi_len:
                    return -float("inf")
                return math.log(0.4)  # encourage transition once min length met
            else:
                return -float("inf")

        elif prev_s == STATE_RPF:
            if curr_s == STATE_RPF:
                return math.log(0.95)  # RPFs are long
            elif curr_s == STATE_ADAPTER:
                # Soft encouragement to remain in RPF until emissions demand adapter
                return math.log(0.05)
            else:
                return -float("inf")

        elif prev_s == STATE_ADAPTER:
            # No transitions out of adapter; terminal region
            if curr_s == STATE_ADAPTER:
                return math.log(0.999)
            else:
                return -float("inf")

        elif prev_s == STATE_START:
            # Should rely on initialization
            return -float("inf")

        return -float("inf")

    def _path_to_segments(self, path: List[int]) -> List[SegmentInfo]:
        """Convert state path to SegmentInfo objects."""
        segments: List[SegmentInfo] = []
        if not path:
            return segments

        current_state = path[0]
        start_pos = 0

        for t, state in enumerate(path):
            if state != current_state:
                # End of segment
                if current_state in STATE_NAMES:  # Valid states
                    name = STATE_NAMES[current_state].lower()
                    if name != "end" and name != "start":
                        # TODO: Add confidence score per segment
                        seg = SegmentInfo(name, start_pos, t, 0.9)
                        segments.append(seg)

                current_state = state
                start_pos = t

        # Add last segment
        if current_state in STATE_NAMES:
            name = STATE_NAMES[current_state].lower()
            if name != "end" and name != "start":
                seg = SegmentInfo(name, start_pos, len(path), 0.9)
                segments.append(seg)

        return segments


def segment_reads(
    stats: SignalStats, config: Optional[SegmenterConfig] = None
) -> List[SegmentInfo]:
    """Infer read segments from signal statistics using the HMM segmenter."""
    return ProbabilisticSegmenter(config or SegmenterConfig()).segment(stats)


def decode_segments_with_posteriors(
    stats: SignalStats, config: Optional[SegmenterConfig] = None
):
    """Return HMM state posteriors and decoded segments for visualization."""
    return ProbabilisticSegmenter(config or SegmenterConfig()).decode_with_posteriors(
        stats
    )
