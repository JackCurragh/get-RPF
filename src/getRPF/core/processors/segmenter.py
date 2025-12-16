"""Probabilistic read segmentation using Viterbi algorithm.

This module implements an HMM-based segmenter to identify read components
(UMI, Barcode, RPF, Adapter) by decoding observed signal statistics.
"""

import math
import logging
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
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
    STATE_END: "END"
}

@dataclass
class SegmenterConfig:
    """Configuration for HMM transition and emission probabilities."""
    # Priors
    p_umi_start: float = 0.3    # Probability read starts with UMI
    p_rpf_start: float = 0.6    # Probability read starts with RPF (no UMI)
    p_umi_len_mu: float = 6.0   # Expected UMI length
    p_rpf_len_mu: float = 30.0  # Expected RPF length
    
    # Emission profiles (Entropy means)
    mu_entropy_umi: float = 1.8   # High entropy
    mu_entropy_rpf: float = 1.2   # Medium/Variable entropy
    mu_entropy_adapter: float = 0.2  # Low entropy (consensus)
    
    # Composition markers
    adapter_consensus_threshold: float = 0.8 # Min freq to be considered 'consed' adapter base


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
            
        # Viterbi Algorithm
        # T = sequence length
        # N = number of states
        T = len(obs_seq)
        N = 6 
        
        # log_prob matrix: T x N
        # path matrix: T x N (stores backpointers)
        viterbi = [[-float('inf')] * N for _ in range(T)]
        backpointer = [[0] * N for _ in range(T)]
        
        # Initialization (t=0)
        # We can start in UMI, RPF, or rarely ADAPTER (if very short read)
        viterbi[0][STATE_UMI] = math.log(self.config.p_umi_start) + self._log_emission(STATE_UMI, obs_seq[0])
        viterbi[0][STATE_RPF] = math.log(self.config.p_rpf_start) + self._log_emission(STATE_RPF, obs_seq[0])
        # Allow barcode start? Let's treat Barcode same as UMI for now structure-wise or just after UMI.
        
        # Iteration
        for t in range(1, T):
            for s in range(N):
                # Find best transition to state s
                max_tr_prob = -float('inf')
                best_prev_s = -1
                
                for prev_s in range(N):
                    # Transition probability prev_s -> s
                    tr_prob = self._log_transition(prev_s, s, t)
                    
                    prob = viterbi[t-1][prev_s] + tr_prob
                    if prob > max_tr_prob:
                        max_tr_prob = prob
                        best_prev_s = prev_s
                
                # Multiply by emission probability
                emission_prob = self._log_emission(s, obs_seq[t])
                viterbi[t][s] = max_tr_prob + emission_prob
                backpointer[t][s] = best_prev_s

        # Termination
        # Best final state (likely ADAPTER or RPF or END)
        max_final_prob = -float('inf')
        best_last_state = -1
        
        for s in range(N):
            if viterbi[T-1][s] > max_final_prob:
                max_final_prob = viterbi[T-1][s]
                best_last_state = s
                
        # Backtracking
        best_path = [best_last_state]
        for t in range(T-1, 0, -1):
            prev_s = backpointer[t][best_path[-1]]
            best_path.append(prev_s)
            
        best_path = best_path[::-1] # Reverse
        
        # Convert path to SegmentInfo
        return self._path_to_segments(best_path)

    def _prepare_observations(self, stats: SignalStats) -> List[Dict]:
        """Convert stats to observation sequence."""
        obs = []
        for i in range(len(stats.entropy_5p)):
            o = {
                "entropy": stats.entropy_5p[i],
                # Add composition info? 
                # e.g. "max_freq" to detect low entropy conservation
                "max_freq": max(stats.composition_5p[i].values()) if stats.composition_5p[i] else 0.0
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
            return -10.0 # Unlikely observation for placeholder states
            
        # Log Gaussian PDF: -0.5 * ((x-mu)/sigma)^2 - log(sigma * sqrt(2pi))
        # Ignore constant terms for comparison
        log_prob = -0.5 * ((entropy - mu) / sigma) ** 2
        
        # Boost adapter probability if high single-nucleotide frequency (composition consensus)
        if state == STATE_ADAPTER and max_freq > self.config.adapter_consensus_threshold:
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
            if curr_s == STATE_UMI:
                return math.log(0.8) # Tend to stay
            elif curr_s == STATE_RPF:
                return math.log(0.2)
            else:
                return -float('inf')
                
        elif prev_s == STATE_RPF:
            if curr_s == STATE_RPF:
                return math.log(0.95) # RPFs are long
            elif curr_s == STATE_ADAPTER:
                return math.log(0.05)
            else:
                return -float('inf')
                
        elif prev_s == STATE_ADAPTER:
            if curr_s == STATE_ADAPTER:
                return math.log(0.99) # Stay in adapter till end
            else:
                return -float('inf')
                
        elif prev_s == STATE_START:
            # Should rely on initialization
            return -float('inf')
            
        return -float('inf')

    def _path_to_segments(self, path: List[int]) -> List[SegmentInfo]:
        """Convert state path to SegmentInfo objects."""
        segments = []
        if not path:
            return segments
            
        current_state = path[0]
        start_pos = 0
        
        for t, state in enumerate(path):
            if state != current_state:
                # End of segment
                if current_state in STATE_NAMES: # Valid states
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
