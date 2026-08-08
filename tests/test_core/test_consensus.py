from getRPF.core.processors.consensus import TrimDecider, decide_trim_consensus


def test_trim_consensus_function_preserves_legacy_wrapper_behavior():
    architecture = {
        "architecture_match": "known_protocol",
        "trim_recommendations": {
            "recommended_5prime_trim": 4,
            "recommended_3prime_trim": 0,
            "three_prime_adapter": "AGATC",
        },
    }
    alignment = {
        "trim_recommendations": {
            "recommended_5prime_trim": 4,
            "recommended_3prime_trim": 0,
            "consensus_level": 0.9,
            "global_pattern_detected": True,
        }
    }

    functional_result = decide_trim_consensus(architecture, alignment)
    wrapper_result = TrimDecider().decide(architecture, alignment)

    assert functional_result == wrapper_result
    assert functional_result.trim_5p == 4
    assert functional_result.adapter_sequence == "AGATC"
