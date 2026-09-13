"""Observable structure scoring does not turn protocol unknowns into misses."""

from getRPF.core.structure.benchmark import prediction, score_reports


def _report(blocks, emit=True):
    return {"architecture": {"blocks": blocks}, "transform": {"emit": emit}}


def _block(kind, length=(0, 0), **extra):
    return {"type": kind, "length": list(length), **extra}


def test_protocol_unobservable_umi_is_excluded_but_observed_umi_is_scored():
    truth = {
        "runs": [
            {
                "run": "TRIMMED",
                "benchmark": {
                    "layout": "insert",
                    "adapter": False,
                    "tail": None,
                    "transform": "emit",
                    "umi": {"state": "not_observable"},
                },
            },
            {
                "run": "UMI",
                "benchmark": {
                    "layout": "random:5|insert|adapter",
                    "adapter": True,
                    "tail": None,
                    "transform": "emit",
                    "umi": {"state": "observed", "five_prime": [5]},
                },
            },
        ]
    }
    reports = {
        "TRIMMED": _report([_block("insert", (20, 40))]),
        "UMI": _report(
            [
                _block("random", (5, 5), keep_as_umi=True),
                _block("insert", (20, 40)),
                _block("adapter", (34, 34)),
            ]
        ),
    }
    summaries = {
        "TRIMMED": {"input_reads": 10, "accepted": 10},
        "UMI": {"input_reads": 10, "accepted": 10},
    }

    score = score_reports(truth, reports, summaries)

    assert score["metrics"]["umi"]["evaluated"] == 1
    assert score["metrics"]["umi"]["correct"] == 1
    assert score["metrics"]["umi"]["not_observable_runs"] == ["TRIMMED"]
    assert not [f for f in score["failures"] if f["type"] == "umi"]


def test_signature_ignores_an_interval_nta_and_preserves_umi_coordinates():
    report = _report(
        [
            _block("random", (3, 3), keep_as_umi=True),
            _block("insert", (20, 40)),
            _block("random", (4, 4), keep_as_umi=True),
            _block("nta", (0, 1)),
            _block("adapter", (17, 17)),
        ]
    )

    got = prediction(report)

    assert got["layout"] == "random:3|insert|random:4|adapter"
    assert got["five_prime_umi"] == (3,)
    assert got["three_prime_umi"] == (4,)
