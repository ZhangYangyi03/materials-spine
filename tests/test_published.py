"""The numbers the documents quote, re-derived -- so a doc cannot drift from the code.

The old repositories each had one of these and it caught real drift. The rule is
narrow on purpose: a test may re-derive a number that appears verbatim in README.md
or docs/paper.md, and nothing else. If a refactor changes the result, the failure
should name the document that is now lying rather than the function that moved.

Every one of these runs whole campaigns on the 4599-cell Buchwald table, so they
are marked `slow` and opt-in (`pytest -m slow`), and CI runs both sets. They are
also the reason the campaigns are cached per test session: this host's allocator
dies when the same Buchwald family is re-materialised a dozen times in one
process, and a test suite that only passes on a big machine is a test suite that
gets annotated away.

The first version of this file asserted on best-so-far/ceiling, because that is
what the documents quoted. When the documents moved to held-out gain -- the
metric that asks about rows the campaign did not choose -- the assertions moved
with them, and both metrics turn out to show the same thing: one seed carries
the refit gap. That agreement is the reason the claim survived a change of
metric, and it is worth stating so that neither metric looks like the safe one.
"""
import gc
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mspine import campaign, data

pytestmark = pytest.mark.slow

_CACHE = {}


@pytest.fixture(scope="module")
def buchwald():
    space, codes, y, meta = data.buchwald()
    _CACHE.clear()
    return space, codes, y


def _gain(space, codes, y, arm, seed, refit, budget=64):
    """Held-out gain: the number docs/paper.md quotes. It calls the same function
    the benchmark calls, so a refactor that moves the metric moves this too."""
    key = ("gain", arm, seed, refit, budget)
    if key not in _CACHE:
        _CACHE[key] = campaign.heldout_gain(space, codes, y, "buchwald", arm=arm,
                                            budget=budget, seed=seed, refit=refit)
        gc.collect()
    return _CACHE[key]


@pytest.mark.parametrize("refit,mean,sd", [(0, 22.20, 7.10), (8, 26.06, 8.94)])
def test_buchwald_qei_gain_at_both_fitting_policies(buchwald, refit, mean, sd):
    """docs/paper.md section 2 quotes these for qei, 5 seeds, budget 64, batch 8."""
    space, codes, y = buchwald
    vals = [_gain(space, codes, y, "qei", s, refit) for s in range(5)]
    assert np.mean(vals) == pytest.approx(mean, abs=0.05)
    assert np.std(vals, ddof=1) == pytest.approx(sd, abs=0.05)


def test_the_refit_gap_is_one_seed_and_its_spread_is_wider(buchwald):
    """docs/paper.md section 3, and the sentence in section 2 that says the effect
    is smaller than its own spread.

    The claim is narrow on purpose: one seed out of five recovers from a stall by
    a large margin, three of the other four move the wrong way, and the sd of the
    per-seed difference is larger than its mean. The first version of this file
    asserted the same thing on best-so-far/ceiling and it was true there too --
    two metrics, one seed. If either test starts failing, the document is what is
    wrong, not the test.
    """
    space, codes, y = buchwald
    a = np.array([_gain(space, codes, y, "qei", s, 0) for s in range(5)])
    b = np.array([_gain(space, codes, y, "qei", s, 8) for s in range(5)])
    d = b - a
    i = int(np.argmax(d))
    assert d[i] > 25.0, "the stall-and-recover seed is gone: %r" % (d.tolist(),)
    assert int((d < 0).sum()) >= 2, "re-tuning now helps most seeds: %r" % (d.tolist(),)
    assert d.std(ddof=1) > abs(d.mean()), "the effect is no longer smaller than its spread"
    # and the mean is what the section quotes
    assert b.mean() - a.mean() == pytest.approx(3.86, abs=0.1)


def test_no_bayesian_arm_beats_ofat_on_buchwald(buchwald):
    """The README's uncomfortable claim, at the metric it is stated in."""
    space, codes, y = buchwald
    recs = campaign.compare(space, codes, y, "buchwald",
                            arms=("ofat", "ei", "qei"),
                            seeds=(0, 1), controls=False, refits=(0, 8))
    gc.collect()
    rows = campaign.aggregate(recs)
    ofat = [r for r in rows if r["arm"] == "ofat"][0]["gain_mean"]
    for r in rows:
        if r["arm"] == "ofat":
            continue
        assert r["gain_mean"] < ofat, "%s (refit=%d) now beats ofat (%+0.2f vs %+0.2f)" % (
            r["arm"], r["refit"], r["gain_mean"], ofat)


def test_the_control_collapses():
    """Every table, one arm: the permuted run must be near zero, not merely smaller."""
    for name in ("buchwald", "concrete"):
        space, codes, y, meta = data.load(name)
        rows = campaign.aggregate(campaign.compare(
            space, codes, y, name, arms=("ei",), seeds=(0,), controls=True,
            budget=40, refits=(0,)))
        gc.collect()
        real = [r for r in rows if not r["control"]][0]
        ctrl = [r for r in rows if r["control"]][0]
        assert real["gain_mean"] > 3.0
        assert abs(ctrl["gain_mean"]) < 3.0, "%s control leaked: %+0.2f" % (
            name, ctrl["gain_mean"])
