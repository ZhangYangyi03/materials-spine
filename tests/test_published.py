"""The numbers the documents quote, re-derived -- so a doc cannot drift from the code.

The old repositories each had one of these and it caught real drift. The rule is
narrow on purpose: a test may re-derive a number that appears verbatim in README.md
or docs/paper.md, and nothing else. If a refactor changes the result, the failure
should name the document that is now lying rather than the function that moved.

Every one of these runs whole campaigns on the 4599-cell Buchwald table, so they
are marked `slow` and opt-in (`pytest -m slow`). They are also the reason the
campaigns are run once and shared: this host's allocator dies when the same
Buchwald family is re-materialised a dozen times in one process, and a test suite
that only passes on a big machine is a test suite that gets annotated away.
"""
import gc
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mspine import campaign, data
from mspine.loop import Campaign

pytestmark = pytest.mark.slow

_CACHE = {}


def _best_frac(space, codes, y, arm, seed, refit, budget=64):
    key = (arm, seed, refit, budget)
    if key not in _CACHE:
        c = Campaign(space=space, codes=codes, y=y, batch=8, seed=seed, init=8,
                     refit_every=refit)
        obs, tr = c.run(budget, arm=arm)
        _CACHE[key] = tr[-1]["best"] / c.ceiling
        del c
        gc.collect()
    return _CACHE[key]


@pytest.fixture(scope="module")
def buchwald():
    space, codes, y, meta = data.buchwald()
    _CACHE.clear()
    return space, codes, y


def test_buchwald_qei_at_refit8_is_the_published_0_9896(buchwald):
    """docs/paper.md section 2 quotes 0.9896 +- 0.010 for qei at refit=8."""
    space, codes, y = buchwald
    vals = [_best_frac(space, codes, y, "qei", s, 8) for s in range(5)]
    assert np.mean(vals) == pytest.approx(0.9896, abs=2e-3)
    assert np.std(vals) == pytest.approx(0.010, abs=3e-3)


def test_one_seed_carries_the_refit_gap(buchwald):
    """docs/paper.md section 3: three seeds untouched, one seed is the gap.

    The first version asserted four untouched seeds and failed on seed 4, which
    moves by -0.0038. That disagreement between a document and a run is exactly
    what this file exists to surface: the claim is "one seed carries the gap",
    not "four seeds are bit-identical".
    """
    space, codes, y = buchwald
    same, moved = 0, []
    for s in range(5):
        a = _best_frac(space, codes, y, "qei", s, 0)
        b = _best_frac(space, codes, y, "qei", s, 8)
        if abs(a - b) < 1e-9:
            same += 1
        else:
            moved.append((s, a, b))
    assert same == 3, "unchanged seeds: %d, moved: %r" % (same, moved)
    big = [m for m in moved if m[2] - m[1] > 0.1]
    small = [m for m in moved if abs(m[2] - m[1]) <= 0.1]
    assert len(big) == 1, "the gap is no longer one seed: %r" % (moved,)
    assert len(small) == 1 and small[0][1] > small[0][2], moved


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
