"""The invariants a campaign must not be able to break.

Each of these is a mistake that was made once in this project's history, or one
that a plausible-looking implementation makes silently:

  the metered API refuses cells nobody measured -- otherwise "budget" is a claim
  about the optimiser and not about experiments;

  the initial design is identical across arms -- otherwise a gain can come from a
  luckier start rather than from the surrogate;

  replicates are folded into one cell -- otherwise a table where the good recipe
  was run twice gets a budget advantage over one where it was run once;

  permuting the responses destroys the held-out gain -- otherwise the gain comes
  from the features' correlation with themselves.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mspine import campaign, data
from mspine.loop import Campaign, cells_from_table


def test_query_refuses_unmeasured_cell():
    space, codes, y, meta = data.concrete()
    c = Campaign(space=space, codes=codes, y=y, seed=0)
    with pytest.raises(KeyError):
        c.query(c.n_cells + 5)


def test_replicates_fold_to_one_cell():
    space, codes, y, meta = data.buchwald()
    cc, cy, cn = cells_from_table(space, codes, y)
    assert len(cy) == len(set(map(tuple, codes.tolist())))
    assert cn.sum() == len(y)
    assert len(cy) <= len(y)


def test_same_initial_design_across_arms():
    space, codes, y, meta = data.buchwald()
    a = Campaign(space=space, codes=codes, y=y, seed=3, init=8)
    b = Campaign(space=space, codes=codes, y=y, seed=3, init=8)
    assert a._initial(np.ones(a.n_cells, bool)) == b._initial(np.ones(b.n_cells, bool))


def test_permuted_responses_kill_the_gain():
    """The negative control has to actually destroy the signal."""
    space, codes, y, meta = data.concrete()
    real = campaign.run_one(space, codes, y, "concrete", arm="ei", budget=40, seed=0)
    ctrl = campaign.run_one(space, codes, y, "concrete", arm="ei", budget=40, seed=0,
                            control=True, ctrl_seed=0)
    assert real["heldout"]["gain"] > ctrl["heldout"]["gain"]


def test_budget_is_respected():
    space, codes, y, meta = data.buchwald()
    for arm in ("random", "ofat", "ei"):
        c = Campaign(space=space, codes=codes, y=y, seed=1, batch=8)
        obs, trace = c.run(32, arm=arm)
        assert len(obs) == 32
        assert trace[-1]["spent"] == 32
        assert len({i for i, _ in obs}) == 32      # no cell spent twice


def test_fitting_policy_is_part_of_the_record():
    """Two policies are two arms; neither may overwrite the other in the aggregate."""
    space, codes, y, meta = data.concrete()
    recs = []
    for refit in (0, 8):
        r = campaign.run_one(space, codes, y, "concrete", arm="ei", budget=24,
                             seed=0, refit=refit)
        assert r["refit"] == refit
        recs.append(r)
    rows = campaign.aggregate(recs)
    assert len(rows) == 2
    assert {row["refit"] for row in rows} == {0, 8}
