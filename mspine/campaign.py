"""What a campaign is allowed to claim, and the control that has to agree.

The campaign's own best-so-far curve cannot answer the question a screening tool
is bought for. That curve is measured on rows the method chose and then looked
at; a method that picks nothing would also look fine on it if the table were
kind. The testable claim is about cells the method never saw:

    fit on exactly what the campaign spent
    rank the unmeasured cells by predicted value, take the top k
    compare their measured mean with the measured mean of everything unmeasured

That gap is the screening gain, it is computed from rows outside the training
set, and it is zero when the method has learned nothing. Which is what the
negative control is for: repeat the whole pipeline with the responses permuted
across cells. The marginal distribution of the target is preserved and the
recipe-response relation is destroyed, so a gain that survives permutation is a
gain coming from the features' correlation with themselves, and every number in
the report is void.

Two more quantities are reported because the ranking above can be right for the
wrong reason: `top10_recovery` (share of the table's true best ten cells that
the model put in its own top ten) and `heldout_rank_corr` (Spearman over the
unmeasured cells). A surrogate can order cells coarsely and still never point at
the best ten, and that difference is the finding, not a footnote.
"""
from __future__ import annotations

import numpy as np

from .loop import Campaign
from .gp import MixedGP


def _permute(y, rng):
    return np.asarray(y, float)[rng.permutation(len(y))]


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / den) if den else float("nan")


def assess(campaign, observed, k_frac=0.1, k_top=10):
    """The held-out verdict for one finished campaign."""
    seen = np.zeros(campaign.n_cells, bool)
    seen[[i for i, _ in observed]] = True
    unseen = np.flatnonzero(~seen)
    if len(unseen) < 20:
        return None
    X = campaign.cell_codes[[i for i, _ in observed]]
    y = np.asarray([v for _, v in observed], float)
    gp = MixedGP.for_space(campaign.space).fit(X, y)
    mu, sd = gp.predict(campaign.cell_codes[unseen])
    truth = campaign.cell_y[unseen]

    k = max(1, int(len(unseen) * k_frac))
    top = np.argsort(mu)[::-1][:k]
    rest_mean = float(truth.mean())
    out = {
        "n_unseen": int(len(unseen)),
        "k": int(k),
        "top_mean": float(truth[top].mean()),
        "unseen_mean": rest_mean,
        "gain": float(truth[top].mean() - rest_mean),
        "gain_frac_of_range": float((truth[top].mean() - rest_mean) /
                                    (campaign.ceiling - campaign.cell_y.min() or 1.0)),
        "rank_corr": _spearman(mu, truth),
        "best_seen": float(max(v for _, v in observed)),
        "best_seen_frac": float(max(v for _, v in observed) / campaign.ceiling),
    }
    # how much of the table's own best-k the model put in its own best-k
    true_top = set(np.argsort(truth)[::-1][:k_top].tolist())
    pred_top = set(np.argsort(mu)[::-1][:k_top].tolist())
    out["top%d_recovery" % k_top] = len(true_top & pred_top) / float(k_top)
    return out


def run_one(space, codes, y, name="", arm="ei", budget=64, seed=0, batch=8,
            k_frac=0.1, control=False, ctrl_seed=None, refit=0):
    """One campaign plus its held-out assessment, optionally with destroyed targets.

    `refit` is the fitting policy, and it is part of the arm rather than a detail:
    0 fits the hyperparameters once on the initial design and then only conditions
    the posterior, while `refit=8` re-runs the marginal-likelihood fit after every
    batch. That single choice moves the Buchwald result by more than the gap
    between the arms, so a record that does not carry it is not comparable with
    one that does.
    """
    yy = np.asarray(y, float)
    if control:
        yy = _permute(yy, np.random.default_rng(ctrl_seed if ctrl_seed is not None else seed + 991))
    c = Campaign(space=space, codes=codes, y=yy, table=name, batch=batch, seed=seed,
                 refit_every=refit)
    observed, trace = c.run(budget, arm=arm)
    a = assess(c, observed, k_frac=k_frac)
    return {"table": name, "arm": arm, "budget": budget, "seed": seed,
            "control": bool(control), "refit": int(refit), "trace": trace, "heldout": a,
            "n_cells": c.n_cells, "coverage": c.coverage}


ARMS = ("random", "ofat", "fill", "ei", "logei", "ucb", "pi", "qei")


def compare(space, codes, y, name="", arms=("random", "ofat", "ei", "logei", "ucb"),
            budget=64, seeds=(0, 1, 2), batch=8, k_frac=0.1, controls=True,
            refits=(0,)):
    """Every arm, every seed, with the permuted control for each arm.

    Returned as a flat list of records: one row per (arm, seed, control) so the
    report can aggregate without re-deriving anything, and the ledger can store
    exactly what was computed.
    """
    out = []
    for refit in refits:
        for arm in arms:
            for seed in seeds:
                out.append(run_one(space, codes, y, name, arm, budget, seed, batch,
                                   k_frac, refit=refit))
                if controls:
                    out.append(run_one(space, codes, y, name, arm, budget, seed, batch,
                                       k_frac, control=True, refit=refit))
    return out


def aggregate(records):
    """Mean +- sd per arm, real and control, on the two numbers that matter."""
    keys = {}
    for r in records:
        if r["heldout"] is None:
            continue
        k = (r["arm"], r["control"], r.get("refit", 0))
        keys.setdefault(k, []).append(r["heldout"])
    rows = []
    for (arm, ctrl, refit), hs in sorted(keys.items()):
        g = np.asarray([h["gain"] for h in hs], float)
        rc = np.asarray([h["rank_corr"] for h in hs], float)
        rec = np.asarray([h["top10_recovery"] for h in hs], float)
        bf = np.asarray([h["best_seen_frac"] for h in hs], float)
        rows.append({
            "arm": arm, "control": ctrl, "refit": refit, "n": len(hs),
            "gain_mean": float(g.mean()), "gain_sd": float(g.std(ddof=1)) if len(g) > 1 else 0.0,
            "rank_corr_mean": float(np.nanmean(rc)),
            "top10_recovery_mean": float(rec.mean()),
            "best_seen_frac_mean": float(bf.mean()),
        })
    return rows
