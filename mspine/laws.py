"""Engineering laws, tested so that they can fail.

Two exercises, both with the same shape: name a law that predicts something the
data does not contain, fit its one parameter, and then run the control that would
destroy the claim if the claim were an artefact.

    hydration   given ONE measured value of a concrete mix (its 28-day strength)
                and the ages it was tested at, which one-parameter law of time
                predicts its other ages? Two separate questions are kept apart
                on purpose: the FORM of the law (leave-one-age-out, same mix)
                and the CARRIER of its parameter (does composition predict it on
                mixes the fit never saw, against permuted targets?).

    weibull     weakest-link says a single Weibull modulus governs both the
                scatter within a gauge length and the scale shift across gauge
                lengths. On 255 single-filament carbon-fibre tests the two
                disagree -- and the shuffle null, which contains no length effect
                at all, disagrees by more. That negative is the result: the test
                as designed cannot support the claim either way, and what is
                reusable is the instrument plus its null.

Both are reported the same way: the number, the control, and the sentence the
evidence does not support.
"""
from __future__ import annotations

import math
import statistics as st

import numpy as np

from . import _law_core as L
from . import _weibull_core as W


# ---------------------------------------------------------------------------
# hydration law
# ---------------------------------------------------------------------------

LAW_FORMS = L.KINDS          # power / ceb (bounded) / aci


def hydration_form():
    """Which one-parameter law predicts an untested age of the SAME mix."""
    ms = L.mixes()
    scores, rows = L.form_scores(ms)
    return {"mixes": len(ms),
            "mixes_with_28d": sum(1 for m in ms if m["s28"]),
            "forms": {k: scores[k] for k in ("power", "ceb", "aci")},
            "ceb_beats_power": scores["ceb_beats_power"],
            "paired_n": scores["paired_n"],
            "long_ages": {str(a): {k: float(st.mean([abs(r[k] - r["truth"]) for r in rows
                                                    if r["age"] == a and r[k] is not None]))
                                   for k in LAW_FORMS}
                          for a in (56, 90, 180, 365)}}


def hydration_carrier(seed=42, nperm=300):
    """Is the law's parameter carried by composition, or by strength level?"""
    ms = L.mixes()
    out, n = L.carrier(ms, seed=seed, nperm=nperm)
    return {"n_mixes": n, "carriers": out}


def hydration_stability(n_splits=20):
    """Held-out R2 across splits: one split is an anecdote."""
    ms = L.mixes()
    agg = {}
    for seed in range(n_splits):
        c, _ = L.carrier(ms, seed=seed, nperm=0)
        for name, d in c.items():
            agg.setdefault(name, []).append(d["r2_heldout"])
    return {k: {"mean": float(st.mean(v)), "sd": float(st.pstdev(v)),
                "min": float(min(v)), "max": float(max(v)), "n_splits": len(v)}
            for k, v in agg.items()}


def hydration_controls(nrep=40, nrep_end=200, seed=1):
    """Two controls that would each kill the claim above on their own.

    1 time-order scramble: keep each mix's own strengths, give them to the wrong
      ages. A law that wins only because it can sit on any increasing data loses.
    2 composition shuffle: same fit, composition labels permuted between mixes.
      The spread of end-member numbers that comes out is what "no relation" looks
      like at this sample size.
    """
    ms = L.mixes()
    scr = L.form_control(ms, nrep=nrep)
    end_ctrl = L.end_member_control(ms, seed=seed, nrep=nrep_end)
    em = L.end_members(ms)
    real = {k: [c[k] for c in end_ctrl] for k in ("k_cement", "k_slag", "k_flyash")}
    return {
        "scrambled_mae": {k: {"mean": float(st.mean(v)), "sd": float(st.pstdev(v)),
                              "n_rep": len(v)} for k, v in scr.items()},
        "end_members": em,
        "end_member_null": {k: {"mean": float(np.mean(v)),
                                "p05": float(np.percentile(v, 5)),
                                "p95": float(np.percentile(v, 95))}
                            for k, v in real.items()},
    }


# ---------------------------------------------------------------------------
# weakest-link / Weibull
# ---------------------------------------------------------------------------

def weibull_test(field="sigma"):
    """The two estimates of the same modulus, and their ratio."""
    rows = W.load()
    r = W.scatter_vs_scaling(rows, field=field)
    return {"field": field, "n_specimens": len(rows),
            "per_length": {str(k): v for k, v in r["per_length"].items()},
            "slope": r["slope"], "slope_se": r["slope_se"],
            "m_scaling": r["m_scaling"], "m_scatter": r["m_scatter"],
            "ratio": r["ratio"], "required": 1.0}


def weibull_null(nrep=300, seed=7, field="sigma"):
    """The pool-and-redistribute null: no length effect by construction."""
    rows = W.load()
    null = W.shuffle_null(rows, field=field, nrep=nrep, seed=seed)
    base = W.scatter_vs_scaling(rows, field=field)
    arr = np.asarray(null, float)
    return {"n_rep": len(arr),
            "m_scaling_median": float(np.median(arr)),
            "m_scaling_p05": float(np.percentile(arr, 5)),
            "m_scaling_min": float(arr.min()),
            "share_at_or_above_observed": float((arr >= base["m_scaling"]).mean()),
            "share_plausible_range_3_25": float(((arr >= 3) & (arr <= 25)).mean()),
            "observed_ratio": base["ratio"]}


def weibull_batches():
    rows = W.load()
    pb = W.per_batch_moduli(rows)
    return {f"{k[0]}mm/{k[1]}": float(v) for k, v in pb.items()}
