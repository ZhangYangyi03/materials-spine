"""The campaign: a budgeted sequence of experiments drawn from a measured table.

The unit of cost is an *experiment*, not a model call. That is the whole reason
this package exists: in a simulation the objective can be called a million times,
in a lab it cannot, so every arm is scored in queries against a fixed budget and
nothing else.

Three rules make the comparison mean something:

  * a candidate is a row somebody already ran. A cell nobody measured is not a
    cheap observation, it is an experiment that has not happened, so `query()`
    refuses it. This is the mistake the module exists not to make.
  * replicates are folded: a cell run twice is one design measured twice, and
    the budget counts the design once.
  * the surrogate's hyperparameters are fitted once, on the initial design, and
    then held while the posterior is updated by conditioning. Re-tuning a prior
    after every plate is not what a campaign does, and letting the optimiser
    re-search it per round would smuggle a second budget into the comparison.
    `refit_every` exists for the study where you *want* that, and off is default.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import baselines as B
from .gp import MixedGP


def cells_from_table(space, codes, y):
    """Fold a table of rows into its distinct measured cells.

    Returns (cell_codes, cell_y, cell_n). A 1030-row concrete table defines a
    65536-cell product grid and a 4599-row plate defines 4608, and in both cases
    the grid is the product of the levels while the table is a sparse
    observation of it -- so anything that assumes a full factorial silently
    looks at the wrong thing.
    """
    codes = np.asarray(codes, dtype=int)
    y = np.asarray(y, dtype=float)
    if codes.ndim == 1:
        codes = codes.reshape(-1, 1)
    card = tuple(space.cardinality)
    if codes.shape[1] != len(card):
        raise ValueError(f"table has {codes.shape[1]} factor columns, space has {len(card)}")
    key = np.ravel_multi_index(tuple(codes.T), card)
    uniq, inv, counts = np.unique(key, return_inverse=True, return_counts=True)
    sums = np.zeros(len(uniq))
    np.add.at(sums, inv, y)
    cell_y = sums / counts
    first = np.unique(inv, return_index=True)[1]
    return codes[first], cell_y, counts


@dataclass
class Campaign:
    """One table, exposed as a metered screening campaign."""

    space: object
    codes: object
    y: object
    table: str = ""
    batch: int = 8
    seed: int = 0
    init: int = 8
    refit_every: int = 0          # 0 = fit once, then condition
    acq_kw: dict = field(default_factory=dict)
    objective_log: bool = False   # for targets with a floor of censored zeros

    def __post_init__(self):
        self.cell_codes, self.cell_y, self.replicates = cells_from_table(
            self.space, self.codes, self.y)
        self.n_cells = len(self.cell_y)
        self.rng = np.random.default_rng(self.seed)
        self.ceiling = float(self.cell_y.max())
        self.n_factors = len(self.space.cardinality)

    # ---- shape -----------------------------------------------------------
    @property
    def coverage(self):
        """Fraction of the product grid the table actually observes. Reported,
        never assumed: 99.8% for the amination plate, 0.4% for concrete, and the
        difficulty of a campaign is not the same thing in the two."""
        return self.n_cells / self.space.size

    def query(self, cell):
        if not (0 <= cell < self.n_cells):
            raise KeyError(f"cell {cell} out of range; the table has {self.n_cells} measured cells")
        return float(self.cell_y[cell])

    def label(self, cell):
        """A recommendation a chemist can read without a lookup table."""
        return self.space.decode(self.cell_codes[cell]) if hasattr(self.space, "decode") else None

    # ---- public numbers --------------------------------------------------
    def optimal(self, k=1):
        order = np.argsort(self.cell_y)[::-1][:k]
        return [{"cell": int(i), "value": float(self.cell_y[i]),
                 "recipe": self.label(int(i))} for i in order]

    # ---- the loop --------------------------------------------------------
    def run(self, budget, arm="ei", trace_every_batch=True):
        """Spend `budget` experiments. Returns (observed, trace).

        `observed` is a list of (cell, value) in the order they were spent --
        the only rows a method is allowed to have seen.
        """
        available = np.ones(self.n_cells, bool)
        observed = []

        for i in self._initial(available):
            available[i] = False
            observed.append((i, self.query(i)))
        trace = [self._snap(observed, arm)]

        gp = None
        spent_at_fit = 0
        while len(observed) < budget:
            nb = min(self.batch, budget - len(observed))
            if gp is None:
                gp = self.fit_gp(observed)          # the only marginal-likelihood fit
                spent_at_fit = len(observed)
            elif self.refit_every and len(observed) - spent_at_fit >= self.refit_every:
                gp = self.fit_gp(observed)          # explicit re-tuning study
                spent_at_fit = len(observed)
            else:
                self._fold(gp, observed, spent_at_fit)
                spent_at_fit = len(observed)
            picks = self._pick(arm, available, nb, observed, gp)
            picks = [p for p in picks if available[p]]
            if not picks:
                break
            for i in picks:
                available[i] = False
                observed.append((i, self.query(i)))
            trace.append(self._snap(observed, arm))
        return observed, trace

    def _fold(self, gp, observed, since):
        """Fold everything measured since the last fit into the held-prior posterior."""
        if since >= len(observed):
            return
        X = self.cell_codes[[i for i, _ in observed]]
        y = np.asarray([v for _, v in observed], dtype=float)
        gp.refit_posterior(X, y)

    def _initial(self, available):
        idx = np.flatnonzero(available)
        k = min(self.init, len(idx))
        return list(self.rng.choice(idx, size=k, replace=False))

    def _pick(self, arm, available, nb, observed, gp):
        if arm == "random":
            return B.random_arm(available, nb, self.rng)
        if arm == "fill":
            return B.fill_arm(self.cell_codes, available, nb)
        if arm == "ofat":
            return B.ofat_arm(self.space, self.cell_codes, available, observed, nb, self.rng)
        return self._propose(arm, available, nb, observed, gp)

    def fit_gp(self, observed):
        X = self.cell_codes[[o[0] for o in observed]]
        y = np.asarray([o[1] for o in observed], dtype=float)
        return MixedGP.for_space(self.space).fit(X, y)

    def _propose(self, arm, available, nb, observed, gp):
        from . import acq
        cand = np.flatnonzero(available)
        best = max(v for _, v in observed)
        if arm in ("qei", "qei_pess"):
            pessimism = 1.0 if arm == "qei" else 0.5
            picks = acq.q_ei_batch(gp, self.cell_codes[cand], nb, best,
                                   pessimism=pessimism,
                                   chunk=self.acq_kw.get("chunk", 192))
            return [int(cand[p]) for p in picks]
        mu, sd = gp.predict(self.cell_codes[cand])
        picks = acq.greedy_batch(arm, mu, sd, np.ones(len(cand), bool), nb, best,
                                 **self.acq_kw)
        return [int(cand[p]) for p in picks]

    def _snap(self, observed, arm):
        vals = [v for _, v in observed]
        best = max(vals) if vals else float("nan")
        return {"arm": arm, "spent": len(observed), "best": best,
                "best_frac": best / self.ceiling if self.ceiling else float("nan"),
                "remaining": int(self.n_cells - len(observed))}

    def regret_curve(self, trace):
        return [(t["spent"], t["best_frac"]) for t in trace]
