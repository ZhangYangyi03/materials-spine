"""A reaction network as arrays, and the four stages run over it.

Ported from kinetics-spine (Apache-2.0, https://github.com/ZhangYangyi03/kinetics-spine,
which was itself the four-stage skeleton of eda-spine pointed at a chemical
mechanism): search with a gate, re-ask a proved property of the candidate, break
one thing and require the gate to say no, measure what moved. The upstream
licence travels in LICENSE-kinetics-spine. What is added here is that every
verdict is handed to the ledger, and that the gate's own blindness is measured
across observation windows instead of assumed away.
"""
from __future__ import annotations

import importlib.util
import os

import numpy as np

R_GAS = 8.314462618  # J/mol/K


def load_mechanism(path):
    spec = importlib.util.spec_from_file_location("mech", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Mechanism:
    """Species, reactions, elements, stoichiometry -- and nothing else.

    A reduction is a Mechanism too, so the property stage does not need to know
    which one it was handed.
    """

    def __init__(self, species, reactions, elements):
        self.species = list(species)
        self.index = {s: i for i, s in enumerate(self.species)}
        self.elements = list(elements)
        self.reactions = list(reactions)   # (name, {reactant: n}, {product: n}, (A, b, Ea))
        self.n = len(self.species)
        self.m = len(self.reactions)
        self.nu = np.zeros((self.m, self.n))
        for j, (_name, react, prod, _a) in enumerate(self.reactions):
            for s, c in react.items():
                self.nu[j, self.index[s]] -= c
            for s, c in prod.items():
                self.nu[j, self.index[s]] += c
        self.comp = np.zeros((self.n, len(self.elements)))
        for i, s in enumerate(self.species):
            for e, c in species[s].items():
                self.comp[i, self.elements.index(e)] = c
        self.arrh = np.array([r[3] for r in self.reactions], dtype=float)
        # A multiplicative factor applied at evaluation time, kept OUT of arrh on
        # purpose: a reduction is rebuilt from the raw reaction tuples, and when
        # the factor lived inside arrh the reduced mechanism silently lost it and
        # overflowed while the full one did not (measured: nan for every candidate).
        self.rate_scale = 1.0

    # ---- the gate's ground truth -----------------------------------------
    def balance_residual(self):
        """Max element imbalance per reaction.

        Zero is the only acceptable value: this is an identity, not a tolerance.
        """
        return float(np.max(np.abs(self.nu @ self.comp)))

    # ---- integration ------------------------------------------------------
    def rates(self, T):
        A, b, Ea = self.arrh[:, 0], self.arrh[:, 1], self.arrh[:, 2]
        return self.rate_scale * A * (T ** b) * np.exp(-Ea / (R_GAS * T))

    def normalize_rates(self, T, target=1.0e3):
        """Rescale the prefactors so the fastest reaction is `target` per unit time.

        The prefactors are GRI values in cm3/mol/s while y is mole fractions;
        together they make k*t_end of order 1e15 and RK4 overflows in a few steps.
        Dividing A by one constant fixes the time scale without touching the
        stoichiometry, so every invariant is exactly what it was.
        """
        peak = float(np.max(self.rates(T))) / self.rate_scale
        if peak > 0:
            self.rate_scale = target / peak
        return self.rate_scale

    def rhs(self, y, T):
        q = self.rates(T).copy()
        for j, (_name, react, _prod, _a) in enumerate(self.reactions):
            for s, c in react.items():
                q[j] *= max(y[self.index[s]], 0.0) ** c
        return self.nu.T @ q

    def integrate(self, y0, T, t_end, steps):
        y = np.asarray(y0, dtype=float)
        h = t_end / steps
        traj = np.empty((steps + 1, self.n))
        traj[0] = y
        for k in range(steps):
            k1 = self.rhs(y, T)
            k2 = self.rhs(y + 0.5 * h * k1, T)
            k3 = self.rhs(y + 0.5 * h * k2, T)
            k4 = self.rhs(y + h * k3, T)
            y = np.maximum(y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4), 0.0)
            traj[k + 1] = y
        return traj

    def integrate_stiff(self, y0, T, t_end, steps):
        """An implicit adaptive solver over the same equations.

        Without a second family, a step-size bug looks exactly like a physical
        result and the gate agrees with it.
        """
        from scipy.integrate import solve_ivp
        t_eval = np.linspace(0.0, t_end, steps + 1)
        # LSODA passes (t, y) or (y) depending on the signature; wrap explicitly
        # rather than rely on `args`, which does not match this rhs's shape.
        f = lambda _t, y: self.rhs(y, T)
        sol = solve_ivp(f, (0.0, t_end), np.asarray(y0, float), method="LSODA",
                        t_eval=t_eval, rtol=1e-9, atol=1e-16)
        return sol.y.T

    # ---- edits ------------------------------------------------------------
    def reduce_dropping(self, names):
        keep = [r for r in self.reactions if r[0] not in names]
        red = rebuild(self.species, self.comp, self.elements, keep)
        red.rate_scale = self.rate_scale
        return red

    def with_scaled_rate(self, rname, factor):
        reac = [(n, r, p, (a[0] * factor, a[1], a[2]) if n == rname else a)
                for (n, r, p, a) in self.reactions]
        out = rebuild(self.species, self.comp, self.elements, reac)
        out.rate_scale = self.rate_scale
        return out

    def fluxes(self, y0, T, t_end, steps):
        """Integrated flux per reaction: which reactions a change to is visible in.

        Perturbing a dead reaction comes back PROVED, and that is gate blindness,
        not a property. Measuring the flux first is what makes the choice of
        perturb an experiment instead of a guess.
        """
        traj = self.integrate(y0, T, t_end, steps)
        h = t_end / steps
        tot = np.zeros(self.m)
        for k in range(steps):
            q = self.rates(T).copy()
            for j, (_n, react, _p, _a) in enumerate(self.reactions):
                for sp, c in react.items():
                    q[j] *= max(traj[k][self.index[sp]], 0.0) ** c
            tot += q
        return tot * h


def rebuild(species, comp, elements, reactions):
    m = object.__new__(Mechanism)
    m.species = list(species)
    m.index = {s: i for i, s in enumerate(m.species)}
    m.elements = list(elements)
    m.reactions = list(reactions)
    m.n = len(m.species)
    m.m = len(m.reactions)
    m.nu = np.zeros((m.m, m.n))
    for j, (_n, react, prod, _a) in enumerate(m.reactions):
        for s, c in react.items():
            m.nu[j, m.index[s]] -= c
        for s, c in prod.items():
            m.nu[j, m.index[s]] += c
    m.comp = comp
    m.arrh = np.array([r[3] for r in m.reactions], dtype=float)
    m.rate_scale = 1.0
    return m


# ---------------------------------------------------------------------------
# stage 1: search, with a gate that can only reject
# ---------------------------------------------------------------------------

def initial_vector(mech, initial):
    y = np.zeros(mech.n)
    for s, x in initial.items():
        if s in mech.index:
            y[mech.index[s]] = x
    return y


def integration_floor(path, t_end=None, steps=None, refine=16):
    """The gate's own noise floor at one observation window.

    RK4 at `steps` steps is not exact, so "the trajectories agree to 1e-6" is a
    statement that cannot be checked more finely than the integrator's own error
    at that resolution. Measured on this mechanism: at the configured 400 steps the
    RK4 error against a 16x refined run is 9e-4 -- three orders of magnitude
    coarser than the gate's tolerance. Two consequences, both of them the reason
    this function exists:

      a reduction is accepted only when the dropped steps carry less flux than
      the integrator can resolve, and saying "equivalent to 1e-6" without that
      floor next to it is a claim the numerics do not support;

      the second integrator can only be used as a bug check at a *resolved* step
      count, because at 400 steps LSODA and RK4 disagree by their discretisation
      error, not by a defect.
    """
    mod = load_mechanism(path)
    full = Mechanism(mod.SPECIES, mod.REACTIONS, mod.ELEMENTS)
    cond = dict(mod.CONDITIONS)
    if t_end is not None:
        cond["t_end"] = float(t_end)
    if steps is not None:
        cond["steps"] = int(steps)
    full.normalize_rates(cond["T"])
    y0 = initial_vector(full, mod.INITIAL)
    coarse = full.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    fine = full.integrate(y0, cond["T"], cond["t_end"], cond["steps"] * refine)
    sub = fine[::refine]
    return {"t_end": cond["t_end"], "steps": cond["steps"], "refine": refine,
            "rk4_self_error": float(np.max(np.abs(coarse - sub)) /
                                    max(float(np.max(np.abs(fine))), 1e-30))}


def equivalence_gate(full, red, y0, cond, tol=1e-6):
    """Two mandatory questions: does the reduction conserve every element, and does
    it reproduce the full trajectory under the same conditions?"""
    T, t_end, steps = cond["T"], cond["t_end"], cond["steps"]
    tf = full.integrate(y0, T, t_end, steps)
    tr = red.integrate(y0, T, t_end, steps)
    scale = max(float(np.max(np.abs(tf))), 1e-30)
    err = float(np.max(np.abs(tf - tr))) / scale
    resid = red.balance_residual()
    return {"balance": resid, "trajectory": err,
            "ok": bool(resid < 1e-12 and err < tol)}


def search(full, y0, cond, candidates, tol=1e-6):
    """Try each reduction; keep only what the gate passes. The gate can only
    reject -- without it "fewer reactions" is not a result, it is a smaller file."""
    rows, kept = [], []
    for label, names in candidates:
        red = full.reduce_dropping(names)
        g = equivalence_gate(full, red, y0, cond, tol=tol)
        rows.append({"label": label, "dropped": list(names), "reactions": red.m, **g})
        if g["ok"]:
            kept.append((label, names, red, g))
    return rows, kept


# ---------------------------------------------------------------------------
# stage 2: the property, re-proved of the reduction
# ---------------------------------------------------------------------------

def conservation_property(mech, y0, T, t_end, steps, rtol=1e-8):
    """PROVED iff every element total is constant along the trajectory."""
    traj = mech.integrate(y0, T, t_end, steps)
    totals = np.array([mech.comp.T @ y for y in traj])
    spread = totals.max(axis=0) - totals.min(axis=0)
    scale = np.maximum(totals.max(axis=0), 1e-30)
    return {"status": "PROVED" if float(np.max(spread / scale)) < rtol else "REFUTED",
            "max_relative_drift": float(np.max(spread / scale))}


def h2_ceiling_property(mech, y0, cond):
    """H2 can never exceed what the total H in the mixture could supply.

    Derived from element conservation instead of restated from the integrator:
    total H is invariant, so H2 <= total_H/2 holds of the full mechanism by
    construction and proving it of a reduction is a real check of the reduction.
    The first version of this had a loose ceiling (initial H2 + 1.0) and a
    mechanism that manufactured H2 came back PROVED. A property whose slack is
    larger than the effect is not a property.
    """
    traj = mech.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    i = mech.index.get("H2")
    if i is None or "H" not in mech.elements:
        return {"status": "REFUTED", "why": "H2 or element H absent"}
    h = mech.elements.index("H")
    ceiling = float((mech.comp.T @ np.asarray(y0, float))[h]) / 2.0
    col = traj[:, i]
    return {"status": "PROVED" if float(col.max()) <= ceiling + 1e-12 else "REFUTED",
            "min": float(col.min()), "max": float(col.max()), "ceiling": ceiling}


# ---------------------------------------------------------------------------
# stage 4: coverage -- measurement, no property involved
# ---------------------------------------------------------------------------

def coverage(mech, y0, cond, atol=1e-12):
    """Which species ever move. The one stage here that would still be true if
    every other stage were wrong."""
    traj = mech.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    moved, still = {}, []
    for s in mech.species:
        col = traj[:, mech.index[s]]
        d = float(col.max() - col.min())
        if d > atol:
            moved[s] = d
        else:
            still.append(s)
    return {"moved": moved, "never_moved": still}


# Reductions written by hand, exactly like a mutation table, and for the same
# reason: if the mechanism's shape changes and no pattern matches, the run must
# stop rather than invent a verdict.
CANDIDATES = [
    ("drop_R11_HO2_H",           ["R11_HO2_H"]),
    ("drop_R01_and_R11",         ["R1_H2_O2", "R11_HO2_H"]),
    ("drop_hydroperoxyl_branch", ["R10_HO2_H", "R11_HO2_H", "R12_HO2_O", "R13_HO2_OH"]),
    ("drop_third_body_H",        ["R6_H_H_M", "R7_H_OH_M"]),
]


WINDOW_GRID = (1e-3, 1e-2, 1e-1, 1.0, 10.0)


def window_sweep(path, candidates=None, tolerance=1e-6, grid=WINDOW_GRID):
    """Does the gate discriminate, and from which observation window on?

    This is the measurement that decided the module's own conclusion. Over an
    observation window short enough that the network has barely moved, every
    reduction passes: at t_end = 1 ms two of fourteen reactions carry any flux,
    H2 conversion is 0.1%, and a dropped reaction is invisible because there was
    nothing to see. Widening the window is the only way to find out whether the
    gate is measuring the reduction or its own blindness, and the answer is in
    this table rather than in a sentence.
    """
    mod = load_mechanism(path)
    candidates = candidates or CANDIDATES
    rows = []
    for t_end in grid:
        full = Mechanism(mod.SPECIES, mod.REACTIONS, mod.ELEMENTS)
        cond = dict(mod.CONDITIONS)
        cond["t_end"] = float(t_end)
        full.normalize_rates(cond["T"])
        y0 = initial_vector(full, mod.INITIAL)
        traj = full.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
        fl = full.fluxes(y0, cond["T"], cond["t_end"], cond["steps"])
        tot = float(np.max(np.abs(fl))) or 1.0
        active = int((np.abs(fl) / tot > 1e-4).sum())
        converted = float((mod.INITIAL.get("H2", 1.0) -
                           traj[-1][full.index["H2"]]) / mod.INITIAL.get("H2", 1.0))
        row = {"t_end": float(t_end), "active_reactions": active,
               "h2_converted": converted, "candidates": {}}
        for label, names in candidates:
            red = full.reduce_dropping(names)
            g = equivalence_gate(full, red, y0, cond, tol=tolerance)
            row["candidates"][label] = {"trajectory": g["trajectory"],
                                        "ok": g["ok"],
                                        "balance": g["balance"]}
        row["accepted"] = sum(1 for c in row["candidates"].values() if c["ok"])
        row["rk4_self_error"] = integration_floor(
            path, t_end, cond["steps"])["rk4_self_error"]
        rows.append(row)
    return {"tolerance": tolerance, "grid": rows,
            "finding": ("a reduction that passes at the shortest window is not "
                        "evidence about the reduction unless the window is long "
                        "enough for the dropped step to matter")}


def stages(path, t_end=None, perturb="R3_H2_O", factor=1e6, tolerance=1e-6):
    """All four stages at one observation window. Returns data, never prints."""
    mod = load_mechanism(path)
    full = Mechanism(mod.SPECIES, mod.REACTIONS, mod.ELEMENTS)
    cond = dict(mod.CONDITIONS)
    if t_end is not None:
        cond["t_end"] = float(t_end)
    full.normalize_rates(cond["T"])
    y0 = initial_vector(full, mod.INITIAL)

    rows, kept = search(full, y0, cond, CANDIDATES, tol=tolerance)

    # A reduction can pass the gate for the wrong reason: if the reactions it
    # drops carry no flux, the trajectory agrees because nothing was removed in
    # any observable sense. That is a vacuous pass and it is reported as one --
    # the share of total flux carried by the dropped steps makes it a number
    # instead of an opinion.
    base_flux = full.fluxes(y0, cond["T"], cond["t_end"], cond["steps"])
    tot_flux = float(np.sum(np.abs(base_flux))) or 1.0
    idx = {r[0]: i for i, r in enumerate(full.reactions)}
    for row in rows:
        dropped = sum(abs(base_flux[idx[n]]) for n in row["dropped"] if n in idx)
        row["dropped_flux_share"] = dropped / tot_flux
        row["vacuous"] = bool(row["ok"] and row["dropped_flux_share"] < 1e-3)

    out = {"mechanism": os.path.basename(path), "species": full.n,
           "reactions": full.m, "search": rows,
           "full_balance": full.balance_residual(),
           "conditions": {"T": cond["T"], "t_end": cond["t_end"],
                          "steps": cond["steps"], "tolerance": tolerance}}
    if not kept:
        out["verdict"] = "NO_WIN"
        return out

    label, names, red, g = kept[0]
    out["kept"] = {"label": label, "dropped": list(names), "reactions": red.m,
                   "balance": g["balance"], "trajectory": g["trajectory"],
                   "dropped_flux_share": next(r["dropped_flux_share"] for r in rows
                                              if r["label"] == label),
                   "vacuous": next(r["vacuous"] for r in rows if r["label"] == label)}

    # Compared at 16x the configured resolution, on purpose: at 400 steps the two
    # families differ by their own discretisation error (~1e-3 here), which reads
    # like a defect and is not one. Two integrators agree at the resolution they
    # are run at, not at the resolution the tolerance was written for.
    try:
        resolved = dict(cond)
        resolved["steps"] = cond["steps"] * 16
        tf = full.integrate(y0, cond["T"], cond["t_end"], resolved["steps"])
        ts = full.integrate_stiff(y0, cond["T"], cond["t_end"], resolved["steps"])
        out["integrator_agreement"] = float(
            np.max(np.abs(tf - ts)) / max(float(np.max(np.abs(tf))), 1e-30))
        out["integrator_steps"] = resolved["steps"]
    except Exception as exc:
        out["integrator_agreement"] = None
        out["integrator_note"] = "second integrator unavailable: %s" % exc
    out["error_floor"] = integration_floor(path, cond["t_end"], cond["steps"])

    out["property"] = {"full": h2_ceiling_property(full, y0, cond),
                       "reduced": h2_ceiling_property(red, y0, cond)}

    fl = red.fluxes(y0, cond["T"], cond["t_end"], cond["steps"])
    top = red.reactions[int(np.argmax(fl))][0]
    perturbed = {}
    for rname in (top, perturb):
        broken = red.with_scaled_rate(rname, factor)
        gb = equivalence_gate(red, broken, y0, cond, tol=tolerance)
        perturbed[rname] = {"factor": factor, **gb,
                            "control": "REFUTED" if not gb["ok"] else "PROVED",
                            "flux": float(fl[[r[0] for r in red.reactions].index(rname)])}
    out["perturb"] = perturbed
    out["perturb_observable"] = top
    out["coverage"] = coverage(red, y0, cond)
    ok = (out["property"]["reduced"]["status"] == "PROVED"
          and out["property"]["full"]["status"] == "PROVED"
          and perturbed[top]["control"] == "REFUTED"
          and full.balance_residual() < 1e-12)
    out["verdict"] = "SPINE_OK" if ok else "SPINE_FAIL"
    return out


def run_spine(path, perturb="R3_H2_O", factor=1e6, tolerance=1e-6, sweep=True):
    """The headline run at the configured window, plus the window sweep, plus --
    when the configured window accepts nothing -- the same four stages at the
    shortest window that accepts something.

    That last part is not a convenience. It is the finding: agreement is a
    property of a *horizon*, so a report that says "the reduction is equivalent"
    without naming the window has not said anything yet.
    """
    out = stages(path, None, perturb, factor, tolerance)
    if sweep:
        out["window_sweep"] = window_sweep(path, tolerance=tolerance)
    if out["verdict"] == "NO_WIN" and sweep:
        for row in out["window_sweep"]["grid"]:
            if row["accepted"]:
                cond = stages(path, row["t_end"], perturb, factor, tolerance)
                cond["note"] = ("accepted only at t_end = %g s, where the network "
                                "has converted %.1f%% of the H2 and %d of %d "
                                "reactions carry flux: a horizon-conditional result, "
                                "not an unconditional one"
                                % (row["t_end"], 100 * row["h2_converted"],
                                   row["active_reactions"], out["reactions"]))
                out["conditional"] = cond
                break
    return out
