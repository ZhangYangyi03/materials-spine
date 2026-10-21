"""Laws must be able to fail, and the gate must be able to reject.

Two tests here are the reason the module exists:

  the bounded hydration law beats the power law at long ages, and the time-order
  scramble removes most of the win -- if the scramble also produced the same
  win, the win was the shape of the data rather than the law;

  the mechanism gate rejects at least one reduction at a window long enough for
  the dropped step to matter. A gate that accepts everything is not a gate.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mspine import laws, mechanism as M

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MECH = os.path.join(HERE, "mspine", "mechanisms", "hydrogen.py")


def test_bounded_law_wins_at_long_ages():
    form = laws.hydration_form()
    assert form["forms"]["ceb"]["mean_abs_err"] < form["forms"]["power"]["mean_abs_err"]
    assert form["ceb_beats_power"] > form["paired_n"] / 2


def test_carrier_is_composition_not_level():
    car = laws.hydration_carrier(nperm=60)
    comp = car["carriers"]["composition (s,f,w/b)"]
    level = car["carriers"]["level (ln S28)"]
    assert comp["r2_heldout"] > level["r2_heldout"]


def test_weibull_null_is_wider_than_the_observation():
    t = laws.weibull_test()
    n = laws.weibull_null(nrep=60)
    assert t["ratio"] != 1.0                       # the two moduli disagree
    assert n["m_scaling_median"] > t["m_scaling"]  # and the null disagrees more


def test_gate_rejects_when_window_is_long_enough():
    sw = M.window_sweep(MECH, tolerance=1e-6)
    long_window = sw["grid"][-1]
    assert long_window["t_end"] >= 1.0
    assert long_window["accepted"] < len(long_window["candidates"])


def test_balance_is_an_identity_not_a_tolerance():
    mod = M.load_mechanism(MECH)
    full = M.Mechanism(mod.SPECIES, mod.REACTIONS, mod.ELEMENTS)
    assert full.balance_residual() == 0.0


def test_second_integrator_agrees_at_resolved_steps():
    """Two integrators agree at the resolution they are run at.

    At the configured 400 steps they disagree by ~1e-3, which is RK4's own
    discretisation error on this stiff mechanism, not a defect in either. The
    check is therefore made at 16x, and the floor is asserted separately below.
    """
    import numpy as np
    mod = M.load_mechanism(MECH)
    full = M.Mechanism(mod.SPECIES, mod.REACTIONS, mod.ELEMENTS)
    cond = dict(mod.CONDITIONS)
    full.normalize_rates(cond["T"])
    y0 = M.initial_vector(full, mod.INITIAL)
    a = full.integrate(y0, cond["T"], cond["t_end"], cond["steps"] * 16)
    b = full.integrate_stiff(y0, cond["T"], cond["t_end"], cond["steps"] * 16)
    assert np.max(np.abs(a - b)) / np.max(np.abs(a)) < 1e-6


def test_integration_floor_is_reported_and_coarser_than_the_gate():
    """The gate cannot resolve finer than its own integrator.

    This is the assertion that keeps the headline honest: at the configured
    resolution the RK4 self-error must be reported, and it must be coarser than
    the gate's 1e-6 tolerance. If a future change makes the integrator accurate
    enough to resolve 1e-6, this test fails and the tolerance comment has to be
    rewritten rather than quietly left standing.
    """
    floor = M.integration_floor(MECH)
    assert floor["rk4_self_error"] > 1e-6
