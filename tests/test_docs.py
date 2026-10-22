"""The documents checked against the chain, not against a re-run.

`test_published.py` re-derives numbers by running campaigns, which costs minutes
and only covers what it was told to cover. This file is the cheaper and stricter
check: every row of the README's screening table and every row of its law table is
parsed out of the prose and compared with `results/ledger.jsonl`, the file that is
committed. A number edited in the README to look better than the run that produced
it fails here in a second, and no amount of re-running hides it.

That pairing is the point of the merge: prose without a chain is a story, and a
chain nobody reads is a log file.
"""
import json
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from mspine import bench, ledger as L

LEDGER = os.path.join(ROOT, "results", "ledger.jsonl")
pytestmark = pytest.mark.skipif(not os.path.exists(LEDGER),
                                reason="no results/ledger.jsonl; run `python -m mspine bench`")


@pytest.fixture(scope="module")
def chain():
    led = L.Ledger(LEDGER)
    doc = bench.summary_from_ledger(None)
    return {(n, a["arm"], a["refit"], a["control"]): a
            for n, t in doc["tables"].items() for a in t["arms"]}, led


def _readme():
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as f:
        return f.read()


def test_the_readme_screening_table_is_the_chain(chain):
    rows, led = chain
    text = _readme()
    block = text[text.index("    table      arm"):text.index("(abridged")]
    # derived from the tables the block actually names, so adding a table to the
    # README cannot silently stop being checked (it did, for gasturbine).
    tables_named = {l.split()[0] for l in block.splitlines()
                    if l.strip() and l.split()[0].isalpha()
                    and l.split()[0] in {k[0] for k in rows}}
    seen, bad = 0, []
    table = None
    for line in block.splitlines()[1:]:
        f = line.split()
        if not f:
            continue
        if f[0] in tables_named:
            table = f[0]
            f = f[1:]
        if not table or len(f) < 8 or not re.fullmatch(r"[+-][\d.]+", f[2]):
            continue
        arm, fit = f[0], int(f[1])
        gain, sd, ctrl = float(f[2]), float(f[3]), float(f[6])
        real = rows[(table, arm, fit, False)]
        c = rows[(table, arm, fit, True)]
        seen += 1
        for got, want, what in ((gain, real["gain_mean"], "gain"),
                                (sd, real["gain_sd"], "sd"),
                                (ctrl, c["gain_mean"], "control")):
            if abs(got - want) > 0.006:
                bad.append("%s/%s/fit%d %s: README %.2f, chain %.2f"
                           % (table, arm, fit, what, got, want))
    assert seen >= 12, "the screening table stopped being parseable (%d rows)" % seen
    assert not bad, "\n".join(bad)


def test_the_readme_law_numbers_are_the_chain(chain):
    _, led = chain
    text = _readme()
    forms = [r["payload"] for r in led.of_kind("law") if r["subject"] == "hydration/form"][-1]
    carrier = [r["payload"] for r in led.of_kind("law") if r["subject"] == "hydration/carrier"][-1]
    weibull = [r["payload"] for r in led.of_kind("law") if r["subject"] == "weibull/test"][-1]
    null = [r["payload"] for r in led.of_kind("law") if r["subject"] == "weibull/null"][-1]
    assert "%d mixes" % forms["mixes"] in text
    assert "%.2f MPa" % forms["forms"]["ceb"]["mean_abs_err"] in text
    assert "%.2f MPa" % forms["forms"]["power"]["mean_abs_err"] in text
    assert "closer than power on %d of %d" % (forms["ceb_beats_power"], forms["paired_n"]) in text
    comp = carrier["carriers"]["composition (s,f,w/b)"]
    assert "+%.2f," % comp["r2_heldout"] in text or "+%.2f " % comp["r2_heldout"] in text
    assert "factor of %.2f" % weibull["ratio"] in text
    assert "%.1f%% of %d" % (100 * null["share_at_or_above_observed"], null["n_rep"]) in text
