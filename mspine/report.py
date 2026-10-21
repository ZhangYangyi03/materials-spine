"""`python -m mspine report` -- print what the ledger holds, and nothing else.

Every line here is read back from results/ledger.jsonl, which is the append-only
hash chain the benchmark wrote. If a number is not on the chain, it does not
appear in this report. That is the whole difference between a result and a story:
the report has no independent memory of its own.
"""
from __future__ import annotations

import json
import os
import sys

from . import ledger as L

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(ROOT, "results", "ledger.jsonl")


def _latest(led, kind, prefix=None):
    rows = [r for r in led.of_kind(kind)
            if prefix is None or str(r["subject"]).startswith(prefix)]
    return rows


def _last(rows):
    """The report shows the latest run of each law, not every run ever made."""
    return rows[-1:]


def main(argv=None):
    path = LEDGER
    if argv:
        path = argv[0]
    if not os.path.exists(path):
        print("no ledger at %s -- run `python -m mspine bench` first" % path)
        return 1
    led = L.Ledger(path)
    ok, bad = led.verify()
    s = led.summary()
    print("ledger %s" % os.path.relpath(path, ROOT))
    print("  rows %d   head %s   chain %s"
          % (s["rows"], s["head"], "OK" if ok else "BROKEN at seq %d" % bad))
    print("  %s" % ", ".join("%s=%d" % kv for kv in sorted(s["by_kind"].items())))

    tables = {}
    for r in _latest(led, "table"):
        tables[r["subject"]] = r["payload"]
    if tables:
        print("\nTABLES (every one SHA-256 pinned and verified on load)")
        print("  %-11s %6s %8s %8s %8s  %s" %
              ("table", "rows", "grid", "cells", "cover", "source"))
        for name, m in tables.items():
            grid = m.get("grid") or 1
            cells = m.get("cells_observed") or 0
            print("  %-11s %6d %8d %8d %7.2f%%  %s"
                  % (name, m.get("rows", 0), grid, cells, 100.0 * cells / grid,
                     (m.get("source") or "")[:52]))

    rows = _latest(led, "campaign")
    if rows:
        # The ledger is append-only and a table gets benched more than once, so
        # the report keeps the LAST row per (table, arm, control) rather than
        # printing the whole history. The first version of this loop printed
        # every row and every law twice for exactly that reason.
        # A row is keyed by (arm, fitting policy, control). The policy is part of
        # the key on purpose: two runs of the same arm under different fitting
        # policies are two different arms, and collapsing them would let the
        # later row silently overwrite the earlier one.
        # A campaign row is keyed by (table, arm, fitting policy, control). The
        # policy is part of the key on purpose: one arm under two fitting policies
        # is two arms, and collapsing them would let the later row overwrite the
        # earlier one silently. Within a key the LAST row wins, because the ledger
        # is append-only and a re-run is a newer measurement of the same thing --
        # superseding, not contradicting. The seed count is printed so a 1-seed
        # check can never be mistaken for a study.
        by_table = {}
        for r in rows:
            table = str(r["subject"]).split("/")[0]
            p = r["payload"]
            key = (p["arm"], int(p.get("refit", 0)))
            by_table.setdefault(table, {}).setdefault(key, {})[p["control"]] = p
        for table in sorted(by_table):
            by = by_table[table]
            counts = {}
            for slot in by.values():
                for p in slot.values():
                    counts[p["n"]] = counts.get(p["n"], 0) + 1
            n_mode = max(counts, key=lambda k: (counts[k], -k))
            hidden = sum(1 for slot in by.values()
                         if not any(p["n"] == n_mode for p in slot.values()))
            note = ("" if not hidden else
                    "\n(%d arm/policy rows come from runs with a different seed"
                    " count, are on the chain, and are not shown)" % hidden)
            print("\nSCREENING %s -- held-out gain, %d seeds; control = the same"
                  % (table, n_mode))
            print("           pipeline on permuted responses, which must collapse to zero")
            print("  %-7s %3s %7s %8s %8s %8s %8s %7s"
                  % ("arm", "rfi", "n", "gain", "sd", "rank", "top10", "control"))
            order = sorted(by, key=lambda a: -by[a].get(False, {}).get("gain_mean", -1e9))
            for key in order:
                arm, refit = key
                p = by[key].get(False)
                if p is None or p["n"] != n_mode:
                    continue
                c = by[key].get(True)
                print("  %-7s %3d %7d %+8.2f %8.2f %+8.3f %8.2f %7s"
                      % (arm, refit, p["n"], p["gain_mean"], p["gain_sd"],
                         p["rank_corr_mean"], p["top10_recovery_mean"],
                         ("%+.2f" % c["gain_mean"]) if c else "-"))
            shown = [k for k in order if by[k].get(False)
                     and by[k][False]["n"] == n_mode]
            if shown:
                b = shown[0]
                print("  best held-out gain: %s (refit=%d, %+.2f)"
                      % (b[0], b[1], by[b][False]["gain_mean"]))
            print(note)

    for r in _last(_latest(led, "law", "hydration/form")):
        p = r["payload"]
        print("\nHYDRATION LAW (one measured S28 + the ages tested -> the other ages)")
        for k in ("power", "ceb", "aci"):
            if k in p["forms"]:
                d = p["forms"][k]
                print("  %-6s mean |err| %6.2f MPa   median %5.2f   p90 %6.2f   n=%d"
                      % (k, d["mean_abs_err"], d["median_abs_err"], d["p90"], d["n"]))
        print("  CEB-FIP (bounded) closer than the power law on %d of %d paired rows"
              % (p["ceb_beats_power"], p["paired_n"]))
    for r in _last(_latest(led, "law", "hydration/carrier")):
        p = r["payload"]
        print("  carrier of k, held-out R2 with a permuted-target control:")
        for name, d in p["carriers"].items():
            print("     %-24s %+0.3f   control95 %+0.3f   p=%.3f"
                  % (name, d["r2_heldout"], d.get("perm_p95", float("nan")),
                     d.get("p", float("nan"))))
    for r in _last(_latest(led, "law", "weibull/test")):
        p = r["payload"]
        print("\nWEAKEST-LINK (one modulus must govern scatter AND length scaling)")
        print("  m_scaling %.2f   m_scatter %.2f   ratio %.2f   required 1.00"
              % (p["m_scaling"], p["m_scatter"], p["ratio"]))
    for r in _last(_latest(led, "law", "weibull/null")):
        p = r["payload"]
        print("  null (pool and redistribute, no length effect by construction):")
        print("     median m_scaling %.1f, %.0f%% of nulls at or above the observed value"
              % (p["m_scaling_median"], 100 * p["share_at_or_above_observed"]))
        print("     -> the test as designed cannot separate law from noise; the")
        print("        reusable product is the instrument plus this null.")

    for r in _last(_latest(led, "mechanism")):
        p = r["payload"]
        print("\nMECHANISM SPINE (%s, %d species / %d reactions)"
              % (p["mechanism"], p["species"], p["reactions"]))
        for row in p["search"]:
            print("  %-26s -> %2d rxns  balance %.1e  traj %.1e  %s"
                  % (row["label"], row["reactions"], row["balance"],
                     row["trajectory"], "KEPT" if row["ok"] else "refused"))
        if "property" in p:
            print("  property H2 ceiling: full %s, reduced %s"
                  % (p["property"]["full"]["status"], p["property"]["reduced"]["status"]))
            print("  independent integrator agreement: %.2e"
                  % p.get("integrator_agreement", float("nan")))
        print("  verdict at t_end = %g s: %s"
              % (p["conditions"]["t_end"], p["verdict"]))
        sw = p.get("window_sweep")
        if sw:
            print("  gate discrimination against the observation window:")
            print("    %-9s %-7s %-9s %s" % ("t_end", "active", "H2 conv", "accepted"))
            for row in sw["grid"]:
                print("    %-9.0e %-7d %-9.3f %d of %d"
                      % (row["t_end"], row["active_reactions"], row["h2_converted"],
                         row["accepted"], len(row["candidates"])))
        c = p.get("conditional")
        if c:
            print("  only an earlier window accepts anything (%g s); the kept %s"
                  % (c["conditions"]["t_end"], c["kept"]["label"]))
            print("  drops a step carrying %.1e of the flux, so it is vacuous: %s"
                  % (c["kept"]["dropped_flux_share"], c["kept"]["vacuous"]))
            print("  %s" % c["note"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
