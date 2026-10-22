"""The benchmark that produces every number the README quotes.

One entry point, one ledger, no hand-written results. `python -m mspine bench`
runs it; `python -m mspine report` prints what came out. Nothing in docs/ is
typed by hand from a run that is not on the chain.

Design of the comparison, stated here so a reader can attack it:

  * every arm gets the same initial design, the same budget, the same batch size;
  * the fitting policy is stated, because it is a variable and not a detail: by
    default the surrogate's hyperparameters are fitted once on the initial design
    and then held, and `--refits 0,8` runs both policies side by side so the cost
    of re-tuning after every batch is a number in the report rather than an
    assumption in a docstring;
  * each (arm, seed) is run twice: once on the real responses and once on
    responses permuted across cells. The permuted run is the negative control and
    it is not optional;
  * the number quoted is the held-out screening gain (see campaign.assess), not
    the campaign's own best-so-far, because the latter is measured on rows the
    method chose.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from . import campaign, data, ledger as L

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def _ledger(path=None):
    return L.Ledger(path or os.path.join(RESULTS, "ledger.jsonl"))


def bench_campaign(tables=("buchwald", "concrete"),
                   arms=("random", "ofat", "fill", "ei", "logei", "ucb", "qei"),
                   budget=64, seeds=(0, 1, 2), batch=8, control=True, quiet=False,
                   refits=(0,)):
    led = _ledger()
    env = led.environment()
    led.append("environment", "host", env)
    out = {"budget": budget, "batch": batch, "seeds": list(seeds),
           "env": env, "tables": {}}
    for name in tables:
        space, codes, y, meta = data.load(name)
        led.append("table", name, {k: (v if not isinstance(v, dict) else v)
                                   for k, v in meta.items()})
        t0 = time.time()
        recs = campaign.compare(space, codes, y, name, arms=arms, budget=budget,
                                seeds=seeds, batch=batch, controls=control,
                                refits=tuple(refits))
        agg = campaign.aggregate(recs)
        out["tables"][name] = {
            "meta": {k: v for k, v in meta.items() if k != "snap_max_frac_of_range"},
            "arms": agg,
            "seconds": round(time.time() - t0, 1),
            "records": [{k: r.get(k) for k in ("arm", "seed", "control", "refit",
                                               "heldout")} for r in recs],
        }
        for a in agg:
            led.append("campaign",
                       "%s/%s/refit%d/%s" % (name, a["arm"], a["refit"],
                                             "control" if a["control"] else "real"),
                       a, data_sha256=meta["sha256"])
        if not quiet:
            print(f"[campaign] {name}: {len(recs)} runs in {out['tables'][name]['seconds']}s",
                  flush=True)
            for a in agg:
                if a["control"]:
                    continue
                print("   %-7s refit=%-2d gain %+6.2f +- %5.2f   rank %+0.3f   top10 %.2f"
                      % (a["arm"], a["refit"], a["gain_mean"], a["gain_sd"],
                         a["rank_corr_mean"], a["top10_recovery_mean"]), flush=True)
    return out


def bench_laws(quiet=False):
    from . import laws
    led = _ledger()
    out = {}
    t0 = time.time()
    form = laws.hydration_form()
    car = laws.hydration_carrier()
    stab = laws.hydration_stability()
    ctrl = laws.hydration_controls()
    out["hydration"] = {"form": form, "carrier": car, "stability": stab, "controls": ctrl}
    led.append("law", "hydration/form", form)
    led.append("law", "hydration/carrier", {k: v for k, v in car.items()})
    led.append("law", "hydration/stability", stab)
    led.append("law", "hydration/controls", ctrl)
    if not quiet:
        print("[laws] hydration in %.1fs" % (time.time() - t0), flush=True)
        for k, v in form["forms"].items():
            print("   %-6s mean |err| %6.2f MPa (n=%d)" % (k, v["mean_abs_err"], v["n"]))
        for n, d in car["carriers"].items():
            print("   carrier %-22s R2 %+0.3f  p=%.3f" % (n, d["r2_heldout"], d["p"]))
    t0 = time.time()
    w = laws.weibull_test()
    wn = laws.weibull_null()
    led.append("law", "weibull/test", w)
    led.append("law", "weibull/null", wn)
    out["weibull"] = {"test": w, "null": wn}
    if not quiet:
        print("[laws] weibull in %.1fs" % (time.time() - t0), flush=True)
        print("   m_scaling %.2f  m_scatter %.2f  ratio %.2f (required 1.00)"
              % (w["m_scaling"], w["m_scatter"], w["ratio"]))
        print("   null: median m_scaling %.1f, %.0f%% of nulls at or above observed"
              % (wn["m_scaling_median"], 100 * wn["share_at_or_above_observed"]))
    return out


def bench_spine(quiet=False):
    from . import mechanism as M
    led = _ledger()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "mechanisms", "hydrogen.py")
    out = M.run_spine(path)
    led.append("mechanism", out.get("mechanism", "hydrogen"), out)
    if not quiet:
        print("[spine] %s: %d species, %d reactions, window t_end = %g s"
              % (out["mechanism"], out["species"], out["reactions"],
                 out["conditions"]["t_end"]))
        print("   %-26s %-6s %-9s %-9s %s"
              % ("reduction", "verdict", "balance", "trajectory", "dropped flux"))
        for r in out["search"]:
            print("   %-26s %-6s %-9.1e %-9.1e %.1e"
                  % (r["label"], "KEPT" if r["ok"] else "refused", r["balance"],
                     r["trajectory"], r["dropped_flux_share"]))
        print("   verdict at this window: %s" % out["verdict"])
        if "window_sweep" in out:
            print("   gate discrimination against the observation window:")
            print("   %-9s %-7s %-9s %s" % ("t_end", "active", "H2 conv", "accepted"))
            for row in out["window_sweep"]["grid"]:
                print("   %-9.0e %-7d %-9.3f %d of %d"
                      % (row["t_end"], row["active_reactions"], row["h2_converted"],
                         row["accepted"], len(row["candidates"])))
        c = out.get("conditional")
        if c:
            print("   -> only an earlier window accepts anything (%g s); kept %s,"
                  % (c["conditions"]["t_end"], c["kept"]["label"]))
            print("      which drops a step carrying %.1e of the flux: vacuous = %s"
                  % (c["kept"]["dropped_flux_share"], c["kept"]["vacuous"]))
            print("      %s" % c["note"])
    return out


def summary_from_ledger(out=None, need=3):
    """Rebuild the summary from the chain alone.

    The part files are a convenience; the chain is the record. This reads the last
    row per (table, arm, refit, control) out of results/ledger.jsonl and prints
    the table -- so a benchmark run in ten pieces across ten processes, or across
    ten restarts, is one result. `need` is the number of seeds below which a group
    is reported as underpowered instead of being quietly averaged.
    """
    led = _ledger()
    ok, bad = led.verify()
    tables, envs = {}, []
    for r in led.of_kind("environment"):
        envs.append(r["payload"])
    for r in led.of_kind("table"):
        # the newest table row wins: a re-run re-measured the shape, and an older
        # row may predate a field (response_range) that the report now prints.
        tables.setdefault(r["subject"], {"meta": {}, "arms": {}})["meta"] = r["payload"]
    for r in led.of_kind("campaign"):
        name = r["subject"].split("/")[0]
        p = r["payload"]
        key = (p["arm"], p.get("refit", 0), bool(p["control"]))
        # append-only: a later row for the same key supersedes, and superseding is
        # visible in the history rather than an overwrite.
        tables.setdefault(name, {"meta": {}, "arms": {}})["arms"][key] = p
    out_doc = {"chain": {"ok": ok, "first_bad": bad, "rows": len(led.rows),
                         "head": led.head},
               "env": envs[-1] if envs else None, "tables": {}}
    for name, d in tables.items():
        rows = [{"arm": k[0], "refit": k[1], "control": k[2], **v}
                for k, v in d["arms"].items()]
        out_doc["tables"][name] = {
            "meta": {k: v for k, v in (d["meta"] or {}).items()
                     if k != "snap_max_frac_of_range"},
            "arms": sorted(rows, key=lambda a: (a["refit"], a["arm"], a["control"])),
            "underpowered": [("%s/%s/refit%d/%s"
                              % (name, a["arm"], a["refit"],
                                 "control" if a["control"] else "real"))
                             for a in rows if a["n"] < need],
        }
    if out:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(out_doc, f, indent=2, sort_keys=True)
    return out_doc


def merge(paths, out):
    """Combine several bench runs into one summary.

    Not a convenience. The campaigns are run in pieces (arms one group at a time)
    because a single process long enough to hold them all is long enough to be
    killed by whatever watchdog is watching, and a benchmark that only completes
    on an idle machine is a benchmark nobody can reproduce. Each piece writes its
    own rows to the same hash-chained ledger, so the combined report is stitched
    from the chain rather than from an in-memory dict that a crash would erase.
    """
    merged = {"tables": {}, "budget": None, "batch": None, "seeds": [], "env": None}
    for path in paths:
        with open(path, encoding="utf-8") as f:
            part = json.load(f)
        merged["env"] = merged["env"] or part.get("env")
        for k in ("budget", "batch"):
            if merged[k] is None:
                merged[k] = part.get(k)
        if part.get("seeds") and len(part["seeds"]) > len(merged["seeds"]):
            merged["seeds"] = part["seeds"]
        for name, d in part.get("tables", {}).items():
            cur = merged["tables"].setdefault(name, {"meta": d["meta"], "records": [],
                                                     "seconds": 0.0})
            cur["records"].extend(d["records"])
            cur["seconds"] = round(cur["seconds"] + d.get("seconds", 0.0), 1)
        for section in ("laws", "spine", "_smoke"):
            if section in part:
                merged[section] = part[section]
    for name, d in merged["tables"].items():
        d["arms"] = campaign.aggregate(d["records"])

    # Cross-check against the chain: every (table, arm, refit, control) group in
    # the merged summary must exist on the ledger with the same gain. A merged
    # file that disagrees with the chain is a merged file nobody can audit.
    led = _ledger()
    chain = {}
    for r in led.of_kind("campaign"):
        p = r["payload"]
        chain["%s/%s/refit%d/%s" % (r["subject"].split("/")[0], p["arm"],
                                    p.get("refit", 0),
                                    "control" if p["control"] else "real")] = p
    merged["chain_crosscheck"] = {"rows": len(chain), "checked": 0, "mismatch": []}
    for name, d in merged["tables"].items():
        for a in d["arms"]:
            key = "%s/%s/refit%d/%s" % (name, a["arm"], a["refit"],
                                        "control" if a["control"] else "real")
            ref = chain.get(key)
            if ref is None:
                merged["chain_crosscheck"]["mismatch"].append([key, "absent from chain"])
                continue
            merged["chain_crosscheck"]["checked"] += 1
            if abs(ref["gain_mean"] - a["gain_mean"]) > 1e-9:
                merged["chain_crosscheck"]["mismatch"].append(
                    [key, ref["gain_mean"], a["gain_mean"]])
    if out:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(merged, f, indent=2, sort_keys=True)
    return merged


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="bench")
    ap.add_argument("--tables", default="buchwald,concrete")
    ap.add_argument("--arms", default="random,ofat,fill,ei,logei,ucb,qei")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--budget", type=int, default=64)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(RESULTS, "summary.json"))
    ap.add_argument("--skip", default="")
    ap.add_argument("--refits", default="0",
                    help="fitting policies to run side by side, e.g. 0,8")
    a = ap.parse_args(argv)
    skip = set(s for s in a.skip.split(",") if s)
    out = {}
    if "campaign" not in skip:
        out["campaign"] = bench_campaign(
            tables=tuple(a.tables.split(",")),
            arms=tuple(a.arms.split(",")),
            budget=a.budget, seeds=tuple(range(a.seeds)), batch=a.batch,
            refits=tuple(int(x) for x in a.refits.split(",") if x.strip()))
    if "laws" not in skip:
        out["laws"] = bench_laws()
    if "spine" not in skip:
        out["spine"] = bench_spine()
    led = _ledger()
    ok, bad = led.verify()
    out["ledger"] = led.summary()
    out["ledger"]["verify"] = {"ok": ok, "first_bad": bad}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print("\nledger: %d rows, chain %s, head %s"
          % (out["ledger"]["rows"], "OK" if ok else f"BROKEN at {bad}", out["ledger"]["head"]))
    print("wrote %s" % a.out)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
