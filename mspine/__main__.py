"""`python -m mspine <command>` -- the four entry points.

    tables    what ships: shape, coverage, provenance hash per table
    bench     run the campaigns, the laws and the mechanism spine; write the ledger
    report    print the numbers the ledger holds, grouped by arm and by control
    verify    re-derive the hash chain and every pinned table digest
    anchor    write results/anchor.json -- the head digest, published outside the
              ledger; verify then checks the chain against it, which is the only
              way an edit to the last row becomes detectable

Nothing here needs a notebook, and nothing here writes a number that is not also
a row on the chain.
"""
from __future__ import annotations

import json
import os
import sys

from . import __version__


def _tables(argv=None):
    from . import data
    rows = data.describe()
    print(f"mspine {__version__} -- {len(rows)} measured tables, all SHA-256 pinned\n")
    hdr = ("table", "factors", "levels", "rows", "grid", "cells", "cover", "target")
    print("%-11s %7s %-18s %7s %9s %7s %7s  %s" % hdr)
    for r in rows:
        print("%-11s %7d %-18s %7d %9d %7d %6.2f%%  %s"
              % (r["table"], r["factors"], "x".join(str(x) for x in r["levels"]),
                 r["rows"], r["grid"], r["cells_observed"], 100 * r["coverage"],
                 r["target"]))
    print("\ncoverage = measured cells / product grid. It is the honest measure of")
    print("how much of a design space a table actually explores, and it ranges from")
    print("99.8% (a near-complete plate) to 0.4% (mix designs nobody crossed).")


def _bench(argv):
    from . import bench
    return bench.main(argv)


def _report(argv):
    from . import report
    return report.main(argv)


def _verify(argv):
    from . import data, ledger as L
    ok = True
    for name in data.REGISTRY:
        try:
            space, codes, y, meta = data.load(name)
            print("table  %-11s sha256 %s  verified" % (name, meta["sha256"][:16]))
        except Exception as exc:
            ok = False
            print("table  %-11s FAILED: %s" % (name, exc))
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "results", "ledger.jsonl")
    if os.path.exists(path):
        led = L.Ledger(path)
        anchor_path = os.path.join(os.path.dirname(path), "anchor.json")
        anchor_head = None
        if os.path.exists(anchor_path) and "--chain-only" not in (argv or []):
            with open(anchor_path, encoding="utf-8") as f:
                anchor = json.load(f)
            anchor_head = anchor.get("head")
            print("anchor %-11s %d rows, head %s (%s)"
                  % ("", anchor.get("rows", 0), (anchor_head or "")[:16], anchor.get("at")))
        good, bad = led.verify(anchor_head=anchor_head)
        print("ledger %-11s %d rows, chain %s" % ("", len(led.rows),
              "OK" if good else "BROKEN at seq %d" % bad))
        ok = ok and good
    else:
        print("ledger %-11s absent (run `python -m mspine bench`)" % "")
    return 0 if ok else 1


def _anchor(argv):
    from . import ledger as L
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "results", "ledger.jsonl")
    if not os.path.exists(path):
        print("no ledger at %s -- run `python -m mspine bench` first" % path)
        return 1
    led = L.Ledger(path)
    good, bad = led.verify()
    if not good:
        print("refusing to anchor a broken chain (first bad row: %s)" % bad)
        return 1
    out = os.path.join(root, "results", "anchor.json")
    a = led.anchor(out)
    print("anchored %d rows, head %s" % (a["rows"], a["head"]))
    print("wrote %s -- commit it next to the ledger it names." % os.path.relpath(out, root))
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    cmds = {"tables": _tables, "bench": _bench, "report": _report,
            "verify": _verify, "anchor": _anchor}
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in cmds:
        print("unknown command %r; have %s" % (cmd, sorted(cmds)))
        return 2
    return cmds[cmd](rest) or 0


if __name__ == "__main__":
    sys.exit(main())
