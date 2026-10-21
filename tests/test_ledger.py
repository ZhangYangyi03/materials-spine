"""The ledger must verify itself on integers for keys, and detect a rewritten row.

The first version of mspine.ledger failed this: json.dumps writes an int key as a
string in *numeric* order, and reloads it as a string that sorts *lexically*, so
the bytes on disk stopped matching the bytes that were hashed and an honest chain
reported itself as tampered. This test is the regression.
"""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mspine.ledger import Ledger, _canon, _digest


def test_int_keys_round_trip(tmp_path):
    led = Ledger(str(tmp_path / "l.jsonl"))
    led.append("k", "s", {"long_ages": {56: 1.0, 90: 2.0, 180: 3.0}})
    again = Ledger(str(tmp_path / "l.jsonl"))
    ok, bad = again.verify()
    assert ok, bad
    assert again.rows[0]["payload"]["long_ages"] == {"56": 1.0, "90": 2.0, "180": 3.0}


def test_chain_detects_an_interior_rewrite(tmp_path):
    """Three rows, edit the middle one's payload without touching its hash.

    This is the attack a hash chain does catch: the edited row's own digest no
    longer matches, so verification stops at that sequence number.
    """
    path = str(tmp_path / "l.jsonl")
    led = Ledger(path)
    for i in range(3):
        led.append("a", "s%d" % i, {"v": i})
    lines = open(path, encoding="utf-8").read().splitlines()
    row = json.loads(lines[1])
    row["payload"]["v"] = 999                 # payload edited, hash left alone
    lines[1] = _canon(row)
    open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    ok, bad = Ledger(path).verify()
    assert not ok and bad == 1


def test_rewrite_of_a_tail_needs_an_anchor(tmp_path):
    """A re-hashed tail is invisible to the chain, and that is what anchors are for.

    Rewriting the last row AND recomputing its digest leaves a chain that verifies
    on its own -- a property of hash chains, not a defect here. The defence is an
    external anchor: the head digest published where the ledger cannot reach it.
    """
    path = str(tmp_path / "l.jsonl")
    led = Ledger(path)
    led.append("a", "s1", {"v": 1})
    led.append("a", "s2", {"v": 2})
    good_head = led.anchor()["head"]

    lines = open(path, encoding="utf-8").read().splitlines()
    row = json.loads(lines[1])
    row["payload"]["v"] = 999
    body = {k: v for k, v in row.items() if k != "hash"}
    row["hash"] = _digest(body)
    lines[1] = _canon(row)
    open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")

    # self-consistently rewritten: the chain alone cannot tell
    assert Ledger(path).verify()[0] is True
    # the published head can
    ok, bad = Ledger(path).verify(anchor_head=good_head)
    assert not ok


def test_deleted_row_breaks_chain(tmp_path):
    path = str(tmp_path / "l.jsonl")
    led = Ledger(path)
    for i in range(3):
        led.append("a", "s%d" % i, {"v": i})
    lines = open(path, encoding="utf-8").read().splitlines()
    open(path, "w", encoding="utf-8").write("\n".join(lines[:1] + lines[2:]) + "\n")
    ok, bad = Ledger(path).verify()
    assert not ok


def test_pinned_tables_verify():
    from mspine import data
    for name in data.REGISTRY:
        space, codes, y, meta = data.load(name)
        assert meta["sha256"] and len(meta["sha256"]) == 64
        assert len(y) == meta["rows"]
