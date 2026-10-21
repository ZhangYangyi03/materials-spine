"""A hash-chained ledger: every number this package emits, and what it rests on.

A benchmark that prints numbers is a claim. A benchmark whose every number is a
row in an append-only chain, each row hashing the row before it, is a claim you
can check years later without trusting the author -- and a run that is edited
afterwards cannot hide, because the chain stops matching.

The shape of a row is fixed and small:

    seq, at, kind, subject, payload, data_sha256, prev, hash

`hash` is SHA-256 over the canonical JSON of the row without `hash`, and `prev`
is the previous row's hash. `verify_chain` re-derives all of it. This is the same
mechanism the tool ledger in my own runtime uses, applied to experiment results.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time

GENESIS = "0" * 64


def _stringify_keys(obj):
    """Recursively force dict keys to str before hashing.

    This is not defensive tidiness -- it is the bug that made the first version of
    this file's own chain fail to verify. `json.dumps` silently writes an int key
    as a string, so `{56: ..., 180: ...}` dumps as `"56"` before `"180"` (numeric
    order) and reloads as strings that sort `"180"` before `"56"`. The bytes on
    disk are then not the bytes the hash was computed over, and a perfectly honest
    ledger reports itself as tampered. A canonical form has to be canonical for
    keys too.
    """
    if isinstance(obj, dict):
        return {str(k): _stringify_keys(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_stringify_keys(v) for v in obj]
    return obj


def _canon(obj) -> str:
    return json.dumps(_stringify_keys(obj), sort_keys=True,
                      separators=(",", ":"), default=float)


def _digest(obj) -> str:
    return hashlib.sha256(_canon(obj).encode("utf-8")).hexdigest()


class Ledger:
    """Append-only JSONL chain. One file, one line per result."""

    def __init__(self, path):
        self.path = path
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self.rows = []
        if os.path.exists(path):
            with open(path, encoding="utf-8") as f:
                self.rows = [json.loads(line) for line in f if line.strip()]

    # ---- write -----------------------------------------------------------
    @property
    def head(self):
        return self.rows[-1]["hash"] if self.rows else GENESIS

    def append(self, kind, subject, payload, data_sha256=None, note=None):
        row = {
            "seq": len(self.rows),
            "at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "kind": kind,
            "subject": subject,
            "payload": payload,
            "data_sha256": data_sha256,
            "note": note,
            "prev": self.head,
        }
        row["hash"] = _digest(row)
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(_canon(row) + "\n")
        self.rows.append(row)
        return row

    def environment(self):
        return {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": _version("numpy"),
            "scipy": _version("scipy"),
            "machine": platform.machine(),
        }

    def anchor(self, path=None):
        """Publish the head digest somewhere the ledger cannot rewrite.

        An append-only chain detects an edit to any *interior* row, because the
        next row's `prev` no longer matches. It cannot detect an edit to the last
        row, or to the last N rows together -- that is a property of hash chains
        in general, not a bug here, and the fix is the one every transparency log
        uses: the head is published elsewhere (a commit, a paper, a tweet) and
        checked against later. `verify(anchor_head=...)` does the checking.
        """
        out = {"rows": len(self.rows), "head": self.head,
               "at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        if path:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(out, f, indent=2, sort_keys=True)
        return out

    # ---- read ------------------------------------------------------------
    def verify(self, anchor_head=None):
        """Re-derive the whole chain. Returns (ok, first_bad_seq or None).

        A rewritten row, an inserted row or a deleted row all show up here: the
        hash of a row covers its own content and the hash of its parent, so
        changing one line breaks every line after it.
        """
        prev = GENESIS
        for i, row in enumerate(self.rows):
            body = {k: v for k, v in row.items() if k != "hash"}
            if body.get("prev") != prev:
                return False, i
            if _digest(body) != row.get("hash"):
                return False, i
            prev = row["hash"]
        if anchor_head is not None and self.head != anchor_head:
            return False, len(self.rows) - 1 if self.rows else 0
        return True, None

    def of_kind(self, kind):
        return [r for r in self.rows if r["kind"] == kind]

    def summary(self):
        by_kind = {}
        for r in self.rows:
            by_kind[r["kind"]] = by_kind.get(r["kind"], 0) + 1
        ok, bad = self.verify()
        return {"rows": len(self.rows), "head": self.head[:16],
                "chain_ok": ok, "first_bad": bad, "by_kind": by_kind}


def _version(name):
    try:
        return __import__(name).__version__
    except Exception:
        return None
