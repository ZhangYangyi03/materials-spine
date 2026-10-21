"""The measured tables, and the rules for touching them.

Two rules are enforced here rather than documented and hoped for:

  * a table is verified against a pinned SHA-256 on load. A mismatch means the
    bytes are not the bytes every published number was computed from, and the
    honest move is to stop, not to produce a plausible result from unknown data.
  * a zero is not a measurement of zero. In the amination table 0% yield means
    "no product detected", which is a censored reading: it is kept as a row (it
    cost real time) and flagged, never smoothed away.

Shapes. A categorical table (named ligands, named additives) becomes a
DesignSpace of named levels. A continuous table (kiln settings, ambient
conditions) becomes a grid of real settings snapped to n levels, with the
worst-case quantisation reported as a fraction of each factor's range -- a
number, not an impression.
"""
from __future__ import annotations

import csv
import hashlib
import os
from array import array
from collections import OrderedDict

import numpy as np

from .space import DesignSpace, Factor

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(HERE, "tables")

# Pinned on the bytes that ship. src/ is where the upstream licence and the
# one-step conversion live; nothing here is downloaded at run time.
PINNED = {
    "buchwald/data_table.csv":
        "fe310b50897e97578078558909efc2edaa7b98ad8939b1caad740a718398ed62",
    "ccpp/ccpp.csv":
        "79e1c4524022de9468fc84289ab2c8f9ca5030790fa25e79e13443ddce4a24f3",
    "gasturbine/gasturbine.csv":
        "8d7e812a5e964bba5e09d3e2de9baab0a494e7b60bd11ef3b64ddf9e39874203",
    "concrete/concrete.csv":
        "98072ff035fc27116079be5c0075c3144677c66877fb5d95a026f1404af80aa5",
    "carbon_fiber/carbon_fiber_single_filament.csv":
        "3640fc6415e588da03ebf015fd426a5744c138b548303ad3423a25dac767a316",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(path):
    """Check a bundled table against its pin; return the digest.

    Raises on any mismatch. An unverified table is a table whose numbers mean
    nothing, and silently skipping the check is how a benchmark becomes a story.
    """
    key = "/".join((os.path.basename(os.path.dirname(path)),
                    os.path.basename(path)))
    got = sha256(path)
    want = PINNED.get(key)
    if want is None:
        raise ValueError(f"{key} has no pinned hash; refusing to load it")
    if want == "PENDING":
        return got
    if got != want:
        raise ValueError(f"{key} hash mismatch: {got} != {want}")
    return got


def _columns(path, needed):
    """Stream a CSV into one C-double array per column.

    Not a stylistic choice. This host caps a single allocation at roughly 10 MB, so
    `list(csv.DictReader(...))` on the 36733-row turbine log dies with a
    MemoryError that reads exactly like a bug in the analysis. Reading column-wise
    into `array('d')` keeps the peak at the file's own size and makes the loader
    independent of how many rows a table has.
    """
    out = {c: array("d") for c in needed}
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        rdr = csv.reader(f)
        hdr = next(rdr)
        idx = {}
        for c in needed:
            if c not in hdr:
                raise KeyError(f"column {c!r} absent from {os.path.basename(path)}")
            idx[c] = hdr.index(c)
        for row in rdr:
            try:
                vals = [float(row[idx[c]]) for c in needed]
            except (ValueError, IndexError):
                continue
            for c, v in zip(needed, vals):
                out[c].append(v)
    return out


# ---------------------------------------------------------------------------
# categorical: the factors are named levels, so no order may be invented
# ---------------------------------------------------------------------------

def buchwald():
    """Doyle-group HTE amination: 3 bases x 4 ligands x 16 halides x 24 additives.

    4599 of the 4608 cells were run. The standard real-data benchmark for
    categorical Bayesian optimisation -- which is why it is here: the objective
    is not a function of this repo's own design.
    """
    path = os.path.join(TABLES_DIR, "buchwald", "data_table.csv")
    digest = verify(path)
    cols = ["base", "ligand", "aryl_halide", "additive"]
    levels = OrderedDict((c, []) for c in cols)
    ys = array("d")
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        for r in csv.DictReader(f):
            yv = (r.get("yield") or "").strip()
            if not yv:
                continue
            for c in cols:
                if r[c] not in levels[c]:
                    levels[c].append(r[c])
            ys.append(float(yv))
    space = DesignSpace.of([Factor.of(c, levels[c]) for c in cols])
    # second pass: with the levels fixed, encode without materialising dict rows
    codes = np.empty((len(ys), len(cols)), dtype=np.int8)
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        k = 0
        for r in csv.DictReader(f):
            if not (r.get("yield") or "").strip():
                continue
            for j, c in enumerate(cols):
                codes[k, j] = levels[c].index(r[c])
            k += 1
    y = np.asarray(ys, dtype=float)

    meta = {
        "name": "buchwald", "kind": "categorical",
        "source": "doylelab/rxnpredict data_table.csv (MIT); Ahneman et al., Science 2018",
        "sha256": digest, "rows": int(len(y)), "grid": space.size,
        "cells_observed": int(len({tuple(c) for c in codes.tolist()})),
        "target": "isolated yield (%)", "categorical": True,
        "censored": int((y == 0).sum()),
        "censored_note": "yield == 0 means no product detected, not a measured 0.0",
    }
    return space, codes, y, meta


# ---------------------------------------------------------------------------
# continuous: the factors are settings, so the metric must be kept
# ---------------------------------------------------------------------------

def continuous_table(path, factor_cols, target_col, n_levels=12, q=(0.01, 0.99),
                     source="", rename=None, positive=False):
    """A table of measured settings as a continuous design space.

    Each factor is cut into `n_levels` equal-width levels across [q_lo, q_hi] of
    the observed range, each row snaps to the nearest level, and a cell is
    answerable if at least one real row landed in it. The half-step worst case
    is reported per factor as `snap` so the campaign's resolution is a number.
    """
    rename = rename or {}
    digest = verify(path)
    cols_d = _columns(path, list(factor_cols) + [target_col])
    y = np.frombuffer(cols_d[target_col], dtype=float)
    if not len(y):
        raise ValueError(f"no usable rows in {path}")
    if positive:
        keep = y > 0
    else:
        keep = np.ones(len(y), bool)
    y = np.ascontiguousarray(y[keep])

    factors, cols, snap = [], [], {}
    for j, name in enumerate(factor_cols):
        col = np.frombuffer(cols_d[name], dtype=float)[keep]
        lo, hi = float(np.quantile(col, q[0])), float(np.quantile(col, q[1]))
        if not (hi > lo):
            lo, hi = float(col.min()), float(col.max())
        label = rename.get(name, name)
        f = Factor.continuous_of(label, lo, hi, n=n_levels)
        # nearest level by arithmetic, not by an (n_rows x n_levels) distance
        # matrix: the levels are equal-width, so rounding is exact and the
        # broadcast copy is what a memory-capped run dies on.
        step = (hi - lo) / (n_levels - 1)
        cols.append(np.clip(np.rint((col - lo) / step), 0, n_levels - 1).astype(int))
        snap[label] = ((hi - lo) / (n_levels - 1)) / 2.0 / (hi - lo)
        factors.append(f)
    space = DesignSpace.of(factors)
    codes = np.asarray(cols, dtype=np.int16).T
    # fold codes to a flat cell id by Horner's rule in int64 arithmetic we never
    # materialise as a matrix: the canonical index for a mixed-radix grid
    card = np.asarray(space.cardinality, dtype=np.int64)
    flat = np.zeros(len(codes), dtype=np.int64)
    for j in range(codes.shape[1]):
        flat = flat * card[j] + codes[:, j]
    cells = int(np.unique(flat).size)
    del flat
    meta = {
        "kind": "continuous", "source": source, "sha256": digest,
        "rows": int(len(y)), "factors": len(factor_cols),
        "levels_per_factor": n_levels, "grid": int(space.size),
        "cells_observed": cells, "coverage": cells / space.size,
        "target": target_col, "categorical": False,
        "quantisation": "equal-width levels across the [q_lo, q_hi] range, "
                        "rows snapped to the nearest level",
        "snap_max_frac_of_range": snap,
    }
    return space, codes, y, meta


def ccpp():
    """Combined-cycle power plant: four ambient conditions -> net hourly output."""
    return continuous_table(
        os.path.join(TABLES_DIR, "ccpp", "ccpp.csv"),
        ["AT", "V", "AP", "RH"], "PE", n_levels=16, q=(0.01, 0.99),
        source="UCI #294 Combined Cycle Power Plant (Kaya/Tufekci 2012, CC BY 4.0); "
               "9568 hourly records",
        rename={"AT": "ambient_temp", "V": "exhaust_vacuum",
                "AP": "ambient_pressure", "RH": "rel_humidity"})


def gasturbine():
    """Gas turbine: eight operating variables -> NOx emission."""
    return continuous_table(
        os.path.join(TABLES_DIR, "gasturbine", "gasturbine.csv"),
        ["AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP"], "NOX",
        n_levels=5, q=(0.01, 0.99), positive=True,
        source="UCI #551 Gas Turbine CO and NOx Emission (CC BY 4.0); 36733 records",
        rename={"AT": "ambient_temp", "AP": "ambient_pressure",
                "AH": "ambient_humidity", "AFDP": "air_filter_dp",
                "GTEP": "gt_exhaust_pressure", "TIT": "turbine_inlet_temp",
                "TAT": "turbine_after_temp", "CDP": "compressor_discharge_pressure"})


def concrete():
    """Concrete mix design: eight mix variables and curing age -> strength."""
    return continuous_table(
        os.path.join(TABLES_DIR, "concrete", "concrete.csv"),
        ["cement", "slag", "fly_ash", "water", "superplasticizer",
         "coarse_agg", "fine_agg", "age"], "strength",
        n_levels=4, q=(0.0, 1.0),
        source="UCI #165 Concrete Compressive Strength (Yeh 2007, CC BY 4.0); "
               "1030 mix designs")


REGISTRY = {
    "buchwald": buchwald,
    "ccpp": ccpp,
    "gasturbine": gasturbine,
    "concrete": concrete,
}


def load(name):
    if name not in REGISTRY:
        raise KeyError(f"unknown table {name!r}; have {sorted(REGISTRY)}")
    return REGISTRY[name]()


def describe():
    """One row per table: shape, coverage, provenance. The first thing to read."""
    out = []
    for name in REGISTRY:
        space, codes, y, meta = REGISTRY[name]()
        out.append({
            "table": name,
            "factors": space.dim,
            "levels": list(space.cardinality),
            "rows": meta["rows"],
            "grid": meta["grid"],
            "cells_observed": meta["cells_observed"],
            "coverage": meta["cells_observed"] / meta["grid"],
            "target": meta["target"],
            "response_range": [float(y.min()), float(y.max())],
            "sha256": meta["sha256"][:16],
        })
    return out
