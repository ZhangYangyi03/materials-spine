
"""A testable law of strength growth, and the control that could have killed it.

The question is narrow on purpose: given one measured value of a concrete mix --
its 28-day strength -- and the ages at which that same mix was tested, which
one-parameter law of time predicts its OTHER ages, and is the parameter a
property of the material or of the mix?

Two things are separated here that are usually reported together:

    the FORM of the age law     tested by leaving one age out of the same mix
                                (the mix's own remaining ages are the only input)
    the CARRIER of its parameter tested on mixes the fit never saw, against a
                                permuted-target control

The candidate forms, each carrying exactly one fitted number and anchored
through the mix's own S(28):

    power law      S = S28 * (t/28)^n
    CEB-FIP        S = S28 * exp(k * (1 - sqrt(28/t)))          <- bounded axis
    ACI-209 style  S = S28 / (1 + b * (t-28)/28)

The bounded form is the one that can saturate. That distinction is the whole
claim: an unbounded power law extrapolates past what the data does.
"""
import math, csv, os, collections, statistics as st, random

from . import data as _data
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(_data.TABLES_DIR, "concrete")

def mixes(path=None):
    """Group the table by mix proportions. One mix = one physical material."""
    path = path or os.path.join(DATA, "concrete.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    key = lambda r: (r["cement"], r["slag"], r["fly_ash"], r["water"],
                     r["superplasticizer"], r["coarse_agg"], r["fine_agg"])
    g = collections.defaultdict(list)
    for r in rows:
        g[key(r)].append((int(r["age"]), float(r["strength"])))
    out = []
    for k, v in g.items():
        c, s, fa, w, sp, ca, fine = [float(x) for x in k]
        b = c + s + fa
        if b <= 0:
            continue
        out.append({"pts": sorted(v), "s28": dict(v).get(28),
                    "sfrac": s / b, "ffrac": fa / b, "wb": w / b})
    return out

# ---- the three one-parameter forms -------------------------------------------------
def fit(pts, s28, kind):
    """Least squares on the mean of the other points; parameter only, anchor fixed."""
    num = den = 0.0
    for a, y in pts:
        if a == 28 or y <= 0:
            continue
        if kind == "power":
            x, z = math.log(a / 28.0), math.log(y / s28)
        elif kind == "ceb":
            x, z = 1 - math.sqrt(28.0 / a), math.log(y / s28)
        else:
            x, z = (a - 28.0) / 28.0, s28 / y - 1.0
        num += x * z
        den += x * x
    return num / den if den > 0 else None

def predict(a, p, s28, kind):
    if p is None:
        return None
    if kind == "power":
        v = s28 * (a / 28.0) ** p
    elif kind == "ceb":
        v = s28 * math.exp(p * (1 - math.sqrt(28.0 / a)))
    else:
        d = 1.0 + p * (a - 28.0) / 28.0
        v = s28 / d if d > 1e-6 else None
    return v if (v is not None and 0 < v < 400) else None

KINDS = ("power", "ceb", "aci")

def leave_one_age_out(ms, min_other=2):
    """Every (mix, age) whose age can be predicted from the mix's other ages."""
    for m in ms:
        s28 = m["s28"]
        if not s28:
            continue
        for a, y in m["pts"]:
            if a == 28:
                continue
            rest = [(aa, yy) for aa, yy in m["pts"] if aa != a]
            if len([1 for aa, _ in rest if aa != 28]) < min_other:
                continue
            row = {"age": a, "truth": y, "mix": m}
            for kind in KINDS:
                row[kind] = predict(a, fit(rest, s28, kind), s28, kind)
            yield row

def form_scores(ms):
    """The FORM test: which law predicts an untested age of the SAME mix."""
    rows = list(leave_one_age_out(ms))
    out = {}
    for kind in KINDS:
        e = [abs(r[kind] - r["truth"]) for r in rows if r[kind] is not None]
        out[kind] = {"n": len(e), "mean_abs_err": st.mean(e), "median_abs_err": st.median(e),
                     "p90": sorted(e)[int(.9 * len(e))]}
    # paired, on the rows where both exist
    both = [r for r in rows if r["ceb"] is not None and r["power"] is not None]
    out["ceb_beats_power"] = sum(1 for r in both if abs(r["ceb"] - r["truth"]) < abs(r["power"] - r["truth"]))
    out["paired_n"] = len(both)
    return out, rows

# ---- the CARRIER test ---------------------------------------------------------------
def _ols(X, y, lam=1e-9):
    k = len(X[0]); A = [[0.0] * k for _ in range(k)]; b = [0.0] * k
    for x, yy in zip(X, y):
        for i in range(k):
            b[i] += x[i] * yy
            for j in range(k):
                A[i][j] += x[i] * x[j]
    for i in range(k):
        A[i][i] += lam
    for i in range(k):
        p = max(range(i, k), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]; b[i], b[p] = b[p], b[i]
        for r in range(i + 1, k):
            f = A[r][i] / A[i][i]
            for c in range(i, k):
                A[r][c] -= f * A[i][c]
            b[r] -= f * b[i]
    x = [0.0] * k
    for i in reversed(range(k)):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, k))) / A[i][i]
    return x

def _r2(X, y, c):
    yh = [sum(a * b for a, b in zip(x, c)) for x in X]
    my = st.mean(y)
    return 1 - sum((v - a) ** 2 for v, a in zip(y, yh)) / sum((v - my) ** 2 for v in y)

def exponent(m):
    """CEB-FIP k for this mix, from all its ages. The parameter we try to explain."""
    s28 = m["s28"]
    if not s28:
        return None
    num = den = 0.0
    for a, y in m["pts"]:
        if a == 28 or y <= 0:
            continue
        x = 1 - math.sqrt(28.0 / a)
        num += x * math.log(y / s28)
        den += x * x
    return num / den if den > 0 else None

FEATURES = {
    "level (ln S28)":        lambda m: [1.0, math.log(m["s28"])],
    "water/binder":          lambda m: [1.0, m["wb"]],
    "composition (s,f,w/b)": lambda m: [1.0, m["sfrac"], m["ffrac"], m["wb"]],
    "level + composition":   lambda m: [1.0, math.log(m["s28"]), m["sfrac"], m["ffrac"], m["wb"]],
}

def carrier(ms, seed=42, nperm=300, cut=0.75):
    """nperm=0 skips the permutation control (for stability studies)."""
    """Held-out R2 for each candidate carrier of k, with a permuted-target control."""
    rec = [m for m in ms if m["s28"] and len(m["pts"]) >= 3 and exponent(m) is not None]
    rnd = random.Random(seed)
    idx = list(range(len(rec))); rnd.shuffle(idx)
    c = int(cut * len(idx)); tr, te = idx[:c], idx[c:]
    out = {}
    for name, F in FEATURES.items():
        Xtr = [F(rec[i]) for i in tr]; ytr = [exponent(rec[i]) for i in tr]
        Xte = [F(rec[i]) for i in te]; yte = [exponent(rec[i]) for i in te]
        coef = _ols(Xtr, ytr)
        obs = _r2(Xte, yte, coef)
        if nperm:
            perm = []
            for _ in range(nperm):
                y2 = ytr[:]; rnd.shuffle(y2)
                perm.append(_r2(Xte, yte, _ols(Xtr, y2)))
            perm.sort()
            out[name] = {"r2_heldout": obs, "coef": coef,
                         "perm_mean": st.mean(perm), "perm_p95": perm[int(.95 * len(perm))],
                         "p": sum(1 for q in perm if q >= obs) / len(perm)}
        else:
            out[name] = {"r2_heldout": obs, "coef": coef}
    return out, len(rec)

def end_members(ms):
    """k at the pure end members, from the composition fit on every mix."""
    rec = [m for m in ms if m["s28"] and len(m["pts"]) >= 3 and exponent(m) is not None]
    X = [FEATURES["composition (s,f,w/b)"](m) for m in rec]
    y = [exponent(m) for m in rec]
    c = _ols(X, y)
    return {"k_cement": c[0], "k_slag": c[0] + c[1], "k_flyash": c[0] + c[2],
            "k_water_coef": c[3], "n_mixes": len(rec)}


# ---- the controls -------------------------------------------------------------------
def form_control(ms, seed=0, nrep=60):
    """Destroy the time-ORDER of each mix's own measurements, keep the multiset.

    If a law wins only because it has a shape that can sit on any increasing data,
    or because mixes with more ages dominate, this control removes the ordering and
    the winner must disappear. The anchor S(28) stays where it is -- only the other
    ages get each other's strengths.
    """
    out = {k: [] for k in KINDS}
    for rep in range(nrep):
        rnd = random.Random(seed * 1000 + rep)
        scr = []
        for m in ms:
            pts = list(m["pts"])
            if not m["s28"]:
                continue
            others = [(a, y) for a, y in pts if a != 28]
            if len(others) < 3:
                continue
            ys = [y for _, y in others]
            rnd.shuffle(ys)
            scr.append({"pts": sorted(others_i for others_i in zip([a for a, _ in others], ys)) + [(28, m["s28"])],
                        "s28": m["s28"], "sfrac": m["sfrac"], "ffrac": m["ffrac"], "wb": m["wb"]})
        r2 = list(leave_one_age_out(scr))
        for k in KINDS:
            e = [abs(r[k] - r["truth"]) for r in r2 if r[k] is not None]
            if e:
                out[k].append(st.mean(e))
    return out

def end_member_control(ms, seed=1, nrep=300):
    """Same fit, composition labels shuffled between mixes. Report the spread of the
    end-member numbers that comes out -- that spread is what 'no relation' looks like."""
    rec = [m for m in ms if m["s28"] and len(m["pts"]) >= 3 and exponent(m) is not None]
    F = FEATURES["composition (s,f,w/b)"]
    X = [F(m) for m in rec]
    y = [exponent(m) for m in rec]
    rnd = random.Random(seed)
    out = []
    for _ in range(nrep):
        y2 = y[:]
        rnd.shuffle(y2)
        c = _ols(X, y2)
        out.append((c[0], c[0] + c[1], c[0] + c[2]))
    return [{"k_cement": a, "k_slag": b, "k_flyash": c} for a, b, c in out]
