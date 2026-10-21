"""Test a real engineering law against real carbon-fibre data, and report honestly.

The law is weakest-link (Weibull 1939): a fibre breaks where its worst flaw is, the
strength is extreme-value distributed, and the scatter WITHIN one gauge length and
the SHIFT of the distribution ACROSS gauge lengths are two views of ONE parameter.
If the law holds, the Weibull modulus fitted to within-length scatter and the
modulus implied by the length scaling of the characteristic strength must agree.

They do not agree here, and the direction of the disagreement is itself the result.

Data: 255 single-filament tensile tests at 30 / 50 / 70 mm gauge length, 16 fibre
batches, from the open dataset of unidirectional carbon-fibre pultruded composites
(Zenodo 18435237, CC-BY-4.0). The .xlsx workbooks were converted to one CSV; the
conversion, not the analysis, is where the risk of error sits, so both the source
SHA-256 and the derived CSV are recorded in results/manifest.json.
"""
import csv, math, os, random, statistics as st, collections

from . import data as _data
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(_data.TABLES_DIR, "carbon_fiber", "carbon_fiber_single_filament.csv")

def load():
    return [{"L": int(r["gauge_length_mm"]), "batch": r["batch"],
             "specimen": r["specimen"], "sigma": float(r["max_stress_MPa"]),
             "area": float(r["area_um2"]), "E": float(r["e_modulus_GPa"]),
             "force": float(r["force_mN"])}
            for r in csv.DictReader(open(CSV, encoding="utf-8"))]

def weibull_mle(x, iters=40):
    """2-parameter Weibull by maximum likelihood, profiled on m alone."""
    lv = [math.log(v) for v in x]
    n = len(lv); lx = sum(lv)
    def nll(m):
        mx = sum(math.exp(m * l) for l in lv) / n
        return -(n * math.log(m) - n * math.log(mx) + (m - 1) * lx - n)
    lo, hi = 0.2, 40.0
    for _ in range(iters):
        a = lo + (hi - lo) / 3; b = hi - (hi - lo) / 3
        if nll(a) < nll(b):
            hi = b
        else:
            lo = a
    m = (lo + hi) / 2
    return m, (sum(math.exp(m * l) for l in lv) / n) ** (1.0 / m)

def _slope(xs, ys):
    mx, my = st.mean(xs), st.mean(ys)
    den = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    se = math.sqrt(sum((y - (my + b * (x - mx))) ** 2 for x, y in zip(xs, ys)) / (len(xs) - 2) / den)
    return b, se

LENGTHS = (30, 50, 70)

def fits(rows, field="sigma"):
    out = {}
    for L in LENGTHS:
        x = [r[field] for r in rows if r["L"] == L]
        if len(x) < 5:
            continue
        m, s0 = weibull_mle(x)
        out[L] = {"m": m, "s0": s0, "n": len(x)}
    return out

def scatter_vs_scaling(rows, field="sigma"):
    """The two estimates of the SAME modulus, and their ratio."""
    f = fits(rows, field)
    Ls = [L for L in LENGTHS if L in f]
    b, se = _slope([math.log(L) for L in Ls], [math.log(f[L]["s0"]) for L in Ls])
    m_scaling = -1.0 / b
    m_scatter = st.mean([f[L]["m"] for L in Ls])
    return {"per_length": f, "slope": b, "slope_se": se, "m_scaling": m_scaling,
            "m_scatter": m_scatter, "ratio": m_scaling / m_scatter}

def shuffle_null(rows, field="sigma", nrep=400, seed=7):
    """Pool every strength, redistribute at random, refit. By construction there is no
    length effect in these surrogates -- so whatever ratio comes out is what 'no law'
    looks like at this sample size."""
    rnd = random.Random(seed)
    pool = sorted(r[field] for r in rows)
    sizes = [len([1 for r in rows if r["L"] == L]) for L in LENGTHS]
    out = []
    for _ in range(nrep):
        p = pool[:]; rnd.shuffle(p)
        cuts = [0, sizes[0], sizes[0] + sizes[1], len(p)]
        v = [(math.log(LENGTHS[i]), weibull_mle(p[cuts[i]:cuts[i + 1]])[1]) for i in range(3)]
        b, _ = _slope([x for x, _ in v], [math.log(y) for _, y in v])
        if b < 0:
            out.append(-1.0 / b)
    out.sort()
    return out

def per_batch_moduli(rows, min_n=9):
    g = collections.defaultdict(list)
    for r in rows:
        g[(r["L"], r["batch"])].append(r["sigma"])
    return {k: weibull_mle(v)[0] for k, v in g.items() if len(v) >= min_n}
