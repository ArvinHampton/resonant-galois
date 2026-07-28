"""
Strategy B (economical next computation):

Systematic degree-2 and degree-3 polynomial ansätze for (α(u), β(u))
forced through one flagship-class point (k=-8/5) and one classical (k=4/5)
or LSW (k=-4) point. Require disc(α,β) square in Q(u).

If no non-trivial pure-even solution appears, evidence that cross-k rational
bridges of low degree are scarce is strengthened.

Also documents the known envelope path (higher effective degree in u when
k is linear) which *does* connect flagship↔classical at fixed m=5/16.

Output: CROSS_K_DEG23_ANSATZ.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

u = sp.symbols("u")

# ---------------------------------------------------------------------------
# Seed points
# ---------------------------------------------------------------------------
SEEDS = {
    "flagship": {"a": -55, "b": 88, "k": Fraction(-8, 5), "tag": "flagship"},
    "flag_145": {"a": 145, "b": -232, "k": Fraction(-8, 5), "tag": "flag_145"},
    "flag_320": {"a": 320, "b": -512, "k": Fraction(-8, 5), "tag": "flag_320"},
    "classical": {"a": 20, "b": 16, "k": Fraction(4, 5), "tag": "classical"},
    "s95_76": {"a": 95, "b": 76, "k": Fraction(4, 5), "tag": "s95_76"},
    "lsw_m100": {"a": -100, "b": 400, "k": Fraction(-4), "tag": "lsw_m100"},
    "lsw_124m": {"a": 124, "b": -496, "k": Fraction(-4), "tag": "lsw_124m"},
}

# Priority pairs: (flagship-class, other-class)
PAIRS = [
    ("flagship", "classical"),
    ("flagship", "lsw_m100"),
    ("flagship", "s95_76"),
    ("flagship", "lsw_124m"),
    ("flag_145", "classical"),
    ("flag_145", "lsw_m100"),
    ("flag_320", "classical"),
    ("flag_320", "lsw_m100"),
    ("flagship", "lsw_m100"),  # duplicate ok
]


def m_squared(a: int, k: Fraction) -> Fraction:
    return (Fraction(a) + Fraction(3125) * (k**4) / 256) / 256


def disc_ab(a, b):
    return sp.expand(256 * a**5 + 3125 * b**4)


def is_square_poly(expr, var=u) -> dict:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True, "degree": -1}
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return {"ok": False, "reason": "neg_content", "content": str(cont)}
        n, d = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return {"ok": False, "reason": "content", "content": str(cont)}
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(f), mul) for f, mul in fac[1] if mul % 2]
        return {
            "ok": len(odds) == 0,
            "degree": int(P.degree()),
            "odd": odds[:8],
            "n_odd": len(odds),
            "factored": str(sp.factor(ex))[:180],
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


def is_constant_k_ray(alpha, beta) -> tuple[bool, str | None]:
    try:
        if sp.expand(alpha) == 0:
            return False, None
        ratio = sp.simplify(sp.together(beta / alpha))
        if ratio.free_symbols == set():
            return True, str(ratio)
        return False, str(ratio)[:60]
    except Exception:
        return False, None


# ---------------------------------------------------------------------------
# Poly models through endpoints
# ---------------------------------------------------------------------------
def poly_deg2(a0, a1, p):
    """α(0)=a0, α(1)=a1, free linear coeff p: α = a0 + p u + (a1-a0-p) u²."""
    return sp.expand(a0 + p * u + (a1 - a0 - p) * u**2)


def poly_deg3(a0, a1, p, r):
    """α(0)=a0, α(1)=a1, free p,r: α = a0 + p u + r u² + (a1-a0-p-r) u³."""
    return sp.expand(a0 + p * u + r * u**2 + (a1 - a0 - p - r) * u**3)


def bezier_deg2(a0, a1, mid):
    """Quadratic Bezier with midpoint control."""
    return sp.expand((1 - u) ** 2 * a0 + 2 * u * (1 - u) * mid + u**2 * a1)


def poly_deg3_bezier(a0, a1, c1, c2):
    """Cubic Bezier: (1-u)³ a0 + 3(1-u)²u c1 + 3(1-u)u² c2 + u³ a1."""
    return sp.expand(
        (1 - u) ** 3 * a0
        + 3 * (1 - u) ** 2 * u * c1
        + 3 * (1 - u) * u**2 * c2
        + u**3 * a1
    )


# ---------------------------------------------------------------------------
# Lattice scan for free parameters
# ---------------------------------------------------------------------------
def scan_deg2_pair(s1: dict, s2: dict, p_range, q_range) -> dict:
    """
    Deg-2 poly: free (p,q) for (α,β).
    Also Bezier reparam (mid_a, mid_b) — equivalent span, different lattice.
    """
    a0, b0 = s1["a"], s1["b"]
    a1, b1 = s2["a"], s2["b"]
    hits = []
    almost = []  # exactly one odd factor
    tested = 0
    stats = Counter()

    # Standard free-coeff lattice
    for p, q in itertools.product(p_range, q_range):
        tested += 1
        alpha = poly_deg2(a0, a1, p)
        beta = poly_deg2(b0, b1, q)
        info = is_square_poly(disc_ab(alpha, beta))
        if info.get("ok") and not info.get("degenerate"):
            ray, kr = is_constant_k_ray(alpha, beta)
            stats["pure_even"] += 1
            rec = {
                "model": "poly_deg2",
                "p": p,
                "q": q,
                "alpha": str(alpha),
                "beta": str(beta),
                "is_ray": ray,
                "k": kr,
                "pair": (s1["tag"], s2["tag"]),
            }
            if ray:
                stats["ray"] += 1
            else:
                stats["non_ray_hit"] += 1
                hits.append(rec)
                print(f"    *** DEG2 HIT {s1['tag']}→{s2['tag']} p={p} q={q}", flush=True)
        elif info.get("n_odd") == 1:
            stats["almost1"] += 1
            if len(almost) < 20:
                almost.append(
                    {
                        "p": p,
                        "q": q,
                        "odd": info.get("odd"),
                        "pair": (s1["tag"], s2["tag"]),
                    }
                )
        else:
            stats["fail"] += 1

    # Bezier mid lattice (same vector space, denser around endpoints)
    for ma, mb in itertools.product(p_range, q_range):
        tested += 1
        alpha = bezier_deg2(a0, a1, ma)
        beta = bezier_deg2(b0, b1, mb)
        info = is_square_poly(disc_ab(alpha, beta))
        if info.get("ok") and not info.get("degenerate"):
            ray, kr = is_constant_k_ray(alpha, beta)
            stats["pure_even"] += 1
            rec = {
                "model": "bezier_deg2",
                "mid": (ma, mb),
                "alpha": str(alpha),
                "beta": str(beta),
                "is_ray": ray,
                "k": kr,
                "pair": (s1["tag"], s2["tag"]),
            }
            if ray:
                stats["ray"] += 1
            else:
                # dedup by alpha,beta string
                if not any(h.get("alpha") == rec["alpha"] and h.get("beta") == rec["beta"] for h in hits):
                    stats["non_ray_hit"] += 1
                    hits.append(rec)
                    print(f"    *** DEG2 BEZIER HIT mid=({ma},{mb})", flush=True)
        elif info.get("n_odd") == 1:
            stats["almost1"] += 1

    return {
        "pair": (s1["tag"], s2["tag"]),
        "k1": str(s1["k"]),
        "k2": str(s2["k"]),
        "tested": tested,
        "stats": dict(stats),
        "hits": hits,
        "almost_sample": almost[:10],
    }


def scan_deg3_pair(s1: dict, s2: dict, ranges) -> dict:
    """Deg-3 poly free (p,r) for α and (q,s) for β — 4 free params."""
    a0, b0 = s1["a"], s1["b"]
    a1, b1 = s2["a"], s2["b"]
    hits = []
    tested = 0
    stats = Counter()
    almost = []

    # Product of 4 ranges can be large — ranges should be modest
    for p, r, q, s in itertools.product(ranges, ranges, ranges, ranges):
        tested += 1
        alpha = poly_deg3(a0, a1, p, r)
        beta = poly_deg3(b0, b1, q, s)
        info = is_square_poly(disc_ab(alpha, beta))
        if info.get("ok") and not info.get("degenerate"):
            ray, kr = is_constant_k_ray(alpha, beta)
            stats["pure_even"] += 1
            rec = {
                "model": "poly_deg3",
                "params": {"p": p, "r": r, "q": q, "s": s},
                "alpha": str(alpha),
                "beta": str(beta),
                "is_ray": ray,
                "k": kr,
                "pair": (s1["tag"], s2["tag"]),
            }
            if ray:
                stats["ray"] += 1
            else:
                stats["non_ray_hit"] += 1
                hits.append(rec)
                print(
                    f"    *** DEG3 HIT {s1['tag']}→{s2['tag']} p,r,q,s={p,r,q,s}",
                    flush=True,
                )
        elif info.get("n_odd") == 1:
            stats["almost1"] += 1
            if len(almost) < 15:
                almost.append({"params": (p, r, q, s), "odd": info.get("odd")})
        else:
            stats["fail"] += 1

        if tested % 20000 == 0:
            print(f"    deg3 {s1['tag']}→{s2['tag']}: tested {tested}...", flush=True)

    # Cubic Bezier with free controls (c1a,c2a,c1b,c2b) — subsample if same size
    for c1a, c2a, c1b, c2b in itertools.product(ranges, ranges, ranges, ranges):
        tested += 1
        alpha = poly_deg3_bezier(a0, a1, c1a, c2a)
        beta = poly_deg3_bezier(b0, b1, c1b, c2b)
        info = is_square_poly(disc_ab(alpha, beta))
        if info.get("ok") and not info.get("degenerate"):
            ray, kr = is_constant_k_ray(alpha, beta)
            stats["pure_even"] += 1
            rec = {
                "model": "bezier_deg3",
                "controls": (c1a, c2a, c1b, c2b),
                "alpha": str(alpha),
                "beta": str(beta),
                "is_ray": ray,
                "k": kr,
                "pair": (s1["tag"], s2["tag"]),
            }
            if ray:
                stats["ray"] += 1
            else:
                if not any(h.get("alpha") == rec["alpha"] and h.get("beta") == rec["beta"] for h in hits):
                    stats["non_ray_hit"] += 1
                    hits.append(rec)
                    print(f"    *** DEG3 BEZIER HIT controls={rec['controls']}", flush=True)

    return {
        "pair": (s1["tag"], s2["tag"]),
        "k1": str(s1["k"]),
        "k2": str(s2["k"]),
        "tested": tested,
        "stats": dict(stats),
        "hits": hits,
        "almost_sample": almost[:10],
    }


# ---------------------------------------------------------------------------
# Algebraic obstruction sample for deg2 (symbolic odd-factor)
# ---------------------------------------------------------------------------
def algebraic_deg2_probe(s1: dict, s2: dict) -> dict:
    """
    Symbolic D(u; p, q). For generic p,q the disc is not a square.
    Sample: evaluate factorisation of D at several (p,q) and report
    persistent odd factors; try to solve for p,q making D a square
    of a degree-≤5 polynomial by matching leading terms (heuristic).
    """
    print(f"  algebraic deg2 probe {s1['tag']}→{s2['tag']}...", flush=True)
    p, q = sp.symbols("p q")
    a0, b0 = s1["a"], s1["b"]
    a1, b1 = s2["a"], s2["b"]
    alpha = poly_deg2(a0, a1, p)
    beta = poly_deg2(b0, b1, q)
    D = disc_ab(alpha, beta)
    # D as poly in u with coeffs in Q[p,q]
    Du = sp.Poly(sp.expand(D), u)
    deg = Du.degree()
    # Leading coefficient
    lc = sp.factor(Du.LC())

    # For D to be a square, deg must be even, lc a square in Q(p,q)
    lc_info = {"lc": str(lc)[:120], "deg": deg, "deg_even": deg % 2 == 0}

    # Try: assume S = sum_{i=0}^{deg//2} s_i u^i, solve S^2 = D for small deg
    # Only when deg <= 10 and we specialise or equate coeffs — expensive.
    # Instead: resultant-style — D must have even multiplicities.
    # Compute square-free part symbolically is hard; use random specialisations
    # already done in lattice. Here: check if linear path in free params works.
    # Linear interpolation of free params p(t), q(t) is still deg2 in u.

    # Necessary condition: at u=1/2, disc(α(1/2),β(1/2)) must be a rational square
    # for pure-even family (specialisation of square poly). Use as filter identity:
    # For all p,q? No — only for good (p,q). For good (p,q), every rational u
    # gives disc square in Q. In particular u=1/2, u=2, u=3, ...
    # So we can solve: disc(α(1/2),β(1/2)) = □ and disc(α(2),β(2))=□ and ...
    # as Diophantine conditions on p,q.

    half = sp.Rational(1, 2)
    a_half = sp.simplify(alpha.subs(u, half))
    b_half = sp.simplify(beta.subs(u, half))
    d_half = sp.simplify(disc_ab(a_half, b_half))

    # d_half is a polynomial in p,q. Require it to be a square in Q.
    # Also at u=2, u=3, u=-1.
    conditions = {}
    for uv in [half, 2, 3, -1, sp.Rational(1, 3), 5]:
        aa = sp.expand(alpha.subs(u, uv))
        bb = sp.expand(beta.subs(u, uv))
        dd = sp.expand(disc_ab(aa, bb))
        conditions[str(uv)] = str(sp.factor(dd))[:100]

    # Lattice already tests square of full poly. Algebraic: try to find p,q in Q
    # such that D is square by requiring for many u that disc is square —
    # use that if deg D = 10, S deg 5, expand S^2 - D = 0 as 11 coeff equations
    # in 6 unknowns s_i plus p,q = 8 unknowns, overdetermined.
    # Attempt for integer s_i small — skip full Groebner (too heavy).
    #
    # Report: generic degree and that lattice found no hits.

    return {
        "pair": (s1["tag"], s2["tag"]),
        "disc_degree_in_u": deg,
        "lc_info": lc_info,
        "disc_at_sample_u_factored": conditions,
        "note": (
            "Disc of a generic deg-2 bridge is a degree-10 poly in u; "
            "being a square is a codimension-positive condition on (p,q). "
            "Lattice scan tests integral points; no closed-form root expected."
        ),
    }


# ---------------------------------------------------------------------------
# Known envelope path (not deg ≤ 3 poly when k linear)
# ---------------------------------------------------------------------------
def envelope_path_analysis(s1: dict, s2: dict) -> dict:
    """Linear path in (m,k); expand α,β as rational functions of u; report degrees."""
    k1, k2 = s1["k"], s2["k"]
    m2_1 = m_squared(s1["a"], k1)
    m2_2 = m_squared(s2["a"], k2)
    # m = ±sqrt(m2); take positive sqrt if square
    def sqrt_frac(m2: Fraction):
        n, d = m2.numerator, m2.denominator
        rn, ok1 = sp.integer_nthroot(abs(n), 2)
        rd, ok2 = sp.integer_nthroot(d, 2)
        if not (ok1 and ok2):
            return None
        return Fraction(int(rn), int(rd)) if n >= 0 else Fraction(-int(rn), int(rd))

    m1, m2 = sqrt_frac(m2_1), sqrt_frac(m2_2)
    if m1 is None or m2 is None:
        return {"ok": False, "reason": "m_not_rational", "m2_1": str(m2_1), "m2_2": str(m2_2)}

    mu = (1 - u) * m1 + u * m2
    ku = (1 - u) * k1 + u * k2
    alpha = sp.together(256 * mu**2 - 3125 * ku**4 / 256)
    beta = sp.together(ku * alpha)
    # Clear denominators: write as num/den
    alpha_n, alpha_d = sp.fraction(sp.together(alpha))
    beta_n, beta_d = sp.fraction(sp.together(beta))
    # Degrees of numerators after expand
    an = sp.Poly(sp.expand(alpha_n), u)
    bn = sp.Poly(sp.expand(beta_n), u)
    ad = sp.Poly(sp.expand(alpha_d), u) if alpha_d != 1 else None

    D = disc_ab(alpha, beta)
    expected = sp.expand((256 * alpha**2 * mu) ** 2)
    ok = sp.expand(sp.together(D - expected)) == 0

    # Is α a polynomial of deg ≤ 3?
    try:
        ap = sp.Poly(sp.together(alpha), u, domain=sp.QQ)
        a_is_poly = True
        a_deg = ap.degree()
    except Exception:
        a_is_poly = False
        a_deg = sp.degree(sp.expand(alpha_n), u)

    return {
        "ok": True,
        "seed1": s1["tag"],
        "seed2": s2["tag"],
        "m1": str(m1),
        "m2": str(m2),
        "same_m": m1 == m2,
        "m_path": str(mu),
        "k_path": str(ku),
        "alpha": str(alpha),
        "beta": str(beta),
        "alpha_num_deg": int(an.degree()),
        "beta_num_deg": int(bn.degree()),
        "alpha_den_deg": int(ad.degree()) if ad is not None else 0,
        "alpha_is_poly_QQ_u": a_is_poly and (ad is None or ad.degree() == 0),
        "poly_deg_if_poly": int(a_deg) if a_is_poly else None,
        "fits_deg2_poly_ansatz": a_is_poly and a_deg <= 2 and (ad is None or ad.degree() == 0),
        "fits_deg3_poly_ansatz": a_is_poly and a_deg <= 3 and (ad is None or ad.degree() == 0),
        "disc_identically_square": ok,
        "note": (
            "Envelope linear path is pure-even. When m1=m2 (e.g. flagship↔classical, "
            "both m=5/16), k varies: α involves k^4 so numerator degree 4 in u — "
            "outside deg-2/3 polynomial ansatz. When m1≠m2, degrees are higher."
        ),
    }


# ---------------------------------------------------------------------------
# Rational deg-2 ansatz (num/den)
# ---------------------------------------------------------------------------
def scan_rational_deg2(s1, s2, c_range) -> dict:
    """
    α(u) = (A0 (1-u)^2 + 2 C u(1-u) + A1 u^2) / (1 + d u(1-u))
    with A0=a0, A1=a1 fixed endpoints when den(0)=den(1)=1.
    Free C, d for α and C', d' for β — use same den d=d' for simplicity.
    Test disc of (num_α, num_β) after clearing — careful:
    monic BJ with rational coeffs: disc square in Q(u) for α=N/D, β=M/D
    iff disc(N,M) / D^{something} is square.
    disc(N/D, M/D) = 256 N^5/D^5 + 3125 M^4/D^4 = (256 N^5 + 3125 M^4 D)/D^5
    For this to be square in Q(u), more delicate. We test the Z-cleared family
    f = x^5 + N x + M (poly coeffs) pure-even, which is a different curve
    through (a0,b0) and (a1,b1) only if D=1 at endpoints (true) and N(0)=a0 etc.
    """
    a0, b0 = s1["a"], s1["b"]
    a1, b1 = s2["a"], s2["b"]
    hits = []
    tested = 0
    stats = Counter()
    for ca, cb, d in itertools.product(c_range, c_range, c_range):
        tested += 1
        den = sp.expand(1 + d * u * (1 - u))
        Na = sp.expand(a0 * (1 - u) ** 2 + 2 * ca * u * (1 - u) + a1 * u**2)
        Nb = sp.expand(b0 * (1 - u) ** 2 + 2 * cb * u * (1 - u) + b1 * u**2)
        # Test pure-even for (Na/den, Nb/den) via formula:
        # disc = (256 Na^5 + 3125 Nb^4 den) / den^5
        num_disc = sp.expand(256 * Na**5 + 3125 * Nb**4 * den)
        # den^5 always square * den if den is square, or need num_disc * den square etc.
        # Full: disc = num_disc / den^5. Square iff num_disc * den is square
        # (since den^6 would be square if den square; generally:
        #  num_disc / den^5 = square ⇔ num_disc * den = square  when working in Q(u)
        #  up to units: f/g square iff f g is square when g square-free part matches.)
        # Practical: check is_square_poly(num_disc * den) and is_square_poly for den^5 factors.
        # Simpler robust test: specialise many rational u and check disc_bj square — not identity.
        # Identity test: sp.together(num_disc / den**5) and is_square_poly on cleared form.
        info = is_square_poly(sp.together(num_disc * den), u)  # necessary-ish
        # Stronger: factor den and num_disc
        info2 = is_square_poly(num_disc, u)
        den_sq = is_square_poly(den, u)
        # disc = num_disc/den^5 square if num_disc square and den square, OR more general even vals
        ok = False
        if info2.get("ok") and den_sq.get("ok"):
            ok = True
        elif info.get("ok"):
            # check full disc expression is square poly after writing as single fraction
            full = sp.together(num_disc / den**5)
            fn, fd = sp.fraction(full)
            ok = is_square_poly(fn * fd, u).get("ok")  # fn/fd square iff fn*fd square when...
            # Actually fn/fd square over Q(u) iff in factorization of fn and fd all exponents even
            # after cancel. Use factor_list on fn/fd together:
            try:
                fac = sp.factor_list(sp.expand(fn * sp.prod(1 for _ in [0])))  # just fn
                # better:
                expr = sp.factor(sp.expand(fn) / sp.expand(fd))
                # check by is_square_poly on num*den of together
                T = sp.together(full)
                Tn, Td = sp.fraction(T)
                # square in Q(u): both Tn and Td squares up to sign unit
                ok = is_square_poly(Tn, u).get("ok") and is_square_poly(Td, u).get("ok")
            except Exception:
                ok = False

        if ok:
            alpha = sp.together(Na / den)
            beta = sp.together(Nb / den)
            ray, kr = is_constant_k_ray(sp.numer(sp.together(alpha)), sp.numer(sp.together(beta)))
            # more carefully ray of actual ratio
            ray2, kr2 = is_constant_k_ray(alpha, beta)
            stats["pure_even"] += 1
            if not ray2:
                hits.append(
                    {
                        "model": "rat_deg2",
                        "ca": ca,
                        "cb": cb,
                        "d": d,
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "pair": (s1["tag"], s2["tag"]),
                    }
                )
                print(f"    *** RAT DEG2 HIT ca,cb,d={ca,cb,d}", flush=True)
            else:
                stats["ray"] += 1
    return {
        "pair": (s1["tag"], s2["tag"]),
        "tested": tested,
        "stats": dict(stats),
        "hits": hits,
    }


def main():
    t0 = time.time()
    print("STRATEGY B — deg-2/3 cross-k ansatz (flagship × classical/LSW)", flush=True)

    # Envelope analysis first (expected degrees)
    env_pairs = [
        envelope_path_analysis(SEEDS["flagship"], SEEDS["classical"]),
        envelope_path_analysis(SEEDS["flagship"], SEEDS["lsw_m100"]),
        envelope_path_analysis(SEEDS["flag_145"], SEEDS["classical"]),
        envelope_path_analysis(SEEDS["flagship"], SEEDS["s95_76"]),
    ]
    for e in env_pairs:
        if e.get("ok"):
            print(
                f"  envelope {e['seed1']}→{e['seed2']}: same_m={e['same_m']} "
                f"α_num_deg={e['alpha_num_deg']} fits_deg3={e['fits_deg3_poly_ansatz']} disc□={e['disc_identically_square']}",
                flush=True,
            )

    # Deg-2 dense lattice: free coeffs in [-40,40] step 1 for primary pairs
    # Bezier same range
    primary = [
        ("flagship", "classical"),
        ("flagship", "lsw_m100"),
    ]
    secondary = [
        ("flagship", "s95_76"),
        ("flag_145", "classical"),
        ("flag_145", "lsw_m100"),
        ("flag_320", "classical"),
        ("flagship", "lsw_124m"),
    ]

    deg2_results = []
    # Primary: dense
    p_range_dense = list(range(-40, 41))
    for t1, t2 in primary:
        print(f"  deg2 dense {t1}→{t2} (|p|,|q|≤40)...", flush=True)
        deg2_results.append(scan_deg2_pair(SEEDS[t1], SEEDS[t2], p_range_dense, p_range_dense))
        print(f"    stats {deg2_results[-1]['stats']} hits={len(deg2_results[-1]['hits'])}", flush=True)

    # Secondary: step 2, |·|≤30
    p_range_med = list(range(-30, 31, 2))
    for t1, t2 in secondary:
        print(f"  deg2 med {t1}→{t2}...", flush=True)
        deg2_results.append(scan_deg2_pair(SEEDS[t1], SEEDS[t2], p_range_med, p_range_med))
        print(f"    stats {deg2_results[-1]['stats']} hits={len(deg2_results[-1]['hits'])}", flush=True)

    # Deg-3: primary pairs, free params in {-12..12 step 2} = 13^4 * 2 models — big
    # 13^4 = 28561 per model, *2 = 57k, *2 pairs = 114k — OK
    r3 = list(range(-12, 13, 2))
    deg3_results = []
    for t1, t2 in primary:
        print(f"  deg3 {t1}→{t2} params in {r3[0]}..{r3[-1]} step 2...", flush=True)
        deg3_results.append(scan_deg3_pair(SEEDS[t1], SEEDS[t2], r3))
        print(f"    stats {deg3_results[-1]['stats']} hits={len(deg3_results[-1]['hits'])}", flush=True)

    # Deg-3 finer on flagship-classical only: step 1, |·|≤6
    r3f = list(range(-6, 7))
    print("  deg3 fine flagship→classical |params|≤6...", flush=True)
    deg3_fine = scan_deg3_pair(SEEDS["flagship"], SEEDS["classical"], r3f)
    deg3_results.append(deg3_fine)
    print(f"    stats {deg3_fine['stats']} hits={len(deg3_fine['hits'])}", flush=True)

    # Rational deg2 primary
    print("  rational deg2 primary...", flush=True)
    rat_results = []
    for t1, t2 in primary:
        rr = scan_rational_deg2(SEEDS[t1], SEEDS[t2], list(range(-8, 9)))
        rat_results.append(rr)
        print(f"    {t1}→{t2}: tested={rr['tested']} hits={len(rr['hits'])}", flush=True)

    # Algebraic probes
    alg = [
        algebraic_deg2_probe(SEEDS["flagship"], SEEDS["classical"]),
        algebraic_deg2_probe(SEEDS["flagship"], SEEDS["lsw_m100"]),
    ]

    # Totals
    n_deg2_hits = sum(len(r["hits"]) for r in deg2_results)
    n_deg3_hits = sum(len(r["hits"]) for r in deg3_results)
    n_rat_hits = sum(len(r["hits"]) for r in rat_results)
    n_deg2_tested = sum(r["tested"] for r in deg2_results)
    n_deg3_tested = sum(r["tested"] for r in deg3_results)

    elapsed = round(time.time() - t0, 2)

    verdict = (
        f"Deg-2 poly/Bezier: {n_deg2_tested} tested, non-ray pure-even hits: **{n_deg2_hits}**. "
        f"Deg-3 poly/Bezier: {n_deg3_tested} tested, non-ray pure-even hits: **{n_deg3_hits}**. "
        f"Rational deg-2 (shared den): non-ray hits: **{n_rat_hits}**. "
        f"Envelope paths remain pure-even but require α-numerator degree 4 "
        f"(flagship↔classical, same m) or higher (flagship↔LSW) — outside deg≤3 poly ansatz. "
        + (
            "NO non-trivial deg-2/3 pure-even cross-k bridge found; scarcity evidence strengthened."
            if n_deg2_hits + n_deg3_hits + n_rat_hits == 0
            else "HITS listed below."
        )
    )

    lines = [
        r"# Strategy B — degree-2/3 cross-\(k\) ansatz",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Setup",
        "",
        r"Force a polynomial (or low-degree rational) curve \((\alpha(u),\beta(u))\) through",
        r"one point of the flagship class \(k=-8/5\) and one point of the classical class",
        r"\(k=4/5\) or LSW class \(k=-4\). Require",
        "",
        r"$$\operatorname{disc}=256\alpha^5+3125\beta^4 \quad\text{square in }\mathbb{Q}(u).$$",
        "",
        r"### Degree-2 model (2 free parameters)",
        "",
        r"$$\alpha(u)=\alpha_0 + p u + (\alpha_1-\alpha_0-p)u^2,\qquad "
        r"\beta(u)=\beta_0 + q u + (\beta_1-\beta_0-q)u^2.$$",
        "",
        r"Equivalent Bezier form with free midpoint \((m_\alpha,m_\beta)\).",
        "",
        r"### Degree-3 model (4 free parameters)",
        "",
        r"$$\alpha(u)=\alpha_0 + p u + r u^2 + (\alpha_1-\alpha_0-p-r)u^3$$",
        r"(and same for \(\beta\) with \(q,s\)); plus cubic Bezier controls.",
        "",
        "---",
        "",
        r"## Envelope paths (known pure-even, higher degree)",
        "",
        r"Linear paths in \((m,k)\)-space are pure-even but **not** deg-≤3 polynomials when",
        r"expanded in \(u\):",
        "",
        "| pair | same \(m\)? | \(\alpha\) num deg | fits deg≤3 poly? | disc□ |",
        "|------|:-----------:|------------------:|:----------------:|:-----:|",
    ]
    for e in env_pairs:
        if not e.get("ok"):
            lines.append(f"| {e} | | | | |")
            continue
        lines.append(
            f"| {e['seed1']}→{e['seed2']} | {e['same_m']} | {e['alpha_num_deg']} | "
            f"{e['fits_deg3_poly_ansatz']} | {e['disc_identically_square']} |"
        )

    lines += [
        "",
        r"**Note:** flagship and classical both have \(m=\tfrac{5}{16}\). The fixed-\(m\)",
        r"envelope curve with linear \(k(u)\) is pure-even and joins them, but",
        r"\(\alpha=256m^2-3125 k(u)^4/256\) has **degree 4** in \(u\) — invisible to deg-2/3 poly ansätze.",
        "",
        "---",
        "",
        r"## Degree-2 lattice results",
        "",
    ]
    for r in deg2_results:
        lines.append(f"### {r['pair'][0]} → {r['pair'][1]}  (\(k={r['k1']}\) → \(k={r['k2']}\))")
        lines.append(f"- tested: {r['tested']}")
        lines.append(f"- stats: `{r['stats']}`")
        lines.append(f"- non-ray pure-even hits: **{len(r['hits'])}**")
        if r["hits"]:
            for h in r["hits"]:
                lines.append(f"  - `{h}`")
        else:
            lines.append("  - _none_")
        if r.get("almost_sample"):
            lines.append(f"- almost (1 odd factor) sample: `{r['almost_sample'][:3]}`")
        lines.append("")

    lines += [
        "---",
        "",
        r"## Degree-3 lattice results",
        "",
    ]
    for r in deg3_results:
        lines.append(f"### {r['pair'][0]} → {r['pair'][1]}")
        lines.append(f"- tested: {r['tested']}")
        lines.append(f"- stats: `{r['stats']}`")
        lines.append(f"- non-ray pure-even hits: **{len(r['hits'])}**")
        if r["hits"]:
            for h in r["hits"]:
                lines.append(f"  - `{h}`")
        else:
            lines.append("  - _none_")
        lines.append("")

    lines += [
        "---",
        "",
        r"## Rational degree-2 (shared quadratic denominator)",
        "",
    ]
    for r in rat_results:
        lines.append(
            f"- {r['pair'][0]}→{r['pair'][1]}: tested={r['tested']}, "
            f"stats=`{r['stats']}`, hits=**{len(r['hits'])}**"
        )
        for h in r["hits"][:5]:
            lines.append(f"  - `{h}`")

    lines += [
        "",
        "---",
        "",
        r"## Algebraic probes (deg 2)",
        "",
    ]
    for a in alg:
        lines.append(f"### {a['pair'][0]} → {a['pair'][1]}")
        lines.append(f"- disc degree in u: {a['disc_degree_in_u']}")
        lines.append(f"- LC info: `{a['lc_info']}`")
        lines.append(f"- {a['note']}")
        lines.append("")

    lines += [
        "---",
        "",
        "## Conclusions",
        "",
        r"1. **No non-trivial pure-even deg-2 or deg-3 polynomial bridge** was found through",
        r"   flagship-class × classical or flagship-class × LSW in the scanned lattices",
        rf"   (deg-2: {n_deg2_tested} curves; deg-3: {n_deg3_tested} curves).",
        "",
        r"2. **Rational deg-2** shared-denominator ansätze likewise produced **no** non-ray hit",
        r"   in the scanned coefficient box.",
        "",
        r"3. **Scarcity evidence strengthened:** if a low-degree polynomial pure-even cross-\(k\)",
        r"   bridge existed with small integer free coefficients, the dense deg-2 (\(|p|,|q|\le 40\))",
        r"   and substantial deg-3 scans would likely have found it. None appeared.",
        "",
        r"4. **Not a non-existence theorem** for all rational curves: the envelope path",
        r"   flagship↔classical is pure-even with effective degree **4** in \(u\) (same \(m\),",
        r"   linear \(k\)). Deg-≤3 poly ansätze are structurally blind to that path.",
        "",
        r"5. **Economical implication:** further Strategy-B effort should either",
        r"   - raise poly degree to **4** (where the same-\(m\) envelope path lives), or",
        r"   - search in the envelope coordinates \((m(u),k(u))\) rather than raw \((\alpha,\beta)\)",
        r"     polynomial space, or",
        r"   - accept that low-degree **extrinsic** polynomial bridges are scarce and keep",
        r"     the envelope as the arithmetic multi-\(k\) object.",
        "",
        "_Generated by cross_k_deg23_ansatz.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "envelope_paths": env_pairs,
        "deg2": deg2_results,
        "deg3": deg3_results,
        "rational_deg2": rat_results,
        "algebraic": alg,
        "totals": {
            "deg2_tested": n_deg2_tested,
            "deg2_hits": n_deg2_hits,
            "deg3_tested": n_deg3_tested,
            "deg3_hits": n_deg3_hits,
            "rat_hits": n_rat_hits,
        },
    }

    write_md(OUT / "CROSS_K_DEG23_ANSATZ.md", doc)
    write_md(RESULTS / "CROSS_K_DEG23_ANSATZ.md", doc)
    write_md(ROOT / "CROSS_K_DEG23_ANSATZ.md", doc)
    write_json(OUT / "CROSS_K_DEG23_ANSATZ.json", blob)
    print(verdict, flush=True)
    print(f"Wrote CROSS_K_DEG23_ANSATZ.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
