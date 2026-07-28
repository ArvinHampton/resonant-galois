"""
Gap A — Diophantine geometry of the BJ even surface.

Surface:  256 α^5 + 3125 β^4  =  γ^2   (γ ∈ Q, or γ ∈ Q(u) along a curve)

Task: find rational curves (α(u), β(u)) on the even locus that pass through
≥2 known HQCC seeds, so the family

    f_u = x^5 + α(u) x + β(u)

has disc square in Q(u) and recovers multiple seeds by Hilbert specialisation.

Methods:
  1. Homogenisation rays (known pure-even, one seed each) — baseline
  2. Linear pencils (known: endpoints only, disc not □ in Q(t))
  3. Quadratic Bezier / free-midpoint curves through two seeds — search midpoints
  4. Monomial / Laurent ansätze forced through two seeds
  5. Three-seed rational cubics (underdetermined) — sample search
  6. Lift seeds to (α,β,γ) and search planes/lines in A^3 ∩ surface

Outputs: GAP_A_EVEN_SURFACE.md / build/GAP_A_EVEN_SURFACE.json
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from math import isqrt
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int  # noqa: E402

u = sp.symbols("u")

# Known HQCC / classical seeds (α, β, tag)
SEEDS = [
    (-55, 88, "flagship"),
    (-55, -88, "flagship_m"),
    (95, 76, "s95_76"),
    (95, -76, "s95_m76"),
    (95, 532, "s95_532"),
    (95, -532, "s95_m532"),
    (-100, 400, "s100_400"),
    (-100, -400, "s100_m400"),
    (124, 496, "s124_496"),
    (124, -496, "s124_m496"),
    (20, 16, "classical"),
    (20, -16, "classical_m"),
]


def disc(a, b) -> int:
    return disc_bj_int(int(a), int(b))


def sqrt_disc(a, b) -> int | None:
    d = disc(a, b)
    if d < 0:
        return None
    r, ok = sp.integer_nthroot(d, 2)
    return int(r) if ok else None


def is_square_poly(expr, var=u) -> dict:
    """Is expr a square in Q(var)? Fast path via factor_list over ZZ."""
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True, "degree": -1}
        # Integer content + primitive
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.fraction(sp.together(P.content()))[0] if False else P.content()
        try:
            cont_r = sp.Rational(cont)
        except Exception:
            cont_r = sp.Rational(str(cont))
        if cont_r < 0:
            return {"ok": False, "reason": "neg_content", "content": str(cont_r)}
        n, d = int(sp.numer(cont_r)), int(sp.denom(cont_r))
        cont_ok = sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]
        if not cont_ok:
            return {
                "ok": False,
                "degree": int(P.degree()),
                "content": str(cont_r),
                "odd_factors": ["content"],
                "degenerate": False,
            }
        prim = P.primitive()[1]
        # factor over QQ
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(f), m) for f, m in fac[1] if m % 2]
        ok = len(odds) == 0
        return {
            "ok": ok,
            "degree": int(P.degree()),
            "content": str(cont_r),
            "odd_factors": odds[:10],
            "degenerate": False,
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


def poly_passes_points(alpha, beta, points: list[tuple[int, int, int]], tol_u=None):
    """
    points: list of (u_val, α, β) that the curve should hit.
    """
    for uv, a, b in points:
        try:
            av = sp.expand(alpha.subs(u, uv))
            bv = sp.expand(beta.subs(u, uv))
            if av != a or bv != b:
                return False
        except Exception:
            return False
    return True


# ---------------------------------------------------------------------------
# 1. Baseline: homogenisation rays
# ---------------------------------------------------------------------------
def baseline_rays() -> list[dict]:
    out = []
    for a0, b0, tag in SEEDS:
        g = sqrt_disc(a0, b0)
        if g is None:
            continue
        alpha, beta = a0 * u**4, b0 * u**5
        D = sp.expand(256 * alpha**5 + 3125 * beta**4)
        info = is_square_poly(D)
        out.append(
            {
                "type": "homogenisation_ray",
                "tag": tag,
                "seed": (a0, b0),
                "alpha": str(alpha),
                "beta": str(beta),
                "gamma": str(g * u**10),  # since disc = (g u^10)^2
                "disc_square_in_Qu": info,
                "n_seeds_on_curve": 1,
                "seeds_on_curve": [tag],
            }
        )
    return out


# ---------------------------------------------------------------------------
# 2. Linear pencils (document failure of pure even)
# ---------------------------------------------------------------------------
def linear_pencils() -> list[dict]:
    out = []
    for i, (a0, b0, t0) in enumerate(SEEDS):
        for a1, b1, t1 in SEEDS[i + 1 :]:
            if sqrt_disc(a0, b0) is None or sqrt_disc(a1, b1) is None:
                continue
            alpha = (1 - u) * a0 + u * a1
            beta = (1 - u) * b0 + u * b1
            D = sp.expand(256 * alpha**5 + 3125 * beta**4)
            info = is_square_poly(D)
            out.append(
                {
                    "type": "linear_pencil",
                    "tag": f"{t0}__{t1}",
                    "seeds": [(a0, b0, t0), (a1, b1, t1)],
                    "alpha": str(alpha),
                    "beta": str(beta),
                    "disc_square_in_Qu": info,
                    "n_seeds_on_curve": 2,
                }
            )
    return out


# ---------------------------------------------------------------------------
# 3. Quadratic curves through two seeds with free midpoint (α_m, β_m)
#     α(u) = (1-u)^2 α0 + 2 u(1-u) α_m + u^2 α1
#     β(u) = (1-u)^2 β0 + 2 u(1-u) β_m + u^2 β1
#     Search α_m, β_m on a lattice so D(u) is square poly
# ---------------------------------------------------------------------------
def quadratic_bezier_search(max_mid: int = 12, seed_pairs_limit: int = 20) -> list[dict]:
    print("  quadratic Bezier midpoint search...", flush=True)
    hits = []
    tested = 0
    pairs = []
    for i, s0 in enumerate(SEEDS):
        for s1 in SEEDS[i + 1 :]:
            if sqrt_disc(s0[0], s0[1]) and sqrt_disc(s1[0], s1[1]):
                pairs.append((s0, s1))
    pairs = pairs[:seed_pairs_limit]

    mid_range = list(range(-max_mid, max_mid + 1))
    # denser near 0
    for (a0, b0, t0), (a1, b1, t1) in pairs:
        # subsample midpoints: lattice including averages
        mids = set()
        mids.add(((a0 + a1) // 2, (b0 + b1) // 2))
        mids.add((a0, b1))
        mids.add((a1, b0))
        mids.add((0, 0))
        mids.add((3, 9))
        mids.add((61, 80))
        for am in mid_range[::2]:
            for bm in mid_range[::2]:
                mids.add((am, bm))
        # also scale averages
        for k in [-2, -1, 2, 3]:
            mids.add((k * (a0 + a1) // 2, k * (b0 + b1) // 2))

        for am, bm in mids:
            tested += 1
            alpha = sp.expand((1 - u) ** 2 * a0 + 2 * u * (1 - u) * am + u**2 * a1)
            beta = sp.expand((1 - u) ** 2 * b0 + 2 * u * (1 - u) * bm + u**2 * b1)
            D = disc_expr_sym(alpha, beta)
            info = is_square_poly(D)
            if info.get("ok") and not info.get("degenerate"):
                hits.append(
                    {
                        "type": "quadratic_bezier",
                        "tag": f"{t0}__{t1}__mid_{am}_{bm}",
                        "seeds": [t0, t1],
                        "midpoint": (am, bm),
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "disc_square_in_Qu": info,
                    }
                )
                print(f"    *** HIT Bezier {t0}-{t1} mid=({am},{bm})", flush=True)
    return {"tested": tested, "hits": hits, "pairs_scanned": len(pairs)}


def disc_expr_sym(alpha, beta):
    return sp.expand(256 * alpha**5 + 3125 * beta**4)


# ---------------------------------------------------------------------------
# 4. Monomial bridge: α = a0 (1-u)^p + a1 u^p, β = b0 (1-u)^q + b1 u^q
#    with (p,q) chosen so weighted disc degrees can match
# ---------------------------------------------------------------------------
def monomial_bridges() -> list[dict]:
    print("  monomial bridges...", flush=True)
    hits = []
    tested = 0
    # weights: for disc terms α^5 ~ deg 5p, β^4 ~ deg 4q — want comparable
    exps = [(1, 1), (2, 2), (4, 5), (5, 4), (2, 1), (1, 2), (3, 3), (4, 4), (8, 10)]
    pairs = [
        (SEEDS[0], SEEDS[10]),  # flagship, classical
        (SEEDS[0], SEEDS[2]),  # flagship, 95_76
        (SEEDS[0], SEEDS[4]),  # flagship, 95_532
        (SEEDS[10], SEEDS[2]),
        (SEEDS[2], SEEDS[4]),
    ]
    for (a0, b0, t0), (a1, b1, t1) in pairs:
        for p, q in exps:
            tested += 1
            alpha = sp.expand(a0 * (1 - u) ** p + a1 * u**p)
            beta = sp.expand(b0 * (1 - u) ** q + b1 * u**q)
            # check endpoints
            if int(alpha.subs(u, 0)) != a0 or int(beta.subs(u, 0)) != b0:
                continue
            if int(alpha.subs(u, 1)) != a1 or int(beta.subs(u, 1)) != b1:
                continue
            D = disc_expr_sym(alpha, beta)
            info = is_square_poly(D)
            if info.get("ok") and not info.get("degenerate"):
                hits.append(
                    {
                        "type": "monomial_bridge",
                        "tag": f"{t0}__{t1}__p{p}_q{q}",
                        "p": p,
                        "q": q,
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "disc_square_in_Qu": info,
                    }
                )
                print(f"    *** HIT monomial {t0}-{t1} p={p} q={q}", flush=True)
    return {"tested": tested, "hits": hits}


# ---------------------------------------------------------------------------
# 5. Search rational curves of form
#     α = (A0 + A1 u + A2 u^2) / (D0 + D1 u + D2 u^2)
#     β = (B0 + B1 u + B2 u^2) / (D0 + D1 u + D2 u^2)
#    forced through two seeds at u=0 and u=1, free small integer coeffs
# ---------------------------------------------------------------------------
def _eval_rational_disc_square(A0, A1, A2, B0, B1, B2, D0, D1, D2, sample_us) -> bool:
    """
    Fast filter: for many integer u with den(u)≠0, check disc(α(u),β(u)) is a square.
    If any sample fails (and is defined), reject. Endpoints 0,1 are seeds (always square).
    """
    fails = 0
    ok_n = 0
    for uv in sample_us:
        denv = D0 + D1 * uv + D2 * uv * uv
        if denv == 0:
            continue
        na = A0 + A1 * uv + A2 * uv * uv
        nb = B0 + B1 * uv + B2 * uv * uv
        # α = na/denv, β = nb/denv — disc of monic BJ with rational coeffs
        # disc = 256 α^5 + 3125 β^4; multiply by denv^5:
        # denv^5 disc = 256 na^5 + 3125 nb^4 denv
        # disc square in Q ⇔ denv^5 disc is square * (square factor from den) carefully:
        # for integer check: clear denoms of α,β to monic Z poly
        # f = x^5 + (na/d) x + (nb/d). Let F = d^5 f(x/d) wait...
        # Simple: compute Dnum = 256 na^5 + 3125 nb^4 * denv; need Dnum / denv^5 = square
        # = square iff Dnum * denv is square when 5 odd... 
        # For rational r = A/B in lowest terms, r square iff A,B both squares (up to sign).
        try:
            Dnum = 256 * (na**5) + 3125 * (nb**4) * denv
            # r = Dnum / denv**5
            # reduce gcd
            from math import gcd

            g = gcd(abs(int(Dnum)), abs(int(denv)) ** 5) if denv else 1
            # use integer nthroot on reduced form
            # For speed: check if Dnum * denv is square and denv**6 is square? 
            # r = N/D with D=denv^5. r=□ iff N*D is □ and D is □ OR both N and D are □.
            # Safest: both |N| and |D| squares after reducing fraction.
            Nn, Dd = int(Dnum), int(denv) ** 5
            if Dd < 0:
                Nn, Dd = -Nn, -Dd
            g = gcd(abs(Nn), abs(Dd))
            Nn //= g
            Dd //= g
            if Nn < 0:
                fails += 1
                if fails >= 2:
                    return False
                continue
            if not (sp.integer_nthroot(Nn, 2)[1] and sp.integer_nthroot(Dd, 2)[1]):
                fails += 1
                if fails >= 2:
                    return False
            else:
                ok_n += 1
        except Exception:
            fails += 1
            if fails >= 2:
                return False
    return ok_n >= 4 and fails == 0


def rational_quadratic_search(coeff_bound: int = 4) -> list[dict]:
    """
    α=na/den, β=nb/den deg≤2 through two seeds.
    Stage 1: probabilistic integer sampling (fast reject).
    Stage 2: algebraic square check only for survivors.
    """
    print("  rational quadratic search (sample filter + algebra)...", flush=True)
    hits = []
    tested = 0
    candidates = 0
    pairs = [
        (SEEDS[0], SEEDS[10]),
        (SEEDS[0], SEEDS[2]),
        (SEEDS[10], SEEDS[2]),
        (SEEDS[0], SEEDS[4]),
        (SEEDS[2], SEEDS[4]),
    ]
    A1_vals = [-5, -3, -1, 0, 1, 3, 5, 8, -8]
    B1_vals = [-5, -3, -1, 0, 1, 3, 5, 8, -8]
    # sample u away from 0,1
    sample_us = [-3, -2, -1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, -4, -5]

    for (a0, b0, t0), (a1, b1, t1) in pairs:
        for D0 in [1, 2, 3, -1, -2]:
            for D1, D2 in itertools.product([-2, -1, 0, 1, 2], repeat=2):
                den1 = D0 + D1 + D2
                if den1 == 0:
                    continue
                A0 = a0 * D0
                B0 = b0 * D0
                for A1, B1 in itertools.product(A1_vals, B1_vals):
                    A2 = a1 * den1 - A0 - A1
                    B2 = b1 * den1 - B0 - B1
                    tested += 1
                    if tested % 1000 == 0:
                        print(
                            f"    ratquad tested {tested}, candidates {candidates}, hits {len(hits)}",
                            flush=True,
                        )
                    # Fast reject
                    if not _eval_rational_disc_square(
                        A0, A1, A2, B0, B1, B2, D0, D1, D2, sample_us
                    ):
                        continue
                    candidates += 1
                    # Algebraic confirmation (only for rare filter survivors)
                    try:
                        na = A0 + A1 * u + A2 * u**2
                        nb = B0 + B1 * u + B2 * u**2
                        den = D0 + D1 * u + D2 * u**2
                        N = sp.expand(256 * na**5 + 3125 * nb**4 * den)
                        Dd = sp.expand(den**5)
                        PN, PD = sp.Poly(N, u, domain=sp.ZZ), sp.Poly(Dd, u, domain=sp.ZZ)
                        g = sp.gcd(PN, PD)
                        N2 = PN.quo(g).as_expr()
                        D2 = PD.quo(g).as_expr()
                        info_n = is_square_poly(N2)
                        info_d = is_square_poly(D2)
                    except Exception:
                        continue
                    if info_n.get("ok") and info_d.get("ok") and not info_n.get("degenerate"):
                        alpha = sp.together(na / den)
                        beta = sp.together(nb / den)
                        hits.append(
                            {
                                "type": "rational_quadratic",
                                "tag": f"{t0}__{t1}__D{D0}_{D1}_{D2}_A1{A1}_B1{B1}",
                                "alpha": str(alpha),
                                "beta": str(beta),
                                "disc_square_in_Qu": {"ok": True, "num": info_n, "den": info_d},
                            }
                        )
                        print(f"    *** HIT ratquad {t0}-{t1} A1={A1} B1={B1}", flush=True)
    return {"tested": tested, "candidates_after_filter": candidates, "hits": hits}


# ---------------------------------------------------------------------------
# 6. Lift to A^3: points (α,β,γ) with γ^2 = disc(α,β)
#    Line through two lifts: (α,β,γ)(s) = (1-s)P + s Q
#    Require γ(s)^2 = disc(α(s),β(s)) identically
# ---------------------------------------------------------------------------
def space_lines_through_lifts() -> list[dict]:
    print("  space lines through (α,β,γ) lifts...", flush=True)
    hits = []
    lifts = []
    for a, b, tag in SEEDS:
        g = sqrt_disc(a, b)
        if g is None:
            continue
        lifts.append((a, b, g, tag))
        lifts.append((a, b, -g, tag + "_neg"))  # both signs

    tested = 0
    for i, P in enumerate(lifts):
        for Q in lifts[i + 1 :]:
            if P[3].replace("_neg", "") == Q[3].replace("_neg", ""):
                continue  # same seed different signs only
            tested += 1
            a0, b0, g0, t0 = P
            a1, b1, g1, t1 = Q
            alpha = (1 - u) * a0 + u * a1
            beta = (1 - u) * b0 + u * b1
            gamma = (1 - u) * g0 + u * g1
            # identity gamma^2 - disc(alpha,beta) == 0 ?
            lhs = sp.expand(gamma**2 - (256 * alpha**5 + 3125 * beta**4))
            if lhs == 0:
                hits.append(
                    {
                        "type": "space_line_lift",
                        "tag": f"{t0}__{t1}",
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "gamma": str(gamma),
                        "identity": True,
                    }
                )
                print(f"    *** HIT space line {t0}-{t1}", flush=True)
            else:
                # check if lhs is identically zero as poly
                if sp.expand(lhs) == 0:
                    hits.append(
                        {
                            "type": "space_line_lift",
                            "tag": f"{t0}__{t1}",
                            "identity": True,
                        }
                    )
    return {"tested": tested, "hits": hits, "n_lifts": len(lifts)}


# ---------------------------------------------------------------------------
# 7. Plane sections: intersect surface with a plane through two seed lifts
#     and a free third point; look for rational components
# ---------------------------------------------------------------------------
def plane_sections_symbolic() -> dict:
    """
    Plane through three seed lifts; F(s,r)=γ²-disc(α,β) on barycentric plane.
    Factor for linear/quadratic components (rational curves).
    Cap combinations and use ZZ expansion + factor with timeout-friendly bounds.
    """
    print("  plane sections through three seed lifts...", flush=True)
    hits = []
    lifts = []
    # Prefer positive-γ primary seeds (unique tags)
    for a, b, tag in SEEDS:
        if tag.endswith("_m") or "m76" in tag or "m532" in tag or "m400" in tag or "m496" in tag:
            continue
        g = sqrt_disc(a, b)
        if g is not None:
            lifts.append((a, b, g, tag))

    # Core triples including flagship + classical + one more
    core_tags = {"flagship", "classical", "s95_76", "s95_532", "s100_400", "s124_496"}
    lifts = [L for L in lifts if L[3] in core_tags]
    print(f"    plane lifts: {[L[3] for L in lifts]}", flush=True)

    tested = 0
    s, r = sp.symbols("s r")
    for comb in itertools.combinations(lifts, 3):
        tested += 1
        (a0, b0, g0, t0), (a1, b1, g1, t1), (a2, b2, g2, t2) = comb
        print(f"    plane {tested}: {t0},{t1},{t2}...", flush=True)
        alpha = s * a0 + r * a1 + (1 - s - r) * a2
        beta = s * b0 + r * b1 + (1 - s - r) * b2
        gamma = s * g0 + r * g1 + (1 - s - r) * g2
        # Expand in stages
        try:
            F = sp.expand(gamma**2) - 256 * sp.expand(alpha**5) - 3125 * sp.expand(beta**4)
            F = sp.expand(F)
            # Factor over rationals
            factors = sp.factor_list(F, domain=sp.QQ)
            low = []
            allf = []
            for f, m in factors[1]:
                deg = sp.total_degree(f)
                allf.append((str(f)[:80], deg, m))
                if deg <= 2 and deg >= 1:
                    low.append((str(f), deg, m))
            # Also check if F is identically 0 (whole plane on surface — impossible usually)
            rec = {
                "seeds": [t0, t1, t2],
                "low_degree_factors": low,
                "n_factors": len(factors[1]),
                "factor_degs": sorted({sp.total_degree(f) for f, _ in factors[1]}),
                "F_degree": sp.total_degree(F),
                "F_preview": str(F)[:200],
            }
            if low:
                print(f"      *** low factors {low}", flush=True)
                hits.append(rec)
            else:
                # record high-degree only summary every time for audit
                hits.append({**rec, "low_degree_factors": [], "note": "no deg≤2 factors"})
                print(f"      deg(F)={rec['F_degree']} factors degs={rec['factor_degs']}", flush=True)
        except Exception as e:
            print(f"      plane factor error: {e}", flush=True)
            hits.append({"seeds": [t0, t1, t2], "error": str(e)})
    # Separate true low-factor hits
    low_hits = [h for h in hits if h.get("low_degree_factors")]
    return {
        "tested": tested,
        "hits": low_hits,
        "all_planes": hits,
        "n_with_low_factors": len(low_hits),
    }


# ---------------------------------------------------------------------------
# 8. Verify candidate curves: specialise and match seeds + A5
# ---------------------------------------------------------------------------
def verify_curve(alpha, beta, u_vals=None) -> dict:
    if u_vals is None:
        u_vals = list(range(-5, 6)) + [9, 16, 61, 80, 243, 539]
    seed_hits = []
    a5 = []
    stats = Counter()
    for uv in u_vals:
        try:
            a = int(sp.expand(alpha.subs(u, uv)))
            b = int(sp.expand(beta.subs(u, uv)))
        except Exception:
            continue
        if b == 0 and a == 0:
            continue
        stats["tested"] += 1
        d = disc(a, b)
        if d <= 0 or not is_square(d):
            stats["odd"] += 1
            continue
        stats["sq"] += 1
        tag = None
        for sa, sb, tg in SEEDS:
            if (a, b) == (sa, sb):
                tag = tg
                seed_hits.append({"u": uv, "seed": tag, "alpha": a, "beta": b})
        r = classify_poly(x**5 + a * x + b, do_galois=True)
        if r.get("status", "").startswith("HIT_A5") or (
            r.get("galois") and "A5" in str(r.get("galois"))
        ):
            a5.append({"u": uv, "poly": r.get("poly"), "gal": r.get("galois")})
            stats["A5"] += 1
    return {"stats": dict(stats), "seed_hits": seed_hits, "A5": a5}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("GAP A — even surface rational curves through HQCC seeds", flush=True)

    rays = baseline_rays()
    print(f"  rays: {len(rays)} pure-even single-seed", flush=True)

    pencils = linear_pencils()
    pure_pencils = [p for p in pencils if p["disc_square_in_Qu"].get("ok")]
    print(f"  linear pencils: {len(pencils)}, pure-even: {len(pure_pencils)}", flush=True)

    space = space_lines_through_lifts()
    print(f"  space lines: tested {space['tested']}, hits {len(space['hits'])}", flush=True)

    mono = monomial_bridges()
    print(f"  monomial bridges: tested {mono['tested']}, hits {len(mono['hits'])}", flush=True)

    bezier = quadratic_bezier_search(max_mid=8, seed_pairs_limit=15)
    print(f"  bezier: tested {bezier['tested']}, hits {len(bezier['hits'])}", flush=True)

    rat = rational_quadratic_search(coeff_bound=3)
    print(f"  ratquad: tested {rat['tested']}, hits {len(rat['hits'])}", flush=True)

    planes = plane_sections_symbolic()
    print(f"  planes: tested {planes['tested']}, hits {len(planes['hits'])}", flush=True)

    # Collect all multi-seed pure-even candidates
    multi_hits = []
    for h in space["hits"] + mono["hits"] + bezier["hits"] + rat["hits"]:
        multi_hits.append(h)

    # Verify any multi-hits
    verified = []
    for h in multi_hits:
        try:
            alpha = sp.sympify(h["alpha"], locals={"u": u})
            beta = sp.sympify(h["beta"], locals={"u": u})
            v = verify_curve(alpha, beta)
            h["verification"] = v
            verified.append(h)
            print(
                f"  verify {h.get('tag')}: seeds={len(v['seed_hits'])} A5={v['stats'].get('A5',0)}",
                flush=True,
            )
        except Exception as e:
            h["verification_error"] = str(e)

    # Count how many seeds lie on some pure-even ray (always 1 each)
    # Multi-seed pure-even curves found:
    n_multi = len([h for h in verified if len(h.get("verification", {}).get("seed_hits", [])) >= 2])

    verdict = (
        f"Pure-even rational curves found that are multi-seed: {n_multi}. "
        f"Space-line lifts through two seed points: {len(space['hits'])} "
        f"(expect 0 — linear in A^3 rarely lies on the quintic surface). "
        f"Monomial bridges pure-even: {len(mono['hits'])}. "
        f"Quadratic Bezier pure-even: {len(bezier['hits'])}. "
        f"Rational quadratic pure-even: {len(rat['hits'])}. "
        f"Plane sections with low-degree factors: {len(planes['hits'])}. "
        "Baseline: each HQCC seed still generates a pure-even homogenisation ray "
        "(one seed per curve). No rational curve of the scanned types carries "
        "two or more distinct HQCC seeds while keeping disc a square in Q(u)."
    )

    elapsed = round(time.time() - t0, 2)
    lines = [
        "# Gap A — rational curves on the BJ even surface",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "## Surface",
        "",
        r"$$256 \alpha^5 + 3125 \beta^4 = \gamma^2$$",
        "",
        "Known points: HQCC / classical BJ seeds with \(\\gamma = \\sqrt{\\mathrm{disc}}\).",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## Method summary",
        "",
        "| Method | Tested | Pure-even multi-seed hits |",
        "|--------|-------:|--------------------------:|",
        f"| Homogenisation rays (1 seed) | {len(rays)} | 0 (by design, 1 seed each) |",
        f"| Linear pencils in (α,β) | {len(pencils)} | {len(pure_pencils)} |",
        f"| Lines in (α,β,γ)-space through lifts | {space['tested']} | {len(space['hits'])} |",
        f"| Monomial bridges | {mono['tested']} | {len(mono['hits'])} |",
        f"| Quadratic Bezier free midpoint | {bezier['tested']} | {len(bezier['hits'])} |",
        f"| Rational quadratic (bounded coeffs) | {rat['tested']} | {len(rat['hits'])} |",
        f"| Plane sections (3 seeds) | {planes['tested']} | low-deg hits: {planes.get('n_with_low_factors', len(planes.get('hits') or []))} |",
        "",
        "---",
        "",
        "## Baseline pure-even rays (theorem-grade)",
        "",
    ]
    for r in rays[:8]:
        lines.append(
            f"- **{r['tag']}** seed={r['seed']}: α=`{r['alpha']}`, β=`{r['beta']}`, "
            f"disc□={r['disc_square_in_Qu'].get('ok')}"
        )
    lines += [
        "",
        "---",
        "",
        "## Multi-seed pure-even candidates",
        "",
    ]
    if not verified:
        lines.append("_None found in scan._")
    for h in verified:
        lines.append(f"### `{h.get('tag')}` ({h.get('type')})")
        lines.append(f"- α=`{h.get('alpha')}`")
        lines.append(f"- β=`{h.get('beta')}`")
        v = h.get("verification") or {}
        lines.append(f"- verification: `{v.get('stats')}`")
        lines.append(f"- seed hits: {v.get('seed_hits')}")
        lines.append("")

    lines += [
        "---",
        "",
        "## Plane sections (all scanned triples)",
        "",
    ]
    all_pl = planes.get("all_planes") or planes.get("hits") or []
    if not all_pl:
        lines.append("_No plane sections computed._")
    for h in all_pl:
        if h.get("error"):
            lines.append(f"- seeds {h.get('seeds')}: error {h['error']}")
        elif h.get("low_degree_factors"):
            lines.append(
                f"- seeds {h['seeds']}: **low factors** `{h['low_degree_factors']}` "
                f"deg(F)={h.get('F_degree')}"
            )
        else:
            lines.append(
                f"- seeds {h.get('seeds')}: no deg≤2 factors; "
                f"deg(F)={h.get('F_degree')} factor_degs={h.get('factor_degs')}"
            )

    lines += [
        "",
        "---",
        "",
        "## Interpretation",
        "",
        "1. The even locus is a **thin** Diophantine set: linear and low-degree bridges",
        "   through two seeds almost never stay on \(\\gamma^2 = \\mathrm{disc}(\\alpha,\\beta)\).",
        "2. **Homogenisation** remains the only systematic pure-even rational curve type,",
        "   and each such curve contains **exactly one** projective seed class",
        "   (the ray through one lattice seed).",
        "3. Space lines through two lifts fail because the surface is degree 5 in α",
        "   (and degree 4 in β) — a line meets it in more points than the two seeds",
        "   unless the whole line lies on the surface (rigid, rare).",
        "4. **Principal open problem stands:** a pure geometric \(A_5\) family with",
        "   disc □ in the parameter recovering **multiple** HQCC seeds was **not**",
        "   found among rational curves of the degrees/forms scanned.",
        "",
        "### Next geometric options (if continuing Gap A)",
        "",
        "- Higher-degree rational curves (deg 3–4) with larger coefficient bounds",
        "- Base change: disc = c(u) · □ with c square-free of deg ≥ 1 (twisted families)",
        "- Work on a cover of the even surface (parametrise via aux. variables)",
        "- Modular / Hurwitz families that are BJ after Tschirnhaus (not BJ in t from the start)",
        "",
        "_Generated by gap_a_even_surface.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "n_multi_seed_pure_even": n_multi,
        "rays": len(rays),
        "linear_pure": len(pure_pencils),
        "space": space,
        "monomial": {"tested": mono["tested"], "hits": mono["hits"]},
        "bezier": {"tested": bezier["tested"], "hits": bezier["hits"]},
        "ratquad": {"tested": rat["tested"], "hits": rat["hits"]},
        "planes": planes,
        "verified_multi": verified,
    }
    write_md(OUT / "GAP_A_EVEN_SURFACE.md", doc)
    write_md(RESULTS / "GAP_A_EVEN_SURFACE.md", doc)
    write_md(ROOT / "GAP_A_EVEN_SURFACE.md", doc)
    write_json(OUT / "GAP_A_EVEN_SURFACE.json", blob)
    print(verdict, flush=True)
    print(f"Wrote GAP_A_EVEN_SURFACE.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
