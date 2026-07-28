"""
Abandon rigid φ. Search for a non-rigid pure-even family that specialises
onto several fixed-k arithmetic pure-even slices (Route 1, ambitious form).

Layers of search:

  A. Arithmetic envelope (known structure, not Hurwitz):
       α(m,s) = 256 m² - 3125 k(s)⁴/256,  β = k(s)·α
     pure-even over Q(m,s); recovers every fixed-k slice by freezing s.

  B. Cross-k pure-even rational curves on the BJ even surface
       256 α⁵ + 3125 β⁴ = γ²
     through seeds from *different* multi-seed k-families
     (monomial / Bezier / bi-degree ansätze).

  C. Non-rigid geometric ansätze (Hurwitz-adjacent, explicit equations):
       - Mestre-style f - t r with disc □ in Q(t)
       - two-parameter BJ with disc □ in Q(u,v)
       - icosahedral / principal quintic deformations staying pure-even
       - pencil of k-slices with geometric monodromy test

  D. Verdict against programme goal:
       pure geometric A5 + Hilbert recovery of several fixed-k families.

Output: NONRIGID_HURWITZ_SEARCH.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
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

u, v, t, m, s = sp.symbols("u v t m s")


# ---------------------------------------------------------------------------
# Multi-seed k-families (from enlarged catalogue)
# ---------------------------------------------------------------------------
# Representative seeds per multi-seed pure-even k (tag, α, β, k_str)
MULTI_K_SEEDS = {
    "-4": [
        ("lsw_m100", -100, 400),
        ("lsw_124m", 124, -496),
        ("lsw_m209", -209, 836),
        ("lsw_239", 239, -956),
    ],
    "4": [
        ("lsw4_m100", -100, -400),
        ("lsw4_124", 124, 496),
    ],
    "-8/5": [
        ("flagship", -55, 88),
        ("flag_145", 145, -232),
        ("flag_320", 320, -512),
        ("flag_1145", 1145, -1832),
    ],
    "8/5": [
        ("flagship_m", -55, -88),
        ("flag_145p", 145, 232),
    ],
    "4/5": [
        ("classical", 20, 16),
        ("s95_76", 95, 76),
        ("s220_176", 220, 176),
    ],
    "-4/5": [
        ("classical_m", 20, -16),
        ("s95_m76", 95, -76),
    ],
    "-12/5": [
        ("s180", -180, 432),
        ("s220m", 220, -528),
        ("s380", -380, 912),
    ],
    "12/5": [
        ("s180m", -180, -432),
        ("s220", 220, 528),
    ],
    "-16/5": [
        ("s55_176", -55, 176),
        ("s655", -655, 2096),
    ],
    "16/5": [
        ("s55_m176", -55, -176),
        ("s655m", -655, -2096),
    ],
}


def is_square_poly(expr, var) -> dict:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True, "degree": -1}
        # Multivariate: check as poly in primary var after expand
        if isinstance(var, (list, tuple)):
            # square in Q(var) if factored form has even exponents and content square
            fac = sp.factor_list(ex)
            cont = sp.Rational(fac[0])
            if cont < 0:
                return {"ok": False, "reason": "neg"}
            n, d = int(sp.numer(cont)), int(sp.denom(cont))
            if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
                return {"ok": False, "reason": "content", "content": str(cont)}
            odds = [(str(f), mul) for f, mul in fac[1] if mul % 2]
            return {"ok": len(odds) == 0, "odd": odds[:6], "factored": str(sp.factor(ex))[:200]}
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return {"ok": False, "reason": "neg"}
        n, d = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return {"ok": False, "reason": "content", "content": str(cont)}
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(f), mul) for f, mul in fac[1] if mul % 2]
        return {
            "ok": len(odds) == 0,
            "degree": int(P.degree()),
            "odd": odds[:6],
            "factored": str(sp.factor(ex))[:200],
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


def disc_ab(a, b):
    return sp.expand(256 * a**5 + 3125 * b**4)


# ---------------------------------------------------------------------------
# A. Arithmetic multi-k envelope (pure-even 2-param, not Hurwitz)
# ---------------------------------------------------------------------------
def arithmetic_envelope() -> dict:
    """
    Universal pure-even envelope over Q(m,s):

      k = s  (or any nonconst rational k(s))
      α = 256 m² - 3125 s⁴ / 256
      β = s · α

    disc = (256 α² m)² identically.

    Freezing s = k0 recovers the fixed-k pure-even family.
    Freezing (m,s) to catalogue values recovers multi-seed lists.
    """
    print("  A. arithmetic multi-k envelope...", flush=True)
    k = s  # parameter
    alpha = sp.together(256 * m**2 - 3125 * k**4 / 256)
    beta = sp.together(k * alpha)
    D = disc_ab(alpha, beta)
    expected = sp.expand((256 * alpha**2 * m) ** 2)
    ok = sp.expand(sp.together(D - expected)) == 0

    # Recover each multi-seed k by setting s = k0, match seeds
    recoveries = []
    for k_str, seeds in MULTI_K_SEEDS.items():
        k0 = Fraction(k_str)
        hits = []
        for tag, a0, b0 in seeds:
            # need α(m,k0)=a0, β=k0*a0=b0 (already on ray)
            # a0 = 256 m² - 3125 k0⁴/256
            # 256 m² = a0 + 3125 k0⁴/256
            rhs = Fraction(a0) + Fraction(3125) * (k0**4) / Fraction(256)
            m2 = rhs / Fraction(256)
            if m2 < 0:
                hits.append({"tag": tag, "on": False, "reason": "m2_neg"})
                continue
            num, den = m2.numerator, m2.denominator
            n_ok = sp.integer_nthroot(abs(num), 2)[1]
            d_ok = sp.integer_nthroot(den, 2)[1]
            if n_ok and d_ok:
                rn = int(sp.integer_nthroot(abs(num), 2)[0])
                rd = int(sp.integer_nthroot(den, 2)[0])
                mm = Fraction(rn, rd) if num >= 0 else Fraction(-rn, rd)
                hits.append({"tag": tag, "on": True, "m": str(mm), "a": a0, "b": b0})
            else:
                hits.append({"tag": tag, "on": False, "reason": "m2_not_square", "m2": str(m2)})
        recoveries.append(
            {
                "k": k_str,
                "n_seeds": len(seeds),
                "n_on_envelope": sum(1 for h in hits if h.get("on")),
                "hits": hits,
            }
        )

    total_on = sum(r["n_on_envelope"] for r in recoveries)
    total_seeds = sum(r["n_seeds"] for r in recoveries)

    # Sample Gal at random (m,s) integer points
    gal_samples = []
    for mv, sv in [(3, -4), (5, Fraction(-8, 5)), (5, Fraction(4, 5)), (2, Fraction(-12, 5)), (7, 4)]:
        aa = sp.simplify(alpha.subs({m: mv, s: sv}))
        bb = sp.simplify(beta.subs({m: mv, s: sv}))
        try:
            ai, bi = int(aa), int(bb)
        except Exception:
            # clear denominators if rational
            aa_r, bb_r = sp.Rational(aa), sp.Rational(bb)
            if aa_r.denominator == 1 and bb_r.denominator == 1:
                ai, bi = int(aa_r), int(bb_r)
            else:
                gal_samples.append({"m": str(mv), "s": str(sv), "alpha": str(aa), "beta": str(bb), "status": "non_Z"})
                continue
        if ai == 0:
            continue
        r = classify_poly(x**5 + ai * x + bi, do_galois=True)
        gal_samples.append(
            {
                "m": str(mv),
                "s": str(sv),
                "alpha": ai,
                "beta": bi,
                "status": r.get("status"),
                "gal": r.get("galois"),
            }
        )

    return {
        "name": "arithmetic_multi_k_envelope",
        "alpha": str(alpha),
        "beta": str(beta),
        "disc_identically_square": ok,
        "family": "x^5 + α(m,s) x + β(m,s) over Q(m,s)",
        "recovers_all_fixed_k_slices": True,
        "seed_recoveries": recoveries,
        "seeds_on_envelope": f"{total_on}/{total_seeds}",
        "gal_samples": gal_samples,
        "is_hurwitz_geometric": False,
        "note": (
            "This is the universal LSW-type parametrization of pure-even points "
            "on rays β=kα. It is pure-even and multi-k by construction, but it is "
            "an arithmetic parametrization of a subvariety of the BJ even surface — "
            "not a positive-dimensional Hurwitz space of branched covers with "
            "prescribed geometric monodromy data independent of the BJ embedding."
        ),
    }


# ---------------------------------------------------------------------------
# B. Cross-k pure-even rational curves on the even surface
# ---------------------------------------------------------------------------
def cross_k_curve_search() -> dict:
    """
    Seek α(u), β(u) ∈ Q(u) with:
      disc(α,β) square in Q(u)
      curve not contained in a single ray β=kα
      passes through ≥1 seed from k1 and ≥1 seed from k2 ≠ k1
    """
    print("  B. cross-k pure-even rational curves...", flush=True)
    # Build seed pairs from different multi-seed k's
    k_list = list(MULTI_K_SEEDS.keys())
    pairs = []
    for i, k1 in enumerate(k_list):
        for k2 in k_list[i + 1 :]:
            # skip pure sign-flip pairs as "same geometry" optional — still try a few
            s1 = MULTI_K_SEEDS[k1][0]
            s2 = MULTI_K_SEEDS[k2][0]
            pairs.append((k1, s1, k2, s2))
            # also flagship vs LSW, flagship vs classical explicitly with more seeds
            if len(MULTI_K_SEEDS[k1]) > 1 and len(MULTI_K_SEEDS[k2]) > 1:
                pairs.append((k1, MULTI_K_SEEDS[k1][1], k2, MULTI_K_SEEDS[k2][1]))

    # Priority pairs: different |k| families (not just sign flips)
    priority = [
        ("-8/5", "flagship", -55, 88, "-4", "lsw_m100", -100, 400),
        ("-8/5", "flagship", -55, 88, "4/5", "classical", 20, 16),
        ("-8/5", "flagship", -55, 88, "-12/5", "s180", -180, 432),
        ("-8/5", "flag_145", 145, -232, "-4", "lsw_124m", 124, -496),
        ("4/5", "classical", 20, 16, "-4", "lsw_m100", -100, 400),
        ("4/5", "s95_76", 95, 76, "-12/5", "s180", -180, 432),
        ("-4", "lsw_m100", -100, 400, "-12/5", "s220m", 220, -528),
        ("-8/5", "flagship", -55, 88, "-16/5", "s55_176", -55, 176),  # same α, different k!
        ("-8/5", "flag_320", 320, -512, "4/5", "classical", 20, 16),
        ("-12/5", "s180", -180, 432, "4", "lsw4_124", 124, 496),
    ]

    exp_pairs = [
        (1, 1), (2, 2), (3, 3), (4, 4), (5, 5),
        (4, 5), (5, 4), (2, 1), (1, 2), (3, 2), (2, 3),
        (8, 10), (6, 6), (3, 4), (4, 3), (2, 5), (5, 2),
        (1, 4), (4, 1), (6, 5), (5, 6),
    ]

    hits = []
    tested = 0
    near = []  # disc almost square (one odd factor)

    def test_curve(alpha, beta, meta):
        nonlocal tested
        tested += 1
        info = is_square_poly(disc_ab(alpha, beta), u)
        if info.get("ok") and not info.get("degenerate"):
            # check not single ray: β/α constant?
            ratio = sp.simplify(sp.together(beta / alpha)) if alpha != 0 else None
            is_ray = ratio is not None and ratio.free_symbols == set()
            rec = {
                **meta,
                "alpha": str(alpha),
                "beta": str(beta),
                "disc_ok": True,
                "is_constant_k_ray": bool(is_ray),
                "k_ratio": str(ratio) if ratio is not None else None,
                "info": {k: info[k] for k in info if k != "factored"},
            }
            if not is_ray:
                hits.append(rec)
                print(f"    *** CROSS-K PURE-EVEN *** {meta}", flush=True)
            else:
                near.append({**rec, "note": "pure-even but single-k ray"})
        elif info.get("odd") and len(info.get("odd") or []) == 1:
            near.append({**meta, "almost": info.get("odd"), "alpha": str(alpha)[:80]})

    for k1, tag1, a1, b1, k2, tag2, a2, b2 in priority:
        meta_base = {"k1": k1, "seed1": tag1, "k2": k2, "seed2": tag2, "ab1": (a1, b1), "ab2": (a2, b2)}
        # monomial bridges
        for p, q in exp_pairs:
            alpha = sp.expand(a1 * (1 - u) ** p + a2 * u**p)
            beta = sp.expand(b1 * (1 - u) ** q + b2 * u**q)
            test_curve(alpha, beta, {**meta_base, "type": "monomial", "exps": (p, q)})

        # weighted homo-style
        alpha = sp.expand(a1 * (1 - u) ** 4 + a2 * u**4)
        beta = sp.expand(b1 * (1 - u) ** 5 + b2 * u**5)
        test_curve(alpha, beta, {**meta_base, "type": "homo_4_5"})

        # Bezier midpoints on small lattice
        for am, bm in itertools.product(range(-20, 21, 4), range(-20, 21, 4)):
            alpha = sp.expand((1 - u) ** 2 * a1 + 2 * u * (1 - u) * am + u**2 * a2)
            beta = sp.expand((1 - u) ** 2 * b1 + 2 * u * (1 - u) * bm + u**2 * b2)
            test_curve(alpha, beta, {**meta_base, "type": "bezier", "mid": (am, bm)})

        # rational quadratic: α = (a1(1-u)² + c u(1-u) + a2 u²)/(1+d u(1-u)) etc — simple
        for c, d in itertools.product(range(-10, 11, 5), range(-6, 7, 3)):
            if d == -4:  # avoid denom issues sometimes
                pass
            den = sp.expand(1 + d * u * (1 - u))
            if den == 0:
                continue
            alpha = sp.together(
                (a1 * (1 - u) ** 2 + c * u * (1 - u) + a2 * u**2) / den
            )
            beta = sp.together(
                (b1 * (1 - u) ** 2 + c * u * (1 - u) + b2 * u**2) / den
            )
            # clear: work with numerators as coeffs of monic after scale — test disc of (A,B) = num
            A = sp.expand(a1 * (1 - u) ** 2 + c * u * (1 - u) + a2 * u**2)
            B = sp.expand(b1 * (1 - u) ** 2 + c * u * (1 - u) + b2 * u**2)
            # disc of (A/den, B/den) = den^{-something} * poly; square condition on num disc careful
            # Use A,B as BJ coeffs directly (different curve through same endpoints if den=1 at 0,1)
            test_curve(A, B, {**meta_base, "type": "quad_num", "c": c, "d": d})

    # Also: space curves α(u),β(u),γ(u) with γ² = disc — already testing disc square

    return {
        "tested": tested,
        "cross_k_pure_even_hits": hits,
        "n_hits": len(hits),
        "single_k_ray_pure_even": [n for n in near if n.get("is_constant_k_ray")],
        "n_almost": len([n for n in near if "almost" in n]),
        "priority_pairs": len(priority),
    }


# ---------------------------------------------------------------------------
# C. Geometric / Hurwitz-adjacent ansätze
# ---------------------------------------------------------------------------
def mestre_two_seed_search() -> dict:
    """
    Mestre: f - t r for f a seed poly, r low degree.
    Seek disc(f - t r) square in Q(t) AND specialization hitting another k-family seed.
    """
    print("  C1. Mestre f-tr cross-k...", flush=True)
    hits = []
    tested = 0
    # f = flagship, r degree ≤ 2, small coeffs
    f0 = x**5 - 55 * x + 88
    for deg in (0, 1, 2):
        ranges = range(-3, 4)
        coeff_lists = itertools.product(ranges, repeat=deg + 1)
        for coeffs in coeff_lists:
            if all(c == 0 for c in coeffs):
                continue
            r = sum(c * x**i for i, c in enumerate(coeffs))
            # f - t r
            ft = sp.expand(f0 - t * r)
            pol = sp.Poly(ft, x)
            if pol.degree() != 5 or pol.LC() != 1:
                continue
            tested += 1
            try:
                D = sp.factor(sp.expand(pol.discriminant()))
            except Exception:
                continue
            info = is_square_poly(D, t)
            if not info.get("ok") or info.get("degenerate"):
                continue
            # pure-even family found — check if any other multi-k seed appears at rational t
            seed_hits = []
            for k_str, seeds in MULTI_K_SEEDS.items():
                for tag, a, b in seeds:
                    # does x^5+a x+b equal f - t0 r for some t0?
                    # coeffs must match: only possible if r is BJ-shaped (deg≤1 in x for BJ stay)
                    target = x**5 + a * x + b
                    diff = sp.expand(f0 - target)
                    # diff = t0 * r  ⇒ r must be proportional to diff
                    if r == 0:
                        continue
                    # equate: f0 - target = t0 * r as polynomials
                    # so r must be c * (f0 - target) and t0 = 1/c, or match coeffs
                    try:
                        q, rem = sp.div(diff, r, domain=sp.QQ)
                        if rem == 0 and q.free_symbols <= set():
                            t0 = sp.Rational(q)
                            seed_hits.append({"tag": tag, "k": k_str, "t": str(t0), "a": a, "b": b})
                    except Exception:
                        # coeff solve for t: compare one coeff
                        pass
            rec = {
                "r": str(r),
                "disc_square": True,
                "seed_hits": seed_hits,
                "n_seeds": len({h["tag"] for h in seed_hits}),
            }
            if rec["n_seeds"] >= 2:
                hits.append(rec)
                print(f"    *** Mestre multi-seed pure-even r={r} seeds={seed_hits}", flush=True)
            elif rec["n_seeds"] == 1 and seed_hits[0]["tag"] != "flagship":
                hits.append({**rec, "note": "single other seed"})

    return {"tested": tested, "hits": hits, "n_multi": sum(1 for h in hits if h.get("n_seeds", 0) >= 2)}


def biparam_bj_ansatz() -> dict:
    """
    Two-parameter BJ ansätze with disc □ in Q(u,v):
      α = p(u,v), β = q(u,v) low degree
    Known pure-even: α = 256 m² - 3125 k⁴/256 with m,k rational in (u,v).
    Search other forms: α = a u² + b v² + c, β = d u v, etc.
    """
    print("  C2. biparam BJ pure-even ansätze...", flush=True)
    hits = []
    tested = 0
    # grid of simple biparam forms
    forms = []
    # form: α = A u^2 + B v^2 + C, β = D u v + E u + F v
    for A, B, C, D, E, F in itertools.product(
        [-2, -1, 0, 1, 2, 4, 5, 16, 256],
        [-2, -1, 0, 1, 2, 4],
        [-3125, -80, -5, 0, 1, 5, 16, 20],
        [-4, -2, -1, 0, 1, 2, 4],
        [-2, -1, 0, 1, 2],
        [-2, -1, 0, 1, 2],
    ):
        if A == 0 and B == 0:
            continue
        if D == 0 and E == 0 and F == 0:
            continue
        forms.append((A, B, C, D, E, F))

    # subsample for runtime
    step = max(1, len(forms) // 800)
    for A, B, C, D, E, F in forms[::step]:
        tested += 1
        alpha = A * u**2 + B * v**2 + C
        beta = D * u * v + E * u + F * v
        info = is_square_poly(disc_ab(alpha, beta), (u, v))
        if info.get("ok") and not info.get("degenerate"):
            # which k-slices / seeds recovered?
            # freeze v=0: β = E u, α = A u² + C → k = β/α not const generally
            seed_hits = []
            for k_str, seeds in MULTI_K_SEEDS.items():
                for tag, a0, b0 in seeds:
                    # solve A u² + B v² + C = a0, D u v + E u + F v = b0
                    sols = sp.solve(
                        [alpha - a0, beta - b0],
                        [u, v],
                        dict=True,
                    )
                    rat_sols = []
                    for sol in sols:
                        try:
                            uu, vv = sol[u], sol[v]
                            if uu.is_rational and vv.is_rational:
                                rat_sols.append((str(uu), str(vv)))
                        except Exception:
                            pass
                    if rat_sols:
                        seed_hits.append({"tag": tag, "k": k_str, "sols": rat_sols[:3]})
            n_k = len({h["k"] for h in seed_hits})
            rec = {
                "A": A, "B": B, "C": C, "D": D, "E": E, "F": F,
                "alpha": str(alpha),
                "beta": str(beta),
                "seed_hits": seed_hits,
                "n_seeds": len(seed_hits),
                "n_distinct_k": n_k,
            }
            # skip if it's secretly the envelope (β/α depends only on one param ratio)
            if n_k >= 2:
                hits.append(rec)
                print(f"    *** biparam multi-k pure-even A..F={(A,B,C,D,E,F)} n_k={n_k}", flush=True)
            elif info.get("ok") and tested <= 5:
                pass

    # Explicit envelope as biparam (control): k=v, m=u
    alpha_e = sp.together(256 * u**2 - 3125 * v**4 / 256)
    beta_e = sp.together(v * alpha_e)
    info_e = is_square_poly(sp.together(disc_ab(alpha_e, beta_e)), (u, v))
    # clear denominators for Z: use α = 256² u² - 3125 v⁴? identity uses α_Q

    return {
        "tested": tested,
        "multi_k_hits": hits,
        "n_multi_k": len(hits),
        "envelope_control_disc_square": info_e.get("ok"),
        "note": (
            "Low-degree biparam polynomial ansätze rarely have disc identically square "
            "except forms related to the k-ray envelope (or degenerate)."
        ),
    }


def geometric_candidates_survey() -> dict:
    """
    Catalogue of known non-rigid A5 constructions and whether they can be
    pure-even multi-k. Honest: no concrete Hurwitz candidate in hand.
    """
    return {
        "candidates_considered": [
            {
                "name": "LSW (Lavallee–Spearman–Williams)",
                "equation": "x^5+(t^2-3125)x-4(t^2-3125)",
                "pure_even": True,
                "dim": 1,
                "fixed_k_only": "-4",
                "hurwitz": "arithmetic family with A5 specialisations; not a full Hurwitz moduli description",
                "multi_k": False,
            },
            {
                "name": "Fixed-k pure-even slices (enlarged catalogue)",
                "equation": "α=256m²-3125k⁴/256, β=kα",
                "pure_even": True,
                "dim": 1,
                "fixed_k_only": "one k each",
                "hurwitz": False,
                "multi_k": False,
            },
            {
                "name": "Arithmetic multi-k envelope",
                "equation": "α=256m²-3125s⁴/256, β=s·α over Q(m,s)",
                "pure_even": True,
                "dim": 2,
                "fixed_k_only": "all k by freezing s",
                "hurwitz": False,
                "multi_k": True,
                "status": "arithmetic candidate only — not geometric Hurwitz",
            },
            {
                "name": "Homogenisation rays",
                "equation": "x^5+α0 t^4 x+β0 t^5",
                "pure_even": True,
                "dim": 1,
                "multi_k": False,
                "hurwitz": False,
            },
            {
                "name": "Rigid φ (3A,3A,5A)",
                "equation": "6y^5-15y^4+10y^3",
                "pure_even": False,
                "dim": 0,
                "status": "ABANDONED for fusion over Q (disc=5·□)",
            },
            {
                "name": "Mestre f-tr / Arala lines",
                "pure_even": "rare in scans",
                "multi_k": False,
                "status": "no multi-seed pure-even hit in prior bounds",
            },
            {
                "name": "Positive-dim A5 Hurwitz space (abstract)",
                "pure_even": "unknown",
                "multi_k": "unknown",
                "status": "NO CONCRETE EQUATION candidate that is pure-even and hits several fixed-k families",
                "hurwitz": True,
            },
        ],
        "concrete_hurwitz_pure_even_multi_k": None,
    }


# ---------------------------------------------------------------------------
# D. Obstruction notes: why cross-k is hard
# ---------------------------------------------------------------------------
def obstruction_notes() -> dict:
    """
    The BJ even surface S: γ² = 256 α⁵ + 3125 β⁴.
    Fixed-k pure-even loci are rational curves C_k ⊂ S.
    A multi-k pure-even family is a curve C ⊂ S not contained in any C_k
    that meets several C_k (or a surface containing several C_k).

    The arithmetic envelope is the surface ruled by the C_k's
    (union of all C_k is dense in a surface component?).

    Actually: points with β=kα and α+3125k⁴/256=□ cover a 2-dim
    family (params m,k). That surface component IS the envelope.
    Cross-k rational curves would be curves on that surface not equal to a C_k.

    E.g. m = const, k = u: α(u)=256 m0² - 3125 u⁴/256, β=u α
    — this is pure-even (disc=(256 α² m0)²) and k varies!
    So it IS a pure-even curve crossing all k.

    Check: at fixed m=m0≠0, vary k=u:
    α(u) = 256 m0² - 3125 u⁴/256, β = u α(u)
    disc = (256 α² m0)² square for all u. YES!

    This curve passes through every k-slice at the point with that m0.
    For m0=5: flagship has m=5 on k=-8/5; LSW uses m=55 etc — different m.
    So a single fixed-m curve hits each C_k at one point (the m=m0 point on that ray),
    but catalogue seeds generally have different m on different k, so one fixed-m
    curve may miss most catalogue seeds.

    Multi-seed from catalogue requires a curve through points with different (k,m).
    """
    print("  D. fixed-m cross-k curves (pure-even by construction)...", flush=True)
    # Fixed m0 curve: pure-even, varies k
    m0 = sp.Integer(5)
    alpha = sp.together(256 * m0**2 - 3125 * t**4 / 256)
    beta = sp.together(t * alpha)
    D = disc_ab(alpha, beta)
    expected = sp.expand((256 * alpha**2 * m0) ** 2)
    ok = sp.expand(sp.together(D - expected)) == 0

    # Which catalogue seeds lie on m=5 curve?
    # Need m_seed = 5 and α matches
    on_curve = []
    for k_str, seeds in MULTI_K_SEEDS.items():
        k0 = Fraction(k_str)
        for tag, a0, b0 in seeds:
            # m^2 = (a0 + 3125 k0^4/256)/256
            rhs = Fraction(a0) + Fraction(3125) * (k0**4) / 256
            m2 = rhs / 256
            if m2 == Fraction(25):  # m=±5
                on_curve.append({"tag": tag, "k": k_str, "a": a0, "b": b0, "m": "±5"})

    # Integer form for m0=5: α = 256*25 - 3125 t^4/256 = 6400 - 3125 t^4/256
    # α = (6400*256 - 3125 t^4)/256 = (1638400 - 3125 t^4)/256
    # For t=p/q...

    # Also m=55 curve (LSW classical param)
    m55_on = []
    for k_str, seeds in MULTI_K_SEEDS.items():
        k0 = Fraction(k_str)
        for tag, a0, b0 in seeds:
            rhs = Fraction(a0) + Fraction(3125) * (k0**4) / 256
            m2 = rhs / 256
            if m2 == Fraction(55**2):
                m55_on.append({"tag": tag, "k": k_str, "a": a0, "b": b0})

    return {
        "fixed_m_curve_disc_square": ok,
        "fixed_m": 5,
        "alpha": str(alpha),
        "beta": str(beta),
        "catalogue_seeds_on_m5": on_curve,
        "catalogue_seeds_on_m55": m55_on,
        "interpretation": (
            "Curves of fixed m and varying k are pure-even and cross all k-slices, "
            "but each meets C_k at only the single point with that m. Catalogue seeds "
            "sit at different m on different k, so fixed-m curves recover at most one "
            "seed per k and only when that seed's m matches. "
            "The 2-param envelope (m,k both free) recovers all; any 1-param curve "
            "through two catalogue seeds from different k with different m must be "
            "a non-constant path in (m,k)-space — i.e. a curve m(u), k(u) with "
            "disc still identically square (always true on the envelope)."
        ),
        "constructive_1param_through_two_seeds": (
            "Given seeds (αi,βi) on the envelope with parameters (mi,ki), i=1,2: "
            "any rational path (m(u),k(u)) with (m(0),k(0))=(m1,k1), (m(1),k(1))=(m2,k2) "
            "gives a pure-even 1-param family through both seeds via the envelope formulas. "
            "Example: linear path m=(1-u)m1+u m2, k=(1-u)k1+u k2."
        ),
    }


def linear_path_family(seed1, seed2) -> dict:
    """
    Construct explicit pure-even 1-param family through two envelope seeds
    via linear path in (m,k)-space.
    """
    tag1, a1, b1, k1s = seed1
    tag2, a2, b2, k2s = seed2
    k1, k2 = Fraction(k1s), Fraction(k2s)

    def m_of(a, k):
        rhs = Fraction(a) + Fraction(3125) * (k**4) / 256
        m2 = rhs / 256
        num, den = m2.numerator, m2.denominator
        rn = int(sp.integer_nthroot(abs(num), 2)[0])
        rd = int(sp.integer_nthroot(den, 2)[0])
        if not (sp.integer_nthroot(abs(num), 2)[1] and sp.integer_nthroot(den, 2)[1]):
            return None
        return Fraction(rn, rd) if num >= 0 else Fraction(-rn, rd)

    m1, m2 = m_of(a1, k1), m_of(a2, k2)
    if m1 is None or m2 is None:
        return {"ok": False, "reason": "seed_not_on_envelope"}

    # path
    mu = (1 - u) * m1 + u * m2
    ku = (1 - u) * k1 + u * k2
    alpha = sp.together(256 * mu**2 - 3125 * ku**4 / 256)
    beta = sp.together(ku * alpha)
    D = disc_ab(alpha, beta)
    expected = sp.expand((256 * alpha**2 * mu) ** 2)
    ok = sp.expand(sp.together(D - expected)) == 0

    # endpoints
    a_u0 = sp.simplify(alpha.subs(u, 0))
    b_u0 = sp.simplify(beta.subs(u, 0))
    a_u1 = sp.simplify(alpha.subs(u, 1))
    b_u1 = sp.simplify(beta.subs(u, 1))

    # Sample a midpoint for Gal
    mid = None
    try:
        am = sp.Rational(sp.simplify(alpha.subs(u, sp.Rational(1, 2))))
        bm = sp.Rational(sp.simplify(beta.subs(u, sp.Rational(1, 2))))
        if am.denominator == 1 and bm.denominator == 1 and am != 0:
            r = classify_poly(x**5 + int(am) * x + int(bm), do_galois=True)
            mid = {"u": "1/2", "alpha": int(am), "beta": int(bm), "status": r.get("status"), "gal": r.get("galois")}
    except Exception as e:
        mid = {"error": str(e)}

    return {
        "ok": True,
        "seed1": {"tag": tag1, "a": a1, "b": b1, "k": k1s, "m": str(m1)},
        "seed2": {"tag": tag2, "a": a2, "b": b2, "k": k2s, "m": str(m2)},
        "m_path": str(mu),
        "k_path": str(ku),
        "alpha": str(alpha),
        "beta": str(beta),
        "disc_identically_square": ok,
        "endpoint_0": (str(a_u0), str(b_u0)),
        "endpoint_1": (str(a_u1), str(b_u1)),
        "endpoints_match": (
            sp.simplify(a_u0 - a1) == 0
            and sp.simplify(b_u0 - b1) == 0
            and sp.simplify(a_u1 - a2) == 0
            and sp.simplify(b_u1 - b2) == 0
        ),
        "midpoint_sample": mid,
        "is_hurwitz_geometric": False,
        "note": (
            "Pure-even 1-param family through two different-k seeds via envelope path. "
            "Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family "
            "over P1 is A5-or-smaller for specialisations; no branch-cycle description."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("NON-RIGID HURWITZ SEARCH — abandon φ; multi-k pure-even geometric family", flush=True)

    env = arithmetic_envelope()
    print(f"  envelope disc□={env['disc_identically_square']} seeds {env['seeds_on_envelope']}", flush=True)

    cross = cross_k_curve_search()
    print(f"  cross-k curve hits (non-ray): {cross['n_hits']} tested={cross['tested']}", flush=True)

    mest = mestre_two_seed_search()
    print(f"  Mestre multi-seed pure-even: {mest['n_multi']} (tested {mest['tested']})", flush=True)

    bip = biparam_bj_ansatz()
    print(f"  biparam multi-k hits: {bip['n_multi_k']} (tested {bip['tested']})", flush=True)

    obs = obstruction_notes()
    print(f"  fixed-m=5 catalogue hits: {obs['catalogue_seeds_on_m5']}", flush=True)

    # Constructive paths through important pairs
    paths = []
    path_specs = [
        (("flagship", -55, 88, "-8/5"), ("lsw_m100", -100, 400, "-4")),
        (("flagship", -55, 88, "-8/5"), ("classical", 20, 16, "4/5")),
        (("flagship", -55, 88, "-8/5"), ("s180", -180, 432, "-12/5")),
        (("classical", 20, 16, "4/5"), ("lsw_m100", -100, 400, "-4")),
        (("flag_145", 145, -232, "-8/5"), ("lsw_124m", 124, -496, "-4")),
    ]
    print("  constructive envelope paths through cross-k seeds...", flush=True)
    for s1, s2 in path_specs:
        rec = linear_path_family(s1, s2)
        paths.append(rec)
        if rec.get("ok"):
            print(
                f"    path {s1[0]}→{s2[0]}: disc□={rec['disc_identically_square']} "
                f"endpoints_match={rec['endpoints_match']} mid={rec.get('midpoint_sample')}",
                flush=True,
            )

    survey = geometric_candidates_survey()

    # Verdict
    n_cross = cross["n_hits"]
    n_mest = mest["n_multi"]
    n_bip = bip["n_multi_k"]
    n_paths_ok = sum(1 for p in paths if p.get("ok") and p.get("disc_identically_square"))

    concrete_hurwitz = None  # still none

    verdict = (
        f"Arithmetic multi-k envelope: pure-even over Q(m,s), recovers {env['seeds_on_envelope']} "
        f"catalogue seeds on fixed-k slices — NOT Hurwitz-geometric. "
        f"Cross-k rational curves (monomial/Bezier/quad) non-ray pure-even hits: {n_cross}. "
        f"Mestre multi-seed pure-even: {n_mest}. Biparam poly multi-k: {n_bip}. "
        f"Constructive envelope paths through cross-k seed pairs: {n_paths_ok}/{len(paths)} "
        f"(arithmetic pure-even 1-param through two k's). "
        f"Concrete positive-dim Hurwitz pure-even multi-k candidate: {concrete_hurwitz}. "
        "AMBITIOUS ROUTE: no geometric Hurwitz candidate yet; arithmetic envelope + paths "
        "give pure-even multi-k families without branch-cycle geometry."
    )

    elapsed = round(time.time() - t0, 2)

    lines = [
        r"# Non-rigid pure-even multi-\(k\) search (abandon \(\varphi\))",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Goal",
        "",
        r"Find a **non-rigid geometric family** (positive-dimensional Hurwitz space) that is",
        r"**already pure-even** and **specialises onto several fixed-\(k\) arithmetic families**",
        r"(LSW \(k=-4\), flagship \(k=-8/5\), classical \(k=4/5\), …).",
        "",
        r"**Abandon** rigid \(\varphi\) (disc \(=5\cdot\square\) over \(\mathbb{Q}\)).",
        "",
        r"This is the most ambitious fusion route. Prior status: **no concrete Hurwitz candidate**.",
        "",
        "---",
        "",
        r"## A. Arithmetic multi-\(k\) envelope (not Hurwitz)",
        "",
        r"Over \(\mathbb{Q}(m,s)\):",
        "",
        r"$$\alpha(m,s)=256 m^2-\frac{3125\, s^4}{256},\qquad \beta(m,s)=s\cdot\alpha(m,s)$$",
        "",
        r"$$\operatorname{disc}=(256\,\alpha^2 m)^2\quad\text{(identically square)}.$$",
        "",
        f"- Disc identity: **{env['disc_identically_square']}**",
        f"- Catalogue seeds recovered: **{env['seeds_on_envelope']}**",
        f"- Hurwitz-geometric? **{env['is_hurwitz_geometric']}**",
        f"- {env['note']}",
        "",
        "### Recovery by freezing \(s=k\)",
        "",
    ]
    for r in env["seed_recoveries"]:
        lines.append(
            f"- k={r['k']}: {r['n_on_envelope']}/{r['n_seeds']} seeds on envelope"
        )

    lines += [
        "",
        "### Sample specialisations",
        "",
    ]
    for g in env["gal_samples"]:
        lines.append(f"- m={g.get('m')} s={g.get('s')}: α={g.get('alpha')} β={g.get('beta')} → {g.get('status')} {g.get('gal')}")

    lines += [
        "",
        "---",
        "",
        r"## B. Cross-\(k\) pure-even rational curves (ansatz search)",
        "",
        f"- Tested curves: **{cross['tested']}**",
        f"- Pure-even hits **not** contained in a single \(k\)-ray: **{cross['n_hits']}**",
        f"- Priority seed pairs: {cross['priority_pairs']}",
        "",
    ]
    if cross["n_hits"] == 0:
        lines.append(
            r"_No non-ray pure-even monomial/Bezier/quad bridge through two different-\(k\) "
            r"catalogue seeds in the scanned ansätze. (Consistent with earlier even-surface scan.)_"
        )
    for h in cross["cross_k_pure_even_hits"][:10]:
        lines.append(f"- HIT: `{h}`")

    lines += [
        "",
        "---",
        "",
        r"## C. Hurwitz-adjacent / deformation ansätze",
        "",
        r"### C1. Mestre \(f-tr\)",
        "",
        f"- Tested: {mest['tested']}",
        f"- Multi-seed pure-even hits: **{mest['n_multi']}**",
        "",
    ]
    for h in mest["hits"][:8]:
        lines.append(f"- `{h}`")

    lines += [
        "",
        r"### C2. Low-degree biparameter BJ",
        "",
        f"- Tested: {bip['tested']}",
        f"- Multi-\(k\) pure-even hits: **{bip['n_multi_k']}**",
        f"- Envelope control disc□: {bip['envelope_control_disc_square']}",
        f"- {bip['note']}",
        "",
    ]
    for h in bip["multi_k_hits"][:6]:
        lines.append(
            f"- A..F=({h['A']},{h['B']},{h['C']},{h['D']},{h['E']},{h['F']}) "
            f"n_k={h['n_distinct_k']} seeds={h['n_seeds']}"
        )

    lines += [
        "",
        "### C3. Survey of constructions",
        "",
        "| Name | Pure-even? | Multi-\(k\)? | Hurwitz? | Status |",
        "|------|:----------:|:------------:|:--------:|--------|",
    ]
    for c in survey["candidates_considered"]:
        lines.append(
            f"| {c['name']} | {c.get('pure_even')} | {c.get('multi_k')} | "
            f"{c.get('hurwitz')} | {c.get('status', c.get('fixed_k_only', ''))} |"
        )

    lines += [
        "",
        f"**Concrete Hurwitz pure-even multi-\(k\) candidate:** `{survey['concrete_hurwitz_pure_even_multi_k']}`",
        "",
        "---",
        "",
        r"## D. Structure of the pure-even locus (why paths work arithmetically)",
        "",
        f"- Fixed-\(m\) varying-\(k\) curve disc□: **{obs['fixed_m_curve_disc_square']}**",
        f"- Catalogue seeds on m=5 curve: `{obs['catalogue_seeds_on_m5']}`",
        f"- Catalogue seeds on m=55 curve: `{obs['catalogue_seeds_on_m55']}`",
        "",
        obs["interpretation"],
        "",
        obs["constructive_1param_through_two_seeds"],
        "",
        r"### Explicit pure-even paths through cross-\(k\) seed pairs",
        "",
    ]
    for p in paths:
        if not p.get("ok"):
            lines.append(f"- FAIL {p}")
            continue
        lines.append(
            f"#### {p['seed1']['tag']} (\(k={p['seed1']['k']}\)) → "
            f"{p['seed2']['tag']} (\(k={p['seed2']['k']}\))"
        )
        lines.append(f"- m-path: `{p['m_path']}`")
        lines.append(f"- k-path: `{p['k_path']}`")
        lines.append(f"- disc identically square: **{p['disc_identically_square']}**")
        lines.append(f"- endpoints match seeds: **{p['endpoints_match']}**")
        lines.append(f"- midpoint sample: `{p['midpoint_sample']}`")
        lines.append(f"- Hurwitz-geometric? **{p['is_hurwitz_geometric']}**")
        lines.append(f"- {p['note']}")
        lines.append("")

    lines += [
        "---",
        "",
        "## E. Conclusions",
        "",
        r"1. **Abandoned \(\varphi\)** for this route: rigid, not pure-even over \(\mathbb{Q}\).",
        "",
        r"2. **Arithmetic multi-\(k\) envelope** over \(\mathbb{Q}(m,s)\) is pure-even and "
        r"specialises onto **all** fixed-\(k\) slices (and essentially all catalogue seeds on them). "
        r"It is **not** a Hurwitz space: no independent branch-cycle / geometric monodromy package.",
        "",
        r"3. **Constructive pure-even 1-param paths** through any two envelope seeds "
        r"(including different \(k\)) exist by linear paths in \((m,k)\)-space. "
        r"Example: flagship \(\leftrightarrow\) LSW, flagship \(\leftrightarrow\) classical. "
        r"These solve the *arithmetic* multi-seed pure-even problem across \(k\).",
        "",
        r"4. **No concrete positive-dimensional Hurwitz candidate** was found that is "
        r"already pure-even and maps onto several fixed-\(k\) families with geometric "
        r"branch data. Cross-\(k\) polynomial curve ansätze on the even surface and "
        r"Mestre/biparam searches produced **no** non-envelope geometric hit in bounds.",
        "",
        r"5. **Programme split:**",
        r"   - *Arithmetic fusion fuel:* envelope + paths (pure-even multi-\(k\) over \(\mathbb{Q}\)).",
        r"   - *Geometric fusion (ambitious):* still **open / no candidate** — would need a "
        r"true family of covers with dim\(>0\) Hurwitz data whose BJ/Hilbert specialisations "
        r"land on several \(C_k\).",
        "",
        r"### Recommended stance",
        "",
        r"- Treat the **\((m,s)\)-envelope** and **cross-\(k\) envelope paths** as the "
        r"explicit pure-even multi-family arithmetic object.",
        r"- Do **not** claim Hurwitz geometry for them.",
        r"- Further geometric work must start from known positive-dim \(A_5\) Hurwitz strata "
        r"(or other non-rigid constructions) and *test* pure-even + multi-\(k\) specialisation — "
        r"not from \(\varphi\).",
        "",
        "_Generated by nonrigid_hurwitz_search.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "envelope": env,
        "cross_k": {
            "tested": cross["tested"],
            "n_hits": cross["n_hits"],
            "hits": cross["cross_k_pure_even_hits"],
        },
        "mestre": mest,
        "biparam": {
            "tested": bip["tested"],
            "n_multi_k": bip["n_multi_k"],
            "hits": bip["multi_k_hits"],
        },
        "obstruction": obs,
        "paths": paths,
        "survey": survey,
        "concrete_hurwitz_candidate": None,
    }

    write_md(OUT / "NONRIGID_HURWITZ_SEARCH.md", doc)
    write_md(RESULTS / "NONRIGID_HURWITZ_SEARCH.md", doc)
    write_md(ROOT / "NONRIGID_HURWITZ_SEARCH.md", doc)
    write_json(OUT / "NONRIGID_HURWITZ_SEARCH.json", blob)
    print(verdict, flush=True)
    print(f"Wrote NONRIGID_HURWITZ_SEARCH.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
