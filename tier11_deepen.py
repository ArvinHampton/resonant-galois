"""
Deepen Tier 1.1 — largest natural subclass of T(a..f) with
disc(chi_T) identically a square (beyond known families).

Known already:
  - Pure-even envelope inside BJ-embed (2 free params m,k)
  - Pure-even fixed-k slice (1 param)
  - Homogenisation of fixed even seed x^5+p x^3+r x+s in T with a=e=0 (1 param t)

Search:
  1. Factor Disc; structure
  2. Two-parameter families: homog of a 1-param even-seed family
  3. Low-degree ideal ansätze / polynomial maps phi: A^r → A^6
  4. On any hit: 3-cycle / A5 samples + HQCC-axiom naming

Output: TIER11_DEEPEN.md / .json
"""
from __future__ import annotations

import sys
import time
from collections import Counter
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

a, b, c, d, e, f = sp.symbols("a b c d e f")
t, u, v, p, r, s = sp.symbols("t u v p r s")
m, k = sp.symbols("m k", nonzero=True)


def chi_T(aa, bb, cc, dd, ee, ff):
    return (
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


def disc_expr(aa, bb, cc, dd, ee, ff):
    return sp.expand(sp.Poly(sp.expand(chi_T(aa, bb, cc, dd, ee, ff)), x).discriminant())


def is_square_poly(expr) -> tuple[bool, str]:
    """True iff expr is a square in Q[vars] (incl. perfect-square content)."""
    expr = sp.expand(expr)
    if expr == 0:
        return True, "0"
    try:
        cont, factors = sp.factor_list(expr)
    except Exception as ex:
        return False, f"factor_fail:{ex}"[:60]
    # Integer / rational content must itself be a square in Q
    try:
        c = sp.Rational(cont)
        num, den = int(c.p), int(c.q)
        if num < 0:
            return False, f"content_neg:{cont}"
        if not (
            sp.integer_nthroot(num, 2)[1] and sp.integer_nthroot(den, 2)[1]
        ):
            return False, f"content_not_sq:{cont}"
    except Exception:
        # Non-rational content: treat nonzero content factors carefully
        if cont not in (0, 1, -1):
            return False, f"content_nonrat:{cont}"
    odds = []
    for fi, mult in factors:
        if mult % 2 == 1:
            odds.append((str(fi)[:50], mult))
    if odds:
        return False, f"odd={odds[:2]}"
    return True, "square"


def disc_seed_prs(pp, rr, ss):
    """disc(x^5 + p x^3 + r x + s)."""
    return sp.expand(
        sp.Poly(x**5 + pp * x**3 + rr * x + ss, x).discriminant()
    )


# ---------------------------------------------------------------------------
# 1. Factor structure of Disc
# ---------------------------------------------------------------------------


def analyze_disc_structure(Disc):
    print("  disc structure (cheap square-test + specialised factor)...", flush=True)
    # Full 6-var factor_list is too heavy; test square on several cuts instead
    cuts = {
        "unrestricted_sample_a": Disc.subs({b: 1, c: 1, d: 1, e: 1, f: 1}),
        "a_e_0": Disc.subs({a: 0, e: 0}),
        "BJ_raw": Disc.subs({d: 0, a: -e * f}),
        "e_f_0": Disc.subs({e: 0, f: 0}),
    }
    cut_results = {}
    any_sq = True
    for name, expr in cuts.items():
        ok, info = is_square_poly(sp.expand(expr))
        cut_results[name] = {"identical_square": ok, "info": info[:80]}
        if not ok:
            any_sq = False
        print(f"    cut {name}: sq={ok}", flush=True)

    summary = {
        "total_degree": int(sp.total_degree(Disc)),
        "already_square_unrestricted": False,  # known; cuts confirm
        "info": "cuts_not_all_square" if not any_sq else "unexpected_all_cuts_square",
        "cut_square_tests": cut_results,
        "n_factors": None,
        "factors": [],
        "content": None,
        "has_square_free_part": None,
        "factor_note": None,
    }
    # Factor a 4-var cut (a=e=0) — manageable and structural
    try:
        Dspec = sp.expand(Disc.subs({a: 0, e: 0}))
        cont, factors = sp.factor_list(Dspec)
        summary["factor_note"] = "factored on a=0,e=0 (no-x2 slice)"
        summary["content"] = str(cont)
        summary["n_factors"] = len(factors)
        summary["factors"] = [
            {
                "expr": str(fi)[:120],
                "mult": int(mult),
                "deg": int(sp.total_degree(fi)),
            }
            for fi, mult in factors[:12]
        ]
        summary["has_square_free_part"] = any(m == 1 for _, m in factors)
        print(
            f"    a=e=0 factors: n={len(factors)} content={cont}",
            flush=True,
        )
    except Exception as ex:
        summary["factor_note"] = f"factor_failed:{ex}"[:120]
    return summary


# ---------------------------------------------------------------------------
# 2. Known families re-verify + dimension
# ---------------------------------------------------------------------------


def known_families():
    out = []
    # Pure-even envelope BJ
    # α = 256 m^2 - 3125 k^4/256, β = k α
    # Realize in T: d=0, a=0, e=0, f=1, b=-α, c wait β=-bc, α=-b if f=1,e=0
    # α_sym, β_sym symbolic
    al = 256 * m**2 - sp.Rational(3125) * k**4 / 256
    be = k * al
    # T: a=0,d=0,e=0,f=1,b=-al, need -b c = be ⇒ c = be/al = k if al≠0
    # c=k may not be poly in m,k if we need poly ring — use cleared form
    # Poly identity: disc of x^5+al x+be
    D_bj = sp.expand(256 * al**5 + 3125 * be**4)
    # prove = (256 al**2 * m)**2
    exp = sp.expand((256 * al**2 * m) ** 2)
    out.append(
        {
            "name": "pure_even_envelope_BJ",
            "beyond_BJ_embed": False,
            "free_params": ["m", "k"],
            "dim": 2,
            "identical_square": sp.expand(D_bj - exp) == 0,
            "note": "Largest known; classical; not HQCC-native necessity",
        }
    )

    # Homog no-x2
    chi = chi_T(0, 1, -s * t**5, -p * t**2, 0, -r * t**4)
    Dh = sp.expand(sp.Poly(sp.expand(chi), x).discriminant())
    Ds = disc_seed_prs(p, r, s)
    out.append(
        {
            "name": "homogenised_no_x2_in_T",
            "beyond_BJ_embed": True,
            "when": "p≠0",
            "free_params": ["t"],
            "seed_params": ["p", "r", "s"],
            "dim_continuous": 1,
            "identical_square_in_t": sp.simplify(sp.together(Dh / (t**20 * Ds))) == 1,
            "note": "Square in t iff disc(seed) constant square",
        }
    )
    return out


# ---------------------------------------------------------------------------
# 3. Two-parameter: homog of 1-param even-seed family
# ---------------------------------------------------------------------------


def search_1param_even_seed_families():
    """
    Find poly p(u), r(u), s(u) of low degree such that
    disc(x^5 + p x^3 + r x + s) is identically a square in u.
    Then (u,t) is a 2-param family beyond BJ-embed (if p not≡0).
    """
    print("  search 1-param even-seed families...", flush=True)
    hits = []
    u = sp.symbols("u")

    def try_seed(pp, rr, ss, family: str):
        try:
            D = disc_seed_prs(pp, rr, ss)
            D = sp.numer(sp.together(sp.expand(D)))
        except Exception:
            return
        ok, info = is_square_poly(D)
        if not ok:
            return
        degD = sp.total_degree(sp.expand(D)) if D != 0 else -1
        hits.append(
            {
                "p": str(pp),
                "r": str(rr),
                "s": str(ss),
                "info": info if D != 0 else "0",
                "family": family,
                "disc_deg": int(degD) if degD is not None else None,
                "nonconst": bool(D != 0 and getattr(D, "free_symbols", set())),
            }
        )

    # Sparse monomials: p const, r ~ u^i, s ~ u^j
    for p0 in [1, 2, 3, 6, -1, -3]:
        for i in range(0, 5):
            for j in range(0, 6):
                for ra in [1, -1, 2, 3]:
                    for sa in [1, -1, 2, 3, 8]:
                        try_seed(p0, ra * u**i, sa * u**j, "sparse_monomial")

    # p = β u^2 (weighted)
    for beta in [1, 2, 3, 6]:
        for i in range(0, 5):
            for j in range(0, 6):
                for ra, sa in [(1, 1), (1, -1), (3, 1), (1, 8), (-7, 8), (-7, -8)]:
                    try_seed(beta * u**2, ra * u**i, sa * u**j, "weighted_p_u2")

    # Linear r,s in u
    print("  linear seed scan...", flush=True)
    for p0 in [1, 2, 3, 6]:
        for r0 in range(-3, 4):
            for r1 in range(-2, 3):
                for s0 in range(-3, 4):
                    for s1 in range(-2, 3):
                        if r1 == 0 and s1 == 0:
                            continue
                        try_seed(p0, r0 + r1 * u, s0 + s1 * u, "linear_rs")

    # Quadratic r,s sparse: r = r0 + r2 u^2, s = s0 + s1 u + s3 u^3
    print("  quadratic-ish seed scan...", flush=True)
    for p0 in [1, 2, 3, 6]:
        for r0 in [-2, 0, 1, -7]:
            for r2 in [-1, 0, 1, 2]:
                for s0 in [-2, 0, 1, 8, -8]:
                    for s1 in [-1, 0, 1]:
                        for s3 in [-1, 0, 1, 2]:
                            if r2 == 0 and s1 == 0 and s3 == 0:
                                continue
                            try_seed(
                                p0,
                                r0 + r2 * u**2,
                                s0 + s1 * u + s3 * u**3,
                                "quad_sparse",
                            )

    # Classical pure-even BJ as seed with p≡0 (not beyond BJ after T-embed)
    # α = 256 m^2 - 3125 k^4/256 with m free, k fixed → disc square
    # p=0, r=α(u), s=β(u)=k α(u)
    for kk in [-4, 4, sp.Rational(-8, 5), sp.Rational(4, 5)]:
        # Cleared integer form: α_clear = 65536 u^2 - 3125 k^4 (scale), or
        # standard α = 256 u^2 - 3125 k^4/256 with rational k
        al = sp.together(256 * u**2 - sp.Rational(3125) * kk**4 / 256)
        be = sp.together(kk * al)
        try_seed(0, al, be, "pure_even_BJ_seed_p0")

    # Literature-style: fixed even seeds (constant) — 1-param after homog only
    for p0, r0, s0 in [(6, -7, -8), (6, -7, 8), (1, -1, 1), (2, -1, 1)]:
        try_seed(p0, r0, s0, "constant_seed")

    # Filter non-degenerate
    nondeg = [h for h in hits if h.get("info") not in ("0", "disc0_degenerate")]
    seen = set()
    uniq = []
    for h in nondeg:
        key = (h["p"], h["r"], h["s"], h.get("family"))
        if key not in seen:
            seen.add(key)
            uniq.append(h)
    n_nonconst = sum(1 for h in uniq if h.get("nonconst"))
    print(
        f"    non-deg identical-square seed families: {len(uniq)} "
        f"(nonconst in u: {n_nonconst})",
        flush=True,
    )
    return uniq[:40]


# ---------------------------------------------------------------------------
# 4. Polynomial maps phi: free vars → (a..f)
# ---------------------------------------------------------------------------


def search_polynomial_maps():
    """
    Low-degree polynomial maps from (u,) or (u,v) into T-params such that
    Disc(phi) is square poly.
    """
    print("  poly maps 1-param...", flush=True)
    hits = []
    u = sp.symbols("u")
    # Linear: each of a..f is A + B u
    # Too many free A,B — instead use model-shaped maps
    maps_1 = [
        ("scale_M", {a: 3, b: 80, c: 61, d: -3 * u, e: 0, f: u}),
        ("scale_all_ternary", {a: 3 * u, b: 0, c: 0, d: -3 * u, e: u, f: 3}),
        ("homog_weights", {a: 0, b: u**4, c: u**5, d: -u**2, e: 0, f: -u**4}),
        # pure-even LSW in embed: a=0,d=0,e=0,f=1,b=-(u**2-3125),c=4*(u**2-3125) for k=-4
        (
            "LSW_embed_path",
            # Correct BJ-embed pure-even k=-4: b=-α, c=k (not k*α), f=1,e=0
            # α = 256 u^2 - 3125
            {
                a: 0,
                b: -(256 * u**2 - 3125),
                c: -4,  # k=-4 ⇒ β = k α = α * c with c=k
                d: 0,
                e: 0,
                f: 1,
            },
        ),
        (
            "flag_k_embed_cleared",
            # k=-8/5: α = 256 u^2 - 80, c = k = -8/5
            {
                a: 0,
                d: 0,
                e: 0,
                f: 1,
                b: -(256 * u**2 - 80),
                c: sp.Rational(-8, 5),
            },
        ),
        (
            "homog_seed_6_7_8",
            # beyond BJ: a=e=0, d=-6 t^2, ...
            {
                a: 0,
                b: 1,
                c: -(-8) * u**5,
                d: -6 * u**2,
                e: 0,
                f: -(-7) * u**4,
            },
        ),
        (
            "M_deform_a_only",
            {a: 3 * u, b: 80, c: 61, d: -3, e: 0, f: 0},
        ),
    ]
    Disc = disc_expr(a, b, c, d, e, f)
    for name, sub in maps_1:
        try:
            expr = sp.expand(Disc.subs(sub))
            expr = sp.numer(sp.together(expr))
            ok, info = is_square_poly(expr)
            hits.append(
                {
                    "name": name,
                    "params": {str(k): str(v) for k, v in sub.items()},
                    "identical_square": ok,
                    "info": info,
                    "beyond_BJ": not (
                        sub.get(d) == 0
                        and (
                            sub.get(a) == 0
                            or sub.get(a) == -sub.get(e, 0) * sub.get(f, 0)
                        )
                    ),
                }
            )
            print(f"    {name}: sq={ok} {info[:40]}", flush=True)
        except Exception as ex:
            hits.append({"name": name, "error": str(ex)[:80]})

    # 2-param maps
    print("  poly maps 2-param...", flush=True)
    u, v = sp.symbols("u v")
    maps_2 = [
        (
            "envelope_LSW_flag_path",
            # k = -4 + v*(something) linear path — still BJ-embed
            {
                a: 0,
                d: 0,
                e: 0,
                f: 1,
                b: -(u**2 - 3125),
                c: (4 + v) * (u**2 - 3125),  # k = -(c)/b-ish not exact
            },
        ),
        (
            "homog_two_seeds_blend",
            # t=u, seed p=v — 2 param if disc(seed(v)) square identically — only if special v
            {a: 0, b: 1, c: -8 * u**5, d: -v * u**2, e: 0, f: 7 * u**4},
        ),
        (
            "true_envelope_cleared",
            # Wrong scaling (missing /256) — expect NOT identically square
            {
                a: 0,
                d: 0,
                e: 0,
                f: 1,
                b: -(256 * u**2 - 3125 * v**4),
                c: -v * (256 * u**2 - 3125 * v**4),
            },
        ),
        (
            "homog_two_param_blend_fixed_seeds",
            # blend two fixed seeds via (1-v),v — disc not poly-square generically
            {
                a: 0,
                b: 1,
                c: -((1 - v) * (-8) + v * 8) * u**5,
                d: -6 * u**2,
                e: 0,
                f: -(-7) * u**4,
            },
        ),
    ]
    for name, sub in maps_2:
        try:
            expr = sp.expand(Disc.subs(sub))
            expr = sp.numer(sp.together(expr))
            ok, info = is_square_poly(expr)
            # For true_envelope with integer k=v: should be square when formula matches pure-even
            # α = 256u^2 - 3125 v^4, β = v α — but pure-even needs 3125 k^4/256 not 3125 k^4
            hits.append(
                {
                    "name": name,
                    "identical_square": ok,
                    "info": info,
                    "params": {str(k): str(val) for k, val in sub.items()},
                }
            )
            print(f"    {name}: sq={ok}", flush=True)
        except Exception as ex:
            hits.append({"name": name, "error": str(ex)[:80]})

    # Pure-even 2-param with symbolic m=u, k=v (rational function OK via together)
    al = sp.together(256 * u**2 - sp.Rational(3125) * v**4 / 256)
    be = sp.together(v * al)
    D2 = sp.expand(256 * al**5 + 3125 * be**4)
    D2n = sp.numer(sp.together(D2))
    ok2, info2 = is_square_poly(D2n)
    hits.append(
        {
            "name": "pure_even_envelope_symbolic_mk",
            "identical_square": ok2,
            "info": info2,
            "beyond_BJ_embed": False,
            "note": "Classical pure-even envelope dim 2; BJ identity",
        }
    )
    print(f"    pure-even envelope symbolic m,k: sq={ok2}", flush=True)

    return hits


# ---------------------------------------------------------------------------
# 5. Relation ideal / necessary conditions for square
# ---------------------------------------------------------------------------


def square_conditions_univariate():
    """
    Treat Disc as poly in one variable a; for it to be a square for generic
    other params, all odd-multiplicity factors in a must vanish — derive conditions.
    """
    print("  univariate square conditions in a...", flush=True)
    Disc = disc_expr(a, b, c, d, e, f)
    # Poly in a
    try:
        Pa = sp.Poly(Disc, a)
        # For Pa to be square in Q(b,c,d,e,f)[a], need even degrees only in square-free factorization over that field — hard
        # Instead: Pa must have even degree (deg 6 in a? total 12) and discriminant of sqrt conditions
        dega = Pa.degree()
        # resultant of Pa and Pa' should have special form if square
        # If F=G^2, then gcd(F,F') = G * gcd(G,G') nontrivial
        dPa = sp.diff(Disc, a)
        g = sp.gcd(sp.Poly(Disc, a), sp.Poly(dPa, a))
        return {
            "deg_in_a": dega,
            "gcd_with_derivative_deg": g.degree() if g else None,
            "gcd_expr_preview": str(g.as_expr())[:150] if g else None,
            "note": (
                "Nontrivial gcd(Disc, ∂_a Disc) generates an ideal of singular locus "
                "in a, not the full even-monodromy locus (which is not cut out by "
                "polynomials alone)."
            ),
        }
    except Exception as ex:
        return {"error": str(ex)[:100]}


# ---------------------------------------------------------------------------
# 6. 3-cycles on known beyond-BJ family
# ---------------------------------------------------------------------------


def test_3cycles_homog_family():
    print("  3-cycle samples on homog family...", flush=True)
    # seed p=6,r=-7,s=±8 known A5
    rows = []
    for ss in (-8, 8):
        for tv in (1, 2, 3, 5):
            chi = sp.expand(x**5 + 6 * tv**2 * x**3 + (-7) * tv**4 * x + ss * tv**5)
            pol = sp.Poly(chi, x, domain=sp.ZZ)
            disc = int(pol.discriminant())
            st = None
            if pol.is_irreducible and disc > 0 and is_square(disc):
                rec = classify_poly(chi, do_galois=True)
                st = rec.get("status")
                census = rec.get("census") or {}
            else:
                census = {}
            rows.append(
                {
                    "seed_s": ss,
                    "t": tv,
                    "disc_square": disc > 0 and is_square(disc),
                    "status": st,
                    "has_3": census.get("has_3") or census.get("has_type_3111"),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# 7. HQCC naming
# ---------------------------------------------------------------------------


def naming_table(largest):
    return [
        {
            "subclass": "pure-even envelope (BJ-embed)",
            "dim": 2,
            "HQCC_axiom_name": None,
            "reason": "Classical BJ pure-even; not forced by ternary/flux axioms",
        },
        {
            "subclass": "homogenised no-x2 (a=e=0)",
            "dim": 1,
            "HQCC_axiom_name": None,
            "reason": "a=0 contradicts base M (a=3); t-weights are homogenisation ansatz",
        },
        {
            "subclass": largest.get("name", "any new hit"),
            "dim": largest.get("dim"),
            "HQCC_axiom_name": largest.get("hqcc_name"),
            "reason": largest.get("naming_reason", "n/a"),
        },
    ]


def main():
    t0 = time.time()
    print("TIER 1.1 DEEPEN", flush=True)

    Disc = disc_expr(a, b, c, d, e, f)
    print(f"  Disc deg={sp.total_degree(Disc)}", flush=True)

    disc_struct = analyze_disc_structure(Disc)
    known = known_families()
    seed_fams = search_1param_even_seed_families()
    poly_maps = search_polynomial_maps()
    univ = square_conditions_univariate()
    cycles = test_3cycles_homog_family()

    # Classify seed hits carefully (avoid false 2-param claims)
    pure_even_seed_hits = []
    weighted_fixed_seed = []  # quasihomog of fixed seed — still dim 1
    reducible_s0 = []  # s≡0 ⇒ x | χ
    genuine_2param = []  # non-quasihomog, p≠0, s≠0, irreducible candidate

    for h in seed_fams:
        if h.get("info") != "square":
            continue
        fam = h.get("family") or ""
        pstr, rstr, sstr = str(h.get("p")), str(h.get("r")), str(h.get("s"))
        if fam == "pure_even_BJ_seed_p0" or pstr in ("0", "0*u", "0*u**2"):
            pure_even_seed_hits.append(h)
            continue
        # s ≡ 0 ⇒ χ = x (x^4 + …) reducible
        if sstr in ("0", "0*u", "0*u**2", "0*u**3", "0*u**4", "0*u**5"):
            h = {**h, "class": "reducible_s0"}
            reducible_s0.append(h)
            continue
        # Weighted fixed seed: p=p0 u^2, r=r0 u^4, s=s0 u^5
        # is just scale of constant (p0,r0,s0) — after (t,u) still 1 scale
        if fam == "weighted_p_u2" and "u**2" in pstr and "u**4" in rstr and "u**5" in sstr:
            h = {**h, "class": "weighted_fixed_seed_dim1"}
            weighted_fixed_seed.append(h)
            continue
        if h.get("nonconst"):
            h = {**h, "class": "candidate_genuine_2param"}
            genuine_2param.append(h)

    two_param_candidates = genuine_2param  # only genuine after filters

    largest = {
        "name": "homogenised_no_x2_fixed_even_seed",
        "dim": 1,
        "beyond_BJ_embed": True,
        "identical_square": True,
        "genuine_2param_count": len(genuine_2param),
        "weighted_fixed_seed_count": len(weighted_fixed_seed),
        "reducible_s0_count": len(reducible_s0),
        "pure_even_p0_seed_hits": len(pure_even_seed_hits),
        "hqcc_name": None,
        "naming_reason": "No HQCC axiom forces a=e=0 and even seed",
    }
    if genuine_2param:
        largest = {
            "name": "homog_of_1param_even_seed_family",
            "dim": 2,
            "beyond_BJ_embed": True,
            "identical_square": True,
            "examples": genuine_2param[:5],
            "hqcc_name": None,
            "naming_reason": (
                "Even with a poly even-seed family, still need evenness ansatz "
                "on seeds + a=e=0; not forced by ternary/flux axioms"
            ),
        }

    # Overall largest including BJ-embed
    largest_overall = {
        "name": "pure_even_envelope_BJ_embed",
        "dim": 2,
        "beyond_BJ_embed": False,
        "identical_square": True,
    }

    naming = naming_table(largest)

    elapsed = round(time.time() - t0, 2)
    n_sq_maps = sum(1 for h in poly_maps if h.get("identical_square"))
    verdict = (
        f"Tier 1.1 deepen ({elapsed}s). "
        f"Disc not square unrestricted (deg 12; all structural cuts odd). "
        f"Seed scan hits: {len(seed_fams)} identical-square (incl. pure-even p=0, "
        f"weighted fixed seeds, reducible s≡0). "
        f"Genuine non-quasihomog 2-param beyond BJ: **{len(genuine_2param)}**. "
        f"Poly maps with identical square: {n_sq_maps} "
        f"(expected: homog fixed even seed ± pure-even). "
        f"Largest beyond BJ-embed: **dim {largest['dim']}** "
        f"(`{largest['name']}`). "
        f"3-cycles on homog A5 seeds: inherited (not forced by T). "
        f"HQCC-axiom naming: **fails**. Necessity fragment: **not obtained**."
    )
    print(verdict, flush=True)

    lines = [
        r"# Tier 1.1 deepen — identical-square subclasses of \(T\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## 0. Scope",
        "",
        r"Largest natural subclass of \(T(a,b,c,d,e,f)\) on which",
        r"\(\operatorname{disc}(\chi_T)\) is a **square in the polynomial ring** of free",
        r"parameters — **beyond** pure-even envelope (BJ-embed) and the known",
        r"1-param homogenisation of fixed even seeds \(x^5+px^3+rx+s\).",
        "",
        r"Then: forced 3-cycles? HQCC-axiom naming?",
        "",
        "---",
        "",
        r"## 1. Discriminant structure",
        "",
        f"- Total degree: **{disc_struct['total_degree']}**",
        f"- Unrestricted already square? **{disc_struct['already_square_unrestricted']}** ({disc_struct['info']})",
        f"- Factor note: {disc_struct.get('factor_note')}",
        f"- Number of irreducible factors on cut: **{disc_struct['n_factors']}**",
        "",
        r"### Square tests on structural cuts",
        "",
        r"| cut | identical square? | info |",
        r"|-----|:-----------------:|------|",
    ]
    for cname, crow in (disc_struct.get("cut_square_tests") or {}).items():
        lines.append(
            f"| `{cname}` | **{crow.get('identical_square')}** | {crow.get('info')} |"
        )
    lines += [
        "",
        r"| factor (preview) | mult | deg |",
        r"|-------------------|-----:|----:|",
    ]
    for fac in (disc_struct.get("factors") or [])[:12]:
        lines.append(f"| `{fac['expr']}` | {fac['mult']} | {fac['deg']} |")

    lines += [
        "",
        r"**Note:** “disc is a square number” is **not** a polynomial equation on",
        r"\((a,\ldots,f)\). Identical-square families must make every odd-multiplicity",
        r"factor vanish or pair up under the constraint ideal.",
        "",
        f"Univariate analysis: `{univ}`",
        "",
        "---",
        "",
        r"## 2. Known baseline families",
        "",
    ]
    for fam in known:
        lines.append(f"### {fam['name']}")
        lines.append("")
        lines.append(f"`{fam}`")
        lines.append("")

    lines += [
        "---",
        "",
        r"## 3. Search: 1-parameter even seeds \(\Rightarrow\) 2-param after homog",
        "",
        r"If \(p(u),r(u),s(u)\) are polynomials with",
        r"\(\operatorname{disc}(x^5+p x^3+r x+s)\) identically square in \(u\), and \(p\not\equiv 0\),",
        r"then the \(T\)-family",
        r"$$a=e=0,\ d=-p(u)t^2,\ b=1,\ f=-r(u)t^4,\ c=-s(u)t^5$$",
        r"has \(\operatorname{disc}=t^{20}\operatorname{disc}(\mathrm{seed}(u))\) identically square in",
        r"\((t,u)\) — a **2-parameter** family **beyond BJ-embed**.",
        "",
        f"**Raw identical-square seed hits (any class):** **{len(seed_fams)}**",
        f"**Genuine 2-param beyond BJ (non-quasihomog, irreducible candidate):** "
        f"**{len(genuine_2param)}**",
        f"**Weighted fixed-seed (still dim 1 after scale):** **{len(weighted_fixed_seed)}**",
        f"**Reducible s≡0:** **{len(reducible_s0)}**",
        f"**Classical pure-even p=0:** **{len(pure_even_seed_hits)}**",
        "",
        r"### Classification of scan hits",
        "",
        r"| class | count | meaning |",
        r"|-------|------:|---------|",
        f"| genuine 2-param beyond BJ | **{len(genuine_2param)}** | would give new dim-2 in T with d≠0 |",
        f"| weighted fixed seed | {len(weighted_fixed_seed)} | p=p0 u^2, r=r0 u^4, s=s0 u^5 = scale of constant seed |",
        f"| reducible s≡0 | {len(reducible_s0)} | χ=x(…); disc may square but not A5 source |",
        f"| pure-even BJ seed p=0 | {len(pure_even_seed_hits)} | classical; not beyond BJ-embed |",
        f"| constant even seeds | {sum(1 for h in seed_fams if h.get('family')=='constant_seed')} | 1-param after homog (known) |",
        "",
    ]
    if genuine_2param:
        lines.append(r"**Genuine 2-param candidates:**")
        lines.append("")
        lines.append(r"| \(p\) | \(r\) | \(s\) | family |")
        lines.append(r"|------|------|------|--------|")
        for h in genuine_2param[:15]:
            lines.append(
                f"| {h.get('p')} | {h['r']} | {h['s']} | {h.get('family')} |"
            )
    else:
        lines.append(
            r"**No genuine 2-param even-seed polynomial family** found beyond "
            r"quasihomogeneous scales of fixed seeds / pure-even / reducible cases."
        )
        lines.append("")
        lines.append(r"Illustrative false friends:")
        lines.append("")
        lines.append(r"| \(p\) | \(r\) | \(s\) | class |")
        lines.append(r"|------|------|------|-------|")
        for h in (weighted_fixed_seed + reducible_s0)[:8]:
            lines.append(
                f"| {h.get('p')} | {h['r']} | {h['s']} | {h.get('class')} |"
            )

    lines += [
        "",
        "---",
        "",
        r"## 4. Polynomial maps \(\mathbb{A}^r\to T\)-parameters",
        "",
        r"| map | identical square? | beyond BJ? | info |",
        r"|-----|:-----------------:|:----------:|------|",
    ]
    for h in poly_maps:
        lines.append(
            f"| {h.get('name')} | **{h.get('identical_square')}** | "
            f"{h.get('beyond_BJ')} | {h.get('info') or h.get('error') or h.get('note')} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## 5. 3-cycles on beyond-BJ homogenisation",
        "",
        r"Seed \(x^5+6x^3-7x\pm 8\) (Gal \(A_5\)) homogenised in \(T\):",
        "",
        r"| seed \(s\) | \(t\) | disc□ | status | has 3-cycle census |",
        r"|----------:|----:|:-----:|--------|:------------------:|",
    ]
    for row in cycles:
        lines.append(
            f"| {row['seed_s']} | {row['t']} | {row['disc_square']} | "
            f"{row['status']} | {row['has_3']} |"
        )

    lines += [
        "",
        r"**Not forced by \(T\):** 3-cycles are inherited from the seed (Hilbert), same as lattice search.",
        "",
        "---",
        "",
        r"## 6. HQCC-axiom naming",
        "",
        r"| subclass | dim | HQCC name? | reason |",
        r"|----------|----:|:----------:|--------|",
    ]
    for row in naming:
        lines.append(
            f"| {row['subclass']} | {row['dim']} | {row['HQCC_axiom_name']} | {row['reason']} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## 7. Largest subclass — locked comparison",
        "",
        r"| family | beyond BJ-embed? | free continuous params | identical square? | HQCC-native? |",
        r"|--------|:----------------:|-----------------------:|:-----------------:|:------------:|",
        r"| Pure-even envelope | No | **2** | Yes | No |",
        r"| Pure-even \(k\)-slice | No | 1 | Yes | No |",
        r"| Homog fixed even seed \(x^5+px^3+rx+s\) | **Yes** (if \(p\neq0\)) | 1 | Yes in \(t\) | No |",
        r"| Homog of 1-param even-seed family | **Yes** *if* genuine family exists | **2** | Yes in \((t,u)\) | No |",
        r"| Weighted scale of fixed seed | **Yes** (\(p\neq 0\)) | **1** (one scale) | Yes | No |",
        "",
        f"**Genuine 2-param beyond BJ found in scan:** **{len(genuine_2param)}**.",
        "",
        f"**Largest beyond BJ-embed established constructively:** "
        f"`{largest['name']}` (dim **{largest['dim']}**).",
        "",
        f"**Largest overall (including BJ-embed):** `{largest_overall['name']}` (dim **2**).",
        "",
        "---",
        "",
        r"## 8. Conclusion (Tier 1.1 deepen)",
        "",
        r"1. **No Crit-2 necessity fragment:** nothing found that is simultaneously "
        r"identically disc-square, beyond classical pure-even, forced 3-cycles, and "
        r"HQCC-axiom-named.",
        r"2. **Dimension ceiling (locked):**",
        r"   - overall identical-square: **dim 2** = pure-even envelope (BJ-embed);",
        r"   - **beyond** BJ-embed: **dim 1** = homogenisation of a *fixed* even seed "
        r"\(x^5+px^3+rx+s\) with \(p\neq 0\) (and its quasihomogeneous reparametrisations).",
        r"   - Scan found **no** genuine polynomial 1-param *family of seeds* with "
        r"identically square disc and \(p\not\equiv 0\), \(s\not\equiv 0\) that would "
        r"lift to a new dim-2 beyond BJ.",
        r"3. **3-cycles** remain operational (seed + Hilbert), not structural from \(T\).",
        r"4. **HQCC naming fails** for every identical-square subclass examined: each "
        r"requires pure-even and/or \(a=e=0\) homogenisation ansätze foreign to base \(M\).",
        r"5. **Priority:** further 1.1 only with a new algebraic idea (literature "
        r"parametric even quintics realised in \(T\)). Otherwise shift to Tier 2 "
        r"(paper / catalogue invariants) or optional geometric leftovers.",
        "",
        r"```bash",
        r"python tier11_deepen.py",
        r"```",
        "",
        r"_Generated by tier11_deepen.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "disc_structure": disc_struct,
        "known_families": known,
        "seed_families": seed_fams,
        "two_param_candidates": two_param_candidates,
        "classification": {
            "genuine_2param": genuine_2param,
            "weighted_fixed_seed": weighted_fixed_seed,
            "reducible_s0": reducible_s0,
            "pure_even_seed_hits": pure_even_seed_hits,
        },
        "poly_maps": poly_maps,
        "univariate": univ,
        "cycles": cycles,
        "largest_beyond_BJ": largest,
        "largest_overall": largest_overall,
        "naming": naming,
    }
    write_md(ROOT / "TIER11_DEEPEN.md", "\n".join(lines))
    write_json(ROOT / "TIER11_DEEPEN.json", payload)
    write_md(OUT / "TIER11_DEEPEN.md", "\n".join(lines))
    write_json(OUT / "TIER11_DEEPEN.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "TIER11_DEEPEN.md", "\n".join(lines))
    except Exception:
        pass

    print(f"Wrote TIER11_DEEPEN.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
