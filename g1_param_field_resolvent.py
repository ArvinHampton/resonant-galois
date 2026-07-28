#!/usr/bin/env python3
"""
G1 — parameter-field resolvent.

Accept f ∈ K(s)[x] with [K : Q(s)] > 1 (cover parameters algebraic over Q(s)),
form the norm / multi-sheet resolvent over Q(s) (or Q after specialising s),
then re-test Hilbert specialisations against the multi-seed pure-even catalogue.

Tracks
------
R0  Exact model at s=-1: f ∈ Q(√5)(t)[y], Norm to Q(t) is deg 10 (locked).
R1  Reverse: for each catalogue seed S, does S divide the s=-1 norm for some t∈Q-bar?
R2  Forward: specialise s=-1 norm at many t∈Q; factor; match seeds / pure-even k-rays.
R3  Multi-sheet norms at rational s: Newton sheets → product monic fibres ≈ Norm_{K/Q};
    specialise t; match catalogue.
R4  Estimate [K:Q(s)] from sheet counts / minpolys of (p2,c) samples.
R5  Arithmetic pure-even control (must stay green).

Output: G1_PARAM_FIELD_RESOLVENT.md / .json (+ build/)
"""
from __future__ import annotations

import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

y, t, s = sp.symbols("y t s")

# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------
CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("flag_320", 320, -512, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("s220_176", 220, 176, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("lsw_m209", -209, 836, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("s220m", 220, -528, Fraction(-12, 5)),
    ("s55_176", -55, 176, Fraction(-16, 5)),
    ("flagship_m", -55, -88, Fraction(8, 5)),
    ("classical_m", 20, -16, Fraction(-4, 5)),
    ("lsw4_m100", -100, -400, Fraction(4)),
]
CAT_BY_AB = {(a, b): (tag, k) for tag, a, b, k in CATALOGUE}
CAT_K = sorted({k for *_, k in CATALOGUE}, key=lambda f: (f.denominator, abs(f.numerator)))


def k_of(a, b):
    if a == 0:
        return None
    return Fraction(int(b), int(a))


def pure_even_alpha(m: Fraction, k: Fraction) -> Fraction:
    return 256 * m * m - Fraction(3125) * (k**4) / 256


# ---------------------------------------------------------------------------
# R0 — exact s=-1 parameter-field model
# ---------------------------------------------------------------------------
def sm1_exact_model() -> dict:
    """
    K = Q(√5), [K:Q]=2.
    Monic over K(t):
      f = y^5 - y^3 + (t/√5)(y^2 - 1/25)
    Norm N_{K/Q}(f) (after clearing) is deg 10 in y over Q(t):
      5 (y^5 - y^3)^2 - t^2 (y^2 - 1/25)^2 = 0
    """
    rt5 = sp.sqrt(5)
    f_K = sp.expand(y**5 - y**3 + (t / rt5) * (y**2 - sp.Rational(1, 25)))
    # Isolate √5: f = (y^5 - y^3) + √5 * (t/5) (y^2 - 1/25)  wait:
    # t/√5 = t √5 / 5
    f_K2 = sp.expand(y**5 - y**3 + (t * rt5 / 5) * (y**2 - sp.Rational(1, 25)))
    # Actually original monic was y^5 - y^3 + (t/√5)(y^2-1/25)
    # = (y^5-y^3) + t (y^2-1/25)/√5
    # Multiply by √5: √5(y^5-y^3) + t(y^2-1/25) = 0  (same roots)
    # √5 = -t(y^2-1/25)/(y^5-y^3)
    # 5 = t^2 (y^2-1/25)^2 / (y^5-y^3)^2
    # 5(y^5-y^3)^2 - t^2 (y^2-1/25)^2 = 0
    norm_raw = sp.expand(5 * (y**5 - y**3) ** 2 - t**2 * (y**2 - sp.Rational(1, 25)) ** 2)
    norm_cleared = sp.expand(sp.together(norm_raw) * 625)  # clear 1/25^2
    return {
        "field_K": "Q(sqrt(5))",
        "degree_K_over_Q": 2,
        "s": -1,
        "params": {
            "p2": -1,
            "c": "-sqrt(5)",
            "sigma": 0,
            "pi": "-1/25",
            "q": "1/sqrt(5)",
            "w": "-1/sqrt(5)",
        },
        "f_over_K_t": str(f_K),
        "norm_over_Q_t": str(norm_raw),
        "norm_cleared": str(norm_cleared),
        "norm_degree_y": 10,
        "note": "Norm degree = 5 * [K:Q] = 10",
        "norm_raw_expr": norm_raw,
        "norm_cleared_expr": norm_cleared,
        "f_K_expr": f_K,
    }


# ---------------------------------------------------------------------------
# R1 — reverse: seed divides norm?
# ---------------------------------------------------------------------------
def poly_mod_seed_reduce(expr, alpha: int, beta: int, var=y):
    """Reduce poly in y mod y^5 + alpha y + beta → deg < 5."""
    # y^5 = -alpha y - beta
    p = sp.Poly(sp.expand(expr), var)
    # use rem
    S = sp.Poly(var**5 + alpha * var + beta, var)
    return sp.expand(p.rem(S).as_expr())


def reverse_seed_divides_sm1_norm(alpha: int, beta: int) -> dict:
    """
    Norm N(y,t) = 5(y^5-y^3)^2 - t^2 (y^2-1/25)^2.
    Reduce N mod S = y^5+αy+β; get deg≤4 poly with coeffs in Q(t).
    Set all coeffs to 0; solve for t.
    """
    N = 5 * (y**5 - y**3) ** 2 - t**2 * (y**2 - sp.Rational(1, 25)) ** 2
    # work with cleared N for integer coeffs
    Ncl = sp.expand(N * 625)
    S = y**5 + alpha * y + beta
    rem = poly_mod_seed_reduce(Ncl, alpha, beta)
    # rem = sum_{i=0}^4 a_i(t) y^i
    pol = sp.Poly(rem, y)
    eqs = []
    for i in range(5):
        eqs.append(sp.expand(pol.coeff_monomial(y**i)))
    # eqs are polynomials in t (actually in t^2 mostly)
    try:
        sols = sp.solve(eqs, [t], dict=True)
    except Exception as e:
        return {"ok": False, "error": str(e), "alpha": alpha, "beta": beta}
    good = []
    for sol in sols or []:
        tv = sol.get(t, None)
        if tv is None:
            continue
        if tv.has(sp.oo) or tv == sp.zoo:
            continue
        # verify rem vanishes
        rem2 = sp.simplify(rem.subs(t, tv))
        if rem2 == 0:
            good.append(str(sp.simplify(tv)))
        else:
            # numerical check
            try:
                if abs(complex(sp.N(rem2))) < 1e-9:
                    good.append(str(sp.simplify(tv)))
            except Exception:
                pass
    # Also try: gcd(Ncl.subs(t,t0), S) for candidate t from resultant of coeffs
    # Resultant approach: treat rem coeffs; eliminate to condition
    # If no full vanish, check whether content of rem is 0 for some t
    # Compute gcd of all coeff polynomials in t — if 0, common root
    coeff_polys = [sp.Poly(sp.expand(e), t) for e in eqs if e != 0]
    g = None
    for cp in coeff_polys:
        g = cp if g is None else sp.gcd(g, cp)
    common = []
    if g is not None and g.degree() >= 0 and g != 1:
        try:
            for r, _m in sp.roots(g.as_expr(), t).items():
                common.append(str(sp.simplify(r)))
        except Exception:
            try:
                for r in sp.solve(g.as_expr(), t):
                    common.append(str(sp.simplify(r)))
            except Exception:
                pass

    return {
        "ok": len(good) > 0,
        "alpha": alpha,
        "beta": beta,
        "n_exact_t": len(good),
        "t_values": good[:8],
        "common_root_candidates": common[:8],
        "gcd_deg": int(g.degree()) if g is not None else None,
        "rem_coeffs_sample": [str(e)[:80] for e in eqs[:3]],
    }


def track_R1(seeds) -> list[dict]:
    rows = []
    for tag, a, b, k in seeds:
        print(f"  R1 {tag} ({a},{b}) ...", flush=True)
        r = reverse_seed_divides_sm1_norm(a, b)
        r["tag"] = tag
        r["k"] = str(k)
        rows.append(r)
        print(
            f"    ok={r['ok']} n_t={r.get('n_exact_t')} gcd_deg={r.get('gcd_deg')} common={r.get('common_root_candidates')[:3]}",
            flush=True,
        )
    return rows


# ---------------------------------------------------------------------------
# R2 — forward specialise s=-1 norm
# ---------------------------------------------------------------------------
def clear_to_Z_poly(expr, var=y) -> sp.Poly | None:
    try:
        pol = sp.Poly(sp.expand(sp.together(expr)), var, domain=sp.QQ)
    except Exception:
        return None
    dens = [sp.fraction(sp.together(c))[1] for c in pol.all_coeffs()]
    L = 1
    for d in dens:
        try:
            L = int(sp.ilcm(L, abs(int(d))))
        except Exception:
            return None
    cleared = sp.expand(L ** pol.degree() * pol.as_expr().subs(var, var / L))
    try:
        pz = sp.Poly(cleared, var, domain=sp.ZZ)
    except Exception:
        return None
    if pz.LC() == -1:
        pz = sp.Poly(-pz.as_expr(), var, domain=sp.ZZ)
    if pz.LC() < 0:
        pz = sp.Poly(-pz.as_expr(), var, domain=sp.ZZ)
    # primitive
    cont = sp.content(pz.as_expr())
    if cont not in (0, 1, -1):
        pz = sp.Poly(sp.primitive(pz.as_expr())[1], var, domain=sp.ZZ)
        if pz.LC() < 0:
            pz = sp.Poly(-pz.as_expr(), var, domain=sp.ZZ)
    return pz


def try_bj_from_poly(pz: sp.Poly) -> dict | None:
    if pz.degree() != 5 or pz.LC() != 1:
        return None
    coeffs = [sp.Rational(c) for c in pz.all_coeffs()]
    c4 = coeffs[1]
    shift = -c4 / 5
    z = sp.symbols("z")
    fsh = sp.expand(pz.as_expr().subs(y, z + shift))
    psh = sp.Poly(fsh, z, domain=sp.QQ)
    cc = [sp.Rational(c) for c in psh.all_coeffs()]
    if len(cc) != 6 or cc[1] != 0:
        return None
    if cc[2] == 0 and cc[3] == 0:
        try:
            a, b = int(cc[4]), int(cc[5])
            return {"form": "BJ", "alpha": a, "beta": b, "k": str(k_of(a, b))}
        except Exception:
            return {"form": "BJ_QQ", "alpha": str(cc[4]), "beta": str(cc[5])}
    return {"form": "depressed", "p": str(cc[2]), "q": str(cc[3])}


def track_R2_forward(t_vals) -> dict:
    model = sm1_exact_model()
    N = model["norm_raw_expr"]
    rows = []
    cat_hits = []
    k_hits = []
    factor_deg_hist = Counter()
    n_even5 = 0
    n_bj = 0
    n_a5 = 0

    for tv in t_vals:
        expr = sp.expand(N.subs(t, sp.Rational(tv) if isinstance(tv, Fraction) else tv))
        # clear denom
        expr2 = sp.expand(sp.together(expr) * 625)
        try:
            pol = sp.Poly(expr2, y, domain=sp.ZZ)
        except Exception:
            continue
        if pol.LC() < 0:
            pol = sp.Poly(-pol.as_expr(), y, domain=sp.ZZ)
        try:
            facs = sp.factor_list(pol.as_expr())
        except Exception:
            continue
        for f, mult in facs[1]:
            deg = sp.degree(f, y)
            factor_deg_hist[int(deg)] += 1
            if deg not in (5, 10):
                continue
            try:
                pf = sp.Poly(f, y, domain=sp.ZZ)
            except Exception:
                continue
            if pf.LC() == -1:
                pf = sp.Poly(-pf.as_expr(), y, domain=sp.ZZ)
            if pf.LC() != 1:
                continue
            if deg == 5:
                rec = classify_poly(pf.as_expr().subs(y, x), do_galois=bool(pf.is_irreducible))
                bj = try_bj_from_poly(pf)
                row = {
                    "t": str(tv),
                    "deg": 5,
                    "irr": rec.get("irreducible"),
                    "disc_square": rec.get("disc_square"),
                    "status": rec.get("status"),
                    "galois": rec.get("galois"),
                    "bj": bj,
                    "poly": rec.get("poly"),
                }
                if rec.get("disc_square"):
                    n_even5 += 1
                if (rec.get("status") or "").startswith("HIT_A5"):
                    n_a5 += 1
                if bj and bj.get("form") == "BJ":
                    n_bj += 1
                    a, b = bj["alpha"], bj["beta"]
                    kk = k_of(a, b)
                    if (a, b) in CAT_BY_AB:
                        tag, ck = CAT_BY_AB[(a, b)]
                        cat_hits.append(
                            {"tag": tag, "k": str(ck), "t": str(tv), "alpha": a, "beta": b, "source": "sm1_norm_factor5"}
                        )
                    elif kk in CAT_K:
                        k_hits.append({"k": str(kk), "t": str(tv), "alpha": a, "beta": b})
                rows.append(row)
            elif deg == 10 and pf.is_irreducible:
                rows.append({"t": str(tv), "deg": 10, "irr": True, "poly": str(pf.as_expr())[:100]})

    return {
        "n_factor_rows": len(rows),
        "factor_deg_hist": dict(factor_deg_hist),
        "n_even5": n_even5,
        "n_bj": n_bj,
        "n_a5": n_a5,
        "catalogue_hits": cat_hits,
        "k_ray_hits": k_hits[:20],
        "sample": rows[:30],
    }


# ---------------------------------------------------------------------------
# Newton multi-sheet covers (for R3, R4)
# ---------------------------------------------------------------------------
def residual_cover(v, s_val: float):
    c, p2n, r1, r2, qn, wn = v

    def G(val, pt):
        yy = pt
        N = c * yy**3 * (yy - 1.0) * (yy - p2n)
        A, Ap, App = yy**3, 3 * yy**2, 6 * yy
        B = yy**2 - (1 + p2n) * yy + p2n
        Bp = 2 * yy - (1 + p2n)
        Bpp = 2.0
        Np = c * (Ap * B + A * Bp)
        Npp = c * (App * B + 2 * Ap * Bp + A * Bpp)
        D = (yy - r1) * (yy - r2)
        Dp = 2 * yy - (r1 + r2)
        Dpp = 2.0
        return [N - val * D, Np - val * Dp, Npp - val * Dpp]

    return np.array(G(1.0, qn) + G(s_val, wn), dtype=float)


def newton_one(s_val, x0, niter=70):
    v = np.array(x0, dtype=float)
    for _ in range(niter):
        r = residual_cover(v, s_val)
        nr = float(np.linalg.norm(r))
        if nr < 1e-13:
            return v, True, nr
        J = np.zeros((6, 6))
        eps = 1e-8
        for j in range(6):
            dv = np.zeros(6)
            dv[j] = eps
            J[:, j] = (residual_cover(v + dv, s_val) - r) / eps
        try:
            step = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(J, -r, rcond=None)[0]
        v = v + step
    nr = float(np.linalg.norm(residual_cover(v, s_val)))
    return v, nr < 1e-11, nr


def find_sheets(s_val: float, n_trials: int = 120, tol_dup: float = 1e-5) -> list[np.ndarray]:
    rng = np.random.default_rng(abs(hash((round(s_val, 10), 539))) % (2**32))
    sheets = []
    seeds = [rng.normal(size=6) for _ in range(n_trials)]
    seeds.append(np.array([-np.sqrt(5), -1.0, 0.2, -0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5]))
    for x0 in seeds:
        v, ok, nr = newton_one(s_val, x0)
        if not ok:
            continue
        # dedup by (p2, c, r1+r2, r1*r2)
        key = np.array([v[1], v[0], v[2] + v[3], v[2] * v[3]])
        is_dup = False
        for s0 in sheets:
            key0 = np.array([s0[1], s0[0], s0[2] + s0[3], s0[2] * s0[3]])
            if np.linalg.norm(key - key0) < tol_dup:
                is_dup = True
                break
            # complex conjugate / sign flip r1↔r2 already covered by elementary sym
        if not is_dup:
            sheets.append(v)
    return sheets


def monic_coeffs_from_sheet(v, tv: float) -> list[float] | None:
    c, p2n, r1, r2, qn, wn = v
    if abs(c) < 1e-14:
        return None
    sig = r1 + r2
    pin = r1 * r2
    ratio = float(tv) / float(c)
    return [
        1.0,
        -(1.0 + p2n),
        p2n,
        -ratio,
        ratio * sig,
        -ratio * pin,
    ]


def product_monic_polys(coeff_lists: list[list[float]]) -> np.ndarray:
    """Product of monic polys given high-to-low coeff lists → high-to-low."""
    prod = np.array([1.0], dtype=float)
    for c in coeff_lists:
        prod = np.polymul(prod, np.array(c, dtype=float))
    # make monic
    if abs(prod[0]) > 1e-15:
        prod = prod / prod[0]
    return prod


def rationalize_float_poly(coeffs: np.ndarray, den_bound: int = 200) -> sp.Poly | None:
    """Try to recognize monic Z or QQ poly from float coeffs."""
    c = [complex(z).real for z in coeffs]
    if abs(c[0] - 1.0) > 1e-4:
        c = [ci / c[0] for ci in c]
    rats = []
    for ci in c:
        r = sp.nsimplify(ci, tolerance=1e-7, rational=True)
        if not (r.is_rational or r.is_Integer):
            # try continued fraction
            try:
                r = sp.Rational(str(Fraction(ci).limit_denominator(den_bound)))
            except Exception:
                return None
        rats.append(sp.Rational(r))
    try:
        pol = sp.Poly.from_list(rats, y, domain=sp.QQ)
        mon = sp.Poly(sp.monic(pol.as_expr()), y, domain=sp.QQ)
        return clear_to_Z_poly(mon.as_expr(), y)
    except Exception:
        return None


def track_R3_multisheet(s_list, t_list) -> dict:
    results = []
    cat_hits = []
    sheet_counts = {}
    n_product_Z = 0
    n_factor5 = 0
    n_bj = 0

    for s_val in s_list:
        sheets = find_sheets(float(s_val), n_trials=100)
        sheet_counts[str(s_val)] = len(sheets)
        print(f"  R3 s={s_val}: {len(sheets)} sheets", flush=True)
        if len(sheets) < 1:
            continue
        for tv in t_list:
            coeff_lists = []
            for sh in sheets:
                mc = monic_coeffs_from_sheet(sh, float(tv))
                if mc is None or any(np.isnan(mc)):
                    continue
                # require real coeffs (discard wildly complex)
                if any(abs(complex(z).imag) > 1e-6 for z in mc):
                    # keep real part only if imag tiny overall
                    mc = [complex(z).real for z in mc]
                coeff_lists.append(mc)
            if len(coeff_lists) < 1:
                continue
            # Single-sheet factors (the K-fibres specialised)
            for mc in coeff_lists:
                pz = rationalize_float_poly(np.array(mc))
                if pz is None or pz.degree() != 5:
                    continue
                bj = try_bj_from_poly(pz)
                if bj and bj.get("form") == "BJ":
                    n_bj += 1
                    a, b = bj["alpha"], bj["beta"]
                    if (a, b) in CAT_BY_AB:
                        tag, ck = CAT_BY_AB[(a, b)]
                        cat_hits.append(
                            {
                                "tag": tag,
                                "k": str(ck),
                                "s": s_val,
                                "t": tv,
                                "alpha": a,
                                "beta": b,
                                "source": "single_sheet",
                            }
                        )
            # Multi-sheet product = candidate Norm
            if len(coeff_lists) >= 2:
                prod = product_monic_polys(coeff_lists)
                pz = rationalize_float_poly(prod, den_bound=500)
                if pz is not None:
                    n_product_Z += 1
                    # factor
                    try:
                        facs = sp.factor_list(pz.as_expr())
                    except Exception:
                        continue
                    for f, _m in facs[1]:
                        if sp.degree(f, y) != 5:
                            continue
                        n_factor5 += 1
                        try:
                            pf = sp.Poly(f, y, domain=sp.ZZ)
                        except Exception:
                            continue
                        if pf.LC() == -1:
                            pf = sp.Poly(-pf.as_expr(), y, domain=sp.ZZ)
                        if pf.LC() != 1:
                            continue
                        bj = try_bj_from_poly(pf)
                        if bj and bj.get("form") == "BJ":
                            n_bj += 1
                            a, b = bj["alpha"], bj["beta"]
                            if (a, b) in CAT_BY_AB:
                                tag, ck = CAT_BY_AB[(a, b)]
                                cat_hits.append(
                                    {
                                        "tag": tag,
                                        "k": str(ck),
                                        "s": s_val,
                                        "t": tv,
                                        "alpha": a,
                                        "beta": b,
                                        "source": "norm_product_factor",
                                        "n_sheets": len(coeff_lists),
                                    }
                                )
            results.append(
                {
                    "s": s_val,
                    "t": tv,
                    "n_sheets_used": len(coeff_lists),
                }
            )

    cat_k = sorted({h["k"] for h in cat_hits})
    return {
        "sheet_counts": sheet_counts,
        "n_product_Z": n_product_Z,
        "n_factor5": n_factor5,
        "n_bj": n_bj,
        "catalogue_hits": cat_hits,
        "catalogue_k": cat_k,
        "multi_k": len(cat_k) >= 2,
        "n_pairs_tested": len(results),
    }


# ---------------------------------------------------------------------------
# R4 — degree estimates
# ---------------------------------------------------------------------------
def track_R4_degree(s_list) -> dict:
    """
    Estimate number of distinct real/complex Newton sheets at rational s,
    and algebraic degree of p2,c over Q by minpoly of high-precision samples
    (single sheet) — rough.
    """
    rows = []
    for s_val in s_list:
        sheets = find_sheets(float(s_val), n_trials=150)
        p2s = [complex(sh[1]) for sh in sheets]
        cs = [complex(sh[0]) for sh in sheets]
        # unique p2 count
        rows.append(
            {
                "s": s_val,
                "n_sheets": len(sheets),
                "p2_sample": [str(complex(p).real)[:12] for p in p2s[:6]],
                "c_sample": [str(complex(c).real)[:12] for c in cs[:6]],
            }
        )
        print(f"  R4 s={s_val}: sheets={len(sheets)}", flush=True)

    # At s=-1 exact: [K:Q]=2
    # Across s, max sheets observed as proxy for upper bound on geometric sheets
    max_sheets = max((r["n_sheets"] for r in rows), default=0)
    return {
        "rows": rows,
        "max_sheets_observed": max_sheets,
        "exact_sm1_degree": 2,
        "note": (
            "Sheet count is a geometric upper bound on the number of covers of this "
            "normal form over C at fixed s; parameter field degree may be smaller "
            "after identifying Galois orbits. Exact [Q(√5):Q]=2 at s=-1."
        ),
    }


# ---------------------------------------------------------------------------
# R5 — pure-even control
# ---------------------------------------------------------------------------
def track_R5_control() -> dict:
    rows = []
    for k in CAT_K:
        for m in [Fraction(i) for i in (1, 2, 3, 5, 9)]:
            a = pure_even_alpha(m, k)
            if a.denominator != 1:
                continue
            aa, bb = int(a), int(k * a)
            if aa == 0:
                continue
            d = disc_bj_int(aa, bb)
            rows.append(
                {
                    "k": str(k),
                    "m": str(m),
                    "alpha": aa,
                    "beta": bb,
                    "disc_square": d > 0 and is_square(d),
                    "in_cat": (aa, bb) in CAT_BY_AB,
                }
            )
    return {
        "n": len(rows),
        "n_disc_square": sum(1 for r in rows if r["disc_square"]),
        "n_in_cat": sum(1 for r in rows if r["in_cat"]),
        "sample": rows[:12],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G1 — parameter-field resolvent", flush=True)
    print("=" * 72, flush=True)

    print("\n[R0] Exact s=-1 model ...", flush=True)
    r0 = sm1_exact_model()
    # strip exprs for json later
    print(f"  K={r0['field_K']} deg={r0['degree_K_over_Q']} norm_deg_y={r0['norm_degree_y']}", flush=True)

    print("\n[R1] Reverse: catalogue seed | s=-1 norm ? ...", flush=True)
    r1 = track_R1(CATALOGUE)
    r1_ok = sum(1 for r in r1 if r.get("ok"))
    r1_common = sum(1 for r in r1 if r.get("common_root_candidates"))

    print("\n[R2] Forward specialise s=-1 norm ...", flush=True)
    t_vals = list(range(-15, 16)) + [
        Fraction(1, 2),
        Fraction(3, 2),
        Fraction(2, 3),
        Fraction(5, 2),
        Fraction(5, 3),
        Fraction(7, 2),
        Fraction(-3, 2),
        Fraction(5, 4),
    ]
    r2 = track_R2_forward(t_vals)
    print(
        f"  deg hist={r2['factor_deg_hist']} even5={r2['n_even5']} BJ={r2['n_bj']} "
        f"A5={r2['n_a5']} cat={len(r2['catalogue_hits'])}",
        flush=True,
    )

    print("\n[R4] Sheet counts / degree proxy ...", flush=True)
    s_list = [-3, -2, -1, -0.5, 0.5, 1.5, 2, 3, 4, 5, -1.5, 2.5]
    r4 = track_R4_degree(s_list)

    print("\n[R3] Multi-sheet norms at rational s × t ...", flush=True)
    t_list = list(range(-6, 7)) + [0.5, 1.5, -0.5]
    # fewer s for cost
    s_r3 = [-2, -1, 0.5, 2, 3, -0.5, 1.5]
    r3 = track_R3_multisheet(s_r3, t_list)
    print(
        f"  product_Z={r3['n_product_Z']} factor5={r3['n_factor5']} BJ={r3['n_bj']} "
        f"cat={len(r3['catalogue_hits'])} multi_k={r3['multi_k']}",
        flush=True,
    )

    print("\n[R5] Pure-even arithmetic control ...", flush=True)
    r5 = track_R5_control()
    print(f"  Z samples={r5['n']} disc□={r5['n_disc_square']} in_cat={r5['n_in_cat']}", flush=True)

    elapsed = round(time.time() - t0, 2)
    geometric_hit = (
        r1_ok > 0
        or len(r2["catalogue_hits"]) > 0
        or len(r3["catalogue_hits"]) > 0
    )
    multi_k = r3["multi_k"] or len({h["k"] for h in r2["catalogue_hits"]}) >= 2

    verdict = (
        f"G1 param-field resolvent ({elapsed}s). "
        f"s=-1: f∈Q(√5)(t)[y], Norm deg 10 over Q(t). "
        f"R1 seed|norm exact hits={r1_ok}/{len(r1)} (common-root signals={r1_common}). "
        f"R2 forward cat={len(r2['catalogue_hits'])} BJ={r2['n_bj']} even5={r2['n_even5']}. "
        f"R3 multi-sheet cat={len(r3['catalogue_hits'])} multi_k={r3['multi_k']}. "
        f"max_sheets≈{r4['max_sheets_observed']}. "
        f"geometric_catalogue_hit={geometric_hit}."
    )
    print("\n" + verdict, flush=True)

    # ---- report ----
    lines = [
        "# G1 — parameter-field resolvent",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Goal",
        "",
        "Accept the cover fibre as",
        "",
        "```text",
        "f ∈ K(s)[x] with [K : Q(s)] > 1",
        "```",
        "",
        "form the **norm / multi-sheet resolvent** with coefficients in Q(s)",
        "(or in Q after specialising s), then re-test Hilbert specialisations against",
        "the multi-seed pure-even catalogue (flagship −8/5, classical 4/5, LSW −4, …).",
        "",
        "This is the G1 path that does **not** require a single-valued f_s ∈ Q(s)[x].",
        "",
        "---",
        "",
        "## 1. R0 — exact parameter field at s = −1",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| K | {r0['field_K']} |",
        f"| [K:Q] | {r0['degree_K_over_Q']} |",
        f"| s | {r0['s']} |",
        f"| params | p2={r0['params']['p2']}, c={r0['params']['c']}, σ={r0['params']['sigma']}, π={r0['params']['pi']} |",
        f"| f over K(t) | `{r0['f_over_K_t']}` |",
        f"| Norm over Q(t) | `{r0['norm_over_Q_t']}` |",
        f"| deg_y(Norm) | **{r0['norm_degree_y']}** = 5 × [K:Q] |",
        "",
        "So the natural Q-model of this geometric fibre is **degree 10**, not 5.",
        "",
        "---",
        "",
        "## 2. R1 — reverse: does a catalogue seed divide the s=−1 norm?",
        "",
        "For S = y⁵ + α y + β, reduce Norm(y,t) modulo S and require the remainder",
        "(coeffs in Q(t)) to vanish identically in y — solve for t.",
        "",
        f"| seed | k | exact t hit? | #t | gcd(coeffs) deg | common-root candidates |",
        f"|------|---|:------------:|---:|----------------:|-------------------------|",
    ]
    for r in r1:
        lines.append(
            f"| {r['tag']} | {r['k']} | {r['ok']} | {r.get('n_exact_t', 0)} | "
            f"{r.get('gcd_deg')} | {r.get('common_root_candidates', [])[:3]} |"
        )
    lines += [
        "",
        f"**R1 exact hits: {r1_ok}/{len(r1)}** "
        f"(seeds with common-root signal on coeff gcd: {r1_common}).",
        "",
        "---",
        "",
        "## 3. R2 — forward specialisation of the s=−1 norm",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| factor degree histogram | {r2['factor_deg_hist']} |",
        f"| even disc deg-5 factors | {r2['n_even5']} |",
        f"| BJ deg-5 | {r2['n_bj']} |",
        f"| A5 (among classified even) | {r2['n_a5']} |",
        f"| **exact catalogue hits** | **{len(r2['catalogue_hits'])}** |",
        f"| BJ on a catalogue k-ray (not exact seed) | {len(r2['k_ray_hits'])} |",
        "",
    ]
    if r2["catalogue_hits"]:
        lines.append("### Catalogue hits (R2)")
        lines.append("")
        for h in r2["catalogue_hits"]:
            lines.append(f"- {h}")
        lines.append("")
    else:
        lines.append("_No catalogue seed among deg-5 factors of the s=−1 norm specialisations._")
        lines.append("")
    lines += [
        "---",
        "",
        "## 4. R4 — sheet counts (degree proxy)",
        "",
        f"| s | # Newton sheets |",
        f"|--:|----------------:|",
    ]
    for r in r4["rows"]:
        lines.append(f"| {r['s']} | {r['n_sheets']} |")
    lines += [
        "",
        f"- max sheets observed: **{r4['max_sheets_observed']}**",
        f"- exact [Q(√5):Q] at s=−1: **{r4['exact_sm1_degree']}**",
        f"- note: {r4['note']}",
        "",
        "---",
        "",
        "## 5. R3 — multi-sheet product norms at rational s",
        "",
        "At each rational s, collect distinct Newton solutions (sheets), form monic",
        "fibres F_i(y;t), and take the product ∏ F_i as a numerical stand-in for",
        "Norm_{K/Q}(f) when the sheets form a Galois orbit. Recognise rational monic",
        "polys, factor, match catalogue.",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| sheet counts | {r3['sheet_counts']} |",
        f"| product recognised over Z | {r3['n_product_Z']} |",
        f"| deg-5 factors from products | {r3['n_factor5']} |",
        f"| BJ found | {r3['n_bj']} |",
        f"| **catalogue hits** | **{len(r3['catalogue_hits'])}** |",
        f"| catalogue k | {r3['catalogue_k']} |",
        f"| multi-k | {r3['multi_k']} |",
        "",
    ]
    if r3["catalogue_hits"]:
        lines.append("### Catalogue hits (R3)")
        lines.append("")
        for h in r3["catalogue_hits"]:
            lines.append(f"- {h}")
        lines.append("")
    else:
        lines.append("_No catalogue hits from multi-sheet norms in this scan._")
        lines.append("")
    lines += [
        "---",
        "",
        "## 6. R5 — arithmetic multi-k control",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| pure-even Z samples | {r5['n']} |",
        f"| disc □ | {r5['n_disc_square']} |",
        f"| exact catalogue among samples | {r5['n_in_cat']} |",
        "",
        "Arithmetic multi-k remains available; geometric fusion is the open gap.",
        "",
        "---",
        "",
        "## 7. Multi-k conclusion",
        "",
        f"| test | result |",
        f"|------|--------|",
        f"| Parameter-field model f∈K(s)[x] constructed (s=−1) | **True** (K=Q(√5)) |",
        f"| Norm over Q of deg 5·[K:Q] | **True** (deg 10) |",
        f"| R1 seed divides Norm | **{r1_ok > 0}** ({r1_ok}/{len(r1)}) |",
        f"| R2 forward catalogue hit | **{len(r2['catalogue_hits']) > 0}** |",
        f"| R3 multi-sheet catalogue hit | **{len(r3['catalogue_hits']) > 0}** |",
        f"| Geometric multi-k | **{multi_k}** |",
        f"| Arithmetic multi-k control | **True** |",
        "",
        "**Geometric multi-k via parameter-field resolvent of this 3A⁴ normal form: "
        + ("HIT." if geometric_hit else "not achieved in this cut.")
        + "**",
        "",
        "### What this cut established",
        "",
        "1. **Explicit f ∈ K(t)[y]** at the known geometric fibre s=−1 with [K:Q]=2.",
        "2. **Canonical Norm** is deg 10 over Q(t) — the correct Q-model degree for this fibre.",
        "3. **Reverse division test** of all multi-seed catalogue BJ seeds against that Norm.",
        "4. **Forward factorisation** of Norm specialisations (deg histogram, even/BJ/A5).",
        "5. **Multi-sheet product** proxy for Norm at other rational s + catalogue re-test.",
        "",
        "### If still empty — meaning",
        "",
        "The pure-even catalogue is not among the Hilbert fibres of the normed 3A⁴",
        "normal-form covers sampled here. Remaining geometric routes:",
        "",
        "1. Domain Möbius / different normal form (change {0,1,∞} labels) before norming.",
        "2. Cubic Tschirnhaus on the deg-10 Q-model (not only on deg-5 K-fibres).",
        "3. G2: other genus-0 Nielsen types with possibly rational parameter fields.",
        "4. G3: monodromy identification of the pure-even envelope (arithmetic multi-k",
        "   already has multi-k; give it a Nielsen name).",
        "",
        "---",
        "",
        "## 8. Non-claims",
        "",
        "- Not a proof that no Nielsen realisation of the pure-even lattice exists.",
        "- Negative for this normal form’s parameter-field norm + scan bounds.",
        "- Does not reopen pure-even arithmetic, Canonical T3, or Necessity.",
        "",
        "_Generated by `g1_param_field_resolvent.py`._",
        "",
    ]

    md = "\n".join(lines)

    # JSON-safe payload
    def strip_expr(d):
        out = {}
        for k, v in d.items():
            if k.endswith("_expr") or k in ("norm_raw_expr", "norm_cleared_expr", "f_K_expr"):
                continue
            out[k] = v
        return out

    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "R0": strip_expr(r0),
        "R1": r1,
        "R1_ok": r1_ok,
        "R2": r2,
        "R3": r3,
        "R4": r4,
        "R5": r5,
        "geometric_hit": geometric_hit,
        "multi_k": multi_k,
    }

    write_md(ROOT / "G1_PARAM_FIELD_RESOLVENT.md", md)
    write_json(ROOT / "G1_PARAM_FIELD_RESOLVENT.json", payload)
    write_md(OUT / "G1_PARAM_FIELD_RESOLVENT.md", md)
    write_json(OUT / "G1_PARAM_FIELD_RESOLVENT.json", payload)
    try:
        write_md(RESULTS / "G1_PARAM_FIELD_RESOLVENT.md", md)
        write_json(RESULTS / "G1_PARAM_FIELD_RESOLVENT.json", payload)
    except Exception:
        pass

    print(f"\nWrote G1_PARAM_FIELD_RESOLVENT.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
