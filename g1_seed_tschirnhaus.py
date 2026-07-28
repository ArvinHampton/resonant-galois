#!/usr/bin/env python3
"""
G1 next: seed-first + domain Tschirnhaus freedom.

G1 baseline (g1_3a4_triple_root.py):
  - Catalogue seeds are reverse-compatible with N−tD after *translation only*
  - Those reverse (p2,σ,π) never sit on the triple-root locus (0/16)

This module enlarges the reverse by a polynomial Tschirnhaus on the fibre
coordinate:
  τ(y) = c0 + c1 y + c2 y²
  Res_y( F(y), X − τ(y) )  =?=  X⁵ + α X + β   (catalogue seed)

where F is the monic N−tD fibre
  F = y⁵ − (1+p2)y⁴ + p2 y³ − λ(y² − σ y + π)

and (p2,σ,π) are required to lie on the triple-root locus
  (chart from (q,w) with P(q,w)=0).

Tracks
------
T1  Fixed geometric fibre s=−1 (p2,σ,π locked): free (λ,c0,c1,c2); exact solve per seed.
T2  Locus-constrained reverse: free (q,w,λ,c0,c2) with P=0, c1=1; numeric nsolve.
T3  Forward: numeric (3,1,1)^4 covers at rational s × t-grid → Tschirnhaus→BJ → catalogue.

Output: G1_SEED_TSCHIRNHAUS.md / .json (+ build/)
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, is_square, write_json, write_md  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

y, X = sp.symbols("y X")
p2, sig, pi, lam = sp.symbols("p2 sigma pi lam")
c0, c1, c2 = sp.symbols("c0 c1 c2")
q, w = sp.symbols("q w")

# ---------------------------------------------------------------------------
# Catalogue (multi-seed pure-even ratios)
# ---------------------------------------------------------------------------
CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("s55_176", -55, 176, Fraction(-16, 5)),
    ("flagship_m", -55, -88, Fraction(8, 5)),
    ("classical_m", 20, -16, Fraction(-4, 5)),
]
CAT_BY_AB = {(a, b): (tag, k) for tag, a, b, k in CATALOGUE}


def k_of(a, b):
    if a == 0:
        return None
    return Fraction(int(b), int(a))


def least_squares_numpy(fun, x0, max_nfev=400, step=1e-6, xtol=1e-12, ftol=1e-12, damp0=1e-3):
    """
    Pure-NumPy Levenberg–Marquardt-style least squares (no scipy).
    fun: R^n -> R^m residual vector.
    """
    x = np.array(x0, dtype=float).copy()
    n = x.size
    fx = np.asarray(fun(x), dtype=float).ravel()
    best_x, best_f = x.copy(), fx.copy()
    best_n = float(np.linalg.norm(best_f))
    damp = damp0
    nfev = 1
    for _ in range(max_nfev):
        if best_n < ftol:
            break
        # Jacobian by forward differences
        m = fx.size
        J = np.zeros((m, n), dtype=float)
        for j in range(n):
            xj = x.copy()
            h = step * max(1.0, abs(x[j]))
            xj[j] += h
            fj = np.asarray(fun(xj), dtype=float).ravel()
            nfev += 1
            J[:, j] = (fj - fx) / h
        # LM step: (J^T J + λ I) δ = -J^T f
        JTJ = J.T @ J
        g = J.T @ fx
        improved = False
        for _try in range(8):
            try:
                delta = np.linalg.solve(JTJ + damp * np.eye(n), -g)
            except np.linalg.LinAlgError:
                delta = -g / (np.linalg.norm(g) + 1e-15)
            x_new = x + delta
            f_new = np.asarray(fun(x_new), dtype=float).ravel()
            nfev += 1
            n_new = float(np.linalg.norm(f_new))
            if n_new < best_n * (1.0 - 1e-8) or n_new < best_n - ftol:
                x, fx = x_new, f_new
                best_x, best_f, best_n = x.copy(), f_new.copy(), n_new
                damp = max(damp / 3.0, 1e-12)
                improved = True
                break
            damp = min(damp * 3.0, 1e8)
        if not improved:
            # random restart kick
            if nfev > max_nfev:
                break
            x = best_x + np.random.default_rng(nfev).normal(scale=0.05, size=n)
            fx = np.asarray(fun(x), dtype=float).ravel()
            nfev += 1
            if float(np.linalg.norm(fx)) < best_n:
                best_x, best_f, best_n = x.copy(), fx.copy(), float(np.linalg.norm(fx))
        if float(np.linalg.norm(delta if improved else 0.0)) < xtol and improved:
            break
    class _Sol:
        pass
    sol = _Sol()
    sol.x = best_x
    sol.fun = best_f
    sol.cost = 0.5 * best_n**2
    sol.nfev = nfev
    sol.success = best_n < 1e-6
    return sol


# ---------------------------------------------------------------------------
# Eliminant / chart (locked)
# ---------------------------------------------------------------------------
P = (
    20 * q**3 * w**3
    - 40 * q**3 * w**2
    + 27 * q**3 * w
    - 6 * q**3
    - 40 * q**2 * w**3
    + 73 * q**2 * w**2
    - 45 * q**2 * w
    + 9 * q**2
    + 27 * q * w**3
    - 45 * q * w**2
    + 26 * q * w
    - 5 * q
    - 6 * w**3
    + 9 * w**2
    - 5 * w
    + 1
)

F1 = (
    16 * p2**2 * q * w
    - 8 * p2**2 * q
    - 8 * p2**2 * w
    + 3 * p2**2
    - 30 * p2 * q**2 * w
    + 15 * p2 * q**2
    - 30 * p2 * q * w**2
    + 37 * p2 * q * w
    - 8 * p2 * q
    + 15 * p2 * w**2
    - 8 * p2 * w
    + 50 * q**2 * w**2
    - 30 * q**2 * w
    - 30 * q * w**2
    + 16 * q * w
)


def DEN(p2v, qv):
    return 6 * p2v * qv - 3 * p2v - 10 * qv**2 + 6 * qv


def chart(p2v, qv):
    den = DEN(p2v, qv)
    if den == 0:
        return None
    return {
        "sigma": sp.together(qv * (8 * p2v * qv - 3 * p2v - 15 * qv**2 + 8 * qv) / den),
        "pi": sp.together(qv**2 * (3 * p2v * qv - p2v - 6 * qv**2 + 3 * qv) / den),
        "c": sp.together(-1 / (qv * den)),
    }


def s_form(p2v, qv, wv):
    num = 6 * p2v * wv**2 - 3 * p2v * wv - 10 * wv**3 + 6 * wv**2
    den = 6 * p2v * qv**2 - 3 * p2v * qv - 10 * qv**3 + 6 * qv**2
    if den == 0:
        return None
    return sp.together(num / den)


def fibre_expr(p2v, sigv, piv, lamv):
    return (
        y**5
        - (1 + p2v) * y**4
        + p2v * y**3
        - lamv * (y**2 - sigv * y + piv)
    )


def physical_p2_at(qv, wv):
    f1 = sp.simplify(F1.subs({q: qv, w: wv}))
    sols = sp.solve(sp.Eq(f1, 0), p2)
    if not sols:
        return None
    best, best_key = None, None
    for sol in sols:
        ss = sp.simplify(sol)
        try:
            n = complex(sp.N(ss))
            key = (abs(n.imag), abs(n.real + 1.0), sp.count_ops(ss))
        except Exception:
            key = (99.0, 99.0, sp.count_ops(ss))
        if best is None or key < best_key:
            best, best_key = ss, key
    return best


# ---------------------------------------------------------------------------
# Tschirnhaus resultant machinery
# ---------------------------------------------------------------------------
def tschirnhaus_resultant(F, tau, x_sym=X, y_sym=y):
    """Res_y(F(y), x - τ(y)), monic in x if possible."""
    res = sp.resultant(F, x_sym - tau, y_sym)
    res = sp.expand(res)
    pol = sp.Poly(res, x_sym)
    if pol.degree() < 1:
        return None
    lc = pol.LC()
    if lc == 0:
        return None
    mon = sp.expand(sp.together(res / lc))
    return sp.Poly(mon, x_sym)


def match_seed_coeffs(monic_poly: sp.Poly, alpha: int, beta: int) -> list:
    """
    monic poly = X^5 + a4 X^4 + a3 X^3 + a2 X^2 + a1 X + a0
    want a4=a3=a2=0, a1=alpha, a0=beta.
    """
    if monic_poly is None or monic_poly.degree() != 5:
        return [sp.Integer(1)]  # impossible sentinel
    coeffs = monic_poly.all_coeffs()  # [1, a4, a3, a2, a1, a0]
    a4, a3, a2, a1, a0 = coeffs[1:]
    return [
        sp.expand(a4),
        sp.expand(a3),
        sp.expand(a2),
        sp.expand(a1 - alpha),
        sp.expand(a0 - beta),
    ]


def translation_only_reverse(alpha: int, beta: int) -> list[dict]:
    """G1 baseline: τ = y + δ only ⇔ monic N-tD after y^4-kill = seed."""
    p2v, sigv, piv, lamv = sp.symbols("p2v sigv piv lamv")
    delta = (1 + p2v) / 5
    z = sp.symbols("z")
    f = sp.expand(
        (z + delta) ** 5
        - (1 + p2v) * (z + delta) ** 4
        + p2v * (z + delta) ** 3
        - lamv * ((z + delta) ** 2 - sigv * (z + delta) + piv)
    )
    pol = sp.Poly(f, z)
    eqs = [
        pol.coeff_monomial(z**3),
        pol.coeff_monomial(z**2),
        pol.coeff_monomial(z) - alpha,
        pol.coeff_monomial(1) - beta,
    ]
    try:
        sols = sp.solve(eqs, [p2v, sigv, piv, lamv], dict=True)
    except Exception:
        return []
    out = []
    for sol in sols:
        if any(v.has(sp.oo) or v == sp.zoo for v in sol.values()):
            continue
        out.append(
            {
                "p2": sp.simplify(sol[p2v]),
                "sigma": sp.simplify(sol[sigv]),
                "pi": sp.simplify(sol[piv]),
                "lam": sp.simplify(sol[lamv]),
            }
        )
    return out


# ---------------------------------------------------------------------------
# T1: fixed s=-1 fibre + Tschirnhaus
# ---------------------------------------------------------------------------
def track_T1_sm1(seeds: list) -> list[dict]:
    """
    Known geometric fibre: p2=-1, σ=0, π=-1/25, c=-√5.
    Free: λ, c0, c1, c2. Solve Res = seed.
    """
    p2v, sigv, piv = -1, 0, sp.Rational(-1, 25)
    rows = []
    for tag, alpha, beta, k in seeds:
        print(f"  T1 {tag} ...", flush=True)
        F = fibre_expr(p2v, sigv, piv, lam)
        tau = c0 + c1 * y + c2 * y**2
        # resultant is polynomial in X with coeffs rational in (lam,c0,c1,c2)
        # degree 5 in y for F; τ deg 2 → Res deg 5 in X typically
        try:
            res_poly = tschirnhaus_resultant(F, tau)
        except Exception as e:
            rows.append({"tag": tag, "ok": False, "error": f"res:{e}"})
            continue
        if res_poly is None or res_poly.degree() != 5:
            rows.append(
                {
                    "tag": tag,
                    "ok": False,
                    "error": f"bad_deg:{None if res_poly is None else res_poly.degree()}",
                }
            )
            continue
        eqs = match_seed_coeffs(res_poly, alpha, beta)
        # Prefer c1 = 1 (affine scaling of τ target can be absorbed into X-scaling
        # only if seed scales — BJ monic forbids X-scale). Keep c1 free but try c1=1 first.
        attempts = []
        for ansatz_name, subs0, free in [
            ("c1=1", {c1: 1}, [lam, c0, c2]),
            ("c2=0_c1=1", {c1: 1, c2: 0}, [lam, c0]),  # translation only on fixed fibre
            ("full", {}, [lam, c0, c1, c2]),
        ]:
            eqs_a = [sp.simplify(e.subs(subs0)) for e in eqs]
            # drop identical zeros
            eqs_a = [e for e in eqs_a if e != 0]
            try:
                sols = sp.solve(eqs_a, free, dict=True, simplify=False)
            except Exception as e:
                attempts.append({"ansatz": ansatz_name, "ok": False, "error": str(e)[:120]})
                continue
            good = []
            for sol in sols or []:
                if any(v.has(sp.oo) or v == sp.zoo for v in sol.values()):
                    continue
                # merge ansatz
                full = dict(subs0)
                full.update({k: sp.simplify(v) for k, v in sol.items()})
                # residual check
                chk = [sp.simplify(e.subs(full)) for e in eqs]
                if all(c == 0 for c in chk):
                    good.append({str(k): str(v) for k, v in full.items()})
            attempts.append(
                {
                    "ansatz": ansatz_name,
                    "ok": len(good) > 0,
                    "n_sols": len(good),
                    "sample": good[:3],
                }
            )
            if good:
                break  # first success sufficient
        ok = any(a.get("ok") for a in attempts)
        rows.append(
            {
                "tag": tag,
                "k": str(k),
                "alpha": alpha,
                "beta": beta,
                "ok": ok,
                "fibre": "s=-1 fixed (p2,σ,π)=(-1,0,-1/25)",
                "attempts": attempts,
            }
        )
        print(f"    ok={ok} attempts={[a['ansatz']+':'+str(a.get('ok')) for a in attempts]}", flush=True)
    return rows


# ---------------------------------------------------------------------------
# T2: locus reverse with Tschirnhaus (numeric)
# ---------------------------------------------------------------------------
def cover_params_from_qw(qv, wv):
    p2v = physical_p2_at(qv, wv)
    if p2v is None:
        return None
    ch = chart(p2v, qv)
    if ch is None:
        return None
    chw = chart(p2v, wv)
    if chw is None:
        return None
    if sp.simplify(ch["sigma"] - chw["sigma"]) != 0 or sp.simplify(ch["pi"] - chw["pi"]) != 0:
        # try other p2 branch
        f1 = sp.simplify(F1.subs({q: qv, w: wv}))
        for sol in sp.solve(sp.Eq(f1, 0), p2):
            p2v = sp.simplify(sol)
            ch = chart(p2v, qv)
            chw = chart(p2v, wv)
            if ch is None or chw is None:
                continue
            if sp.simplify(ch["sigma"] - chw["sigma"]) == 0 and sp.simplify(ch["pi"] - chw["pi"]) == 0:
                break
        else:
            return None
    sv = s_form(p2v, qv, wv)
    return {
        "p2": p2v,
        "sigma": ch["sigma"],
        "pi": ch["pi"],
        "c": ch["c"],
        "s": sv,
        "q": qv,
        "w": wv,
    }


def track_T2_numeric(seeds: list, n_starts: int = 40) -> list[dict]:
    """
    Residual system on real variables:
      P(q,w)=0
      chart consistency already via physical p2
      monic Tschirnhaus Res coeffs match seed
    Free: q,w,lam,c0,c2 with c1=1.
    """
    rows = []
    # Build residual function numerically for each seed
    # Pre-build symbolic residual with symbols
    qv, wv, lamv, c0v, c2v = sp.symbols("qv wv lamv c0v c2v", real=True)
    # Use complex-friendly symbols without real assumption for broader nsolve
    qs, ws, lams, c0s, c2s = sp.symbols("qs ws lams c0s c2s")

    for tag, alpha, beta, k in seeds:
        print(f"  T2 {tag} ...", flush=True)
        # For speed: sample points near known geometric point and random
        hits = []
        best_res = None
        # Known geometric start
        rt5 = float(np.sqrt(5))
        starts = [
            (1 / rt5, -1 / rt5, 1.0, 0.0, 0.0),
            (-1 / rt5, 1 / rt5, 1.0, 0.0, 0.0),
            (1 / rt5, -1 / rt5, -1.0, 0.0, 0.1),
            (1 / rt5, -1 / rt5, 2.0, 0.5, -0.1),
        ]
        rng = np.random.default_rng(539 + abs(alpha) + abs(beta))
        for _ in range(n_starts):
            starts.append(
                (
                    rng.normal(0, 1.2),
                    rng.normal(0, 1.2),
                    rng.normal(0, 3),
                    rng.normal(0, 1),
                    rng.normal(0, 0.5),
                )
            )

        def residual(vec, alpha=alpha, beta=beta):
            qn, wn, ln, c0n, c2n = [float(v) for v in vec]
            # P
            Pn = float(
                20 * qn**3 * wn**3
                - 40 * qn**3 * wn**2
                + 27 * qn**3 * wn
                - 6 * qn**3
                - 40 * qn**2 * wn**3
                + 73 * qn**2 * wn**2
                - 45 * qn**2 * wn
                + 9 * qn**2
                + 27 * qn * wn**3
                - 45 * qn * wn**2
                + 26 * qn * wn
                - 5 * qn
                - 6 * wn**3
                + 9 * wn**2
                - 5 * wn
                + 1
            )
            # p2 from F1 quadratic
            # F1 = A p2^2 + B p2 + C
            A = 16 * qn * wn - 8 * qn - 8 * wn + 3
            B = (
                -30 * qn**2 * wn
                + 15 * qn**2
                - 30 * qn * wn**2
                + 37 * qn * wn
                - 8 * qn
                + 15 * wn**2
                - 8 * wn
            )
            C = 50 * qn**2 * wn**2 - 30 * qn**2 * wn - 30 * qn * wn**2 + 16 * qn * wn
            disc = B * B - 4 * A * C
            if abs(A) < 1e-14:
                return np.ones(6) * 1e3
            if disc < -1e-8:
                # allow tiny negative as 0
                if disc < -1e-3:
                    return np.array([Pn, 1e2, 1e2, 1e2, 1e2, 1e2])
                disc = 0.0
            sqrtD = np.sqrt(max(disc, 0.0))
            # pick branch nearer -1
            p2a = (-B + sqrtD) / (2 * A)
            p2b = (-B - sqrtD) / (2 * A)
            p2n = p2a if abs(p2a + 1) <= abs(p2b + 1) else p2b
            den_q = 6 * p2n * qn - 3 * p2n - 10 * qn**2 + 6 * qn
            den_w = 6 * p2n * wn - 3 * p2n - 10 * wn**2 + 6 * wn
            if abs(den_q) < 1e-12 or abs(den_w) < 1e-12 or abs(qn) < 1e-12:
                return np.ones(6) * 1e3
            sign = (
                qn * (8 * p2n * qn - 3 * p2n - 15 * qn**2 + 8 * qn) / den_q
            )
            pin = (
                qn**2 * (3 * p2n * qn - p2n - 6 * qn**2 + 3 * qn) / den_q
            )
            # chart consistency residual with w
            sigw = wn * (8 * p2n * wn - 3 * p2n - 15 * wn**2 + 8 * wn) / den_w
            piw = wn**2 * (3 * p2n * wn - p2n - 6 * wn**2 + 3 * wn) / den_w
            # Build F coeffs and Tschirnhaus numerically via companion resultant
            # F = y^5 - (1+p2)y^4 + p2 y^3 - ln y^2 + ln*sig y - ln*pi
            # τ = c0 + y + c2 y^2  (c1=1)
            # Res via sympy on floats is slow; use numpy poly resultant
            try:
                Fy = [
                    1.0,
                    -(1.0 + p2n),
                    p2n,
                    -ln,
                    ln * sign,
                    -ln * pin,
                ]  # high to low
                # resultant of F(y) and (c2 y^2 + y + c0) - X  as poly in y, then monic in X
                # For fixed numeric params, compute monic poly in X by evaluating
                # Sylvester resultant coefficients via sympy once... use np
                # X is parameter: g(y) = c2 y^2 + y + c0 - X
                # Res(F,g) is deg 5 in X for deg g = 2? Actually deg_y g = 2, deg F = 5,
                # Res has degree 2 in coeffs of F and degree 5 in coeffs of g.
                # Leading in X: g = -X + ..., Res degree 5 in X.
                # Compute via eigenvalue / companion of quadratic factors — use sympy N
            except Exception:
                return np.ones(6) * 1e3

            # Use complex root method: for monic F, roots r_i; image τ(r_i); elementary symmetric → poly
            try:
                roots = np.roots(Fy)
                images = c0n + roots + c2n * roots**2
                # monic poly with those roots
                mon = np.poly(images)  # length 6: X^5 + ...
                # want mon = [1,0,0,0,alpha,beta]
                target = np.array([1.0, 0.0, 0.0, 0.0, float(alpha), float(beta)])
                # residual of coeffs (skip leading 1)
                coeff_res = mon[1:] - target[1:]
                # also P and chart mismatch
                return np.array(
                    [
                        Pn,
                        sign - sigw,
                        pin - piw,
                        coeff_res[0],
                        coeff_res[1],
                        coeff_res[2],
                        coeff_res[3],
                        coeff_res[4],
                    ],
                    dtype=float,
                )
            except Exception:
                return np.ones(8) * 1e3

        for x0 in starts:
            try:
                sol = least_squares_numpy(
                    residual,
                    np.array(x0, dtype=float),
                    max_nfev=250,
                )
                rnorm = float(np.linalg.norm(sol.fun))
                if best_res is None or rnorm < best_res:
                    best_res = rnorm
                if rnorm < 1e-6:
                    qn, wn, ln, c0n, c2n = sol.x
                    hits.append(
                        {
                            "q": float(qn),
                            "w": float(wn),
                            "lam": float(ln),
                            "c0": float(c0n),
                            "c1": 1.0,
                            "c2": float(c2n),
                            "residual_norm": rnorm,
                        }
                    )
                    break
            except Exception:
                continue

        rows.append(
            {
                "tag": tag,
                "k": str(k),
                "alpha": alpha,
                "beta": beta,
                "ok": len(hits) > 0,
                "n_hits": len(hits),
                "best_residual": best_res,
                "hits": hits[:5],
                "method": "numpy_LM",
            }
        )
        print(
            f"    ok={len(hits)>0} hits={len(hits)} best_res={best_res}",
            flush=True,
        )
    return rows


# ---------------------------------------------------------------------------
# T3: forward numeric covers + Tschirnhaus → BJ → catalogue
# ---------------------------------------------------------------------------
def newton_cover(s_val: float, n_trials: int = 50):
    """Solve triple-root system for cover params at fixed s (from build_3a4_resolvent)."""

    def residual(v):
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

    def newton(x0, niter=60):
        v = np.array(x0, dtype=float)
        for _ in range(niter):
            r = residual(v)
            if np.linalg.norm(r) < 1e-12:
                return v, True, float(np.linalg.norm(r))
            J = np.zeros((6, 6))
            eps = 1e-8
            for j in range(6):
                dv = np.zeros(6)
                dv[j] = eps
                J[:, j] = (residual(v + dv) - r) / eps
            try:
                step = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -r, rcond=None)[0]
            v = v + step
        nr = float(np.linalg.norm(residual(v)))
        return v, nr < 1e-10, nr

    rng = np.random.default_rng(abs(hash(round(s_val, 8))) % (2**32))
    candidates = [rng.normal(size=6) for _ in range(n_trials)]
    candidates.append(np.array([-np.sqrt(5), -1.0, 0.2, -0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5]))
    best = None
    for x0 in candidates:
        v, ok, nr = newton(x0)
        if best is None or nr < best[2]:
            best = (v, ok, nr)
        if ok:
            return v, True, nr
    return best[0], False, best[2]


def fibre_to_bj_translation(coeffs_high_to_low):
    """Kill y^4 by translation; return (form, data)."""
    # monic [1,c4,c3,c2,c1,c0]
    c = [float(x) for x in coeffs_high_to_low]
    if abs(c[0] - 1.0) > 1e-8:
        c = [ci / c[0] for ci in c]
    shift = -c[1] / 5.0
    # expand (z+shift)^5 + c4 (z+shift)^4 + ...
    # use numpy poly composition
    # f(y)=sum c[i] y^{5-i}; y=z+s
    p = np.poly1d(c)
    psh = np.poly1d(p.c)  # copy
    # evaluate poly at z+shift: use poly transform
    z = np.poly1d([1, shift])  # y = z + shift
    # Horner
    res = np.poly1d([0.0])
    for coef in p.c:
        res = np.polymul(res, z) + np.poly1d([coef])
    cc = list(res.c)
    while len(cc) < 6:
        cc = [0.0] + cc
    cc = cc[-6:]
    # expect c4 ~ 0
    if abs(cc[1]) > 1e-6:
        return None
    if abs(cc[2]) < 1e-7 and abs(cc[3]) < 1e-7:
        a, b = cc[4], cc[5]
        # round near integers
        if abs(a - round(a)) < 1e-5 and abs(b - round(b)) < 1e-5:
            return ("BJ", int(round(a)), int(round(b)))
        return ("BJ_float", float(a), float(b))
    return ("depressed", float(cc[2]), float(cc[3]), float(cc[4]), float(cc[5]))


def tschirnhaus_numeric_to_bj(coeffs_high_to_low, n_starts=25):
    """
    Find c0,c2 (c1=1) so images of roots under τ=c0+y+c2 y^2 form a BJ poly
    (no X^4,X^3,X^2 terms). Return best BJ-like integer match if any.
    """
    c = [float(x) for x in coeffs_high_to_low]
    if abs(c[0] - 1.0) > 1e-8:
        c = [ci / c[0] for ci in c]
    try:
        roots = np.roots(c)
    except Exception:
        return None, None
    if np.any(np.isnan(roots)):
        return None, None

    def obj(v):
        c0n, c2n = v
        imgs = c0n + roots + c2n * roots**2
        mon = np.poly(imgs)
        # mon[1], mon[2], mon[3] should vanish for BJ
        return np.array([mon[1], mon[2], mon[3]], dtype=float)

    best = None
    rng = np.random.default_rng(7)
    starts = [(0.0, 0.0), (0.0, 0.1), (0.5, 0.0), (-0.5, 0.1), (0.0, -0.2), (1.0, 0.05)]
    starts += [(rng.normal(), rng.normal() * 0.3) for _ in range(n_starts)]

    for x0 in starts:
        try:
            sol = least_squares_numpy(obj, np.array(x0, float), max_nfev=100)
            rnorm = float(np.linalg.norm(sol.fun))
            mon = np.poly(sol.x[0] + roots + sol.x[1] * roots**2)
            if best is None or rnorm < best[0]:
                best = (rnorm, sol.x, mon)
            if rnorm < 1e-8:
                break
        except Exception:
            continue

    if best is None:
        return None, None
    rnorm, xv, mon = best
    if rnorm > 1e-6:
        return {"ok": False, "residual": rnorm}, None
    # mon ~ [1,0,0,0,a,b]
    a, b = float(mon[4]), float(mon[5])
    if abs(a - round(a)) < 1e-4 and abs(b - round(b)) < 1e-4:
        return {
            "ok": True,
            "residual": rnorm,
            "c0": float(xv[0]),
            "c2": float(xv[1]),
            "alpha": int(round(a)),
            "beta": int(round(b)),
        }, (int(round(a)), int(round(b)))
    return {
        "ok": True,
        "residual": rnorm,
        "c0": float(xv[0]),
        "c2": float(xv[1]),
        "alpha": a,
        "beta": b,
        "integer": False,
    }, (a, b)


def track_T3_forward() -> dict:
    s_list = [
        -3, -2, -1, -0.5, 0.5, 1.5, 2, 2.5, 3, 4, 5, -2.5, 0.333, 0.666, 1.333, -1.5
    ]
    t_list = list(range(-8, 9)) + [0.5, 1.5, 2.5, -0.5]
    covers_ok = 0
    bj_trans = []
    bj_tsch = []
    cat_hits_trans = []
    cat_hits_tsch = []

    for s_val in s_list:
        v, ok, nr = newton_cover(float(s_val), n_trials=40)
        if not ok:
            continue
        covers_ok += 1
        c, p2n, r1, r2, qn, wn = v
        # monic fibre: y^5-(1+p2)y^4+p2 y^3 - (t/c)(y^2-(r1+r2)y+r1 r2)
        sig = r1 + r2
        pin = r1 * r2
        for tv in t_list:
            if abs(c) < 1e-12:
                continue
            ratio = float(tv) / float(c)
            coeffs = [
                1.0,
                -(1.0 + p2n),
                p2n,
                -ratio,
                ratio * sig,
                -ratio * pin,
            ]
            # translation BJ
            tr = fibre_to_bj_translation(coeffs)
            if tr and tr[0] == "BJ":
                _, a, b = tr
                bj_trans.append({"s": s_val, "t": tv, "alpha": a, "beta": b})
                if (a, b) in CAT_BY_AB:
                    tag, kk = CAT_BY_AB[(a, b)]
                    cat_hits_trans.append(
                        {"tag": tag, "k": str(kk), "s": s_val, "t": tv, "alpha": a, "beta": b, "via": "translation"}
                    )
            # Tschirnhaus BJ
            info, ab = tschirnhaus_numeric_to_bj(coeffs, n_starts=12)
            if info and info.get("ok") and ab and info.get("integer", True) is not False:
                if isinstance(ab[0], int):
                    a, b = ab
                    bj_tsch.append({"s": s_val, "t": tv, "alpha": a, "beta": b, "tsch": info})
                    if (a, b) in CAT_BY_AB:
                        tag, kk = CAT_BY_AB[(a, b)]
                        cat_hits_tsch.append(
                            {
                                "tag": tag,
                                "k": str(kk),
                                "s": s_val,
                                "t": tv,
                                "alpha": a,
                                "beta": b,
                                "via": "tschirnhaus",
                                "c0": info.get("c0"),
                                "c2": info.get("c2"),
                            }
                        )
        print(f"  T3 cover s={s_val}: ok res={nr:.1e}", flush=True)

    cat_k = sorted({h["k"] for h in cat_hits_trans + cat_hits_tsch})
    return {
        "covers_ok": covers_ok,
        "covers_attempted": len(s_list),
        "bj_translation": len(bj_trans),
        "bj_tschirnhaus": len(bj_tsch),
        "cat_hits_translation": cat_hits_trans,
        "cat_hits_tschirnhaus": cat_hits_tsch,
        "catalogue_k": cat_k,
        "multi_k": len(cat_k) >= 2,
        "bj_trans_sample": bj_trans[:15],
        "bj_tsch_sample": bj_tsch[:15],
    }


# ---------------------------------------------------------------------------
# Exact symbolic T1 lite for translation on s=-1 (sanity)
# ---------------------------------------------------------------------------
def track_T1_exact_light(seeds) -> list[dict]:
    """
    On s=-1 fibre, try quadratic Tschirnhaus with c1=1 using sympy solve
    on a reduced system — only for first few seeds if full T1 is heavy.
    Uses numerical verification of symbolic residual structure.
    """
    # Already covered by track_T1_sm1; this is a backup using root-image integer search
    rows = []
    p2v, sigv, piv = -1.0, 0.0, -0.04
    for tag, alpha, beta, k in seeds[:4]:
        hits = []
        for lamn in np.linspace(-5, 5, 41):
            if abs(lamn) < 1e-12:
                continue
            coeffs = [1.0, -(1 + p2v), p2v, -lamn, lamn * sigv, -lamn * piv]
            info, ab = tschirnhaus_numeric_to_bj(coeffs, n_starts=20)
            if info and info.get("ok") and ab == (alpha, beta):
                hits.append({"lam": float(lamn), **info})
        rows.append(
            {
                "tag": tag,
                "ok": len(hits) > 0,
                "n_hits": len(hits),
                "hits": hits[:3],
                "note": "numeric lam-scan on s=-1 + Tschirnhaus",
            }
        )
        print(f"  T1b {tag}: hits={len(hits)}", flush=True)
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G1 — seed-first + domain Tschirnhaus", flush=True)
    print("=" * 72, flush=True)

    seeds = CATALOGUE  # 10 priority multi-ratio seeds

    print("\n[T1] Fixed s=-1 fibre + symbolic Tschirnhaus solve ...", flush=True)
    t1 = track_T1_sm1(seeds)
    t1_ok = sum(1 for r in t1 if r.get("ok"))

    print("\n[T1b] Numeric lam-scan on s=-1 + Tschirnhaus (first 4 seeds) ...", flush=True)
    t1b = track_T1_exact_light(seeds)

    print("\n[T2] Locus-constrained numeric reverse (P=0 + Tschirnhaus) ...", flush=True)
    t2 = track_T2_numeric(seeds, n_starts=30)
    t2_ok = sum(1 for r in t2 if r.get("ok"))

    print("\n[T3] Forward Newton covers + Tschirnhaus → catalogue ...", flush=True)
    t3 = track_T3_forward()

    # Baseline reminder: translation reverse always works, locus never (G1)
    print("\n[Base] Translation-only reverse count (sanity) ...", flush=True)
    base = []
    for tag, a, b, k in seeds:
        rev = translation_only_reverse(a, b)
        base.append({"tag": tag, "n_rev": len(rev), "ok": len(rev) > 0})
    print(f"  translation reverse ok: {sum(1 for r in base if r['ok'])}/{len(base)}", flush=True)

    elapsed = round(time.time() - t0, 2)
    multi = bool(t3.get("multi_k")) or (
        len({r["k"] for r in t2 if r.get("ok")}) >= 2
    ) or (
        len({r["k"] for r in t1 if r.get("ok")}) >= 2
    )
    geometric_hit = t1_ok > 0 or t2_ok > 0 or len(t3.get("cat_hits_tschirnhaus", [])) > 0 or len(
        t3.get("cat_hits_translation", [])
    ) > 0

    verdict = (
        f"G1 Tschirnhaus cut ({elapsed}s). "
        f"T1 s=-1 exact/symbolic hits={t1_ok}/{len(t1)}. "
        f"T2 locus numeric hits={t2_ok}/{len(t2)}. "
        f"T3 covers_ok={t3['covers_ok']}/{t3['covers_attempted']}; "
        f"BJ_trans={t3['bj_translation']}, BJ_tsch={t3['bj_tschirnhaus']}; "
        f"cat_trans={len(t3['cat_hits_translation'])}, cat_tsch={len(t3['cat_hits_tschirnhaus'])}; "
        f"multi_k={multi}; geometric_catalogue_hit={geometric_hit}."
    )
    print("\n" + verdict, flush=True)

    # Report
    lines = [
        "# G1 — seed-first + domain Tschirnhaus",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Goal",
        "",
        "Enlarge the reverse of catalogue BJ seeds through the N−tD fibre by a domain / fibre",
        "Tschirnhaus transform",
        "",
        "```",
        "τ(y) = c0 + c1 y + c2 y²",
        "Res_y(F(y), X − τ(y)) = X⁵ + α X + β",
        "```",
        "",
        "with F the monic N−tD fibre and (p2,σ,π) on the triple-root locus P(q,w)=0.",
        "",
        "G1 baseline: translation-only reverse always works; locus hit rate 0/16.",
        "This cut tests whether quadratic Tschirnhaus freedom creates locus hits.",
        "Numeric tracks use pure-NumPy Levenberg–Marquardt least squares (no scipy).",
        "",
        "---",
        "",
        "## 1. Track T1 — fixed geometric fibre s=−1",
        "",
        "Locked params: (p2,σ,π)=(−1, 0, −1/25). Free: (λ,c0,c1,c2).",
        "",
        f"| seed | k | hit? | notes |",
        f"|------|---|:----:|-------|",
    ]
    for r in t1:
        notes = ", ".join(
            f"{a['ansatz']}:{a.get('ok')}" for a in r.get("attempts", [])
        )
        lines.append(f"| {r['tag']} | {r.get('k')} | {r.get('ok')} | {notes} |")
    lines += [
        "",
        f"**T1 hits: {t1_ok}/{len(t1)}**",
        "",
        "### T1b — numeric λ-scan on s=−1 (first 4 seeds)",
        "",
    ]
    for r in t1b:
        lines.append(f"- {r['tag']}: hits={r['n_hits']} ok={r['ok']}")
    lines += [
        "",
        "---",
        "",
        "## 2. Track T2 — locus-constrained numeric reverse",
        "",
        "Free (q,w,λ,c0,c2) with c1=1, residual enforces P(q,w)=0, chart consistency,",
        "and monic Tschirnhaus image = seed. Solver: scipy.least_squares multi-start",
        "(falls back to residual scan if scipy missing).",
        "",
        f"| seed | k | hit? | best residual | #hits |",
        f"|------|---|:----:|--------------:|------:|",
    ]
    for r in t2:
        br = r.get("best_residual")
        brs = f"{br:.2e}" if isinstance(br, float) else str(br)
        lines.append(
            f"| {r['tag']} | {r['k']} | {r['ok']} | {brs} | {r['n_hits']} |"
        )
    lines += [
        "",
        f"**T2 hits: {t2_ok}/{len(t2)}**",
        "",
    ]
    if t2_ok:
        lines.append("### T2 hit details")
        lines.append("")
        for r in t2:
            if not r["ok"]:
                continue
            lines.append(f"- **{r['tag']}**: {r['hits'][:2]}")
        lines.append("")
    lines += [
        "---",
        "",
        "## 3. Track T3 — forward covers + Tschirnhaus",
        "",
        "Newton-solve (3,1,1)⁴ covers at rational s; specialise t; reduce by translation",
        "and by quadratic Tschirnhaus; match multi-seed catalogue.",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| covers ok | {t3['covers_ok']}/{t3['covers_attempted']} |",
        f"| BJ via translation | {t3['bj_translation']} |",
        f"| BJ via Tschirnhaus | {t3['bj_tschirnhaus']} |",
        f"| catalogue hits (translation) | {len(t3['cat_hits_translation'])} |",
        f"| catalogue hits (Tschirnhaus) | {len(t3['cat_hits_tschirnhaus'])} |",
        f"| catalogue k set | {t3['catalogue_k']} |",
        f"| multi-k | {t3['multi_k']} |",
        "",
    ]
    if t3["cat_hits_translation"] or t3["cat_hits_tschirnhaus"]:
        lines.append("### Catalogue hits")
        lines.append("")
        for h in t3["cat_hits_translation"] + t3["cat_hits_tschirnhaus"]:
            lines.append(f"- {h}")
        lines.append("")
    else:
        lines.append("_No catalogue seeds recovered from forward T3 fibres._")
        lines.append("")
    lines += [
        "---",
        "",
        "## 4. Multi-k conclusion",
        "",
        f"| test | result |",
        f"|------|--------|",
        f"| T1 s=−1 + Tschirnhaus seed hit | **{t1_ok > 0}** ({t1_ok}/{len(t1)}) |",
        f"| T2 locus + Tschirnhaus seed hit | **{t2_ok > 0}** ({t2_ok}/{len(t2)}) |",
        f"| T3 forward catalogue hit | **{geometric_hit and (len(t3['cat_hits_translation'])+len(t3['cat_hits_tschirnhaus']))>0}** |",
        f"| Geometric multi-k | **{multi}** |",
        f"| Translation-only reverse (baseline) | **{sum(1 for r in base if r['ok'])}/{len(base)}** |",
        "",
        "**Interpretation.**",
        "",
        "- If T1/T2/T3 stay empty: quadratic fibre Tschirnhaus is **not** enough to place",
        "  pure-even catalogue seeds on this 3A⁴ normal-form locus — the obstruction is",
        "  deeper than the y⁴-translation reverse used in G1.",
        "- If hits appear: record (q,w,s,λ,τ) as the first geometric multi-k candidates",
        "  and verify exactly over Q or Q(√5).",
        "",
        "### Next if still empty",
        "",
        "1. Full cubic Tschirnhaus τ = c0+c1 y+c2 y²+c3 y³ (classical quintic freedom).",
        "2. Domain Möbius before normal form (move {0,1,∞} labels) — true domain freedom.",
        "3. Parameter-field resolvent over multi-sheeted K/Q(s) (G1 remaining item 2).",
        "4. G2 other Nielsen types.",
        "",
        "---",
        "",
        "## 5. Non-claims",
        "",
        "- Does not alter pure-even multi-k arithmetic or Canonical T3.",
        "- Negative results are for quadratic Tschirnhaus + this normal form, not a proof",
        "  that geometric multi-k is impossible.",
        "",
        "_Generated by `g1_seed_tschirnhaus.py`._",
        "",
    ]

    md = "\n".join(lines)
    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "T1": t1,
        "T1b": t1b,
        "T2": t2,
        "T3": t3,
        "baseline_translation_reverse": base,
        "t1_ok": t1_ok,
        "t2_ok": t2_ok,
        "geometric_hit": geometric_hit,
        "multi_k": multi,
    }

    write_md(ROOT / "G1_SEED_TSCHIRNHAUS.md", md)
    write_json(ROOT / "G1_SEED_TSCHIRNHAUS.json", payload)
    write_md(OUT / "G1_SEED_TSCHIRNHAUS.md", md)
    write_json(OUT / "G1_SEED_TSCHIRNHAUS.json", payload)
    try:
        write_md(RESULTS / "G1_SEED_TSCHIRNHAUS.md", md)
        write_json(RESULTS / "G1_SEED_TSCHIRNHAUS.json", payload)
    except Exception:
        pass

    print(f"\nWrote G1_SEED_TSCHIRNHAUS.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
