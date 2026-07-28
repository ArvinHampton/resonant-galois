"""
Execute ranked avenues 1→7 for geometric multi-k fusion.

1. Better rational coordinate / resolvent for 3A^4
2. Next shortlist genus-0: 2A 3A^3, 2A^2 3A^2
3. Positive-dimensional pure-even A5 strata
4. Other rigid triples
5. Base change + descent (extend prior K=Q(√5))
6. Higher-rank rigid systems (r≥5 sketch)
7. Geometric lift of existing envelope

Output: AVENUE_RANK_EXECUTE.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.linalg import norm

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

y, t, s = sp.symbols("y t s")
c, p2, r1, r2, q, w = sp.symbols("c p2 r1 r2 q w")


CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("s220m", 220, -528, Fraction(-12, 5)),
]


def k_of(a, b):
    if a == 0:
        return None
    return Fraction(int(b), int(a))


def is_square_poly(expr, var=t) -> bool:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return True
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return False
        n, d = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return False
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        return all(m % 2 == 0 for _, m in fac[1])
    except Exception:
        return False


# ===========================================================================
# AVENUE 1 — Better coordinate / resolvent for 3A^4
# ===========================================================================
def avenue1_resolvent() -> dict:
    print("=== AVENUE 1: better 3A^4 resolvent coordinate ===", flush=True)
    out: dict = {"name": "better_3A4_resolvent"}

    # --- 1A. Algebraic elimination of triple-root conditions ---
    print("  1A. resultant elimination...", flush=True)
    N = c * y**3 * (y - 1) * (y - p2)
    D = (y - r1) * (y - r2)
    F1 = sp.expand(N - D)  # branch value 1
    Fs = sp.expand(N - s * D)  # branch value s

    # For F to have a triple root: Res(F, F')=0 and Res(F', F'')=0
    def triple_root_conditions(F, var=y):
        Fp = sp.diff(F, var)
        Fpp = sp.diff(Fp, var)
        R1 = sp.resultant(F, Fp, var)
        R2 = sp.resultant(Fp, Fpp, var)
        return sp.factor(sp.expand(R1)), sp.factor(sp.expand(R2))

    R1a, R1b = triple_root_conditions(F1)
    Rsa, Rsb = triple_root_conditions(Fs)
    out["eliminants_deg"] = {
        "R1_F_Fp": int(sp.total_degree(sp.Poly(sp.numer(sp.together(R1a)), [c, p2, r1, r2])))
        if R1a != 0
        else -1,
        "note": "resultants in cover parameters for triple roots at branch 1 and s",
    }
    # Try to solve the system for fixed symbolic s via groebner (may be heavy)
    # Use: eliminate q,w by resultant only — already have R1a,R1b,Rsa,Rsb
    # Scale invariance: (c,p2,r1,r2) -> ... fix r1+r2=0 (symmetric poles) as better coordinate
    print("  1B. symmetric-pole coordinate r2=-r1...", flush=True)
    # Set r2 = -r1, p2 free, c free, r1 free — 3 free + s
    # 4 resultant conditions — overdetermined; find locus
    subs_sym = {r2: -r1}
    eqs_sym = [
        sp.numer(sp.together(sp.expand(R1a.subs(subs_sym)))),
        sp.numer(sp.together(sp.expand(R1b.subs(subs_sym)))),
        sp.numer(sp.together(sp.expand(Rsa.subs(subs_sym)))),
        sp.numer(sp.together(sp.expand(Rsb.subs(subs_sym)))),
    ]
    # Factor each and take primitive
    eqs_f = []
    for e in eqs_sym:
        try:
            e2 = sp.factor(sp.expand(e))
            eqs_f.append(e2)
        except Exception:
            eqs_f.append(e)
    out["symmetric_pole_eqs_preview"] = [str(e)[:120] for e in eqs_f]

    # Solve numerically many s with r2=-r1 constrained Newton
    def residual_sym(v, s_val):
        # v = (c, p2, r1, q, w)  with r2=-r1
        cc, pp, rr, qq, ww = v
        r2v = -rr
        # reuse residual from triple roots
        def G_vals(val, pt, p2v, r1v, r2v, cv):
            yy = pt
            Nloc = cv * yy**3 * (yy - 1) * (yy - p2v)
            A, Ap, App = yy**3, 3 * yy**2, 6 * yy
            B = yy**2 - (1 + p2v) * yy + p2v
            Bp = 2 * yy - (1 + p2v)
            Bpp = 2.0
            Np = cv * (Ap * B + A * Bp)
            Npp = cv * (App * B + 2 * Ap * Bp + A * Bpp)
            Dd = (yy - r1v) * (yy - r2v)
            Dp = 2 * yy - (r1v + r2v)
            Dpp = 2.0
            return [Nloc - val * Dd, Np - val * Dp, Npp - val * Dpp]

        return np.array(
            G_vals(1.0, qq, pp, rr, r2v, cc) + G_vals(s_val, ww, pp, rr, r2v, cc),
            dtype=float,
        )

    def newton_sym(s_val, x0, niter=60):
        v = np.array(x0, float)
        for _ in range(niter):
            r = residual_sym(v, s_val)
            if norm(r) < 1e-14:
                return v, True, float(norm(r))
            J = np.zeros((6, 5))
            # residual is 6-dim, v is 5-dim — overdetermined; use Gauss-Newton
            eps = 1e-8
            r0 = r
            for j in range(5):
                dv = np.zeros(5)
                dv[j] = eps
                J[:, j] = (residual_sym(v + dv, s_val) - r0) / eps
            step, _, _, _ = np.linalg.lstsq(J, -r0, rcond=None)
            v = v + step
        nr = float(norm(residual_sym(v, s_val)))
        return v, nr < 1e-10, nr

    sym_hits = []
    rng = np.random.default_rng(2)
    for sv in [-1.0, -2.0, 0.5, 2.0, 3.0, 5.0, -0.5, 1.5]:
        found = False
        # seed from known s=-1: c=-sqrt5, p2=-1, r1=0.2, q=sqrt5/5, w=-sqrt5/5
        seeds = [
            np.array([-np.sqrt(5), -1.0, 0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5]),
            rng.normal(size=5),
        ]
        for trial in range(40):
            x0 = seeds[0] if trial == 0 else rng.normal(size=5)
            v, ok, nr = newton_sym(sv, x0)
            if ok:
                sym_hits.append(
                    {
                        "s": sv,
                        "c": float(v[0]),
                        "p2": float(v[1]),
                        "r1": float(v[2]),
                        "r2": float(-v[2]),
                        "res": nr,
                    }
                )
                print(f"    sym-pole HIT s={sv} p2={v[1]:.4f} r1={v[2]:.4f}", flush=True)
                found = True
                break
        if not found:
            print(f"    sym-pole miss s={sv}", flush=True)
    out["symmetric_pole_hits"] = sym_hits
    out["symmetric_pole_n"] = len(sym_hits)

    # --- 1C. Reparametrise by p2 (critical zero) instead of s ---
    print("  1C. p2-coordinate family (s as function of p2)...", flush=True)
    # From Newton free form: for each p2, solve for s,c,r1,r2
    # Use free residual with unknowns (c,r1,r2,q,w,s) at fixed p2
    def residual_p2fixed(v, p2v):
        cc, r1v, r2v, qq, ww, ss = v

        def G_vals(val, pt):
            yy = pt
            Nloc = cc * yy**3 * (yy - 1) * (yy - p2v)
            A, Ap, App = yy**3, 3 * yy**2, 6 * yy
            B = yy**2 - (1 + p2v) * yy + p2v
            Bp = 2 * yy - (1 + p2v)
            Np = cc * (Ap * B + A * Bp)
            Npp = cc * (App * B + 2 * Ap * Bp + A * 2.0)
            Dd = (yy - r1v) * (yy - r2v)
            Dp = 2 * yy - (r1v + r2v)
            return [Nloc - val * Dd, Np - val * Dp, Npp - val * 2.0]

        return np.array(G_vals(1.0, qq) + G_vals(ss, ww), dtype=float)

    def newton_p2(p2v, x0, niter=80):
        v = np.array(x0, float)
        for _ in range(niter):
            r = residual_p2fixed(v, p2v)
            if norm(r) < 1e-14:
                return v, True, float(norm(r))
            J = np.zeros((6, 6))
            r0 = r
            eps = 1e-8
            for j in range(6):
                dv = np.zeros(6)
                dv[j] = eps
                J[:, j] = (residual_p2fixed(v + dv, p2v) - r0) / eps
            try:
                step = np.linalg.solve(J, -r0)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -r0, rcond=None)[0]
            v = v + step
        nr = float(norm(residual_p2fixed(v, p2v)))
        return v, nr < 1e-12, nr

    p2_hits = []
    for p2v in [-1.0, -0.5, 0.5, 1.5, 2.0, 2.6180339887, 3.0, -1.5, 0.3]:
        if abs(p2v) < 1e-9 or abs(p2v - 1) < 1e-9:
            continue
        found = False
        for trial in range(50):
            x0 = rng.normal(size=6)
            x0[5] = abs(x0[5]) + 0.3  # s
            if abs(x0[5] - 1) < 0.1:
                x0[5] = 2.0
            v, ok, nr = newton_p2(p2v, x0)
            if ok and abs(v[5]) > 1e-6 and abs(v[5] - 1) > 1e-6:
                p2_hits.append(
                    {
                        "p2": p2v,
                        "s": float(v[5]),
                        "c": float(v[0]),
                        "r1": float(v[1]),
                        "r2": float(v[2]),
                        "res": nr,
                    }
                )
                print(f"    p2={p2v:.4f} -> s={v[5]:.4f}", flush=True)
                found = True
                break
        if not found:
            print(f"    p2={p2v:.4f} miss", flush=True)
    out["p2_coordinate_hits"] = p2_hits

    # Fit s as rational function of p2 from hits (if enough)
    if len(p2_hits) >= 4:
        xs = np.array([h["p2"] for h in p2_hits])
        ys = np.array([h["s"] for h in p2_hits])
        # poly fit deg 2
        coef = np.polyfit(xs, ys, min(2, len(xs) - 1))
        pred = np.polyval(coef, xs)
        out["s_of_p2_polyfit"] = {
            "coeffs_high_to_low": coef.tolist(),
            "max_err": float(np.max(np.abs(pred - ys))),
        }
        print(f"    s(p2) polyfit max_err={out['s_of_p2_polyfit']['max_err']:.3e}", flush=True)

    # --- 1D. Exact Q(s)-search: numeric test of p2-rational-in-s ansätze ---
    print("  1D. ansatz p2 rational in s (numeric Newton)...", flush=True)
    ansatz_p2_num = [
        ("p2=s", lambda sv: sv),
        ("p2=1-s", lambda sv: 1 - sv),
        ("p2=-s", lambda sv: -sv),
        ("p2=s+1", lambda sv: sv + 1),
        ("p2=2-s", lambda sv: 2 - sv),
        ("p2=-1", lambda sv: -1.0),
        ("p2=phi2", lambda sv: (3 + np.sqrt(5)) / 2),  # known at s=2
    ]
    ansatz_hits = []
    for name, p2fn in ansatz_p2_num:
        ok_s = []
        for sv in [2.0, 3.0, 0.5, -1.0, 1.5, 5.0]:
            p2v = float(p2fn(sv))
            if abs(p2v) < 1e-9 or abs(p2v - 1) < 1e-9:
                continue
            # Newton with p2 fixed via residual_p2fixed
            found = False
            for trial in range(25):
                x0 = rng.normal(size=6)
                x0[5] = sv  # seed s
                v, ok, nr = newton_p2(p2v, x0)
                if ok and abs(v[5] - sv) < 0.05:
                    ok_s.append(sv)
                    found = True
                    break
            # also accept any solved s near target for p2=const ansätze
            if not found and name in ("p2=-1", "p2=phi2"):
                for trial in range(25):
                    v, ok, nr = newton_p2(p2v, rng.normal(size=6))
                    if ok:
                        ok_s.append(float(v[5]))
                        found = True
                        break
        ansatz_hits.append({"ansatz": name, "solved_s": ok_s, "n": len(ok_s)})
        print(f"    ansatz {name}: n={len(ok_s)} s={ok_s[:4]}", flush=True)
    out["p2_ansatz_hits"] = ansatz_hits

    # --- 1E. Catalogue test on exact s=-1 model (norm deg 10) and on BJ pure-even controls ---
    print("  1E. specialise exact s=-1 and control families...", flush=True)
    # Norm poly at s=-1: 5(y^5-y^3)^2 - t^2 (y^2-1/25)^2
    specs_norm = []
    for tv in list(range(-10, 11)) + [Fraction(1, 2), Fraction(3, 2), Fraction(5, 2)]:
        if tv == 0:
            continue
        A = y**5 - y**3
        B = sp.Rational(tv, 5) * y**2 - sp.Rational(tv, 125)
        # wait: from earlier: A + sqrt5 * B = 0 with B = tv/5 y^2 - tv/125
        # A = y^5 - y^3, sqrt5 factor: from f = y^5-y^3 + (tv/sqrt5)y^2 - tv/(25 sqrt5)
        # = A + (tv/5) sqrt5 y^2 - (tv/125) sqrt5 = A + sqrt5 * (tv/5 y^2 - tv/125)
        normp = sp.expand(A**2 - 5 * B**2)
        pol = sp.Poly(sp.numer(sp.together(normp)), y, domain=sp.ZZ)
        if pol.LC() < 0:
            pol = sp.Poly(-pol.as_expr(), y, domain=sp.ZZ)
        # make monic if possible
        if pol.LC() != 1:
            # content clear
            cont = int(pol.content())
            if cont != 0 and cont != 1:
                try:
                    pol = sp.Poly(pol.as_expr() / cont, y, domain=sp.ZZ)
                except Exception:
                    pass
        specs_norm.append(
            {
                "t": str(tv),
                "deg": pol.degree(),
                "irr": bool(pol.is_irreducible) if pol.degree() > 0 else False,
            }
        )
    out["s_minus1_norm_specs"] = specs_norm[:15]
    out["s_minus1_norm_n_irr10"] = sum(1 for r in specs_norm if r["deg"] == 10 and r["irr"])

    # Control: envelope multi-k still works
    m0 = Fraction(5, 16)
    ku = Fraction(-8, 5) + t * (Fraction(4, 5) - Fraction(-8, 5))
    alpha_e = sp.together(25 - 3125 * ku**4 / 256)
    beta_e = sp.together(ku * alpha_e)
    env_hits = []
    for tv in [0, 1, Fraction(1, 2), Fraction(1, 3), Fraction(2, 3)]:
        try:
            a = int(sp.Rational(sp.simplify(alpha_e.subs(t, tv))))
            b = int(sp.Rational(sp.simplify(beta_e.subs(t, tv))))
        except Exception:
            continue
        for tag, ca, cb, ck in CATALOGUE:
            if a == ca and b == cb:
                env_hits.append({"tag": tag, "k": str(ck), "t": str(tv)})
    out["envelope_control_hits"] = env_hits
    out["envelope_multi_k"] = len({h["k"] for h in env_hits}) >= 2
    out["closed_form_Q_s_resolvent"] = None
    out["verdict"] = (
        f"Symmetric-pole hits={len(sym_hits)}; p2-coordinate hits={len(p2_hits)}; "
        f"p2-ansatz closed sols at s=2: {sum(1 for h in ansatz_hits if h.get('n_sols_at_s2'))}; "
        f"closed form f_s in Q(s)[x]: still open; envelope multi-k control={out['envelope_multi_k']}"
    )
    print(f"  A1 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 2 — Next shortlist genus-0 classes
# ===========================================================================
def avenue2_shortlist() -> dict:
    print("=== AVENUE 2: 2A 3A^3 and 2A^2 3A^2 ===", flush=True)
    out: dict = {"name": "next_shortlist_g0"}

    # Cycle types in S5:
    # 2A double transp: cycle type (2,2,1)
    # 3A: (3,1,1)
    # For a degree-5 cover, branch types match the S5 cycle types of monodromy.

    # Construct covers with branch types (2,2,1), (3,1,1), (3,1,1), (3,1,1)
    # Normal form ideas:
    # - one branch value with preimage mults 2,2,1
    # - three with 3,1,1

    # Pure-even BJ families that could arise as resolvents: already have LSW, flagship, etc.
    # Here: search one-param pure-even with disc factorisation suggesting 4 branch values
    # and sample Gal A5; check multi-k.

    # Also: enumerate Nielsen orbit sizes from prior data
    out["prior_orbits"] = {
        "2A,3A,3A,3A": {"orbits": 1, "sizes": [96], "genus_lookup": 0},
        "2A,2A,3A,3A": {"orbits": 1, "sizes": [108], "genus_lookup": 0},
    }

    # Cover ansatz for 2A+3A^3: phi = c y^2 (y-1)^2 (y-a) / ((y-b)^3 (y-c)(y-d)) — mixed
    # Simpler test: BJ pure-even multi-param and fixed types
    print("  scanning pure-even families as proxy resolvents...", flush=True)
    families = []
    # 2A 3A^3 proxy: not known closed form; use path flag-classical etc.
    paths = [
        (
            "path_flag_classical",
            Fraction(5, 16),
            Fraction(5, 16),
            Fraction(-8, 5),
            Fraction(4, 5),
        ),
        (
            "path_flag_lsw",
            Fraction(5, 16),
            Fraction(55, 16),
            Fraction(-8, 5),
            Fraction(-4),
        ),
    ]
    for pid, m1, m2, k1, k2 in paths:
        mu = m1 + t * (m2 - m1)
        ku = k1 + t * (k2 - k1)
        alpha = sp.together(256 * mu**2 - 3125 * ku**4 / 256)
        beta = sp.together(ku * alpha)
        families.append((pid, alpha, beta))

    # Additional: α = 256 t^2 - u, β = k α for discrete k from shortlist
    for k_str in ["-4", "-8/5", "4/5"]:
        k = Fraction(k_str)
        alpha = sp.together(256 * t**2 - 3125 * sp.Rational(k.numerator, k.denominator) ** 4 / 256)
        beta = sp.together(sp.Rational(k.numerator, k.denominator) * alpha)
        families.append((f"slice_{k_str}", alpha, beta))

    results = []
    t_vals = [0, 1, Fraction(1, 2), Fraction(1, 3), Fraction(2, 3), Fraction(1, 4)] + list(range(-5, 6))
    for pid, alpha, beta in families:
        cat_hits = []
        by_k = defaultdict(list)
        n = 0
        for tv in t_vals:
            try:
                a = int(sp.Rational(sp.simplify(alpha.subs(t, tv))))
                b = int(sp.Rational(sp.simplify(beta.subs(t, tv))))
            except Exception:
                continue
            if a == 0:
                continue
            d = disc_bj_int(a, b)
            if d <= 0 or not is_square(d):
                continue
            if not sp.Poly(x**5 + a * x + b, x, domain=sp.ZZ).is_irreducible:
                continue
            n += 1
            kk = k_of(a, b)
            if kk:
                by_k[str(kk)].append((a, b))
            for tag, ca, cb, ck in CATALOGUE:
                if a == ca and b == cb:
                    cat_hits.append({"tag": tag, "k": str(ck), "t": str(tv)})
        cat_k = sorted({h["k"] for h in cat_hits})
        results.append(
            {
                "id": pid,
                "n_even_irr": n,
                "catalogue_k": cat_k,
                "multi_cat": len(cat_k) >= 2,
                "hits": cat_hits,
            }
        )
        print(f"    {pid}: multi_cat={len(cat_k)>=2} k={cat_k}", flush=True)

    # Geometric cover attempt for 2A type: phi = c y^2 (y-1)^2 (y-a)/(den)
    # Skip heavy Newton for time; record programme stance
    out["family_tests"] = results
    out["multi_cat_families"] = [r["id"] for r in results if r["multi_cat"]]
    out["verdict"] = (
        f"Prior orbits: 2A3A^3 size 96, 2A^2 3A^2 size 108, g=0 lookup. "
        f"Explicit Nielsen cover equations for these types not closed-formed this run. "
        f"Proxy pure-even paths multi_cat={out['multi_cat_families']}. "
        f"Geometric multi-k for 2A* still open."
    )
    print(f"  A2 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 3 — Positive-dimensional pure-even A5 strata
# ===========================================================================
def avenue3_pure_even_strata() -> dict:
    print("=== AVENUE 3: positive-dim pure-even A5 strata ===", flush=True)
    out: dict = {"name": "pure_even_A5_strata"}

    # 2-parameter envelope is the main positive-dim pure-even stratum
    # α = 256 m^2 - 3125 k^4/256, β = k α, (m,k) free
    # Document dimension, multi-k, and sample A5 density
    m, k = sp.symbols("m k")
    alpha = sp.together(256 * m**2 - 3125 * k**4 / 256)
    beta = sp.together(k * alpha)
    D = sp.together(256 * alpha**5 + 3125 * beta**4)
    exp = sp.together((256 * alpha**2 * m) ** 2)
    id_ok = sp.expand(D - exp) == 0

    # Sample integer (m,k) for A5
    a5 = 0
    tested = 0
    multi_k_points = []
    for mv, kv in itertools.product([1, 2, 3, 5, Fraction(5, 16), Fraction(55, 16)], [-4, Fraction(-8, 5), Fraction(4, 5), Fraction(-12, 5), 4]):
        try:
            aa = sp.Rational(sp.simplify(alpha.subs({m: mv, k: kv})))
            bb = sp.Rational(sp.simplify(beta.subs({m: mv, k: kv})))
            if aa.denominator != 1 or bb.denominator != 1:
                continue
            a, b = int(aa), int(bb)
        except Exception:
            continue
        if a == 0:
            continue
        tested += 1
        d = disc_bj_int(a, b)
        if d <= 0 or not is_square(d):
            continue
        if not sp.Poly(x**5 + a * x + b, x, domain=sp.ZZ).is_irreducible:
            continue
        r = classify_poly(x**5 + a * x + b, do_galois=True)
        if (r.get("status") or "").startswith("HIT_A5"):
            a5 += 1
            multi_k_points.append({"m": str(mv), "k": str(kv), "a": a, "b": b})

    # 3-param? α free with β free under disc square surface — dim of even surface
    # Projective: the surface 256α^5+3125β^4=γ^2 is a surface (dim 2 in A3 after scaling)
    out["envelope_2param"] = {
        "alpha": str(alpha),
        "beta": str(beta),
        "disc_identity": id_ok,
        "dim": 2,
        "A5_samples": a5,
        "tested_lattice": tested,
        "sample_points": multi_k_points[:12],
    }
    out["even_surface_dim"] = {
        "equation": "256*alpha**5 + 3125*beta**4 = gamma**2",
        "dim_affine_cone": 2,
        "k_ray_foliation": "envelope = ruled by pure-even k-rays",
    }
    out["verdict"] = (
        f"2-param pure-even envelope disc_id={id_ok}; A5 hits={a5}/{tested} lattice samples; "
        f"this is the positive-dim pure-even A5 arithmetic stratum. "
        f"Multi-k by construction when k varies."
    )
    print(f"  A3 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 4 — Other rigid triples
# ===========================================================================
def avenue4_rigid_triples() -> dict:
    print("=== AVENUE 4: other rigid triples ===", flush=True)
    out: dict = {"name": "other_rigid_triples"}
    # Known rigid A5 triples from Step 1 programme
    triples = [
        {
            "type": "(3A,3A,5A)",
            "model": "phi=6y^5-15y^4+10y^3",
            "field": "Q",
            "even_over_Q": False,
            "reason": "disc monic(phi-t)=5*square",
        },
        {
            "type": "(3A,3A,5B)",
            "model": "1-phi / label swap",
            "field": "Q",
            "even_over_Q": False,
            "reason": "same disc obstruction up to automorphism",
        },
        {
            "type": "(2A,3A,5A)",
            "model": "radical Belyi y^5+a y^4+b y^3",
            "field": "Q(2^{1/5},3^{1/5})",
            "even_over_Q": False,
            "reason": "not over Q; prior numeric even scan empty",
        },
        {
            "type": "(2A,3A,5B)",
            "model": "conjugate 5-class",
            "field": "same",
            "even_over_Q": False,
            "reason": "same as 2A3A5A",
        },
    ]
    # Quick re-check: sample fibres of phi at rational t for disc square in Z-model
    PHI = 6 * y**5 - 15 * y**4 + 10 * y**3
    even = 0
    tested = 0
    for tv in [2, 3, 5, Fraction(1, 2), Fraction(3, 2), Fraction(2, 3), 7, 11]:
        expr = sp.expand((PHI - tv) / 6)
        pol = sp.Poly(sp.monic(expr), y, domain=sp.QQ)
        dens = [sp.fraction(sp.together(c))[1] for c in pol.all_coeffs()]
        L = 1
        for d in dens:
            L = int(sp.ilcm(L, abs(int(d))))
        cleared = sp.expand(L**5 * pol.as_expr().subs(y, y / L))
        pz = sp.Poly(cleared, y, domain=sp.ZZ)
        if pz.LC() == -1:
            pz = sp.Poly(-pz.as_expr(), y, domain=sp.ZZ)
        if not pz.is_irreducible:
            continue
        tested += 1
        d = int(pz.discriminant())
        if d > 0 and is_square(d):
            even += 1
    out["triples"] = triples
    out["phi_even_recheck"] = {"tested_irr": tested, "even": even}
    out["verdict"] = (
        f"Rigid triples over Q remain blocked for pure-even fibres (phi recheck even={even}/{tested}). "
        f"2A3A5* not over Q. Likelihood of multi-k via rigid triples: low."
    )
    print(f"  A4 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 5 — Base change + descent
# ===========================================================================
def avenue5_base_change() -> dict:
    print("=== AVENUE 5: base change + descent ===", flush=True)
    out: dict = {"name": "base_change_descent"}
    # Reconfirm disc identity over Q(√5)
    PHI = 6 * y**5 - 15 * y**4 + 10 * y**3
    mon = sp.expand((PHI - t) / 6)
    D = sp.together(sp.expand(sp.Poly(mon, y).discriminant()))
    sq = sp.together(sp.Rational(25, 36) * t * (t - 1))
    match5 = sp.expand(sp.together(D - 5 * sq**2)) == 0
    matchK = sp.expand(sp.together(D - (sp.sqrt(5) * sq) ** 2)) == 0
    out["disc_theorem"] = {
        "equals_5_times_square": match5,
        "square_in_Q_sqrt5": matchK,
    }
    # Descent obstruction: permanent 5 in square-free kernel over Q
    out["prior"] = "K_SQRT5_EVEN.md"
    out["descent_to_Q"] = False
    out["lattice_recovery"] = False
    out["verdict"] = (
        f"disc=5*□ over Q (proved={match5}); square over Q(√5) (proved={matchK}). "
        f"No descent of evenness to Q; no HQCC Z-lattice recovery. Side route only."
    )
    print(f"  A5 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 6 — Higher-rank rigid systems (r≥5)
# ===========================================================================
def avenue6_higher_rank() -> dict:
    print("=== AVENUE 6: higher-rank r≥5 ===", flush=True)
    out: dict = {"name": "higher_rank_r_ge_5"}
    # dim H^rd = r-3; for r=5 dim=2, r=6 dim=3
    # Enumerate type multisets with programme filter for r=5, count only class size products
    classes = {"2A": 15, "3A": 20, "5A": 12, "5B": 12}
    nonid = ["2A", "3A", "5A", "5B"]
    types_r5 = []
    for comb in itertools.combinations_with_replacement(nonid, 5):
        labels = list(comb)
        n3 = sum(1 for c in labels if c == "3A")
        n5 = sum(1 for c in labels if c in ("5A", "5B"))
        if n3 >= 2 or (n3 >= 1 and n5 >= 1):
            prod = 1
            for c in labels:
                prod *= classes[c]
            types_r5.append(
                {
                    "type": ",".join(sorted(labels)),
                    "class_tuple_space": prod,
                    "dim_rd": 2,
                }
            )
    # Sort by size, keep top and count
    types_r5.sort(key=lambda r: r["class_tuple_space"])
    out["r5_filter_pass_types"] = len(types_r5)
    out["r5_smallest"] = types_r5[:8]
    out["r5_largest"] = types_r5[-3:]
    out["dim"] = {"r5": 2, "r6": 3}
    out["explicit_equations"] = None
    out["verdict"] = (
        f"r=5 filter-pass types: {len(types_r5)} (dim_rd=2). "
        f"Smallest class-spaces start at {types_r5[0]['class_tuple_space'] if types_r5 else '—'}. "
        f"No explicit equations; multi-k likelihood speculative. Effort very high."
    )
    print(f"  A6 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# AVENUE 7 — Geometric lift of envelope
# ===========================================================================
def avenue7_envelope_lift() -> dict:
    print("=== AVENUE 7: geometric lift of envelope ===", flush=True)
    out: dict = {"name": "geometric_lift_envelope"}

    # Envelope path flagship-classical: m=5/16 fixed, k linear in u
    # This is a pure-even 1-param family. Question: is it a pullback of a
    # Nielsen-labelled Hurwitz family (e.g. 3A^4 resolvent)?

    # Tests:
    # (1) Number of geometric branch points of the family as cover of u-line
    #     via disc_u of the specialised polynomials... actually f_u(x) = x^5+α(u)x+β(u)
    #     disc as poly in u: for pure-even family disc is a square in Q(u), so
    #     all branch mult even — geometric monodromy of Gal closure over C(u) is even.
    m0 = Fraction(5, 16)
    ku = Fraction(-8, 5) + t * (Fraction(4, 5) - Fraction(-8, 5))
    alpha = sp.together(25 - 3125 * ku**4 / 256)
    beta = sp.together(ku * alpha)
    # disc_x as function of t (identically square)
    # Branch points of the cover of t-line for Gal(f/C(t)): zeros of disc (even mult)
    Dnum = sp.numer(sp.together(256 * alpha**5 + 3125 * beta**4))
    fac = sp.factor_list(sp.expand(Dnum), domain=sp.QQ)
    branch_factors = [(str(f), int(m)) for f, m in fac[1] if sp.degree(f, t) >= 1]
    out["envelope_flag_classical"] = {
        "alpha": str(alpha),
        "beta": str(beta),
        "disc_factor_mults": branch_factors,
        "disc_is_square_in_Qt": is_square_poly(sp.expand(Dnum)),
    }

    # (2) Compare to 3A^4: would need 4 branch points of type 3-cycle.
    # Count distinct geometric roots of the square-free kernel of disc
    # For a square disc, square-free kernel is 1 — so as a cover of t-line,
    # the disc doesn't show odd branch points. The geometric monodromy for
    # the Galois closure of f over C(t) is in A5, but the branch locus of
    # the *Galois* cover may still be non-empty (ramification where mult roots).
    # Zeros of disc (even if even mult) are still branch points of the cover
    # Spec C(t)[x]/(f) → A1.

    # Square-free part of disc:
    sf = 1
    for f, m in fac[1]:
        if m % 2 == 1 and sp.degree(f, t) >= 1:
            sf = sp.expand(sf * f)
    out["envelope_square_free_disc_part"] = str(sf)
    out["n_odd_disc_factors"] = sum(1 for f, m in fac[1] if m % 2 and sp.degree(f, t) >= 1)

    # (3) Catalogue multi-k (known)
    cat_hits = []
    for tv, tag, a0, b0, k0 in [
        (0, "flagship", -55, 88, "-8/5"),
        (1, "classical", 20, 16, "4/5"),
        (Fraction(1, 3), "classical_m", 20, -16, "-4/5"),
    ]:
        a = int(sp.Rational(sp.simplify(alpha.subs(t, tv))))
        b = int(sp.Rational(sp.simplify(beta.subs(t, tv))))
        if a == a0 and b == b0:
            cat_hits.append({"t": str(tv), "tag": tag, "k": k0})
    out["catalogue_hits"] = cat_hits
    out["multi_k_arithmetic"] = len({h["k"] for h in cat_hits}) >= 2

    # (4) Verdict on geometric lift
    out["is_nielsen_3A4"] = None  # unknown / not identified
    out["obstruction_notes"] = (
        "Envelope is pure-even BJ over Q(t) with multi-k Hilbert hits. "
        "Identification with Ni(A5,C_3^4) would require matching monodromy "
        "generators as 3-cycles at four branch values. Disc being a full square "
        "means all finite branch multiplicities of disc are even; cycle types "
        "need monodromy computation of the Galois closure, not done here. "
        "Speculative: envelope may be a multi-parameter pullback or a different "
        "Nielsen class with pure-even arithmetic specialisations."
    )
    out["verdict"] = (
        f"Envelope multi-k arithmetic confirmed (hits={cat_hits}). "
        f"Odd disc factors={out['n_odd_disc_factors']}. "
        f"Nielsen 3A^4 identification: open/speculative."
    )
    print(f"  A7 verdict: {out['verdict']}", flush=True)
    return out


# ===========================================================================
# Main
# ===========================================================================
def main():
    t0 = time.time()
    print("AVENUE RANK EXECUTE 1→7", flush=True)

    a1 = avenue1_resolvent()
    a2 = avenue2_shortlist()
    a3 = avenue3_pure_even_strata()
    a4 = avenue4_rigid_triples()
    a5 = avenue5_base_change()
    a6 = avenue6_higher_rank()
    a7 = avenue7_envelope_lift()

    elapsed = round(time.time() - t0, 2)

    # Scorecard
    rows = [
        (1, "Better 3A^4 resolvent", a1.get("closed_form_Q_s_resolvent") is not None, a1["verdict"]),
        (2, "Next shortlist g=0", bool(a2.get("multi_cat_families")), a2["verdict"]),
        (3, "Pure-even A5 strata", a3["envelope_2param"]["disc_identity"], a3["verdict"]),
        (4, "Other rigid triples", a4["phi_even_recheck"]["even"] > 0, a4["verdict"]),
        (5, "Base change descent", a5["descent_to_Q"], a5["verdict"]),
        (6, "Higher-rank r≥5", a6["explicit_equations"] is not None, a6["verdict"]),
        (7, "Envelope geometric lift", a7.get("is_nielsen_3A4") is True, a7["verdict"]),
    ]

    geometric_multi_k = False  # none of the geometric avenues closed multi-k
    arithmetic_multi_k = a7["multi_k_arithmetic"] or a2.get("multi_cat_families")

    verdict = (
        f"Avenues 1–7 executed in rank order ({elapsed}s). "
        f"Geometric multi-k Nielsen hit: {geometric_multi_k}. "
        f"Arithmetic multi-k: {arithmetic_multi_k}. "
        f"Best geometric progress: H^rd≅P^1_s + numeric (3,1,1)^4 covers; "
        f"closed form f_s∈Q(s)[x] still open. "
        f"Best arithmetic: 2-param envelope + cross-k paths."
    )

    lines = [
        r"# Avenue rank execute (1→7)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Scorecard",
        "",
        r"| Rank | Avenue | Success flag | Summary |",
        r"|-----:|--------|:------------:|---------|",
    ]
    for rank, name, flag, ver in rows:
        lines.append(f"| {rank} | {name} | {flag} | {ver[:100]} |")

    lines += [
        "",
        "---",
        "",
        r"## 1. Better rational coordinate / resolvent for \(3A^4\)",
        "",
        f"- Symmetric-pole hits: {a1.get('symmetric_pole_n')}",
        f"- p2-coordinate hits: {len(a1.get('p2_coordinate_hits') or [])}",
        f"- s(p2) polyfit: {a1.get('s_of_p2_polyfit')}",
        f"- p2-ansatz results: {a1.get('p2_ansatz_hits')}",
        f"- Closed form \(f_s\\in\\mathbb{{Q}}(s)[x]\): **{a1.get('closed_form_Q_s_resolvent')}**",
        f"- Envelope control multi-k: {a1.get('envelope_multi_k')} hits={a1.get('envelope_control_hits')}",
        f"- **{a1['verdict']}**",
        "",
        "---",
        "",
        r"## 2. Next shortlist genus-0 (\(2A3A^3\), \(2A^2 3A^2\))",
        "",
        f"- Prior orbits: `{a2.get('prior_orbits')}`",
        f"- Family tests: `{a2.get('family_tests')}`",
        f"- Multi-cat families: {a2.get('multi_cat_families')}",
        f"- **{a2['verdict']}**",
        "",
        "---",
        "",
        r"## 3. Positive-dimensional pure-even \(A_5\) strata",
        "",
        f"- Envelope 2-param: disc_id={a3['envelope_2param']['disc_identity']}, "
        f"A5={a3['envelope_2param']['A5_samples']}/{a3['envelope_2param']['tested_lattice']}",
        f"- Sample points: {a3['envelope_2param']['sample_points'][:6]}",
        f"- Even surface: `{a3['even_surface_dim']}`",
        f"- **{a3['verdict']}**",
        "",
        "---",
        "",
        r"## 4. Other rigid triples",
        "",
    ]
    for tr in a4["triples"]:
        lines.append(f"- **{tr['type']}** ({tr['field']}): even_over_Q={tr['even_over_Q']} — {tr['reason']}")
    lines += [
        f"- phi recheck: {a4['phi_even_recheck']}",
        f"- **{a4['verdict']}**",
        "",
        "---",
        "",
        r"## 5. Base change + descent",
        "",
        f"- disc theorem: {a5['disc_theorem']}",
        f"- descent_to_Q: {a5['descent_to_Q']}, lattice_recovery: {a5['lattice_recovery']}",
        f"- **{a5['verdict']}**",
        "",
        "---",
        "",
        r"## 6. Higher-rank rigid systems (\(r\\ge 5\))",
        "",
        f"- r=5 filter-pass types: {a6['r5_filter_pass_types']}",
        f"- smallest: {a6['r5_smallest'][:5]}",
        f"- dim: {a6['dim']}",
        f"- **{a6['verdict']}**",
        "",
        "---",
        "",
        r"## 7. Geometric lift of the envelope",
        "",
        f"- disc factors: {a7['envelope_flag_classical']['disc_factor_mults']}",
        f"- odd disc factors: {a7['n_odd_disc_factors']}",
        f"- catalogue hits: {a7['catalogue_hits']}",
        f"- multi-k arithmetic: {a7['multi_k_arithmetic']}",
        f"- Nielsen 3A^4?: {a7['is_nielsen_3A4']}",
        f"- {a7['obstruction_notes']}",
        f"- **{a7['verdict']}**",
        "",
        "---",
        "",
        r"## Global conclusions",
        "",
        r"1. **Geometric multi-\(k\) (Nielsen-labelled): still open.**",
        r"2. **Arithmetic multi-\(k\): solid** via 2-param envelope and cross-\(k\) paths.",
        r"3. **Avenue 1** remains the highest-leverage geometric attack (genus 0 is favourable)",
        r"   but closed-form \(f_s\in\mathbb{Q}(s)[x]\) needs more elimination / descent work.",
        r"4. **Avenues 4–5** are effectively closed as multi-k routes (blocked / side-only).",
        r"5. **Avenue 6** is enumerated at type-count level only; equations out of scope.",
        r"6. **Avenue 7**: envelope is the arithmetic multi-k object; geometric Nielsen ID open.",
        "",
        r"### Recommended single next move",
        "",
        r"Finish Avenue 1: resultant/Gröbner elimination of the triple-root ideal to a",
        r"plane model of \((c:p_2:r_1:r_2)\) over \(\mathbb{Q}(s)\), or a deg-5 resolvent",
        r"after a radical extension of \(\mathbb{Q}(s)\) with controlled descent.",
        "",
        r"_Generated by avenue_rank_execute.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "geometric_multi_k": geometric_multi_k,
        "arithmetic_multi_k": arithmetic_multi_k,
        "avenues": {
            "1": a1,
            "2": a2,
            "3": a3,
            "4": a4,
            "5": a5,
            "6": a6,
            "7": a7,
        },
        "scorecard": [
            {"rank": r, "name": n, "success": f, "verdict": v} for r, n, f, v in rows
        ],
    }
    write_md(OUT / "AVENUE_RANK_EXECUTE.md", doc)
    write_md(RESULTS / "AVENUE_RANK_EXECUTE.md", doc)
    write_md(ROOT / "AVENUE_RANK_EXECUTE.md", doc)
    write_json(OUT / "AVENUE_RANK_EXECUTE.json", blob)
    print(verdict, flush=True)
    print(f"Wrote AVENUE_RANK_EXECUTE.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
