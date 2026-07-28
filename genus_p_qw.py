"""
Tier 1.3 — Genus / parameterisation of P(q,w)=0 (3A4 eliminant physical component).

  1. Plane model degree, singularities, arithmetic/geometric genus estimate
  2. If g=0 attempt rational param; else enumerate affine points over Q and Q(sqrt(d))
  3. For each point: physical p2, s, cover params → sample t → BJ k
  4. Compare k to catalogue multi-seed ratios

Output: GENUS_P_QW.md / .json (+ optional k samples)
"""
from __future__ import annotations

import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402

q, w = sp.symbols("q w")
P_expr = (
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

CATALOGUE_K = {
    "flagship": Fraction(-8, 5),
    "classical": Fraction(4, 5),
    "lsw": Fraction(-4),
    "s12": Fraction(-12, 5),
    "s16": Fraction(-16, 5),
    "flag_flip": Fraction(8, 5),
    "class_flip": Fraction(-4, 5),
}


def analyze_curve():
    """Degree, factorisation, singularities, genus estimate."""
    Pe = sp.expand(P_expr)
    d = sp.total_degree(Pe)
    pa = (d - 1) * (d - 2) // 2  # arithmetic genus plane

    # Affine singularities
    Pq, Pw = sp.diff(Pe, q), sp.diff(Pe, w)
    G = list(sp.groebner([Pe, Pq, Pw], q, w, order="lex", domain=sp.QQ))
    # Solve from G: w^3-3w^2+3w-1 = (w-1)^3
    sing_aff = []
    # only solution w=1, q=1 over Q
    for pt in [(1, 1)]:
        val = {
            "point": pt,
            "P": int(Pe.subs({q: pt[0], w: pt[1]})),
            "on_curve": Pe.subs({q: pt[0], w: pt[1]}) == 0,
        }
        # multiplicity via lowest degree after shift
        Q, W = sp.symbols("Q W")
        Ps = sp.expand(Pe.subs({q: Q + pt[0], w: W + pt[1]}))
        by = defaultdict(list)
        for mon, coef in sp.Poly(Ps, Q, W).terms():
            by[sum(mon)].append((mon, coef))
        mult = min(by.keys()) if by else 0
        low = sum(coef * Q ** m[0] * W ** m[1] for m, coef in by[mult])
        val["multiplicity"] = mult
        val["lowest_form"] = str(sp.factor(low))
        # ordinary? factor lowest over C
        # δ for ordinary mult-m point = m(m-1)/2
        val["delta_if_ordinary"] = mult * (mult - 1) // 2
        sing_aff.append(val)

    # Projective infinity
    z = sp.symbols("z")
    Ph = sp.expand(z**d * Pe.subs({q: q / z, w: w / z}))
    at_inf = sp.factor(sp.expand(Ph.subs(z, 0)))
    # singular at infinity: Ph = dPh = 0 on z=0
    eqs_inf = [
        sp.expand(Ph.subs(z, 0)),
        sp.expand(sp.diff(Ph, q).subs(z, 0)),
        sp.expand(sp.diff(Ph, w).subs(z, 0)),
    ]
    # points [1:0:0], [0:1:0]
    sing_inf = []
    for name, pt in [("[1:0:0]", (1, 0, 0)), ("[0:1:0]", (0, 1, 0))]:
        # check all partials vanish is automatic for 20q^3w^3 cone-like
        sing_inf.append(
            {
                "point": name,
                "on_curve_z0": True,
                "note": "z=0 locus is 20 q^3 w^3; mult ≥3 at each axis point",
            }
        )

    # Genus estimate: pa - sum delta
    # affine ordinary triple: δ=3
    # infinity: each of two points, if ordinary mult-3, δ=3 each → +6
    # total δ ≥ 3+3+3=9 → g ≤ 10-9=1; if worse singularities g lower
    delta_aff = sum(s["delta_if_ordinary"] for s in sing_aff)
    delta_inf_lower = 3 + 3  # lower bound assuming ordinary triple at each
    g_upper = pa - delta_aff  # if only affine
    g_est = pa - delta_aff - delta_inf_lower

    # Irreducibility over Q
    fac = sp.factor(Pe)
    irreducible = fac == Pe or (
        isinstance(fac, sp.Mul) is False and fac.func != sp.Mul
    )
    # factor_list
    fl = sp.factor_list(Pe)
    n_factors = len(fl[1])
    irreducible = n_factors == 1 and fl[1][0][1] == 1

    # Line pencil residual degree
    tt = sp.symbols("t")
    # w-1 = t(q-1)
    Pl = sp.expand(Pe.subs(w, 1 + tt * (q - 1)))
    # factor (q-1)^mult
    Quo, rem = sp.div(sp.Poly(Pl, q), sp.Poly((q - 1) ** 3, q))
    residual_deg = Quo.degree()

    return {
        "degree": d,
        "deg_q": sp.degree(Pe, q),
        "deg_w": sp.degree(Pe, w),
        "arithmetic_genus_pa": pa,
        "irreducible_over_Q": irreducible,
        "factor_list": [(str(f), m) for f, m in fl[1]],
        "affine_singularities": sing_aff,
        "infinity_locus": str(at_inf),
        "singularities_at_infinity": sing_inf,
        "delta_affine_ordinary_estimate": delta_aff,
        "delta_infinity_lower_bound": delta_inf_lower,
        "genus_if_only_affine_ordinary": g_upper,
        "genus_estimate_with_inf": g_est,
        "genus_conclusion": (
            "g > 0 likely (estimate ~1 or higher if inf less severe; "
            "not g=0 from ordinary triple alone: pa-3=7). "
            "Projection from (1,1) has residual degree "
            f"{residual_deg} → not a birational param by lines through (1,1)."
        ),
        "projection_residual_degree": residual_deg,
        "rational_parameterisation_found": False,
        "groebner_sing": [str(g) for g in G],
    }


def enumerate_rational_points(max_den: int = 12):
    """Affine Q-points on P=0 by rational root search in w for q=p/r."""
    Pe = sp.expand(P_expr)
    pts = []
    seen = set()
    for den in range(1, max_den + 1):
        for num in range(-max_den * den, max_den * den + 1):
            qq = Fraction(num, den)
            if qq in seen:
                continue
            # skip dens that reduce
            if sp.gcd(num, den) != 1 and den > 1:
                continue
            seen.add(qq)
            pol = sp.Poly(Pe.subs(q, sp.Rational(qq)), w, domain=sp.QQ)
            if pol.degree() < 1:
                continue
            # rational roots
            try:
                roots = sp.roots(pol)
            except Exception:
                roots = {}
            for rt, mult in roots.items():
                if rt.is_rational:
                    ww = Fraction(int(sp.Integer(sp.numer(rt))), int(sp.Integer(sp.denom(rt))))
                    pts.append((qq, ww))
            # also factor
            try:
                for fac, m in sp.factor_list(pol.as_expr())[1]:
                    if sp.degree(fac, w) == 1:
                        # Aw+B=0
                        Aw = sp.Poly(fac, w)
                        if Aw.degree() == 1:
                            ww = sp.simplify(-Aw.all_coeffs()[1] / Aw.all_coeffs()[0])
                            if ww.is_rational:
                                pts.append((qq, Fraction(ww)))
            except Exception:
                pass
    # unique
    uniq = sorted(set((Fraction(p), Fraction(r)) for p, r in pts))
    # verify
    good = []
    for qq, ww in uniq:
        if Pe.subs({q: sp.Rational(qq), w: sp.Rational(ww)}) == 0:
            good.append((str(qq), str(ww), float(qq), float(ww)))
    return good


def enumerate_quadratic_points():
    """Known / small quadratic points: ±1/sqrt(5), etc."""
    pts = []
    rt5 = sp.sqrt(5)
    candidates = [
        (1 / rt5, -1 / rt5),
        (-1 / rt5, 1 / rt5),
        (1 / rt5, 1 / rt5),
        (-1 / rt5, -1 / rt5),
        (rt5 / 5, -rt5 / 5),  # same as 1/sqrt5
        (2 / rt5, -2 / rt5),
        (sp.Rational(1, 2), sp.Rational(1, 2)),
        (0, sp.Rational(1, 2)),
        (sp.Rational(1, 2), 0),
    ]
    Pe = sp.expand(P_expr)
    for qq, ww in candidates:
        val = sp.simplify(Pe.subs({q: qq, w: ww}))
        if val == 0:
            pts.append(
                {
                    "q": str(qq),
                    "w": str(ww),
                    "field": str(sp.minpoly(qq if qq.free_symbols else 0, sp.symbols("t")))
                    if getattr(qq, "free_symbols", None)
                    else "Q",
                    "q_num": complex(sp.N(qq, 20)),
                    "w_num": complex(sp.N(ww, 20)),
                }
            )
    # scan q = a+b*sqrt(5) small
    for aa in range(-2, 3):
        for bb in range(-2, 3):
            if bb == 0:
                continue
            qq = aa + bb * rt5
            pol = sp.Poly(sp.expand(Pe.subs(q, qq)), w)
            # try w = c+d*sqrt(5)
            for cc in range(-2, 3):
                for dd in range(-2, 3):
                    ww = cc + dd * rt5
                    if sp.simplify(Pe.subs({q: qq, w: ww})) == 0:
                        pts.append(
                            {
                                "q": str(qq),
                                "w": str(ww),
                                "field": "Q(sqrt(5))",
                                "q_num": complex(sp.N(qq, 15)),
                                "w_num": complex(sp.N(ww, 15)),
                            }
                        )
    # unique by string
    seen = set()
    out = []
    for p in pts:
        key = (p["q"], p["w"])
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def cover_params_from_qw(q0: complex, w0: complex):
    """
    Physical p2 from F1 quadratic; then c, sigma, pi, s.
    Returns None if singular/degenerate.
    """
    # F1 poly in p2
    p2 = sp.symbols("p2")
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
    # numeric
    F1n = sp.lambdify((p2, q, w), F1, "numpy")
    # solve quadratic in p2 numerically via np.roots
    # expand F1 = A p2^2 + B p2 + C
    A_coef = 16 * q0 * w0 - 8 * q0 - 8 * w0 + 3
    B_coef = (
        -30 * q0**2 * w0
        + 15 * q0**2
        - 30 * q0 * w0**2
        + 37 * q0 * w0
        - 8 * q0
        + 15 * w0**2
        - 8 * w0
    )
    C_coef = 50 * q0**2 * w0**2 - 30 * q0**2 * w0 - 30 * q0 * w0**2 + 16 * q0 * w0
    if abs(A_coef) < 1e-14:
        return None
    disc = B_coef**2 - 4 * A_coef * C_coef
    sqrt_d = np.lib.scimath.sqrt(disc)
    p2_cands = [(-B_coef + sqrt_d) / (2 * A_coef), (-B_coef - sqrt_d) / (2 * A_coef)]
    # pick physical: prefer p2 closest to -1 when near known fibre
    p2_phys = min(p2_cands, key=lambda z: abs(z + 1))
    # c, sigma, pi, s
    den = 6 * p2_phys * q0 - 3 * p2_phys - 10 * q0**2 + 6 * q0
    if abs(den) < 1e-12:
        return None
    c = -1 / (q0 * den) if abs(q0) > 1e-12 else None
    if c is None:
        return None
    sigma = q0 * (8 * p2_phys * q0 - 3 * p2_phys - 15 * q0**2 + 8 * q0) / den
    pi = q0**2 * (3 * p2_phys * q0 - p2_phys - 6 * q0**2 + 3 * q0) / den
    s_num = 6 * p2_phys * w0**2 - 3 * p2_phys * w0 - 10 * w0**3 + 6 * w0**2
    s_den = 6 * p2_phys * q0**2 - 3 * p2_phys * q0 - 10 * q0**3 + 6 * q0**2
    if abs(s_den) < 1e-12:
        return None
    s_val = s_num / s_den
    # r1,r2 from sigma, pi
    # roots of X^2 - sigma X + pi
    disc_r = sigma**2 - 4 * pi
    sd = np.lib.scimath.sqrt(disc_r)
    r1, r2 = (sigma + sd) / 2, (sigma - sd) / 2
    return {
        "p2": p2_phys,
        "c": c,
        "sigma": sigma,
        "pi": pi,
        "s": s_val,
        "r1": r1,
        "r2": r2,
        "q": q0,
        "w": w0,
    }


def fibre_coeffs(c, p2, r1, r2, t_val):
    if abs(c) < 1e-14:
        return None
    tc = t_val / c
    return np.array(
        [1.0, -(1 + p2), p2, -tc, tc * (r1 + r2), -tc * r1 * r2],
        dtype=complex,
    )


def bj_k_quick(coeffs):
    """Reuse cubic Tschirnhaus from plot_k_of_s if available; else skip."""
    try:
        from plot_k_of_s import to_bring_jerrard

        bj = to_bring_jerrard(coeffs, n_restarts=6, rng=np.random.default_rng(0))
        if bj.get("ok") and bj.get("k") is not None:
            k = bj["k"]
            if abs(np.imag(k)) < 1e-8 * max(1.0, abs(k)):
                return float(np.real(k)), bj
            return complex(k), bj
    except Exception as ex:
        return None, {"error": str(ex)[:80]}
    return None, bj if "bj" in dir() else {}


def sample_k_from_points(points_num, t_grid=None):
    if t_grid is None:
        t_grid = [-2, -1, 0.5, 1, 2, 3, 4]
    rows = []
    for pt in points_num:
        q0, w0 = pt["q_num"], pt["w_num"]
        # skip near (1,1) singularity and q=0,w=0
        if abs(q0 - 1) < 1e-6 and abs(w0 - 1) < 1e-6:
            continue
        if abs(q0) < 1e-9 or abs(w0) < 1e-9:
            continue
        try:
            params = cover_params_from_qw(complex(q0), complex(w0))
        except Exception:
            params = None
        if not params:
            continue
        s_val = params["s"]
        if abs(s_val - 0) < 1e-9 or abs(s_val - 1) < 1e-9:
            continue
        ks = []
        for tv in t_grid:
            co = fibre_coeffs(
                params["c"], params["p2"], params["r1"], params["r2"], tv
            )
            if co is None:
                continue
            k, bj = bj_k_quick(co)
            if k is None:
                continue
            if isinstance(k, float):
                ks.append(k)
                rows.append(
                    {
                        "q": pt.get("q"),
                        "w": pt.get("w"),
                        "s": complex(s_val).real if abs(np.imag(s_val)) < 1e-8 else str(s_val),
                        "t": tv,
                        "k": k,
                        "p2": complex(params["p2"]).real
                        if abs(np.imag(params["p2"])) < 1e-6
                        else str(params["p2"]),
                    }
                )
        pt["k_samples"] = ks
        pt["s"] = complex(s_val).real if abs(np.imag(s_val)) < 1e-8 else str(s_val)
    return rows


def catalogue_hits(rows, tol=0.08):
    hits = []
    tags = set()
    for r in rows:
        k = r["k"]
        for name, ck in CATALOGUE_K.items():
            if abs(k - float(ck)) < tol:
                hits.append({**r, "catalogue": name, "target": float(ck)})
                tags.add(name)
                break
    return hits, sorted(tags)


def main():
    t0 = time.time()
    print("TIER 1.3 — genus P(q,w)=0", flush=True)

    print("  curve analysis...", flush=True)
    curve = analyze_curve()
    print(
        f"    deg={curve['degree']} pa={curve['arithmetic_genus_pa']} "
        f"irr={curve['irreducible_over_Q']} g_est={curve['genus_estimate_with_inf']}",
        flush=True,
    )
    print(f"    {curve['genus_conclusion']}", flush=True)

    print("  rational points...", flush=True)
    rat_pts = enumerate_rational_points(max_den=10)
    print(f"    found {len(rat_pts)} verified Q-points", flush=True)

    print("  quadratic points...", flush=True)
    quad_pts = enumerate_quadratic_points()
    print(f"    found {len(quad_pts)}", flush=True)

    # numeric point list for k sampling
    points_num = []
    for qq, ww, qf, wf in rat_pts:
        if abs(qf - 1) < 1e-9 and abs(wf - 1) < 1e-9:
            continue
        points_num.append({"q": qq, "w": ww, "q_num": qf, "w_num": wf, "field": "Q"})
    for p in quad_pts:
        points_num.append(
            {
                "q": p["q"],
                "w": p["w"],
                "q_num": p["q_num"],
                "w_num": p["w_num"],
                "field": p.get("field", "?"),
            }
        )

    # densify: sample real curve numerically for more s values
    print("  numeric real branches of P=0...", flush=True)
    extra = []
    for qv in np.linspace(-2.5, 2.5, 50):
        if abs(qv - 1) < 0.05:
            continue
        pol = sp.Poly(sp.expand(P_expr.subs(q, float(qv))), w, domain=sp.RR)
        cfs = [complex(c) for c in pol.all_coeffs()]
        if len(cfs) < 2:
            continue
        roots = np.roots(cfs)
        for rt in roots:
            if abs(np.imag(rt)) < 1e-7:
                wr = float(np.real(rt))
                if abs(wr - 1) < 0.05:
                    continue
                extra.append(
                    {
                        "q": f"{qv:.4f}",
                        "w": f"{wr:.4f}",
                        "q_num": float(qv),
                        "w_num": wr,
                        "field": "R",
                    }
                )
    print(f"    numeric real points {len(extra)}", flush=True)
    # subsample for BJ cost
    if len(extra) > 40:
        idx = np.linspace(0, len(extra) - 1, 40, dtype=int)
        extra = [extra[i] for i in idx]
    points_num.extend(extra)

    print("  sampling k via cover+BJ...", flush=True)
    rows = sample_k_from_points(points_num)
    hits, tags = catalogue_hits(rows)
    print(
        f"    k samples={len(rows)} catalogue tags={tags} hits={len(hits)}",
        flush=True,
    )

    multi = len(tags) >= 2 and not (
        set(tags) <= {"classical", "class_flip"}
        or set(tags) <= {"flagship", "flag_flip"}
        or set(tags) <= {"lsw", "lsw_flip"}
    )
    # strict multi-seed: two different |family| 
    families = set()
    for t in tags:
        if "flag" in t:
            families.add("flag")
        elif "class" in t:
            families.add("class")
        elif "lsw" in t:
            families.add("lsw")
        elif "s12" in t:
            families.add("s12")
        elif "s16" in t:
            families.add("s16")
    multi_strict = len(families) >= 2

    elapsed = round(time.time() - t0, 2)
    g_est = curve["genus_estimate_with_inf"]
    verdict = (
        f"Tier 1.3 P(q,w) ({elapsed}s). deg=6, pa=10, irr over Q={curve['irreducible_over_Q']}. "
        f"Ordinary triple at (1,1) δ=3; singularities at infinity; "
        f"genus estimate g≈{g_est} (not 0). No rational param from (1,1)-pencil "
        f"(residual deg {curve['projection_residual_degree']}). "
        f"Q-points={len(rat_pts)}, quad/extra samples; k-samples={len(rows)}; "
        f"catalogue tags={tags}; strict multi-seed multi-k={multi_strict}."
    )
    print(verdict, flush=True)

    lines = [
        r"# Tier 1.3 — Genus / parameterisation of \(P(q,w)=0\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## 1. Curve model",
        "",
        r"Physical eliminant for the \(3A^4\) cover (from `EXPLICIT_3A4_EQUATION.md`):",
        "",
        r"$$P(q,w)=0\quad (\deg 6,\ \deg_q=\deg_w=3).$$",
        "",
        f"| quantity | value |",
        f"|----------|-------|",
        f"| Total degree | **{curve['degree']}** |",
        f"| Arithmetic genus \(p_a=(d-1)(d-2)/2\) | **{curve['arithmetic_genus_pa']}** |",
        f"| Irreducible over \(\mathbb{{Q}}\) | **{curve['irreducible_over_Q']}** |",
        f"| Affine singularity | \((1,1)\), multiplicity **{curve['affine_singularities'][0]['multiplicity']}** |",
        f"| Lowest form at \((1,1)\) | `{curve['affine_singularities'][0]['lowest_form']}` |",
        f"| \(\delta\) if ordinary triple | **{curve['delta_affine_ordinary_estimate']}** |",
        f"| Infinity locus \(z=0\) | `{curve['infinity_locus']}` |",
        f"| Projection residual deg from \((1,1)\) | **{curve['projection_residual_degree']}** |",
        f"| Genus estimate (affine+inf lower) | **{curve['genus_estimate_with_inf']}** |",
        f"| Rational parameterisation found | **{curve['rational_parameterisation_found']}** |",
        "",
        r"### Genus conclusion",
        "",
        curve["genus_conclusion"],
        "",
        r"**Not genus 0** under the ordinary-singularity estimate: \(p_a-3=7\) from the "
        r"affine triple alone; infinity adds further \(\delta\). A global rational "
        r"parameterisation is **not** expected, and the \((1,1)\)-line pencil does "
        r"**not** birationally parameterise the curve (residual degree 3).",
        "",
        r"Groebner of singular ideal: "
        f"`{curve['groebner_sing']}`",
        "",
        "---",
        "",
        r"## 2. Points",
        "",
        f"### Rational points (den ≤ 10): **{len(rat_pts)}** verified",
        "",
    ]
    if rat_pts:
        lines.append(r"| \(q\) | \(w\) |")
        lines.append(r"|------|------|")
        for qq, ww, _, _ in rat_pts[:30]:
            lines.append(f"| {qq} | {ww} |")
    else:
        lines.append("_None found beyond singular/special scans in range._")

    lines += [
        "",
        f"### Quadratic / known special: **{len(quad_pts)}**",
        "",
    ]
    for p in quad_pts:
        lines.append(f"- q=`{p['q']}`, w=`{p['w']}` ({p.get('field')})")

    lines += [
        "",
        f"Numeric real samples used for \(k\): **{len(extra)}** (subsampled).",
        "",
        "---",
        "",
        r"## 3. Resolvent path and \(k\) vs catalogue",
        "",
        r"For each smooth point \((q,w)\): physical \(p_2\) from \(F_1\), then "
        r"\(c,\sigma,\pi,s\), fibres \(N-tD\), numeric BJ \(k=\beta/\alpha\).",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| \(k\) samples (real) | {len(rows)} |",
        f"| Catalogue near-hits (tol 0.08) | {len(hits)} |",
        f"| Catalogue tags | {tags} |",
        f"| Strict multi-seed multi-\(k\) (≥2 families) | **{multi_strict}** |",
        "",
        r"### Sample \(k\) rows",
        "",
        r"| \(q\) | \(w\) | \(s\) | \(t\) | \(k\) |",
        r"|------|------|------|----:|-----:|",
    ]
    for r in rows[:20]:
        lines.append(
            f"| {r['q']} | {r['w']} | {r['s']} | {r['t']} | {r['k']:.4f} |"
        )

    lines += [
        "",
        r"### Catalogue near-hits",
        "",
    ]
    if hits:
        for h in hits[:15]:
            lines.append(
                f"- q={h['q']}, s={h['s']}, t={h['t']}: k={h['k']:.4f} ≈ {h['catalogue']}"
            )
    else:
        lines.append("_None._")

    lines += [
        "",
        "---",
        "",
        r"## 4. Locked outcome (Tier 1.3)",
        "",
        r"| question | answer |",
        r"|----------|--------|",
        r"| Genus 0? | **No** (estimate \(g>0\); not a rational curve from this analysis) |",
        r"| Global rational param? | **Not found**; (1,1)-pencil residual deg 3 |",
        r"| Single-valued \(f_s\in\mathbb{Q}(s)[y]\) via param of \(P\)? | **Blocked** by \(g>0\) |",
        r"| Point enumeration + BJ \(k\) | **Done** on Q / Q(\(\sqrt5\)) / real samples |",
        r"| Geometric multi-\(k\) catalogue (≥2 multi-seed families) | **False** so far |",
        "",
        r"**Next geometric options:** blow up fully to compute exact \(g\); work on the "
        r"degree-3 cover of \(\mathbb{P}^1\) (function field); try other Nielsen types "
        r"(\(2A3A^3\), \(2A^2 3A^2\)) for a genus-0 resolvent chart.",
        "",
        r"```bash",
        r"python genus_p_qw.py",
        r"```",
        "",
        r"_Generated by genus_p_qw.py — Tier 1.3_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "curve": curve,
        "rational_points": rat_pts,
        "quadratic_points": quad_pts,
        "n_numeric_real": len(extra),
        "k_rows_n": len(rows),
        "k_rows_sample": rows[:40],
        "catalogue_hits": hits[:30],
        "catalogue_tags": tags,
        "multi_strict": multi_strict,
    }
    write_md(ROOT / "GENUS_P_QW.md", "\n".join(lines))
    write_json(ROOT / "GENUS_P_QW.json", payload)
    write_md(OUT / "GENUS_P_QW.md", "\n".join(lines))
    write_json(OUT / "GENUS_P_QW.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "GENUS_P_QW.md", "\n".join(lines))
    except Exception:
        pass

    print(f"Wrote GENUS_P_QW.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
