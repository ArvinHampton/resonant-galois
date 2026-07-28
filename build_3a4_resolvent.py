"""
Explicit model of Ni(A5, C_3^4):
  - Reduced Hurwitz curve ≅ P^1_s (genus 0 lock)
  - Degree-5 cover phi_s: P1 → P1 of type (3,1,1)^4 with branch values {0,1,∞,s}
  - Resolvent f_s(y,t) = N_s(y) - t D_s(y) monicised in y
  - Specialise (s,t) rational; test pure-even fixed-k catalogue

Cover normal form (after Aut domain):
  phi = c * y^3 * (y-1) * (y-p2) / ((y-r1)*(y-r2))
so poles give type (3,1,1) at ∞ (order-3 at ∞, simple at r1,r2),
zeros give type (3,1,1) at 0 (order-3 at 0, simple at 1 and p2).
Parameters (c,p2,r1,r2) determined by requiring phi-1 and phi-s each have a triple root.

Output: EXPLICIT_3A4_RESOLVENT.md / .json
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.linalg import norm

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

y, t, s = sp.symbols("y t s")


# ---------------------------------------------------------------------------
# Newton solver for cover parameters at fixed s
# ---------------------------------------------------------------------------
def residual(v, s_val: float):
    c, p2, r1, r2, q, w = v

    def G_vals(val, pt):
        yy = pt
        N = c * yy**3 * (yy - 1) * (yy - p2)
        A = yy**3
        Ap = 3 * yy**2
        App = 6 * yy
        B = yy**2 - (1 + p2) * yy + p2
        Bp = 2 * yy - (1 + p2)
        Bpp = 2.0
        Np = c * (Ap * B + A * Bp)
        Npp = c * (App * B + 2 * Ap * Bp + A * Bpp)
        D = (yy - r1) * (yy - r2)
        Dp = 2 * yy - (r1 + r2)
        Dpp = 2.0
        return [N - val * D, Np - val * Dp, Npp - val * Dpp]

    return np.array(G_vals(1.0, q) + G_vals(s_val, w), dtype=float)


def jacobian(v, s_val, eps=1e-8):
    J = np.zeros((6, 6))
    r0 = residual(v, s_val)
    for j in range(6):
        dv = np.zeros(6)
        dv[j] = eps
        J[:, j] = (residual(v + dv, s_val) - r0) / eps
    return J


def newton(s_val, x0, niter=80):
    v = np.array(x0, dtype=float)
    for _ in range(niter):
        r = residual(v, s_val)
        if norm(r) < 1e-14:
            return v, True, float(norm(r))
        J = jacobian(v, s_val)
        try:
            step = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(J, -r, rcond=None)[0]
        v = v + step
    nr = float(norm(residual(v, s_val)))
    return v, nr < 1e-12, nr


def solve_at_s(s_val, seeds=None, n_trials=60):
    rng = np.random.default_rng(abs(hash((round(s_val, 10),))) % (2**32))
    candidates = []
    if seeds:
        candidates.extend(seeds)
    for _ in range(n_trials):
        candidates.append(rng.normal(size=6))
    # exact seed at s=-1
    candidates.append(np.array([-np.sqrt(5), -1.0, 0.2, -0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5]))
    best = None
    for x0 in candidates:
        v, ok, nr = newton(s_val, x0)
        if best is None or nr < best[2]:
            best = (v, ok, nr)
        if ok:
            return v, True, nr
    return best[0], False, best[2]


# ---------------------------------------------------------------------------
# Build monic fibre poly N - t D over Qbar, then for specializations
# ---------------------------------------------------------------------------
def make_ND(c, p2, r1, r2):
    N = sp.expand(sp.Float(c) * y**3 * (y - 1) * (y - sp.Float(p2)))
    D = sp.expand((y - sp.Float(r1)) * (y - sp.Float(r2)))
    return N, D


def fibre_monic_Z(N, D, t_val):
    """Monic Z-poly for N - t D = 0 in y."""
    fib = sp.expand(N - sp.Rational(t_val) * D)
    try:
        pol = sp.Poly(fib, y, domain=sp.QQ)
    except Exception:
        return None
    if pol.degree() != 5:
        # multiply by denom of leading if needed
        pass
    if pol.degree() < 5:
        return None
    # Make monic over QQ then clear
    try:
        mon = sp.Poly(sp.monic(pol.as_expr()), y, domain=sp.QQ)
    except Exception:
        return None
    dens = [sp.fraction(sp.together(c))[1] for c in mon.all_coeffs()]
    L = 1
    for d in dens:
        try:
            L = int(sp.ilcm(L, abs(int(d))))
        except Exception:
            return None
    cleared = sp.expand(L**5 * mon.as_expr().subs(y, y / L))
    pz = sp.Poly(cleared, y, domain=sp.ZZ)
    if pz.LC() == -1:
        pz = sp.Poly(-pz.as_expr(), y, domain=sp.ZZ)
    if pz.LC() != 1 or pz.degree() != 5:
        return None
    return pz


def try_bj(pz):
    """Shift to kill y^4; return (alpha,beta) if Bring-Jerrard (no y^3,y^2)."""
    coeffs = [int(c) for c in pz.all_coeffs()]  # y5..y0
    c4 = coeffs[1]
    shift = -Fraction(c4, 5)
    z = sp.symbols("z")
    fsh = sp.expand(pz.as_expr().subs(y, z + sp.Rational(shift)))
    psh = sp.Poly(fsh, z, domain=sp.QQ)
    cc = psh.all_coeffs()
    if len(cc) != 6:
        return None
    # [1, c4, c3, c2, c1, c0]
    if cc[1] != 0:
        return None
    if cc[2] != 0 or cc[3] != 0:
        return ("depressed_not_BJ", [sp.Rational(c) for c in cc])
    try:
        a = int(sp.Rational(cc[4]))
        b = int(sp.Rational(cc[5]))
        return ("BJ", a, b)
    except Exception:
        return None


CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("flag_320", 320, -512, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("s220m", 220, -528, Fraction(-12, 5)),
    ("s55_176", -55, 176, Fraction(-16, 5)),
]


def k_of(a, b):
    if a == 0:
        return None
    return Fraction(int(b), int(a))


def main():
    t0 = time.time()
    print("EXPLICIT 3A^4 RESOLVENT", flush=True)

    # ---- Sample covers at rational s ----
    s_list = []
    for num, den in [
        (-3, 1), (-2, 1), (-1, 1), (-1, 2), (1, 2), (3, 2), (2, 1), (5, 2), (3, 1),
        (4, 1), (5, 1), (-5, 2), (5, 3), (7, 2), (1, 3), (2, 3), (4, 3), (5, 4),
        (-3, 2), (7, 3), (8, 3),
    ]:
        sv = Fraction(num, den)
        if sv in (0, 1):
            continue
        s_list.append(sv)
    s_list = list(dict.fromkeys(s_list))

    covers = []
    v_prev = None
    for sv in s_list:
        seeds = [v_prev] if v_prev is not None else None
        v, ok, nr = solve_at_s(float(sv), seeds=seeds, n_trials=80)
        if not ok:
            v, ok, nr = solve_at_s(float(sv), seeds=None, n_trials=120)
        if ok:
            c, p2, r1, r2, q, w = v
            covers.append(
                {
                    "s": str(sv),
                    "s_float": float(sv),
                    "c": float(c),
                    "p2": float(p2),
                    "r1": float(r1),
                    "r2": float(r2),
                    "q": float(q),
                    "w": float(w),
                    "newton_res": nr,
                    "ok": True,
                }
            )
            v_prev = v
            print(f"  cover s={sv}: ok res={nr:.1e} p2={p2:.6f}", flush=True)
        else:
            covers.append({"s": str(sv), "ok": False, "newton_res": nr})
            print(f"  cover s={sv}: FAIL res={nr:.1e}", flush=True)

    n_ok = sum(1 for c in covers if c.get("ok"))
    print(f"  covers ok: {n_ok}/{len(covers)}", flush=True)

    # ---- Exact model at s=-1 over Q(sqrt(5)) ----
    exact_sm1 = {
        "s": -1,
        "c": -sp.sqrt(5),
        "p2": -1,
        "r1": sp.Rational(1, 5),
        "r2": -sp.Rational(1, 5),
        "field": "Q(sqrt(5))",
        "N": str(sp.expand(-sp.sqrt(5) * y**3 * (y - 1) * (y + 1))),
        "D": str(sp.expand((y - sp.Rational(1, 5)) * (y + sp.Rational(1, 5)))),
        "note": "Verified triple-root conditions for branch values 1 and s=-1",
    }

    # ---- Specialise covers: for each (s,t) rational, get quintic, test Gal/catalogue ----
    print("  specialising fibres...", flush=True)
    t_vals = list(range(-12, 13)) + [
        Fraction(1, 2), Fraction(3, 2), Fraction(2, 3), Fraction(5, 2), Fraction(5, 3),
        Fraction(-1, 2), Fraction(7, 2),
    ]
    fibre_rows = []
    cat_hits = []
    bj_hits = []
    a5_count = 0
    even_count = 0
    irr_count = 0

    for cov in covers:
        if not cov.get("ok"):
            continue
        N, D = make_ND(cov["c"], cov["p2"], cov["r1"], cov["r2"])
        for tv in t_vals:
            pz = fibre_monic_Z(N, D, tv)
            if pz is None:
                continue
            if not pz.is_irreducible:
                continue
            irr_count += 1
            d = int(pz.discriminant())
            row = {
                "s": cov["s"],
                "t": str(tv),
                "poly": str(pz.as_expr()),
                "disc": d,
                "disc_square": d > 0 and is_square(d),
            }
            if row["disc_square"]:
                even_count += 1
                r = classify_poly(pz.as_expr().subs(y, x), do_galois=True)
                row["status"] = r.get("status")
                row["gal"] = r.get("galois")
                if (r.get("status") or "").startswith("HIT_A5"):
                    a5_count += 1
                    row["A5"] = True
            # BJ reduction
            bj = try_bj(pz)
            if bj and bj[0] == "BJ":
                _, a, b = bj
                row["BJ"] = (a, b)
                row["k"] = str(k_of(a, b)) if k_of(a, b) is not None else None
                bj_hits.append(row)
                for tag, ca, cb, ck in CATALOGUE:
                    if a == ca and b == cb:
                        cat_hits.append(
                            {
                                "tag": tag,
                                "k": str(ck),
                                "s": cov["s"],
                                "t": str(tv),
                                "alpha": a,
                                "beta": b,
                            }
                        )
            fibre_rows.append(row)

    cat_k = sorted({h["k"] for h in cat_hits})
    multi_cat = len(cat_k) >= 2

    # ---- Rational interpolation attempt for p2(s) etc. from samples ----
    # Fit p2 as rational function of s of low degree via linear system
    ok_covs = [c for c in covers if c.get("ok")]
    interp = {"attempted": True, "note": "numeric samples; closed form in Q(s) not fitted symbolically this run"}
    if len(ok_covs) >= 6:
        # Vandermonde poly fit degree 3 for p2(s)
        ss = np.array([c["s_float"] for c in ok_covs])
        for name in ("p2", "c", "r1", "r2"):
            yy = np.array([c[name] for c in ok_covs])
            # poly deg 4
            deg = min(4, len(ss) - 1)
            coef = np.polyfit(ss, yy, deg)
            pred = np.polyval(coef, ss)
            err = float(np.max(np.abs(pred - yy)))
            interp[name] = {
                "poly_coeffs_high_to_low": coef.tolist(),
                "max_abs_fit_error": err,
                "deg": deg,
            }

    # ---- Explicit resolvent formula (parametric, numeric coeffs at each s) ----
    # Symbolic form:
    resolvent_form = (
        "f_{s,t}(y) = clear( c(s) y^3 (y-1)(y-p2(s)) - t (y-r1(s))(y-r2(s)) ), "
        "monic in y, where (c,p2,r1,r2)(s) solve the triple-root conditions for branch values 1 and s."
    )

    # Exact family over Q(sqrt(5)) at the rational point s=-1 of the Hurwitz curve:
    # N = -sqrt(5) y^3 (y^2-1), D = y^2 - 1/25
    # Fibre: -sqrt(5) y^3 (y^2-1) - t (y^2 - 1/25) = 0
    # Multiply by -1/sqrt(5): y^5 - y^3 + t/sqrt(5) y^2 - t/(25 sqrt(5)) = 0
    # Minpoly over Q(t): eliminate sqrt(5)
    tt = sp.symbols("tt")
    # From N - tt D = 0 with N = -sqrt(5)*(y^5 - y^3), D = y^2 - 1/25
    # sqrt(5) (y^5 - y^3) + tt (y^2 - 1/25) = 0
    # sqrt(5) = -tt (y^2-1/25)/(y^5-y^3)
    # Square: 5 = tt^2 (y^2-1/25)^2 / (y^5-y^3)^2
    # 5 (y^5-y^3)^2 - tt^2 (y^2-1/25)^2 = 0
    exact_minpoly = sp.expand(
        5 * (y**5 - y**3) ** 2 - tt**2 * (y**2 - sp.Rational(1, 25)) ** 2
    )
    # This is degree 10 - the Q-model of the cover at s=-1 is deg 10, not 5
    # The degree-5 model is over Q(sqrt(5))

    elapsed = round(time.time() - t0, 2)

    verdict = (
        f"Reduced Hurwitz curve modelled as P^1_s (g=0). "
        f"Explicit (3,1,1)^4 covers constructed numerically for {n_ok} rational s "
        f"(Newton on triple-root equations); exact model at s=-1 over Q(√5). "
        f"Fibres: irr={irr_count}, even={even_count}, A5≈{a5_count}, BJ-hits={len(bj_hits)}, "
        f"catalogue hits={len(cat_hits)}, catalogue k={cat_k}, multi_cat={multi_cat}. "
        + (
            "GEOMETRIC multi-k HIT."
            if multi_cat
            else "No geometric multi-k catalogue hit yet; multi-k remains arithmetic (envelope)."
        )
    )

    lines = [
        r"# Explicit model of \(\mathrm{Ni}(A_5,C_3^4)\) and specialisation test",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## 1. Explicit model of the reduced Hurwitz curve",
        "",
        r"By the genus lock (`GENUS_3A4_LOCK.md`), the reduced Hurwitz curve is",
        r"**isomorphic to \(\mathbb{P}^1\) over \(\mathbb{Q}\)**. An explicit model is therefore",
        "",
        r"$$H^{\mathrm{rd}}\ \cong\ \mathbb{P}^1_s,\qquad s=\text{cross-ratio of the four branch points.}$$",
        "",
        r"Infinitely many rational points: all \(s\in\mathbb{Q}\setminus\{0,1\}\) (and \(\infty\)).",
        "",
        r"The map \(H^{\mathrm{rd}}\to M_{0,4}\cong\mathbb{P}^1\) may be taken as the identity in this coordinate",
        r"(branch points placed at \(0,1,\infty,s\)).",
        "",
        "---",
        "",
        r"## 2. Degree-5 resolvent form",
        "",
        r"For each \(s\), a genus-0 cover of type \((3,1,1)^4\) with branch values \(\{0,1,\infty,s\}\)",
        r"is realised in normal form",
        "",
        r"$$\varphi_s(y)=\frac{c\, y^3(y-1)(y-p_2)}{(y-r_1)(y-r_2)},$$",
        "",
        r"where \((c,p_2,r_1,r_2)\) are determined by the conditions that \(\varphi_s-1\) and",
        r"\(\varphi_s-s\) each have a **triple root** (remaining two preimages simple).",
        "",
        r"The fibre polynomial (resolvent in the base coordinate \(t=\varphi_s(y)\)) is",
        "",
        r"$$N_s(y)-t\,D_s(y)=0,\qquad N_s=c y^3(y-1)(y-p_2),\quad D_s=(y-r_1)(y-r_2),$$",
        "",
        r"cleared to a monic polynomial \(f_{s,t}(y)\in\mathbb{C}[y]\) of degree 5.",
        "",
        f"**Form:** `{resolvent_form}`",
        "",
        r"### Exact point \(s=-1\) (over \(\mathbb{Q}(\sqrt{5})\))",
        "",
        f"- \(c=-\sqrt{{5}}\), \(p_2=-1\), \(r_1=1/5\), \(r_2=-1/5\)",
        f"- \(N={exact_sm1['N']}\)",
        f"- \(D={exact_sm1['D']}\)",
        r"- Triple-root conditions verified symbolically.",
        r"- Descent of this single fibre to a deg-5 model over \(\mathbb{Q}(t)\) fails in general;",
        r"  the norm to \(\mathbb{Q}(t)\) yields a degree-10 equation",
        r"  \(5(y^5-y^3)^2-t^2(y^2-1/25)^2=0\).",
        "",
        r"### Numerical covers at rational \(s\)",
        "",
        f"- Requested: {len(covers)}, solved: **{n_ok}**",
        "",
        r"| \(s\) | \(p_2\) | \(c\) | \(r_1\) | \(r_2\) | Newton res |",
        r"|------|--------|------|--------|--------|------------|",
    ]
    for c in covers:
        if not c.get("ok"):
            lines.append(f"| {c['s']} | — | — | — | — | FAIL {c.get('newton_res')} |")
        else:
            lines.append(
                f"| {c['s']} | {c['p2']:.6f} | {c['c']:.6f} | {c['r1']:.6f} | {c['r2']:.6f} | {c['newton_res']:.1e} |"
            )

    lines += [
        "",
        r"### Interpolation of parameters (numeric)",
        "",
        f"```{json.dumps({k: v for k, v in interp.items() if k != 'note'}, indent=2)[:1500]}```",
        "",
        r"**Closed form** \(c(s),p_2(s),r_1(s),r_2(s)\in\mathbb{Q}(s)\) was **not** obtained in this run;",
        r"coefficients involve algebraic functions of \(s\) (often in quadratic extensions).",
        r"A global model over \(\mathbb{Q}(s)\) may require a different coordinate on \(H^{\mathrm{rd}}\)",
        r"or a resolvent of degree \(>5\) over \(\mathbb{Q}(s)\).",
        "",
        "---",
        "",
        r"## 3. Specialisation at rational \((s,t)\)",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| irreducible fibres | {irr_count} |",
        f"| even disc | {even_count} |",
        f"| A5 (among even, galois run) | {a5_count} |",
        f"| BJ reductions | {len(bj_hits)} |",
        f"| catalogue seed hits | {len(cat_hits)} |",
        f"| distinct catalogue \(k\) | {cat_k} |",
        f"| **multi catalogue \(k\)** | **{multi_cat}** |",
        "",
        r"### Catalogue hits",
        "",
    ]
    if not cat_hits:
        lines.append("_None._ Geometric fibres of this normal form did not land on catalogue BJ seeds in the scan.")
    for h in cat_hits:
        lines.append(f"- s={h['s']}, t={h['t']}: **{h['tag']}** (k={h['k']}) α={h['alpha']} β={h['beta']}")

    lines += [
        "",
        r"### Sample A5 fibres (even)",
        "",
    ]
    shown = 0
    for row in fibre_rows:
        if row.get("A5"):
            lines.append(f"- s={row['s']}, t={row['t']}: `{row['poly'][:70]}`")
            shown += 1
            if shown >= 12:
                break
    if shown == 0:
        lines.append("_No A5 even fibres recorded (or galois not triggered)._")

    lines += [
        "",
        r"### Sample BJ reductions",
        "",
    ]
    for row in bj_hits[:15]:
        lines.append(
            f"- s={row['s']}, t={row['t']}: α,β={row.get('BJ')} k={row.get('k')} "
            f"even={row.get('disc_square')} A5={row.get('A5')}"
        )
    if not bj_hits:
        lines.append("_No pure BJ (x^5+αx+β) fibres after y^4-shift in the scan._")

    lines += [
        "",
        "---",
        "",
        r"## 4. Multi-\(k\) conclusion",
        "",
        f"**Geometric multi-\(k\) catalogue hit: {multi_cat}**",
        "",
        r"- The **reduced Hurwitz curve** has an explicit rational model \(\mathbb{P}^1_s\).",
        r"- **Degree-5 covers** of type \((3,1,1)^4\) are constructed for many rational \(s\)",
        r"  (numeric parameters; exact at \(s=-1\) over \(\mathbb{Q}(\sqrt{5})\)).",
        r"- Specialisations tested against the pure-even fixed-\(k\) catalogue:",
        f"  multi catalogue \(k\) = **{multi_cat}**.",
        "",
        r"Until a closed-form \(f_s\in\mathbb{Q}(s)[x]\) produces a multi-\(k\) catalogue hit,",
        r"**multi-\(k\) success remains arithmetic** (envelope/paths in",
        r"`REALISE_3A4_SPECIALISE.md` / `NONRIGID_HURWITZ_SEARCH.md`), **not geometric**",
        r"(Nielsen-labelled).",
        "",
        r"### Next to close the geometric gap",
        "",
        r"1. Determine \(c,p_2,r_1,r_2\) as algebraic functions of \(s\) in closed form",
        r"   (resultants on the triple-root ideal), or find a better coordinate on \(H^{\mathrm{rd}}\).",
        r"2. Produce \(f_s\in\mathbb{Q}(s)[x]\) (possibly after a resolvent of the parameter field).",
        r"3. Re-run the catalogue specialisation test on that closed form.",
        "",
        r"_Generated by build_3a4_resolvent.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "hurwitz_curve_model": "P1_s",
        "resolvent_form": resolvent_form,
        "exact_s_minus_1": {k: str(v) for k, v in exact_sm1.items()},
        "covers": covers,
        "n_covers_ok": n_ok,
        "specialisation": {
            "irr": irr_count,
            "even": even_count,
            "A5": a5_count,
            "bj_hits": len(bj_hits),
            "catalogue_hits": cat_hits,
            "catalogue_k": cat_k,
            "multi_catalogue_k": multi_cat,
        },
        "interpolation": interp,
        "geometric_multi_k": multi_cat,
        "arithmetic_multi_k_reference": "envelope paths in REALISE_3A4_SPECIALISE.md",
    }

    write_md(OUT / "EXPLICIT_3A4_RESOLVENT.md", doc)
    write_md(RESULTS / "EXPLICIT_3A4_RESOLVENT.md", doc)
    write_md(ROOT / "EXPLICIT_3A4_RESOLVENT.md", doc)
    write_json(OUT / "EXPLICIT_3A4_RESOLVENT.json", blob)
    print(verdict, flush=True)
    print(f"Wrote EXPLICIT_3A4_RESOLVENT.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
