"""
Recommended single next move: eliminate triple-root ideal for 3A^4 cover.

Use direct triple-root conditions (lower degree than pure resultants):
  F1(q)=F1'(q)=F1''(q)=0,  Fs(w)=Fs'(w)=Fs''(w)=0
with N=c y^3(y-1)(y-p2), D=(y-r1)(y-r2).

Scale fixes + groebner / successive resultants over Q(s).
Also bivariate algebraic fit from Newton samples.

Output: TRIPLE_ROOT_ELIMINATE.md / .json
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

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402

y, s, t = sp.symbols("y s t")
c, p2, r1, r2, q, w = sp.symbols("c p2 r1 r2 q w")
ra = sp.symbols("ra")


def F_and_derivs(N, D, val, pt):
    G = sp.expand(N - val * D)
    return (
        sp.expand(G.subs(y, pt)),
        sp.expand(sp.diff(G, y).subs(y, pt)),
        sp.expand(sp.diff(G, y, 2).subs(y, pt)),
    )


def main():
    t0 = time.time()
    print("TRIPLE-ROOT ELIMINATION (direct conditions)", flush=True)

    N = c * y**3 * (y - 1) * (y - p2)
    D = (y - r1) * (y - r2)

    e1, e2, e3 = F_and_derivs(N, D, 1, q)
    e4, e5, e6 = F_and_derivs(N, D, s, w)
    eqs = [e1, e2, e3, e4, e5, e6]
    for i, e in enumerate(eqs):
        print(f"  eq{i+1} deg={sp.total_degree(e)} size={len(str(e))}", flush=True)

    blob_models = []

    # ----- Model A: r2=-ra, r1=ra (symmetric poles), known to work at s=-1 -----
    print("  Model A: r2=-r1 ...", flush=True)
    eqs_A = [sp.expand(e.subs({r1: ra, r2: -ra})) for e in eqs]
    # Eliminate q from e1,e2,e3 via successive resultants
    R12 = sp.resultant(eqs_A[0], eqs_A[1], q)
    R23 = sp.resultant(eqs_A[1], eqs_A[2], q)
    Ra = sp.resultant(R12, R23, q) if R12 != 0 and R23 != 0 else sp.resultant(R12, eqs_A[2], q)
    # Eliminate w from e4,e5,e6
    R45 = sp.resultant(eqs_A[3], eqs_A[4], w)
    R56 = sp.resultant(eqs_A[4], eqs_A[5], w)
    Rb = sp.resultant(R45, R56, w) if R45 != 0 and R56 != 0 else 0
    print(f"    Ra size={len(str(Ra))} Rb size={len(str(Rb))}", flush=True)

    # Now Ra, Rb in c, p2, ra, s
    gens = [c, p2, ra, s]
    Ra_n = sp.numer(sp.together(sp.expand(Ra)))
    Rb_n = sp.numer(sp.together(sp.expand(Rb)))
    try:
        Ra_n = sp.Poly(Ra_n, *gens, domain=sp.QQ).primitive()[1].as_expr()
        Rb_n = sp.Poly(Rb_n, *gens, domain=sp.QQ).primitive()[1].as_expr()
    except Exception:
        pass
    print(
        f"    Ra deg={sp.total_degree(Ra_n)} Rb deg={sp.total_degree(Rb_n)}",
        flush=True,
    )

    # Eliminate c between Ra, Rb
    try:
        Rc = sp.resultant(Ra_n, Rb_n, c)
        Rc = sp.factor(sp.expand(Rc))
        print(f"    Res_c(Ra,Rb) size={len(str(Rc))} factor={str(Rc)[:200]}", flush=True)
        # Rc in p2, ra, s — eliminate ra
        Rc_n = sp.numer(sp.together(Rc))
        # Factor and take non-trivial factors
        fac = sp.factor_list(Rc_n)
        factors = [(str(f), m) for f, m in fac[1]]
        print(f"    factors of Rc: {factors[:12]}", flush=True)
        blob_models.append(
            {
                "name": "A_symmetric",
                "Rc_preview": str(Rc)[:400],
                "factors": factors[:20],
            }
        )
        # For each factor involving p2,s only (or p2,ra,s), try to solve
        p2_s_relations = []
        for f, m in fac[1]:
            free = f.free_symbols if hasattr(f, "free_symbols") else set()
            if free <= {p2, s} and f != 0:
                p2_s_relations.append(str(f))
                print(f"    p2-s factor: {f}", flush=True)
            elif free <= {p2, ra, s}:
                # eliminate ra: treat as poly in ra
                try:
                    pr = sp.Poly(f, ra)
                    # discriminant in ra or content as p2,s relation
                    if pr.degree() >= 1:
                        disc_ra = sp.factor(sp.discriminant(pr))
                        p2_s_relations.append(f"disc_ra({f})={disc_ra}")
                        print(f"    disc_ra: {str(disc_ra)[:150]}", flush=True)
                except Exception:
                    pass
        blob_models[-1]["p2_s_relations"] = p2_s_relations[:15]
    except Exception as e:
        print(f"    Model A eliminate c failed: {e}", flush=True)
        blob_models.append({"name": "A_symmetric", "error": str(e)[:200]})

    # Groebner on the two eliminants Ra_n, Rb_n only (lighter)
    print("  Model A groebner on Ra,Rb...", flush=True)
    try:
        GA = sp.groebner([Ra_n, Rb_n], c, p2, ra, s, order="lex", domain=sp.QQ)
        gA = list(GA)
        print(f"    |G|={len(gA)}", flush=True)
        low = []
        for g in gA:
            if len(g.free_symbols) <= 2:
                low.append(str(g)[:250])
                print(f"    low: {str(g)[:120]}", flush=True)
        blob_models.append(
            {
                "name": "A_groebner_RaRb",
                "n": len(gA),
                "preview": [str(g)[:200] for g in gA[:10]],
                "low_var": low,
            }
        )
        for g in gA:
            if g.free_symbols <= {p2, s} and g != 0:
                try:
                    sols = sp.solve(g, p2)
                    blob_models[-1]["p2_sols"] = [str(z) for z in sols[:8]]
                    print(f"    p2 sols: {sols[:4]}", flush=True)
                except Exception as ex:
                    print(f"    solve p2: {ex}", flush=True)
    except Exception as e:
        print(f"    groebner A light failed: {e}", flush=True)
        blob_models.append({"name": "A_groebner_RaRb", "error": str(e)[:200]})

    # ----- Model B: c=1, full 6 eqs, eliminate q,w first -----
    print("  Model B: c=1, eliminate q,w ...", flush=True)
    eqs_B = [sp.expand(e.subs(c, 1)) for e in eqs]
    R12 = sp.resultant(eqs_B[0], eqs_B[1], q)
    R23 = sp.resultant(eqs_B[1], eqs_B[2], q)
    Rq = sp.resultant(R12, R23, q)
    R45 = sp.resultant(eqs_B[3], eqs_B[4], w)
    R56 = sp.resultant(eqs_B[4], eqs_B[5], w)
    Rw = sp.resultant(R45, R56, w)
    Rq = sp.numer(sp.together(sp.expand(Rq)))
    Rw = sp.numer(sp.together(sp.expand(Rw)))
    print(f"    Rq size={len(str(Rq))} Rw size={len(str(Rw))}", flush=True)
    # Eliminate r1 between Rq, Rw
    try:
        R1 = sp.resultant(Rq, Rw, r1)
        R1 = sp.numer(sp.together(sp.expand(R1)))
        print(f"    Res_r1 size={len(str(R1))}", flush=True)
        # Eliminate r2
        R2 = sp.resultant(R1, sp.diff(R1, r2) if R1.has(r2) else R1, r2)
        # Better: Rq,Rw still have r2 — resultant Rq,Rw wrt r2 after r1
        # R1 already eliminated r1, in p2,r2,s
        fac = sp.factor_list(R1)
        factors = [(str(f), m) for f, m in fac[1] if sp.total_degree(f) > 0]
        print(f"    R1 factors ({len(factors)}):", flush=True)
        p2s = []
        for f, m in fac[1]:
            free = f.free_symbols
            if free <= {p2, s}:
                p2s.append(str(f))
                print(f"      p2-s: {f}", flush=True)
            elif free <= {p2, r2, s} and sp.degree(f, r2) >= 1:
                try:
                    d = sp.factor(sp.discriminant(sp.Poly(f, r2)))
                    if d != 0:
                        p2s.append(f"disc_r2={d}")
                        print(f"      disc_r2: {str(d)[:120]}", flush=True)
                except Exception:
                    pass
        blob_models.append(
            {
                "name": "B_c1_eliminate",
                "R1_factors": factors[:25],
                "p2_s_relations": p2s[:15],
            }
        )
        # Try groebner of Rq, Rw only
        GB = sp.groebner([Rq, Rw], p2, r1, r2, s, order="lex", domain=sp.QQ)
        gB = list(GB)
        print(f"    groebner B: {len(gB)}", flush=True)
        lowB = [str(g)[:200] for g in gB if len(g.free_symbols) <= 2]
        for g in lowB[:8]:
            print(f"    lowB: {g[:120]}", flush=True)
        blob_models[-1]["groebner_n"] = len(gB)
        blob_models[-1]["low_var"] = lowB[:15]
        blob_models[-1]["preview"] = [str(g)[:180] for g in gB[:8]]
        for g in gB:
            if g.free_symbols <= {p2, s} and g != 0:
                try:
                    sols = sp.solve(g, p2)
                    blob_models[-1]["p2_sols"] = [str(z) for z in sols[:8]]
                    print(f"    p2 sols: {sols[:4]}", flush=True)
                except Exception as ex:
                    print(f"    p2 solve err {ex}", flush=True)
    except Exception as e:
        print(f"    Model B failed: {e}", flush=True)
        blob_models.append({"name": "B_c1_eliminate", "error": str(e)[:250]})

    # ----- Numeric bivariate fit F(s,p2)=0 -----
    print("  Numeric samples + algebraic curve fit...", flush=True)
    import numpy as np
    from numpy.linalg import norm as npnorm

    def residual(vec, s_val):
        cc, pp, rr1, rr2, qq, ww = vec

        def G_vals(val, pt):
            yy = pt
            Nloc = cc * yy**3 * (yy - 1) * (yy - pp)
            A, Ap, App = yy**3, 3 * yy**2, 6 * yy
            B = yy**2 - (1 + pp) * yy + pp
            Bp = 2 * yy - (1 + pp)
            Np = cc * (Ap * B + A * Bp)
            Npp = cc * (App * B + 2 * Ap * Bp + A * 2.0)
            Dd = (yy - rr1) * (yy - rr2)
            Dp = 2 * yy - (rr1 + rr2)
            return [Nloc - val * Dd, Np - val * Dp, Npp - val * 2.0]

        return np.array(G_vals(1.0, qq) + G_vals(s_val, ww), dtype=float)

    def newton(s_val, x0, niter=80):
        v = np.array(x0, dtype=float)
        for _ in range(niter):
            r = residual(v, s_val)
            if npnorm(r) < 1e-14:
                return v, True
            J = np.zeros((6, 6))
            eps = 1e-8
            for j in range(6):
                dv = np.zeros(6)
                dv[j] = eps
                J[:, j] = (residual(v + dv, s_val) - r) / eps
            try:
                step = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -r, rcond=None)[0]
            v = v + step
        return v, npnorm(residual(v, s_val)) < 1e-12

    rng = np.random.default_rng(1)
    samples = []
    v = np.array([-np.sqrt(5), -1.0, 0.2, -0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5])
    for num, den in [
        (-4, 1), (-3, 1), (-2, 1), (-1, 1), (-1, 2), (1, 2), (3, 2), (2, 1),
        (5, 2), (3, 1), (4, 1), (5, 1), (-5, 2), (5, 3), (7, 2), (2, 3), (4, 3),
        (7, 3), (8, 5), (-3, 2), (5, 4), (9, 2),
    ]:
        sv = float(Fraction(num, den))
        if abs(sv) < 1e-12 or abs(sv - 1) < 1e-12:
            continue
        v, ok = newton(sv, v)
        if not ok:
            for _ in range(50):
                v2, ok = newton(sv, rng.normal(size=6))
                if ok:
                    v = v2
                    break
        if ok:
            samples.append({"s": sv, "p2": float(v[1]), "c": float(v[0]), "r1": float(v[2]), "r2": float(v[3])})
    print(f"  samples: {len(samples)}", flush=True)

    fit_info = {}
    if len(samples) >= 8:
        di, dj = 4, 4
        rows = []
        for sm in samples:
            sv, pv = sm["s"], sm["p2"]
            row = [sv**i * pv**j for i in range(di + 1) for j in range(dj + 1)]
            rows.append(row)
        Mf = np.array(rows, dtype=float)
        _, S, Vt = np.linalg.svd(Mf, full_matrices=True)
        null_vec = Vt[-1, :]
        res = float(np.max(np.abs(Mf @ null_vec)))
        # rationalize
        scale = null_vec[np.argmax(np.abs(null_vec))]
        nv = null_vec / scale
        ss, pp = sp.symbols("s p2")
        poly = 0
        terms = []
        idx = 0
        for i in range(di + 1):
            for j in range(dj + 1):
                aij = sp.nsimplify(nv[idx], tolerance=2e-7, rational=True)
                idx += 1
                if aij != 0:
                    terms.append((i, j, str(aij)))
                    poly += aij * ss**i * pp**j
        poly = sp.expand(poly)
        errs = [abs(complex(poly.subs({ss: sm["s"], pp: sm["p2"]}))) for sm in samples]
        fit_info = {
            "bideg": [di, dj],
            "svd_tail": [float(x) for x in S[-5:]],
            "max_residual_float_null": res,
            "rational_poly": str(poly),
            "factored": str(sp.factor(poly)),
            "terms": terms,
            "max_eval_err": float(max(errs)) if errs else None,
        }
        print(f"  fit poly={sp.factor(poly)}", flush=True)
        print(f"  max eval err={max(errs):.3e}", flush=True)

        # If good fit, solve p2(s) and try to recover c,r1,r2 symbolically at generic s
        if max(errs) < 1e-5 and poly != 0:
            print("  attempting p2(s) from fit + reconstruct cover...", flush=True)
            try:
                p2sols = sp.solve(poly, pp)
                fit_info["p2_branches"] = [str(z) for z in p2sols[:6]]
                print(f"    p2 branches: {p2sols[:4]}", flush=True)
            except Exception as e:
                fit_info["p2_solve_error"] = str(e)[:100]

    # Exact s=-1 resolvent
    s5 = sp.sqrt(5)
    exact = {
        "s": -1,
        "c": str(-s5),
        "p2": -1,
        "r1": "1/5",
        "r2": "-1/5",
        "f_over_Qsqrt5": str(sp.expand(y**5 - y**3 + (t / s5) * y**2 - t / (25 * s5))),
        "norm_deg10": str(sp.expand(5 * (y**5 - y**3) ** 2 - t**2 * (y**2 - sp.Rational(1, 25)) ** 2)),
    }

    # Envelope control
    ku = Fraction(-8, 5) + t * (Fraction(4, 5) - Fraction(-8, 5))
    alpha = sp.together(25 - 3125 * ku**4 / 256)
    beta = sp.together(ku * alpha)
    env = []
    for tv, tag, a0, b0, k0 in [
        (0, "flagship", -55, 88, "-8/5"),
        (1, "classical", 20, 16, "4/5"),
        (Fraction(1, 3), "classical_m", 20, -16, "-4/5"),
    ]:
        a = int(sp.Rational(sp.simplify(alpha.subs(t, tv))))
        b = int(sp.Rational(sp.simplify(beta.subs(t, tv))))
        if (a, b) == (a0, b0):
            env.append({"tag": tag, "k": k0, "t": str(tv)})

    elapsed = round(time.time() - t0, 2)
    closed = None
    if fit_info.get("max_eval_err") is not None and fit_info["max_eval_err"] < 1e-5:
        closed = {
            "type": "plane_model_F_s_p2",
            "poly": fit_info.get("factored") or fit_info.get("rational_poly"),
            "err": fit_info["max_eval_err"],
            "p2_branches": fit_info.get("p2_branches"),
        }
    for m in blob_models:
        if m.get("p2_sols"):
            closed = {"type": "groebner_p2_s", "sols": m["p2_sols"], "model": m["name"]}
        if m.get("p2_s_relations"):
            if closed is None:
                closed = {"type": "relations", "rels": m["p2_s_relations"][:5], "model": m["name"]}

    verdict = (
        f"Elimination complete ({elapsed}s). "
        f"Models: {len(blob_models)}. "
        f"Closed-form candidate: {closed is not None} — {closed}. "
        f"Numeric samples={len(samples)}, fit_err={fit_info.get('max_eval_err')}. "
        f"Exact s=-1 over Q(√5) retained. Envelope multi-k: {env}."
    )

    lines = [
        r"# Triple-root elimination — recommended next move",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Method",
        "",
        r"Direct triple-root conditions (degree \(\le 5\) in auxiliaries \(q,w\)):",
        r"\(F_1(q)=F_1'(q)=F_1''(q)=0\), \(F_s(w)=F_s'(w)=F_s''(w)=0\),",
        r"then successive resultants / Gröbner to eliminate \(q,w,c,r_i\).",
        "",
        "---",
        "",
        r"## Algebraic models",
        "",
    ]
    for m in blob_models:
        lines.append(f"### `{m.get('name')}`")
        for k, val in m.items():
            if k == "name":
                continue
            lines.append(f"- **{k}**: `{val}`")
        lines.append("")

    lines += [
        "---",
        "",
        r"## Plane model fit \(F(s,p_2)=0\)",
        "",
        f"```{fit_info}```",
        "",
        "---",
        "",
        r"## Exact \(s=-1\)",
        "",
        f"```{exact}```",
        "",
        r"## Envelope multi-\(k\) control",
        "",
        f"{env}",
        "",
        "---",
        "",
        r"## Closed-form status",
        "",
        f"**Candidate:** `{closed}`",
        "",
        r"### If plane model \(F(s,p_2)=0\) is reliable",
        r"1. Solve \(p_2=p_2(s)\) on each branch.",
        r"2. Back-substitute into triple-root eqs for \(c,r_1,r_2\).",
        r"3. Form \(f_{s,t}(y)=\mathrm{monic}(N_s-t D_s)\).",
        r"4. Specialise and test fixed-\(k\) catalogue.",
        "",
        r"### Geometric multi-\(k\)",
        r"Still open until a \(\mathbb{Q}(s)\)-resolvent is assembled and catalogue-tested.",
        r"Arithmetic multi-\(k\) remains the envelope.",
        "",
        r"_Generated by triple_root_eliminate.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "models": blob_models,
        "fit": fit_info,
        "exact_s_minus1": exact,
        "closed_form_candidate": closed,
        "envelope_hits": env,
        "n_samples": len(samples),
        "samples": samples[:15],
    }
    write_md(OUT / "TRIPLE_ROOT_ELIMINATE.md", doc)
    write_md(RESULTS / "TRIPLE_ROOT_ELIMINATE.md", doc)
    write_md(ROOT / "TRIPLE_ROOT_ELIMINATE.md", doc)
    write_json(OUT / "TRIPLE_ROOT_ELIMINATE.json", blob)
    print(verdict, flush=True)
    print(f"Wrote TRIPLE_ROOT_ELIMINATE.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
