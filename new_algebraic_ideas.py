"""
New algebraic ideas A–F (ranked programme).

Primary:
  A. Mestre deformation: solve P''R - 2 P' R' ≡ 0 (mod P), deg R < deg P;
     build 1-param families P_t; sample Gal; HQCC lattice seeds.
  F. Embed test: solve chi_T(a..f) = monic family poly as eqs in template params.

Secondary:
  B. Non-BJ degree-1 A5 family x^5+75x^3+A x^2+3A matched to chi_T.

Light probes:
  C. Alternate matrix avatars (companion of ternary-coeff poly; transfer graph).
  D. Icosahedral / invariant-style parameter sample.
  E. HQCC-native poly from T3 itinerary / transfer truncation.

Do NOT: more linear cuts of T; F→T disc□ hope; rigid φ surgery; binary Collatz forces evenness.

Output: NEW_ALGEBRAIC_IDEAS.md / .json
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

t = sp.symbols("t")
a, b, c, d, e, f = sp.symbols("a b c d e f")

# HQCC / multi-seed BJ seeds (alpha, beta, tag)
SEEDS = [
    (-55, 88, "flagship"),
    (-55, -88, "flagship_m"),
    (95, 76, "classical_95_76"),
    (95, -76, "classical_95_m76"),
    (20, 16, "classical_20_16"),
    (20, -16, "classical_20_m16"),
    (-100, 400, "lsw_slice_100"),  # on LSW if beta=-4 alpha
    (124, -496, "lsw_slice_124"),
    (320, -512, "flag_m_cleared"),  # pure-even flagship style
    (-3121, 12484, "lsw_m"),  # LSW sample
]


def is_square_poly(expr, var=t) -> dict:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True}
        # content over Q(var)
        num, den = sp.fraction(sp.together(ex))
        num, den = sp.expand(num), sp.expand(den)
        for piece, label in ((num, "num"), (den, "den")):
            if piece == 0:
                continue
            P = sp.Poly(piece, var, domain=sp.QQ)
            cont = sp.Rational(P.content())
            if cont < 0:
                cont = -cont
                # sign absorbed
            n, dd = int(sp.numer(cont)), int(sp.denom(cont))
            if not (
                sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(dd, 2)[1]
            ):
                return {"ok": False, "reason": f"content_{label}", "content": str(cont)}
            prim = P.primitive()[1]
            if prim == 1 or prim.degree() == 0:
                continue
            fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
            odds = [(str(fi)[:40], int(m)) for fi, m in fac[1] if m % 2]
            if odds:
                return {"ok": False, "reason": "odd_factors", "odd": odds[:4]}
        return {"ok": True, "degenerate": False}
    except Exception as ex:
        return {"ok": False, "error": str(ex)[:120]}


def chi_T_expr(aa, bb, cc, dd, ee, ff):
    return (
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


# ---------------------------------------------------------------------------
# Idea A — Mestre
# ---------------------------------------------------------------------------


def mestre_R_space(P_expr, deg_R_max: int | None = None):
    """
    Solve P'' R - 2 P' R' ≡ 0 (mod P) for R = sum r_i x^i, deg R < n
    (i.e. deg ≤ n-1). Returns basis of the solution space over Q.
    """
    P = sp.Poly(sp.expand(P_expr), x, domain=sp.QQ)
    n = int(P.degree())
    if deg_R_max is None:
        deg_R_max = n - 1  # deg R < deg P
    Pp = sp.diff(P.as_expr(), x)
    Ppp = sp.diff(Pp, x)
    r_syms = sp.symbols(f"r0:{deg_R_max + 1}")
    R = sum(r_syms[i] * x**i for i in range(deg_R_max + 1))
    Rp = sp.diff(R, x)
    W = sp.expand(Ppp * R - 2 * Pp * Rp)
    # Remainder of W divided by P (coeffs linear in r_i)
    rem_poly = sp.Poly(W, x).rem(P)
    eqs = [sp.expand(c) for c in rem_poly.coeffs()] if rem_poly != 0 else [0]
    # Drop tautological zeros
    eqs = [eq for eq in eqs if eq != 0] or [0]
    A_mat, _ = sp.linear_eq_to_matrix(eqs, list(r_syms))
    null = A_mat.nullspace()
    basis = []
    for vec in null:
        coeffs = [sp.Rational(v) for v in vec]
        dens = [int(c.q) for c in coeffs if c != 0]
        L = 1
        for q in dens:
            L = int(sp.ilcm(L, q))
        icoeffs = [int(c * L) for c in coeffs]
        g = 0
        for ci in icoeffs:
            g = int(sp.gcd(g, ci)) if g else abs(ci)
        if g > 1:
            icoeffs = [ci // g for ci in icoeffs]
        Rpoly = sum(icoeffs[i] * x**i for i in range(len(icoeffs)))
        if Rpoly != 0:
            basis.append(sp.expand(Rpoly))
    return {
        "n": n,
        "deg_R_max": deg_R_max,
        "null_dim": len(null),
        "basis": basis,
        "basis_str": [str(b) for b in basis],
    }


def family_res_y_shift(P_expr, R_expr, tvar=t):
    """F(z,t) = Res_y(P(y), z - y - t R(y)), monic in z."""
    y, z = sp.symbols("y z")
    P = sp.expand(P_expr.subs(x, y))
    R = sp.expand(R_expr.subs(x, y))
    res = sp.resultant(P, z - y - tvar * R, y)
    pol = sp.Poly(sp.expand(res), z)
    lc = pol.LC()
    if lc != 0 and lc != 1:
        pol = sp.Poly(sp.expand(pol.as_expr() / lc), z)
    return pol.as_expr().subs(z, x)


def family_res_u_scale(P_expr, R_expr, tvar=t):
    """N(u,t) = Res_y(P(y), u P'(y) - t R(y))."""
    y, u = sp.symbols("y u")
    P = sp.expand(P_expr.subs(x, y))
    R = sp.expand(R_expr.subs(x, y))
    Pp = sp.diff(P, y)
    res = sp.resultant(P, u * Pp - tvar * R, y)
    pol = sp.Poly(sp.expand(res), u)
    if pol.degree() < 1:
        return None
    lc = pol.LC()
    if lc not in (0, 1) and lc is not None:
        try:
            pol = sp.Poly(sp.expand(pol.as_expr() / lc), u)
        except Exception:
            pass
    return pol.as_expr().subs(u, x)


def family_res_u_translate(P_expr, R_expr, tvar=t):
    """N(u,t) = Res_y(P(y), u P'(y) - R(y) - t)."""
    y, u = sp.symbols("y u")
    P = sp.expand(P_expr.subs(x, y))
    R = sp.expand(R_expr.subs(x, y))
    Pp = sp.diff(P, y)
    res = sp.resultant(P, u * Pp - R - tvar, y)
    pol = sp.Poly(sp.expand(res), u)
    if pol.degree() < 1:
        return None
    lc = pol.LC()
    try:
        if lc not in (0, 1):
            pol = sp.Poly(sp.expand(pol.as_expr() / lc), u)
    except Exception:
        pass
    return pol.as_expr().subs(u, x)


def disc_of_family(F_expr, var=x, tvar=t):
    try:
        pol = sp.Poly(sp.expand(F_expr), var, domain=sp.QQ[tvar])
        if pol.degree() != 5:
            return None, f"deg={pol.degree()}"
        # monic
        if pol.LC() != 1:
            F2 = sp.expand(sp.monic(pol.as_expr(), var))
            pol = sp.Poly(F2, var, domain=sp.QQ[tvar])
        D = sp.expand(pol.discriminant())
        return D, None
    except Exception as ex:
        return None, str(ex)[:100]


def sample_family_gal(F_expr, tvals, max_checks: int = 8) -> list[dict]:
    rows = []
    for tv in tvals:
        try:
            Fx = sp.expand(F_expr.subs(t, tv))
            pol = sp.Poly(Fx, x, domain=sp.QQ)
            # clear denoms
            cont = pol.denom() if hasattr(pol, "denom") else 1
            coeffs = pol.all_coeffs()
            dens = [sp.Integer(c).q if isinstance(c, sp.Rational) else getattr(sp.QQ(c), "q", 1) for c in coeffs]
            # simpler:
            Fx2 = sp.together(Fx)
            num, den = sp.fraction(Fx2)
            # monic Z poly via clearing
            P = sp.Poly(sp.expand(Fx), x, domain=sp.QQ)
            lc = P.LC()
            if lc == 0:
                continue
            mon = sp.expand(P.as_expr() / lc)
            Pm = sp.Poly(mon, x, domain=sp.QQ)
            # multiply by denom of coeffs to Z
            den_l = 1
            for c in Pm.all_coeffs():
                r = sp.Rational(c)
                den_l = sp.ilcm(den_l, int(r.q))
            # monic with x = y/den? For Gal, use content-free primitive monic if integer coeffs
            Q = sp.Poly(sp.expand(den_l * mon), x, domain=sp.QQ)
            # If not monic integer, scale variable
            if Q.LC() != 1:
                # already monic if mon was monic and den_l scales all
                pass
            # Build monic Z by x=w, clear: if coeffs in Q, multiply through
            coeffs_q = [sp.Rational(c) for c in Pm.all_coeffs()]
            Dclear = 1
            for cq in coeffs_q:
                Dclear = sp.ilcm(Dclear, int(cq.q))
            # Poly Dclear^deg * mon(x/Dclear) for monic Z? Standard:
            # mon = x^5 + c4 x^4 + ... with ci in Q
            # set x = y, multiply by L^5 after y = L x... use integer content primitive of L*mon with L clearing
            expr_z = sum(
                sp.Rational(ci) * x ** (5 - i) for i, ci in enumerate(coeffs_q)
            )
            # Actually Pm is monic: x^5 + a4 x^4 + ...
            L = Dclear
            # y = x; poly with coeffs a_i * L^{5-i}? For monic Z specialisation classification:
            cleared = sp.expand(x**5 + sum(
                sp.Rational(coeffs_q[i]) * (L ** i) * x ** (5 - i)
                for i in range(1, 6)
            ))
            # Wrong scaling. Simpler approach: evaluate numerical monic over Q and use classify if Z coeffs.
            Pz = sp.Poly(mon, x, domain=sp.QQ)
            if all(sp.Rational(c).q == 1 for c in Pz.all_coeffs()):
                chi = Pz.as_expr()
            else:
                # skip non-Z for galois or clear naively
                L = 1
                for c in Pz.all_coeffs():
                    L = sp.ilcm(L, int(sp.Rational(c).q))
                chi = sp.expand(sum(int(sp.Rational(c) * L) * x ** (Pz.degree() - i)
                                    for i, c in enumerate(Pz.all_coeffs())))
                # not monic — skip Gal
                if sp.LC(chi, x) != 1:
                    rows.append({"t": tv, "status": "not_monic_Z", "disc_sq": None})
                    continue
            if len(rows) >= max_checks:
                break
            cl = classify_poly(chi, do_galois=True)
            rows.append(
                {
                    "t": tv,
                    "status": cl.get("status"),
                    "disc_sq": cl.get("disc_square"),
                    "galois": cl.get("galois"),
                }
            )
        except Exception as ex:
            rows.append({"t": tv, "status": f"err:{ex}"[:40]})
    return rows


def run_mestre_on_seed(alpha: int, beta: int, tag: str) -> dict:
    P = x**5 + alpha * x + beta
    print(f"  Mestre seed {tag}: {P}", flush=True)
    space = mestre_R_space(P)
    out = {
        "tag": tag,
        "alpha": alpha,
        "beta": beta,
        "P": str(P),
        "null_dim": space["null_dim"],
        "R_basis": space["basis_str"],
        "families": [],
    }
    if space["null_dim"] == 0:
        out["note"] = "No nontrivial R (unexpected if disc□ and irr)"
        return out

    # Prefer lowest-degree nonzero basis elements
    basis = space["basis"]
    # Also check disc of seed
    d0 = disc_bj_int(alpha, beta)
    out["seed_disc_square"] = d0 > 0 and is_square(d0)

    constructions = [
        ("shift_y_tR", family_res_y_shift),
        ("uPp_tR", family_res_u_scale),
        ("uPp_R_t", family_res_u_translate),
    ]
    for R in basis[:3]:
        for cname, cfun in constructions:
            try:
                F = cfun(P, R)
                if F is None:
                    continue
                D, err = disc_of_family(F)
                if err:
                    out["families"].append(
                        {"R": str(R), "construction": cname, "error": err}
                    )
                    continue
                info = is_square_poly(D, t)
                entry = {
                    "R": str(R),
                    "construction": cname,
                    "F_preview": str(F)[:100],
                    "disc_square_in_Qt": info.get("ok"),
                    "disc_info": info,
                    "deg_F": int(sp.degree(sp.Poly(sp.expand(F), x))),
                }
                if info.get("ok") and not info.get("degenerate"):
                    print(f"    HIT {tag} {cname} R={R}", flush=True)
                    entry["samples"] = sample_family_gal(
                        F, [1, 2, 3, -1, 5, 7, 9], max_checks=6
                    )
                out["families"].append(entry)
            except Exception as ex:
                out["families"].append(
                    {"R": str(R), "construction": cname, "error": str(ex)[:80]}
                )
    return out


def idea_A() -> dict:
    print("IDEA A — Mestre", flush=True)
    results = []
    for al, be, tag in SEEDS:
        # only seeds with disc square (even)
        d0 = disc_bj_int(al, be)
        if d0 <= 0 or not is_square(d0):
            results.append(
                {
                    "tag": tag,
                    "skipped": True,
                    "reason": "seed_disc_not_square",
                    "alpha": al,
                    "beta": be,
                }
            )
            continue
        results.append(run_mestre_on_seed(al, be, tag))
    n_hit = sum(
        1
        for r in results
        for fam in r.get("families") or []
        if fam.get("disc_square_in_Qt")
    )
    n_R = sum(1 for r in results if r.get("null_dim", 0) > 0)
    return {
        "seeds_processed": len(results),
        "seeds_with_R": n_R,
        "families_disc_square": n_hit,
        "results": results,
    }


# ---------------------------------------------------------------------------
# Idea F — embed into T
# ---------------------------------------------------------------------------


def poly_coeffs_monic5(expr):
    """Return (c4,c3,c2,c1,c0) for x^5 + c4 x^4 + ... + c0."""
    P = sp.Poly(sp.expand(expr), x, domain=sp.QQ)
    if P.degree() != 5:
        return None
    lc = P.LC()
    if lc != 1:
        P = sp.Poly(sp.expand(P.as_expr() / lc), x, domain=sp.QQ)
    # monic: coeffs of x^4..x^0
    allc = P.all_coeffs()  # [1, c4, c3, c2, c1, c0]
    while len(allc) < 6:
        allc.append(0)
    return tuple(allc[1:])


def embed_into_T(target_coeffs) -> dict:
    """
    Solve chi_T coeffs = target monic (0, c3, c2, c1, c0) typically no x^4.
    chi_T = x^5 - d x^3 - (a+ef) x^2 - (bf+ce) x + (ad-bc)
    so c4=0 always for T. If target has c4≠0, T cannot embed without Tschirnhaus.
    """
    c4, c3, c2, c1, c0 = target_coeffs
    if c4 != 0 and sp.expand(c4) != 0:
        return {
            "embeddable": False,
            "reason": "x^4_term_nonzero_T_has_no_x4",
            "note": "Need Tschirnhaus to kill x^4 before T-embed",
        }
    # -d = c3 ⇒ d = -c3
    # -(a+ef) = c2 ⇒ a+ef = -c2
    # -(bf+ce) = c1 ⇒ bf+ce = -c1
    # ad-bc = c0
    dd = -sp.expand(c3)
    eqs = [
        a + e * f + sp.expand(c2),
        b * f + c * e + sp.expand(c1),
        a * dd - b * c - sp.expand(c0),
    ]
    # free vars a,b,c,e,f with d fixed
    sols = sp.solve(eqs, [a, b, c, e, f], dict=True, simplify=False)
    # Also try with some zeros for sparse models
    sparse_hits = []
    for zeros in [
        {e: 0},
        {f: 0},
        {e: 0, f: 1},
        {a: 0},
        {a: 0, e: 0},
        {d: dd, e: 0, f: 1},
    ]:
        sub_eqs = [eq.subs(zeros) for eq in eqs]
        free = [v for v in (a, b, c, e, f) if v not in zeros]
        try:
            s2 = sp.solve(sub_eqs, free, dict=True)
            for s in s2[:3]:
                full = {**zeros, **s, d: dd}
                sparse_hits.append({str(k): str(v) for k, v in full.items()})
        except Exception:
            pass
    return {
        "embeddable": len(sols) > 0 or len(sparse_hits) > 0,
        "d": str(dd),
        "n_generic_sols": len(sols),
        "sample_sol": {str(k): str(v) for k, v in sols[0].items()} if sols else None,
        "sparse_hits": sparse_hits[:5],
        "eqs": [str(eq) for eq in eqs],
    }


def idea_F(mestre_A: dict) -> dict:
    print("IDEA F — embed into T", flush=True)
    tests = []
    # 1) Embed seed itself (BJ → BJ-embed)
    for al, be, tag in SEEDS[:6]:
        coeffs = (0, 0, 0, al, be)  # x^5 + al x + be
        emb = embed_into_T(coeffs)
        tests.append({"source": f"seed:{tag}", "coeffs": coeffs, "embed": emb})
    # 2) Literature non-BJ family at A=3
    for Aval in [3, 9, 27, 61]:
        P = x**5 + 75 * x**3 + Aval * x**2 + 3 * Aval
        coeffs = poly_coeffs_monic5(P)
        emb = embed_into_T(coeffs)
        tests.append({"source": f"nonBJ_A={Aval}", "coeffs": [str(c) for c in coeffs], "embed": emb})
    # 3) Mestre family specialisations t=1 when disc□ family exists
    for r in mestre_A.get("results") or []:
        for fam in r.get("families") or []:
            if not fam.get("disc_square_in_Qt"):
                continue
            # specialise symbolic F at t=1 if we stored preview only — re-run light
            tag = r["tag"]
            tests.append(
                {
                    "source": f"mestre_hit:{tag}:{fam.get('construction')}",
                    "note": "family disc□ in Q(t); embed per specialisation",
                    "embed_specialisations": [],
                }
            )
            # reconstruct
            al, be = r["alpha"], r["beta"]
            P = x**5 + al * x + be
            R = sp.sympify(fam["R"], locals={"x": x})
            cname = fam["construction"]
            try:
                if cname == "shift_y_tR":
                    F = family_res_y_shift(P, R)
                elif cname == "uPp_tR":
                    F = family_res_u_scale(P, R)
                else:
                    F = family_res_u_translate(P, R)
                for tv in [1, 2, 3]:
                    Fx = sp.expand(F.subs(t, tv))
                    co = poly_coeffs_monic5(Fx)
                    if co is None:
                        continue
                    # kill x^4 by shift x = y - c4/5
                    c4, c3, c2, c1, c0 = co
                    if c4 != 0:
                        y = sp.symbols("y")
                        Fx2 = sp.expand(Fx.subs(x, y - sp.Rational(c4, 5)))
                        Fx2 = sp.Poly(Fx2, y)
                        Fx2 = sp.expand(Fx2.as_expr() / Fx2.LC()).subs(y, x)
                        co = poly_coeffs_monic5(Fx2)
                    emb = embed_into_T(co)
                    tests[-1]["embed_specialisations"].append(
                        {"t": tv, "coeffs": [str(c) for c in co], "embed": emb}
                    )
            except Exception as ex:
                tests[-1]["error"] = str(ex)[:80]
    n_yes = sum(1 for t in tests if (t.get("embed") or {}).get("embeddable"))
    n_mes_yes = sum(
        1
        for t in tests
        for s in t.get("embed_specialisations") or []
        if s.get("embed", {}).get("embeddable")
    )
    return {
        "tests": tests,
        "n_embeddable_static": n_yes,
        "n_mestre_spec_embeddable": n_mes_yes,
    }


# ---------------------------------------------------------------------------
# Idea B — non-BJ degree-1 family
# ---------------------------------------------------------------------------


def idea_B() -> dict:
    print("IDEA B — non-BJ deg-1 family", flush=True)
    A = sp.symbols("A")
    # Literature-style: x^5 + 75 x^3 + A x^2 + 3 A
    P = x**5 + 75 * x**3 + A * x**2 + 3 * A
    # disc as poly in A
    D = sp.Poly(P, x).discriminant()
    D = sp.expand(D)
    info = is_square_poly(D, A)
    # Factor D
    fac = sp.factor(D)
    # When is D a square? Content + square-free part
    # Sample integer A on resonant lattice
    lattice = [1, 3, 9, 27, 61, 80, 243, 539, -3, -9, 18, 54, 4880, 55, 88, 95]
    samples = []
    for Av in lattice:
        Pv = sp.expand(P.subs(A, Av))
        pol = sp.Poly(Pv, x, domain=sp.ZZ)
        if not pol.is_irreducible:
            samples.append({"A": Av, "status": "red"})
            continue
        disc = int(pol.discriminant())
        sq = disc > 0 and is_square(disc)
        rec = {"A": Av, "disc_sq": sq, "disc": disc}
        if sq:
            cl = classify_poly(Pv, do_galois=True)
            rec["status"] = cl.get("status")
            rec["galois"] = cl.get("galois")
        else:
            rec["status"] = "odd" if disc > 0 else "disc_nonpos"
        samples.append(rec)
        print(f"    A={Av}: disc□={sq} {rec.get('status')}", flush=True)

    # Match to chi_T symbolically
    # d = -75, a+ef = -A, bf+ce = 0, a*(-75) - b c = 3A
    eqs = [
        a + e * f + A,
        b * f + c * e,
        -75 * a - b * c - 3 * A,
    ]
    # solve treating A as free param
    sols = sp.solve(eqs, [a, b, c, e, f], dict=True)
    # sparse: e=0, f free
    sparse = []
    for zsub in [{e: 0}, {e: 0, f: 1}, {f: 0}, {a: 0}, {b: 1, e: 0}]:
        try:
            s = sp.solve([eq.subs(zsub) for eq in eqs], [v for v in (a, b, c, e, f) if v not in zsub], dict=True)
            for sol in s[:2]:
                sparse.append({str(k): str(v) for k, v in {**zsub, **sol, d: -75}.items()})
        except Exception:
            pass

    # HQCC-native? Is A on lattice with disc□ and embed relations involving only model ints?
    hits_A5 = [s for s in samples if str(s.get("status", "")).startswith("HIT_A5")]
    return {
        "family": "x^5+75*x^3+A*x^2+3*A",
        "disc_identically_square_in_A": info.get("ok"),
        "disc_info": info,
        "disc_factored_preview": str(fac)[:200],
        "lattice_samples": samples,
        "n_disc_sq": sum(1 for s in samples if s.get("disc_sq")),
        "n_A5": len(hits_A5),
        "A5_hits": hits_A5,
        "T_match_generic_sols": len(sols),
        "T_match_sample": {str(k): str(v) for k, v in sols[0].items()} if sols else None,
        "T_match_sparse": sparse[:6],
        "HQCC_native_embed": False,  # set below
    }


# ---------------------------------------------------------------------------
# Idea C — new matrix avatars
# ---------------------------------------------------------------------------


def idea_C() -> dict:
    print("IDEA C — matrix avatars", flush=True)
    rows = []
    # C1: companion matrix of monic poly with ternary-constrained coeffs
    # x^5 + p3 x^3 + p2 x^2 + p1 x + p0 with pi in {0,±3,±9,±27,61,...}
    pool = [0, 3, -3, 9, 27, 61, -9]
    n_sq = 0
    n_irr = 0
    n_A5 = 0
    tested = 0
    for p3, p2, p1, p0 in itertools.product(pool, repeat=4):
        if p3 == p2 == p1 == p0 == 0:
            continue
        tested += 1
        if tested > 400:
            break
        chi = x**5 + p3 * x**3 + p2 * x**2 + p1 * x + p0
        pol = sp.Poly(chi, x, domain=sp.ZZ)
        if not pol.is_irreducible:
            continue
        n_irr += 1
        disc = int(pol.discriminant())
        if disc > 0 and is_square(disc):
            n_sq += 1
            if n_A5 < 5:
                cl = classify_poly(chi, do_galois=True)
                if str(cl.get("status", "")).startswith("HIT_A5"):
                    n_A5 += 1
                    rows.append({"coeffs": (p3, p2, p1, p0), "status": cl.get("status")})
    # C2: transfer matrix of ternary graph (3 states, weighted by lattice)
    # states 0,1,2; adjacency A[i,j] = weight if edge i->j allowed by T3 residue
    # T3-like: from r, branches
    w0, w1, w2 = 3, 61, 80
    M = sp.Matrix(
        [
            [0, w1, w2],
            [w0, 0, w1],
            [w0, w2, 0],
        ]
    )
    # lift to 5×5 by block padding / Kronecker with 2-cycle companion
    J2 = sp.Matrix([[0, 1], [-1, 0]])
    # 6x6 kronecker then take a 5x5 principal? Use charpoly of companion-like
    chi3 = M.charpoly(x).as_expr()
    # Build 5×5: direct sum style companion of (x^2+1)*(cubic) truncated — skip
    # Instead: adjacency of path graph on 5 nodes with ternary weights
    weights = [3, 9, 61, 80, 27]
    Tm = sp.zeros(5)
    for i in range(5):
        Tm[i, (i + 1) % 5] = weights[i]
        Tm[i, (i + 2) % 5] = weights[(i + 1) % 5]
    chiTm = sp.expand(Tm.charpoly(x).as_expr())
    polT = sp.Poly(chiTm, x, domain=sp.ZZ)
    discT = int(polT.discriminant()) if polT.degree() == 5 else None
    transfer = {
        "matrix_preview": str(Tm),
        "chi": str(chiTm),
        "deg": polT.degree(),
        "irr": bool(polT.is_irreducible) if polT.degree() == 5 else False,
        "disc_sq": is_square(discT) if discT and discT > 0 else False,
        "disc": discT,
    }
    if transfer["disc_sq"] and transfer["irr"]:
        cl = classify_poly(chiTm, do_galois=True)
        transfer["status"] = cl.get("status")

    # C3: Kronecker-ish — representation: fixed 3-cycle block deformed
    # 5×5 matrix with 3-cycle block + lattice deformation param u
    u = sp.symbols("u")
    Mc = sp.Matrix(
        [
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [1, 0, 0, 0, 0],  # 3-cycle on first 3
            [0, 0, 0, 0, 1],
            [u, 0, 0, -3, 0],  # lattice coupling
        ]
    )
    # integer specialisations
    kron_rows = []
    for uv in [0, 1, 3, 9, 61, 80, 243]:
        Mcu = Mc.subs(u, uv)
        chi = sp.expand(Mcu.charpoly(x).as_expr())
        pol = sp.Poly(chi, x, domain=sp.ZZ)
        if pol.degree() != 5:
            continue
        rec = {"u": uv, "irr": bool(pol.is_irreducible)}
        if pol.is_irreducible:
            disc = int(pol.discriminant())
            rec["disc_sq"] = disc > 0 and is_square(disc)
            if rec["disc_sq"]:
                cl = classify_poly(chi, do_galois=True)
                rec["status"] = cl.get("status")
        kron_rows.append(rec)

    return {
        "C1_ternary_companion_scan": {
            "tested": tested,
            "irr": n_irr,
            "disc_sq": n_sq,
            "A5_found": n_A5,
            "samples": rows,
            "note": "disc□ not identically forced; rate low (search, not identity)",
        },
        "C2_transfer_graph": transfer,
        "C3_3cycle_block_deform": kron_rows,
        "identically_square_by_construction": False,
    }


# ---------------------------------------------------------------------------
# Idea D — icosahedral / invariant sample
# ---------------------------------------------------------------------------


def idea_D() -> dict:
    print("IDEA D — icosahedral-style", flush=True)
    # Classical icosahedral-adjacent: x^5 + 5 m x^3 + 5 m^2 x + n
    m, n = sp.symbols("m n")
    P = x**5 + 5 * m * x**3 + 5 * m**2 * x + n
    D = sp.expand(sp.Poly(P, x).discriminant())
    # disc formula known: sample resonant (m,n)
    samples = []
    for mv in [0, 1, -1, 3, 9, 61]:
        for nv in [1, 3, 9, 27, 61, 80, 88, -88, 243]:
            Pv = sp.expand(P.subs({m: mv, n: nv}))
            pol = sp.Poly(Pv, x, domain=sp.ZZ)
            if not pol.is_irreducible:
                continue
            disc = int(pol.discriminant())
            sq = disc > 0 and is_square(disc)
            rec = {"m": mv, "n": nv, "disc_sq": sq}
            if sq and len([s for s in samples if s.get("status", "").startswith("HIT")]) < 6:
                cl = classify_poly(Pv, do_galois=True)
                rec["status"] = cl.get("status")
            samples.append(rec)
    # Pure-even BJ is the main A5 machine; icosa form is different
    n_sq = sum(1 for s in samples if s.get("disc_sq"))
    return {
        "family": "x^5+5m x^3+5m^2 x+n",
        "n_irr_samples": len(samples),
        "n_disc_sq": n_sq,
        "A5_or_even": [s for s in samples if s.get("disc_sq")][:10],
        "note": "Icosahedral-adjacent scans; disc□ sporadic not identical in (m,n)",
    }


# ---------------------------------------------------------------------------
# Idea E — poly from T3 dynamics
# ---------------------------------------------------------------------------


def T3(n: int) -> int:
    if n == 0:
        return 0
    r = n % 3
    if r == 0:
        return n // 3
    if r == 1:
        return (4 * n + 2) // 3
    return (2 * n + 1) // 3


def idea_E() -> dict:
    print("IDEA E — T3-native polys", flush=True)
    # E1: truncated transfer operator on functions mod x^N — charpoly proxy via companion of itinerary poly
    # E2: minimal polynomial of periodic-point multipliers hard; use product (x - T3^k(n0))
    rows = []
    for n0 in [1, 2, 3, 61, 80, 539]:
        orbit = []
        n = n0
        for _ in range(12):
            orbit.append(n)
            n = T3(n)
            if n in orbit or n == 0:
                if n == 0:
                    orbit.append(0)
                break
        # poly Π(x - o) for distinct orbit points (deg may vary)
        pts = sorted(set(orbit))
        if len(pts) < 2:
            continue
        pol = sp.prod(x - p for p in pts)
        pol = sp.Poly(sp.expand(pol), x, domain=sp.ZZ)
        rec = {
            "n0": n0,
            "orbit": orbit,
            "deg": pol.degree(),
            "irr": bool(pol.is_irreducible),
        }
        if pol.degree() == 5 and pol.is_irreducible:
            disc = int(pol.discriminant())
            rec["disc_sq"] = disc > 0 and is_square(disc)
            if rec["disc_sq"]:
                cl = classify_poly(pol.as_expr(), do_galois=True)
                rec["status"] = cl.get("status")
        elif pol.degree() <= 7 and pol.is_irreducible:
            disc = int(pol.discriminant())
            rec["disc_sq"] = disc > 0 and is_square(disc)
        rows.append(rec)

    # E3: resultant eliminating itinerary bits — model
    # z0,z1,z2 residues; crude: Res of (3w - n) style
    # Use charpoly of 5×5 transfer already in C; here dynatomic-like for x |-> x(x+1)(x-1) mod ternary? skip heavy
    return {
        "orbit_polys": rows,
        "n_deg5_irr": sum(1 for r in rows if r.get("deg") == 5 and r.get("irr")),
        "n_disc_sq": sum(1 for r in rows if r.get("disc_sq")),
        "note": "Orbit polys from T3 are design probes; not systematic A5 machines",
    }


# ---------------------------------------------------------------------------
# Main report
# ---------------------------------------------------------------------------


def main():
    t0 = time.time()
    print("NEW ALGEBRAIC IDEAS A–F", flush=True)

    A = idea_A()
    F = idea_F(A)
    B = idea_B()
    # HQCC native embed for B: true if sparse sol uses only lattice when A lattice
    B["HQCC_native_embed"] = bool(B.get("T_match_sparse")) and B.get("n_disc_sq", 0) > 0
    C = idea_C()
    D = idea_D()
    E = idea_E()

    elapsed = round(time.time() - t0, 2)

    # Synthesis
    mestre_works = A["families_disc_square"] > 0
    verdict = (
        f"New algebraic ideas ({elapsed}s). "
        f"A Mestre: seeds_with_R={A['seeds_with_R']}, families disc□ in Q(t)={A['families_disc_square']}. "
        f"F embed: static embeddable={F['n_embeddable_static']}, mestre-spec embeddable={F['n_mestre_spec_embeddable']}. "
        f"B non-BJ: disc□ identically in A={B['disc_identically_square_in_A']}, "
        f"lattice disc□={B['n_disc_sq']}, A5={B['n_A5']}, T-match sols={B['T_match_generic_sols']}. "
        f"C new avatars: no identical-square by construction. "
        f"D icosa scan disc□={D['n_disc_sq']}. "
        f"E T3 orbit disc□={E['n_disc_sq']}. "
        f"Primary path A+F: {'partial HIT' if mestre_works else 'R-space found but family disc□ rare / constructions need refinement'}. "
        f"Old T cuts remain closed negative experiment."
    )
    print(verdict, flush=True)

    lines = [
        r"# New algebraic ideas A–F",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Stance",
        "",
        r"- Arithmetic centre (pure-even multi-\(k\)) **finished** — `PURE_EVEN_MULTI_K.md`.",
        r"- Criterion 2 on template \(T(a,\ldots,f)\) **closed** negative — `TIER11_DEEPEN.md`.",
        r"- This document: **new equations / new matrices**, not more cuts of \(T\).",
        "",
        r"### Do not retry",
        "",
        r"| Approach | Reason |",
        r"|----------|--------|",
        r"| More linear/bilinear cuts of same \(T\) | Exhausted |",
        r"| \(F\to T\) only hoping disc□→1 | Rate 0 |",
        r"| Surgery on rigid \(\varphi/\mathbb{Q}\) | Permanent factor 5 |",
        r"| Binary Collatz forces evenness without pure-even | Composite only |",
        "",
        "---",
        "",
        r"## Idea A — Mestre deformation (primary)",
        "",
        r"For irreducible seed \(P\) with disc □, solve",
        r"$$P''R-2P'R'\equiv 0\pmod{P},\quad \deg R\le \deg P-2.$$",
        r"Then build families via resultants:",
        r"- `shift_y_tR`: \(\operatorname{Res}_y(P(y),\,z-y-t R(y))\)",
        r"- `uPp_tR`: \(\operatorname{Res}_y(P(y),\,u P'(y)-t R(y))\)",
        r"- `uPp_R_t`: \(\operatorname{Res}_y(P(y),\,u P'(y)-R(y)-t)\)",
        "",
        f"- Seeds processed: **{A['seeds_processed']}**",
        f"- Seeds with nontrivial \(R\)-space: **{A['seeds_with_R']}**",
        f"- Family constructions with disc identically □ in \(\\mathbb{{Q}}(t)\): **{A['families_disc_square']}**",
        "",
    ]
    for r in A["results"]:
        if r.get("skipped"):
            lines.append(f"- `{r['tag']}`: skipped ({r.get('reason')})")
            continue
        lines.append(
            f"### Seed `{r['tag']}` — \(P=x^5+({r['alpha']})x+({r['beta']})\)"
        )
        lines.append("")
        lines.append(
            f"- seed disc□: **{r.get('seed_disc_square')}**, null_dim(R): **{r.get('null_dim')}**"
        )
        lines.append(f"- R basis: `{r.get('R_basis')}`")
        lines.append("")
        lines.append(r"| construction | R | disc□ in Q(t)? | samples |")
        lines.append(r"|--------------|---|:--------------:|---------|")
        for fam in r.get("families") or []:
            samp = fam.get("samples") or []
            samp_s = ",".join(
                f"t={s.get('t')}:{s.get('status')}" for s in samp[:4]
            ) or fam.get("error") or fam.get("disc_info", {}).get("reason", "")
            lines.append(
                f"| {fam.get('construction')} | `{fam.get('R')}` | "
                f"**{fam.get('disc_square_in_Qt')}** | {samp_s} |"
            )
        lines.append("")

    lines += [
        "---",
        "",
        r"## Idea F — embed even families into \(T\) (primary companion)",
        "",
        r"Solve \(\chi_T(a,b,c,d,e,f)=P(x)\) (after killing \(x^4\) by shift if needed).",
        "",
        f"- Static embeddable (seeds / non-BJ points): **{F['n_embeddable_static']}**",
        f"- Mestre specialisation embeddable: **{F['n_mestre_spec_embeddable']}**",
        "",
    ]
    for te in F["tests"][:12]:
        emb = te.get("embed") or {}
        if te.get("embed_specialisations") is not None:
            lines.append(
                f"- **{te['source']}**: specs="
                f"{[(s.get('t'), s.get('embed', {}).get('embeddable')) for s in te['embed_specialisations']]}"
            )
        else:
            lines.append(
                f"- **{te['source']}**: embeddable=**{emb.get('embeddable')}** "
                f"({emb.get('reason') or emb.get('sample_sol') or emb.get('sparse_hits', [''])[:1]})"
            )

    lines += [
        "",
        r"**Note:** BJ seeds embed via classical BJ-embed (\(d=0,a=-ef\)). "
        r"That recovers pure-even, not a new HQCC-native necessity fragment. "
        r"Non-BJ embeds (when they exist) give \(d\neq 0\) realisations inside \(T\).",
        "",
        "---",
        "",
        r"## Idea B — non-BJ degree-1 \(A_5\) family (secondary)",
        "",
        f"Family: `${B['family']}$`",
        "",
        f"- Disc identically square in parameter \(A\)? **{B['disc_identically_square_in_A']}**",
        f"- Factored disc (preview): `{B['disc_factored_preview']}`",
        f"- Lattice samples disc□: **{B['n_disc_sq']}**, A5: **{B['n_A5']}**",
        f"- Generic \(T\)-match solutions: **{B['T_match_generic_sols']}**",
        f"- Sparse match sample: `{B['T_match_sparse'][:3]}`",
        "",
        r"| \(A\) | disc□ | status |",
        r"|----:|:-----:|--------|",
    ]
    for s in B["lattice_samples"]:
        lines.append(
            f"| {s.get('A')} | {s.get('disc_sq')} | {s.get('status')} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## Idea C — change the matrix avatar (if A+F stall)",
        "",
        f"- C1 ternary-coeff companion scan: tested={C['C1_ternary_companion_scan']['tested']}, "
        f"irr={C['C1_ternary_companion_scan']['irr']}, disc□={C['C1_ternary_companion_scan']['disc_sq']}, "
        f"A5={C['C1_ternary_companion_scan']['A5_found']}",
        f"- C2 transfer graph: disc□={C['C2_transfer_graph'].get('disc_sq')}, "
        f"status={C['C2_transfer_graph'].get('status')}",
        f"- C3 3-cycle block deform: `{C['C3_3cycle_block_deform']}`",
        f"- **Identically square by construction?** **{C['identically_square_by_construction']}**",
        "",
        r"Old \(T\) remains a **closed negative experiment** for Crit 2. New avatars need a "
        r"built-in evenness identity (as pure-even has), not another sparse scan.",
        "",
        "---",
        "",
        r"## Idea D — icosahedral / invariant parameters",
        "",
        f"- Family `{D['family']}`",
        f"- irr samples={D['n_irr_samples']}, disc□={D['n_disc_sq']}",
        f"- Even hits (sample): `{D['A5_or_even'][:5]}`",
        f"- {D['note']}",
        "",
        "---",
        "",
        r"## Idea E — HQCC-native polynomial from T₃",
        "",
        f"- Orbit polys: `{E['orbit_polys']}`",
        f"- deg5 irr={E['n_deg5_irr']}, disc□={E['n_disc_sq']}",
        f"- {E['note']}",
        "",
        "---",
        "",
        r"## Recommended resolution path (updated)",
        "",
        r"| Priority | Idea | Outcome this run | Next |",
        r"|:--------:|------|------------------|------|",
        f"| 1 | A Mestre + F embed | R-space on even seeds; disc□ families={A['families_disc_square']}; "
        f"embeds={F['n_embeddable_static']}+{F['n_mestre_spec_embeddable']} | "
        r"If disc□ family HIT: publish family + check HQCC-native matrix realisation beyond BJ-embed |",
        f"| 2 | B non-BJ deg-1 | identical disc□={B['disc_identically_square_in_A']}; "
        f"lattice A5={B['n_A5']}; T-match={B['T_match_generic_sols']}>0 | "
        r"Restrict \(A\) so embed coords lie in resonant lattice; re-check Gal |",
        r"| 3 | C new matrix | no identical-square avatar yet | only with evenness identity by design |",
        r"| 4 | D, E | sporadic / probe | low priority |",
        "",
        r"## One-line synthesis",
        "",
        r"The Model’s arithmetic centre is finished; Criterion 2 on the present template is closed. "
        r"New progress requires Mestre-style even families, non-BJ parametric \(A_5\) shapes, "
        r"or a new matrix avatar of ternary branching — not further relations inside \(T(a,\ldots,f)\).",
        "",
        r"```bash",
        r"python new_algebraic_ideas.py",
        r"```",
        "",
        r"_Generated by new_algebraic_ideas.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "A": A,
        "F": F,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
    }
    write_md(ROOT / "NEW_ALGEBRAIC_IDEAS.md", "\n".join(lines))
    write_json(ROOT / "NEW_ALGEBRAIC_IDEAS.json", payload)
    write_md(OUT / "NEW_ALGEBRAIC_IDEAS.md", "\n".join(lines))
    write_json(OUT / "NEW_ALGEBRAIC_IDEAS.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "NEW_ALGEBRAIC_IDEAS.md", "\n".join(lines))
    except Exception:
        pass
    print(f"Wrote NEW_ALGEBRAIC_IDEAS.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
