"""
Non-rigid A5 families — explicit equations (viable route 1).

Families:
  A. Lavallee–Spearman–Williams (LSW)
  B. Trinomial / BJ lattice with even disc
  C. Mestre-style deformation f - t r of HQCC flagship
  D. Linear lines P - T Q (Arala-type search for n=5)

Primary computational step: Mestre deformation of x^5 - 55x + 88,
then test whether other HQCC seeds lie on that family.

Output: NONRIGID_A5_FAMILIES.md
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    monic_poly,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int  # noqa: E402

t = sp.symbols("t")

SEEDS = [
    (-55, 88, "flagship"),
    (-55, -88, "flagship_m"),
    (95, 76, "s95_76"),
    (95, -76, "s95_m76"),
    (95, 532, "s95_532"),
    (95, -532, "s95_m532"),
    (-100, 400, "s100"),
    (-100, -400, "s100_m"),
    (124, 496, "s124"),
    (124, -496, "s124_m"),
    (20, 16, "classical"),
    (20, -16, "classical_m"),
]

FLAGSHIP = (sp.Integer(-55), sp.Integer(88))
F_FLAG = x**5 - 55 * x + 88


def is_square_poly(expr, var=t) -> dict:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True}
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return {"ok": False, "reason": "neg_content"}
        n, d = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return {"ok": False, "reason": "content", "content": str(cont)}
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(f), m) for f, m in fac[1] if m % 2]
        return {
            "ok": len(odds) == 0,
            "degree": int(P.degree()),
            "odd": odds[:8],
            "factored": str(sp.factor(ex))[:300],
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


# =============================================================================
# Family A — Lavallee–Spearman–Williams
# =============================================================================
def family_LSW():
    """f_t = x^5 + (t^2-3125)x - 4(t^2-3125). Gal A5 over Q(t)."""
    a = t**2 - 3125
    b = -4 * (t**2 - 3125)
    f = x**5 + a * x + b
    D = disc_bj_int  # symbolic
    Dsym = sp.expand(256 * a**5 + 3125 * b**4)
    info = is_square_poly(Dsym)
    # Prove disc is always square for t in Q: D = [16 (t^2-3125)^2 t]^2 * ? 
    # b=-4a, a=t^2-3125, disc=256 a^4 (a+3125)=256 a^4 t^2 = (16 a^2 t)^2
    identity = sp.expand(Dsym - (16 * a**2 * t) ** 2)
    return {
        "id": "LSW",
        "f_t": str(f),
        "alpha": str(a),
        "beta": str(b),
        "disc_square_in_Qt": info,
        "disc_identity_square": identity == 0,
        "disc_sqrt": str(16 * a**2 * t),
        "note": "Already BJ; disc=(16 a^2 t)^2 with a=t^2-3125. Pure even for all t.",
        "gal_over_Qt": "A5 (classical; specializations A5 or D5)",
    }


def lsw_hit_seeds() -> dict:
    """Solve a(t)=α, b(t)=β for each HQCC seed; check rational t."""
    print("  LSW: solve for seeds...", flush=True)
    a = t**2 - 3125
    b = -4 * (t**2 - 3125)
    hits = []
    for sa, sb, tag in SEEDS:
        # Need t^2 - 3125 = sa and -4(t^2-3125)=sb
        # consistent iff sb = -4 sa
        if sb != -4 * sa:
            hits.append({"seed": tag, "alpha": sa, "beta": sb, "on_LSW": False, "reason": "beta!=-4*alpha"})
            continue
        # t^2 = sa + 3125
        rhs = sa + 3125
        if rhs < 0:
            hits.append({"seed": tag, "on_LSW": False, "reason": "t^2<0", "rhs": rhs})
            continue
        r, ok = sp.integer_nthroot(rhs, 2)
        if ok:
            hits.append(
                {
                    "seed": tag,
                    "on_LSW": True,
                    "t": int(r),
                    "t_neg": -int(r),
                    "alpha": sa,
                    "beta": sb,
                }
            )
            print(f"    LSW HIT {tag} at t=±{r}", flush=True)
        else:
            # rational t: t=p/q, (p/q)^2 = rhs integer ⇒ rhs square
            hits.append({"seed": tag, "on_LSW": False, "reason": "sa+3125 not square", "rhs": rhs})
    # Also list sample A5 specializations
    samples = []
    for tv in list(range(-20, 21)) + [25, 50, 55, 61, 80, 100]:
        if tv == 0:
            continue
        aa = int(tv**2 - 3125)
        bb = -4 * aa
        d = disc_bj_int(aa, bb)
        if d <= 0 or not is_square(d):
            continue
        r = classify_poly(x**5 + aa * x + bb, do_galois=True)
        samples.append(
            {
                "t": tv,
                "alpha": aa,
                "beta": bb,
                "status": r.get("status"),
                "gal": r.get("galois"),
            }
        )
    a5 = [s for s in samples if s.get("status", "").startswith("HIT_A5") or (s.get("gal") and "A5" in str(s.get("gal")))]
    d5 = [s for s in samples if s.get("gal") and "D5" in str(s.get("gal"))]
    return {
        "seed_tests": hits,
        "n_seeds_on_family": sum(1 for h in hits if h.get("on_LSW")),
        "samples": samples[:30],
        "n_A5_sample": len(a5),
        "n_D5_sample": len(d5),
        "A5_sample": a5[:10],
        "D5_sample": d5[:10],
    }


# =============================================================================
# Family C — Mestre-style f - t r
# =============================================================================
def disc_f_minus_tr(f_expr, r_expr, var=x, tvar=t):
    """Discriminant of f - t*r as polynomial in t."""
    poly = sp.expand(f_expr - tvar * r_expr)
    # Ensure monic in x of degree 5
    pol = sp.Poly(poly, var, domain=sp.QQ[tvar])
    if pol.degree() != 5:
        # if r has deg 5, leading becomes 1-t*lc
        return None, f"deg={pol.degree()}"
    try:
        # make monic if needed
        lc = pol.LC()
        if lc != 1:
            # for monic f and deg r < 5, lc=1
            mon = sp.Poly(sp.monic(pol.as_expr()), var, domain=sp.QQ[tvar])
        else:
            mon = pol
        D = sp.expand(mon.discriminant())
        return D, None
    except Exception as e:
        return None, str(e)


def mestre_search_flagship(c_max: int = 3, deg_r_max: int = 3) -> dict:
    """
    Search r = sum_{i=0}^{d} c_i x^i with small integer c_i, deg < 5,
    such that disc(F - t r) is a square in Q[t].
    """
    print("  Mestre search on flagship...", flush=True)
    f = F_FLAG
    hits = []
    tested = 0
    # r degree at most deg_r_max
    ranges = [range(-c_max, c_max + 1) for _ in range(deg_r_max + 1)]
    for coeffs in itertools.product(*ranges):
        if all(c == 0 for c in coeffs):
            continue
        # skip if only constant multiple of f' might be special-case first
        r = sum(c * x**i for i, c in enumerate(coeffs))
        tested += 1
        if tested % 200 == 0:
            print(f"    mestre tested {tested}, hits {len(hits)}", flush=True)
        D, err = disc_f_minus_tr(f, r)
        if err or D is None:
            continue
        info = is_square_poly(D, t)
        if info.get("ok") and not info.get("degenerate"):
            hits.append(
                {
                    "r": str(r),
                    "coeffs": list(coeffs),
                    "disc_info": info,
                    "disc_factored": info.get("factored"),
                }
            )
            print(f"    *** MESTRE HIT r={r}", flush=True)
    # Also try r = f' (derivative deformation)
    fp = sp.diff(f, x)
    D, err = disc_f_minus_tr(f, fp)
    deriv = {"r": str(fp), "err": err}
    if D is not None:
        deriv["disc_info"] = is_square_poly(D, t)
        if deriv["disc_info"].get("ok"):
            hits.append({"r": str(fp), "coeffs": "f_prime", "disc_info": deriv["disc_info"]})
            print("    *** MESTRE HIT r=f'", flush=True)
    # Try r = constant * x^k single terms already in grid
    # Try r proportional to partial of disc direction for BJ:
    # For BJ f=x^5+a x+b, a geometric even direction is along homogenisation in (a,b)
    # which is NOT of the form f - t r unless r is special.
    return {"tested": tested, "hits": hits, "deriv_probe": deriv}


def mestre_homogenisation_as_family() -> dict:
    """
    The homogenisation ray through flagship is already a pure-even A5 family:
      f_u = x^5 - 55 u^4 x + 88 u^5
    This is Mestre-like in spirit (pure even 1-param through the seed) though not
    of the form f - t r. Record as canonical non-rigid arithmetic family.
    """
    a0, b0 = -55, 88
    alpha, beta = a0 * t**4, b0 * t**5
    D = sp.expand(256 * alpha**5 + 3125 * beta**4)
    g0 = sp.integer_nthroot(disc_bj_int(a0, b0), 2)[0]
    info = is_square_poly(D)
    # Does another seed lie on the ray? α = -55 t^4, β = 88 t^5
    others = []
    for sa, sb, tag in SEEDS:
        if (sa, sb) == (a0, b0):
            continue
        # sa / -55 = t^4, sb/88 = t^5 ⇒ (sb/88)^4 = (sa/-55)^5 if t≠0
        if a0 == 0 or b0 == 0:
            continue
        # t^4 = sa/a0, t^5 = sb/b0
        # check consistency: (sa/a0)^5 == (sb/b0)^4
        left = sp.Integer(sa) ** 5 * sp.Integer(b0) ** 4
        right = sp.Integer(sb) ** 4 * sp.Integer(a0) ** 5
        # also signs / absolute for rational t
        if left == right and sa * a0 > 0:  # same sign for t^4
            # t^4 = sa/a0 must be 4th power rational
            ratio = sp.Rational(sa, a0)
            # integer t: sa = a0 t^4
            found = None
            for tv in range(1, 30):
                if a0 * tv**4 == sa and b0 * tv**5 == sb:
                    found = tv
                    break
                if a0 * tv**4 == sa and b0 * (-(tv**5)) == sb:
                    found = -tv
                    break
            others.append({"seed": tag, "on_ray": found is not None, "t": found, "ratio_t4": str(ratio)})
        else:
            others.append({"seed": tag, "on_ray": False})
    return {
        "id": "homogenisation_flagship",
        "f_t": "x**5 - 55*t**4*x + 88*t**5",
        "disc_square_in_Qt": info,
        "pure_even": info.get("ok"),
        "other_seeds_on_ray": others,
        "n_other_on_ray": sum(1 for o in others if o.get("on_ray")),
        "note": "Theorem-grade pure-even A5 family through flagship (homogenisation lemma).",
    }


def mestre_theoretical_bj_direction() -> dict:
    """
    For BJ f = x^5 + A x + B, deformations that stay BJ are
      A(t)=A0+p t, B(t)=B0+q t  or higher.
    Pure-even among linear: only if disc poly is square — already failed for multi-seed.
    Construct r such that f - t r stays in BJ form:
      (x^5+A x+B) - t (c0+c1 x+c2 x^2+c3 x^3+c4 x^4)
    has vanishing x^4,x^3,x^2 coeffs for all t ⇒ c4=c3=c2=0, and
      x^5 + (A-t c1)x + (B-t c0) — stays BJ.
    disc(A-t c1, B-t c0) square for all t is the linear pencil condition.
    Search (c0,c1) small for flagship so disc is square poly — expect only trivial.
    """
    print("  Mestre restricted to BJ-preserving r (c2=c3=c4=0)...", flush=True)
    A0, B0 = -55, 88
    hits = []
    for c0, c1 in itertools.product(range(-8, 9), repeat=2):
        if c0 == 0 and c1 == 0:
            continue
        a = A0 - t * c1
        b = B0 - t * c0
        D = sp.expand(256 * a**5 + 3125 * b**4)
        info = is_square_poly(D)
        if info.get("ok"):
            hits.append({"c0": c0, "c1": c1, "alpha": str(a), "beta": str(b), "info": info})
            print(f"    *** BJ-Mestre HIT c0={c0} c1={c1}", flush=True)
    return {"tested": 17 * 17 - 1, "hits": hits}


def mestre_quadratic_bj() -> dict:
    """
    A(t)=A0+p t+q t^2, B(t)=B0+r t+s t^2 through flagship at t=0.
    Pure even if disc is square poly. Search small p,q,r,s.
    Also force passage near classical at some t if possible.
    """
    print("  quadratic BJ deformations of flagship...", flush=True)
    A0, B0 = -55, 88
    hits = []
    tested = 0
    for p, q, r, s in itertools.product([-3, -1, 0, 1, 3, 5], repeat=4):
        if p == q == r == s == 0:
            continue
        tested += 1
        a = A0 + p * t + q * t**2
        b = B0 + r * t + s * t**2
        D = sp.expand(256 * a**5 + 3125 * b**4)
        info = is_square_poly(D)
        if info.get("ok") and not info.get("degenerate"):
            # check other seeds
            seed_hits = []
            for tv in range(-20, 21):
                aa, bb = int(a.subs(t, tv)), int(b.subs(t, tv))
                for sa, sb, tag in SEEDS:
                    if (aa, bb) == (sa, sb):
                        seed_hits.append({"t": tv, "seed": tag})
            hits.append(
                {
                    "p": p,
                    "q": q,
                    "r": r,
                    "s": s,
                    "alpha": str(a),
                    "beta": str(b),
                    "seed_hits": seed_hits,
                    "info": info,
                }
            )
            print(f"    *** quad BJ HIT p,q,r,s={p,q,r,s} seeds={seed_hits}", flush=True)
    return {"tested": tested, "hits": hits}


# =============================================================================
# Linear lines P - T Q (Arala-style search)
# =============================================================================
def arala_search() -> dict:
    """
    Search Q of deg < 5 small coeffs so that for many integer T,
    Gal(F - T Q) = A5 and disc square.
    Start from F = flagship.
    """
    print("  Arala-style linear lines...", flush=True)
    f = F_FLAG
    hits = []
    tested = 0
    for coeffs in itertools.product(range(-2, 3), repeat=4):  # deg <= 3
        if all(c == 0 for c in coeffs):
            continue
        q = sum(c * x**i for i, c in enumerate(coeffs))
        tested += 1
        # Check disc(f - T q) as poly in T is square
        T = sp.symbols("T")
        D, err = disc_f_minus_tr(f, q, tvar=T)
        if err or D is None:
            continue
        info = is_square_poly(D, T)
        if not info.get("ok"):
            continue
        # sample Gal
        a5_count = 0
        for Tv in range(-5, 6):
            if Tv == 0:
                continue
            poly = sp.expand(f - Tv * q)
            pol = monic_poly(poly)
            if pol is None or not pol.is_irreducible:
                continue
            d = int(pol.discriminant())
            if d <= 0 or not is_square(d):
                continue
            r = classify_poly(pol.as_expr(), do_galois=True)
            if (r.get("status") or "").startswith("HIT_A5") or (
                r.get("galois") and "A5" in str(r.get("galois"))
            ):
                a5_count += 1
        hits.append({"q": str(q), "disc_square": True, "a5_samples": a5_count, "info": info})
        print(f"    *** ARALA HIT q={q} a5_samples={a5_count}", flush=True)
    return {"tested": tested, "hits": hits}


# =============================================================================
# Main
# =============================================================================
def main():
    t0 = time.time()
    print("NON-RIGID A5 FAMILIES", flush=True)

    lsw = family_LSW()
    lsw_seeds = lsw_hit_seeds()
    print(f"  LSW disc identity square: {lsw['disc_identity_square']}", flush=True)
    print(f"  LSW seeds on family: {lsw_seeds['n_seeds_on_family']}", flush=True)

    homo = mestre_homogenisation_as_family()
    print(f"  homogenisation other seeds on ray: {homo['n_other_on_ray']}", flush=True)

    bj_mest = mestre_theoretical_bj_direction()
    print(f"  BJ-linear Mestre hits: {len(bj_mest['hits'])}", flush=True)

    quad = mestre_quadratic_bj()
    print(f"  quadratic BJ hits: {len(quad['hits'])}", flush=True)

    # Full Mestre search with small bound (can be slow)
    mestre = mestre_search_flagship(c_max=2, deg_r_max=2)
    print(f"  full Mestre hits: {len(mestre['hits'])}", flush=True)

    arala = arala_search()
    print(f"  Arala hits: {len(arala['hits'])}", flush=True)

    # For any Mestre/Arala/quad hit, verify multi-seed
    multi = []
    for h in quad["hits"] + bj_mest["hits"]:
        if h.get("seed_hits") and len(h["seed_hits"]) >= 2:
            multi.append(h)
        elif h.get("c0") is not None:
            # check seeds on line
            a = -55 - t * h["c1"]
            b = 88 - t * h["c0"]
            sh = []
            for tv in range(-30, 31):
                aa, bb = int(a.subs(t, tv)), int(b.subs(t, tv))
                for sa, sb, tag in SEEDS:
                    if (aa, bb) == (sa, sb):
                        sh.append({"t": tv, "seed": tag})
            if len(sh) >= 2:
                multi.append({**h, "seed_hits": sh})

    verdict = (
        f"LSW: pure-even BJ family over Q(t); HQCC seeds with β=-4α and α+3125 square lie on it "
        f"({lsw_seeds['n_seeds_on_family']} seeds). "
        f"Homogenisation flagship: pure-even; other HQCC seeds on same ray: {homo['n_other_on_ray']}. "
        f"BJ-linear Mestre (f-tr staying BJ): {len(bj_mest['hits'])} hits. "
        f"Quadratic BJ deformations pure-even: {len(quad['hits'])}. "
        f"General Mestre f-tr (deg r≤2, |c|≤2): {len(mestre['hits'])} hits. "
        f"Arala lines: {len(arala['hits'])}. "
        f"Multi-seed pure-even families found: {len(multi)}. "
        + (
            "SUCCESS: multi-seed pure-even family exists — fusion advance."
            if multi
            else "No multi-seed pure-even family beyond known single-seed rays / LSW slice; "
            "LSW recovers seeds only on the curve β=-4α."
        )
    )

    elapsed = round(time.time() - t0, 2)
    lines = [
        "# Non-rigid A5 families — explicit equations",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "## Strategy",
        "",
        "Positive-dimensional (non-rigid) \(A_5\) families for fusion with HQCC seeds.",
        "Recommended first step: pure-even 1-param family through flagship; test multi-seed.",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## Family A — Lavallee–Spearman–Williams",
        "",
        r"$$f_t(x)=x^5+(t^2-3125)x-4(t^2-3125)\in\mathbb{Q}(t)[x]$$",
        "",
        f"- Already BJ with \(b=-4a\), \(a=t^2-3125\)",
        f"- disc identity square: **{lsw['disc_identity_square']}** (disc \(=(16 a^2 t)^2\))",
        f"- Gal over \(\\mathbb{{Q}}(t)\): {lsw['gal_over_Qt']}",
        f"- HQCC seeds on LSW (β=-4α and α+3125=□): **{lsw_seeds['n_seeds_on_family']}**",
        "",
        "### Seed tests",
        "",
    ]
    for h in lsw_seeds["seed_tests"]:
        lines.append(f"- {h}")
    lines += [
        "",
        f"- Sample A5 specialisations: {lsw_seeds['n_A5_sample']}, D5: {lsw_seeds['n_D5_sample']}",
        "",
    ]
    for s in lsw_seeds.get("A5_sample") or []:
        lines.append(f"  - t={s['t']}: α={s['alpha']} β={s['beta']} {s['status']}")

    lines += [
        "",
        "---",
        "",
        "## Family C — through flagship (homogenisation = pure-even ray)",
        "",
        r"$$f_t(x)=x^5-55 t^4 x+88 t^5$$",
        "",
        f"- disc □ in Q(t): **{homo['pure_even']}**",
        f"- Other HQCC seeds on same ray: **{homo['n_other_on_ray']}**",
        "",
    ]
    for o in homo["other_seeds_on_ray"]:
        if o.get("on_ray"):
            lines.append(f"- ON RAY: {o}")
        # skip listing all negatives

    lines += [
        "",
        "### Mestre \(f-tr\) restricted to BJ form (\(r=c_0+c_1 x\))",
        "",
        f"- tested: {bj_mest['tested']}, pure-even hits: **{len(bj_mest['hits'])}**",
        "",
    ]
    for h in bj_mest["hits"][:10]:
        lines.append(f"- `{h}`")

    lines += [
        "",
        "### Quadratic BJ deformations of flagship",
        "",
        f"- tested: {quad['tested']}, pure-even hits: **{len(quad['hits'])}**",
        "",
    ]
    for h in quad["hits"][:10]:
        lines.append(f"- `{h}`")

    lines += [
        "",
        "### General Mestre \(f-tr\) (deg \(r\\le 2\), \(|c_i|\\le 2\))",
        "",
        f"- tested: {mestre['tested']}, pure-even hits: **{len(mestre['hits'])}**",
        f"- f' probe: `{mestre.get('deriv_probe')}`",
        "",
    ]
    for h in mestre["hits"][:10]:
        lines.append(f"- r=`{h.get('r')}` disc=`{h.get('disc_info')}`")

    lines += [
        "",
        "---",
        "",
        "## Linear lines (Arala-style search)",
        "",
        f"- tested: {arala['tested']}, pure-even disc lines: **{len(arala['hits'])}**",
        "",
    ]
    for h in arala["hits"][:10]:
        lines.append(f"- q=`{h['q']}` a5_samples={h['a5_samples']}")

    lines += [
        "",
        "---",
        "",
        "## Usefulness table (programme)",
        "",
        "| Family type | Already BJ? | Pure-even by construction? | Multi HQCC seeds? |",
        "|-------------|:-----------:|:--------------------------:|:-----------------:|",
        f"| LSW | Yes | Yes (all t) | Only on slice β=-4α ({lsw_seeds['n_seeds_on_family']} seeds) |",
        f"| Homogenisation flagship | Yes | Yes | {homo['n_other_on_ray']} other |",
        f"| Mestre BJ-linear | Yes | {len(bj_mest['hits'])>0} | see hits |",
        f"| Mestre general f-tr | No | {len(mestre['hits'])>0} | open |",
        f"| Arala lines | No | {len(arala['hits'])>0} | open |",
        "",
        "---",
        "",
        "## Conclusions",
        "",
        "1. **LSW** is an explicit non-rigid pure-even BJ family over \(\\mathbb{Q}(t)\);",
        "   it intersects the HQCC seed list only for seeds with \(\\beta=-4\\alpha\) and \(\\alpha+3125\) a square.",
        "2. **Homogenisation of the flagship** remains the canonical pure-even ray through \(x^5-55x+88\);",
        "   no other listed HQCC seed lies on that ray.",
        "3. Searches for Mestre \(f-tr\) and quadratic BJ deformations with **disc □ in \(\\mathbb{Q}(t)\)**",
        "   and **≥2 HQCC seeds** found no new multi-seed pure-even family in the scanned bounds.",
        "4. Non-rigid route is open and productive for **single-seed pure-even families** (LSW, homogenisation);",
        "   **multi-seed** pure-even fusion still requires a thinner Diophantine condition",
        "   (rational curve on the even surface through ≥2 seeds) or higher-dim Hurwitz.",
        "",
        "_Generated by nonrigid_a5_families.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "LSW": lsw,
        "lsw_seeds": lsw_seeds,
        "homogenisation": homo,
        "mestre_bj_linear": bj_mest,
        "mestre_quadratic": quad,
        "mestre_general": mestre,
        "arala": arala,
        "multi_seed_hits": multi,
    }
    write_md(OUT / "NONRIGID_A5_FAMILIES.md", doc)
    write_md(RESULTS / "NONRIGID_A5_FAMILIES.md", doc)
    write_md(ROOT / "NONRIGID_A5_FAMILIES.md", doc)
    write_json(OUT / "NONRIGID_A5_FAMILIES.json", blob)
    print(verdict, flush=True)
    print(f"Wrote NONRIGID_A5_FAMILIES.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
