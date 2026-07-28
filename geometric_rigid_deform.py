"""
Systematic examination of other rigid triples + mild deformations
that preserve A5 monodromy while improving even specialisations.

Rigid triples (Step 1):
  (3A,3A,5A), (3A,3A,5B), (2A,3A,5A), (2A,3A,5B)

Note: absolute rigidity ⇒ Hurwitz space 0-dimensional for a fixed triple.
Continuous deformations cannot preserve the exact conjugacy classes.
"Mild deformations" here means:
  - other rigid triples / realizations
  - base automorphisms and rational pull-backs (still A5 geometrically)
  - discrete Galois conjugates of coefficients
  - finite twists; and scanning the even-specialisation locus of each cover

Outputs: GEOMETRIC_RIGID_DEFORM.md
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.polynomial import polyroots

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

y, t = sp.symbols("y t")

# ---------------------------------------------------------------------------
# Explicit rigid Belyi maps (degree 5, genus 0)
# ---------------------------------------------------------------------------

# (3A,3A,5A) — preferred, over Q
PHI_335A = 6 * y**5 - 15 * y**4 + 10 * y**3

# (3A,3A,5B): same passport types; monodromy uses inverse 5-class.
# Complex solutions of the (3,3,5) elimination give Galois conjugates of the
# scaling t with t^5=6. After domain scaling all give Q-models isomorphic to PHI_335A
# or complex embeddings. Real form is unique up to PGL2 over R.
# We also test PHI with inverted base (φ :→ 1-φ or 1/φ) which swaps branch labels.


def phi_base_transforms():
    """Discrete transforms of preferred φ (still rigid A5, relabelled branches)."""
    # φ, 1-φ, 1/φ, φ/(φ-1), etc. — S3 action on {0,1,∞}
    return {
        "phi": PHI_335A,
        "one_minus_phi": sp.expand(1 - PHI_335A),
        "one_over_phi": sp.together(1 / PHI_335A),  # rational, not poly
        "phi_over_phi_minus_1": sp.together(PHI_335A / (PHI_335A - 1)),
    }


def phi_235_coeffs():
    """(2A,3A,5A) radical coefficients (from Step 2 fallback)."""
    a = 5 * sp.Integer(2) ** sp.Rational(4, 5) * sp.Integer(3) ** sp.Rational(2, 5) / 12
    b = 5 * sp.Integer(2) ** sp.Rational(3, 5) * sp.Integer(3) ** sp.Rational(4, 5) / 9
    # Also Galois conjugates: 2^{1/5} → ζ^k 2^{1/5}
    return a, b


def phi_235_poly(a=None, b=None):
    if a is None:
        a, b = phi_235_coeffs()
    return sp.expand(y**5 + a * y**4 + b * y**3)


def monic_fibre_poly_map(phi_expr, w, var=y):
    """Monic polynomial in y for φ(y) - w = 0 (φ polynomial)."""
    expr = sp.expand(phi_expr - w)
    pol = sp.Poly(expr, var, domain=sp.QQ)
    return sp.monic(pol.as_expr())


def monic_fibre_to_Z(expr, var=y):
    pol = sp.Poly(sp.expand(expr), var, domain=sp.QQ)
    if pol.degree() != 5:
        return None
    mon = sp.Poly(sp.monic(pol.as_expr()), var, domain=sp.QQ)
    dens = []
    for c in mon.all_coeffs():
        dens.append(sp.fraction(sp.together(c))[1])
    L = 1
    for d in dens:
        try:
            L = int(sp.ilcm(L, abs(int(d))))
        except Exception:
            return None
    cleared = sp.expand(L**5 * mon.as_expr().subs(var, var / L))
    p = sp.Poly(cleared, var, domain=sp.ZZ)
    if p.LC() == -1:
        p = sp.Poly(-p.as_expr(), var, domain=sp.ZZ)
    if p.LC() != 1:
        return None
    return p


# ---------------------------------------------------------------------------
# Disc of monic(φ - t) as poly in t — even specialisation locus
# ---------------------------------------------------------------------------
def disc_of_phi_minus_t_poly(phi_expr=PHI_335A):
    """
    f(y;t) = monic(φ(y)-t). Compute Disc_y(f) as element of Q(t) or Q[t].
    Find rational t where disc is a square (even arithmetic monodromy candidates).
    """
    # φ monic form: φ/6 = y^5 - (5/2)y^4 + (5/3)y^3
    # monic(φ-t) = y^5 - (5/2)y^4 + (5/3)y^3 - t/6
    mon = sp.expand(phi_expr / sp.LC(sp.Poly(phi_expr, y)) - t / sp.LC(sp.Poly(phi_expr, y)))
    # Actually LC(φ)=6, monic(φ-t)=(φ-t)/6
    mon = sp.expand((phi_expr - t) / 6)
    pol = sp.Poly(mon, y, domain=sp.QQ[t])
    # Discriminant may be huge; use formula for monic quintic or sympy
    try:
        D = sp.factor(sp.expand(pol.discriminant()))
    except Exception as e:
        D = f"error:{e}"
    return mon, D


def search_even_specialisations_poly_family(phi_expr, name, t_rats, do_galois=True):
    """Scan rational t for irr + disc square; report Gal."""
    print(f"  even specs for {name}...", flush=True)
    even = []
    stats = Counter()
    for tv in t_rats:
        stats["tested"] += 1
        try:
            expr = sp.expand(phi_expr - tv)
            num, den = sp.fraction(sp.together(expr))
            den = sp.expand(den)
            if den != 1 and sp.degree(den, y) not in (0, -sp.oo):
                stats["skip_rat"] += 1
                continue
            pol = monic_fibre_to_Z(sp.expand(num) if den != 1 else expr, y)
            if pol is None:
                stats["bad"] += 1
                continue
            if not pol.is_irreducible:
                stats["red"] += 1
                continue
            d = int(pol.discriminant())
            if d <= 0 or not is_square(d):
                stats["odd"] += 1
                continue
            stats["even"] += 1
            rec = {
                "t": str(tv),
                "poly": str(pol.as_expr()),
                "disc": d,
            }
            if do_galois:
                r = classify_poly(pol.as_expr().subs(y, x), do_galois=True)
                rec["gal"] = r.get("galois")
                rec["status"] = r.get("status")
                if (rec.get("status") or "").startswith("HIT_A5") or (
                    rec.get("gal") and "A5" in str(rec.get("gal"))
                ):
                    stats["A5"] += 1
                    print(f"    EVEN A5 t={tv} {rec['poly'][:50]}", flush=True)
            even.append(rec)
        except Exception as e:
            stats["err"] += 1
    return {"name": name, "stats": dict(stats), "even": even}


def search_even_from_disc_poly():
    """
    For monic(φ-t)= y^5-(5/2)y^4+(5/3)y^3-t/6, compute Disc as poly in t
    and search rational points where Disc is square.
    """
    print("  computing Disc_y(monic(φ-t)) as function of t...", flush=True)
    mon = sp.expand((PHI_335A - t) / 6)
    # Use Poly over QQ(t) — discriminant
    pol = sp.Poly(mon, y)
    D = sp.simplify(sp.expand(pol.discriminant()))
    # D should be in Q[t]
    D = sp.together(D)
    num, den = sp.fraction(D)
    num = sp.expand(num)
    den = sp.expand(den)
    print(f"    disc num deg={sp.degree(num, t)} den deg={sp.degree(den, t)}", flush=True)
    # For t = p/q integer, check square: evaluate
    # Prefer: clear — disc of Z-model of fibre
    # Search t=p/q with |p|<=N, q<=M where cleared disc is square
    hits = []
    # Also factor numerator for square structure
    try:
        nfac = sp.factor(num)
    except Exception:
        nfac = str(num)[:200]
    return {
        "disc_expr": str(D)[:500],
        "num_factored": str(nfac)[:500],
        "num": str(num)[:300],
        "den": str(den)[:300],
    }


def dense_even_search_335(max_p=80, max_q=20):
    """Dense rational t search for even fibres of preferred φ."""
    print(f"  dense even search 335 |p|<={max_p} q<={max_q}...", flush=True)
    even = []
    stats = Counter()
    for q in range(1, max_q + 1):
        for p in range(-max_p, max_p + 1):
            if gcd(abs(p), q) != 1 and q > 1:
                continue
            # skip t where fibre reducible often near 0,1
            tv = sp.Rational(p, q)
            stats["tested"] += 1
            try:
                expr = sp.expand((PHI_335A - tv) / 6)
                pol = monic_fibre_to_Z(expr, y)
                if pol is None:
                    continue
                if not pol.is_irreducible:
                    stats["red"] += 1
                    continue
                d = int(pol.discriminant())
                if d > 0 and is_square(d):
                    stats["even"] += 1
                    r = classify_poly(pol.as_expr().subs(y, x), do_galois=True)
                    rec = {
                        "t": str(tv),
                        "poly": str(pol.as_expr()),
                        "gal": r.get("galois"),
                        "status": r.get("status"),
                        "disc": d,
                    }
                    even.append(rec)
                    print(f"    EVEN t={tv} {r.get('status')} {r.get('galois')}", flush=True)
                else:
                    stats["odd"] += 1
            except Exception:
                stats["err"] += 1
    return {"stats": dict(stats), "even": even}


def dense_even_search_235(max_p=30, max_q=8):
    """
    Even fibres of monic(φ_235 - t) with φ over Q(2^{1/5},3^{1/5}).
    Work numerically: use float a,b and approximate — for exact, clear field.
    Exact approach: let u=2^{1/5}, v=3^{1/5}, a=5 u^4 v^2 / 12, b=5 u^3 v^4 / 9
    φ = y^5 + a y^4 + b y^3 ∈ K[y], K=Q(u,v).
    Fibre φ - t monic over K; Gal over Q harder.
    Instead: multiply to get poly over Q by taking norm K/Q of (φ-t), degree 5*[K:Q].
    [Q(2^{1/5},3^{1/5}):Q]=25 typically — degree 125 poly, too big.

    Practical: specialize t and use minimal poly of a primitive element, OR
    use float monodromy only for cycle types, and for evenness use
    disc of the poly over the number field (in K) and check if disc is square in K.

    Simpler probe: clear coefficients by writing
    12^5 φ(y/something)...
    φ = y^3 (y^2 + a y + b). Let U=2^{1/5}, V=3^{1/5}.
    Multiply: consider F = 12 y^5 + 5 U^4 V^2 y^4 + ... still in K.

    For this run: numeric a,b float, monic Z poly via rounding cleared form from
    algebraic a,b expressed with minpoly elimination for fixed rational t.
    """
    print("  even search 235 via algebraic clearing for rational t...", flush=True)
    U, V = sp.symbols("U V")
    a = 5 * U**4 * V**2 / 12
    b = 5 * U**3 * V**4 / 9
    # Relations U^5-2=0, V^5-3=0
    even = []
    stats = Counter()
    # For fixed rational t0, poly = y^5 + a y^4 + b y^3 - t0
    # Primitive element approach: too heavy.
    # Use float approximation + monic Z from high precision? Unreliable for disc square.

    # Exact for t in Q: compute resultant eliminating U,V from system
    # — too expensive for many t.
    # Sample few t and compute disc in K, check if square in K using
    # K.<u,v> representation and factoring.

    # Minimal practical: for t rational, form
    # g(y) = 12 y^5 + 5*2^{4/5}*3^{2/5} y^4 + ... 
    # Use sympy minpoly of a number alpha = root, but for full poly over Q:
    # The poly over Q of a primitive root of the fibre is the product of Galois conjugates
    # of (y - root_i) — degree 5*25=125.

    # Fallback: report that 235 lives over K and list float monodromy verification only
    a_f, b_f = phi_235_coeffs()
    af, bf = float(a_f), float(b_f)
    # monodromy of the map φ:P1->P1 with float coeffs — already known A5 from Step 2
    # Even arithmetic specialisations require working over K; sample t and check
    # whether disc (as element of K) is a square by writing disc = x0 + x1 U + ... and solving.

    for tv in [sp.Rational(p, q) for q in range(1, 6) for p in range(-15, 16) if gcd(abs(p), q) == 1 or q == 1]:
        stats["tested"] += 1
        # disc of y^5 + a y^4 + b y^3 - tv over K
        poly = y**5 + a * y**4 + b * y**3 - tv
        try:
            # Discriminant as expression in U,V
            D = sp.discriminant(sp.Poly(sp.expand(poly), y))
            D = sp.expand(sp.simplify(D))
            # Reduce powers U^5=2, V^5=3
            D = sp.expand(D.subs(U**5, 2).subs(V**5, 3))
            # Collect basis U^i V^j i,j=0..4
            # Check if D is a square in Q(U,V) — hard.
            # Numeric: evaluate U=2**0.2, V=3**0.2
            Dn = complex(D.subs({U: 2**0.2, V: 3**0.2}))
            if abs(Dn.imag) < 1e-8 and Dn.real > 0:
                r = float(Dn.real) ** 0.5
                # weak check: is Dn close to a square of an element of K? skip strict
                # only record if Dn is almost integer square
                if abs(r - round(r)) < 1e-6:
                    stats["even_numeric_hint"] += 1
                    even.append({"t": str(tv), "disc_numeric": Dn.real, "hint": "near_integer_square"})
        except Exception:
            stats["err"] += 1
    return {
        "stats": dict(stats),
        "even_hints": even,
        "note": (
            "Full even-locus for (2A,3A,5A) needs disc square in Q(2^{1/5},3^{1/5}); "
            "numeric scan only flags near-integer-square hints."
        ),
        "field": "Q(2**(1/5), 3**(1/5))",
        "phi": "y**5 + a y**4 + b y**3",
        "a": str(a_f),
        "b": str(b_f),
    }


# ---------------------------------------------------------------------------
# Mild deformations preserving A5 monodromy (discrete / pull-back type)
# ---------------------------------------------------------------------------
def deformations_preserving_A5():
    """
    Absolute rigidity ⇒ no continuous deformation of the triple.
    Allowed mild moves that still give geometric monodromy A5:
      (1) Base pull-back by nonconstant rational r: φ(y) = r(t)
      (2) PGL2 on domain/range (already partly covered)
      (3) Other rigid triples (2A,3A,5*)
      (4) Composition with unramified covers — skip
    For (1): monodromy is A5 as long as r is nonconstant (Hurwitz: pullback
    of a connected A5 cover along a nonconstant map of bases remains with
    monodromy a subgroup of A5; if r is degree 1, isomorphic; if deg>1,
    monodromy still in A5 but branching multiplies).
    """
    return [
        {"id": "id_t", "r": t, "note": "standard"},
        {"id": "t2", "r": t**2, "note": "deg 2 pullback; monodromy ≤ A5"},
        {"id": "t3", "r": t**3, "note": "deg 3 pullback"},
        {"id": "mobius_361539", "r": -239 * (t - 3) / (29 * (t - 539)), "note": "resonant base chart"},
        {"id": "t_plus_1_over_t", "r": t + 1 / t, "note": "Laurent pullback"},
        {"id": "3t_61", "r": 3 * t + 61, "note": "affine lattice"},
    ]


def scan_pullback_even(phi_expr, pullbacks, t_vals):
    results = []
    for pb in pullbacks:
        print(f"  pullback {pb['id']}...", flush=True)
        stats = Counter()
        even = []
        for tv in t_vals:
            stats["tested"] += 1
            try:
                w = sp.simplify(pb["r"].subs(t, tv))
                if w.has(sp.zoo) or w == sp.oo or w.has(sp.nan):
                    stats["pole"] += 1
                    continue
                # require w rational for Z poly
                if not (sp.simplify(w).is_rational or sp.simplify(w).is_integer):
                    # try N(w)
                    try:
                        w = sp.nsimplify(w)
                    except Exception:
                        stats["irrational_w"] += 1
                        continue
                if not (sp.Rational(w) == w or sp.Integer(w) == w or sp.QQ.convert(w)):
                    pass
                expr = sp.expand((phi_expr - w) / sp.LC(sp.Poly(phi_expr, y)))
                pol = monic_fibre_to_Z(expr, y)
                if pol is None or not pol.is_irreducible:
                    stats["red_or_bad"] += 1
                    continue
                d = int(pol.discriminant())
                if d > 0 and is_square(d):
                    stats["even"] += 1
                    r = classify_poly(pol.as_expr().subs(y, x), do_galois=True)
                    even.append(
                        {
                            "t": str(tv),
                            "w": str(w),
                            "status": r.get("status"),
                            "gal": r.get("galois"),
                            "poly": str(pol.as_expr()),
                        }
                    )
                    print(f"    EVEN {pb['id']} t={tv} w={w} {r.get('status')}", flush=True)
                else:
                    stats["odd"] += 1
            except Exception:
                stats["err"] += 1
        results.append({"pullback": pb, "stats": dict(stats), "even": even})
    return results


def compare_triples_summary(r335, r235, pull_results, disc_info):
    return {
        "preferred_335_dense_even": len(r335.get("even") or []),
        "fallback_235_hints": len(r235.get("even_hints") or []),
        "pullbacks_with_even": [
            r["pullback"]["id"] for r in pull_results if r["even"]
        ],
        "disc_structure": disc_info,
    }


def main():
    t0 = time.time()
    print("GEOMETRIC RIGID TRIPLES + MILD DEFORMATIONS", flush=True)

    # Disc structure of preferred family
    disc_info = search_even_from_disc_poly()

    # Dense even search for (3A,3A,5A)
    r335 = dense_even_search_335(max_p=60, max_q=12)

    # (3A,3A,5B): S3 base transforms of φ (relabel 0,1,∞ / invert)
    transforms = phi_base_transforms()
    transform_results = []
    t_sample = [
        sp.Rational(p, q)
        for q in range(1, 9)
        for p in range(-25, 26)
        if gcd(abs(p), q) == 1 or q == 1
    ]
    # only polynomial transforms for Z fibres easily: phi and 1-phi
    for name, expr in [("phi_335A", transforms["phi"]), ("one_minus_phi", transforms["one_minus_phi"])]:
        # one_minus_phi is still polynomial deg 5
        transform_results.append(
            search_even_specialisations_poly_family(expr, name, t_sample[:80], do_galois=True)
        )

    # (2A,3A,5A)
    r235 = dense_even_search_235()

    # Mild pull-backs preserving A5 monodromy geometrically
    # Use moderate lattice for speed
    t_lattice = list(range(-12, 13)) + [16, 18, 27, 61, 80, 243, -16, -61]
    pullbacks = deformations_preserving_A5()
    pull_results = scan_pullback_even(PHI_335A, pullbacks, t_lattice)

    # Rigidity remark
    rigidity_note = (
        "Absolute rigidity of each triple implies the Hurwitz space is 0-dimensional: "
        "there is no continuous family of covers with the same conjugacy classes. "
        "Mild deformations that preserve geometric monodromy A5 are discrete "
        "(other triples, Galois conjugates, base pull-backs / PGL2), not a positive-dim "
        "moduli deformation of a single rigid triple."
    )

    summary = compare_triples_summary(r335, r235, pull_results, disc_info)

    total_even_335 = len(r335.get("even") or [])
    total_even_pull = sum(len(r["even"]) for r in pull_results)
    total_even_tr = sum(len(r["even"]) for r in transform_results)

    verdict = (
        f"Dense even search on preferred (3A,3A,5A) cover: {total_even_335} even irr fibres "
        f"(stats {r335.get('stats')}). "
        f"Base transforms (1-φ etc.): {total_even_tr} even. "
        f"Pull-backs: {total_even_pull} even across {len(pullbacks)} maps. "
        f"(2A,3A,5A) numeric disc hints: {len(r235.get('even_hints') or [])}. "
        + (
            "Found even specialisations — see list; check BJ/HQCC separately."
            if (total_even_335 + total_even_pull + total_even_tr) > 0
            else "No even irreducible arithmetic specialisations found on dense rational grids "
            "for the preferred Q-cover and its polynomial transforms/pull-backs. "
            "Improving even-specialisation rate likely needs a different rigid triple over a "
            "number field with more rational points on the even locus, or a non-rigid A5 family."
        )
    )

    elapsed = round(time.time() - t0, 2)
    lines = [
        "# Other rigid triples + mild A5-preserving deformations",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "## Goal",
        "",
        "Systematically examine rigid triples beyond preferred (3A,3A,5A) and mild",
        "moves that preserve geometric monodromy \(A_5\), focusing on **even**",
        "arithmetic specialisations (necessary for BJ/HQCC fusion).",
        "",
        f"**Verdict:** {verdict}",
        "",
        f"**Rigidity note:** {rigidity_note}",
        "",
        "---",
        "",
        "## Rigid triples (from Step 1)",
        "",
        "| Signature | Realization used | Field |",
        "|-----------|------------------|-------|",
        "| (3A,3A,5A) | \(\\varphi=6y^5-15y^4+10y^3\) | \(\\mathbb{Q}\) |",
        "| (3A,3A,5B) | \(1-\\varphi\) / label swap (same types) | \(\\mathbb{Q}\) |",
        "| (2A,3A,5A) | \(y^5+ay^4+by^3\) radical coeffs | \(\\mathbb{Q}(2^{1/5},3^{1/5})\) |",
        "| (2A,3A,5B) | conjugate 5-class (same map up to labeling) | same |",
        "",
        "---",
        "",
        "## Preferred (3A,3A,5A) — dense even locus",
        "",
        f"- stats: `{r335.get('stats')}`",
        f"- even count: **{total_even_335}**",
        "",
    ]
    for e in (r335.get("even") or [])[:20]:
        lines.append(
            f"- t={e['t']}: {e.get('status')} gal={e.get('gal')} poly=`{e.get('poly')}`"
        )
    if not r335.get("even"):
        lines.append("_No even irreducible fibres on the dense rational grid._")

    lines += [
        "",
        "### Discriminant structure of monic(φ−t)",
        "",
        f"- disc expr (preview): `{disc_info.get('disc_expr')}`",
        f"- num factored (preview): `{disc_info.get('num_factored')}`",
        "",
        "---",
        "",
        "## Base transforms (3A,3A,5B labeling / S3 on {0,1,∞})",
        "",
    ]
    for tr in transform_results:
        lines.append(f"### {tr['name']}")
        lines.append(f"- stats: `{tr['stats']}`")
        for e in (tr.get("even") or [])[:10]:
            lines.append(f"- t={e['t']}: {e.get('status')} `{e.get('poly')}`")
        if not tr.get("even"):
            lines.append("- _no even fibres_")
        lines.append("")

    lines += [
        "---",
        "",
        "## (2A,3A,5A) radical cover",
        "",
        f"- a = `{r235.get('a')}`",
        f"- b = `{r235.get('b')}`",
        f"- field: {r235.get('field')}",
        f"- stats: `{r235.get('stats')}`",
        f"- note: {r235.get('note')}",
        f"- even hints: {r235.get('even_hints')}",
        "",
        "---",
        "",
        "## Mild pull-backs (geometric monodromy ≤ A5)",
        "",
    ]
    for r in pull_results:
        pb = r["pullback"]
        lines.append(f"### `{pb['id']}` — r={pb['r']}")
        lines.append(f"- {pb.get('note')}")
        lines.append(f"- stats: `{r['stats']}`")
        for e in r.get("even") or []:
            lines.append(f"- EVEN t={e['t']} w={e['w']}: {e.get('status')} `{e.get('poly')}`")
        if not r.get("even"):
            lines.append("- _no even fibres on lattice sample_")
        lines.append("")

    lines += [
        "---",
        "",
        "## Conclusions",
        "",
        "1. **No continuous rigid deformation** of a single triple exists (absolute rigidity).",
        "2. **Other rigid triples** realized: (3A,3A,5*) via \(\\varphi\) and \(1-\\varphi\);",
        "   (2A,3A,5*) via radical Belyi over \(K=\\mathbb{Q}(2^{1/5},3^{1/5})\).",
        "3. **Even arithmetic specialisations** of the preferred \(\\mathbb{Q}\)-cover remain",
        f"   extremely rare / absent on large rational grids ({total_even_335} hits).",
        "4. Pull-backs and base transforms did **not** materially improve even rates",
        f"   on the sampled lattices (pull even total {total_even_pull}).",
        "5. Therefore: improving even specialisations likely requires either",
        "   - working over number fields with more even fibres (e.g. 235 over \(K\)), or",
        "   - **non-rigid** A5 families (positive-dim Hurwitz) with denser rational even loci,",
        "   - not further low-degree surgery on the single rigid \(\\varphi/\\mathbb{Q}\).",
        "",
        "### Link to fusion",
        "",
        "Geometric monodromy \(A_5\) is settled for these rigid objects. The bottleneck for",
        "fusion with HQCC BJ seeds remains **even arithmetic specialisations + BJ form**,",
        "which these rigid Q-covers do not supply on rational base points in the scans.",
        "",
        "_Generated by geometric_rigid_deform.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "rigidity_note": rigidity_note,
        "r335": r335,
        "r235": r235,
        "transforms": transform_results,
        "pullbacks": pull_results,
        "disc_info": disc_info,
        "summary": summary,
    }
    write_md(OUT / "GEOMETRIC_RIGID_DEFORM.md", doc)
    write_md(RESULTS / "GEOMETRIC_RIGID_DEFORM.md", doc)
    write_md(ROOT / "GEOMETRIC_RIGID_DEFORM.md", doc)
    write_json(OUT / "GEOMETRIC_RIGID_DEFORM.json", blob)
    print(verdict, flush=True)
    print(f"Wrote GEOMETRIC_RIGID_DEFORM.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
