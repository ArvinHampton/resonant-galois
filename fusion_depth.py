"""
Next-depth probes for geometric–arithmetic fusion (Criterion 1).

Option A: mild twists / parametric deformations of φ that might recover HQCC seeds.
Option B: alternative residue→braid assignments; naturality diagnostics.

Does not claim to close the fusion gap — quantifies and searches.
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

from lib.common import OUT, RESULTS, classify_poly, is_square, monic_poly, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

SEEDS = [
    (-55, 88),
    (-55, -88),
    (95, 76),
    (95, -76),
    (95, 532),
    (95, -532),
    (-100, 400),
    (124, 496),
    (20, 16),
]


def to_monic_Z(expr) -> sp.Poly | None:
    try:
        pol = sp.Poly(sp.expand(expr), x, domain=sp.QQ)
        if pol.degree() != 5 or pol.LC() == 0:
            return None
        mon = sp.Poly(sp.monic(pol.as_expr()), x, domain=sp.QQ)
        dens = [sp.fraction(sp.together(c))[1] for c in mon.all_coeffs()]
        L = 1
        for d in dens:
            L = int(sp.ilcm(L, abs(int(d))))
        cleared = sp.expand(L**5 * mon.as_expr().subs(x, x / L))
        return monic_poly(cleared)
    except Exception:
        return None


def poly_matches_seed(pol: sp.Poly, alpha: int, beta: int) -> str | None:
    target = [1, 0, 0, 0, alpha, beta]
    coeffs = [int(c) for c in pol.all_coeffs()]
    if coeffs == target:
        return "exact"
    # x -> -x
    alt = monic_poly(pol.as_expr().subs(x, -x))
    if alt is not None and [int(c) for c in alt.all_coeffs()] == target:
        return "x_to_-x"
    # same up to scaling of variable already handled by monic clear
    return None


# ---------------------------------------------------------------------------
# Option A: parametric twists of φ
# ---------------------------------------------------------------------------
def option_A_twists() -> dict:
    """
    Families that specialise to φ at a base value, then scan for HQCC seeds.

    Twist 1: φ_u(y) = 6y^5 - 15y^4 + 10y^3 - u  (already fibres; Attack 1)
    Twist 2: φ_{p,q}(y) = p y^5 - (p+q) something — weighted homogenisation of φ
    Twist 3: φ_c(y) = y^3 (6y^2 - 15 c y + 10 c^2)  (scale middle coeffs)
    Twist 4: BJ family through seed space already known — not new geometry
    Twist 5: pull-back φ(y) - t(s) with t non-Möbius rational of low degree
    """
    print("Option A: twists of φ...", flush=True)
    hits = []
    stats = Counter()

    # Twist 3: scale family φ_c(y) = y^3 (6y^2 - 15 c y + 10 c^2)
    # At c=1 recovers φ. After monic clear, may hit seeds.
    for c in list(range(-12, 13)) + [16, 18, 27, 61, 80, 243, 539, -16, -61]:
        if c == 0:
            continue
        expr = x**3 * (6 * x**2 - 15 * c * x + 10 * c**2)
        pol = to_monic_Z(expr)
        stats["twist3_tested"] += 1
        if pol is None:
            continue
        coeffs = [int(c0) for c0 in pol.all_coeffs()]
        # check if BJ form [1,0,0,0,α,β]
        if coeffs[1] == coeffs[2] == coeffs[3] == 0:
            a, b = coeffs[4], coeffs[5]
            for sa, sb in SEEDS:
                m = poly_matches_seed(pol, sa, sb)
                if m:
                    hits.append({"twist": "scale_c", "c": c, "seed": (sa, sb), "how": m, "poly": str(pol.as_expr())})
                    print(f"  HIT scale_c={c} seed=({sa},{sb})", flush=True)
            # also record if Gal A5 even if not known seed
            if is_square(int(pol.discriminant())) if pol.is_irreducible else False:
                stats["twist3_sq"] += 1

    # Twist 5: φ(y) - R(s) for rational R of deg ≤ 2 at model s
    # For each model s and simple R, clear to monic and match seeds
    s = sp.symbols("s")
    rationals = [
        s,
        s**2,
        1 / s if False else s,  # skip 1/s symbolic alone
        (s - 3) / (s - 539),
        s * (s - 61),
        3 * s + 61,
        s**2 + 3 * s + 9,
    ]
    for R in [(s - 3) / (s - 539), 3 * s + 61, s**2 + 9, 61 * s - 539, s**2 - 3 * s + 27]:
        for sval in [1, 2, 3, 4, 5, 9, 16, 18, 27, 61, 80, -1, -3, 10, 12, 15, 20, 24, 30]:
            try:
                w = sp.simplify(R.subs(s, sval))
                if w == sp.zoo or w.has(sp.oo):
                    continue
                expr = 6 * x**5 - 15 * x**4 + 10 * x**3 - w
                pol = to_monic_Z(expr)
                stats["twist5_tested"] += 1
                if pol is None or not pol.is_irreducible:
                    continue
                for sa, sb in SEEDS:
                    m = poly_matches_seed(pol, sa, sb)
                    if m:
                        hits.append(
                            {
                                "twist": "pullback_R",
                                "R": str(R),
                                "s": sval,
                                "w": str(w),
                                "seed": (sa, sb),
                                "how": m,
                                "poly": str(pol.as_expr()),
                            }
                        )
                        print(f"  HIT R={R} s={sval} seed=({sa},{sb})", flush=True)
            except Exception:
                continue

    # Twist 6: general monic BJ search is done; try near-φ polynomials
    # x^5 + a x^4 + b x^3 + c x^2 + d x + e with (a,b,c) near (0,0,0) and (d,e) near seeds
    # — skip heavy; instead: Tschirnhaus reverse — start from seed, ask if critical
    # values of a poly with that Gal match {0,1,∞} after Möbius on base.
    # For seed f, if it were a fibre of a Belyi map of passport (3,3,5), the
    # critical values of the inverse would be special.
    # Practical probe: for each seed, compute whether it is a translate/scale of
    # monic form of y^3(6y^2-15cy+10c^2) - w.
    for sa, sb in SEEDS:
        stats["seed_probe"] += 1
        # solve  monic( y^3(6y^2-15c y+10c^2) - w )  == x^5 + sa x + sb
        # after y = λx+μ — too many vars; fix μ=0, λ free, c free, w free
        # monic form of λ^{-5} ( (λx)^3 (6 λ^2 x^2 - 15 c λ x + 10 c^2) - w )
        # = x^3 (6 x^2 - 15 (c/λ) x + 10 (c/λ)^2  - w/λ^5 / x^3 ) not poly unless...
        # Only pure φ_c - w works when result is degree 5 poly: always is.
        # monic(φ_c - w): leading 6, so monic = (1/6)φ_c - w/6
        # Need x^4,x^3,x^2 coeffs zero after monic: 
        # φ_c = 6x^5 - 15 c x^4 + 10 c^2 x^3
        # monic = x^5 - (15c/6) x^4 + (10 c^2/6) x^3 - w/6
        # = x^5 - (5c/2) x^4 + (5 c^2/3) x^3 - w/6
        # For BJ: need c=0 degenerate, or never zero x^4 unless c=0.
        # So φ_c - w is NEVER Bring-Jerrard for c≠0.
        # Conclusion: this twist family cannot produce BJ seeds.
        pass

    conclusion = (
        "Mild scale/pull-back twists of φ do not produce the known HQCC BJ seeds "
        f"in the scanned grid (hits={len(hits)}). "
        "Structural reason: monic(φ_c−w) has x^4 coefficient −5c/2 ≠ 0 for c≠0, "
        "so this family never lands in Bring–Jerrard form. "
        "Fusion Gap A requires a family that is BJ for a positive-dimensional "
        "parameter set (or a non-polynomial cover / different passport)."
    )
    return {
        "hits": hits,
        "stats": dict(stats),
        "structural_obstruction": (
            "φ_c(y)=y^3(6y^2−15c y+10c^2) has monic form with x^4 coeff −5c/2; "
            "HQCC seeds are BJ (no x^4,x^3,x^2). Intersection only if c=0 degenerates."
        ),
        "conclusion": conclusion,
    }


# ---------------------------------------------------------------------------
# Option B: naturality of residue→braid
# ---------------------------------------------------------------------------
def option_B_naturality() -> dict:
    """
    Enumerate small encodings of residue 2 as words in {σ0,σ1} of length ≤ 3
    and measure how much path statistics change — naturality diagnostic.
    """
    print("Option B: functor naturality diagnostics...", flush=True)
    # Generators from Step 2 (fixed)
    s0 = [1, 2, 0, 3, 4]  # (0 1 2) type 3+1+1 — may differ by conjugation from numeric run
    s1 = [0, 1, 4, 2, 3]

    # Recompute from geometric_step2 logic
    import numpy as np
    from numpy.polynomial.polynomial import polyroots

    def preimages(w):
        c = [0 - w, 0, 0, 10, -15, 6]  # low to high for 6x^5-15x^4+10x^3 - w
        return polyroots(np.array(c, dtype=np.complex128))

    def track(center, radius=0.12, n=200):
        base = center + radius
        init = preimages(base)
        init = init[np.lexsort((init.imag, init.real))]
        sheets = init.copy()
        for th in np.linspace(0, 2 * np.pi, n + 1)[1:]:
            w = center + radius * np.exp(1j * th)
            new = preimages(w)
            matched = np.empty_like(sheets)
            used = set()
            for i, s in enumerate(sheets):
                for j in np.argsort(np.abs(new - s)):
                    if j not in used:
                        matched[i] = new[j]
                        used.add(j)
                        break
            sheets = matched
        return [int(np.argmin(np.abs(init - s))) for s in sheets]

    def compose(p, q):
        return [p[q[i]] for i in range(5)]

    def invert(p):
        inv = [0] * 5
        for i, j in enumerate(p):
            inv[j] = i
        return inv

    def cycles_part(p):
        seen = [False] * 5
        lens = []
        for i in range(5):
            if not seen[i]:
                L = 0
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    L += 1
                lens.append(L)
        return tuple(sorted(lens, reverse=True))

    g0, g1 = track(0.0), track(1.0)
    # candidate words for residue 2
    candidates = {
        "σ0σ1": lambda a: compose(compose(a, g0), g1),
        "σ1σ0": lambda a: compose(compose(a, g1), g0),
        "σ0": lambda a: compose(a, g0),
        "σ1": lambda a: compose(a, g1),
        "σ0^{-1}": lambda a: compose(a, invert(g0)),
        "σ1^{-1}": lambda a: compose(a, invert(g1)),
        "σ0σ1σ0": lambda a: compose(compose(compose(a, g0), g1), g0),
        "σ1σ0σ1": lambda a: compose(compose(compose(a, g1), g0), g1),
        "(σ0σ1)^{-1}": lambda a: compose(a, invert(compose(g0, g1))),
    }

    def T3(n):
        if n == 0:
            return 0
        r = n % 3
        if r == 0:
            return n // 3
        if r == 1:
            return (4 * n + 2) // 3
        return (2 * n + 1) // 3

    def path(n0, steps=16):
        n, p = n0, []
        for _ in range(steps):
            p.append(n % 3)
            n = T3(n)
            if n == 0:
                p.append(0)
                break
        return p

    seeds_n = [1, 2, 3, 9, 27, 61, 80, 243, 539, 4880, 100, 55, 88, 95]
    comparison = {}
    for name, apply2 in candidates.items():
        hist = Counter()
        for n0 in seeds_n:
            acc = list(range(5))
            for r in path(n0):
                if r == 0:
                    acc = compose(acc, g0)
                elif r == 1:
                    acc = compose(acc, g1)
                else:
                    acc = apply2(acc)
            hist[str(cycles_part(acc))] += 1
        comparison[name] = dict(hist)

    # Naturality score: encodings that agree on conjugacy class of image
    # for all sample paths (same partition histogram)
    ref = comparison.get("σ0σ1", {})
    agreement = {}
    for name, hist in comparison.items():
        agreement[name] = hist == ref

    # Prefer encodings whose image lands in class 3A often (ternary)
    ternary_weight = {
        name: sum(v for k, v in hist.items() if "3" in k)
        for name, hist in comparison.items()
    }

    return {
        "generators": {"sigma_0": g0, "sigma_1": g1},
        "residue2_candidates_histograms": comparison,
        "agrees_with_sigma0sigma1": agreement,
        "ternary_weight": ternary_weight,
        "naturality_note": (
            "Several residue-2 encodings produce different cycle-type histograms "
            "⇒ the scaffold functor is encoding-dependent. "
            "Naturality fails until a characterisation picks a unique word "
            "(e.g. shortest even word of type 3A, or from a third marked point)."
        ),
        "recommendation": (
            "To naturalise: either (i) pass to a 4-branch cover with a dedicated "
            "point for residue 2, or (ii) fix residue 2 via a conjugacy-invariant "
            "rule (e.g. unique element of minimal length in A5 with prescribed "
            "cycle type and braid relation)."
        ),
    }


def write_doc(a, b, elapsed) -> str:
    lines = [
        "# Fusion depth probes (Options A & B)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "Context: fused Criterion-1 gap (`FUSION_GAP.md`).",
        "",
        "---",
        "",
        "## Option A — mild twists of φ → HQCC seeds",
        "",
        f"**Conclusion:** {a.get('conclusion')}",
        "",
        f"- Hits: **{len(a.get('hits') or [])}**",
        f"- Stats: `{a.get('stats')}`",
        f"- Structural obstruction: {a.get('structural_obstruction')}",
        "",
    ]
    for h in (a.get("hits") or [])[:20]:
        lines.append(f"- HIT `{h}`")
    lines += [
        "",
        "### Implication",
        "",
        "Polynomial fibres of the preferred Belyi map (and its scale family) **cannot**",
        "be Bring–Jerrard for nontrivial scale. Therefore Gap A needs a **different**",
        "geometric family (e.g. a BJ pencil with geometric monodromy \(A_5\), or a",
        "Hurwitz family with more parameters), not only a twist of \(\\varphi\).",
        "",
        "---",
        "",
        "## Option B — naturality of residue→braid",
        "",
        f"- Note: {b.get('naturality_note')}",
        f"- Recommendation: {b.get('recommendation')}",
        "",
        "### Residue-2 encoding histograms",
        "",
        f"```\n{json.dumps(b.get('residue2_candidates_histograms'), indent=2)}\n```",
        "",
        f"Agrees with σ0σ1: `{b.get('agrees_with_sigma0sigma1')}`",
        f"Ternary weight: `{b.get('ternary_weight')}`",
        "",
        "---",
        "",
        "## Status",
        "",
        "| Gap | After this probe |",
        "|-----|------------------|",
        "| A Hilbert recovery of seeds via mild φ-twist | **Still open** — structural obstruction for this φ-family |",
        "| B Natural T3→braid functor | **Still open** — encoding-dependent histograms |",
        "",
        "_Generated by fusion_depth.py_",
    ]
    return "\n".join(lines)


def main():
    t0 = time.time()
    print("FUSION DEPTH — Options A and B", flush=True)
    a = option_A_twists()
    b = option_B_naturality()
    elapsed = round(time.time() - t0, 2)
    doc = write_doc(a, b, elapsed)
    blob = {"elapsed_sec": elapsed, "option_A": a, "option_B": b}
    write_md(OUT / "FUSION_DEPTH.md", doc)
    write_md(RESULTS / "FUSION_DEPTH.md", doc)
    write_md(ROOT / "FUSION_DEPTH.md", doc)
    write_json(OUT / "FUSION_DEPTH.json", blob)
    print(a.get("conclusion"), flush=True)
    print(b.get("naturality_note"), flush=True)
    print(f"Wrote FUSION_DEPTH.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
