"""
Criterion 1 — Natural algebraic objects linked to HQCC / resonant data.

Constructs explicit one-parameter families and branch-cycle style polynomials
inspired by HQCC ternary branches, then specialises at model integers and
runs the standard Gal pipeline.

This is constructive scaffolding toward a monodromy theorem, not the theorem itself.
"""
from __future__ import annotations

import itertools
import sys
from collections import Counter

import sympy as sp

from lib.common import OUT, RESULTS, classify_poly, write_json, write_md

x, t = sp.symbols("x t")

# HQCC / Ad branches as Möbius maps n |-> (A n + B)/C
BRANCHES = {
    "div3": (1, 0, 3),
    "Ad_plus": (3, 1, 1),   # 3n+1
    "Ad_minus": (3, -1, 1),  # 3n-1
    "Syr1": (4, 2, 3),
    "Syr2": (2, 1, 3),
}

MODEL_T = [1, 3, 9, 16, 18, 27, 61, 80, 243, 539, 4880, -3, -9, -61]


def family_BJ_model() -> list[tuple[str, sp.Expr]]:
    """Bring–Jerrard style x^n + a x + b with model a,b — geometric monodromy often S_n/A_n."""
    out = []
    for n in (5, 6):
        for a, b in itertools.product(MODEL_T, MODEL_T):
            if b == 0:
                continue
            out.append((f"BJ_n{n}_a{a}_b{b}", x**n + a * x + b))
    return out


def family_hqcc_resultant_style() -> list[tuple[str, sp.Expr]]:
    """
    Primitive-element style: eliminate y from
      y^3 - 3*s*y - t = 0  (trig/cubic related to Z/3)
      x = y + m/y
    Resultant gives a deg-6 poly in x for each (s,m,t).
    Classic path to dihedral/A_n specialisations.
    """
    y = sp.symbols("y")
    out = []
    for s in [1, 3, 9, 61]:
        for m in [1, 3, 9]:
            for tv in [1, 3, 16, 61, 80, 539]:
                # resultant_y( y^3 - 3 s y - t,  y^2 - x y + m )
                f = y**3 - 3 * s * y - tv
                g = y**2 - x * y + m
                r = sp.resultant(f, g, y)
                out.append((f"res_s{s}_m{m}_t{tv}", sp.expand(r)))
    return out


def family_branch_composition_charpoly() -> list[tuple[str, sp.Expr]]:
    """
    2×2 Möbius lifts of HQCC branches, block-assembled to 6×6 with parameter t,
    take char poly, specialise t at model values.
    """
    out = []

    def mob(A, B, C):
        return sp.Matrix([[A, B], [0, C]])

    B0 = mob(*BRANCHES["div3"])
    B1 = mob(*BRANCHES["Ad_plus"])
    B2 = mob(*BRANCHES["Ad_minus"])
    for tval in MODEL_T:
        M = sp.zeros(6)
        M[0:2, 0:2] = B0
        M[2:4, 2:4] = B1
        M[4:6, 4:6] = B2
        # couple with t
        M[1, 2] = tval
        M[2, 1] = 3
        M[3, 4] = 1
        M[4, 3] = tval
        chi = sp.expand(M.charpoly(x).as_expr())
        out.append((f"mob_block_t{tval}", chi))
    return out


def family_rigid_An_specialisations() -> list[tuple[str, sp.Expr]]:
    """
    Known rigid / standard An seeds deformed by model integers.
    Deg 5: x^5 + 20x + 16 is A5; deform 20,16 toward model.
    Deg 6: x^6 + x + 1 often S6; try disc-square deformations.
    """
    out = []
    for p in [16, 20, 3, 61, 80, 243, 539, -16, -20, 1, -1, 9]:
        for q in [16, 20, 3, 61, 80, 1, -1, 9, -3, 539]:
            out.append((f"nearA5_p{p}_q{q}", x**5 + p * x + q))
    # icosahedral-ish
    for m in MODEL_T:
        out.append((f"x5_5m_n{m}", x**5 + 5 * m * x**3 + 5 * m**2 * x + 3))
    return out


def run_family(name: str, pairs: list[tuple[str, sp.Expr]], max_items: int = 500) -> dict:
    print(f"  Family {name}: {min(len(pairs), max_items)} items", flush=True)
    stats = Counter()
    hits = []
    even = []
    for label, expr in pairs[:max_items]:
        stats["tested"] += 1
        try:
            rec = classify_poly(expr, do_galois=True)
        except Exception as e:
            stats["error"] += 1
            continue
        rec["label"] = label
        rec["family"] = name
        st = rec.get("status") or ""
        stats[st] += 1
        if st.startswith("HIT_A"):
            hits.append(rec)
            print(f"    *** {st} *** {rec['poly']} [{label}]", flush=True)
        elif rec.get("disc_square") and rec.get("irreducible"):
            even.append(rec)
    return {
        "family": name,
        "stats": dict(stats),
        "A_hits": hits,
        "even_sample": even[:40],
    }


def theory_doc(results: list[dict]) -> str:
    lines = [
        "# Criterion 1 — Canonical / HQCC-linked algebraic objects",
        "",
        "## Goal",
        "",
        "A natural object (cover, representation, or moduli space) **canonically**",
        "associated with HQCC / the resonant ring, whose monodromy is **proved**",
        "to be alternating (or to contain \(A_n\)).",
        "",
        "## Constructive scaffolds implemented",
        "",
        "### 1. HQCC Möbius block char polys",
        "2×2 lifts of branches \(\\{n/3,\\,3n\\pm1\\}\) assembled into 6×6 with model coupling \(t\).",
        "",
        "### 2. Cubic resultant families (Z/3 + \(x=y+m/y\))",
        "Classical path from cyclic cubic data to deg-6 polys; specialise \((s,m,t)\) at model ints.",
        "",
        "### 3. Bring–Jerrard / near-rigid An seeds",
        "Standard geometric monodromy sources; specialise coefficients in the model lattice.",
        "",
        "## What would count as a theorem (not yet)",
        "",
        "> There exists a finite cover \(X\\to Y\) (or family \(f_t\\in\\mathbb{Q}(t)[x]\))",
        "> built **only** from HQCC branch data such that the geometric monodromy is \(A_n\).",
        "> Model integers arise as specialisations via Hilbert irreducibility.",
        "",
        "## Experimental specialisation results",
        "",
    ]
    total_hits = []
    for r in results:
        lines.append(f"### {r['family']}")
        lines.append(f"- stats: `{r['stats']}`")
        lines.append(f"- An hits: {len(r['A_hits'])}")
        for h in r["A_hits"]:
            lines.append(f"  - **{h['status']}** `{h['poly']}` label={h.get('label')}")
            total_hits.append(h)
        lines.append("")
    lines += [
        f"**Total An hits in Criterion-1 scaffolds: {len(total_hits)}**",
        "",
        "## Status",
        "",
        "- Scaffolds exist and feed the same Gal pipeline.",
        "- **No proof** that monodromy is alternating for a canonical HQCC cover.",
        "- Hits here are still *specialisations*, not a monodromy theorem.",
        "",
        "_Generated by criterion1_hqcc.py_",
    ]
    return "\n".join(lines)


def main():
    print("Criterion 1: HQCC-linked families", flush=True)
    results = []
    results.append(run_family("BJ_model", family_BJ_model(), max_items=400))
    results.append(run_family("hqcc_resultant", family_hqcc_resultant_style(), max_items=200))
    results.append(run_family("mob_block", family_branch_composition_charpoly(), max_items=50))
    results.append(run_family("near_rigid", family_rigid_An_specialisations(), max_items=300))
    doc = theory_doc(results)
    write_md(OUT / "CRITERION1_HQCC.md", doc)
    write_md(RESULTS / "CRITERION1_HQCC.md", doc)
    write_json(OUT / "CRITERION1_HQCC.json", {"results": results})
    print(f"Wrote {OUT / 'CRITERION1_HQCC.md'}", flush=True)


if __name__ == "__main__":
    main()
