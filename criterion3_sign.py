"""
Criterion 3 — Sign character of S_n vs ternary invariants.

Experimentally correlates:
  - disc square (parity / image in A_n)
  - det(M), signature-like matrix invariants
  - ternary content of coefficients
  - whether residual Z/3 structure is present

Goal of a theorem: sgn ∘ ρ = 1 for monodromy ρ arising from HQCC data.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import sympy as sp

from lib.common import (
    OUT,
    RESULTS,
    charpoly_matrix,
    classify_poly,
    is_square,
    write_json,
    write_md,
)

x = sp.symbols("x")


def T5(a, b, c, d, e=0, f=0):
    return sp.Matrix([
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [a, 0, 0, b, e],
        [0, 0, 0, 0, 1],
        [c, f, 0, d, 0],
    ])


def ternary_weight(M: sp.Matrix) -> int:
    """Count entries divisible by 3."""
    w = 0
    for i in range(M.rows):
        for j in range(M.cols):
            if int(M[i, j]) % 3 == 0 and int(M[i, j]) != 0:
                w += 1
    return w


def det_sign(M: sp.Matrix) -> int:
    try:
        d = int(M.det())
        if d > 0:
            return 1
        if d < 0:
            return -1
        return 0
    except Exception:
        return 0


def pfaffian_like_skew(M: sp.Matrix):
    """If M is skew-symmetric odd size, det=0; for even, related to pfaffian^2=det."""
    n = M.rows
    skew = all(int(M[i, j] + M[j, i]) == 0 for i in range(n) for j in range(n))
    return skew


def sample_correlations(max_n: int = 5000) -> dict:
    """
    Sample T5 over small lattice; record (disc_square, det_sign, ternary_weight, has_3cycle).
    """
    pool = [0, 1, -1, 3, -3, 9, 61, 80]
    rows = []
    stats = Counter()
    count = 0
    for a, b, c, d in itertools.product(pool, repeat=4):
        for e, f in itertools.product([0, 1, -1, 3], repeat=2):
            count += 1
            if count > max_n:
                break
            M = T5(a, b, c, d, e, f)
            chi = charpoly_matrix(M)
            pol = sp.Poly(chi, x, domain=sp.ZZ)
            if pol.LC() != 1 and pol.LC() != -1:
                continue
            if pol.LC() == -1:
                pol = sp.Poly(-pol.as_expr(), x, domain=sp.ZZ)
            if pol.degree() != 5:
                continue
            if not pol.is_irreducible:
                stats["red"] += 1
                continue
            disc = int(pol.discriminant())
            sq = is_square(disc)
            # cheap 3-cycle check: few primes
            has3 = False
            for p in [2, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
                if disc % p == 0:
                    continue
                try:
                    facs = sp.factor_list(pol.as_expr(), modulus=p)
                    degs = []
                    for fct, m in facs[1]:
                        degs.extend([int(sp.degree(fct))] * int(m))
                    if 3 in degs:
                        has3 = True
                        break
                except Exception:
                    pass
            tw = ternary_weight(M)
            ds = det_sign(M)
            sk = pfaffian_like_skew(M)
            row = {
                "params": (a, b, c, d, e, f),
                "disc_square": sq,
                "det_sign": ds,
                "ternary_weight": tw,
                "has_3cycle_mod": has3,
                "skew": sk,
                "poly": str(pol.as_expr()),
            }
            rows.append(row)
            stats["irr"] += 1
            if sq:
                stats["sq"] += 1
            if has3:
                stats["has3"] += 1
            if sq and has3:
                stats["sq_and_3"] += 1
            # correlation buckets
            stats[f"det{ds}_sq{sq}"] += 1
            stats[f"tw{tw}_sq{int(sq)}"] += 1
        if count > max_n:
            break

    # correlation summary: P(sq | det_sign), P(sq | ternary_weight>=k)
    corr = {}
    for ds in (-1, 0, 1):
        sub = [r for r in rows if r["det_sign"] == ds]
        if sub:
            corr[f"P(sq|det={ds})"] = sum(1 for r in sub if r["disc_square"]) / len(sub)
            corr[f"n|det={ds}"] = len(sub)
    for k in range(0, 6):
        sub = [r for r in rows if r["ternary_weight"] >= k]
        if sub:
            corr[f"P(sq|ternary_weight>={k})"] = sum(1 for r in sub if r["disc_square"]) / len(sub)
    sub3 = [r for r in rows if r["has_3cycle_mod"]]
    if sub3:
        corr["P(sq|has3)"] = sum(1 for r in sub3 if r["disc_square"]) / len(sub3)
    sub_sq = [r for r in rows if r["disc_square"]]
    if sub_sq:
        corr["P(has3|sq)"] = sum(1 for r in sub_sq if r["has_3cycle_mod"]) / len(sub_sq)

    return {"stats": dict(stats), "correlations": corr, "n_rows": len(rows), "sample": rows[:50]}


def load_catalogue_parity() -> dict:
    """Known A5/A6 vs D5/S5 from results JSON — sign is trivial exactly when disc^2."""
    cats = []
    for fname, label in [
        ("DEFORM_M.json", "deform"),
        ("SUMMARY_unique.json", "lattice"),
        ("A6_T6.json", "a6"),
    ]:
        p = RESULTS / fname
        if not p.exists():
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        if "A5" in d:
            for h in d["A5"]:
                cats.append({"src": label, "gal": "A5", "disc_sq": True, "poly": h.get("poly")})
        if "A6" in d:
            for h in d["A6"]:
                cats.append({"src": label, "gal": "A6", "disc_sq": True, "poly": h.get("poly")})
        if "unique_A5" in d:
            for h in d["unique_A5"]:
                cats.append({"src": label, "gal": "A5", "disc_sq": True, "poly": h.get("poly")})
        if "D5" in d:
            for h in d["D5"]:
                cats.append({"src": label, "gal": "D5", "disc_sq": True, "poly": h.get("poly")})
    # base odd examples
    cats.append({"src": "base", "gal": "S5", "disc_sq": False, "poly": "x**5 + 3*x**3 - 3*x**2 - 4889"})
    cats.append({"src": "base", "gal": "S4xC2", "disc_sq": False, "poly": "x**6 + 3*x**4 - 3*x**2 - 4889"})
    return {"catalogue": cats, "n": len(cats)}


def theory_doc(corr_blob: dict, cat: dict) -> str:
    lines = [
        "# Criterion 3 — Sign character and ternary invariants",
        "",
        "## Goal",
        "",
        "Link the sign character \(\\operatorname{sgn}: S_n\\to\\{\\pm1\\}\) to an invariant",
        "of the ternary / HQCC structure so that monodromy representations \(\\rho\)",
        "satisfy \(\\operatorname{sgn}\\circ\\rho = 1\) (image in \(A_n\)).",
        "",
        "Equivalently for characteristic polynomials: \(\\operatorname{disc}(\\chi)\) is a square.",
        "",
        "## Classical dictionary",
        "",
        "| Arithmetic | Group-theoretic |",
        "|------------|-----------------|",
        "| disc square | \(\\mathrm{Gal}\\le A_n\) |",
        "| transposition Frobenius | odd permutation |",
        "| 3-cycle Frobenius | even; kills pure \(D_5\) at deg 5 |",
        "",
        "## Experimental correlations (T5 sample)",
        "",
        f"Irrreducible samples: {corr_blob.get('n_rows')}",
        "",
        "### Raw stats",
        f"```\n{json.dumps(corr_blob.get('stats'), indent=2)}\n```",
        "",
        "### Conditional probabilities",
        f"```\n{json.dumps(corr_blob.get('correlations'), indent=2)}\n```",
        "",
        "## Catalogue parity check",
        "",
        f"Known hits loaded: {cat.get('n')}",
        "",
        "- All stored \(A_5/A_6/D_5\) hits have disc² by construction of the catalogues.",
        "- Base structural matrices \(M\), \(T_6\) are **odd** (sgn nontrivial) despite ternary entries.",
        "",
        "**Conclusion from data:** ternary weight alone does **not** force trivial sign.",
        "Any Criterion-3 theorem needs a stronger invariant than “has a factor of 3 in \(M\).”",
        "",
        "## Candidate invariants for a future theorem",
        "",
        "| Invariant | Mechanism |",
        "|-----------|-----------|",
        "| Volume form on root space preserved by monodromy | forces \(\\operatorname{sgn}=1\) |",
        "| T-complementarity involution commuting with monodromy | may force evenness or a fixed twist |",
        "| Quadratic character of flux lattice (4880, 61) | \(\\operatorname{sgn}\\cdot\\chi_{\\mathrm{model}}=1\) |",
        "| Norm construction from \(\\mathbb{Q}(\\omega)\) | built-in disc control |",
        "",
        "## Status",
        "",
        "- Sign/evenness is **exactly** the disc² gate used in scans.",
        "- **No ternary invariant yet identified** that implies disc² for all structural \(M\).",
        "- Correlations above quantify how det(M) / ternary_weight relate to disc² empirically.",
        "",
        "_Generated by criterion3_sign.py_",
    ]
    return "\n".join(lines)


def main():
    print("Criterion 3: sign / ternary correlations", flush=True)
    corr = sample_correlations(max_n=4000)
    print(f"  correlations: {corr['correlations']}", flush=True)
    cat = load_catalogue_parity()
    doc = theory_doc(corr, cat)
    write_md(OUT / "CRITERION3_SIGN.md", doc)
    write_md(RESULTS / "CRITERION3_SIGN.md", doc)
    write_json(OUT / "CRITERION3_SIGN.json", {"correlations": corr, "catalogue": cat})
    print(f"Wrote {OUT / 'CRITERION3_SIGN.md'}", flush=True)


if __name__ == "__main__":
    main()
