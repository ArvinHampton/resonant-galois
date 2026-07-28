"""
Criterion 2 — Axioms on structural matrices.

Evenness obstruction: matrices satisfying "ternary structural" axioms
need NOT have square discriminant (base M, base T6).

Also searches for stronger axiom subclasses where disc^2 holds more often.
"""
from __future__ import annotations

import itertools
import sys
from collections import Counter

import sympy as sp

sys.path.insert(0, str(__file__).rsplit("\\", 1)[0] if "\\" in __file__ else str(__file__).rsplit("/", 1)[0])
from lib.common import (  # noqa: E402
    OUT,
    RESULTS,
    charpoly_matrix,
    classify_poly,
    is_square,
    write_json,
    write_md,
)

x = sp.symbols("x")


def T5(a, b, c, d, e=0, f=0) -> sp.Matrix:
    """Successful deg-5 structural template (companion + couplings)."""
    return sp.Matrix([
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [a, 0, 0, b, e],
        [0, 0, 0, 0, 1],
        [c, f, 0, d, 0],
    ])


def T6(a, b, c, d, e=0, f=0, g=0) -> sp.Matrix:
    return sp.Matrix([
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [a, 0, 0, 0, b, e],
        [0, 0, 0, 0, 0, 1],
        [c, f, 0, 0, d, g],
    ])


def G4H_M() -> sp.Matrix:
    return sp.Matrix([
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [3, 0, 0, 80, 0],
        [0, 0, 0, 0, 1],
        [61, 0, 0, -3, 0],
    ])


# ---------------------------------------------------------------------------
# Axiom predicates (boolean properties of integer matrices)
# ---------------------------------------------------------------------------
def axiom_companion_chain(M: sp.Matrix) -> bool:
    """Superdiagonal 1's on a companion-like prefix."""
    n = M.rows
    for i in range(n - 3):
        if M[i, i + 1] != 1:
            return False
        for j in range(n):
            if j != i + 1 and M[i, j] != 0:
                return False
    return True


def axiom_ternary_entry(M: sp.Matrix) -> bool:
    """Some entry is in {±3, ±9, ±27, ±81, ±243}."""
    tern = {3, -3, 9, -9, 27, -27, 81, -81, 243, -243}
    for i in range(M.rows):
        for j in range(M.cols):
            if int(M[i, j]) in tern:
                return True
    return False


def axiom_model_entries(M: sp.Matrix) -> bool:
    """All nonzero entries in model lattice union {1,-1}."""
    allowed = {0, 1, -1, 3, -3, 9, -9, 18, -18, 61, -61, 80, -80, 243, -243, 539, -539, 4880, -4880}
    for i in range(M.rows):
        for j in range(M.cols):
            if int(M[i, j]) not in allowed:
                return False
    return True


def axiom_det_pm1(M: sp.Matrix) -> bool:
    try:
        return abs(int(M.det())) == 1
    except Exception:
        return False


def axiom_trace0(M: sp.Matrix) -> bool:
    return int(M.trace()) == 0


def axiom_integer_charpoly_even_degree_terms_only(M: sp.Matrix) -> bool:
    """Palindromic / reciprocal candidate: only even powers (weak)."""
    chi = charpoly_matrix(M)
    pol = sp.Poly(chi, x)
    coeffs = pol.all_coeffs()
    # monic deg n: coeffs[0]=1; odd-degree terms (from x^{n-1}, x^{n-3},...) zero
    n = pol.degree()
    # coeff of x^k is at index n-k
    for k in range(n):
        if k % 2 == 1:  # odd powers
            c = pol.coeff_monomial(x**k)
            if c != 0:
                return False
    return True


def evaluate_axioms(M: sp.Matrix) -> dict:
    return {
        "companion_chain": axiom_companion_chain(M),
        "ternary_entry": axiom_ternary_entry(M),
        "model_entries": axiom_model_entries(M),
        "det_pm1": axiom_det_pm1(M),
        "trace0": axiom_trace0(M),
        "even_powers_only": axiom_integer_charpoly_even_degree_terms_only(M),
    }


def analyze_M(name: str, M: sp.Matrix) -> dict:
    chi = charpoly_matrix(M)
    rec = classify_poly(chi, do_galois=True)
    rec["name"] = name
    rec["matrix"] = [[int(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    rec["axioms"] = evaluate_axioms(M)
    return rec


def obstruction_examples() -> list[dict]:
    """Matrices that look 'structural' but fail disc^2 — Criterion 2 counterexamples."""
    examples = []
    # Base G4H M
    examples.append(analyze_M("G4H_M", G4H_M()))
    # Base T6
    examples.append(analyze_M("T6_base", T6(3, 80, 61, -3, 0, 0, 0)))
    # Pure block order-3 (reducible)
    C3 = sp.Matrix([[0, 3, 0], [0, 0, 61], [80, 0, 0]])
    F = sp.Matrix([[0, -80], [80, 0]])
    Mb = sp.zeros(5)
    Mb[0:3, 0:3] = C3
    Mb[3:5, 3:5] = F
    examples.append(analyze_M("block_diag_C3_flux", Mb))
    # Coupled but odd
    examples.append(analyze_M("T5_base_like", T5(3, 80, 61, -3, 0, 0)))
    return examples


def search_stronger_subclass():
    """
    Search T5 with extra axioms:
      - det ±1, or
      - even powers only in chi, or
      - M + M^T = 0 (skew), or
      - all entries in {0,±1,±3}
    Measure disc^2 rate among irreducible.
    """
    results = {}
    pools = {
        "tiny_ternary": [0, 1, -1, 3, -3],
        "small_model": [0, 1, -1, 3, -3, 9, 61, 80],
    }

    def gen_T5(pool):
        for a, b, c, d, e, f in itertools.product(pool, repeat=6):
            yield (a, b, c, d, e, f), T5(a, b, c, d, e, f)

    for name, pool in pools.items():
        stats = Counter()
        hits_a5 = []
        disc_sq_irr = []
        # limit
        count = 0
        max_n = 8000 if name == "tiny_ternary" else 4000
        for params, M in gen_T5(pool):
            count += 1
            if count > max_n:
                break
            # filter subclasses
            subclasses = {
                "all": True,
                "det_pm1": axiom_det_pm1(M),
                "trace0": axiom_trace0(M),
                "skew": all(int(M[i, j] + M[j, i]) == 0 for i in range(5) for j in range(5)),
                "even_powers": axiom_integer_charpoly_even_degree_terms_only(M),
            }
            chi = charpoly_matrix(M)
            pol = sp.Poly(chi, x, domain=sp.ZZ)
            if pol.LC() == -1:
                pol = sp.Poly(-pol.as_expr(), x, domain=sp.ZZ)
            if pol.degree() != 5:
                continue
            irr = bool(pol.is_irreducible)
            disc = int(pol.discriminant()) if irr else None
            sq = is_square(disc) if disc is not None else False

            for sc, ok in subclasses.items():
                if not ok:
                    continue
                key = f"{name}/{sc}"
                stats[f"{key}/tested"] += 1
                if not irr:
                    stats[f"{key}/red"] += 1
                    continue
                stats[f"{key}/irr"] += 1
                if sq:
                    stats[f"{key}/disc_sq"] += 1
                    rec = classify_poly(pol.as_expr(), do_galois=True)
                    rec["params"] = params
                    rec["subclass"] = key
                    disc_sq_irr.append(rec)
                    if rec.get("status", "").startswith("HIT_A5") or (
                        rec.get("galois") and "A5" in str(rec.get("galois"))
                    ):
                        hits_a5.append(rec)
                        stats[f"{key}/A5"] += 1

        results[name] = {
            "stats": dict(stats),
            "A5_hits": hits_a5[:20],
            "disc_sq_sample": disc_sq_irr[:30],
            "rates": {},
        }
        # compute rates
        for sc in ["all", "det_pm1", "trace0", "skew", "even_powers"]:
            key = f"{name}/{sc}"
            irr = stats.get(f"{key}/irr", 0)
            sq = stats.get(f"{key}/disc_sq", 0)
            results[name]["rates"][sc] = {
                "irr": irr,
                "disc_sq": sq,
                "disc_sq_rate_among_irr": (sq / irr) if irr else None,
                "A5": stats.get(f"{key}/A5", 0),
            }
    return results


def proposed_axioms_doc(obstructions, subclass) -> str:
    lines = [
        "# Criterion 2 — Structural axioms and the evenness obstruction",
        "",
        "## Goal",
        "",
        "Extract axioms (A1)–(Ak) from the 9 Maths such that any integer matrix \(M\)",
        "satisfying them has:",
        "",
        "1. \(\\chi_M\) irreducible over \(\\mathbb{Q}\),",
        "2. \(\\operatorname{disc}(\\chi_M)\) a square,",
        "3. Gal contains a 3-cycle (hence \(A_n\) under standard generation for small \(n\)).",
        "",
        "## Weak axiom set (current templates)",
        "",
        "| Axiom | Meaning |",
        "|-------|---------|",
        "| A1 companion_chain | Superdiagonal 1-prefix |",
        "| A2 ternary_entry | Some entry in \(\\{\\pm 3^k\\}\) |",
        "| A3 model_entries | Nonzero entries in model lattice |",
        "",
        "**These do NOT force disc².** Counterexamples below satisfy A1–A3 (or close)",
        "but have non-square discriminant.",
        "",
        "## Evenness obstruction — explicit counterexamples",
        "",
    ]
    for r in obstructions:
        lines.append(f"### `{r['name']}`")
        lines.append(f"- poly: `{r['poly']}`")
        lines.append(f"- irr={r['irreducible']} disc_sq={r['disc_square']} status={r['status']} Gal={r.get('galois')}")
        lines.append(f"- axioms: `{r['axioms']}`")
        lines.append(f"- matrix: `{r.get('matrix')}`")
        lines.append("")

    lines += [
        "## Experimental rates under stronger axioms (T5 lattice)",
        "",
        "Disc² rate among irreducibles for various subclasses:",
        "",
    ]
    for pool_name, blob in subclass.items():
        lines.append(f"### Pool `{pool_name}`")
        lines.append("")
        lines.append("| Subclass | irr | disc² | rate | A5 |")
        lines.append("|----------|----:|------:|-----:|---:|")
        for sc, rates in (blob.get("rates") or {}).items():
            rate = rates.get("disc_sq_rate_among_irr")
            rate_s = f"{rate:.4f}" if rate is not None else "—"
            lines.append(
                f"| {sc} | {rates.get('irr')} | {rates.get('disc_sq')} | {rate_s} | {rates.get('A5')} |"
            )
        lines.append("")
        if blob.get("A5_hits"):
            lines.append("Sample A5 under stronger filters:")
            for h in blob["A5_hits"][:5]:
                lines.append(f"- `{h['poly']}` ({h.get('subclass')})")
            lines.append("")

    lines += [
        "## Proposed stronger axioms (open — for theorems)",
        "",
        "| ID | Axiom (candidate) | Intended effect |",
        "|----|-------------------|-----------------|",
        "| B1 | \(M\) preserves a rational volume form (det monodromy = 1 on ambient) | evenness |",
        "| B2 | \(\\chi_M\) is a norm from a cyclic cubic étale algebra in a controlled way | 3-cycles + arithmetic |",
        "| B3 | \(M\) lies in an algebraic group with trivial sign character on eigenvalues | disc² |",
        "| B4 | Self-adjoint w.r.t. a model quadratic form with square disc | disc² |",
        "| B5 | Rigid branch-cycle type fixed by HQCC (geometric, not matrix) | full \(A_n\) |",
        "",
        "## Status",
        "",
        "- **Obstruction documented:** structural-looking matrices can be odd.",
        "- **No axiom list yet proved** to force disc² for all \(M\) in the class.",
        "- Subclass rates above show whether det±1 / skew / even-powers help empirically.",
        "",
        "_Generated by criterion2_axioms.py_",
    ]
    return "\n".join(lines)


def main():
    print("Criterion 2: axioms + evenness obstruction", flush=True)
    obs = obstruction_examples()
    for r in obs:
        print(f"  {r['name']}: disc_sq={r['disc_square']} Gal={r.get('galois')} status={r['status']}", flush=True)
    print("Searching stronger subclasses...", flush=True)
    sub = search_stronger_subclass()
    for name, blob in sub.items():
        print(f"  {name} rates: {blob['rates']}", flush=True)
    doc = proposed_axioms_doc(obs, sub)
    write_md(OUT / "CRITERION2_AXIOMS.md", doc)
    write_json(OUT / "CRITERION2_AXIOMS.json", {"obstructions": obs, "subclass_search": sub})
    # also copy pointer into results
    write_md(RESULTS / "CRITERION2_AXIOMS.md", doc)
    print(f"Wrote {OUT / 'CRITERION2_AXIOMS.md'}", flush=True)


if __name__ == "__main__":
    main()
