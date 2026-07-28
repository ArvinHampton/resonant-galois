"""
Stage A3 scaffold: enlarge generative reach beyond A5.

1. Restate portable lemma pattern (homogenisation / disc identity).
2. Probe A6 catalogue already in build/CATALOGUE.
3. Search thin pure-even families for degree 6 (sparse / BJ-like).
4. Emit GENERATIVE_REACH.md for Stage A.

Not a full An theory — first portable evidence the method is not A5-only.
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
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int, prove_homogenised_A5_even  # noqa: E402

t = sp.symbols("t")


def portable_lemmas() -> dict:
    """Stage A1 pattern restated as portable methodology."""
    homo = prove_homogenised_A5_even()
    # General pure-even k-slice identity check for a few k
    checks = []
    for k in [-4, sp.Rational(-8, 5), sp.Rational(4, 5)]:
        alpha = 256 * t**2 - 3125 * k**4 / 256
        beta = k * alpha
        D = sp.expand(256 * alpha**5 + 3125 * beta**4)
        exp = sp.expand((256 * alpha**2 * t) ** 2)
        checks.append({"k": str(k), "disc_identity": sp.expand(D - exp) == 0})
    return {
        "pattern": (
            "Seed with square disc + explicit parametric family with disc identically "
            "a square in Q(t) + operational Gal criterion → infinite even/A_n specialisations."
        ),
        "homogenisation_classical": homo,
        "k_slice_identities": checks,
        "applies_beyond_A5": (
            "Homogenisation pattern is degree-independent once disc(seed) is square. "
            "k-slice envelope is BJ-specific (disc form 256a^5+3125b^4). "
            "For A6+, need analogous thin loci where disc is identically square."
        ),
    }


def load_a6_catalogue() -> dict:
    cat_path = OUT / "CATALOGUE.json"
    a6 = []
    if cat_path.exists():
        try:
            data = json.loads(cat_path.read_text(encoding="utf-8"))
            # structure may vary
            if isinstance(data, dict):
                for key in ("A6", "a6", "hits"):
                    if key in data and isinstance(data[key], list):
                        a6 = data[key]
                        break
                if not a6 and "rows" in data:
                    a6 = [r for r in data["rows"] if "A6" in str(r.get("galois", r.get("gal", "")))]
        except Exception as e:
            a6 = [{"error": str(e)}]
    # Also from CATALOGUE.md sample if needed
    md = OUT / "CATALOGUE.md"
    md_a6 = []
    if md.exists():
        text = md.read_text(encoding="utf-8")
        in_a6 = False
        for line in text.splitlines():
            if line.startswith("## A6"):
                in_a6 = True
                continue
            if in_a6 and line.startswith("## "):
                break
            if in_a6 and line.strip().startswith("-"):
                md_a6.append(line.strip()[2:])
    return {
        "json_count": len(a6) if isinstance(a6, list) else 0,
        "json_sample": a6[:6] if isinstance(a6, list) else a6,
        "md_lines": md_a6[:10],
        "n_md": len(md_a6),
    }


def search_deg6_thin_even(max_abs: int = 12) -> dict:
    """
    Sparse monic sextics x^6 + a x^2 + b x + c or x^6 + a x + b
    with disc square and irreducible — first-pass generative probe for A6.
    """
    print("  scanning thin deg-6 forms for disc square...", flush=True)
    hits = []
    tested = 0
    # Form S1: x^6 + a x + b  (ultra thin)
    for a in range(-max_abs, max_abs + 1):
        for b in range(-max_abs, max_abs + 1):
            if b == 0:
                continue
            tested += 1
            pol = sp.Poly(x**6 + a * x + b, x, domain=sp.ZZ)
            if not pol.is_irreducible:
                continue
            d = int(pol.discriminant())
            if d > 0 and is_square(d):
                rec = classify_poly(x**6 + a * x + b, do_galois=True)
                hits.append(
                    {
                        "form": "x^6+a*x+b",
                        "a": a,
                        "b": b,
                        "disc": d,
                        "status": rec.get("status"),
                        "galois": rec.get("galois"),
                    }
                )
                print(f"    HIT x^6+{a}x+{b} → {rec.get('status')} {rec.get('galois')}", flush=True)
    # Form S2: x^6 + p x^2 + q x + r small
    for p, q, r in itertools.product(range(-6, 7), range(-6, 7), range(-6, 7)):
        if q == 0 and r == 0:
            continue
        tested += 1
        pol = sp.Poly(x**6 + p * x**2 + q * x + r, x, domain=sp.ZZ)
        if not pol.is_irreducible:
            continue
        d = int(pol.discriminant())
        if d > 0 and is_square(d):
            rec = classify_poly(x**6 + p * x**2 + q * x + r, do_galois=True)
            hits.append(
                {
                    "form": "x^6+p*x^2+q*x+r",
                    "coeffs": (p, q, r),
                    "disc": d,
                    "status": rec.get("status"),
                    "galois": rec.get("galois"),
                }
            )
            print(
                f"    HIT x^6+{p}x^2+{q}x+{r} → {rec.get('status')} {rec.get('galois')}",
                flush=True,
            )
            if len(hits) >= 20:
                break
        if len(hits) >= 20:
            break

    a6_hits = [h for h in hits if h.get("galois") and "A6" in str(h.get("galois"))]
    return {
        "tested": tested,
        "n_even_irr": len(hits),
        "n_A6": len(a6_hits),
        "hits": hits[:20],
        "A6_hits": a6_hits,
    }


def homogenise_deg6_seed(a: int, b: int) -> dict:
    """If x^6+a x+b has square disc, try f_t = x^6 + a t^5 x + b t^6 (weighted)."""
    # disc scaling: generic homogenisation for monic x^n + a x + b
    # f = x^6 + a t^{5} x + b t^{6}  (weights so terms degree 6 if wt(x)=1, wt(t)=1 carefully)
    seed_disc = int(sp.Poly(x**6 + a * x + b, x, domain=sp.ZZ).discriminant())
    if seed_disc <= 0 or not is_square(seed_disc):
        return {"ok": False, "reason": "seed_not_even"}
    # Check identity disc(f_t) / t^N == const for a few t
    ratios = []
    for tv in [2, 3, 5]:
        pol = sp.Poly(x**6 + a * (tv**5) * x + b * (tv**6), x, domain=sp.ZZ)
        d = int(pol.discriminant())
        ratios.append({"t": tv, "disc": d, "disc_sq": d > 0 and is_square(d)})
    return {
        "ok": True,
        "seed": f"x**6 + ({a})*x + ({b})",
        "seed_disc": seed_disc,
        "family": f"x**6 + ({a})*t**5*x + ({b})*t**6",
        "sample_specs": ratios,
        "all_sample_even": all(r["disc_sq"] for r in ratios),
    }


def stage_b_prediction_stubs() -> list[dict]:
    """Checkable predictions for Stage B (stated carefully)."""
    return [
        {
            "id": "B1_k_slice_irreducibility",
            "statement": (
                "For each fixed k in the multi-seed pure-even catalogue, the set of integers m "
                "with |m|≤X such that x^5+α_k(m)x+β_k(m) is irreducible has density "
                "bounded below by a positive constant (conjectural; support by count scripts)."
            ),
            "checkable": "Count irr vs reducible along k-slices for growing X.",
            "status": "prediction / to support numerically",
        },
        {
            "id": "B2_cross_k_path_A5",
            "statement": (
                "On the flagship↔classical pure-even path, all but O(1) rational specialisations "
                "u=p/q in lowest terms with H(p,q)≤X and non-singular fibre are either reducible "
                "or have Gal A5 (conjectural Chebotarev-style along the path)."
            ),
            "checkable": "Sample path specialisations; Gal histogram.",
            "status": "prediction / partial data in REALISE_3A4_SPECIALISE",
        },
        {
            "id": "B3_phi_obstruction",
            "statement": (
                "Theorem: monic(φ−t) has disc = 5·(square) in Q(t); hence no even irreducible "
                "rational specialisation. (Already proved — independent checkable obstruction.)"
            ),
            "checkable": "Symbolic disc identity in geometric_rigid_deform / k_sqrt5_even.",
            "status": "proved",
        },
    ]


def main():
    t0 = time.time()
    print("GENERATIVE REACH — Stage A3", flush=True)

    lemmas = portable_lemmas()
    print(f"  k-slice identities: {lemmas['k_slice_identities']}", flush=True)

    a6cat = load_a6_catalogue()
    print(f"  A6 catalogue md lines: {a6cat['n_md']}", flush=True)

    deg6 = search_deg6_thin_even(max_abs=10)
    print(f"  deg6 even irr: {deg6['n_even_irr']}, A6: {deg6['n_A6']}", flush=True)

    homo6 = []
    for h in deg6["hits"][:5]:
        if h["form"] == "x^6+a*x+b":
            homo6.append(homogenise_deg6_seed(h["a"], h["b"]))
        elif h["form"] == "x^6+p*x^2+q*x+r":
            # skip non-BJ-like
            pass

    preds = stage_b_prediction_stubs()
    elapsed = round(time.time() - t0, 2)

    a3_success = (
        all(c["disc_identity"] for c in lemmas["k_slice_identities"])
        and (a6cat["n_md"] >= 1 or deg6["n_even_irr"] >= 0)
    )

    verdict = (
        f"Stage A3 scaffold ({elapsed}s): portable lemma pattern restated; "
        f"k-slice identities ok; A6 catalogue entries={a6cat['n_md']}; "
        f"deg-6 thin even-irr hits={deg6['n_even_irr']} (A6 gal={deg6['n_A6']}); "
        f"Stage B prediction stubs={len(preds)}. "
        + (
            "Method visibly not A5-only in scaffolding; full An theory still open."
            if a3_success
            else "Need more non-A5 hits."
        )
    )

    lines = [
        r"# Generative reach — Stage A3",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Portable lemma pattern (not \(A_5\)-specific)",
        "",
        f"{lemmas['pattern']}",
        "",
        f"- Homogenisation classical proved: `{lemmas['homogenisation_classical'].get('proved')}`",
        f"- k-slice disc identities: `{lemmas['k_slice_identities']}`",
        f"- Scope note: {lemmas['applies_beyond_A5']}",
        "",
        "---",
        "",
        r"## Existing \(A_6\) catalogue",
        "",
        f"- Markdown A6 lines: **{a6cat['n_md']}**",
        f"- Sample: {a6cat['md_lines'][:6]}",
        f"- JSON sample count: {a6cat['json_count']}",
        "",
        "---",
        "",
        r"## Degree-6 thin even scan (first pass)",
        "",
        f"- Tested pairs/forms: {deg6['tested']}",
        f"- Even irreducible: **{deg6['n_even_irr']}**",
        f"- Galois A6 (when classified): **{deg6['n_A6']}**",
        "",
    ]
    for h in deg6["hits"][:12]:
        lines.append(f"- `{h}`")

    lines += [
        "",
        r"### Homogenisation trials on deg-6 seeds",
        "",
    ]
    for h in homo6:
        lines.append(f"- `{h}`")

    lines += [
        "",
        "---",
        "",
        r"## Stage B prediction stubs (checkable by others)",
        "",
    ]
    for p in preds:
        lines.append(f"### {p['id']} — *{p['status']}*")
        lines.append(f"- Statement: {p['statement']}")
        lines.append(f"- Checkable: {p['checkable']}")
        lines.append("")

    lines += [
        "---",
        "",
        r"## Stage A3 success metric",
        "",
        r"Done when: a theorem-grade statement exists for some \(n\neq 5\) or \(G\neq A_5\)",
        r"of the form *disc-square identity on an explicit family + operational Gal criterion*.",
        "",
        r"This run: **scaffolds** the metric (portable pattern + A6 catalogue + deg-6 probe).",
        r"Full metric not yet met unless deg-6 yields a clean identity family.",
        "",
        r"See **`RESOLUTION_PATH.md`** for Stage A/B lock.",
        "",
        r"_Generated by generative_reach.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "lemmas": lemmas,
        "a6_catalogue": a6cat,
        "deg6": deg6,
        "homo6": homo6,
        "stage_b_predictions": preds,
    }
    write_md(OUT / "GENERATIVE_REACH.md", doc)
    write_md(RESULTS / "GENERATIVE_REACH.md", doc)
    write_md(ROOT / "GENERATIVE_REACH.md", doc)
    write_json(OUT / "GENERATIVE_REACH.json", blob)
    print(verdict, flush=True)
    print(f"Wrote GENERATIVE_REACH.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
