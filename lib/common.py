"""Shared Galois / model utilities for the resonant Galois programme."""
from __future__ import annotations

import json
from collections import Counter
from math import isqrt
from pathlib import Path

import sympy as sp

ROOT = Path(r"C:\Users\bradl\Desktop\resonant_galois")
RESULTS = Path(r"C:\Users\bradl\Desktop\a5_brute_results")
OUT = ROOT / "build"
OUT.mkdir(parents=True, exist_ok=True)

x = sp.symbols("x")

MODEL_CORE = {
    3: "ternary/generations",
    9: "3^2",
    18: "visible_digits",
    61: "punctures",
    80: "flux/61",
    243: "3^5_towers",
    520: "tower_res",
    539: "period",
    4880: "flux_budget",
}


def is_square(n: int) -> bool:
    if n is None or n <= 0:
        return False
    return bool(sp.integer_nthroot(int(n), 2)[1])


def monic_poly(expr) -> sp.Poly | None:
    pol = sp.Poly(sp.expand(expr), x, domain=sp.ZZ)
    if pol.LC() == -1:
        pol = sp.Poly(-pol.as_expr(), x, domain=sp.ZZ)
    if pol.LC() != 1:
        return None
    return pol


def classify_poly(expr, do_galois: bool = True) -> dict:
    pol = monic_poly(expr)
    if pol is None:
        return {"status": "not_monic_Z", "poly": str(expr)}
    rec = {
        "poly": str(pol.as_expr()),
        "coeffs": [int(c) for c in pol.all_coeffs()],
        "degree": pol.degree(),
        "irreducible": bool(pol.is_irreducible),
        "discriminant": None,
        "disc_square": None,
        "galois": None,
        "status": None,
        "census": None,
    }
    if not rec["irreducible"]:
        rec["status"] = "reducible"
        rec["factor_Q"] = str(sp.factor(pol.as_expr()))
        return rec
    disc = int(pol.discriminant())
    rec["discriminant"] = disc
    rec["disc_square"] = is_square(disc)
    rec["census"] = cycle_census(pol)
    if not rec["disc_square"]:
        rec["status"] = "odd_monodromy"
        if do_galois and pol.degree() <= 7:
            try:
                g, alt = pol.galois_group(by_name=True)
                rec["galois"] = str(g)
                rec["galois_alt"] = bool(alt)
            except Exception as e:
                rec["galois_error"] = str(e)
        return rec
    if do_galois and pol.degree() <= 7:
        try:
            g, alt = pol.galois_group(by_name=True)
            rec["galois"] = str(g)
            rec["galois_alt"] = bool(alt)
            short = str(g).split(".")[-1]
            if short.startswith("A") and short[1:].isdigit():
                rec["status"] = f"HIT_{short}"
            else:
                rec["status"] = f"even:{short}"
        except Exception as e:
            rec["status"] = f"sq_disc_gal_err:{type(e).__name__}"
            rec["galois_error"] = str(e)
    else:
        rec["status"] = "sq_disc"
    return rec


def cycle_census(poly: sp.Poly, max_p: int = 40) -> dict:
    counts: Counter = Counter()
    used = 0
    disc = int(poly.discriminant())
    for p in sp.primerange(2, 300):
        if used >= max_p:
            break
        if disc % p == 0:
            continue
        try:
            facs = sp.factor_list(poly.as_expr(), modulus=int(p))
            degs = []
            for f, m in facs[1]:
                degs.extend([int(sp.degree(f))] * int(m))
            counts[tuple(sorted(degs))] += 1
            used += 1
        except Exception:
            continue
    return {
        "primes_used": used,
        "patterns": {str(k): v for k, v in sorted(counts.items(), key=lambda kv: -kv[1])},
        "has_3": any(3 in k for k in counts),
        "has_type_3111": any(sorted(k) == [1, 1, 1, 3] or sorted(k) == [1, 1, 3] for k in counts),
        "has_type_33": any(sorted(k) == [3, 3] for k in counts),
        "has_5": any(5 in k for k in counts),
    }


def charpoly_matrix(M: sp.Matrix) -> sp.Expr:
    return sp.expand(M.charpoly(x).as_expr())


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


def write_md(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
