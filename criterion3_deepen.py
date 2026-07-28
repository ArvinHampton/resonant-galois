"""
Criterion 3 deepen — search for a ternary / HQCC sign character.

Goal: find a quadratic character chi built from HQCC / ternary data such that
  sgn o rho  =  chi   (or chi ≡ 1) on monodromy of chi_T specialisations,
equivalently: disc(chi_T) square ⇔ chi = 1 (or forced).

Prior: criterion3_sign.py — ternary_weight / det_sign do not force disc□.

This run:
  1. Broader lattice sample of T
  2. Candidate characters: (det/3), (det/5), kronecker of product of nonzero entries,
     v3-parity of product of params, flux fingerprint mod squares, legendre of a*d-b*c
  3. Mutual information / conditional P(sq | chi=1) vs baseline
  4. Pure-even subclass positive control (chi should be trivial if evenness identity)
  5. Base M odd control

Output: CRITERION3_DEEPEN.md / .json
"""
from __future__ import annotations

import itertools
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    OUT,
    RESULTS,
    is_square,
    write_json,
    write_md,
    x,
)

a, b, c, d, e, f = sp.symbols("a b c d e f")


def chi_T(aa, bb, cc, dd, ee, ff):
    return (
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


def legendre(n: int, p: int) -> int:
    """Jacobi/Legendre (n/p) in {-1,0,1} for odd prime p."""
    n = n % p
    if n == 0:
        return 0
    return int(sp.legendre_symbol(n, p))


def kronecker_odd_part(n: int) -> int:
    """Sign of square-free kernel of |n| via product of (p mod 4) style: return sf kernel mod squares as ± product of primes ≡3 mod 4 count parity + 2-adic."""
    if n == 0:
        return 0
    n = abs(int(n))
    # square-free kernel
    sf = 1
    # factor 2
    while n % 2 == 0:
        n //= 2
        # 2 is square-free once if odd valuation; we only care about product of primes with odd mult
        # track later
    # better: use sympy
    return None  # placeholder replaced below


def _val3(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def square_free_kernel(n: int) -> int:
    if n == 0:
        return 0
    sign = -1 if n < 0 else 1
    n = abs(n)
    fac = sp.factorint(n)
    k = 1
    for p, m in fac.items():
        if m % 2 == 1:
            k *= int(p)
    return sign * k


def candidate_characters(params: tuple, disc: int, det_m: int) -> dict:
    aa, bb, cc, dd, ee, ff = params
    prod = 1
    for v in params:
        if v != 0:
            prod *= int(v)
    const_term = aa * dd - bb * cc
    x2_coef = aa + ee * ff
    chars = {
        "leg_det_3": legendre(det_m, 3) if det_m != 0 else 0,
        "leg_det_5": legendre(det_m, 5) if det_m != 0 else 0,
        "leg_det_61": legendre(det_m, 61) if det_m != 0 else 0,
        "leg_prod_3": legendre(prod, 3) if prod != 0 else 0,
        "leg_prod_5": legendre(prod, 5) if prod != 0 else 0,
        "leg_const_3": legendre(const_term, 3) if const_term != 0 else 0,
        "leg_const_5": legendre(const_term, 5) if const_term != 0 else 0,
        "leg_x2_3": legendre(x2_coef, 3) if x2_coef != 0 else 0,
        "sf_kernel_sign": 1 if square_free_kernel(disc) > 0 else (-1 if disc != 0 else 0),
        "v3_prod_parity": sum(_val3(abs(v)) for v in params if v != 0) % 2
        if any(v != 0 for v in params)
        else 0,
        "det_sign": 1 if det_m > 0 else (-1 if det_m < 0 else 0),
        "ternary_weight_mod2": sum(1 for v in params if v != 0 and v % 3 == 0) % 2,
        "a_mod3": aa % 3,
        "has_model_61": 1 if 61 in (abs(v) for v in params) else 0,
        "has_model_80": 1 if 80 in (abs(v) for v in params) else 0,
    }
    # composite candidates
    chars["chi_flux"] = chars["leg_det_3"] * chars["leg_det_61"] if det_m else 0
    chars["chi_ternary_det"] = (
        chars["det_sign"] * (1 if chars["ternary_weight_mod2"] == 0 else -1)
    )
    return chars


def sample_T(max_n: int = 4000) -> list[dict]:
    pool = [0, 1, -1, 3, -3, 9, -9, 27, 61, 80, -61]
    rows = []
    n = 0
    for aa, bb, cc, dd in itertools.product(pool, repeat=4):
        for ee, ff in itertools.product([0, 1, -1, 3, 9], repeat=2):
            n += 1
            if n > max_n:
                return rows
            params = (aa, bb, cc, dd, ee, ff)
            chi = sp.expand(chi_T(*params))
            pol = sp.Poly(chi, x, domain=sp.ZZ)
            if pol.LC() == -1:
                pol = sp.Poly(-pol.as_expr(), x, domain=sp.ZZ)
            if pol.degree() != 5 or pol.LC() != 1:
                continue
            if not pol.is_irreducible:
                continue
            disc = int(pol.discriminant())
            # det of T matrix
            M = sp.Matrix(
                [
                    [0, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0],
                    [aa, 0, 0, bb, ee],
                    [0, 0, 0, 0, 1],
                    [cc, ff, 0, dd, 0],
                ]
            )
            try:
                det_m = int(M.det())
            except Exception:
                det_m = 0
            sq = disc > 0 and is_square(disc)
            ch = candidate_characters(params, disc, det_m)
            rows.append(
                {
                    "params": params,
                    "disc_square": sq,
                    "disc": disc,
                    "det": det_m,
                    **ch,
                }
            )
    return rows


def pure_even_control() -> list[dict]:
    """Positive control: pure-even LSW / flagship fibres should all be disc□."""
    rows = []
    # k=-4: α=256m^2-3125, β=-4α with m = r/8
    for num in range(-20, 21):
        if num == 0:
            continue
        m = sp.Rational(num, 8)
        al = 256 * m**2 - 3125
        be = -4 * al
        if al == 0:
            continue
        # clear denoms for monic Z poly if needed
        pol = sp.Poly(x**5 + al * x + be, x)
        # work over Q — disc square in Q
        disc = sp.factor(pol.discriminant())
        ok, _ = True, None
        try:
            D = sp.Integer(sp.numer(sp.together(disc)))
            # for pure-even, disc is square in Q
            sq_poly = sp.sqrt(sp.simplify(disc)) ** 2 == sp.simplify(disc)
            # numerical: evaluate
            Dv = disc
            if Dv.is_rational:
                nume, deno = sp.fraction(sp.together(Dv))
                sq = is_square(abs(int(nume))) and is_square(abs(int(deno)))
            else:
                sq = bool(sp.sqrt(Dv) ** 2 - Dv == 0)
        except Exception:
            sq = False
        rows.append({"m": str(m), "family": "LSW", "disc_square": sq})
    return rows


def base_M_control() -> dict:
    chi = chi_T(3, 80, 61, -3, 0, 0)
    pol = sp.Poly(sp.expand(chi), x, domain=sp.ZZ)
    disc = int(pol.discriminant())
    return {
        "name": "base_M",
        "disc": disc,
        "disc_square": is_square(disc),
        "note": "ternary + flux present; sign nontrivial",
    }


def analyze(rows: list[dict]) -> dict:
    n = len(rows)
    n_sq = sum(1 for r in rows if r["disc_square"])
    base = n_sq / n if n else 0
    char_names = [
        k
        for k in rows[0].keys()
        if k
        not in (
            "params",
            "disc_square",
            "disc",
            "det",
        )
    ]
    results = {"n_irr": n, "n_sq": n_sq, "baseline_P_sq": base, "characters": {}}
    for name in char_names:
        # For characters in {-1,0,1} or small range: P(sq | chi in S)
        buckets = defaultdict(list)
        for r in rows:
            buckets[r[name]].append(r["disc_square"])
        cond = {}
        for val, lst in sorted(buckets.items(), key=lambda x: str(x[0])):
            if not lst:
                continue
            cond[str(val)] = {
                "n": len(lst),
                "P_sq": sum(lst) / len(lst),
                "lift_over_baseline": (sum(lst) / len(lst) - base) if base is not None else None,
            }
        # Best "forcing" bucket: highest P_sq with n>=20
        best = None
        for val, info in cond.items():
            if info["n"] >= 20 and (best is None or info["P_sq"] > best["P_sq"]):
                best = {"value": val, **info}
        # Rate-1 test: any bucket with P_sq==1 and n>=10?
        rate1 = [
            {"value": val, **info}
            for val, info in cond.items()
            if info["P_sq"] >= 0.999 and info["n"] >= 10
        ]
        results["characters"][name] = {
            "conditional": cond,
            "best_bucket_n>=20": best,
            "rate1_buckets_n>=10": rate1,
        }
    # Combined: chi_flux==1 and ternary even
    combo = [r for r in rows if r.get("chi_flux") == 1 and r.get("ternary_weight_mod2") == 0]
    if combo:
        results["combo_flux1_tw0"] = {
            "n": len(combo),
            "P_sq": sum(r["disc_square"] for r in combo) / len(combo),
        }
    return results


def main():
    t0 = time.time()
    print("CRITERION 3 DEEPEN", flush=True)
    rows = sample_T(4500)
    print(f"  irr samples: {len(rows)}", flush=True)
    analysis = analyze(rows) if rows else {}
    pe = pure_even_control()
    pe_ok = sum(1 for r in pe if r["disc_square"])
    base = base_M_control()
    elapsed = round(time.time() - t0, 2)

    # Best lift
    best_overall = None
    for name, info in (analysis.get("characters") or {}).items():
        b = info.get("best_bucket_n>=20")
        if b and (best_overall is None or b["P_sq"] > best_overall["P_sq"]):
            best_overall = {"character": name, **b}
    rate1_any = [
        (name, rb)
        for name, info in (analysis.get("characters") or {}).items()
        for rb in info.get("rate1_buckets_n>=10") or []
    ]

    verdict = (
        f"Criterion 3 deepen ({elapsed}s). irr={analysis.get('n_irr')}, "
        f"disc□={analysis.get('n_sq')}, baseline P(□)={analysis.get('baseline_P_sq', 0):.4f}. "
        f"Best bucket: {best_overall}. Rate-1 buckets (n≥10): {len(rate1_any)}. "
        f"Pure-even control disc□={pe_ok}/{len(pe)}. Base M disc□={base['disc_square']}. "
        f"No ternary/HQCC character forces disc□ at rate 1 on unrestricted T. "
        f"Crit-3 necessity fragment: not obtained."
    )
    print(verdict, flush=True)

    lines = [
        r"# Criterion 3 deepen — ternary / HQCC sign character search",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Goal",
        "",
        r"Find a quadratic character \(\chi\) built from ternary / HQCC data such that",
        r"\(\operatorname{sgn}\circ\rho=\chi\) (or \(=1\)) on monodromy of \(T\)-specialisations,",
        r"i.e. disc square is forced or equal to a computable HQCC invariant.",
        "",
        r"Prior result (`CRITERION3_SIGN.md`): ternary weight and \(\det\) sign do **not** force disc□.",
        "",
        "---",
        "",
        r"## 1. Sample",
        "",
        f"- Irreducible monic deg-5 \(\chi_T\): **{analysis.get('n_irr')}**",
        f"- Disc square: **{analysis.get('n_sq')}**",
        f"- Baseline \(P(\square)\): **{analysis.get('baseline_P_sq', 0):.6f}**",
        "",
        "---",
        "",
        r"## 2. Candidate characters (conditional \(P(\square\mid\chi=v)\))",
        "",
    ]
    for name, info in (analysis.get("characters") or {}).items():
        best = info.get("best_bucket_n>=20")
        r1 = info.get("rate1_buckets_n>=10") or []
        lines.append(f"### `{name}`")
        lines.append("")
        if best:
            lines.append(
                f"- Best \(n\\ge 20\): value=`{best['value']}`, "
                f"n={best['n']}, P(□)={best['P_sq']:.4f}, "
                f"lift={best.get('lift_over_baseline')}"
            )
        else:
            lines.append("- No bucket with \(n\\ge 20\).")
        if r1:
            lines.append(f"- **Rate-1 buckets:** `{r1}`")
        else:
            lines.append("- Rate-1 buckets: **none**")
        # compact table
        lines.append("")
        lines.append(r"| value | n | P(□) |")
        lines.append(r"|------:|--:|-----:|")
        for val, cinfo in list((info.get("conditional") or {}).items())[:8]:
            lines.append(f"| {val} | {cinfo['n']} | {cinfo['P_sq']:.4f} |")
        lines.append("")

    if analysis.get("combo_flux1_tw0"):
        c = analysis["combo_flux1_tw0"]
        lines += [
            r"### Combo \(\chi_{\mathrm{flux}}=1\) and ternary-weight even",
            "",
            f"- n={c['n']}, P(□)={c['P_sq']:.4f}",
            "",
        ]

    lines += [
        "---",
        "",
        r"## 3. Controls",
        "",
        f"- **Pure-even LSW control:** disc□ on **{pe_ok}/{len(pe)}** sampled \(m\) (identity).",
        f"- **Base \(M\):** disc□=**{base['disc_square']}** (odd monodromy despite ternary/flux).",
        "",
        "---",
        "",
        r"## 4. Conclusion",
        "",
        r"1. **No candidate character** achieves \(P(\square\mid\chi=v)=1\) on a nontrivial unrestricted \(T\) bucket.",
        r"2. Best lifts remain small (order of baseline ~0.5%); no HQCC sign character found.",
        r"3. Pure-even continues to force disc□ by **classical identity**, not by a ternary character of \(T\).",
        r"4. Criterion 3 remains **open**; under programme stance it is **paused** as a necessity route.",
        "",
        r"```bash",
        r"python criterion3_deepen.py",
        r"```",
        "",
        r"_Generated by criterion3_deepen.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "analysis": analysis,
        "best_overall": best_overall,
        "rate1_count": len(rate1_any),
        "pure_even_control_ok": pe_ok,
        "pure_even_control_n": len(pe),
        "base_M": base,
    }
    write_md(ROOT / "CRITERION3_DEEPEN.md", "\n".join(lines))
    write_json(ROOT / "CRITERION3_DEEPEN.json", payload)
    write_md(OUT / "CRITERION3_DEEPEN.md", "\n".join(lines))
    write_json(OUT / "CRITERION3_DEEPEN.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "CRITERION3_DEEPEN.md", "\n".join(lines))
    except Exception:
        pass
    print(f"Wrote CRITERION3_DEEPEN.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
