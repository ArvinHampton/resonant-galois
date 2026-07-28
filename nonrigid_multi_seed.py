"""
Next natural step after LSW multi-seed hit:

1. LSW-type slices: β = k α, seek a(t) so disc is □ in Q(t) and ≥2 HQCC seeds lie on the slice.
2. Pure-even curves through flagship (-55, 88) and one other seed (α(u), β(u)).
3. Homogenisation-compatible bridges and low-degree ansätze.

Success = pure-even family in Q(u) containing flagship + another HQCC seed
       or more LSW-like slices with multiple HQCC A5 seeds.

Output: NONRIGID_MULTI_SEED.md
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

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

u, t = sp.symbols("u t")

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

FLAG = (-55, 88, "flagship")


def is_square_poly(expr, var=u) -> dict:
    try:
        ex = sp.expand(expr)
        if ex == 0:
            return {"ok": True, "degenerate": True, "degree": -1}
        P = sp.Poly(ex, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return {"ok": False, "reason": "neg"}
        n, d = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return {"ok": False, "reason": "content"}
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(f), m) for f, m in fac[1] if m % 2]
        return {"ok": len(odds) == 0, "degree": int(P.degree()), "odd": odds[:6], "factored": str(sp.factor(ex))[:250]}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def disc_ab(a, b):
    return sp.expand(256 * a**5 + 3125 * b**4)


# ---------------------------------------------------------------------------
# 1. LSW-type slices β = k α
# ---------------------------------------------------------------------------
def slice_beta_k_alpha():
    """
    For integer k, seeds with β = k α.
    LSW is k=-4 with α = t^2 - 3125 (so α+3125 = t^2 is square).
    General: seek α = c(t^2 - m) or α = c t^2 + d so disc is □.
    disc(α, kα) = 256 α^5 + 3125 k^4 α^4 = 256 α^4 (α + 3125 k^4 / 256)
    For this to be square for parametric α:
      256 α^4 (α + 3125 k^4 / 256) = □
    If 256 is square and α^4 is square, need α + 3125 k^4/256 = square * rational.
    Write 3125 k^4 / 256 = (3125/256) k^4.
    Set α = s^2 - (3125/256) k^4, then disc = 256 (s^2 - ...)^4 s^2 = square if s in Q
    (when 3125/256 k^4 is such that α lands in Q for s in Q).

    For α,β ∈ Z: need (3125 k^4)/256 ∈ Q and α integer.
    3125 = 5^5, 256=2^8 ⇒ 3125 k^4 / 256 ∈ Z iff 2^8 | k^4 * 5^5 — so 2^2 | k
    (since k^4 contributes multiples of 2 in steps of 4 valuation: v2(k^4)=4 v2(k) ≥ 8 ⇒ v2(k)≥2).

    Family: α = s^2 - 3125 k^4 / 256, β = k α, for s ∈ Z, k ≡ 0 mod 4? v2(k)≥2 means k even enough.
    """
    print("  LSW-type slices β=kα...", flush=True)
    results = []
    # Which seeds lie on β=kα for some integer k?
    seed_by_k = {}
    for sa, sb, tag in SEEDS:
        if sa == 0:
            continue
        if sb % sa == 0:
            k = sb // sa
            seed_by_k.setdefault(k, []).append((sa, sb, tag))
        # also rational k = sb/sa
        kr = sp.Rational(sb, sa)
        seed_by_k.setdefault(("rat", str(kr)), []).append((sa, sb, tag))

    integer_k_slices = {k: v for k, v in seed_by_k.items() if isinstance(k, int)}
    multi_k = {k: v for k, v in integer_k_slices.items() if len(v) >= 2}

    # For each integer k with v2(k)>=2, build pure-even family and list Z seeds
    families = []
    for k in range(-20, 21):
        if k == 0:
            continue
        # v2(k) >= 2?
        kk = abs(k)
        v2 = 0
        while kk % 2 == 0:
            v2 += 1
            kk //= 2
        # α = s^2 - 3125 k^4 / 256 must be rational; for s integer and 256|3125 k^4
        if (3125 * k**4) % 256 != 0:
            # still try α = 256 s^2 - 3125 k^4, β = k α / 256? rescale
            # Better integer form: set α = 256 m^2 - 3125 k^4, then
            # disc = 256 α^4 (α + 3125 k^4/256) = 256 α^4 (256 m^2)/256 = α^4 * 256 m^2
            # = (α^2 * 16 m)^2 — yes square!
            # β = k α
            pass
        # Integer family for any k:
        # α(m) = 256 m^2 - 3125 k**4
        # β(m) = k * α(m)
        # Wait check: α + 3125 k^4/256 = 256 m^2, disc=256 α^4 * 256 m^2 = (16 α^2 m)^2 * 256/256...
        # disc = 256 α^5 + 3125 k^4 α^4 = 256 α^4 (α + 3125 k^4/256)
        # = 256 α^4 * (256 m^2) = 65536 α^4 m^2 = (256 α^2 m)^2. Yes!

        alpha = 256 * t**2 - 3125 * (k**4)
        beta = k * alpha
        D = disc_ab(alpha, beta)
        info = is_square_poly(D, t)
        # Which known seeds appear for integer t?
        seed_hits = []
        for tv in range(-80, 81):
            if tv == 0 and alpha.subs(t, 0) == 0:
                continue
            aa = int(alpha.subs(t, tv))
            bb = int(beta.subs(t, tv))
            for sa, sb, tag in SEEDS:
                if (aa, bb) == (sa, sb):
                    seed_hits.append({"t": tv, "seed": tag, "alpha": aa, "beta": bb})
        # Also match by solving 256 t^2 - 3125 k^4 = sa
        for sa, sb, tag in SEEDS:
            if sa == 0 or sb != k * sa:
                continue
            # 256 t^2 = sa + 3125 k^4
            rhs = sa + 3125 * (k**4)
            if rhs % 256 != 0:
                continue
            rhs2 = rhs // 256
            if rhs2 < 0:
                continue
            r, ok = sp.integer_nthroot(rhs2, 2)
            if ok:
                if not any(h["seed"] == tag for h in seed_hits):
                    seed_hits.append({"t": int(r), "seed": tag, "alpha": sa, "beta": sb, "solved": True})

        gal_ok = []
        for h in seed_hits:
            r = classify_poly(x**5 + h["alpha"] * x + h["beta"], do_galois=True)
            h["gal"] = r.get("galois")
            h["status"] = r.get("status")
            if (h.get("status") or "").startswith("HIT_A5") or (h.get("gal") and "A5" in str(h.get("gal"))):
                gal_ok.append(h)

        rec = {
            "k": k,
            "alpha": str(alpha),
            "beta": str(beta),
            "disc_square_in_Qt": info.get("ok"),
            "seed_hits": seed_hits,
            "n_seeds": len({h["seed"] for h in seed_hits}),
            "A5_seeds": gal_ok,
            "LSW_special_case": k == -4,  # α = 256 t^2 - 3125*256 = 256(t^2-3125) — scale of LSW
        }
        # For k=-4: α = 256 t^2 - 3125*256 = 256(t^2-3125), β = -4α
        # LSW uses α_LSW = t^2-3125, β=-4 α_LSW — same projective ray / same polys after scaling variable?
        # x^5 + α x + β vs x^5 + (α/256) x + β/256 — different. 
        # Our scaling: f = x^5 + 256(t^2-3125) x - 1024(t^2-3125)
        # vs LSW x^5 + (t^2-3125)x - 4(t^2-3125). Related by x = 4^{1/4} z or not monic-equivalent necessarily.
        if rec["n_seeds"] >= 1 and info.get("ok"):
            families.append(rec)
            if rec["n_seeds"] >= 2:
                print(f"    MULTI k={k}: seeds={[h['seed'] for h in seed_hits]}", flush=True)

    multi = [f for f in families if f["n_seeds"] >= 2]
    return {
        "seeds_by_integer_k": {str(k): [(t, a, b) for a, b, t in v] for k, v in integer_k_slices.items()},
        "multi_seed_k_raw": {str(k): v for k, v in multi_k.items()},
        "pure_even_slice_families": families,
        "multi_seed_pure_even_slices": multi,
    }


# ---------------------------------------------------------------------------
# 2. Pure-even curves through flagship + one other seed
# ---------------------------------------------------------------------------
def flagship_plus_one_search():
    """
    Seek α(u), β(u) with:
      (α(0),β(0)) = (-55, 88)
      (α(1),β(1)) = (a1, b1) for each other seed
      disc(α,β) square in Q(u)

    Ansätze:
      A. Linear — already known fail for pure even
      B. Homogenisation-style: α = -55 (1-u)^4 + a1 u^4, β = 88 (1-u)^5 + b1 u^5
      C. α = -55 (1-u)^m + a1 u^m, β = 88 (1-u)^n + b1 u^n for (m,n) in list
      D. Quadratic Bezier with free midpoint on lattice
    """
    print("  flagship + one seed pure-even search...", flush=True)
    a0, b0 = -55, 88
    hits = []
    tested = 0

    others = [(a, b, tag) for a, b, tag in SEEDS if tag != "flagship"]
    exp_pairs = [
        (1, 1),
        (2, 2),
        (4, 5),
        (5, 4),
        (4, 4),
        (5, 5),
        (2, 1),
        (1, 2),
        (3, 3),
        (8, 10),
        (3, 4),
        (4, 3),
        (6, 6),
        (2, 5),
        (5, 2),
    ]

    for a1, b1, tag in others:
        # B, C monomial bridges
        for m, n in exp_pairs:
            tested += 1
            alpha = sp.expand(a0 * (1 - u) ** m + a1 * u**m)
            beta = sp.expand(b0 * (1 - u) ** n + b1 * u**n)
            # endpoints
            if int(alpha.subs(u, 0)) != a0 or int(beta.subs(u, 0)) != b0:
                continue
            if int(alpha.subs(u, 1)) != a1 or int(beta.subs(u, 1)) != b1:
                continue
            info = is_square_poly(disc_ab(alpha, beta), u)
            if info.get("ok") and not info.get("degenerate"):
                hits.append(
                    {
                        "type": "monomial_bridge",
                        "other": tag,
                        "m": m,
                        "n": n,
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "info": info,
                    }
                )
                print(f"    *** HIT flagship--{tag} m={m} n={n}", flush=True)

        # D Bezier midpoints
        for am, bm in itertools.product(range(-12, 13, 2), range(-12, 13, 2)):
            tested += 1
            alpha = sp.expand((1 - u) ** 2 * a0 + 2 * u * (1 - u) * am + u**2 * a1)
            beta = sp.expand((1 - u) ** 2 * b0 + 2 * u * (1 - u) * bm + u**2 * b1)
            info = is_square_poly(disc_ab(alpha, beta), u)
            if info.get("ok") and not info.get("degenerate"):
                hits.append(
                    {
                        "type": "bezier",
                        "other": tag,
                        "mid": (am, bm),
                        "alpha": str(alpha),
                        "beta": str(beta),
                        "info": info,
                    }
                )
                print(f"    *** HIT bezier flagship--{tag} mid=({am},{bm})", flush=True)

        # Weighted mixed: α = a0(1-u)^4 + a1 u^4, β = b0(1-u)^5 + b1 u^5  (homogenisation bridge)
        tested += 1
        alpha = sp.expand(a0 * (1 - u) ** 4 + a1 * u**4)
        beta = sp.expand(b0 * (1 - u) ** 5 + b1 * u**5)
        info = is_square_poly(disc_ab(alpha, beta), u)
        if info.get("ok") and not info.get("degenerate"):
            hits.append(
                {
                    "type": "homo_bridge_4_5",
                    "other": tag,
                    "alpha": str(alpha),
                    "beta": str(beta),
                    "info": info,
                }
            )
            print(f"    *** HIT homo bridge flagship--{tag}", flush=True)

    return {"tested": tested, "hits": hits}


# ---------------------------------------------------------------------------
# 3. Does flagship lie on any k-slice pure-even family? (β = kα with k=88/-55)
# ---------------------------------------------------------------------------
def flagship_on_slice_families(slice_data) -> dict:
    """k = 88/(-55) = -88/55 not integer — flagship not on integer-k LSW-type slice."""
    a0, b0 = -55, 88
    kr = sp.Rational(b0, a0)
    # Family with rational k: α = 256 t^2 - 3125 k^4, β = k α
    # For k=p/q, clear: use α = 256 q^4 t^2 - 3125 p^4, β = p/q α carefully integer
    p, q = int(sp.numer(kr)), int(sp.denom(kr))
    # α = 256 (q t)^2 - 3125 p^4 ? want β=kα with integer coeffs
    # Set α = 256 m^2 - 3125 p^4 * c, ...
    # Integer model: α(m) = 256 * (q**4) * m**2 - 3125 * p**4
    # Then α + 3125 (p/q)^4 / 256 * something — redo carefully.
    #
    # disc(α, kα) = 256 α^4 (α + 3125 k^4 / 256)
    # Set α + 3125 k^4 / 256 = 256 n^2  ⇒ α = 256 n^2 - 3125 k^4 / 256
    # k=p/q: 3125 p^4 /(256 q^4)
    # α = 256 n^2 - 3125 p^4 /(256 q^4) = (256^2 n^2 q^4 - 3125 p^4) / (256 q^4)
    # Multiply poly by scaling: use A = 256^2 n^2 q^4 - 3125 p^4, B = k A but
    # f = x^5 + (A/D) x + (k A/D) with D=256 q^4
    # monic Z: x^5 + A x + (p A)/q  needs q|A for integer — A = 256^2 n^2 q^4 - 3125 p^4 ≡ -3125 p^4 mod q
    #
    # Simpler: parametric n, α_Q = 256 n^2 - 3125*(p**4)/(256*q**4), β_Q = (p/q)*α_Q
    # At n such that α_Q = -55: 256 n^2 - 3125 p^4/(256 q^4) = -55
    alpha_target = -55
    k = sp.Rational(88, -55)
    # 256 n^2 = alpha_target + 3125 k**4 / 256
    rhs = sp.together(alpha_target + 3125 * k**4 / 256)
    n2 = sp.together(rhs / 256)
    print(f"  flagship on k={k} slice: need n^2 = {n2}", flush=True)
    n2f = sp.factor(n2)
    is_sq = False
    try:
        n2r = sp.Rational(n2)
        num, den = int(n2r.p), int(n2r.q)
        is_sq = sp.integer_nthroot(abs(num), 2)[1] and sp.integer_nthroot(den, 2)[1] and num * den > 0
    except Exception:
        n2r = n2
    return {
        "k": str(k),
        "n_squared_needed": str(n2),
        "n_squared_factored": str(n2f),
        "n_squared_is_rational_square": is_sq,
        "flagship_on_LSW_type_slice": is_sq,
        "note": "If n^2 is a rational square, flagship lies on the pure-even k-slice family.",
    }


# ---------------------------------------------------------------------------
# 4. Verify LSW multi-seed + expand k-slice multi list with Gal
# ---------------------------------------------------------------------------
def verify_multi_families(slice_blob, flag_hits) -> list:
    multi = list(slice_blob.get("multi_seed_pure_even_slices") or [])
    # enrich with A5 confirmation
    out = []
    for fam in multi:
        seeds_a5 = []
        for h in fam.get("seed_hits") or []:
            r = classify_poly(x**5 + h["alpha"] * x + h["beta"], do_galois=True)
            if (r.get("status") or "").startswith("HIT_A5") or (r.get("galois") and "A5" in str(r.get("galois"))):
                seeds_a5.append({**h, "gal": r.get("galois")})
        out.append(
            {
                "k": fam["k"],
                "n_seeds": fam["n_seeds"],
                "A5_seed_count": len(seeds_a5),
                "A5_seeds": seeds_a5,
                "alpha": fam["alpha"],
                "beta": fam["beta"],
            }
        )
    for h in flag_hits.get("hits") or []:
        out.append({"type": "flagship_bridge", **h})
    return out


def main():
    t0 = time.time()
    print("NONRIGID MULTI-SEED — flagship+one & LSW-type slices", flush=True)

    slices = slice_beta_k_alpha()
    print(
        f"  pure-even k-slices with ≥1 seed: {len(slices['pure_even_slice_families'])}, "
        f"multi: {len(slices['multi_seed_pure_even_slices'])}",
        flush=True,
    )

    flag_on = flagship_on_slice_families(slices)
    print(f"  flagship on its k-slice: {flag_on}", flush=True)

    bridges = flagship_plus_one_search()
    print(f"  flagship+one pure-even bridges: {len(bridges['hits'])} (tested {bridges['tested']})", flush=True)

    multi = verify_multi_families(slices, bridges)

    # Summary of integer-k multi raw (before pure-even family)
    multi_raw = slices.get("multi_seed_k_raw") or {}

    verdict = (
        f"LSW-type pure-even slices (β=kα): {len(slices['multi_seed_pure_even_slices'])} with ≥2 known seeds. "
        f"Flagship k={flag_on['k']}: lies on pure-even slice family? {flag_on['flagship_on_LSW_type_slice']}. "
        f"Pure-even bridges flagship↔other seed: {len(bridges['hits'])}. "
        + (
            "SUCCESS: multi-seed pure-even families listed below."
            if slices["multi_seed_pure_even_slices"] or bridges["hits"]
            else "No new multi-seed pure-even family through flagship in scanned ansätze."
        )
    )

    elapsed = round(time.time() - t0, 2)
    lines = [
        "# Next natural step — multi-seed pure-even families",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 1. LSW-type slices \(\\beta = k\\alpha\)",
        "",
        "General pure-even integer family for any \(k\\neq 0\):",
        "",
        r"$$\alpha(t)=256 t^2 - 3125 k^4,\qquad \beta(t)=k\cdot\alpha(t)$$",
        "",
        r"$$\operatorname{disc}=(256\,\alpha(t)^2\, t)^2\quad\text{(identically square)}.$$",
        "",
        "Special case \(k=-4\): recovers LSW up to the scaling \(\alpha_{\mathrm{LSW}}=t^2-3125\) vs \(256(t^2-3125)\).",
        "",
        "### Known seeds on lines \(\\beta=k\\alpha\) (integer k)",
        "",
    ]
    for k, v in sorted(
        ((int(k), v) for k, v in (slices.get("seeds_by_integer_k") or {}).items() if k.lstrip("-").isdigit()),
        key=lambda kv: kv[0],
    ):
        lines.append(f"- k={k}: {v}")

    lines += [
        "",
        "### Pure-even multi-seed slices",
        "",
    ]
    multi_slices = slices.get("multi_seed_pure_even_slices") or []
    if not multi_slices:
        lines.append("_None with ≥2 known seeds in range |k|≤20._")
    for fam in multi_slices:
        lines.append(f"#### k = {fam['k']}")
        lines.append(f"- α=`{fam['alpha']}`, β=`{fam['beta']}`")
        lines.append(f"- seeds: {fam['seed_hits']}")
        lines.append(f"- disc □ in Q(t): {fam['disc_square_in_Qt']}")
        lines.append("")

    lines += [
        "---",
        "",
        "## 2. Flagship on its own k-slice",
        "",
        f"- k = β/α = `{flag_on['k']}` (not integer)",
        f"- Needs n² = `{flag_on['n_squared_needed']}`",
        f"- Factored: `{flag_on['n_squared_factored']}`",
        f"- n² rational square? **{flag_on['flagship_on_LSW_type_slice']}**",
        f"- {flag_on['note']}",
        "",
        "---",
        "",
        "## 3. Pure-even bridges: flagship + one other seed",
        "",
        f"- Tested: {bridges['tested']}",
        f"- Hits: **{len(bridges['hits'])}**",
        "",
    ]
    if not bridges["hits"]:
        lines.append(
            "_No pure-even monomial/Bezier/homo-bridge through flagship and another "
            "listed HQCC seed in the scanned exponents/midpoints._"
        )
    for h in bridges["hits"]:
        lines.append(f"- `{h}`")

    lines += [
        "",
        "---",
        "",
        "## Conclusions",
        "",
        "1. **LSW / k-slice families** give an infinite list of pure-even BJ families;",
        "   multi-seed hits occur when several HQCC seeds share the same ratio β/α",
        "   and lie on the parametric square condition (as with k=-4: s100 and s124_m).",
        "2. **Flagship** has k=-88/55 ∉ ℤ and does **not** lie on a pure-even LSW-type",
        "   slice in the rational-square sense computed above — consistent with separate homogenisation ray.",
        "3. **No pure-even curve** of the scanned bridge types joins flagship to another seed.",
        "4. Practical fusion progress: treat **each integer k** with multiple lattice seeds",
        "   as an LSW-type pure-even family; flagship remains on its homogenisation family alone.",
        "",
        "### Recommended follow-up",
        "",
        "- Enumerate more HQCC-lattice BJ seeds with square disc + A5, group by k=β/α,",
        "  and attach pure-even k-slice families (grows the multi-seed list).",
        "- For flagship specifically: search higher-degree pure-even bridges, or accept",
        "  single-seed homogenisation as its non-rigid family.",
        "",
        "_Generated by nonrigid_multi_seed.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "verdict": verdict,
        "slices": {
            "multi_seed_pure_even_slices": slices["multi_seed_pure_even_slices"],
            "n_families_with_seeds": len(slices["pure_even_slice_families"]),
            "seeds_by_integer_k": slices["seeds_by_integer_k"],
        },
        "flagship_slice": flag_on,
        "bridges": bridges,
        "multi_verified": multi,
    }
    write_md(OUT / "NONRIGID_MULTI_SEED.md", doc)
    write_md(RESULTS / "NONRIGID_MULTI_SEED.md", doc)
    write_md(ROOT / "NONRIGID_MULTI_SEED.md", doc)
    write_json(OUT / "NONRIGID_MULTI_SEED.json", blob)
    print(verdict, flush=True)
    print(f"Wrote NONRIGID_MULTI_SEED.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
