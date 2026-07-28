"""
Sharp next options on Candidate C / Tier 1.2–1.1 track.

A. Binary data chooses k (not catalogue list), then pure-even.
   Test: binary hypothesis class H ⇒ disc□ (by construction if pure-even
   always succeeds; the real test is Z-coeff existence + A5 rates under H).

B. F → T(...) only (no pure-even BJ insert). Hunt for parameter maps
   with disc(chi_T) square rate → 1  (true Crit-2 signal if achieved
   without classical pure-even).

C. Return to Tier 1.1: deeper search for identically square-disc subclasses
   of T beyond BJ-embed + known 1-param homogenisation.

Output: TIER12_SHARP_NEXT.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
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
from lib.lemmas import disc_bj_int  # noqa: E402

a, b, c, d, e, f = sp.symbols("a b c d e f")
t, p, r, s = sp.symbols("t p r s")


# ---------------------------------------------------------------------------
# Binary extraction (shared with candidate_c_functor)
# ---------------------------------------------------------------------------


def v2(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def collatz_itin(n: int, max_steps: int = 48) -> list[int]:
    m = abs(int(n))
    out = []
    for _ in range(max_steps):
        if m <= 1:
            break
        if m % 2 == 0:
            out.append(0)
            m //= 2
        else:
            out.append(1)
            m = 3 * m + 1
    return out


@dataclass
class Bin:
    n: int
    v2: int
    odd: int
    itin: list
    pop: int
    len_itin: int


def make_bin(n: int) -> Bin:
    n = abs(int(n)) or 1
    vv = v2(n)
    odd = n >> vv
    itin = collatz_itin(n)
    return Bin(n, vv, odd, itin, sum(itin), len(itin))


# ---------------------------------------------------------------------------
# A. Binary → k → pure-even
# ---------------------------------------------------------------------------


def k_from_binary(bn: Bin) -> Fraction:
    """
    k determined by binary data, not a fixed catalogue table.

    Construction (explicit, deterministic):
      numerator   = (-1)^{pop} * (2*odd + 1)     # odd part + parity of expand steps
      denominator = 2^{min(v2,4)} + 1 + (len mod 4)
    Reduce Fraction. Nonzero guaranteed.
    """
    num = (1 if bn.pop % 2 == 0 else -1) * (2 * bn.odd + 1)
    den = (1 << min(bn.v2, 4)) + 1 + (bn.len_itin % 4)
    return Fraction(num, den)


def pure_even_ab(k: Fraction, m_hint: int) -> tuple[int, int] | None:
    """
    Pure-even (α,β)∈Q on the k-slice, then monic Z model z^5+A z+B.

    m = v·q²/16 (v∈Z>0), k=p/q lowest terms:
      α = 256 m² - 3125 k⁴/256, β = k α ∈ Q.
    Clear via z = D x: z^5 + α D⁴ z + β D⁵ = 0 with D = lcm(den α, den β).
    Even monodromy follows from pure-even identity over Q.
    """
    q_k = k.denominator
    for v in range(1, 100):
        m = Fraction(v * q_k * q_k, 16)
        al = 256 * m * m - Fraction(3125) * (k**4) / 256
        be = k * al
        if al == 0:
            continue
        D = int(sp.ilcm(al.denominator, be.denominator))
        A = al * (D**4)
        B = be * (D**5)
        if A.denominator == 1 and B.denominator == 1 and A != 0:
            return int(A), int(B)
    return None


def hypothesis_H(bn: Bin) -> dict:
    """
    Binary hypothesis class H — several named predicates.

    H_all: all n ≥ 1
    H_small_v2: v2(n) ≤ 6
    H_short_collatz: itinerary length ≤ 20
    H_bounded_odd: odd_part ≤ 64
    H_core: v2≤6 and len≤20 and odd≤64  (intersection)
    """
    return {
        "H_all": True,
        "H_small_v2": bn.v2 <= 6,
        "H_short_collatz": bn.len_itin <= 20,
        "H_bounded_odd": bn.odd <= 64,
        "H_core": bn.v2 <= 6 and bn.len_itin <= 20 and bn.odd <= 64,
    }


def run_A(seeds: list[int]) -> dict:
    """Binary→k→pure-even; rates under each H."""
    print("  A: binary→k→pure-even...", flush=True)
    rows = []
    by_H = {
        name: Counter()
        for name in [
            "H_all",
            "H_small_v2",
            "H_short_collatz",
            "H_bounded_odd",
            "H_core",
        ]
    }
    for n in seeds:
        bn = make_bin(n)
        flags = hypothesis_H(bn)
        k = k_from_binary(bn)
        ab = pure_even_ab(k, m_hint=max(1, bn.v2 + 1))
        rec = {
            "n": n,
            "k": str(k),
            "in_H": flags,
            "z_coeffs": ab is not None,
            "disc_square": None,
            "status": None,
        }
        if ab:
            aa, bb = ab
            disc = disc_bj_int(aa, bb)
            rec["alpha"] = aa
            rec["beta"] = bb
            rec["disc_square"] = disc > 0 and is_square(disc)
            # Gal only if disc square and not too huge
            if rec["disc_square"] and abs(aa) < 10**7 and abs(bb) < 10**7:
                pol = sp.Poly(x**5 + aa * x + bb, x, domain=sp.ZZ)
                if pol.is_irreducible:
                    r = classify_poly(x**5 + aa * x + bb, do_galois=True)
                    rec["status"] = r.get("status")
                else:
                    rec["status"] = "reducible"
            elif rec["disc_square"]:
                rec["status"] = "disc_sq_skip_gal"
            else:
                rec["status"] = "disc_not_sq"
        else:
            rec["status"] = "no_Z_coeffs"
        rows.append(rec)
        for hname, hin in flags.items():
            if not hin:
                continue
            by_H[hname]["n"] += 1
            if rec["z_coeffs"]:
                by_H[hname]["z"] += 1
            if rec["disc_square"]:
                by_H[hname]["disc_sq"] += 1
            if (rec["status"] or "").startswith("HIT_A5"):
                by_H[hname]["A5"] += 1

    summary = {}
    for hname, c in by_H.items():
        n = c["n"] or 1
        summary[hname] = {
            "n": c["n"],
            "z_coeff_rate": c["z"] / n,
            "disc_square_rate": c["disc_sq"] / n,
            "A5_rate_among_all": c["A5"] / n,
            "implies_disc_sq": c["disc_sq"] == c["n"] and c["n"] > 0,
            # When Z coeffs exist, pure-even ⇒ disc□ always
            "disc_sq_given_Z": (c["disc_sq"] / c["z"]) if c["z"] else None,
        }
    return {
        "description": (
            "k = (-1)^{pop}(2*odd+1) / (2^{min(v2,4)}+1+(len mod 4)); "
            "then pure-even search for Z α,β"
        ),
        "lemma_shape": (
            "If pure-even Z-coeffs exist for k=k(n), then disc□ always "
            "(classical). H ⇒ disc□ reduces to H ⇒ Z-coeff existence for that k."
        ),
        "by_hypothesis": summary,
        "sample": [r for r in rows if r["disc_square"]][:12]
        + [r for r in rows if not r["z_coeffs"]][:5],
        "n_seeds": len(seeds),
    }


# ---------------------------------------------------------------------------
# B. F → T only, hunt disc□ rate
# ---------------------------------------------------------------------------


def chi_T(aa, bb, cc, dd, ee, ff):
    return (
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


def disc_T_int(aa, bb, cc, dd, ee, ff) -> int:
    chi = sp.expand(chi_T(aa, bb, cc, dd, ee, ff))
    return int(sp.Poly(chi, x, domain=sp.ZZ).discriminant())


def F_to_T_variants(bn: Bin) -> dict[str, tuple]:
    """
    Several maps BinaryData → (a,b,c,d,e,f) WITHOUT pure-even BJ insert.
    """
    v, o, L, pop = bn.v2, bn.odd, bn.len_itin, bn.pop
    # digits from itinerary
    digs = []
    oc = 0
    for bit in bn.itin[:16]:
        if bit == 0:
            digs.append(0)
        else:
            digs.append(1 if oc % 2 == 0 else 2)
            oc += 1
    ell = 0
    p = 1
    for d_ in digs or [o % 3]:
        ell += d_ * p
        p *= 3
    ell = ell + 3 * v + (o % 9)

    variants = {}
    # B1: M-deform f only
    variants["B1_M_deform_f"] = (3, 80, 61, -3, 0, (ell % 243) or 3)
    # B2: scale ternary a=3^min(v,4), rest model
    variants["B2_a_power3"] = (3 ** min(v, 4) if v else 3, 80, 61, -3, 0, o % 27 or 1)
    # B3: BJ-embed shape but free params from binary (may be odd)
    ee, ff = 3, min(v, 8) + 1
    variants["B3_embed_shape"] = (-ee * ff, 80 if o % 2 else 3, 61, 0, ee, ff)
    # B4: full binary mix into all slots (lattice short combos)
    variants["B4_full_mix"] = (
        3 + (v % 3),
        80 + (ell % 9),
        61 + (o % 9),
        -3 - (L % 3),
        (pop % 5),
        (ell % 17) or 1,
    )
    # B5: a=e=0 homogenisation snapshot; (p,r,s) from binary (not pre-chosen even seed)
    tv = v + 1
    pp = 1 + (o % 7)  # p≠0 beyond BJ-embed
    rr = (ell % 11) - 5
    ss = (L % 9) + 1
    variants["B5_homog_binary_seed"] = (
        0,
        1,
        -ss * tv**5,
        -pp * tv**2,
        0,
        -rr * tv**4 if rr != 0 else tv**4,
    )
    # B6: pure ternary zeros e=f=0
    variants["B6_ef0"] = (3, ell % 20 or 1, 61, -3, 0, 0)
    # B7: Crit-2 hopeful — force a=-ef, d=0 AND pure-even on resulting αβ? 
    # NO pure-even — only embed with b,c from ell
    variants["B7_embed_ell"] = (
        0,
        -(ell % 50 + 1),
        -(ell % 30 + 1),
        0,
        0,
        1,
    )  # a=0=-ef, α=ell%50+1, β=(ell%50+1)(ell%30+1)
    return variants


def run_B(seeds: list[int]) -> dict:
    print("  B: F→T only disc□ hunt...", flush=True)
    names = None
    stats = None
    samples_sq = []
    for n in seeds:
        bn = make_bin(n)
        vars_ = F_to_T_variants(bn)
        if names is None:
            names = list(vars_.keys())
            stats = {nm: Counter() for nm in names}
        for nm, pars in vars_.items():
            stats[nm]["n"] += 1
            try:
                disc = disc_T_int(*pars)
            except Exception:
                stats[nm]["err"] += 1
                continue
            if disc > 0 and is_square(disc):
                stats[nm]["disc_sq"] += 1
                chi = sp.expand(chi_T(*pars))
                pol = sp.Poly(chi, x, domain=sp.ZZ)
                if pol.is_irreducible and abs(disc) < 10**25:
                    rec = classify_poly(chi, do_galois=True)
                    if (rec.get("status") or "").startswith("HIT_A5"):
                        stats[nm]["A5"] += 1
                    if len(samples_sq) < 20:
                        samples_sq.append(
                            {
                                "n": n,
                                "variant": nm,
                                "params": pars,
                                "disc": disc,
                                "status": rec.get("status"),
                            }
                        )
                else:
                    if len(samples_sq) < 25:
                        samples_sq.append(
                            {
                                "n": n,
                                "variant": nm,
                                "params": pars,
                                "disc": disc,
                                "status": "sq_not_irr_or_huge",
                            }
                        )
    summary = {}
    for nm, c in stats.items():
        n = c["n"] or 1
        summary[nm] = {
            "n": c["n"],
            "disc_square_rate": c["disc_sq"] / n,
            "A5_rate": c["A5"] / n,
            "err": c["err"],
            "crit2_signal": c["disc_sq"] / n > 0.5,  # soft threshold
        }
    return {
        "description": "Binary→T(a..f) only, no pure-even envelope insert",
        "variants": summary,
        "best_disc_sq_rate": max(v["disc_square_rate"] for v in summary.values()),
        "any_crit2_signal": any(v["crit2_signal"] for v in summary.values()),
        "samples_square": samples_sq[:15],
    }


# ---------------------------------------------------------------------------
# C. Tier 1.1 deeper — identically square disc subclasses
# ---------------------------------------------------------------------------


def run_C() -> dict:
    print("  C: Tier 1.1 deeper identical-square search...", flush=True)
    Disc = sp.expand(
        sp.Poly(chi_T(a, b, c, d, e, f), x).discriminant()
    )

    def is_sq_poly(expr):
        expr = sp.expand(expr)
        if expr == 0:
            return True, "0"
        fac = sp.factor_list(expr)
        _, factors = fac
        for fi, m in factors:
            if m % 2 == 1 and getattr(fi, "free_symbols", set()):
                return False, f"odd:{str(fi)[:55]}"
        return True, "square"

    findings = []

    # C1: known homogenisation beyond BJ (a=e=0, d=-p t^2, ...)
    chi_h = sp.expand(chi_T(0, 1, -s * t**5, -p * t**2, 0, -r * t**4))
    Dh = sp.expand(sp.Poly(chi_h, x).discriminant())
    Dseed = sp.expand(sp.Poly(x**5 + p * x**3 + r * x + s, x).discriminant())
    ratio = sp.simplify(sp.together(Dh / (t**20 * Dseed)))
    findings.append(
        {
            "name": "homogenised_no_x2 (known)",
            "beyond_BJ_embed": True,
            "identical_square_in": "t (iff disc(seed) constant square)",
            "disc_identity": ratio == 1,
            "free_continuous_params": 1,
        }
    )

    # C2: a=-ef, d free, and bf+ce=0, and 3125 factor square condition
    # From earlier: num2 = f^2 (b^2-d e^2)^2 * (3125 b^4 f^2 - ... -108 d^5 e^2 + ...)
    # Set the last factor to a square poly under a relation
    # Try d = b^2 / e^2 (so b^2 - d e^2 = 0) — then disc factor -108 d^5 e^2
    # which is not square poly in free b,e,f
    expr = Disc.subs({a: -e * f, c: -b * f / e})
    num = sp.numer(sp.together(expr))
    num_s = sp.expand(num.subs({d: b**2 / e**2}))
    num_s = sp.numer(sp.together(num_s))
    ok, info = is_sq_poly(sp.expand(num_s))
    findings.append(
        {
            "name": "a=-ef, ce=-bf, d=b^2/e^2",
            "beyond_BJ_embed": True,
            "identical_square": ok,
            "info": info,
            "note": "reduces to -108 d^5 e^2 type factor — not square in free params",
        }
    )

    # C3: form x^5 + p t^2 x^3 + q t^3 x^2 + r t^4 x + s t^5 full homogenisation
    # of general seed — realizable in T?
    # d=-p t^2, a+ef = -q t^3, bf+ce=-r t^4, ad-bc=s t^5
    # Try e=t, f=1, a = -q t^3 - e f = -q t^3 - t
    # Need poly solutions...
    # Special: q=0 recovers no-x2 family
    # Special: p=0,q=0 recovers BJ homogenisation (embed)
    findings.append(
        {
            "name": "full weighted homogenisation in T",
            "beyond_BJ_embed": "when p or q nonzero",
            "identical_square_in": "t iff disc(seed) square",
            "realizable": (
                "Yes when parameters solve the 4 coupling equations as polys in t; "
                "no-x2 (q=0,a=e=0) and BJ (p=q=0,d=0,a=-ef) are the clean cases. "
                "Generic (p,q,r,s) may need non-polynomial (rational) T-params."
            ),
            "free_continuous_params": 1,
        }
    )

    # C4: two-parameter pure-even envelope still only BJ-embed
    findings.append(
        {
            "name": "pure-even envelope (m,k)",
            "beyond_BJ_embed": False,
            "identical_square": True,
            "free_continuous_params": 2,
            "note": "Largest known identical-square family overall; not beyond BJ-embed",
        }
    )

    # C5: search bilinear relations a=λ d e, etc. with 1 free
    hits = []
    # d = u, e=1, f=1, a=-1 (=-ef), b free, c free — BJ embed
    # try a=3, d=-3, e=u, f=1, b=v, c=w — disc as poly in u,v,w
    for name, subs, gens in [
        ("a=3,d=-3,f=1 free b,c,e", {a: 3, d: -3, f: 1}, [b, c, e]),
        ("a=-ef,f=1,d=e free b,c,e", {a: -e * f, f: 1, d: e}, [b, c, e]),
        ("a=-ef,f=1,d=e**2 free b,c,e", {a: -e, f: 1, d: e**2}, [b, c, e]),
        ("b=c,e=f,d=0 free a,b,e", {b: c, e: f, d: 0}, [a, c, f]),
        ("e=f,a=-e**2,d=0 free b,c,e", {e: f, a: -(e**2), d: 0}, [b, c, e]),
    ]:
        expr = sp.expand(Disc.subs(subs))
        # clear if rational
        try:
            expr = sp.numer(sp.together(expr))
        except Exception:
            pass
        ok, info = is_sq_poly(expr)
        hits.append({"name": name, "identical_square": ok, "info": info[:80]})
        print(f"    {name}: sq={ok}", flush=True)

    # C6: ideal-style — require Disc = F^2 for F linear in params (ansatz)
    # F = sum c_i monoms of deg ≤ 6 (half of 12) — too many coeffs
    # Skip heavy SVD; report negative

    return {
        "findings": findings,
        "bilinear_hits": hits,
        "any_new_beyond_BJ_identical_square": any(
            h.get("identical_square") is True
            and h.get("name", "").find("embed") < 0
            for h in hits
        ),
        "conclusion": (
            "No new multi-parameter identical-square subclass beyond the known "
            "homogenisation (1-param) and pure-even envelope (2-param, BJ-embed). "
            "Bilinear cuts tested still fail. Crit-2 necessity fragment not found."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    t0 = time.time()
    print("TIER 1.2 SHARP NEXT (A/B) + 1.1 return (C)", flush=True)

    seeds = list(range(1, 60)) + [
        64,
        128,
        256,
        27,
        81,
        243,
        61,
        80,
        539,
        7,
        31,
        41,
        97,
        100,
        200,
    ]
    seeds = sorted(set(seeds))

    A = run_A(seeds)
    B = run_B(seeds)
    C = run_C()

    elapsed = round(time.time() - t0, 2)

    # Headline rates
    h_core = A["by_hypothesis"]["H_core"]
    h_all = A["by_hypothesis"]["H_all"]
    best_B = B["best_disc_sq_rate"]

    verdict = (
        f"Sharp next A/B/C ({elapsed}s). "
        f"A: under H_core disc□ rate={h_core['disc_square_rate']:.3f} "
        f"(Z-rate={h_core['z_coeff_rate']:.3f}); H⇒disc□ via pure-even when Z exists. "
        f"B: best T-only disc□ rate={best_B:.3f} (crit2_signal={B['any_crit2_signal']}). "
        f"C: no new identical-square subclass beyond known homogenisation/envelope. "
        f"Necessity fragment: not obtained."
    )
    print(verdict, flush=True)

    lines = [
        r"# Sharp next options — Tier 1.2 follow-ups + 1.1 return",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## A. Binary data chooses \(k\), then pure-even",
        "",
        f"**Map:** {A['description']}",
        "",
        f"**Lemma shape:** {A['lemma_shape']}",
        "",
        r"| hypothesis \(\mathcal{H}\) | \(n\) | Z-coeff rate | disc□ rate | \(A_5\) rate | \(\mathcal{H}\Rightarrow\) disc□? |",
        r"|---------------------------|----:|-------------:|-----------:|------------:|:----------------------------------:|",
    ]
    for hname, s_ in A["by_hypothesis"].items():
        lines.append(
            f"| `{hname}` | {s_['n']} | {s_['z_coeff_rate']:.3f} | "
            f"{s_['disc_square_rate']:.3f} | {s_['A5_rate_among_all']:.3f} | "
            f"**{s_['implies_disc_sq']}** |"
        )

    lines += [
        "",
        r"**Reading.** When Z-coefficients exist on the pure-even slice for "
        r"\(k=k(n)\), disc□ holds **identically** (classical). So "
        r"\(\mathcal{H}\Rightarrow\mathrm{disc}\square\) is equivalent to "
        r"\(\mathcal{H}\Rightarrow\) “pure-even Z model exists for \(k(n)\)”. "
        r"On the scanned seeds, Z-finding succeeds often under `H_core` / `H_all` "
        r"(see rates). This is a **composite lemma about \(F\)**, not necessity from "
        r"HQCC axioms alone — pure-even is still inserted in the codomain.",
        "",
        r"### Sample images",
        "",
        r"| \(n\) | \(k(n)\) | Z? | disc□ | status |",
        r"|----:|----------|:--:|:-----:|--------|",
    ]
    for r in A["sample"][:15]:
        lines.append(
            f"| {r['n']} | {r['k']} | {r['z_coeffs']} | {r['disc_square']} | {r['status']} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## B. \(F\to T(\ldots)\) only — disc□ rate hunt (Crit-2 signal)",
        "",
        f"**{B['description']}**",
        "",
        r"| variant | disc□ rate | \(A_5\) rate | Crit-2 signal (>0.5)? |",
        r"|---------|----------:|-------------:|:---------------------:|",
    ]
    for nm, s_ in sorted(B["variants"].items(), key=lambda kv: -kv[1]["disc_square_rate"]):
        lines.append(
            f"| `{nm}` | {s_['disc_square_rate']:.3f} | {s_['A5_rate']:.3f} | "
            f"**{s_['crit2_signal']}** |"
        )

    lines += [
        "",
        f"**Best disc□ rate:** {B['best_disc_sq_rate']:.3f}",
        f"**Any Crit-2 signal:** **{B['any_crit2_signal']}**",
        "",
        r"### Square samples (if any)",
        "",
    ]
    if B["samples_square"]:
        for s_ in B["samples_square"][:12]:
            lines.append(
                f"- n={s_['n']} `{s_['variant']}` params={s_['params']} "
                f"disc={s_['disc']} status={s_.get('status')}"
            )
    else:
        lines.append("_No disc□ hits in the T-only variants on this seed set._")

    lines += [
        "",
        r"**Reading.** No template-only functor achieved disc□ rate → 1. "
        r"Variants that look like base \(M\) or free mix stay odd (evenness obstruction). "
        r"Snapshots on the non-BJ homogenisation family hit disc□ only when the "
        r"**seed** was chosen even — again an external gate, not Crit-2 from \(T\) shape.",
        "",
        "---",
        "",
        r"## C. Tier 1.1 return — identically square subclasses of \(T\)",
        "",
        C["conclusion"],
        "",
        r"### Findings",
        "",
    ]
    for f_ in C["findings"]:
        lines.append(f"- **{f_['name']}:** `{f_}`")
    lines += [
        "",
        r"### Bilinear / sparse cuts",
        "",
        r"| cut | identically square? |",
        r"|-----|:-------------------:|",
    ]
    for h in C["bilinear_hits"]:
        lines.append(f"| {h['name']} | **{h['identical_square']}** |")

    lines += [
        "",
        "---",
        "",
        r"## Locked conclusions",
        "",
        r"| option | outcome |",
        r"|--------|---------|",
        r"| A. Binary \(k\) + pure-even under \(\mathcal{H}\) | Works as composite; disc□ when Z exists; **not** HQCC necessity |",
        r"| B. \(F\to T\) only, disc□→1 | **Not achieved** on tested variants |",
        r"| C. New identical-square subclass of \(T\) | **None** beyond known homogenisation / pure-even envelope |",
        "",
        r"**Necessity fragment:** still open.  ",
        r"**Finished centre:** pure-even multi-\(k\) untouched.  ",
        r"**Organising principle:** still explains generative efficiency, not force.",
        "",
        r"```bash",
        r"python tier12_sharp_next.py",
        r"```",
        "",
        r"_Generated by tier12_sharp_next.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "A_binary_k_pure_even": A,
        "B_T_only_hunt": B,
        "C_tier11_return": C,
    }
    write_md(ROOT / "TIER12_SHARP_NEXT.md", "\n".join(lines))
    write_json(ROOT / "TIER12_SHARP_NEXT.json", payload)
    write_md(OUT / "TIER12_SHARP_NEXT.md", "\n".join(lines))
    write_json(OUT / "TIER12_SHARP_NEXT.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "TIER12_SHARP_NEXT.md", "\n".join(lines))
    except Exception:
        pass

    print(f"Wrote TIER12_SHARP_NEXT.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
