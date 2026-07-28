"""
Tier 1.2 — Explicit functor for Candidate C (mod-2 ↔ mod-3 design mirror).

  F : BinaryData → TernaryOutput

Precise types, three concrete realizations F1–F3, monodromy tests on images
(disc□ rate, Gal A5 rate vs random controls of comparable height).

Output: CANDIDATE_C_FUNCTOR.md / .json
"""
from __future__ import annotations

import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    MODEL_CORE,
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int  # noqa: E402

# ---------------------------------------------------------------------------
# Input / output types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BinaryData:
    """
    Binary height / 2-adic branched data.

    Fields
    ------
    n : positive integer seed (the primary binary object)
    v2 : 2-adic valuation v_2(n)  [contract depth at start]
    odd_part : n / 2^{v2}  [expand residue class odd]
    collatz_itinerary : finite word in {0,1}*  (0=even/contract, 1=odd/expand)
        Classical Collatz step: even → n/2 (label 0); odd → 3n+1 (label 1) then continue.
        Truncated at max_steps or first entry ≤ 1.
    collatz_length : len(itinerary)
    """

    n: int
    v2: int
    odd_part: int
    collatz_itinerary: tuple  # of 0/1
    collatz_length: int


@dataclass(frozen=True)
class TernaryLatticeElement:
    """Integer built from order-3 data (ternary digits / model generators)."""

    value: int
    ternary_digits: tuple  # digits in {0,1,2}, least significant first
    construction: str


@dataclass(frozen=True)
class TemplateParams:
    """Parameters (a,b,c,d,e,f) for structural template T."""

    a: int
    b: int
    c: int
    d: int
    e: int
    f: int
    construction: str


@dataclass(frozen=True)
class BJPair:
    """Bring–Jerrard coefficients."""

    alpha: int
    beta: int
    k: str  # rational as string
    construction: str


@dataclass(frozen=True)
class TernaryOutput:
    """Full image of F: lattice + optional T params + optional BJ pair."""

    lattice: TernaryLatticeElement
    template: TemplateParams | None
    bj: BJPair | None
    source_n: int
    functor_id: str


# ---------------------------------------------------------------------------
# Binary extraction
# ---------------------------------------------------------------------------


def v2(n: int) -> int:
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def extract_binary(n: int, max_steps: int = 40) -> BinaryData:
    n0 = abs(int(n))
    if n0 == 0:
        n0 = 1
    vv = v2(n0)
    odd = n0 >> vv
    # Collatz itinerary from n0
    m = n0
    itin = []
    for _ in range(max_steps):
        if m <= 1:
            break
        if m % 2 == 0:
            itin.append(0)
            m //= 2
        else:
            itin.append(1)
            m = 3 * m + 1
    return BinaryData(
        n=n0,
        v2=vv,
        odd_part=odd,
        collatz_itinerary=tuple(itin),
        collatz_length=len(itin),
    )


# ---------------------------------------------------------------------------
# Ternary maps (T3-style)
# ---------------------------------------------------------------------------

MODEL_LIST = sorted(MODEL_CORE.keys())


def binary_word_to_ternary_digits(itin: tuple, mode: str = "01_to_01_alt2") -> list[int]:
    """
    Map binary itinerary → ternary digits.

    Modes
    -----
    01_to_01_alt2 : 0→0 (contract), 1→1 or 2 alternating by odd-step index
    01_to_012_run : 0→0, 1→1, and double-1 runs use 2
    """
    digits = []
    odd_count = 0
    prev_one = False
    for bit in itin:
        if bit == 0:
            digits.append(0)
            prev_one = False
        else:
            if mode == "01_to_01_alt2":
                digits.append(1 if odd_count % 2 == 0 else 2)
                odd_count += 1
            elif mode == "01_to_012_run":
                if prev_one:
                    digits.append(2)
                else:
                    digits.append(1)
                prev_one = True
            else:
                digits.append(1)
    return digits


def digits_to_int(digits: list[int], base: int = 3) -> int:
    v = 0
    p = 1
    for d in digits:
        v += int(d) * p
        p *= base
    return v


def mix_with_model(val: int, model_idx: int) -> int:
    """Fold in a model generator (order-3 / flux lattice element)."""
    m = MODEL_LIST[model_idx % len(MODEL_LIST)]
    # short combination: val + m or val * 3 + m
    return int(val) + int(m)


# ---------------------------------------------------------------------------
# Functors F1, F2, F3
# ---------------------------------------------------------------------------


def F1_lattice(bin_data: BinaryData) -> TernaryOutput:
    """
    F1 — itinerary evaluation.

    Binary itinerary → ternary digits → integer ℓ ∈ N.
    Template: T(3, 80, 61, -3, 0, ℓ mod 243) sparse deformation of base M.
    BJ: pure-even on k=-4 (LSW) with m related to v2 and ℓ.
    """
    digits = binary_word_to_ternary_digits(bin_data.collatz_itinerary, "01_to_01_alt2")
    if not digits:
        digits = [bin_data.odd_part % 3]
    ell = digits_to_int(digits, 3)
    ell = mix_with_model(ell, bin_data.v2)
    lat = TernaryLatticeElement(
        value=ell,
        ternary_digits=tuple(digits),
        construction="F1: Collatz itinerary→ternary digits→base3 + model[v2]",
    )
    # Template deformation of M: f = ell mod 243 (bounded), keep ternary 3,80,61
    f_par = ell % 243
    if f_par == 0:
        f_par = 3
    tmpl = TemplateParams(
        a=3, b=80, c=61, d=-3, e=0, f=int(f_par), construction="F1: M-deform f=ℓ mod 243"
    )
    # BJ pure-even LSW k=-4: α = m^2 - 3125 after scale... use α=256m^2-3125*256, wait
    # k=-4: α = 256 m^2 - 3125*256 = 256(m^2 - 3125), β = -4α
    # Use m = 3^{v2} + (ell % 16) for variety, m≠0
    m = 3 ** min(bin_data.v2, 8) + (ell % 16)
    if m == 0:
        m = 1
    # integer α,β: use α = m*m - 3125, β = -4*α  (LSW scale) when m^2>3125
    # standard LSW: x^5+(t^2-3125)x-4(t^2-3125)
    alpha = m * m - 3125
    beta = -4 * alpha
    bj = BJPair(
        alpha=int(alpha),
        beta=int(beta),
        k="-4",
        construction=f"F1: LSW pure-even m={m} from 3^v2+(ℓ mod 16)",
    )
    return TernaryOutput(
        lattice=lat, template=tmpl, bj=bj, source_n=bin_data.n, functor_id="F1"
    )


def F2_height_map(bin_data: BinaryData) -> TernaryOutput:
    """
    F2 — valuation transport v2 → v3-style height.

    Map (v2, odd_part) → ternary height h = v2, residue rho = odd_part mod 3^w.
    Lattice: ℓ = 3^h * q + model, q from odd_part.
    Template: BJ-embed style a=-e f, d=0 with e,f,b,c from binary data.
    """
    h = min(bin_data.v2, 10)
    q = bin_data.odd_part % (3**5)
    ell = (3**h) * (q if q else 1) + 61  # puncture generator
    digits = []
    tmp = ell
    for _ in range(12):
        digits.append(tmp % 3)
        tmp //= 3
        if tmp == 0:
            break
    lat = TernaryLatticeElement(
        value=ell,
        ternary_digits=tuple(digits),
        construction="F2: ℓ=3^{v2}*(odd mod 3^5)+61",
    )
    # BJ-embed params: e = 3, f = h+1, a = -e*f = -3(h+1), d=0
    # b = 80 or from odd, c = 61
    ee, ff = 3, h + 1
    aa = -ee * ff
    bb = 80 if bin_data.odd_part % 2 == 1 else 0
    cc = 61
    dd = 0
    tmpl = TemplateParams(
        a=aa, b=bb, c=cc, d=dd, e=ee, f=ff, construction="F2: BJ-embed a=-ef,d=0 from v2"
    )
    # chi BJ: α=-(bf+ce), β=-bc
    alpha = -(bb * ff + cc * ee)
    beta = -bb * cc
    k_str = str(Fraction(beta, alpha)) if alpha else "inf"
    bj = BJPair(
        alpha=int(alpha),
        beta=int(beta),
        k=k_str,
        construction="F2: BJ-embed chi coefficients",
    )
    return TernaryOutput(
        lattice=lat, template=tmpl, bj=bj, source_n=bin_data.n, functor_id="F2"
    )


def F3_pure_even_k_from_itinerary(bin_data: BinaryData) -> TernaryOutput:
    """
    F3 — itinerary selects multi-seed ratio k, pure-even envelope in m.

    k chosen from catalogue by popcount(itinerary) mod #catalogue.
    m from v2 and odd_part. Forces disc□ by construction on BJ image.
    """
    catalogue_k = [
        Fraction(-4),
        Fraction(4),
        Fraction(-8, 5),
        Fraction(8, 5),
        Fraction(4, 5),
        Fraction(-4, 5),
        Fraction(-12, 5),
        Fraction(-16, 5),
    ]
    pop = sum(bin_data.collatz_itinerary) if bin_data.collatz_itinerary else bin_data.odd_part
    kv = catalogue_k[pop % len(catalogue_k)]
    # m rational so α,β ∈ Z: use m = q^2 * t / 16 style — simple: m = Fraction(5,16)*s
    s_int = max(1, bin_data.v2 + (bin_data.odd_part % 7))
    # For general k=p/q, α=256m^2-3125 k^4/256 ∈ Z for suitable m
    # Use m = Fraction(kv.denominator**2 * s_int, 16) wait
    # Practical: use integer parameter s on cleared form
    # α = 256*(den^2 * s)^2 wait simpler LSW-style only for k=-4, else:
    # m = Fraction(s_int * kv.denominator**2, 16)
    m = Fraction(s_int * (kv.denominator**2), 16)
    alpha = 256 * (m**2) - Fraction(3125) * (kv**4) / 256
    beta = kv * alpha
    # clear denominators of alpha,beta for Z poly? may be rational
    # scale: if alpha = A/D, use poly x^5 + A x + B with common denom cleared via x=u
    A, B = alpha, beta
    if A.denominator != 1 or B.denominator != 1:
        # use integer model: multiply by denoms carefully for BJ monic Z
        # x = z / den → not monic BJ with Z coeffs always
        # Fall back: take numerators after common denominator
        den = sp.ilcm(A.denominator, B.denominator)
        A_i = int(A * den)
        B_i = int(B * den)
        # This is not the same Gal as x^5+A x+B unless den=1
        # Better: find m in Z with Z coeffs
        A_i, B_i = None, None
        for s2 in range(1, 80):
            m2 = Fraction(s2 * kv.denominator**2, 16)
            al = 256 * m2**2 - Fraction(3125) * (kv**4) / 256
            be = kv * al
            if al.denominator == 1 and be.denominator == 1 and al != 0:
                A_i, B_i = int(al), int(be)
                m = m2
                break
        if A_i is None:
            # force k=-4 integer path
            kv = Fraction(-4)
            m = Fraction(s_int, 1)
            A_i = int(m * m - 3125)
            B_i = int(-4 * A_i)
    else:
        A_i, B_i = int(A), int(B)

    digits = binary_word_to_ternary_digits(bin_data.collatz_itinerary, "01_to_012_run")
    ell = digits_to_int(digits) if digits else bin_data.n
    lat = TernaryLatticeElement(
        value=ell,
        ternary_digits=tuple(digits) if digits else (0,),
        construction="F3: itinerary→ternary for lattice label; k from popcount",
    )
    # template optional BJ-embed recovering (A_i,B_i) approximately
    tmpl = TemplateParams(
        a=0, b=-A_i if A_i else 0, c=-B_i if B_i else 1, d=0, e=0, f=1,
        construction="F3: rough BJ-embed sketch b~-α,c~-β (may not match chi exactly)",
    )
    # Fix BJ-embed properly: e=0,f=1,d=0,a=0, α=-(b f)= -b so b=-α, β=-b c = α c ⇒ c=β/α=k
    # c must be integer: k integer only. For general k use a=-ef with e,f free
    # Skip exact template match; BJ pair is the main image
    bj = BJPair(
        alpha=int(A_i),
        beta=int(B_i),
        k=str(kv),
        construction=f"F3: pure-even catalogue k={kv}, m={m}",
    )
    return TernaryOutput(
        lattice=lat, template=tmpl, bj=bj, source_n=bin_data.n, functor_id="F3"
    )


FUNCTORS = {
    "F1": F1_lattice,
    "F2": F2_height_map,
    "F3": F3_pure_even_k_from_itinerary,
}


# ---------------------------------------------------------------------------
# Image monodromy tests
# ---------------------------------------------------------------------------


def chi_from_template(tmpl: TemplateParams):
    aa, bb, cc, dd, ee, ff = tmpl.a, tmpl.b, tmpl.c, tmpl.d, tmpl.e, tmpl.f
    return sp.expand(
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


def analyse_poly_expr(expr) -> dict:
    try:
        pol = sp.Poly(sp.expand(expr), x, domain=sp.ZZ)
    except Exception:
        return {"status": "not_Z", "disc_square": False}
    if pol.degree() < 2:
        return {"status": "deg_low", "disc_square": False}
    if pol.LC() != 1:
        # make monic if LC=-1
        if pol.LC() == -1:
            pol = sp.Poly((-pol.as_expr()), x, domain=sp.ZZ)
        else:
            return {"status": "not_monic", "disc_square": False}
    disc = int(pol.discriminant())
    sq = disc > 0 and is_square(disc)
    if not pol.is_irreducible:
        return {
            "status": "reducible",
            "disc": disc,
            "disc_square": sq,
            "irreducible": False,
        }
    if not sq:
        rec = classify_poly(pol.as_expr(), do_galois=True)
        return {
            "status": rec.get("status") or "odd",
            "galois": rec.get("galois"),
            "disc": disc,
            "disc_square": False,
            "irreducible": True,
        }
    rec = classify_poly(pol.as_expr(), do_galois=True)
    return {
        "status": rec.get("status"),
        "galois": rec.get("galois"),
        "disc": disc,
        "disc_square": True,
        "irreducible": True,
    }


def run_functor_battery(functor_id: str, seeds: list[int]) -> dict:
    F = FUNCTORS[functor_id]
    rows = []
    stats = Counter()
    for n in seeds:
        bd = extract_binary(n)
        out = F(bd)
        # Analyse BJ image (primary for monodromy test)
        bj_stat = None
        if out.bj and out.bj.alpha != 0:
            expr = x**5 + out.bj.alpha * x + out.bj.beta
            bj_stat = analyse_poly_expr(expr)
            stats["bj_" + str(bj_stat.get("status"))] += 1
            if bj_stat.get("disc_square"):
                stats["bj_disc_square"] += 1
            if (bj_stat.get("status") or "").startswith("HIT_A5"):
                stats["bj_A5"] += 1
        # Template chi
        tmpl_stat = None
        if out.template:
            expr_t = chi_from_template(out.template)
            tmpl_stat = analyse_poly_expr(expr_t)
            stats["tmpl_" + str(tmpl_stat.get("status"))] += 1
            if tmpl_stat.get("disc_square"):
                stats["tmpl_disc_square"] += 1
            if (tmpl_stat.get("status") or "").startswith("HIT_A5"):
                stats["tmpl_A5"] += 1
        rows.append(
            {
                "n": n,
                "v2": bd.v2,
                "odd": bd.odd_part,
                "itin_len": bd.collatz_length,
                "lattice": out.lattice.value,
                "bj": asdict(out.bj) if out.bj else None,
                "bj_gal": bj_stat,
                "template": asdict(out.template) if out.template else None,
                "tmpl_gal": tmpl_stat,
            }
        )
    n_bj = sum(1 for r in rows if r["bj"])
    n_tmpl = sum(1 for r in rows if r["template"])
    return {
        "functor_id": functor_id,
        "n_seeds": len(seeds),
        "n_bj_images": n_bj,
        "n_tmpl_images": n_tmpl,
        "bj_disc_square_rate": stats.get("bj_disc_square", 0) / max(n_bj, 1),
        "bj_A5_rate": stats.get("bj_A5", 0) / max(n_bj, 1),
        "tmpl_disc_square_rate": stats.get("tmpl_disc_square", 0) / max(n_tmpl, 1),
        "tmpl_A5_rate": stats.get("tmpl_A5", 0) / max(n_tmpl, 1),
        "status_counts": dict(stats),
        "sample_rows": rows[:15],
        "all_rows_summary": [
            {
                "n": r["n"],
                "lattice": r["lattice"],
                "bj_status": (r["bj_gal"] or {}).get("status"),
                "bj_disc_sq": (r["bj_gal"] or {}).get("disc_square"),
                "tmpl_status": (r["tmpl_gal"] or {}).get("status"),
            }
            for r in rows
        ],
    }


def random_control(seeds_height: list[int], n_samples: int = 40) -> dict:
    """Random BJ with |α|,|β| comparable to F-image heights — baseline rates."""
    import random

    rng = random.Random(539)
    # height scale from seeds
    heights = []
    for n in seeds_height[:20]:
        heights.append(max(n, 10))
    H = max(heights) if heights else 100
    stats = Counter()
    n_ok = 0
    for _ in range(n_samples):
        aa = rng.randint(-H, H)
        bb = rng.randint(-H, H)
        if aa == 0 or bb == 0:
            continue
        n_ok += 1
        st = analyse_poly_expr(x**5 + aa * x + bb)
        if st.get("disc_square"):
            stats["disc_square"] += 1
        if (st.get("status") or "").startswith("HIT_A5"):
            stats["A5"] += 1
        stats[st.get("status") or "na"] += 1
    return {
        "n_samples": n_ok,
        "disc_square_rate": stats.get("disc_square", 0) / max(n_ok, 1),
        "A5_rate": stats.get("A5", 0) / max(n_ok, 1),
        "note": "Random monic BJ with |coeffs|≤H; not pure-even constrained",
    }


def main():
    t0 = time.time()
    print("CANDIDATE C FUNCTOR — Tier 1.2", flush=True)

    # Seed set: binary-rich integers
    seeds = []
    for n in range(1, 80):
        seeds.append(n)
    for n in [128, 256, 512, 1024, 27, 81, 243, 61, 80, 539, 4880, 100, 200, 300]:
        if n not in seeds:
            seeds.append(n)
    # Collatz interesting
    for n in [7, 27, 31, 41, 47, 54, 63, 73, 97]:
        if n not in seeds:
            seeds.append(n)

    results = {}
    for fid in FUNCTORS:
        print(f"  running {fid} on {len(seeds)} seeds...", flush=True)
        results[fid] = run_functor_battery(fid, seeds)
        r = results[fid]
        print(
            f"    BJ disc□ rate={r['bj_disc_square_rate']:.3f} A5={r['bj_A5_rate']:.3f} | "
            f"Tmpl disc□={r['tmpl_disc_square_rate']:.3f} A5={r['tmpl_A5_rate']:.3f}",
            flush=True,
        )

    ctrl = random_control(seeds, n_samples=50)
    print(
        f"  random BJ control disc□={ctrl['disc_square_rate']:.3f} A5={ctrl['A5_rate']:.3f}",
        flush=True,
    )

    # Type signatures for the document
    type_sig = {
        "BinaryData": {
            "n": "N>0",
            "v2": "N_0 = v_2(n)",
            "odd_part": "n/2^{v2}",
            "collatz_itinerary": "{0,1}* truncated Collatz word",
            "collatz_length": "N_0",
        },
        "TernaryOutput": {
            "lattice": "TernaryLatticeElement",
            "template": "TemplateParams | None  → T(a,b,c,d,e,f)",
            "bj": "BJPair | None  → x^5+αx+β",
            "functor_id": "F1|F2|F3",
        },
        "functors": {
            "F1": "itinerary→ternary base-3 + M-deform + LSW pure-even",
            "F2": "v2→3^{v2}*residue + BJ-embed template",
            "F3": "popcount→catalogue k + pure-even envelope (disc□ by construction)",
        },
    }

    # Interpretation
    f3 = results["F3"]
    f1 = results["F1"]
    f2 = results["F2"]
    interpretation = {
        "F1": (
            "Design-faithful (itinerary mirror). BJ image uses pure-even LSW so disc□ "
            f"rate={f1['bj_disc_square_rate']:.2f} (by construction on LSW ray). "
            f"Template M-deform disc□ rate={f1['tmpl_disc_square_rate']:.2f} "
            "(not forced even — same obstruction as base M)."
        ),
        "F2": (
            "Valuation transport. BJ-embed image has disc□ only when (α,β) hit even locus; "
            f"observed BJ disc□={f2['bj_disc_square_rate']:.2f}, A5={f2['bj_A5_rate']:.2f}. "
            "Template often non-BJ or odd."
        ),
        "F3": (
            "Strongest even monodromy by design: pure-even catalogue k. "
            f"BJ disc□ rate={f3['bj_disc_square_rate']:.2f}, A5={f3['bj_A5_rate']:.2f}. "
            "Evenness is pure-even theory, not a new force from the binary input."
        ),
        "preservation_vs_enrichment": (
            "Binary Collatz data has no native Gal monodromy to 'preserve'. "
            "Tests are enrichment rates on F(image) vs random BJ control. "
            "F3 enriches disc□ to ~1 by invoking pure-even; that is not necessity from "
            "binary hypotheses alone — pure-even is inserted in the codomain construction."
        ),
        "necessity": (
            "None of F1–F3 proves that HQCC axioms force A_n. They make Candidate C "
            "checkable: explicit maps exist; even monodromy on the image is either "
            "imported (F1/F3 pure-even) or partial (F2)."
        ),
    }

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"Candidate C functors F1–F3 ({elapsed}s). "
        f"F1 BJ disc□={f1['bj_disc_square_rate']:.2f} A5={f1['bj_A5_rate']:.2f}; "
        f"F2 BJ disc□={f2['bj_disc_square_rate']:.2f} A5={f2['bj_A5_rate']:.2f}; "
        f"F3 BJ disc□={f3['bj_disc_square_rate']:.2f} A5={f3['bj_A5_rate']:.2f}; "
        f"random control disc□={ctrl['disc_square_rate']:.2f}. "
        f"Functor explicit; even monodromy on image is imported pure-even or partial — "
        f"not a necessity theorem."
    )
    print(verdict, flush=True)

    lines = [
        r"# Candidate C functor — Tier 1.2",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        r"Turns the mod-2 ↔ mod-3 design mirror into **explicit maps** with typed "
        r"input/output. Does **not** prove necessity.",
        "",
        "---",
        "",
        r"## Type signature",
        "",
        r"### Input: `BinaryData`",
        "",
        r"| field | type | meaning |",
        r"|-------|------|---------|",
        r"| `n` | \(\mathbb{N}_{>0}\) | binary seed |",
        r"| `v2` | \(\mathbb{N}_0\) | \(v_2(n)\) |",
        r"| `odd_part` | odd \(\mathbb{N}\) | \(n/2^{v_2}\) |",
        r"| `collatz_itinerary` | \(\{0,1\}^*\) | Collatz word: 0=even→\(n/2\), 1=odd→\(3n+1\) |",
        r"| `collatz_length` | \(\mathbb{N}_0\) | truncation length |",
        "",
        r"### Output: `TernaryOutput`",
        "",
        r"| field | type | meaning |",
        r"|-------|------|---------|",
        r"| `lattice` | `TernaryLatticeElement` | integer from ternary digits / model mix |",
        r"| `template` | `TemplateParams` | \((a,b,c,d,e,f)\) for \(T\) |",
        r"| `bj` | `BJPair` | \((\alpha,\beta)\) for \(x^5+\alpha x+\beta\) |",
        r"| `functor_id` | string | F1 / F2 / F3 |",
        "",
        r"$$F:\ \mathbf{BinaryData}\ \longrightarrow\ \mathbf{TernaryOutput}.$$",
        "",
        "---",
        "",
        r"## Three concrete functors",
        "",
        r"### F1 — itinerary evaluation (design-faithful)",
        "",
        r"1. Map Collatz word bit \(0\mapsto 0\), \(1\mapsto 1\) or \(2\) (alternating).  ",
        r"2. Evaluate as base-3 integer; add model generator indexed by \(v_2\).  ",
        r"3. Template: deform base \(M\) by \(f=\ell \bmod 243\).  ",
        r"4. BJ: LSW pure-even ray \(k=-4\) with \(m=3^{\min(v_2,8)}+(\ell\bmod 16)\).",
        "",
        f"- BJ disc□ rate: **{f1['bj_disc_square_rate']:.3f}** (pure-even by construction)  ",
        f"- BJ A5 rate: **{f1['bj_A5_rate']:.3f}**  ",
        f"- Template disc□ rate: **{f1['tmpl_disc_square_rate']:.3f}**  ",
        f"- Template A5 rate: **{f1['tmpl_A5_rate']:.3f}**",
        "",
        r"### F2 — valuation transport",
        "",
        r"1. \(\ell = 3^{v_2}\cdot(\mathrm{odd}\bmod 3^5)+61\).  ",
        r"2. BJ-embed template: \(d=0\), \(a=-ef\), \(e=3\), \(f=v_2+1\), \(c=61\), \(b\in\{0,80\}\).  ",
        r"3. BJ pair from embed formulae \(\alpha=-(bf+ce)\), \(\beta=-bc\).",
        "",
        f"- BJ disc□ rate: **{f2['bj_disc_square_rate']:.3f}**  ",
        f"- BJ A5 rate: **{f2['bj_A5_rate']:.3f}**  ",
        f"- Template disc□ rate: **{f2['tmpl_disc_square_rate']:.3f}**  ",
        f"- Template A5 rate: **{f2['tmpl_A5_rate']:.3f}**",
        "",
        r"### F3 — catalogue \(k\) + pure-even envelope",
        "",
        r"1. Choose multi-seed \(k\) by popcount(itinerary) mod 8.  ",
        r"2. Pure-even \((\alpha,\beta)\) with integer coefficients when possible.  ",
        r"3. Lattice label from itinerary base-3.",
        "",
        f"- BJ disc□ rate: **{f3['bj_disc_square_rate']:.3f}** (by pure-even design)  ",
        f"- BJ A5 rate: **{f3['bj_A5_rate']:.3f}**  ",
        f"- Template disc□ rate: **{f3['tmpl_disc_square_rate']:.3f}**",
        "",
        "---",
        "",
        r"## Control and comparison",
        "",
        f"Random BJ control (\(|\alpha|,|\beta|\leq H\)): disc□ rate "
        f"**{ctrl['disc_square_rate']:.3f}**, A5 rate **{ctrl['A5_rate']:.3f}**.",
        "",
        r"| map | BJ disc□ | BJ A5 | vs control disc□ |",
        r"|-----|--------:|------:|:----------------:|",
        f"| F1 | {f1['bj_disc_square_rate']:.3f} | {f1['bj_A5_rate']:.3f} | enriched (pure-even) |",
        f"| F2 | {f2['bj_disc_square_rate']:.3f} | {f2['bj_A5_rate']:.3f} | partial |",
        f"| F3 | {f3['bj_disc_square_rate']:.3f} | {f3['bj_A5_rate']:.3f} | enriched (pure-even) |",
        f"| random BJ | {ctrl['disc_square_rate']:.3f} | {ctrl['A5_rate']:.3f} | baseline |",
        "",
        "---",
        "",
        r"## Sample images",
        "",
    ]
    for fid in ("F1", "F2", "F3"):
        lines.append(f"### {fid} samples")
        lines.append("")
        lines.append(r"| \(n\) | lattice | BJ status | disc□ | tmpl status |")
        lines.append(r"|----:|--------:|-----------|:-----:|-------------|")
        for r in results[fid]["all_rows_summary"][:12]:
            lines.append(
                f"| {r['n']} | {r['lattice']} | {r['bj_status']} | "
                f"{r['bj_disc_sq']} | {r['tmpl_status']} |"
            )
        lines.append("")

    lines += [
        "---",
        "",
        r"## Interpretation",
        "",
        interpretation["F1"],
        "",
        interpretation["F2"],
        "",
        interpretation["F3"],
        "",
        f"**Preservation vs enrichment:** {interpretation['preservation_vs_enrichment']}",
        "",
        f"**Necessity:** {interpretation['necessity']}",
        "",
        "---",
        "",
        r"## Relation to the four-face principle",
        "",
        r"F1–F3 make the **dynamics → lattice → matrices → (optional) Galois** pipeline "
        r"into total functions on binary data. The 3-cycle face is encouraged by "
        r"ternary digits and template couplings; the sign face is supplied by pure-even "
        r"when we choose F1/F3 BJ outputs — i.e. the organising principle is "
        r"**implemented**, not newly forced.",
        "",
        "---",
        "",
        r"## What to do next (functor track)",
        "",
        r"1. Refine F2 so BJ-embed \((\alpha,\beta)\) lie on a pure-even \(k\)-slice "
        r"   determined by binary data (merge F2 height with F3 evenness without "
        r"   hard-coding catalogue \(k\)).",
        r"2. Define a binary **hypothesis** class \(\mathcal{H}\) (e.g. bounded Collatz "
        r"   height) and prove a statement of the form "
        r"   \(n\in\mathcal{H}\Rightarrow \mathrm{disc}(F_3(n))\ \square\) — still uses "
        r"   pure-even in \(F\), but is a real lemma about the composite.",
        r"3. Seek an F that outputs **only** template parameters (no pure-even insert) "
        r"   yet has disc□ rate \(\to 1\) — that would be news for Criterion 2.",
        "",
        r"```bash",
        r"python candidate_c_functor.py",
        r"```",
        "",
        r"_Generated by candidate_c_functor.py — Tier 1.2_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "type_signature": type_sig,
        "results": results,
        "random_control": ctrl,
        "interpretation": interpretation,
        "seeds_n": len(seeds),
    }
    # JSON-serialize: remove non-serializable if any
    write_md(ROOT / "CANDIDATE_C_FUNCTOR.md", "\n".join(lines))
    write_json(ROOT / "CANDIDATE_C_FUNCTOR.json", payload)
    write_md(OUT / "CANDIDATE_C_FUNCTOR.md", "\n".join(lines))
    write_json(OUT / "CANDIDATE_C_FUNCTOR.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "CANDIDATE_C_FUNCTOR.md", "\n".join(lines))
    except Exception:
        pass

    print(f"Wrote CANDIDATE_C_FUNCTOR.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
