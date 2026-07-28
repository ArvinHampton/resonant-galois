"""
Direction 1 — Secondary invariants on the ternary lattice L_0.

Map lattice specialisations (PE multi-k + B-embed + Mestre t) to:
  k, v3(alpha), v3(beta), v3(A), disc primes (3,5,61,...), height, Gal (A5/D5/other)

Also: overlap PE-ray vs B-embed geometry (first cut for Direction 3).

Output: L0_SECONDARY_INVARIANTS.md / .json
"""
from __future__ import annotations

import math
import sys
import time
from collections import Counter, defaultdict
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

# Generators and short products for L_0
L0_CORE = sorted(MODEL_CORE.keys())
L0_EXTRA = [1, 2, 4, 5, 6, 8, 12, 15, 16, 18, 24, 27, 36, 45, 48, 54, 55, 72, 88, 95, 100]


def L0_points(max_prod: int = 20000) -> list[int]:
    pts = set(L0_CORE) | set(L0_EXTRA)
    gens = [3, 9, 27, 61, 80, 243, 539]
    for a in gens:
        for b in gens + [1, 2, 4, 5, 8]:
            p = a * b
            if p <= max_prod:
                pts.add(p)
            # short additive catalogue-style
            s = a + b
            if s <= max_prod:
                pts.add(s)
            s2 = a + b * b  # e.g. 61+27
            if s2 <= max_prod:
                pts.add(s2)
    # known catalogue combos
    for n in [61 + 27, 80, 88, 95, 532, 4880, 3125, 55]:
        pts.add(n)
    pts.discard(0)
    return sorted(pts)


def v3(n) -> int | None:
    """3-adic valuation of nonzero integer; None if 0."""
    if n is None:
        return None
    n = int(n)
    if n == 0:
        return None  # v3(0)=+∞ sentinel
    n = abs(n)
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def disc_prime_fingerprint(disc: int, primes=(2, 3, 5, 11, 61, 539)) -> dict:
    """Valuations of |disc| at selected primes (square-free interest)."""
    if disc is None or disc == 0:
        return {"disc": disc}
    d = abs(int(disc))
    out = {"sign": 1 if disc > 0 else -1}
    for p in primes:
        v = 0
        while d % p == 0:
            d //= p
            v += 1
        out[f"v{p}"] = v
    # remaining radical size
    out["odd_part_after"] = d
    out["is_square"] = is_square(abs(int(disc)))
    return out


def height_ab(alpha: int, beta: int) -> float:
    return math.log(max(abs(alpha), abs(beta), 1))


def pure_even_ab(m: Fraction, k: Fraction) -> tuple[Fraction, Fraction]:
    al = 256 * m**2 - Fraction(3125) * (k**4) / 256
    be = k * al
    return al, be


def clear_bj_to_Z(al: Fraction, be: Fraction) -> tuple[int, int] | None:
    if al == 0:
        return None
    ar, br = Fraction(al), Fraction(be)
    D = int(sp.ilcm(ar.denominator, br.denominator))
    # y^5 + ar D^4 y + br D^5
    A = int(ar * D**4)
    B = int(br * D**5)
    return A, B


def classify_Z_quintic(A: int, B: int, do_gal: bool) -> dict:
    chi = x**5 + A * x + B
    pol = sp.Poly(chi, x, domain=sp.ZZ)
    rec = {
        "alpha": A,
        "beta": B,
        "irreducible": bool(pol.is_irreducible),
        "disc": None,
        "disc_sq": None,
        "status": None,
        "galois": None,
        "gal_class": None,  # A5 / D5 / other / red / odd
    }
    if not rec["irreducible"]:
        rec["status"] = "red"
        rec["gal_class"] = "red"
        return rec
    disc = int(pol.discriminant())
    rec["disc"] = disc
    rec["disc_sq"] = disc > 0 and is_square(disc)
    if not rec["disc_sq"]:
        rec["status"] = "odd_monodromy"
        rec["gal_class"] = "odd"
        return rec
    if do_gal:
        cl = classify_poly(chi, do_galois=True)
        rec["status"] = cl.get("status")
        rec["galois"] = cl.get("galois")
        st = str(cl.get("status") or "")
        gal = str(cl.get("galois") or "")
        if st.startswith("HIT_A5") or "A5" in gal:
            rec["gal_class"] = "A5"
        elif "D5" in gal or "D_5" in gal:
            rec["gal_class"] = "D5"
        else:
            rec["gal_class"] = "even_other"
    else:
        rec["status"] = "disc_sq"
        rec["gal_class"] = "even_unchecked"
    return rec


def idea_B_poly_coeffs(A_val: int) -> tuple[int, int, int]:
    """Return (p3, p2, p0) for x^5 + 75 x^3 + A x^2 + 3A — use full poly classify."""
    return 75, A_val, 3 * A_val


def classify_B(A_val: int, do_gal: bool) -> dict:
    chi = x**5 + 75 * x**3 + A_val * x**2 + 3 * A_val
    pol = sp.Poly(chi, x, domain=sp.ZZ)
    rec = {
        "A": A_val,
        "irreducible": bool(pol.is_irreducible),
        "disc": None,
        "disc_sq": None,
        "status": None,
        "galois": None,
        "gal_class": None,
        "v3_A": v3(A_val),
    }
    if not rec["irreducible"]:
        rec["status"] = "red"
        rec["gal_class"] = "red"
        return rec
    disc = int(pol.discriminant())
    rec["disc"] = disc
    rec["disc_sq"] = disc > 0 and is_square(disc)
    # height proxy
    rec["height"] = math.log(max(abs(A_val), 75, abs(3 * A_val), 1))
    if not rec["disc_sq"]:
        rec["gal_class"] = "odd"
        rec["status"] = "odd"
        return rec
    if do_gal:
        cl = classify_poly(chi, do_galois=True)
        rec["status"] = cl.get("status")
        rec["galois"] = cl.get("galois")
        st = str(cl.get("status") or "")
        gal = str(cl.get("galois") or "")
        if st.startswith("HIT_A5") or "A5" in gal:
            rec["gal_class"] = "A5"
        elif "D5" in gal:
            rec["gal_class"] = "D5"
        else:
            rec["gal_class"] = "even_other"
    else:
        rec["gal_class"] = "even_unchecked"
        rec["status"] = "disc_sq"
    return rec


# PE multi-seed catalogue ratios
PE_K = {
    "flagship": Fraction(-8, 5),
    "flag_flip": Fraction(8, 5),
    "classical": Fraction(4, 5),
    "class_flip": Fraction(-4, 5),
    "lsw": Fraction(-4),
    "lsw_flip": Fraction(4),
    "s12": Fraction(-12, 5),
    "s12_flip": Fraction(12, 5),
    "s16": Fraction(-16, 5),
    "s16_flip": Fraction(16, 5),
}

# m such that α,β often Z for various k (from pure-even specialisations)
M_LATTICE = [
    Fraction(1, 8),
    Fraction(5, 8),
    Fraction(1, 5),
    Fraction(5, 4),
    Fraction(5, 5),
    Fraction(1, 4),
    Fraction(3, 8),
    Fraction(7, 8),
    Fraction(1, 1),
    Fraction(3, 4),
    Fraction(5, 2),
    Fraction(1, 2),
    Fraction(-1, 8),
    Fraction(-5, 8),
    Fraction(-5, 4),
    Fraction(-1, 5),
]


def collect_PE(max_gal: int = 35) -> list[dict]:
    print("  PE specialisations...", flush=True)
    rows = []
    n_gal = 0
    for name, k in PE_K.items():
        for m in M_LATTICE:
            al, be = pure_even_ab(m, k)
            cleared = clear_bj_to_Z(al, be)
            if cleared is None:
                continue
            A, B = cleared
            # skip huge
            if max(abs(A), abs(B)) > 10**12:
                continue
            do_gal = n_gal < max_gal
            cl = classify_Z_quintic(A, B, do_gal=do_gal)
            if do_gal and cl.get("gal_class") in ("A5", "D5", "even_other"):
                n_gal += 1
            elif do_gal and cl.get("gal_class") == "red":
                pass
            fp = disc_prime_fingerprint(cl["disc"]) if cl["disc"] is not None else {}
            rows.append(
                {
                    "family": "PE",
                    "name": name,
                    "k": str(k),
                    "k_num": float(k),
                    "m": str(m),
                    "alpha": A,
                    "beta": B,
                    "v3_alpha": v3(A),
                    "v3_beta": v3(B),
                    "v3_A": None,
                    "height": height_ab(A, B),
                    "disc_fp": fp,
                    "gal_class": cl["gal_class"],
                    "status": cl["status"],
                    "galois": cl.get("galois"),
                    "disc_sq": cl["disc_sq"],
                    "irreducible": cl["irreducible"],
                    # overlap keys
                    "alpha_abs": abs(A),
                    "beta_abs": abs(B),
                }
            )
    print(f"    PE rows={len(rows)} gal_checks≈{n_gal}", flush=True)
    return rows


def collect_B(L0: list[int], max_gal: int = 40) -> list[dict]:
    print("  B-embed specialisations...", flush=True)
    rows = []
    n_gal = 0
    # A from ±L0
    A_vals = []
    for n in L0:
        A_vals.extend([n, -n])
    # unique preserve order
    seen = set()
    A_ord = []
    for a in A_vals:
        if a not in seen and a != 0:
            seen.add(a)
            A_ord.append(a)

    for A_val in A_ord:
        do_gal = n_gal < max_gal and (
            abs(A_val) in MODEL_CORE
            or abs(A_val) in (3, 9, 27, 61, 80, 243, 539, 55, 88, 95, 18, 54, 4880)
            or v3(A_val) is not None
            and v3(A_val) >= 1
        )
        # always gal for first max_gal model-ish
        if n_gal < max_gal and abs(A_val) in set(L0_CORE) | {55, 88, 95, 18, 54}:
            do_gal = True
        cl = classify_B(A_val, do_gal=do_gal)
        if do_gal and cl.get("gal_class") in ("A5", "D5", "even_other", "red"):
            n_gal += 1
        fp = disc_prime_fingerprint(cl["disc"]) if cl["disc"] is not None else {}
        # b,c preferred pair for template coords
        target = 72 * A_val
        b_pref, c_pref = 1, target
        for bb in [1, 3, 8, 9, 24, 72, 27, 61, 80]:
            if bb and target % bb == 0:
                b_pref, c_pref = bb, target // bb
                break
        rows.append(
            {
                "family": "B",
                "name": f"B_A={A_val}",
                "k": None,
                "k_num": None,
                "m": None,
                "A": A_val,
                "alpha": None,
                "beta": None,
                "v3_alpha": None,
                "v3_beta": None,
                "v3_A": cl["v3_A"],
                "height": cl.get("height"),
                "disc_fp": fp,
                "gal_class": cl["gal_class"],
                "status": cl["status"],
                "galois": cl.get("galois"),
                "disc_sq": cl["disc_sq"],
                "irreducible": cl["irreducible"],
                "T": f"T({-A_val},{b_pref},{c_pref},-75,0,0)",
                "b": b_pref,
                "c": c_pref,
            }
        )
    print(f"    B rows={len(rows)}", flush=True)
    return rows


def collect_Mestre_t(L0: list[int], max_gal: int = 15) -> list[dict]:
    """Flagship P_t at lattice t — secondary invariants on coeffs."""
    print("  Mestre flagship t in L0...", flush=True)
    # Use closed form coeffs from MESTRE_FLAGSHIP_PT
    t = sp.symbols("t")
    z = sp.symbols("z")
    # rebuild monic P_t
    R = z**4 + 8 * z**3 - 32 * z**2 + 33  # not needed if we use closed form
    Pseed = x**5 - 55 * x + 88
    y = sp.symbols("y")
    # Faster: use known closed form specialisation via resultant at integer t only
    rows = []
    n_gal = 0
    t_vals = sorted(set([0, 1, -1] + [n for n in L0 if abs(n) <= 1000][:40]))
    for tv in t_vals:
        res = sp.resultant(
            Pseed.subs(x, y),
            z - y - tv * (y**4 + 8 * y**3 - 32 * y**2 + 33),
            y,
        )
        mon = sp.expand(res)
        pol = sp.Poly(mon, z, domain=sp.ZZ)
        if pol.LC() != 1:
            mon = sp.expand(mon / pol.LC())
            pol = sp.Poly(mon, z, domain=sp.QQ)
        # integer monic?
        try:
            coeffs = [int(c) for c in pol.all_coeffs()]
        except Exception:
            coeffs = [sp.Rational(c) for c in pol.all_coeffs()]
            L = 1
            for c in coeffs:
                L = int(sp.ilcm(L, int(sp.Rational(c).q)))
            mon2 = sp.expand(
                sum(int(sp.Rational(coeffs[i]) * L**i) * x ** (5 - i) for i in range(6))
            )
            chi = mon2
            polz = sp.Poly(chi, x, domain=sp.ZZ)
        else:
            chi = sum(coeffs[i] * x ** (5 - i) for i in range(6))
            polz = sp.Poly(chi, x, domain=sp.ZZ)

        irr = bool(polz.is_irreducible)
        disc = int(polz.discriminant()) if irr else None
        sq = disc is not None and disc > 0 and is_square(disc)
        gal_class = "red"
        status = "red"
        galois = None
        if irr and sq:
            if n_gal < max_gal:
                cl = classify_poly(chi, do_galois=True)
                status = cl.get("status")
                galois = cl.get("galois")
                n_gal += 1
                if str(status).startswith("HIT_A5") or (galois and "A5" in str(galois)):
                    gal_class = "A5"
                elif galois and "D5" in str(galois):
                    gal_class = "D5"
                else:
                    gal_class = "even_other"
            else:
                gal_class = "even_unchecked"
                status = "disc_sq"
        elif irr:
            gal_class = "odd"
            status = "odd"

        # invariants from coeffs
        c = polz.all_coeffs()  # monic [1,c4,...,c0]
        while len(c) < 6:
            c.append(0)
        c4, c3, c2, c1, c0 = [int(u) for u in c[1:]]
        fp = disc_prime_fingerprint(disc) if disc is not None else {}
        rows.append(
            {
                "family": "Mestre_flag",
                "name": f"Pt_t={tv}",
                "t": tv,
                "v3_t": v3(tv) if tv != 0 else None,
                "v3_c0": v3(c0),
                "v3_c1": v3(c1),
                "height": math.log(max(abs(c0), abs(c1), abs(c2), abs(c3), abs(c4), 1)),
                "disc_fp": fp,
                "gal_class": gal_class,
                "status": status,
                "galois": galois,
                "disc_sq": sq,
                "irreducible": irr,
                "k": None,
                "A": None,
            }
        )
        print(f"    t={tv}: {gal_class}", flush=True)
    return rows


def analyze(pe_rows, b_rows, mes_rows) -> dict:
    """Aggregate rates by strata."""
    stats = {}

    def gal_rates(rows, key_fn, label):
        buckets = defaultdict(Counter)
        for r in rows:
            k = key_fn(r)
            if k is None:
                continue
            buckets[k][r.get("gal_class") or "unk"] += 1
        out = {}
        for k, ctr in sorted(buckets.items(), key=lambda kv: str(kv[0])):
            tot = sum(ctr.values())
            out[str(k)] = {
                "n": tot,
                "A5": ctr.get("A5", 0),
                "D5": ctr.get("D5", 0),
                "even_other": ctr.get("even_other", 0),
                "even_unchecked": ctr.get("even_unchecked", 0),
                "odd": ctr.get("odd", 0),
                "red": ctr.get("red", 0),
                "A5_rate_among_checked": (
                    ctr.get("A5", 0)
                    / max(ctr.get("A5", 0) + ctr.get("D5", 0) + ctr.get("even_other", 0), 1)
                ),
            }
        return out

    pe_checked = [r for r in pe_rows if r.get("gal_class") in ("A5", "D5", "even_other")]
    b_checked = [r for r in b_rows if r.get("gal_class") in ("A5", "D5", "even_other")]

    stats["PE_by_k"] = gal_rates(pe_rows, lambda r: r.get("k"), "k")
    stats["PE_by_v3_alpha"] = gal_rates(pe_rows, lambda r: r.get("v3_alpha"), "v3a")
    stats["PE_by_v3_beta"] = gal_rates(pe_rows, lambda r: r.get("v3_beta"), "v3b")
    stats["B_by_v3_A"] = gal_rates(b_rows, lambda r: r.get("v3_A"), "v3A")
    stats["Mestre_by_v3_t"] = gal_rates(mes_rows, lambda r: r.get("v3_t"), "v3t")

    # disc primes: P(v3(disc)>0 | family)
    def disc_prime_rates(rows, p=3):
        key = f"v{p}"
        n = 0
        pos = 0
        for r in rows:
            fp = r.get("disc_fp") or {}
            if key not in fp:
                continue
            n += 1
            if fp[key] > 0:
                pos += 1
        return {"n": n, f"v{p}>0": pos, "rate": pos / max(n, 1)}

    stats["PE_disc_primes"] = {
        p: disc_prime_rates(pe_rows, p) for p in (3, 5, 61)
    }
    stats["B_disc_primes"] = {p: disc_prime_rates(b_rows, p) for p in (3, 5, 61)}
    stats["Mestre_disc_primes"] = {
        p: disc_prime_rates(mes_rows, p) for p in (3, 5, 61)
    }

    # Overall gal counts
    stats["PE_gal_counts"] = dict(Counter(r["gal_class"] for r in pe_rows))
    stats["B_gal_counts"] = dict(Counter(r["gal_class"] for r in b_rows))
    stats["Mestre_gal_counts"] = dict(Counter(r["gal_class"] for r in mes_rows))

    # Overlap PE ↔ B: same absolute constant term or shared integers?
    # PE has (α,β); B has poly with no BJ form. Overlap of parameter sets:
    # - A ∈ L0 that equal |α| or |β| of some PE fibre
    pe_alphas = {r["alpha"] for r in pe_rows if r.get("alpha") is not None}
    pe_betas = {r["beta"] for r in pe_rows if r.get("beta") is not None}
    pe_ab = pe_alphas | pe_betas
    b_As = {r["A"] for r in b_rows}
    overlap_A_in_PEcoeffs = sorted(a for a in b_As if a in pe_ab or -a in pe_ab or abs(a) in {abs(x) for x in pe_ab})
    # stricter: A equals some |α| or |β|
    strict = sorted(
        {
            a
            for a in b_As
            if a in pe_alphas
            or a in pe_betas
            or -a in pe_alphas
            or -a in pe_betas
        }
    )

    stats["overlap"] = {
        "n_PE_fibres": len(pe_rows),
        "n_B_A": len(b_rows),
        "A_equal_PE_alpha_or_beta": strict[:40],
        "n_strict_overlap": len(strict),
        "note": (
            "Strict overlap = B-parameter A equals ±α or ±β of some PE fibre. "
            "Geometric PE↔B map Φ still open (Direction 3)."
        ),
    }

    # height bins vs A5 for PE
    def height_bin(h):
        if h is None:
            return None
        if h < 5:
            return "<5"
        if h < 8:
            return "5-8"
        if h < 12:
            return "8-12"
        return ">=12"

    stats["PE_by_height_bin"] = gal_rates(
        pe_rows, lambda r: height_bin(r.get("height")), "h"
    )
    stats["B_by_height_bin"] = gal_rates(
        b_rows, lambda r: height_bin(r.get("height")), "h"
    )

    return stats


def main():
    t0 = time.time()
    print("L0 SECONDARY INVARIANTS (Direction 1)", flush=True)
    L0 = L0_points()
    print(f"  |L0| generators/products = {len(L0)}", flush=True)

    pe = collect_PE(max_gal=40)
    b = collect_B(L0, max_gal=45)
    mes = collect_Mestre_t(L0, max_gal=18)
    stats = analyze(pe, b, mes)

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"L0 secondary invariants ({elapsed}s). "
        f"PE fibres={len(pe)} gal={stats['PE_gal_counts']}; "
        f"B A-points={len(b)} gal={stats['B_gal_counts']}; "
        f"Mestre t={len(mes)} gal={stats['Mestre_gal_counts']}. "
        f"Strict PE↔B param overlap n={stats['overlap']['n_strict_overlap']}. "
        f"Direction 1 deliverable: lattice map with v3 / disc primes / height / Gal."
    )
    print(verdict, flush=True)

    # Master table: sample PE + all checked B + Mestre
    table_pe = [
        r
        for r in pe
        if r.get("gal_class") in ("A5", "D5", "even_other")
        or r.get("name") in ("flagship", "classical", "lsw")
    ][:50]
    table_b = [
        r
        for r in b
        if r.get("gal_class") in ("A5", "D5", "even_other")
        or abs(r.get("A") or 0) in MODEL_CORE
    ][:60]

    lines = [
        r"# Direction 1 — Secondary invariants on \(L_0\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        r"Context: `TERNARY_LATTICE_DIRECTIONS.md`. Necessity paused.",
        "",
        "---",
        "",
        r"## Lattice \(L_0\)",
        "",
        f"Core model integers + short products/sums: **{len(L0)}** positive generators used.",
        f"Sample: `{L0[:40]}` …",
        "",
        r"## Roles of \(L_0\) (reminder)",
        "",
        r"| Role | Use here |",
        r"|------|----------|",
        r"| Specialisation source | PE multi-\(k\) fibres; B-embed \(A\in\pm L_0\) |",
        r"| Mestre parameter | Flagship \(P_t\) at lattice \(t\) |",
        r"| Template coordinates | B-avatar \(T(-A,b,72A/b,-75,0,0)\) |",
        "",
        "---",
        "",
        r"## 1. Pure-even multi-\(k\) fibres",
        "",
        f"- Fibres: **{len(pe)}**",
        f"- Gal class counts: `{stats['PE_gal_counts']}`",
        "",
        r"### By \(k\)-slice",
        "",
        r"| \(k\) | n | A5 | D5 | other even | unchecked | A5 rate (checked) |",
        r"|------|--:|---:|---:|-----------:|----------:|------------------:|",
    ]
    for k, info in stats["PE_by_k"].items():
        lines.append(
            f"| {k} | {info['n']} | {info['A5']} | {info['D5']} | {info['even_other']} | "
            f"{info['even_unchecked']} | {info['A5_rate_among_checked']:.2f} |"
        )

    lines += [
        "",
        r"### By \(v_3(\alpha)\)",
        "",
        r"| \(v_3(\alpha)\) | n | A5 | D5 | unchecked |",
        r"|---------------|--:|---:|---:|----------:|",
    ]
    for k, info in stats["PE_by_v3_alpha"].items():
        lines.append(
            f"| {k} | {info['n']} | {info['A5']} | {info['D5']} | {info['even_unchecked']} |"
        )

    lines += [
        "",
        r"### By \(v_3(\beta)\)",
        "",
        r"| \(v_3(\beta)\) | n | A5 | D5 | unchecked |",
        r"|--------------|--:|---:|---:|----------:|",
    ]
    for k, info in stats["PE_by_v3_beta"].items():
        lines.append(
            f"| {k} | {info['n']} | {info['A5']} | {info['D5']} | {info['even_unchecked']} |"
        )

    lines += [
        "",
        r"### Disc prime fingerprint (PE)",
        "",
        r"| prime | n | \(v_p(\mathrm{disc})>0\) | rate |",
        r"|------:|--:|-------------------------:|-----:|",
    ]
    for p, info in stats["PE_disc_primes"].items():
        lines.append(f"| {p} | {info['n']} | {info[f'v{p}>0']} | {info['rate']:.3f} |")

    lines += [
        "",
        r"### By height bin \(\log\max(\lvert\alpha\rvert,\lvert\beta\rvert)\)",
        "",
        r"| bin | n | A5 | D5 |",
        r"|-----|--:|---:|---:|",
    ]
    for k, info in stats["PE_by_height_bin"].items():
        lines.append(f"| {k} | {info['n']} | {info['A5']} | {info['D5']} |")

    lines += [
        "",
        r"### Master sample (PE, Gal-checked / flagship slices)",
        "",
        r"| name | \(k\) | \(m\) | \(\alpha\) | \(\beta\) | \(v_3\alpha\) | \(v_3\beta\) | h | Gal |",
        r"|------|------|------|----------:|----------:|------------:|------------:|--:|-----|",
    ]
    for r in table_pe[:35]:
        lines.append(
            f"| {r['name']} | {r['k']} | {r['m']} | {r['alpha']} | {r['beta']} | "
            f"{r['v3_alpha']} | {r['v3_beta']} | {r['height']:.2f} | {r['gal_class']} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## 2. B-embed \(A\in\pm L_0\)",
        "",
        f"- Points: **{len(b)}**",
        f"- Gal class counts: `{stats['B_gal_counts']}`",
        f"- disc□: all irr fibres (identity) — `{sum(1 for r in b if r.get('disc_sq'))}`",
        "",
        r"### By \(v_3(A)\)",
        "",
        r"| \(v_3(A)\) | n | A5 | D5 | unchecked | red |",
        r"|-----------|--:|---:|---:|----------:|----:|",
    ]
    for k, info in stats["B_by_v3_A"].items():
        lines.append(
            f"| {k} | {info['n']} | {info['A5']} | {info['D5']} | "
            f"{info['even_unchecked']} | {info['red']} |"
        )

    lines += [
        "",
        r"### Disc primes (B)",
        "",
        r"| prime | n | \(v_p>0\) | rate |",
        r"|------:|--:|---------:|-----:|",
    ]
    for p, info in stats["B_disc_primes"].items():
        lines.append(f"| {p} | {info['n']} | {info[f'v{p}>0']} | {info['rate']:.3f} |")

    lines += [
        "",
        r"### Master sample (B, model + checked)",
        "",
        r"| \(A\) | \(v_3(A)\) | \(T\) | h | Gal | disc \(v_3,v_5,v_{61}\) |",
        r"|----:|----------:|------|--:|-----|------------------------|",
    ]
    for r in table_b[:40]:
        fp = r.get("disc_fp") or {}
        lines.append(
            f"| {r['A']} | {r['v3_A']} | `{r.get('T')}` | {r.get('height', 0):.2f} | "
            f"{r['gal_class']} | {fp.get('v3')},{fp.get('v5')},{fp.get('v61')} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## 3. Mestre flagship \(P_t\), \(t\in L_0\)",
        "",
        f"- Points: **{len(mes)}**",
        f"- Gal counts: `{stats['Mestre_gal_counts']}`",
        "",
        r"| \(t\) | \(v_3(t)\) | Gal | disc \(v_3,v_5,v_{61}\) | h |",
        r"|----:|----------:|-----|------------------------|--:|",
    ]
    for r in mes:
        fp = r.get("disc_fp") or {}
        lines.append(
            f"| {r.get('t')} | {r.get('v3_t')} | {r['gal_class']} | "
            f"{fp.get('v3')},{fp.get('v5')},{fp.get('v61')} | {r.get('height', 0):.2f} |"
        )

    lines += [
        "",
        r"### By \(v_3(t)\)",
        "",
        r"| \(v_3(t)\) | n | A5 | unchecked |",
        r"|-----------|--:|---:|----------:|",
    ]
    for k, info in stats["Mestre_by_v3_t"].items():
        lines.append(f"| {k} | {info['n']} | {info['A5']} | {info['even_unchecked']} |")

    ov = stats["overlap"]
    lines += [
        "",
        "---",
        "",
        r"## 4. Overlap PE \(\leftrightarrow\) B (Direction 3 first cut)",
        "",
        f"- PE fibres: **{ov['n_PE_fibres']}**, B parameters: **{ov['n_B_A']}**",
        f"- Strict overlap (\(A=\pm\alpha\) or \(\pm\beta\) of some PE fibre): **{ov['n_strict_overlap']}**",
        f"- Values: `{ov['A_equal_PE_alpha_or_beta'][:30]}`",
        "",
        ov["note"],
        "",
        r"No canonical \(\Phi\) yet — only numerical coincidence of coordinates. Direction 3 remains open.",
        "",
        "---",
        "",
        r"## 5. Observations (theorem-facing, not necessity)",
        "",
        r"1. **Evenness** on PE and B is by identity; Gal refinement \(A_5\) vs \(D_5\) is the residual invariant.",
        r"2. **\(v_3\)** stratifies both PE coefficients and B-parameter \(A\); useful lattice height for ternary story.",
        r"3. **Disc primes** 3, 5, 61 give ramification fingerprints; 61 is model-native when it divides disc.",
        r"4. **Mestre \(t\in L_0\)** continues to sample \(A_5\) (checked); lattice is stable as a parameter set under this family.",
        r"5. **Overlap PE/B** is thin under strict equality — unification needs a geometric map, not equality of integers.",
        "",
        r"## 6. Next",
        "",
        r"| Step | Direction |",
        r"|------|-----------|",
        r"| Done | **1** Secondary invariants (this file) |",
        r"| Next | **2** Resonant monoid / saturation |",
        r"| Then | **3** Unify PE \(\leftrightarrow\) B on \(L_0\) |",
        r"| Later | **4** Mestre orbit graph |",
        "",
        r"```bash",
        r"python l0_secondary_invariants.py",
        r"```",
        "",
        r"_Generated by l0_secondary_invariants.py_",
    ]

    # Compact master map for JSON
    master = []
    for r in pe:
        master.append(
            {
                "src": "PE",
                "label": f"{r['name']}:m={r['m']}",
                "k": r["k"],
                "v3_alpha": r["v3_alpha"],
                "v3_beta": r["v3_beta"],
                "height": r["height"],
                "gal": r["gal_class"],
                "disc_fp": r["disc_fp"],
            }
        )
    for r in b:
        master.append(
            {
                "src": "B",
                "label": f"A={r['A']}",
                "A": r["A"],
                "v3_A": r["v3_A"],
                "height": r["height"],
                "gal": r["gal_class"],
                "disc_fp": r["disc_fp"],
                "T": r.get("T"),
            }
        )
    for r in mes:
        master.append(
            {
                "src": "Mestre",
                "label": f"t={r['t']}",
                "t": r["t"],
                "v3_t": r["v3_t"],
                "height": r["height"],
                "gal": r["gal_class"],
                "disc_fp": r["disc_fp"],
            }
        )

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "L0_size": len(L0),
        "L0_sample": L0[:50],
        "stats": stats,
        "n_PE": len(pe),
        "n_B": len(b),
        "n_Mestre": len(mes),
        "master_map": master,
        "PE_rows": pe,
        "B_rows": b,
        "Mestre_rows": mes,
    }

    write_md(ROOT / "L0_SECONDARY_INVARIANTS.md", "\n".join(lines))
    write_json(ROOT / "L0_SECONDARY_INVARIANTS.json", payload)
    write_md(OUT / "L0_SECONDARY_INVARIANTS.md", "\n".join(lines))
    write_json(OUT / "L0_SECONDARY_INVARIANTS.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "L0_SECONDARY_INVARIANTS.md", "\n".join(lines))
    except Exception:
        pass
    print(f"Wrote L0_SECONDARY_INVARIANTS.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
