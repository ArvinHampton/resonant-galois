"""
Next steps after NEW_ALGEBRAIC_IDEAS:

  1. Closed-form flagship Mestre family P_t
  2. Systematic HQCC-lattice points on B-embed bc = 72 A
  3. Matrix avatar with built-in evenness identity (not another scan)

Output: MESTRE_FLAGSHIP_PT.md, B_EMBED_LATTICE.md, EVENNESS_AVATAR.md
        + combined NEXT_MESTRE_B_AVATAR.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
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

t, z, y = sp.symbols("t z y")
a, b, c, d, e, f = sp.symbols("a b c d e f")
m, k = sp.symbols("m k")


# ---------------------------------------------------------------------------
# 1. Closed-form flagship P_t
# ---------------------------------------------------------------------------

FLAGSHIP_P = x**5 - 55 * x + 88
FLAGSHIP_R = x**4 + 8 * x**3 - 32 * x**2 + 33


def compute_flagship_Pt():
    """Exact monic P_t via Res_y(P(y), z - y - t R(y))."""
    print("  computing flagship P_t resultant...", flush=True)
    res = sp.resultant(
        FLAGSHIP_P.subs(x, y),
        z - y - t * FLAGSHIP_R.subs(x, y),
        y,
    )
    F = sp.expand(res)
    pol = sp.Poly(F, z)
    assert pol.LC() == 1
    mon = F  # already monic
    # Coefficients of z^4..z^0 as factored polys in t
    coeffs = {}
    powers = ["z^5", "z^4", "z^3", "z^2", "z", "const"]
    allc = pol.all_coeffs()  # [1, c4, c3, c2, c1, c0]
    labels = ["1", "c4", "c3", "c2", "c1", "c0"]
    factored = {}
    for lab, coef in zip(labels, allc):
        factored[lab] = str(sp.factor(sp.expand(coef)))
        coeffs[lab] = sp.expand(coef)

    # Compact closed form (human-readable)
    c4 = coeffs["c4"]
    c3 = coeffs["c3"]
    c2 = coeffs["c2"]
    c1 = coeffs["c1"]
    c0 = coeffs["c0"]

    closed_form = (
        r"P_t(z)=z^5"
        + rf" + ({sp.factor(c4)}) z^4"
        + rf" + ({sp.factor(c3)}) z^3"
        + rf" + ({sp.factor(c2)}) z^2"
        + rf" + ({sp.factor(c1)}) z"
        + rf" + ({sp.factor(c0)})"
    )

    # Disc identity in Q(t)
    print("  discriminant of P_t...", flush=True)
    D = sp.expand(pol.discriminant())
    # Test square in Q[t]
    num, den = sp.fraction(sp.together(D))
    num, den = sp.expand(num), sp.expand(den)

    def square_poly_Q(expr, var=t):
        if expr == 0:
            return True, "0"
        P = sp.Poly(expr, var, domain=sp.QQ)
        cont = sp.Rational(P.content())
        if cont < 0:
            return False, "neg_content"
        n, dd = int(sp.numer(cont)), int(sp.denom(cont))
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(dd, 2)[1]):
            return False, f"content:{cont}"
        prim = P.primitive()[1]
        if prim.degree() == 0:
            return True, "square"
        fac = sp.factor_list(prim.as_expr(), domain=sp.QQ)
        odds = [(str(fi)[:40], int(m)) for fi, m in fac[1] if m % 2]
        return (len(odds) == 0), (f"odd={odds}" if odds else "square")

    ok_n, info_n = square_poly_Q(num)
    ok_d, info_d = square_poly_Q(den)
    disc_square = ok_n and ok_d

    # Explicit sqrt candidate: try factor D completely
    Dfac = sp.factor(D)

    # Specialisations
    samples = []
    for tv in [0, 1, -1, 2, 3, 5, 7, 9, 61, 80]:
        Pt = sp.expand(mon.subs(t, tv)).subs(z, x)
        # monic Q already; integer coeffs for these t?
        polz = sp.Poly(Pt, x, domain=sp.QQ)
        if all(sp.Rational(c).q == 1 for c in polz.all_coeffs()):
            chi = polz.as_expr()
        else:
            L = 1
            coeffs_q = [sp.Rational(c) for c in polz.all_coeffs()]
            for cq in coeffs_q:
                L = int(sp.ilcm(L, int(cq.q)))
            # monic Z via x = w/L clearing: L^5 mon(x/L)
            chi = sp.expand(
                sum(
                    int(coeffs_q[i] * (L**i)) * x ** (5 - i)
                    for i in range(6)
                )
            )
            # leading = 1 * L^0 = 1? mon = x^5 + c4 x^4 + ... with ci rational
            # L^5 mon(x/L) = x^5 + c4 L x^4 + c3 L^2 x^3 + ... yes monic Z if L clears dens of c_i * L^{i-?}
            # c_i * L^i must be integer: L multiple of den(c_i) and we need den(c_i)|L^i — safer L = lcm dens
            chi = sp.expand(
                x**5
                + sum(
                    sp.Integer(sp.Rational(coeffs_q[i]) * (L**i)) * x ** (5 - i)
                    for i in range(1, 6)
                )
            )
        try:
            cl = classify_poly(chi, do_galois=True)
            samples.append(
                {
                    "t": tv,
                    "poly": str(chi)[:80],
                    "status": cl.get("status"),
                    "disc_sq": cl.get("disc_square"),
                    "galois": cl.get("galois"),
                }
            )
            print(f"    t={tv}: {cl.get('status')}", flush=True)
        except Exception as ex:
            samples.append({"t": tv, "error": str(ex)[:60]})

    # t=0 recovers seed
    recovers = sp.expand(mon.subs(t, 0) - (z**5 - 55 * z + 88)) == 0

    return {
        "seed": "x^5 - 55 x + 88",
        "R": str(FLAGSHIP_R),
        "construction": "Res_y(P(y), z - y - t R(y))",
        "coeffs_factored": factored,
        "closed_form_latex": (
            r"P_t(z)=z^5 - 385 t\, z^4 - 440 t(380 t+3)\, z^3 "
            r"+ 3520 t(18150 t^2 + 205 t + 3)\, z^2 "
            r"+ 55(45619200 t^4 - 4364800 t^3 - 21120 t^2 - 256 t - 1)\, z "
            r"+ 11(2269696000 t^5 - 444928000 t^4 + 21120000 t^3 + 70400 t^2 + 165 t + 8)"
        ),
        "closed_form_sympy": str(mon),
        "disc_factored_preview": str(Dfac)[:240],
        "disc_square_in_Qt": disc_square,
        "disc_info": {"num": info_n, "den": info_d},
        "t0_recovers_seed": recovers,
        "samples": samples,
        "n_A5": sum(1 for s in samples if str(s.get("status", "")).startswith("HIT_A5")),
    }


# ---------------------------------------------------------------------------
# 2. B-embed lattice: d=-75, e=f=0, a=-A, b*c = 72 A
# ---------------------------------------------------------------------------

LATTICE_GENS = sorted(
    set(MODEL_CORE.keys())
    | {1, 2, 4, 5, 6, 8, 12, 15, 16, 24, 25, 27, 36, 45, 48, 54, 55, 72, 88, 95, 100, 243, 539, 4880}
)


def idea_B_poly(A_val: int):
    return x**5 + 75 * x**3 + A_val * x**2 + 3 * A_val


def systematic_B_embed(max_A: int = 80, pairs_per_A: int = 6) -> dict:
    """
    Enumerate A from resonant lattice / short products;
    for each A pick a few b | 72A (prefer model divisors); classify Gal once per A.
    """
    print("  B-embed lattice scan...", flush=True)
    A_cands = set()
    for g in LATTICE_GENS:
        A_cands.add(g)
        A_cands.add(-g)
    for g1, g2 in itertools.product([3, 9, 27, 61, 80, 243], repeat=2):
        A_cands.add(g1 * g2)
        A_cands.add(-(g1 * g2))
    for g in [3, 9, 27, 61, 80, 55, 88]:
        A_cands.add(72 * g)
        A_cands.add(-72 * g)
    A_cands.discard(0)

    prefer_b = [1, 2, 3, 4, 6, 8, 9, 12, 18, 24, 27, 36, 61, 72, 80, 243, 539]

    def lattice_score(n: int) -> int:
        n = abs(int(n))
        if n == 0:
            return -1
        score = 0
        for g in sorted(MODEL_CORE.keys(), reverse=True):
            while n % g == 0 and g > 1:
                n //= g
                score += 2 if g in (3, 9, 27, 61, 80, 243, 539) else 1
        if n > 1:
            for p, mlt in sp.factorint(n).items():
                if p not in (2, 5):
                    score -= mlt
        return score

    def pairs_for_A(A_val: int):
        target = 72 * A_val
        out = []
        # preferred model divisors
        for bb in prefer_b + [abs(target)]:
            if bb == 0 or target % bb != 0:
                continue
            cc = target // bb
            out.append((bb, cc))
            out.append((-bb, -cc))
            if len(out) >= pairs_per_A * 2:
                break
        # a few more small positive divisors
        for bb in range(1, min(200, abs(target) + 1)):
            if target % bb == 0:
                cc = target // bb
                if (bb, cc) not in out:
                    out.append((bb, cc))
                if len(out) >= pairs_per_A * 3:
                    break
        # unique preserving order
        seen = set()
        uniq = []
        for p in out:
            if p not in seen and p[0] * p[1] == target:
                seen.add(p)
                uniq.append(p)
        return uniq[: pairs_per_A * 2]

    rows = []
    by_A = {}
    a5 = []
    n_gal_checks = 0
    # Prioritise model A
    ordered_A = sorted(
        A_cands,
        key=lambda u: (
            0 if abs(u) in MODEL_CORE else 1,
            0 if abs(u) in (3, 9, 27, 61, 80, 243, 539, 55, 88, 95) else 1,
            abs(u),
            u,
        ),
    )[:max_A]

    for A_val in ordered_A:
        expected = idea_B_poly(A_val)
        pol = sp.Poly(expected, x, domain=sp.ZZ)
        irr = bool(pol.is_irreducible)
        disc = int(pol.discriminant()) if irr else None
        sq = bool(irr and disc and disc > 0 and is_square(disc))
        status = "red"
        galois = None
        if not irr:
            status = "red"
        elif sq:
            # Gal check for model / primary lattice A
            if n_gal_checks < 50 and (
                abs(A_val) in MODEL_CORE
                or abs(A_val) in (3, 9, 27, 61, 80, 243, 539, 55, 88, 95, 4880, 18, 54)
                or lattice_score(A_val) >= 2
            ):
                cl = classify_poly(expected, do_galois=True)
                status = cl.get("status")
                galois = cl.get("galois")
                n_gal_checks += 1
            else:
                status = "disc_sq"
        else:
            status = "odd"

        pairs = pairs_for_A(A_val)
        best = None
        for bi, ci in pairs:
            rec = {
                "A": A_val,
                "b": bi,
                "c": ci,
                "a": -A_val,
                "d": -75,
                "e": 0,
                "f": 0,
                "T": f"T({-A_val},{bi},{ci},-75,0,0)",
                "disc_sq": sq,
                "irreducible": irr,
                "lattice_score_bc": lattice_score(bi) + lattice_score(ci),
                "lattice_score_A": lattice_score(A_val),
                "A_in_model": abs(A_val) in MODEL_CORE,
                "status": status,
                "galois": galois,
            }
            rows.append(rec)
            if best is None or rec["lattice_score_bc"] > best["lattice_score_bc"]:
                best = rec
        if best:
            by_A[A_val] = best
            if str(status).startswith("HIT_A5"):
                a5.append(best)
        print(f"    A={A_val}: irr={irr} disc□={sq} {status} pairs={len(pairs)}", flush=True)

    model_A5 = [
        r
        for r in a5
        if abs(r["A"]) in MODEL_CORE
        or abs(r["A"]) in (3, 9, 27, 61, 80, 243, 539, 4880, 55, 88, 95)
    ]

    return {
        "relation": "d=-75, e=f=0, a=-A, b*c=72*A",
        "chi": "x^5 + 75 x^3 + A x^2 + 3 A",
        "disc_identity": "324 * A**2 * (A**2 + 84375)**2",
        "n_points_scanned": len(rows),
        "n_unique_A": len(by_A),
        "n_disc_sq": sum(1 for r in by_A.values() if r.get("disc_sq")),
        "n_A5": len(a5),
        "n_gal_checks": n_gal_checks,
        "model_A5_sample": model_A5[:25],
        "best_per_A": sorted(by_A.values(), key=lambda r: (abs(r["A"]), r["A"]))[:50],
        "top_lattice_A5": sorted(a5, key=lambda r: -r.get("lattice_score_bc", 0))[:25],
        "all_A_status": [
            {
                "A": r["A"],
                "status": r.get("status"),
                "disc_sq": r.get("disc_sq"),
                "T": r.get("T"),
            }
            for r in sorted(by_A.values(), key=lambda r: (abs(r["A"]), r["A"]))
        ],
    }


# ---------------------------------------------------------------------------
# 3. Matrix avatar with built-in evenness
# ---------------------------------------------------------------------------


def evenness_avatar() -> dict:
    """
    Matrix model whose characteristic polynomial is *identically* pure-even BJ
    for free parameters (m,k) — evenness by construction, not by search.

    Construction:
      T(a,b,c,d,e,f) with d=0, a=-e f  (BJ-embed)
      Choose e=0, f=1, b=-α, c=k  where
        α = 256 m^2 - 3125 k^4 / 256,  β = k α
      Then chi = x^5 + α x + β with disc = (256 α^2 m)^2 identically.

    Integer/cleared form for lattice k:
      Use cleared α_Z = 65536 m^2 - 3125 k^4   wait careful scaling.
      Standard: α = 256 m^2 - 3125 k^4/256, β = k α.
      Matrix: T(0, -α, k, 0, 0, 1)  — entries may be rational.

    Cleared integer matrix for m = p/q, k fixed rational:
      Work with monic Z poly via homogenisation weights on specialisations.

    Also: non-BJ even avatar from Idea B:
      T(-A, b, 72A/b, -75, 0, 0) — disc identically square in A (for any b | 72A).
    """
    print("  evenness avatar identities...", flush=True)
    # --- Avatar PE: pure-even BJ matrix ---
    al = 256 * m**2 - sp.Rational(3125) * k**4 / 256
    be = k * al
    # Matrix T(0, -α, k, 0, 0, 1)
    chi_pe = (
        x**5
        - 0 * x**3
        - (0 + 0 * 1) * x**2
        - ((-al) * 1 + k * 0) * x
        + (0 * 0 - (-al) * k)
    )
    chi_pe = sp.expand(chi_pe)
    # Should be x^5 + α x + β = x^5 + al x + k*al
    expected = sp.expand(x**5 + al * x + be)
    match_pe = sp.expand(chi_pe - expected) == 0
    D_pe = sp.expand(256 * al**5 + 3125 * be**4)
    sqrt_pe = sp.expand(256 * al**2 * m)
    id_pe = sp.expand(D_pe - sqrt_pe**2) == 0

    # --- Avatar B: non-BJ identically even ---
    A = sp.symbols("A")
    bb = sp.symbols("b", nonzero=True)
    # c = 72 A / b
    chi_B = (
        x**5
        - (-75) * x**3
        - ((-A) + 0) * x**2
        - 0 * x
        + ((-A) * (-75) - bb * (72 * A / bb))
    )
    chi_B = sp.expand(chi_B)
    expected_B = sp.expand(x**5 + 75 * x**3 + A * x**2 + 3 * A)
    match_B = sp.simplify(chi_B - expected_B) == 0
    D_B = sp.expand(sp.Poly(expected_B, x).discriminant())
    # 324 A^2 (A^2+84375)^2
    sqrt_B = 18 * A * (A**2 + 84375)
    id_B = sp.expand(D_B - sqrt_B**2) == 0

    # Sample integer PE specialisations on HQCC lattice k
    pe_samples = []
    for kk, name in [
        (sp.Rational(-8, 5), "flagship"),
        (sp.Rational(4, 5), "classical"),
        (-4, "LSW"),
        (sp.Rational(-12, 5), "s12"),
    ]:
        for mv in [sp.Rational(1, 8), sp.Rational(5, 8), sp.Rational(5, 4), 1]:
            alv = al.subs({m: mv, k: kk})
            bev = be.subs({m: mv, k: kk})
            if alv == 0:
                continue
            # clear to Z monic
            alv, bev = sp.together(alv), sp.together(bev)
            try:
                A_num, A_den = sp.fraction(sp.together(alv))
                B_num, B_den = sp.fraction(sp.together(bev))
                L = sp.ilcm(int(sp.denom(sp.QQ(A_num) / sp.QQ(A_den))), int(sp.denom(sp.QQ(B_num) / sp.QQ(B_den))))
                # simpler: α,β as Rational
                ar, br = sp.Rational(alv), sp.Rational(bev)
                # poly x^5 + ar x + br; clear dens with homog weights
                Da, Db = int(ar.q), int(br.q)
                Dclear = int(sp.ilcm(Da, Db))
                # x = y: want D^5 ( (y/D)^5 + ar (y/D) + br ) = y^5 + ar D^4 y + br D^5
                AA = int(ar * Dclear**4)
                BB = int(br * Dclear**5)
                chi = x**5 + AA * x + BB
                pol = sp.Poly(chi, x, domain=sp.ZZ)
                if not pol.is_irreducible:
                    pe_samples.append(
                        {"k": str(kk), "name": name, "m": str(mv), "status": "red"}
                    )
                    continue
                disc = int(pol.discriminant())
                rec = {
                    "k": str(kk),
                    "name": name,
                    "m": str(mv),
                    "alpha": AA,
                    "beta": BB,
                    "disc_sq": is_square(disc),
                    "matrix": f"T(0,{-ar},{kk},0,0,1) cleared Z-poly",
                }
                if rec["disc_sq"] and len([s for s in pe_samples if str(s.get("status","")).startswith("HIT")]) < 12:
                    cl = classify_poly(chi, do_galois=True)
                    rec["status"] = cl.get("status")
                pe_samples.append(rec)
            except Exception as ex:
                pe_samples.append({"k": str(kk), "m": str(mv), "error": str(ex)[:50]})

    # B-avatar samples with lattice A and b|72A
    b_avatar_samples = []
    for Av in [3, 9, 27, 61, 80, 243, 539, 55, 88, -3, -9]:
        target = 72 * Av
        # pick b dividing target, prefer model
        chosen = None
        for bb in [1, 3, 8, 9, 24, 72, 61, 80, 27]:
            if target % bb == 0:
                chosen = (bb, target // bb)
                break
        if not chosen:
            continue
        bb, cc = chosen
        chi = idea_B_poly(Av)
        cl = classify_poly(chi, do_galois=True)
        b_avatar_samples.append(
            {
                "A": Av,
                "T": f"T({-Av},{bb},{cc},-75,0,0)",
                "status": cl.get("status"),
                "disc_sq": cl.get("disc_square"),
            }
        )

    return {
        "avatar_PE": {
            "name": "pure_even_BJ_matrix",
            "template": "T(0, -α, k, 0, 0, 1)",
            "alpha": "256*m**2 - 3125*k**4/256",
            "beta": "k*alpha",
            "chi_match_identity": match_pe,
            "disc_identity": id_pe,
            "disc_sqrt": "256 * alpha**2 * m",
            "evenness": "identical square by pure-even formula — built-in",
            "HQCC_native": False,
            "note": "Classical pure-even; HQCC only specialises (m,k)",
            "samples": pe_samples[:20],
        },
        "avatar_B": {
            "name": "nonBJ_deg1_matrix",
            "template": "T(-A, b, 72*A/b, -75, 0, 0)",
            "chi_match_identity": match_B,
            "disc_identity": id_B,
            "disc_sqrt": "18*A*(A**2+84375)",
            "evenness": "identical square in A for any b | 72A — built-in",
            "beyond_BJ_embed": True,
            "HQCC_native": False,
            "note": "d=-75≠0; beyond BJ-embed; lattice A specialises",
            "samples": b_avatar_samples,
        },
        "what_is_not_claimed": (
            "Neither avatar is forced by unrestricted ternary matrix axioms alone; "
            "evenness is built into the parameterisation (pure-even or B-family identity)."
        ),
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def main():
    t0 = time.time()
    print("NEXT: Mestre P_t + B-embed lattice + evenness avatar", flush=True)

    Pt = compute_flagship_Pt()
    Bemb = systematic_B_embed(180)
    av = evenness_avatar()

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"Next Mestre/B/avatar ({elapsed}s). "
        f"Flagship P_t closed form: disc□ in Q(t)={Pt['disc_square_in_Qt']}, "
        f"t=0 recovers seed={Pt['t0_recovers_seed']}, sample A5={Pt['n_A5']}. "
        f"B-embed lattice: points={Bemb['n_points_scanned']}, unique A={Bemb['n_unique_A']}, "
        f"disc□={Bemb['n_disc_sq']}, A5={Bemb['n_A5']}. "
        f"Avatars: PE identity={av['avatar_PE']['disc_identity']}, "
        f"B identity={av['avatar_B']['disc_identity']} (beyond BJ). "
        f"Evenness built-in by parameterisation, not Crit-2 forcing."
    )
    print(verdict, flush=True)

    # --- MESTRE_FLAGSHIP_PT.md ---
    lines_pt = [
        r"# Closed-form flagship Mestre family \(P_t\)",
        "",
        f"_Elapsed portion of next-run; total job {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        r"## Seed and Mestre \(R\)",
        "",
        r"$$P(x)=x^5-55x+88,\qquad R(x)=x^4+8x^3-32x^2+33.$$",
        "",
        r"Condition \(P''R-2P'R'\equiv 0\pmod{P}\) holds (1-dimensional solution space).",
        "",
        r"## Construction",
        "",
        r"$$P_t(z)=\operatorname{Res}_y\bigl(P(y),\, z-y-t R(y)\bigr)\in\mathbb{Q}(t)[z].$$",
        "",
        r"Monic of degree 5. At \(t=0\): recovers seed = "
        f"**{Pt['t0_recovers_seed']}**.",
        "",
        r"## Closed form",
        "",
        "$$",
        Pt["closed_form_latex"],
        "$$",
        "",
        r"### Coefficients (factored)",
        "",
        r"| coeff | factored form |",
        r"|-------|---------------|",
    ]
    for lab, fs in Pt["coeffs_factored"].items():
        lines_pt.append(f"| `{lab}` | `{fs}` |")

    lines_pt += [
        "",
        r"## Discriminant",
        "",
        f"- Identically square in \(\\mathbb{{Q}}(t)\)? **{Pt['disc_square_in_Qt']}** "
        f"(`{Pt['disc_info']}`)",
        f"- Factored preview: `{Pt['disc_factored_preview']}`",
        "",
        r"## Specialisations",
        "",
        r"| \(t\) | disc□ | status | Gal |",
        r"|----:|:-----:|--------|-----|",
    ]
    for s in Pt["samples"]:
        lines_pt.append(
            f"| {s.get('t')} | {s.get('disc_sq')} | {s.get('status')} | {s.get('galois')} |"
        )

    lines_pt += [
        "",
        r"## Role",
        "",
        r"- **Generative:** 1-param \(A_5\) family through the HQCC flagship seed.",
        r"- **Not necessity:** Mestre evenness is classical (seed disc□ \(\Rightarrow\) family disc□).",
        r"- HQCC enters by **choice of seed**, not by forcing the differential condition.",
        "",
        r"```bash",
        r"python next_mestre_b_avatar.py",
        r"```",
        "",
        r"_Generated by next_mestre_b_avatar.py_",
    ]

    # --- B_EMBED_LATTICE.md ---
    lines_b = [
        r"# B-embed lattice — systematic \(bc=72A\)",
        "",
        f"_Job elapsed {elapsed}s_",
        "",
        r"## Setup",
        "",
        r"Non-BJ family \(P_A=x^5+75x^3+A x^2+3A\) with",
        r"$$\operatorname{disc}(P_A)=324 A^2(A^2+84375)^2$$",
        r"(identically square).",
        "",
        r"**Matrix realisation** (beyond BJ-embed, \(d\neq 0\)):",
        r"$$d=-75,\quad e=f=0,\quad a=-A,\quad bc=72A.$$",
        r"Then \(\chi_T=P_A\).",
        "",
        f"- Points scanned: **{Bemb['n_points_scanned']}**",
        f"- Unique \(A\): **{Bemb['n_unique_A']}**",
        f"- disc□: **{Bemb['n_disc_sq']}**",
        f"- Gal \(A_5\) among checked: **{Bemb['n_A5']}**",
        "",
        r"## Model / resonant \(A\) with \(A_5\)",
        "",
        r"| \(A\) | \(b\) | \(c\) | \(T\) | status |",
        r"|----:|----:|----:|------|--------|",
    ]
    for r in Bemb["model_A5_sample"]:
        lines_b.append(
            f"| {r['A']} | {r['b']} | {r['c']} | `{r.get('T')}` | {r.get('status')} |"
        )
    if not Bemb["model_A5_sample"]:
        lines_b.append(r"| _(see top lattice A5)_ | | | | |")

    lines_b += [
        "",
        r"## Top lattice-score \(A_5\) points",
        "",
        r"| \(A\) | \(b\) | \(c\) | score | status |",
        r"|----:|----:|----:|------:|--------|",
    ]
    for r in Bemb["top_lattice_A5"]:
        lines_b.append(
            f"| {r['A']} | {r['b']} | {r['c']} | {r.get('lattice_score_bc')} | {r.get('status')} |"
        )

    lines_b += [
        "",
        r"## Best pair per \(A\) (sample)",
        "",
        r"| \(A\) | \(b\) | \(c\) | disc□ | status |",
        r"|----:|----:|----:|:-----:|--------|",
    ]
    for r in Bemb["best_per_A"][:30]:
        lines_b.append(
            f"| {r['A']} | {r['b']} | {r['c']} | {r.get('disc_sq')} | {r.get('status')} |"
        )

    lines_b += [
        "",
        r"## Notes",
        "",
        r"- Every \(A\neq 0\) has disc□ by identity; Gal \(A_5\) when irr + type (3,1,1).",
        r"- Integer \(T\) entries require \(b\mid 72A\). Prefer \(b\) from "
        r"\(\{1,3,8,9,24,72,27,61,80,\ldots\}\).",
        r"- **HQCC-native?** Parameter \(A\) and factors of \(72A\) can be lattice integers; "
        r"the relation \(bc=72A\) and \(d=-75\) are still classical ansatz, not Crit-2 forcing.",
        "",
        r"_Generated by next_mestre_b_avatar.py_",
    ]

    # --- EVENNESS_AVATAR.md ---
    pe, ba = av["avatar_PE"], av["avatar_B"]
    lines_av = [
        r"# Matrix avatars with built-in evenness identity",
        "",
        f"_Job elapsed {elapsed}s_",
        "",
        r"**Design rule:** only admit matrix models whose disc(\(\chi\)) is a "
        r"**square polynomial identity** in free parameters — not sparse search.",
        "",
        "---",
        "",
        r"## Avatar PE — pure-even BJ matrix",
        "",
        r"$$T(0,-\alpha,k,0,0,1),\qquad "
        r"\alpha=256m^2-\frac{3125 k^4}{256},\quad \beta=k\alpha.$$",
        "",
        f"- \(\\chi\) matches \(x^5+\\alpha x+\\beta\)? **{pe['chi_match_identity']}**",
        f"- disc \(=(256\\alpha^2 m)^2\)? **{pe['disc_identity']}**",
        f"- HQCC-native forcing? **{pe['HQCC_native']}** — {pe['note']}",
        "",
        r"| \(k\) | name | \(m\) | disc□ | status |",
        r"|------|------|------|:-----:|--------|",
    ]
    for s in pe["samples"]:
        lines_av.append(
            f"| {s.get('k')} | {s.get('name')} | {s.get('m')} | {s.get('disc_sq')} | {s.get('status')} |"
        )

    lines_av += [
        "",
        "---",
        "",
        r"## Avatar B — non-BJ degree-1 matrix (beyond BJ-embed)",
        "",
        r"$$T(-A,\, b,\, 72A/b,\, -75,\, 0,\, 0)\qquad (b\mid 72A).$$",
        "",
        f"- \(\\chi\) matches \(x^5+75x^3+A x^2+3A\)? **{ba['chi_match_identity']}**",
        f"- disc \(=(18A(A^2+84375))^2\)? **{ba['disc_identity']}**",
        f"- Beyond BJ-embed (\(d\\neq 0\))? **{ba['beyond_BJ_embed']}**",
        f"- HQCC-native forcing? **{ba['HQCC_native']}** — {ba['note']}",
        "",
        r"| \(A\) | \(T\) | disc□ | status |",
        r"|----:|------|:-----:|--------|",
    ]
    for s in ba["samples"]:
        lines_av.append(
            f"| {s['A']} | `{s['T']}` | {s.get('disc_sq')} | {s.get('status')} |"
        )

    lines_av += [
        "",
        "---",
        "",
        r"## What this is / is not",
        "",
        av["what_is_not_claimed"],
        "",
        r"| Claim | Verdict |",
        r"|-------|---------|",
        r"| Built-in disc□ identity | **Yes** (both avatars) |",
        r"| Hosts HQCC lattice specialisations | **Yes** |",
        r"| Unrestricted ternary matrix \(\Rightarrow\) even | **No** (old Crit 2) |",
        r"| Necessity theorem | **No** — evenness parameterised in |",
        "",
        r"_Generated by next_mestre_b_avatar.py_",
    ]

    # Combined
    lines = [
        r"# Next: flagship \(P_t\), B-embed lattice, evenness avatars",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        r"Detail docs: `MESTRE_FLAGSHIP_PT.md`, `B_EMBED_LATTICE.md`, `EVENNESS_AVATAR.md`.",
        "",
        "---",
        "",
        r"## 1. Flagship \(P_t\) (closed form)",
        "",
        "$$",
        Pt["closed_form_latex"],
        "$$",
        "",
        f"- Construction: `{Pt['construction']}` with \(R={Pt['R']}\)",
        f"- disc□ in \(\\mathbb{{Q}}(t)\): **{Pt['disc_square_in_Qt']}**",
        f"- \(t=0\) recovers seed: **{Pt['t0_recovers_seed']}**",
        f"- Sample \(A_5\): **{Pt['n_A5']}** / {len(Pt['samples'])}",
        "",
        r"## 2. B-embed lattice \(bc=72A\)",
        "",
        f"- Scanned points: **{Bemb['n_points_scanned']}**, unique \(A\): **{Bemb['n_unique_A']}**",
        f"- disc□: **{Bemb['n_disc_sq']}**, \(A_5\): **{Bemb['n_A5']}**",
        f"- Model A5 sample size: **{len(Bemb['model_A5_sample'])}**",
        "",
        r"## 3. Evenness avatars",
        "",
        f"- PE matrix identity: **{pe['disc_identity']}**",
        f"- B matrix identity (beyond BJ): **{ba['disc_identity']}**",
        "",
        r"## Synthesis",
        "",
        r"1. Flagship has an explicit Mestre 1-param \(A_5\) family in closed form.",
        r"2. B-family supplies a systematic lattice of integer templates \(T(-A,b,72A/b,-75,0,0)\).",
        r"3. Evenness avatars encode disc□ **by construction**; they package generative success, "
        r"not Crit-2 necessity.",
        "",
        r"```bash",
        r"python next_mestre_b_avatar.py",
        r"```",
        "",
        r"_Generated by next_mestre_b_avatar.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "flagship_Pt": Pt,
        "B_embed": Bemb,
        "avatars": av,
    }

    for name, text in [
        ("MESTRE_FLAGSHIP_PT.md", "\n".join(lines_pt)),
        ("B_EMBED_LATTICE.md", "\n".join(lines_b)),
        ("EVENNESS_AVATAR.md", "\n".join(lines_av)),
        ("NEXT_MESTRE_B_AVATAR.md", "\n".join(lines)),
    ]:
        write_md(ROOT / name, text)
        write_md(OUT / name, text)
        try:
            if RESULTS.exists():
                write_md(RESULTS / name, text)
        except Exception:
            pass

    write_json(ROOT / "NEXT_MESTRE_B_AVATAR.json", payload)
    write_json(OUT / "NEXT_MESTRE_B_AVATAR.json", payload)
    print(f"Wrote MESTRE_FLAGSHIP_PT / B_EMBED_LATTICE / EVENNESS_AVATAR ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
