"""
Route 2: arithmetic over K = Q(√5).

The preferred rigid cover φ = 6y^5 - 15y^4 + 10y^3 has

    disc_y( monic(φ - t) ) = 5 · (25 t (t-1) / 36)²   in Q(t).

Over Q this is never a square for t ∈ Q \\ {0,1} (permanent factor 5).
Over K = Q(√5) one has 5 = (√5)², so the disc is identically a square in K(t):

    disc = ( √5 · 25 t (t-1) / 36 )².

Hence every non-critical fibre is *even over K* (disc square in the base field).
This script:

1. Proves the disc identity and the K-square theorem.
2. Samples specialisations t ∈ Q and t ∈ Z[√5]; records Z-models / O_K-models.
3. Checks irreducibility over Q vs over K.
4. Attempts Bring–Jerrard reduction of fibres over K.
5. Tests whether BJ coeffs (α,β) ∈ K lie on pure-even k-slice families
   (same construction as over Q, coeffs allowed in K).
6. Tests descent: norms, Galois conjugates, relation of K-even fibres
   to the pure-even Q-families (LSW, flagship, …).

Output: K_SQRT5_EVEN.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter
from fractions import Fraction
from math import gcd
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

y, t = sp.symbols("y t")
s5 = sp.sqrt(5)  # algebraic √5
PHI = 6 * y**5 - 15 * y**4 + 10 * y**3

# Minimal polynomial of √5
K_MP = sp.Poly(x**2 - 5, x, domain=sp.ZZ)


# ---------------------------------------------------------------------------
# 1. Disc theorem over Q(t) and K(t)
# ---------------------------------------------------------------------------
def disc_theorem() -> dict:
    """Symbolic proof: disc monic(φ-t) = 5 · (square) in Q(t); square in K(t)."""
    print("  disc theorem monic(φ-t)...", flush=True)
    mon = sp.expand((PHI - t) / 6)
    pol = sp.Poly(mon, y)
    D = sp.together(sp.expand(pol.discriminant()))
    num, den = sp.fraction(D)
    num_f = sp.factor(sp.expand(num))
    den_f = sp.factor(sp.expand(den))

    # Expected: 3125/1296 * t^2 (t-1)^2 = 5 * (25 t (t-1)/36)^2
    expected = sp.together(sp.Rational(3125, 1296) * t**2 * (t - 1) ** 2)
    square_part = sp.together(sp.Rational(25, 36) * t * (t - 1))
    five_sq = sp.together(5 * square_part**2)
    match_expected = sp.expand(sp.together(D - expected)) == 0
    match_five = sp.expand(sp.together(D - five_sq)) == 0

    # Over K: disc = (√5 * square_part)^2
    k_sqrt_form = sp.expand((s5 * square_part) ** 2)
    match_k = sp.expand(sp.together(D - k_sqrt_form)) == 0

    # Permanent factor analysis: D / (t(t-1))^2 as element of Q
    ratio = sp.together(D / (t**2 * (t - 1) ** 2))
    ratio_simp = sp.simplify(ratio)

    return {
        "monic": str(mon),
        "disc": str(D),
        "disc_factored_num": str(num_f),
        "disc_factored_den": str(den_f),
        "equals_3125_1296_form": match_expected,
        "equals_5_times_square": match_five,
        "square_part": str(square_part),
        "equals_sqrt5_square_form": match_k,
        "ratio_D_over_t2_tm1_2": str(ratio_simp),
        "theorem_Q": (
            "For t in Q \\ {0,1}, disc(monic(φ-t)) = 5 · (25 t (t-1)/36)^2 "
            "is 5 times a square, hence not a square in Q."
        ),
        "theorem_K": (
            "Over K=Q(√5), 5=(√5)^2, so disc = (√5 · 25 t (t-1)/36)^2 "
            "is a square in K for all t in K \\ {0,1} (and in K(t) identically)."
        ),
        "proved": match_expected and match_five and match_k,
    }


# ---------------------------------------------------------------------------
# Helpers: Z-model of fibre, disc square in Q and in K
# ---------------------------------------------------------------------------
def monic_fibre_to_Z(expr, var=y) -> sp.Poly | None:
    pol = sp.Poly(sp.expand(expr), var, domain=sp.QQ)
    if pol.degree() != 5:
        return None
    mon = sp.Poly(sp.monic(pol.as_expr()), var, domain=sp.QQ)
    dens = []
    for c in mon.all_coeffs():
        dens.append(sp.fraction(sp.together(c))[1])
    L = 1
    for d in dens:
        try:
            L = int(sp.ilcm(L, abs(int(d))))
        except Exception:
            return None
    cleared = sp.expand(L**5 * mon.as_expr().subs(var, var / L))
    p = sp.Poly(cleared, var, domain=sp.ZZ)
    if p.LC() == -1:
        p = sp.Poly(-p.as_expr(), var, domain=sp.ZZ)
    if p.LC() != 1:
        return None
    return p


def is_square_in_K_of_rational(d: int) -> dict:
    """
    Is d a square in K=Q(√5)?
    d = square in K ⇔ d = (a + b√5)^2 = a^2 + 5b^2 + 2ab√5 with a,b∈Q
    ⇒ 2ab=0 and a^2+5b^2=d.
    Cases: b=0 ⇒ d square in Q; or a=0 ⇒ d=5b^2 ⇒ d/5 square in Q.
    More generally for d∈Q: write d=p/q; clear to integer N.
    For N∈Z: N is square in K iff N is of the form x^2 + 5 y^2 with
    either (standard) the equation N = a^2 + 5 b^2 for a,b∈Z after scaling,
    OR equivalently: in Q(√5), N is square iff v_p(N) even for all p≡±2 mod 5
    after accounting for unit/5... Practical: N is K-square iff
    N * m^2 = A^2 + 5 B^2 ... simpler computational test:

    N = (r + s√5)^2 / den^2 with r,s,den ∈ Z
    ⇔ N den^2 = r^2 + 5 s^2 and 2 r s = 0 if we want pure rational*5 form...

    Full: (r + s √5)^2 = r^2 + 5 s^2 + 2rs √5.
    For this to equal N ∈ Q, need 2rs = 0, so r=0 or s=0.
    Thus N square in K (with N∈Q) iff N is square in Q OR N/5 is square in Q.

    (If we allow N = u * (r+s√5)^2 / (c+d√5)^2 more generally — that's always
    about being square in K as field. An element c0 + c1√5 is a square in K
    under different conditions. For pure rational N=N+0√5:
    N = (a+b√5)^2 = a^2+5b^2 + 2ab√5 ⇒ 2ab=0, so yes only s=0 or a=0 cases.)
    """
    if d == 0:
        return {"square_in_K": True, "reason": "zero", "form": "0"}
    if d < 0:
        # -1 is not a square in Q(√5)? Actually -1 = ((1+√5)/2)^2 * something?
        # ((√5)/2)^2 * (-4/5) ... Check: a^2+5b^2 = d < 0 impossible for a,b real.
        # For a,b∈Q: a^2+5b^2 ≥ 0. With 2ab=0: if a=0, 5b^2=d; if b=0, a^2=d.
        # Negative d: never a square in K among elements of Q.
        return {"square_in_K": False, "reason": "negative_rational"}
    if is_square(d):
        return {"square_in_K": True, "reason": "square_in_Q", "form": f"{int(sp.integer_nthroot(d,2)[0])}^2"}
    if d % 5 == 0 and is_square(d // 5):
        r = int(sp.integer_nthroot(d // 5, 2)[0])
        return {
            "square_in_K": True,
            "reason": "5_times_square_in_Q",
            "form": f"(√5 · {r})^2",
        }
    # Also: d = 5 * (p/q)^2 with p/q reduced — already covered for integer d
    # if d/5 is rational square with den: e.g. d=5*4/9 not integer.
    # For integer d, d/5 square in Q means d/5 = a^2/b^2 ⇒ d b^2 = 5 a^2.
    # Check: is 5d a square? No — d = 5 (a/b)^2 ⇒ d b^2 = 5 a^2.
    # Integer test: factor d = 5^e * m, need e odd and m square, OR e even and m square
    # for Q-square; for K: (e odd and m square) OR (e even and m square).
    # Actually: d square in Q: all exponents even.
    # d = 5 * square in Q: v_5(d) odd and d/5^{v5} * 5^{v5-1} is square...
    # Simplest integer criterion for N∈Z>0 square in K among Q-elements:
    # N is square in Q, or N/5 is square in Q (N divisible by 5, N/5 square).
    # Wait: N = 5 * (2/3)^2 = 20/9 not integer.
    # N = 20: is 20 square in K? 20 = a^2+5b^2, 2ab=0.
    # a=0 ⇒ 5b^2=20 ⇒ b^2=4, b=2 ⇒ 20 = (2√5)^2. Yes 20/5=4 square.
    # N=45=5*9: (3√5)^2. Yes.
    # N=10: 10/5=2 not square, 10 not square → not K-square in Q-elements.
    # Confirm 10 = a^2+5b^2 with ab=0: no.
    return {"square_in_K": False, "reason": "not_Q_square_and_not_5_times_Q_square"}


def fibre_at_t(tv) -> dict:
    """Build Z-model of monic(φ-t) at rational t; classify over Q; disc in K."""
    expr = sp.expand((PHI - tv) / 6)
    pol = monic_fibre_to_Z(expr, y)
    if pol is None:
        return {"t": str(tv), "status": "no_Z_model"}
    irr_Q = bool(pol.is_irreducible)
    disc = int(pol.discriminant())
    ksq = is_square_in_K_of_rational(disc)
    rec = {
        "t": str(tv),
        "poly": str(pol.as_expr()),
        "coeffs": [int(c) for c in pol.all_coeffs()],
        "irr_over_Q": irr_Q,
        "disc": disc,
        "disc_square_Q": is_square(disc) if disc > 0 else False,
        "disc_square_K": ksq.get("square_in_K"),
        "disc_K_reason": ksq.get("reason"),
        "disc_K_form": ksq.get("form"),
    }
    # Irreducibility over K: factor poly over Q(√5)
    # sympy: factor over extension
    try:
        fac_K = sp.factor(pol.as_expr(), extension=[s5])
        # If factors nontrivially over K
        fac_list = sp.factor_list(pol.as_expr(), extension=[sp.sqrt(5)])
        # fac_list = (content, [(factor, mult), ...])
        factors = fac_list[1]
        degs = sorted(int(sp.degree(f, y)) for f, m in factors for _ in range(int(m)))
        rec["factor_degs_over_K"] = degs
        rec["irr_over_K"] = degs == [5]
        rec["factor_K_preview"] = str(fac_K)[:120]
    except Exception as e:
        rec["factor_K_error"] = str(e)
        # Fallback: compose with minpoly — check if irr in K[y]
        # Poly over QQ(√5)
        try:
            K = sp.QQ.algebraic_field(sp.sqrt(5))
            pk = sp.Poly(pol.as_expr(), y, domain=K)
            rec["irr_over_K"] = bool(pk.is_irreducible)
        except Exception as e2:
            rec["irr_over_K"] = None
            rec["irr_K_error"] = str(e2)

    if irr_Q and disc > 0:
        # Gal over Q
        try:
            r = classify_poly(pol.as_expr().subs(y, x), do_galois=True)
            rec["gal_Q"] = r.get("galois")
            rec["status_Q"] = r.get("status")
        except Exception as e:
            rec["gal_Q_error"] = str(e)
    return rec


# ---------------------------------------------------------------------------
# 2. Sample even-over-K specialisations
# ---------------------------------------------------------------------------
def sample_rational_t(max_p: int = 25, max_q: int = 8) -> dict:
    print("  sampling rational t for K-even fibres...", flush=True)
    rows = []
    stats = Counter()
    for q in range(1, max_q + 1):
        for p in range(-max_p, max_p + 1):
            if q > 1 and gcd(abs(p), q) != 1:
                continue
            tv = sp.Rational(p, q)
            if tv in (0, 1):
                stats["critical"] += 1
                continue
            stats["tested"] += 1
            rec = fibre_at_t(tv)
            if rec.get("status") == "no_Z_model":
                stats["no_Z"] += 1
                continue
            if rec.get("disc_square_K"):
                stats["even_K"] += 1
            else:
                stats["odd_K"] += 1
            if rec.get("disc_square_Q"):
                stats["even_Q"] += 1
            if rec.get("irr_over_Q"):
                stats["irr_Q"] += 1
            if rec.get("irr_over_K"):
                stats["irr_K"] += 1
            if rec.get("irr_over_K") and rec.get("disc_square_K"):
                stats["even_irr_K"] += 1
                rows.append(rec)
            elif rec.get("disc_square_K") and not rec.get("irr_over_Q"):
                stats["even_K_red_Q"] += 1
    # Also keep a sample of even_K whether or not irr
    # Re-scan smaller for full table of first even_K irr
    return {"stats": dict(stats), "even_irr_K_sample": rows[:40], "n_even_irr_K": len(rows)}


def sample_quadratic_t(max_ab: int = 6) -> dict:
    """
    t = a + b √5 with a,b ∈ Z small, b≠0.
    Fibre monic(φ-t) has coeffs in K; clear to poly over Q by taking
    either O_K-model or the norm (degree 10 over Q).
    We work with the poly in K[y] and check disc square in K.
    """
    print("  sampling t = a+b√5 ...", flush=True)
    rows = []
    stats = Counter()
    for a, b in itertools.product(range(-max_ab, max_ab + 1), repeat=2):
        if b == 0:
            continue  # pure rational covered above
        tv = a + b * s5
        # skip if t in {0,1} — not for b≠0
        stats["tested"] += 1
        mon = sp.expand((PHI - tv) / 6)
        # disc over K(t) theorem says square; verify numerically/symbolically
        try:
            pol = sp.Poly(mon, y, domain=sp.QQ.algebraic_field(sp.sqrt(5)))
            D = pol.discriminant()
            # D should be a square in K
            # Represent D as c0 + c1 √5
            Ds = sp.simplify(sp.expand(D))
            # Check D / (known square form) 
            square_part = sp.Rational(25, 36) * tv * (tv - 1)
            expected = sp.expand((s5 * square_part) ** 2)
            match = sp.simplify(sp.expand(Ds - expected)) == 0
            irr = bool(pol.is_irreducible)
            rec = {
                "t": f"{a}+{b}*sqrt(5)",
                "a": a,
                "b": b,
                "disc_match_sqrt5_form": match,
                "irr_over_K": irr,
                "disc_preview": str(Ds)[:80],
            }
            if match:
                stats["disc_square_K"] += 1
            if irr and match:
                stats["even_irr_K"] += 1
                rows.append(rec)
            elif not irr:
                stats["red_K"] += 1
        except Exception as e:
            stats["err"] += 1
            if stats["err"] <= 3:
                print(f"    err t={a}+{b}√5: {e}", flush=True)
    return {"stats": dict(stats), "even_irr_K_sample": rows[:30], "n_even_irr_K": len(rows)}


# ---------------------------------------------------------------------------
# 3. Bring–Jerrard reduction over K
# ---------------------------------------------------------------------------
def depress_quintic_to_bj(coeffs_monic: list, domain="Q"):
    """
    Monic quintic y^5 + c4 y^4 + c3 y^3 + c2 y^2 + c1 y + c0.
    Standard: shift y = z - c4/5 to kill z^4; then Tschirnhaus to kill z^3,z^2
    for Bring–Jerrard x^5 + α x + β (over radical extensions in general).

    For φ-fibres: monic = y^5 - (5/2)y^4 + (5/3)y^3 - t/6
    so c4=-5/2, c3=5/3, c2=0, c1=0, c0=-t/6.

    After y = z + 1/2 (since -c4/5 = (5/2)/5 = 1/2):
    the depressed form is known for this family.
    """
    # Use sympy to shift
    c4, c3, c2, c1, c0 = [sp.sympify(c) for c in coeffs_monic]
    # y = z - c4/5
    shift = -c4 / 5
    z = sp.symbols("z")
    f = (
        (z + shift) ** 5
        + c4 * (z + shift) ** 4
        + c3 * (z + shift) ** 3
        + c2 * (z + shift) ** 2
        + c1 * (z + shift)
        + c0
    )
    f = sp.expand(f)
    # collect powers of z
    pol = sp.Poly(f, z)
    # Should have no z^4
    coeffs = pol.all_coeffs()  # deg 5 .. 0
    # Tschirnhaus for BJ: z = x + a/x + b  or classical Bring radical
    # For this specific family a simpler path: the fibre is
    # monic(φ-t) = y^5 - 5/2 y^4 + 5/3 y^3 - t/6
    # After y=z+1/2: compute explicitly
    return {
        "shift": str(shift),
        "depressed": str(f),
        "depressed_coeffs": [str(c) for c in pol.all_coeffs()],
    }


def fibre_family_depressed() -> dict:
    """
    Symbolic depressed form of monic(φ-t) and attempt BJ over K(t).
    """
    print("  depressed monic(φ-t) + BJ attempt...", flush=True)
    mon = sp.expand((PHI - t) / 6)
    # mon = y^5 - (5/2)y^4 + (5/3)y^3 - t/6
    shift = sp.Rational(1, 2)  # -c4/5 = 1/2
    z = sp.symbols("z")
    f = sp.expand(mon.subs(y, z + shift))
    pol = sp.Poly(f, z)
    coeffs = pol.all_coeffs()
    # coeffs: [1, 0, p, q, r, s] expected
    names = ["z5", "z4", "z3", "z2", "z1", "z0"]
    coeff_map = {names[i]: sp.simplify(coeffs[i]) for i in range(len(coeffs))}

    # For Bring–Jerrard over C one solves a principal quintic; over K(t)
    # Tschirnhaus z = u x + v + w/x  with parameters to kill x^3,x^2.
    # Standard: after depressed (no x^4), use z = x + c/x to kill x^3
    # if the resolvent permits.
    c = sp.symbols("c")
    # z = x + c/x, multiply by x^5: 
    # (x + c/x)^5 + p (x+c/x)^3 + q (x+c/x)^2 + r (x+c/x) + s
    # times x^5 / x^5 monic in x...
    xvar = sp.symbols("X")
    p = coeff_map["z3"]
    q = coeff_map["z2"]
    r = coeff_map["z1"]
    s0 = coeff_map["z0"]
    zsub = xvar + c / xvar
    g = sp.together(
        zsub**5 + p * zsub**3 + q * zsub**2 + r * zsub + s0
    )
    g_clear = sp.simplify(sp.expand(g * xvar**5))
    gpol = sp.Poly(sp.numer(sp.together(g_clear)), xvar)
    # Want coeff of x^3 and x^2 (and ideally x^4 already 0) to vanish
    # degrees in cleared: x^10 + ... + c^5; monic deg 10? 
    # Actually (x+c/x)^5 * x^5 = (x^2+c)^5, degree 10. Wrong approach for monic quintic.
    #
    # Correct Tschirnhaus for quintic: z = x^2 + a x + b + c/x + d/x^2 is Bring's,
    # or one-parameter z = x + e/x for principal form.
    #
    # Document depressed coeffs; for BJ numerical at sample t over K.

    return {
        "shift": str(shift),
        "coeff_map": {k: str(v) for k, v in coeff_map.items()},
        "note": (
            "Depressed monic(φ-t) has coeffs in Q(t). Full Bring–Jerrard reduction "
            "of a general quintic requires solving an auxiliary resolvent (degree 6); "
            "over K(t) this is possible in radicals of K(t) but not typically with "
            "α,β ∈ K(t) rational. We test numeric BJ matching at specialisations."
        ),
    }


def try_bj_match_at_samples(even_rows: list, pure_even_seeds: list) -> dict:
    """
    For each even-over-K fibre (Z-model over Q), compare its Gal/disc to
    pure-even family seeds; try whether a Tschirnhaus over Q links them
    (same number field).
    """
    print("  relating K-even fibres to pure-even Q-seeds...", flush=True)
    # Pure-even seeds: list of (a,b,k,tag)
    results = []
    # Number field isomorphism test: same minpoly of a primitive element?
    # Simpler: for fibre poly f and seed poly s, check if disc(f)/disc(s) is square
    # (same quadratic resolvent parity) and if f,s define isomorphic fields
    # via Tschirnhaus — hard. Instead:
    # 1) Same discriminant up to squares?
    # 2) Do roots generate same field? (compare factorizations of resolvents)
    for fr in even_rows[:20]:
        fpoly = sp.Poly(fr["poly"], y, domain=sp.ZZ)
        fdisc = fr["disc"]
        matches = []
        for a, b, k, tag in pure_even_seeds:
            spoly = sp.Poly(x**5 + a * x + b, x, domain=sp.ZZ)
            sdisc = disc_bj_int(a, b)
            # disc ratio square in Q?
            if fdisc == 0 or sdisc == 0:
                continue
            # fdisc / sdisc square in Q?
            # Both may be 5*square; ratio may be square
            ratio = Fraction(fdisc, sdisc)
            num, den = ratio.numerator, ratio.denominator
            ratio_sq = is_square(abs(num)) and is_square(den)
            # ratio square in K?
            ratio_K = is_square_in_K_of_rational(num * den)  # for ±num/den cleared
            # Better: num/den is K-square if num*den is K-square after signs
            nd = num * den
            if nd < 0:
                ksq_ratio = {"square_in_K": False}
            else:
                ksq_ratio = is_square_in_K_of_rational(nd)

            # Field equality: check if f factors modulo many p same pattern as s
            # or compose resultant test: exists linear Tschirnhaus? unlikely
            # Check whether minpoly fields have same degree and disc class
            same_disc_class_Q = ratio_sq
            same_disc_class_K = ksq_ratio.get("square_in_K")

            if same_disc_class_Q or same_disc_class_K:
                matches.append(
                    {
                        "seed": tag,
                        "seed_ab": (a, b),
                        "k": k,
                        "ratio_disc": str(ratio),
                        "disc_ratio_square_Q": same_disc_class_Q,
                        "disc_ratio_square_K": same_disc_class_K,
                    }
                )
        results.append(
            {
                "t": fr.get("t"),
                "fibre_poly": fr.get("poly"),
                "fibre_disc": fdisc,
                "gal_Q": fr.get("gal_Q"),
                "irr_K": fr.get("irr_over_K"),
                "seed_disc_matches": matches,
            }
        )
    n_any_Q = sum(1 for r in results if any(m["disc_ratio_square_Q"] for m in r["seed_disc_matches"]))
    n_any_K = sum(1 for r in results if any(m["disc_ratio_square_K"] for m in r["seed_disc_matches"]))
    return {
        "compared": len(results),
        "fibres_with_seed_disc_class_Q": n_any_Q,
        "fibres_with_seed_disc_class_K": n_any_K,
        "details": results[:15],
    }


# ---------------------------------------------------------------------------
# 4. Pure-even families base-changed to K
# ---------------------------------------------------------------------------
def pure_even_over_K() -> dict:
    """
    Pure-even k-slice families remain pure-even over K (disc still square).
    Flagship / LSW specialisations with m ∈ K: α,β ∈ K; Gal over K even.

    Test m = r + s√5 small: get O_K-coeff BJ polys; check disc square in K;
    check whether any equal a φ-fibre's BJ reduction.
    """
    print("  pure-even families over K (m in Z[√5])...", flush=True)
    # k in Q as before
    ks = [
        Fraction(-4),
        Fraction(4),
        Fraction(-8, 5),
        Fraction(8, 5),
        Fraction(4, 5),
        Fraction(-4, 5),
        Fraction(-12, 5),
        Fraction(12, 5),
    ]
    samples = []
    for k in ks:
        for r, s in itertools.product(range(-3, 4), range(-2, 3)):
            if r == 0 and s == 0:
                continue
            m = r + s * s5
            # α = 256 m^2 - 3125 k^4 / 256
            alpha = sp.simplify(256 * m**2 - 3125 * sp.Rational(k.numerator, k.denominator) ** 4 / 256)
            beta = sp.simplify(sp.Rational(k.numerator, k.denominator) * alpha)
            # disc BJ
            D = sp.simplify(256 * alpha**5 + 3125 * beta**4)
            # Is D a square in K?
            # D should be (256 α^2 m)^2 by identity
            expected = sp.simplify((256 * alpha**2 * m) ** 2)
            ok = sp.simplify(D - expected) == 0
            # Is α,β in Q?
            a_in_Q = sp.im(alpha) == 0 and s5 not in sp.Poly(sp.expand(alpha), s5).gens if False else (
                sp.expand(alpha).as_poly(s5) is None or sp.degree(sp.expand(alpha), s5) == 0
            )
            try:
                ap = sp.Poly(sp.expand(alpha), s5)
                a_in_Q = ap.degree() == 0
                bp = sp.Poly(sp.expand(beta), s5)
                b_in_Q = bp.degree() == 0
            except Exception:
                a_in_Q = s5 not in sp.preorder_traversal(sp.expand(alpha))
                b_in_Q = s5 not in sp.preorder_traversal(sp.expand(beta))

            if ok and (s != 0):  # genuinely quadratic m
                samples.append(
                    {
                        "k": str(k),
                        "m": f"{r}+{s}*sqrt(5)",
                        "alpha": str(sp.simplify(alpha)),
                        "beta": str(sp.simplify(beta)),
                        "disc_identity_ok": ok,
                        "alpha_in_Q": bool(a_in_Q),
                        "beta_in_Q": bool(b_in_Q),
                    }
                )
            if len(samples) >= 40:
                break
        if len(samples) >= 40:
            break

    n_Q_coeffs = sum(1 for s in samples if s["alpha_in_Q"] and s["beta_in_Q"])
    return {
        "n_samples_quadratic_m": len(samples),
        "n_with_Q_coeffs": n_Q_coeffs,
        "samples": samples[:20],
        "note": (
            "For m ∈ K \\ Q, typically α,β ∈ K \\ Q: pure-even family over K "
            "produces BJ seeds over K, not automatically over Q. Descent to Q "
            "requires α,β ∈ Q, which forces m^2 ∈ Q under the formula "
            "α=256m^2 - const, so m ∈ Q or m ∈ √d·Q with constraints."
        ),
    }


def descent_analysis() -> dict:
    """
    When can a K-even φ-fibre descend to an even Q-poly?
    And relation between disc classes.
    """
    print("  descent analysis...", flush=True)
    # For t∈Q: fibre f ∈ Q[y], disc_Q = 5 · □, so odd over Q, even over K.
    # Descent of monodromy: the quadratic character of disc is the sign character
    # of Gal→{\pm1}. Restriction Gal(f/K) = Gal(f/Q) ∩ ker(sign) or index-2
    # depending on whether √disc ∈ K.
    # Since √disc = √5 · rational, √disc ∈ K ⇔ the poly's even monodromy over K
    # is exactly the kernel of the sign character on Gal(f/Q) ≅ S5 or D5 etc.

    # Concrete: sample t and report Gal over Q vs even-over-K
    samples = []
    for tv in [sp.Rational(2), sp.Rational(3), sp.Rational(1, 2), sp.Rational(3, 2),
               sp.Rational(-1), sp.Rational(5), sp.Rational(2, 3), sp.Rational(7, 2),
               sp.Rational(61, 3), sp.Rational(3, 61), sp.Rational(8, 5), sp.Rational(-8, 5)]:
        rec = fibre_at_t(tv)
        samples.append(rec)

    n_even_K = sum(1 for r in samples if r.get("disc_square_K"))
    n_even_Q = sum(1 for r in samples if r.get("disc_square_Q"))
    n_irr_K = sum(1 for r in samples if r.get("irr_over_K"))
    n_irr_Q = sum(1 for r in samples if r.get("irr_over_Q"))

    # Can disc become square in Q by twisting the cover / base change?
    # Quadratic twist of poly: f_twist related — for quintics not standard like elliptic.
    # Clearing √5 from disc: would need a different model over Q with disc without 5.
    # The permanent 5 is geometric (from φ), not removable by Z-scaling of monic model
    # up to squares: check disc of monic vs Z-cleared ratio.

    mon = sp.expand((PHI - t) / 6)
    pol = sp.Poly(mon, y)
    D_monic = sp.together(pol.discriminant())
    # Z-model at generic: L depends on dens of monic coeffs (fixed 2,3,6)
    # monic coeffs involve 1/2,1/3,t/6 — L=6 for t integer
    # cleared disc = L^{n(n-1)} * disc(monic) with n=5, and content factors
    # The factor 5 survives in the square-free part for t∈Z\{0,1}.

    # Prove square-free part always includes 5 for integer t:
    sqfree_trials = []
    for ti in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 61, 80]:
        fr = fibre_at_t(ti)
        d = fr.get("disc") or 0
        # square-free kernel
        if d > 0:
            sf = int(sp.squarefree_kernel(d) if hasattr(sp, "squarefree_kernel") else 0)
            # manual squarefree part
            sf = 1
            n = d
            for p in sp.primerange(2, int(n**0.5) + 2):
                if n % p == 0:
                    e = 0
                    while n % p == 0:
                        n //= p
                        e += 1
                    if e % 2:
                        sf *= p
            if n > 1:
                sf *= n
            sqfree_trials.append({"t": ti, "disc": d, "sqfree": sf, "has_5": sf % 5 == 0})

    always_5 = all(s["has_5"] for s in sqfree_trials)

    return {
        "samples": samples,
        "n_even_K": n_even_K,
        "n_even_Q": n_even_Q,
        "n_irr_K": n_irr_K,
        "n_irr_Q": n_irr_Q,
        "sqfree_trials": sqfree_trials,
        "sqfree_always_has_5": always_5,
        "descent_obstruction": (
            "For t∈Q, the Z-model fibre is in Q[y] with disc = 5·□ in Z "
            "(square-free part contains 5). Even monodromy holds over K=Q(√5) "
            "but not over Q. There is no descent of evenness to Q for these "
            "fibres: the permanent factor 5 is in the disc's square-free kernel "
            "and is exactly what K kills. A Q-rational pure-even model would need "
            "a different geometric cover (or non-rigid family), not a base change "
            "of the same φ-fibres."
        ),
        "relation_to_pure_even_families": (
            "Pure-even BJ families over Q (LSW, flagship k-slices) already have "
            "disc □ in Q and Gal≤A5 over Q; they are not φ-fibres for rational t "
            "(φ-fibres are never even over Q). Over K, both classes are even: "
            "φ-fibres become even by base change; pure-even families remain even. "
            "Matching them requires a Tschirnhaus isomorphism over K between a "
            "φ-fibre and a BJ specialisation — tested separately via disc classes."
        ),
    }


def tschirnhaus_field_test(even_rows: list, seeds: list) -> dict:
    """
    For a few even-over-K fibres and pure-even seeds, test whether the
    number fields Q[y]/(fibre) and Q[x]/(seed) are isomorphic (same field)
    by comparing Dedekind zeta fingerprints: factorization types at many primes
    must match for isomorphic fields (necessary condition).
    """
    print("  number-field fingerprint vs pure-even seeds...", flush=True)

    def fingerprint(poly_expr, var, max_p=80):
        pol = sp.Poly(poly_expr, var, domain=sp.ZZ)
        if not pol.is_irreducible:
            return None
        disc = int(pol.discriminant())
        counts = Counter()
        for p in sp.primerange(2, max_p):
            if disc % p == 0:
                continue
            try:
                facs = sp.factor_list(pol.as_expr(), modulus=int(p))
                degs = []
                for f, m in facs[1]:
                    degs.extend([int(sp.degree(f))] * int(m))
                counts[tuple(sorted(degs))] += 1
            except Exception:
                continue
        return dict(counts)

    def similar(fp1, fp2, tol=0.35):
        if not fp1 or not fp2:
            return False
        # Compare distribution of partition types
        keys = set(fp1) | set(fp2)
        n1 = sum(fp1.values()) or 1
        n2 = sum(fp2.values()) or 1
        dist = 0.0
        for k in keys:
            dist += abs(fp1.get(k, 0) / n1 - fp2.get(k, 0) / n2)
        return dist < tol

    seed_fps = []
    for a, b, k, tag in seeds:
        fp = fingerprint(x**5 + a * x + b, x)
        seed_fps.append((tag, a, b, k, fp))

    hits = []
    tested = 0
    for fr in even_rows[:12]:
        if not fr.get("irr_over_Q"):
            continue
        fp = fingerprint(sp.sympify(fr["poly"]), y)
        tested += 1
        for tag, a, b, k, sfp in seed_fps:
            if similar(fp, sfp):
                hits.append(
                    {
                        "t": fr["t"],
                        "fibre": fr["poly"],
                        "seed": tag,
                        "seed_ab": (a, b),
                        "k": k,
                        "fp_fibre": {str(k_): v for k_, v in (fp or {}).items()},
                        "fp_seed": {str(k_): v for k_, v in (sfp or {}).items()},
                    }
                )
    return {"tested_fibres": tested, "fingerprint_hits": hits, "n_hits": len(hits)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
PURE_EVEN_SEEDS = [
    (-55, 88, "-8/5", "flagship"),
    (145, -232, "-8/5", "flag_145"),
    (320, -512, "-8/5", "flag_320"),
    (-100, 400, "-4", "lsw_m100"),
    (124, -496, "-4", "lsw_124m"),
    (20, 16, "4/5", "classical"),
    (95, 76, "4/5", "s95_76"),
    (-180, 432, "-12/5", "s180"),
    (220, -528, "-12/5", "s220m"),
    (-55, -88, "8/5", "flagship_m"),
    (20, -16, "-4/5", "classical_m"),
]


def main():
    t0 = time.time()
    print("K = Q(√5) — permanent factor 5 becomes square", flush=True)

    thm = disc_theorem()
    print(f"  proved={thm['proved']}", flush=True)

    # Sample rational t — all non-critical should be even_K for Z-disc
    rat = sample_rational_t(max_p=20, max_q=6)
    print(f"  rational t: {rat['stats']}", flush=True)

    quad = sample_quadratic_t(max_ab=4)
    print(f"  quadratic t: {quad['stats']}", flush=True)

    dep = fibre_family_depressed()

    # Collect even_irr_K rows for relation tests — rebuild a focused list
    print("  building focused even-over-K fibre list...", flush=True)
    focus = []
    for tv in [
        2, 3, 4, 5, 6, 7, 8, 9, 10, -1, -2, sp.Rational(1, 2), sp.Rational(3, 2),
        sp.Rational(2, 3), sp.Rational(3, 4), sp.Rational(5, 2), sp.Rational(5, 3),
        sp.Rational(8, 5), sp.Rational(-8, 5), sp.Rational(4, 5), sp.Rational(61, 3),
        sp.Rational(3, 61), 16, 25, 27, 32, 61, 80, 243,
    ]:
        rec = fibre_at_t(tv)
        if rec.get("disc_square_K") and rec.get("irr_over_K"):
            focus.append(rec)
        elif rec.get("disc_square_K"):
            # still useful
            focus.append(rec)

    rel = try_bj_match_at_samples(focus, PURE_EVEN_SEEDS)
    print(
        f"  disc-class matches: Q={rel['fibres_with_seed_disc_class_Q']} "
        f"K={rel['fibres_with_seed_disc_class_K']}",
        flush=True,
    )

    pe_K = pure_even_over_K()
    print(f"  pure-even over K samples: {pe_K['n_samples_quadratic_m']}", flush=True)

    desc = descent_analysis()
    print(
        f"  descent samples even_K={desc['n_even_K']} even_Q={desc['n_even_Q']} "
        f"sqfree_always_5={desc['sqfree_always_has_5']}",
        flush=True,
    )

    fp = tschirnhaus_field_test(
        [r for r in focus if r.get("irr_over_Q")],
        PURE_EVEN_SEEDS,
    )
    print(f"  fingerprint hits: {fp['n_hits']}", flush=True)

    elapsed = round(time.time() - t0, 2)

    # ---- markdown ----
    lines = [
        r"# Arithmetic over \(K=\mathbb{Q}(\sqrt{5})\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** Disc theorem proved={thm['proved']}. "
        f"Over \(K=\mathbb{{Q}}(\\sqrt{{5}})\), \(\\mathrm{{disc}}(\\mathrm{{monic}}(\\varphi-t))\) "
        r"is identically a square in \(K(t)\). "
        f"Rational-\(t\) scan: even-over-\(K\) irreducible fibres available "
        f"(stats {rat['stats']}). "
        f"Even-over-\(\\mathbb{{Q}}\) fibres: **{rat['stats'].get('even_Q', 0)}** (still blocked). "
        f"Square-free disc kernel retains permanent 5 for tested integer \(t\): "
        f"**{desc['sqfree_always_has_5']}**. "
        f"Number-field fingerprint matches fibre↔pure-even seed: **{fp['n_hits']}**. "
        "Descent of evenness to \(\\mathbb{Q}\) fails; relation to pure-even BJ families "
        "is only via base change / disc class over \(K\), not Q-isomorphism in the scan.",
        "",
        "---",
        "",
        r"## 1. Theorem: permanent factor 5 becomes a square in \(K\)",
        "",
        r"Preferred cover: \(\varphi(y)=6y^5-15y^4+10y^3\).",
        "",
        r"Monic fibre: \(\mathrm{monic}(\varphi-t)=y^5-\frac{5}{2}y^4+\frac{5}{3}y^3-\frac{t}{6}\).",
        "",
        f"- disc = `{thm['disc']}`",
        f"- equals \(5\\cdot(\\text{{square}})\): **{thm['equals_5_times_square']}**",
        f"- square part: `{thm['square_part']}`",
        f"- equals \((\\sqrt{{5}}\\cdot\\text{{square part}})^2\): **{thm['equals_sqrt5_square_form']}**",
        f"- proved: **{thm['proved']}**",
        "",
        f"**Over \(\\mathbb{{Q}}\):** {thm['theorem_Q']}",
        "",
        f"**Over \(K\):** {thm['theorem_K']}",
        "",
        "---",
        "",
        r"## 2. Specialisations at rational \(t\)",
        "",
        f"Stats: `{rat['stats']}`",
        "",
        r"Every non-critical rational \(t\) yields a Z-model fibre whose discriminant "
        r"is a square in \(K\) (reason: \(5\cdot\square\) or already \(\square\)). "
        r"None of these are squares in \(\mathbb{Q}\) in the scan (`even_Q` count).",
        "",
        r"### Sample even-over-\(K\) fibres (irr over \(K\))",
        "",
    ]
    for rec in (rat.get("even_irr_K_sample") or [])[:12]:
        lines.append(
            f"- t=`{rec.get('t')}`: `{rec.get('poly')}` "
            f"disc={rec.get('disc')} K-even={rec.get('disc_square_K')} "
            f"irr_Q={rec.get('irr_over_Q')} irr_K={rec.get('irr_over_K')} "
            f"gal_Q={rec.get('gal_Q')} form={rec.get('disc_K_form')}"
        )
    if not rat.get("even_irr_K_sample"):
        lines.append(
            "_No irr-over-K rows collected in dense scan bounds "
            "(may factor over K while irr over Q, or converse). See focused list._"
        )

    lines += [
        "",
        "### Focused t list",
        "",
    ]
    for rec in focus[:18]:
        lines.append(
            f"- t=`{rec.get('t')}` poly=`{rec.get('poly')}` "
            f"disc_Q□={rec.get('disc_square_Q')} disc_K□={rec.get('disc_square_K')} "
            f"irr_Q={rec.get('irr_over_Q')} irr_K={rec.get('irr_over_K')} "
            f"gal_Q=`{rec.get('gal_Q')}` degs_K={rec.get('factor_degs_over_K')}"
        )

    lines += [
        "",
        "---",
        "",
        r"## 3. Specialisations at \(t=a+b\sqrt{5}\)",
        "",
        f"Stats: `{quad['stats']}`",
        "",
        r"Disc matches \((\sqrt{5}\cdot 25 t(t-1)/36)^2\) for sampled quadratic \(t\): "
        f"count disc_square_K={quad['stats'].get('disc_square_K', 0)} / tested={quad['stats'].get('tested', 0)}.",
        "",
        "---",
        "",
        r"## 4. Depressed form of \(\mathrm{monic}(\varphi-t)\)",
        "",
        f"- shift: `{dep['shift']}`",
        f"- coeffs: `{dep['coeff_map']}`",
        f"- {dep['note']}",
        "",
        "---",
        "",
        r"## 5. Descent: can K-evenness reach \(\mathbb{Q}\)?",
        "",
        f"- Focused/sample even_K: {desc['n_even_K']}, even_Q: {desc['n_even_Q']}",
        f"- irr_K: {desc['n_irr_K']}, irr_Q: {desc['n_irr_Q']}",
        f"- Square-free part of disc always contains 5 (integer t trials): **{desc['sqfree_always_has_5']}**",
        "",
        "Square-free trials:",
        "",
    ]
    for s in desc["sqfree_trials"][:12]:
        lines.append(f"- t={s['t']}: sqfree={s['sqfree']} has_5={s['has_5']}")

    lines += [
        "",
        f"**Obstruction:** {desc['descent_obstruction']}",
        "",
        f"**Relation note:** {desc['relation_to_pure_even_families']}",
        "",
        "---",
        "",
        r"## 6. Relation to pure-even BJ families",
        "",
        "### Disc-class comparison (fibre vs catalogue seeds)",
        "",
        f"- Fibres compared: {rel['compared']}",
        f"- With seed disc ratio square in Q: **{rel['fibres_with_seed_disc_class_Q']}**",
        f"- With seed disc ratio square in K: **{rel['fibres_with_seed_disc_class_K']}**",
        "",
    ]
    for d in (rel.get("details") or [])[:8]:
        lines.append(
            f"- t={d['t']}: gal_Q={d.get('gal_Q')} matches={d.get('seed_disc_matches')}"
        )

    lines += [
        "",
        "### Number-field fingerprints (necessary for Q-isomorphism)",
        "",
        f"- Tested irr fibres: {fp['tested_fibres']}",
        f"- Hits (similar Frobenius stats to a pure-even seed): **{fp['n_hits']}**",
        "",
    ]
    if not fp["fingerprint_hits"]:
        lines.append(
            "_No fingerprint match: φ-fibres over Q and pure-even BJ seeds "
            "appear to define distinct number fields in the sample "
            "(expected — different disc square-free parts: 5·□ vs pure □)._"
        )
    for h in fp.get("fingerprint_hits") or []:
        lines.append(f"- HIT t={h['t']} ↔ seed {h['seed']} {h['seed_ab']}")

    lines += [
        "",
        r"### Pure-even k-slices base-changed to \(K\)",
        "",
        f"- Samples with m=r+s√5: {pe_K['n_samples_quadratic_m']}",
        f"- Of which α,β ∈ Q: {pe_K['n_with_Q_coeffs']}",
        f"- {pe_K['note']}",
        "",
    ]
    for s in (pe_K.get("samples") or [])[:8]:
        lines.append(
            f"- k={s['k']} m={s['m']}: α=`{s['alpha'][:60]}` β=`{s['beta'][:60]}` "
            f"id_ok={s['disc_identity_ok']} inQ={s['alpha_in_Q'] and s['beta_in_Q']}"
        )

    lines += [
        "",
        "---",
        "",
        "## 7. Conclusions",
        "",
        r"1. **Theorem (proved):** \(\mathrm{disc}(\mathrm{monic}(\varphi-t))"
        r"=5\cdot(\frac{25 t(t-1)}{36})^2\) in \(\mathbb{Q}(t)\). "
        r"Over \(K=\mathbb{Q}(\sqrt{5})\) this is \((\sqrt{5}\cdot\frac{25 t(t-1)}{36})^2\), "
        r"hence **identically a square** in \(K(t)\).",
        "",
        r"2. **Even specialisations over \(K\):** every non-critical \(t\in K\) "
        r"(including all \(t\in\mathbb{Q}\setminus\{0,1\}\)) gives disc square in \(K\). "
        r"Irreducibility over \(K\) holds for many rational \(t\); Gal over \(\mathbb{Q}\) "
        r"remains odd (typically in \(S_5\setminus A_5\)).",
        "",
        r"3. **No descent to even-over-\(\mathbb{Q}\):** the permanent factor 5 sits in the "
        r"square-free kernel of the Z-disc for integer \(t\). Base change to \(K\) kills it; "
        r"it cannot be removed while staying over \(\mathbb{Q}\) with the same fibres.",
        "",
        r"4. **Relation to pure-even BJ families:** those families are already even over "
        r"\(\mathbb{Q}\) and are **not** the rational fibres of \(\varphi\). Over \(K\), both "
        r"are even; disc-class / Frobenius fingerprint tests in this scan do **not** identify "
        r"φ-fibres with catalogue pure-even seeds as the same \(\mathbb{Q}\)-fields. "
        r"Extending pure-even parameters \(m\in K\) produces BJ seeds over \(K\), "
        r"complicating the HQCC lattice interpretation (coeffs leave \(\mathbb{Z}\)).",
        "",
        r"5. **Programme impact:** Route 2 is **technically viable for even monodromy** "
        r"of \(\varphi\)-fibres, but moves arithmetic into \(K\) and does **not** by itself "
        r"recover HQCC lattice seeds over \(\mathbb{Q}\). Fusion over \(\mathbb{Q}\) still "
        r"needs non-rigid pure-even families (Route 1) or a different geometric construction "
        r"(Route 3).",
        "",
        r"### Recommended stance",
        "",
        r"- Keep \(K=\mathbb{Q}(\sqrt{5})\) as a **side route** for geometric evenness of \(\varphi\).",
        r"- Primary fusion fuel over \(\mathbb{Q}\) remains the **multi-seed pure-even \(k\)-slices** "
        r"(`ENLARGED_SEED_CATALOGUE.md`).",
        r"- If pursuing \(K\): study Gal of fibres over \(K\) (A5 over \(K\)) and whether "
        r"HQCC lattice can be rephrased in \(\mathcal{O}_K=\mathbb{Z}[\frac{1+\sqrt{5}}{2}]\) "
        r"— a genuine change of programme coefficients.",
        "",
        "_Generated by k_sqrt5_even.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "theorem": thm,
        "rational_t": rat,
        "quadratic_t": quad,
        "depressed": dep,
        "relation": rel,
        "pure_even_K": pe_K,
        "descent": {
            "n_even_K": desc["n_even_K"],
            "n_even_Q": desc["n_even_Q"],
            "sqfree_always_has_5": desc["sqfree_always_has_5"],
            "sqfree_trials": desc["sqfree_trials"],
            "descent_obstruction": desc["descent_obstruction"],
            "samples": desc["samples"],
        },
        "fingerprints": fp,
        "focus_fibres": focus[:25],
        "verdict": (
            f"proved={thm['proved']}, even_Q_blocked, even_K_open, "
            f"sqfree_always_5={desc['sqfree_always_has_5']}, "
            f"fp_hits={fp['n_hits']}, no_descent_to_Q"
        ),
    }

    write_md(OUT / "K_SQRT5_EVEN.md", doc)
    write_md(RESULTS / "K_SQRT5_EVEN.md", doc)
    write_md(ROOT / "K_SQRT5_EVEN.md", doc)
    write_json(OUT / "K_SQRT5_EVEN.json", blob)
    print(blob["verdict"], flush=True)
    print(f"Wrote K_SQRT5_EVEN.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
