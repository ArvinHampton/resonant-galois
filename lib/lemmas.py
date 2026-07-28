"""
Closed-form lemmas and thin-class utilities for theorem promotion.

These are the arithmetic facts the Criteria 1–3 attacks rest on.
"""
from __future__ import annotations

import sympy as sp

from .common import is_square, monic_poly, classify_poly, x

# ---------------------------------------------------------------------------
# Lemma: discriminant of Bring–Jerrard quintic x^5 + a x + b
# ---------------------------------------------------------------------------
def disc_bj(a, b) -> sp.Expr:
    """
    disc(x^5 + a x + b) = 256 a^5 + 3125 b^4.

    Verified symbolically in verify_disc_formulas().
    """
    return 256 * a**5 + 3125 * b**4


def disc_bj_int(a: int, b: int) -> int:
    return int(256 * a**5 + 3125 * b**4)


# ---------------------------------------------------------------------------
# Lemma: depressed / icosahedral-adjacent family
#   f = x^5 + p x^3 + q x + r
# ---------------------------------------------------------------------------
def disc_pqr(p, q, r) -> sp.Expr:
    """
    Discriminant of x^5 + p x^3 + q x + r (exact sympy formula, simplified).
    Cached computation via sympy Poly.
    """
    f = sp.Poly(x**5 + p * x**3 + q * x + r, x)
    return sp.factor(sp.simplify(f.discriminant()))


def disc_icosa(m, n) -> sp.Expr:
    """disc of x^5 + 5 m x^3 + 5 m^2 x + n."""
    return disc_pqr(5 * m, 5 * m**2, n)


# ---------------------------------------------------------------------------
# Lemma (group theory, stated): operational A5 criterion
# ---------------------------------------------------------------------------
OPERATIONAL_A5 = (
    "For irreducible monic f in Z[x] of degree 5: "
    "disc(f) square and some unramified p with factorization type (3,1,1) "
    "implies Gal(f/Q) = A5 "
    "(only transitive subgroup of A5 containing a 3-cycle)."
)


def verify_disc_formulas(trials: int = 30) -> dict:
    """Check closed forms against sympy Poly.discriminant on random samples."""
    import random

    rng = random.Random(539)
    bj_ok = 0
    bj_fail = []
    for _ in range(trials):
        a = rng.randint(-20, 20)
        b = rng.randint(-20, 20)
        if b == 0:
            continue
        pol = sp.Poly(x**5 + a * x + b, x, domain=sp.ZZ)
        d_sym = int(pol.discriminant())
        d_form = disc_bj_int(a, b)
        if d_sym == d_form:
            bj_ok += 1
        else:
            bj_fail.append((a, b, d_sym, d_form))

    # Symbolic identity for BJ
    a, b = sp.symbols("a b")
    f = sp.Poly(x**5 + a * x + b, x)
    identity_bj = sp.expand(f.discriminant() - disc_bj(a, b)) == 0

    # Icosa family identity check numerically
    ico_ok = 0
    ico_fail = []
    m, n = sp.symbols("m n")
    fico = sp.Poly(x**5 + 5 * m * x**3 + 5 * m**2 * x + n, x)
    form_ico = sp.factor(fico.discriminant())
    for _ in range(trials):
        mv = rng.randint(-8, 8)
        nv = rng.randint(-20, 20)
        pol = sp.Poly(x**5 + 5 * mv * x**3 + 5 * mv**2 * x + nv, x, domain=sp.ZZ)
        d_sym = int(pol.discriminant())
        d_form = int(form_ico.subs({m: mv, n: nv}))
        if d_sym == d_form:
            ico_ok += 1
        else:
            ico_fail.append((mv, nv, d_sym, d_form))

    return {
        "bj_numeric_ok": bj_ok,
        "bj_numeric_fail": bj_fail[:5],
        "bj_symbolic_identity": bool(identity_bj),
        "bj_formula": "256*a**5 + 3125*b**4",
        "icosa_formula": str(form_ico),
        "icosa_numeric_ok": ico_ok,
        "icosa_numeric_fail": ico_fail[:5],
        "operational_A5": OPERATIONAL_A5,
    }


def bj_evenness_condition(a: int, b: int) -> dict:
    """
    Thin-class decision for BJ: Gal ≤ A5 (among transitive) requires disc square.
    Returns disc, whether square, and classification if irreducible.
    """
    disc = disc_bj_int(a, b)
    rec = {
        "a": a,
        "b": b,
        "poly": f"x**5 + ({a})*x + ({b})",
        "disc": disc,
        "disc_square": is_square(disc) if disc > 0 else False,
    }
    if rec["disc_square"]:
        full = classify_poly(x**5 + a * x + b, do_galois=True)
        rec.update({k: full.get(k) for k in ("irreducible", "galois", "status", "census")})
    else:
        pol = monic_poly(x**5 + a * x + b)
        if pol is not None:
            rec["irreducible"] = bool(pol.is_irreducible)
    return rec


def search_bj_square_disc(values: list[int], max_hits: int = 40) -> list[dict]:
    hits = []
    for a in values:
        for b in values:
            if b == 0:
                continue
            rec = bj_evenness_condition(a, b)
            if rec.get("disc_square") and rec.get("irreducible"):
                hits.append(rec)
                if len(hits) >= max_hits:
                    return hits
    return hits


def disc_homogenised_A5(tval) -> sp.Expr:
    """
    disc(x^5 + 20 t^4 x + 16 t^5) = t^{20} * disc(x^5 + 20 x + 16).

    Proof: BJ formula with a=20 t^4, b=16 t^5:
      256 (20 t^4)^5 + 3125 (16 t^5)^4
      = 256 * 20^5 t^{20} + 3125 * 16^4 t^{20}
      = t^{20} (256*20^5 + 3125*16^4)
      = t^{20} * disc(x^5+20x+16).

    Since disc(x^5+20x+16) is a positive square (known A5 seed),
    disc(f_t) is a square for all t ≠ 0 (and 0 for t=0, where f is not deg 5 monic irr).
    """
    t = sp.symbols("t")
    seed = disc_bj(20, 16)
    return sp.factor(sp.expand(disc_bj(20 * t**4, 16 * t**5)))


def prove_homogenised_A5_even() -> dict:
    """Symbolic proof that homogenised A5 family is always even for t ≠ 0."""
    t = sp.symbols("t")
    formula = disc_bj(20 * t**4, 16 * t**5)
    seed = disc_bj_int(20, 16)
    factored = sp.factor(sp.expand(formula))
    # factored should be t**20 * seed
    quotient = sp.simplify(sp.expand(formula) / (t**20))
    seed_square = is_square(seed)
    # For integer t≠0, t**20 is always a square ((t**10)^2), seed square ⇒ product square
    return {
        "seed_poly": "x**5 + 20*x + 16",
        "seed_disc": seed,
        "seed_disc_is_square": seed_square,
        "family": "x**5 + 20*t**4*x + 16*t**5",
        "disc_factored": str(factored),
        "disc_over_t20": str(quotient),
        "disc_over_t20_equals_seed": int(quotient) == seed if quotient.is_number else sp.expand(quotient - seed) == 0,
        "theorem": (
            "For all t in Z\\{0}, disc(x^5+20 t^4 x+16 t^5) "
            "= t^{20} * disc(x^5+20x+16) is a square in Z, "
            "because t^{20}=(t^{10})^2 and disc(seed) is a square. "
            "Hence Gal(f_t/Q) ≤ A5 whenever f_t is irreducible; "
            "with a (3,1,1) Frobenius, Gal = A5."
        ),
        "proved": seed_square and (sp.expand(quotient - seed) == 0),
    }
