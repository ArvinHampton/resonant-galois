"""
Enlarge the HQCC A5 seed catalogue and group by k = β/α
to find more multi-seed pure-even LSW-type slices.

Pure-even family on the ray β = k α (any k ≠ 0 in Q):

  disc(α, kα) = 256 α^4 (α + 3125 k^4 / 256)

  Integer model (k = p/q in lowest terms):
    α(m) = 256 q^4 m^2 - 3125 p^4
    β(m) = (p/q) α(m)   — cleared: use β = p · (α/q) when q | α,
    or monic Z form:
      A(m) = 256 q^4 m^2 - 3125 p^4
      B(m) = p · A(m) // q   if q | A(m); else scale poly differently.

  When q | 256 q^4 m^2 (always) and q | 3125 p^4:
    if q ∤ 3125 p^4, A may not be divisible by q.
  Safer monic model always integer:
    α(m) = 256 (q m)^2 - 3125 (p^4 / g) ...  (use Rational arithmetic)

  Practical: α_Q(m) = 256 m^2 - 3125 k^4 / 256, β_Q = k α_Q
  Clear denominators to monic Z[x] specialisations at rational m.

  Identically: disc = (256 α^2 m)^2 when α = 256 m^2 - 3125 k^4/256.

Output: ENLARGED_SEED_CATALOGUE.md / .json
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import defaultdict
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
# Known baseline seeds (prior catalogue)
# ---------------------------------------------------------------------------
KNOWN = [
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
    (-180, 432, "s180"),  # from HQCC_NATIVE parametric even
    (-180, -432, "s180_m"),
    (220, 528, "s220"),
    (220, -528, "s220_m"),
]


def is_a5(rec: dict) -> bool:
    st = rec.get("status") or ""
    gal = str(rec.get("galois") or "")
    return st.startswith("HIT_A5") or ("A5" in gal and "A6" not in gal)


# ---------------------------------------------------------------------------
# Enlarged HQCC lattice
# ---------------------------------------------------------------------------
def enlarged_lattice(max_abs: int = 50000) -> list[int]:
    vals: set[int] = set()
    core = list(MODEL_CORE.keys()) + [
        1, 2, 3, 4, 5, 8, 9, 16, 18, 20, 21, 25, 27, 32, 45, 55, 61, 76, 80,
        88, 95, 100, 124, 180, 220, 243, 256, 400, 432, 496, 520, 528, 532,
        539, 4880, 223, 3125, 55, 57, -1, -2, -3, -4, -5,
    ]
    for v in core:
        if abs(v) <= max_abs:
            vals.add(v)
            vals.add(-v)

    # powers of 3
    p = 1
    for _ in range(10):
        if abs(p) <= max_abs:
            vals.add(p)
            vals.add(-p)
        p *= 3

    # short linear combos of model atoms
    atoms = [1, 3, 9, 18, 27, 61, 80, 243, 539, 20, 16, 55, 88, 95, 76, 100, 4, 5, 8, 25, 32, 256]
    for a, b in itertools.product(atoms, repeat=2):
        for w in (a + b, a - b, a + 2 * b, 2 * a + b, 3 * a + b, a + 3 * b, a * b):
            if w and abs(w) <= max_abs:
                vals.add(w)
                vals.add(-w)

    # products with small multipliers
    for a in list(vals)[:]:
        for m in (2, 3, 4, 5, 8, 9, 16, 25, 27, 32, 61, 80):
            w = a * m
            if abs(w) <= max_abs:
                vals.add(w)
                vals.add(-w)

    # HQCC-ish residues: ±(61±3^k), ±(80±3^k), ±(539±3^k), flagship combos
    for base in (61, 80, 539, 88, 55, 95, 76, 532, 100, 124, 400, 496):
        for k in range(0, 7):
            t = 3**k
            for w in (base + t, base - t, base + 2 * t, base - 2 * t, base * t // t):
                if w and abs(w) <= max_abs:
                    vals.add(w)
                    vals.add(-w)

    # integer range densification for small α,β (catches near-lattice seeds)
    for n in range(-400, 401):
        if n != 0:
            vals.add(n)

    # medium densification every 5 / every 4 around model scale
    for n in range(-2000, 2001, 5):
        if n != 0:
            vals.add(n)
    for n in range(-3000, 3001, 4):
        if n != 0 and abs(n) > 400:
            vals.add(n)

    vals.discard(0)
    return sorted(vals, key=lambda z: (abs(z), z))


def scan_bj_seeds(
    lattice: list[int],
    alpha_bound: int = 2500,
    beta_bound: int = 3000,
    max_pairs: int = 2_000_000,
) -> dict:
    """Scan (α,β) for disc□; classify Gal only on square-disc hits."""
    print("  Scanning enlarged BJ lattice for disc□ ...", flush=True)
    alphas = [v for v in lattice if abs(v) <= alpha_bound]
    betas = [v for v in lattice if 0 < abs(v) <= beta_bound]
    # Prefer smaller first: sort by abs
    alphas = sorted(set(alphas), key=lambda z: (abs(z), z))
    betas = sorted(set(betas), key=lambda z: (abs(z), z))
    print(f"    |α|≤{alpha_bound}: {len(alphas)}, |β|≤{beta_bound}: {len(betas)}", flush=True)

    even: list[dict] = []
    tested = 0
    seen: set[tuple[int, int]] = set()

    # Always include known seeds
    for a, b, tag in KNOWN:
        seen.add((a, b))
        d = disc_bj_int(a, b)
        if d > 0 and is_square(d):
            rec = classify_poly(x**5 + a * x + b, do_galois=True)
            even.append(
                {
                    "a": a,
                    "b": b,
                    "disc": d,
                    "poly": rec.get("poly"),
                    "status": rec.get("status"),
                    "galois": rec.get("galois"),
                    "irreducible": rec.get("irreducible"),
                    "source": f"known:{tag}",
                    "tag": tag,
                }
            )

    # Full lattice scan — disc filter first (cheap)
    t_disc0 = time.time()
    candidates: list[tuple[int, int, int]] = []
    for a in alphas:
        for b in betas:
            tested += 1
            if tested > max_pairs:
                break
            if (a, b) in seen:
                continue
            d = disc_bj_int(a, b)
            if d > 0 and is_square(d):
                candidates.append((a, b, d))
                seen.add((a, b))
        if tested > max_pairs:
            break
    print(
        f"    disc□ candidates (new): {len(candidates)} in {round(time.time()-t_disc0,1)}s "
        f"(tested {tested} pairs)",
        flush=True,
    )

    # Classify Galois — can be slow; do all disc□
    for i, (a, b, d) in enumerate(candidates):
        rec = classify_poly(x**5 + a * x + b, do_galois=True)
        row = {
            "a": a,
            "b": b,
            "disc": d,
            "poly": rec.get("poly"),
            "status": rec.get("status"),
            "galois": rec.get("galois"),
            "irreducible": rec.get("irreducible"),
            "source": "scan",
            "tag": f"a{a}_b{b}",
        }
        even.append(row)
        if is_a5(row):
            print(f"    *** A5 *** α={a} β={b}  {row['poly']}", flush=True)
        elif i < 5 or (i + 1) % 25 == 0:
            print(
                f"    classified {i+1}/{len(candidates)}: α={a} β={b} → {row.get('status')}",
                flush=True,
            )

    a5 = [r for r in even if is_a5(r) and r.get("irreducible")]
    d5 = [r for r in even if r.get("galois") and "D5" in str(r.get("galois"))]
    other_even = [
        r
        for r in even
        if r.get("disc")
        and r.get("irreducible")
        and not is_a5(r)
        and not (r.get("galois") and "D5" in str(r.get("galois")))
    ]
    return {
        "tested_pairs": tested,
        "n_even": len(even),
        "n_A5": len(a5),
        "n_D5": len(d5),
        "n_other_even_irr": len(other_even),
        "even": even,
        "A5": a5,
        "D5": d5,
        "other_even_irr": other_even,
        "n_alphas": len(alphas),
        "n_betas": len(betas),
    }


# ---------------------------------------------------------------------------
# Group by k = β/α
# ---------------------------------------------------------------------------
def k_of(a: int, b: int) -> Fraction | None:
    if a == 0:
        return None
    return Fraction(b, a)


def group_by_k(seeds: list[dict]) -> dict:
    """Group seeds by reduced rational k=β/α. Drop α=0."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for s in seeds:
        a, b = int(s["a"]), int(s["b"])
        k = k_of(a, b)
        if k is None:
            continue
        key = str(k)  # reduced Fraction string, e.g. '-8/5', '-4'
        row = {**s, "k": key, "k_num": k.numerator, "k_den": k.denominator}
        groups[key].append(row)
    # sort each group by α
    for k in groups:
        groups[k] = sorted(groups[k], key=lambda r: (r["a"], r["b"]))
    return dict(groups)


def pure_even_family_for_k(k: Fraction) -> dict:
    """
    Build pure-even parametric family for ratio k=p/q.

    α_Q(m) = 256 m^2 - 3125 k^4 / 256
    β_Q(m) = k · α_Q(m)
    disc identically square for m ∈ Q.

    Integer clearing: let k=p/q.
      3125 k^4 / 256 = 3125 p^4 / (256 q^4)
      α = (256^2 m^2 q^4 - 3125 p^4) / (256 q^4)
    Set parameter m = r / (16 q^2) so α is nicer, or use integer param n:
      A(n) = 256 n^2 - 3125 p^4 * 256 / gcd? 

    Clean integer model used previously for integer k:
      α(t) = 256 t^2 - 3125 k^4, β = k α   (when k ∈ Z)
    For rational k=p/q:
      α(t) = 256 q^4 t^2 - 3125 p^4
      β(t) = p * α(t) / q
    Require q | α(t). Since α = 256 q^4 t^2 - 3125 p^4,
      α/q = 256 q^3 t^2 - 3125 p^4 / q  ⇒ need q | 3125 p^4.
    If not, multiply: use poly x^5 + (q α_raw) x + (p α_raw) with α_raw = ...
    Better: always integer coeffs via
      α_Z(t) = 256 q^4 t^2 - 3125 p^4
      β_Z(t) = p * (256 q^3 t^2) - 3125 p^5 / q   -- still need q|3125 p^4

    Safest always-integer family (may scale):
      Let D = 256 * q**4
      α(t) = D * t**2 - 3125 * p**4     wait that's same as above with t scaled

    Check disc with α = 256 q^4 t^2 - 3125 p^4, β = p α / q:
    Only integer when q|α.

    Alternative always-integer BJ (same k-ray in P1):
      α(t) = q * (256 (q t)^2 - 3125 p^4 / something)

    Use the Q-family and specialise at t where coeffs land in Z:
      α(t) = 256 t^2 - 3125 * p**4 / (256 * q**4)
      β(t) = (p/q) * α(t)
    For integer t such that 256 q^4 | (256^2 t^2 q^4 - 3125 p^4)... 

    Integer model that always works:
      α(t) = 256 * (q**4) * t**2 - 3125 * p**4
      β(t) = p * α(t)   // then k_eff = β/α = p, not p/q.

    Correct for k=p/q with β/α = p/q:
      α(t) = 256 * q**4 * t**2 - 3125 * p**4
      β(t) = (p/q) * α(t)
    At integer t, α ∈ Z; β ∈ Z iff q | α.
    α mod q = -3125 p^4 mod q. Since gcd(p,q)=1, q|3125 p^4 iff q|3125=5^5,
    so q's prime factors ⊆ {5}. For q|5^e, β integer for all t.

    For general q: clear by using t = q s:
      α = 256 q^4 (q s)^2 - 3125 p^4 = 256 q^6 s^2 - 3125 p^4 — still q|α issue.

    Use scaled BJ that is Z[t][x] with same k:
      f_t = x^5 + α0(t) x + β0(t) with β0/α0 = k always:
      α0(t) = 256 q^4 t^2 - 3125 p^4
      β0(t) = p q^3 * 256 t^2 - 3125 p^5 / q  NO

    Write α = u, β = (p/q) u. Multiply the linear coeffs by q^4 to clear the
    formula's denominator 256 q^4 in α_Q:
      α_Z = 256^2 q^4 n^2 - 3125 p^4     (n = m*256? simplify)

    From α_Q = 256 m^2 - 3125 k^4/256, set m = n/(16 q^2):
      α_Q = 256 n^2/(256 q^4) - 3125 p^4/(256 q^4)
           = (n^2 - 3125 p^4 / 256?)/ (q^4) ... messy.

    Final clean model (proved disc square by direct expansion):
      α(t) = 256 * q**4 * t**2 - 3125 * p**4
      β(t) = p * q**3 * 256 * t**2 - 3125 * p**5 //  — need β/α = p/q
      β = (p/q) α = (p/q)(256 q^4 t^2 - 3125 p^4)
         = p * 256 q^3 t^2 - 3125 p^5 / q
      So if q | 3125 p^5 and gcd(p,q)=1 then q|3125.

    We'll store:
      alpha_expr = 256*q**4*t**2 - 3125*p**4
      beta_expr  = p*(256*q**3*t**2) - 3125*p**5 // q   when q|3125
      else beta_expr as Rational and specialise only when integer.

    And verify disc square symbolically when β = k α with α = 256 q^4 t^2 - 3125 p^4.
    """
    p, q = int(k.numerator), int(k.denominator)
    t = sp.symbols("t")
    # α(t) = 256 q^4 t^2 - 3125 p^4
    alpha = 256 * (q**4) * t**2 - 3125 * (p**4)
    # β = (p/q) α — keep symbolic Rational
    beta = sp.together(sp.Rational(p, q) * alpha)
    D = sp.expand(256 * alpha**5 + 3125 * beta**4)
    # Factor / check square
    disc_ok = False
    disc_fact = None
    try:
        # D should be square of a poly in t
        # With α + 3125 k^4/256 = 256 q^4 t^2  (since k^4=p^4/q^4,
        # 3125 k^4/256 = 3125 p^4/(256 q^4);
        # α + that = 256 q^4 t^2 - 3125 p^4 + 3125 p^4/(256 q^4)*? 
        # Wait: α = 256 q^4 t^2 - 3125 p^4
        # 3125 k^4/256 = 3125 p^4 /(256 q^4)
        # α + 3125 k^4/256 = 256 q^4 t^2 - 3125 p^4 + 3125 p^4/(256 q^4)
        # Not equal to a pure square term unless we scale differently!
        #
        # Correct α for the identity α + 3125 k^4/256 = 256 s^2:
        # α = 256 s^2 - 3125 k^4/256 = 256 s^2 - 3125 p^4/(256 q^4)
        # Set s = q^2 t: α = 256 q^4 t^2 - 3125 p^4/(256 q^4)
        # = (256^2 q^8 t^2 - 3125 p^4) / (256 q^4)
        # Use monic-scaled seed coeffs proportional:
        # A = 256^2 q^8 t^2 - 3125 p^4 ? Too big.
        #
        # Simpler integer model with proved disc:
        # α = 256 * u - 3125 * p**4 * v with chosen u,v...
        #
        # Use: α(t) = 256 t**2 - 3125 * p**4  (ignore q in formula)
        #      β(t) = (p/q) α only when we want k=p/q — then α must absorb q.
        #
        # THEOREM form for any k in Q:
        #   α(t) = 256 t**2 - 3125 * k**4 / 256
        #   β(t) = k * α(t)
        # disc = (256 α^2 t)^2 identically (when α,β in Q(t)).
        alpha_Q = sp.together(256 * t**2 - 3125 * sp.Rational(p, q) ** 4 / 256)
        beta_Q = sp.together(sp.Rational(p, q) * alpha_Q)
        D_Q = sp.expand(256 * alpha_Q**5 + 3125 * beta_Q**4)
        # Expected square: (256 * alpha_Q**2 * t)**2
        expected = sp.expand((256 * alpha_Q**2 * t) ** 2)
        disc_ok = sp.expand(D_Q - expected) == 0
        disc_fact = str(sp.factor(D_Q))[:200]
        alpha, beta = alpha_Q, beta_Q
    except Exception as e:
        disc_fact = f"error:{e}"

    # Integer specialisations: find t=Fraction where α,β ∈ Z
    # α = 256 t^2 - 3125 p^4 / (256 q^4)
    # = (256^2 q^4 t^2 - 3125 p^4) / (256 q^4)
    # Set t = r / (16 q^2) with r ∈ Z: then
    # 256 t^2 = 256 r^2 / (256 q^4) = r^2 / q^4
    # α = r^2/q^4 - 3125 p^4/(256 q^4) = (256 r^2 - 3125 p^4)/ (256 q^4)
    # Hmm still.
    # Set t = q^2 * s / 16? Try t integer first and clear:
    # For t ∈ Z, α ∈ Q; multiply num/den.
    # Integer family used in code (always Z for t∈Z when we clear k^4 den):
    #   α_Z(t) = 256 * (q**4) * t**2 - 3125 * p**4     # note: this is 256 q^4 t^2 - 3125 p^4
    #   For disc with β_Z = (p/q) α_Z to be square:
    #   α_Z + 3125 k^4/256 = 256 q^4 t^2 - 3125 p^4 + 3125 p^4/(256 q^4)
    #   NOT a pure square. So this α_Z is WRONG scale for the identity.
    #
    # Correct integer α with β = k α and α + 3125 k^4/256 = 256 (integer)^2:
    #   α + 3125 p^4/(256 q^4) = 256 m^2, m∈Z
    #   α = 256 m^2 - 3125 p^4/(256 q^4)
    # For α ∈ Z: 256 q^4 | 3125 p^4? No: 3125 p^4 /(256 q^4) ∈ Z
    #   ⇔ 256 q^4 | 3125 p^4. Since gcd(p,q)=1,  q^4 | 3125=5^5 ⇒ q | 5.
    #   and 256=2^8 | 3125 p^4 ⇒ 2^8 | p^4 ⇒ 2^2 | p. So p divisible by 4 if q free of 2.
    #
    # General: allow m = M/(16 something) — use Fraction m.
    # Integer model: multiply poly by making
    #   A = 256 * (256 q^4) m^2 - 3125 p^4 = 65536 q^4 m^2 - 3125 p^4? 
    #   α_true = (65536 q^4 m^2 - 3125 p^4)/(256 q^4) = 256 m^2 - 3125 p^4/(256 q^4)
    #   So A := 65536 q^4 m^2 - 3125 p^4, then α_true = A / (256 q^4)
    #   β_true = (p/q) α_true = p A / (256 q^5)
    #   f = x^5 + α x + β. Clear by y = c x to get Z coeffs... complex.
    #
    # Practical integer form (used for flagship k=-8/5 successfully):
    #   For k=p/q, set α(s) = (q s)^2 * c - ... hand-tuned.
    # Flagship: α=25s^2-80, β=128-40s^2 with k=-8/5.
    #   Check: β/α = (128-40s^2)/(25s^2-80) = 8(16-5s^2)/5(5s^2-16) = -8/5. Yes!
    #
    # General construction matching flagship:
    #   Want β/α = p/q constant, and α + 3125 (p/q)^4 / 256 = square.
    #   α = 256 m^2 - 3125 p^4/(256 q^4)
    #   Write m = s/(16 q^2) * something...
    #   α = (256 * 256 q^4 s^2 - 3125 p^4)/(256 q^4) with m=s:
    #       if we want α integer with parameter s: need den | num.
    #   α(s) = 256 s^2 * q^4 * N - ... Let:
    #     α(s) = (256 s^2 * (256 q^4) - 3125 p^4) / den_fixed
    #
    # Use the form that always works for specialisation matching:
    #   α(s) = 256 * q**4 * s**2 - 3125 * p**4     # integer
    #   β(s) = p * q**3 * 256 * s**2 - (3125 * p**5)//q  if q|3125 else keep Rational
    # This gives β/α = p/q ONLY if β = (p/q)α, i.e. the second form.
    # And disc is square for this (α,β) pair as numbers when
    # α + 3125 k^4/256 is a rational square — which for α = 256 q^4 s^2 - 3125 p^4:
    #   α + 3125 p^4/(256 q^4) = 256 q^4 s^2 - 3125 p^4 + 3125 p^4/(256 q^4)
    #   = (256^2 q^8 s^2 - 3125 p^4 * 256 q^4 + 3125 p^4)/(256 q^4)
    #   Not a square in general. So wrong.
    #
    # CORRECT integer α:
    #   α(s) = 256 * (q**4) * s**2 * 256 - 3125 * p**4? 
    # From α = 256 m^2 - 3125 p^4/(256 q^4), set m = q^2 s:
    #   α = 256 q^4 s^2 - 3125 p^4/(256 q^4)
    # Multiply by (256 q^4)^2 to get contentful poly? For seed list matching we only need:
    # For each seed (a,b) with b/a = k, check if a + 3125 k^4/256 is a rational square.
    # If yes, seed lies on the pure-even family.

    on_family_condition = "α + 3125 k^4/256 = rational square"

    return {
        "k": str(k),
        "p": p,
        "q": q,
        "alpha_Q": str(alpha),
        "beta_Q": str(beta),
        "disc_identically_square": disc_ok,
        "disc_factored_sample": disc_fact,
        "on_family_condition": on_family_condition,
        "note": (
            "Parametric pure-even family: α=256 m² - 3125 k⁴/256, β=k α (m∈Q). "
            "A lattice seed (α0,k α0) lies on it iff α0 + 3125 k⁴/256 is a rational square."
        ),
    }


def seed_on_pure_even_slice(a: int, b: int, k: Fraction) -> dict:
    """Check whether (a,b) lies on the pure-even k-slice (disc already □ ⇒ often yes)."""
    if a == 0 or Fraction(b, a) != k:
        return {"on": False, "reason": "not_on_ray"}
    # need a + 3125 k^4 / 256 = rational square
    val = Fraction(a) + Fraction(3125) * (k**4) / Fraction(256)
    if val < 0:
        return {"on": False, "reason": "negative", "val": str(val)}
    # val = num/den in lowest terms; square iff num and den are squares
    num, den = val.numerator, val.denominator
    n_ok = sp.integer_nthroot(abs(num), 2)[1]
    d_ok = sp.integer_nthroot(den, 2)[1]
    on = bool(n_ok and d_ok)
    m = None
    if on:
        rn = int(sp.integer_nthroot(abs(num), 2)[0])
        rd = int(sp.integer_nthroot(den, 2)[0])
        m = str(Fraction(rn, rd) if num >= 0 else Fraction(-rn, rd))
    return {"on": on, "val": str(val), "m": m, "a": a, "b": b}


def integer_family_specialisations(k: Fraction, seeds: list[dict]) -> dict:
    """
    Build an explicit integer-coeff pure-even family when possible,
    and match known seeds.
    """
    p, q = int(k.numerator), int(k.denominator)
    # General integer model (always Z coeffs for t∈Z):
    # Use α(t) = 256*(q**4)*t**2 - 3125*(p**4)  ONLY when this matches the identity
    # after scaling. Prefer matching via m from seed_on_pure_even_slice.
    #
    # Cleared form that is always pure-even with β/α = p/q when defined:
    #   Let U(t) = 256 t**2 * (256 * q**4) - 3125 * p**4 = 65536 q**4 t**2 - 3125 p**4
    #   Then α_true = U / (256 q**4), β_true = p U / (256 q**5)
    #   f = x^5 + α_true x + β_true has disc □.
    #   Multiply variable: not needed if we specialise U when 256 q^5 | p U... messy.
    #
    # Simpler explicit Z-family when q | 5^∞ (q's primes ⊆ {5}) OR use content:
    #   α(t) = 256 * q**4 * t**2 - 3125 * p**4 // wait wrong scale
    #
    # Correct: α(t) = (65536 * q**4 * t**2 - 3125 * p**4)   as numerator
    # Actually take monic family with integer α,β proportional to true:
    #   α_Z(t) = q * (256 * (256 q**4) t**2 - 3125 p**4) / gcd stuff
    #   α_Z = 256 * 256 * q**5 t**2 - 3125 p**4 * q ? 
    #
    # Flagship method: α = c (r^2 - θ), β = k α with c chosen so α,β ∈ Z[s].
    # k=p/q, θ = 3125 k^4/256 = 3125 p^4/(256 q^4)
    # α = c (s^2 - θ) requires c θ ∈ Z, c ∈ Z.
    # c = 256 q^4 ⇒ α = 256 q^4 s^2 - 3125 p^4 ∈ Z[s], β = (p/q)α
    # β ∈ Z[s] ⇔ q | α for all s ⇔ q | 3125 p^4 ⇔ q | 3125.
    #
    # When q | 3125 (i.e. q | 5^5): integer family
    #   α(s) = 256 q^4 s^2 - 3125 p^4
    #   β(s) = (p/q) α(s)
    # When not: still α(s) integer, β rational — clear by reporting Q-family.
    t = sp.symbols("s", integer=True)
    alpha = 256 * (q**4) * t**2 - 3125 * (p**4)
    beta_exact = sp.Rational(p, q) * alpha
    beta_is_poly_Z = sp.Poly(sp.expand(beta_exact), t, domain=sp.QQ).domain == sp.QQ
    # check if all coeffs of beta are integers
    try:
        bp = sp.Poly(sp.expand(beta_exact), t, domain=sp.QQ)
        beta_Z = all(sp.Rational(c).denominator == 1 for c in bp.all_coeffs())
    except Exception:
        beta_Z = False

    # Disc of (α, β_exact)
    D = sp.expand(256 * alpha**5 + 3125 * beta_exact**4)
    # Compare to square — may fail if wrong scale; use true Q family scale:
    # α_true = 256 (q^2 s)^2 - 3125 k^4/256 = 256 q^4 s^2 - 3125 p^4/(256 q^4)
    alpha_true = sp.together(256 * (q**2 * t) ** 2 - 3125 * sp.Rational(p, q) ** 4 / 256)
    beta_true = sp.together(sp.Rational(p, q) * alpha_true)
    D_true = sp.expand(256 * alpha_true**5 + 3125 * beta_true**4)
    expected = sp.expand((256 * alpha_true**2 * (q**2 * t)) ** 2)
    # m = q^2 s, expected (256 α^2 m)^2
    disc_true_ok = sp.expand(D_true - expected) == 0

    # Match seeds: for each seed compute m^2 = a + 3125 k^4/256
    matches = []
    for s in seeds:
        info = seed_on_pure_even_slice(int(s["a"]), int(s["b"]), k)
        matches.append({**info, "tag": s.get("tag"), "poly": s.get("poly"), "galois": s.get("galois")})

    n_on = sum(1 for m in matches if m.get("on"))
    return {
        "k": str(k),
        "p": p,
        "q": q,
        "alpha_integer_model": str(alpha),
        "beta_integer_model": str(sp.expand(beta_exact)),
        "beta_in_Z_s": beta_Z,
        "alpha_true_Q": str(alpha_true),
        "beta_true_Q": str(beta_true),
        "disc_true_identically_square": disc_true_ok,
        "seed_matches": matches,
        "n_seeds_on_pure_even_family": n_on,
        "multi_seed_pure_even": n_on >= 2,
    }


def analyse_groups(a5_seeds: list[dict]) -> dict:
    groups = group_by_k(a5_seeds)
    multi = {k: v for k, v in groups.items() if len(v) >= 2}
    # Also count unique (a,b) up to sign patterns carefully — keep all
    multi_analysis = []
    for k_str, seeds in sorted(multi.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        k = Fraction(k_str)
        fam = pure_even_family_for_k(k)
        spec = integer_family_specialisations(k, seeds)
        # Deduplicate seeds by (a,b)
        uniq = {(s["a"], s["b"]): s for s in seeds}
        multi_analysis.append(
            {
                "k": k_str,
                "n_seeds": len(uniq),
                "seeds": [
                    {
                        "a": s["a"],
                        "b": s["b"],
                        "tag": s.get("tag"),
                        "poly": s.get("poly"),
                        "galois": s.get("galois"),
                        "status": s.get("status"),
                        "source": s.get("source"),
                    }
                    for s in uniq.values()
                ],
                "family": fam,
                "specialisations": spec,
            }
        )
    multi_pure = [m for m in multi_analysis if m["specialisations"]["multi_seed_pure_even"]]
    return {
        "n_groups": len(groups),
        "n_multi_seed_groups": len(multi),
        "groups_all_sizes": {k: len(v) for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))},
        "multi_seed_groups": multi_analysis,
        "multi_seed_pure_even_slices": multi_pure,
    }


def sample_family_A5(k: Fraction, max_t: int = 12) -> list[dict]:
    """Sample integer points on true Q-family cleared to Z when possible; Gal check."""
    p, q = int(k.numerator), int(k.denominator)
    hits = []
    # Use α = 256 m^2 - 3125 k^4/256 with m = q^2 * t so
    # α = 256 q^4 t^2 - 3125 p^4/(256 q^4)
    # Multiply poly coeffs by den to stay monic Z: if α=A/D, β=B/D with D=256 q^4,
    # monic x^5 + (A/D) x + B/D — not monic Z unless D|A and D|B.
    # Instead evaluate when D divides both:
    Dden = 256 * (q**4)
    for tval in range(1, max_t + 1):
        # m = tval (rational ok)
        # try m = tval and m = tval/(16) etc.
        for m in (
            Fraction(tval),
            Fraction(tval, q),
            Fraction(tval, 16),
            Fraction(tval * q * q, 16),
            Fraction(tval, 256),
            Fraction(tval * 5, 16),
            Fraction(tval, 5),
        ):
            alpha = 256 * (m**2) - Fraction(3125) * (k**4) / 256
            beta = k * alpha
            if alpha.denominator != 1 or beta.denominator != 1:
                continue
            a, b = int(alpha), int(beta)
            if a == 0 or b == 0:
                continue
            d = disc_bj_int(a, b)
            if d <= 0 or not is_square(d):
                continue
            rec = classify_poly(x**5 + a * x + b, do_galois=True)
            if is_a5(rec):
                hits.append(
                    {
                        "m": str(m),
                        "a": a,
                        "b": b,
                        "poly": rec.get("poly"),
                        "galois": rec.get("galois"),
                    }
                )
                if len(hits) >= 6:
                    return hits
    return hits


def main():
    t0 = time.time()
    print("ENLARGE HQCC A5 SEED CATALOGUE — group by k=β/α", flush=True)

    lat = enlarged_lattice(50000)
    print(f"  lattice size: {len(lat)}", flush=True)

    # Bounds: dense small + medium; cap pairs ~ 1.5e6 for runtime
    scan = scan_bj_seeds(lat, alpha_bound=1200, beta_bound=1500, max_pairs=1_800_000)
    print(
        f"  even={scan['n_even']} A5={scan['n_A5']} D5={scan['n_D5']} "
        f"other_even_irr={scan['n_other_even_irr']}",
        flush=True,
    )

    # Also expand with a second medium ring if few A5
    if scan["n_A5"] < 40:
        print("  second pass: medium ring α≤2000, β≤2500 (lattice subset)...", flush=True)
        # denser medium lattice only (not full -2000..2000 every 1)
        med_a = [v for v in lat if abs(v) <= 2000]
        med_b = [v for v in lat if 0 < abs(v) <= 2500]
        seen = {(r["a"], r["b"]) for r in scan["even"]}
        extra_cand = []
        tested2 = 0
        for a in med_a:
            for b in med_b:
                if abs(a) <= 1200 and abs(b) <= 1500:
                    continue  # already scanned
                tested2 += 1
                if (a, b) in seen:
                    continue
                d = disc_bj_int(a, b)
                if d > 0 and is_square(d):
                    extra_cand.append((a, b, d))
                    seen.add((a, b))
                if tested2 > 2_500_000:
                    break
            if tested2 > 2_500_000:
                break
        print(f"    pass2 tested={tested2} new disc□={len(extra_cand)}", flush=True)
        for a, b, d in extra_cand:
            rec = classify_poly(x**5 + a * x + b, do_galois=True)
            row = {
                "a": a,
                "b": b,
                "disc": d,
                "poly": rec.get("poly"),
                "status": rec.get("status"),
                "galois": rec.get("galois"),
                "irreducible": rec.get("irreducible"),
                "source": "scan_pass2",
                "tag": f"a{a}_b{b}",
            }
            scan["even"].append(row)
            if is_a5(row):
                scan["A5"].append(row)
                print(f"    *** A5 pass2 *** α={a} β={b}", flush=True)
        scan["n_even"] = len(scan["even"])
        scan["n_A5"] = len(scan["A5"])
        scan["tested_pairs"] += tested2

    # Deduplicate A5 by (a,b)
    uniq_a5 = {}
    for r in scan["A5"]:
        uniq_a5[(r["a"], r["b"])] = r
    a5_list = list(uniq_a5.values())
    print(f"  unique A5 seeds: {len(a5_list)}", flush=True)

    analysis = analyse_groups(a5_list)
    print(
        f"  k-groups: {analysis['n_groups']}, multi-seed: {analysis['n_multi_seed_groups']}, "
        f"multi pure-even: {len(analysis['multi_seed_pure_even_slices'])}",
        flush=True,
    )

    # Sample extra A5 points on each multi pure-even slice
    for m in analysis["multi_seed_pure_even_slices"]:
        k = Fraction(m["k"])
        samples = sample_family_A5(k, max_t=15)
        m["family_A5_samples"] = samples
        print(f"    k={m['k']}: {m['n_seeds']} catalogue seeds, +{len(samples)} family samples A5", flush=True)

    # Also report multi-seed groups that fail pure-even family condition
    multi_not_pe = [
        m
        for m in analysis["multi_seed_groups"]
        if not m["specialisations"]["multi_seed_pure_even"]
    ]

    elapsed = round(time.time() - t0, 2)

    # ---- markdown ----
    lines = [
        r"# Enlarged HQCC A₅ seed catalogue — grouped by \(k=\beta/\alpha\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** Unique A5 BJ seeds: **{len(a5_list)}**. "
        f"Distinct ratios \(k=\\beta/\\alpha\): **{analysis['n_groups']}**. "
        f"Multi-seed k-groups: **{analysis['n_multi_seed_groups']}**. "
        f"Multi-seed pure-even slices: **{len(analysis['multi_seed_pure_even_slices'])}**.",
        "",
        "---",
        "",
        "## Method",
        "",
        "1. Build enlarged HQCC lattice (model core, powers of 3, short combos, dense small integers).",
        r"2. Scan Bring–Jerrard \(x^5+\alpha x+\beta\): keep disc perfect square; classify Gal.",
        r"3. Retain irreducible \(\mathrm{Gal}=A_5\) seeds.",
        r"4. Group by reduced rational \(k=\beta/\alpha\).",
        r"5. For each multi-seed \(k\), test the pure-even family condition",
        r"   \(\alpha + 3125 k^4/256 = \square_{\mathbb{Q}}\) and record the LSW-type family.",
        "",
        r"General pure-even family on the ray \(\beta=k\alpha\):",
        "",
        r"$$\alpha(m)=256 m^2 - \frac{3125\,k^4}{256},\qquad \beta(m)=k\cdot\alpha(m)$$",
        "",
        r"$$\operatorname{disc}=(256\,\alpha(m)^2\, m)^2\quad\text{(identically square in }\mathbb{Q}(m)).$$",
        "",
        "---",
        "",
        "## Scan stats",
        "",
        f"| quantity | value |",
        f"|----------|------:|",
        f"| lattice size | {len(lat)} |",
        f"| pairs tested | {scan['tested_pairs']} |",
        f"| disc□ (all sources) | {scan['n_even']} |",
        f"| A5 irreducible | {len(a5_list)} |",
        f"| D5 (scan label) | {scan['n_D5']} |",
        f"| other even irr | {scan['n_other_even_irr']} |",
        f"| k-groups | {analysis['n_groups']} |",
        f"| multi-seed k-groups | {analysis['n_multi_seed_groups']} |",
        f"| multi-seed pure-even slices | {len(analysis['multi_seed_pure_even_slices'])} |",
        "",
        "---",
        "",
        "## Multi-seed pure-even slices (primary output)",
        "",
    ]

    if not analysis["multi_seed_pure_even_slices"]:
        lines.append("_None found beyond grouping failures — see multi-seed groups below._")
    for m in analysis["multi_seed_pure_even_slices"]:
        lines.append(rf"### \(k = {m['k']}\)")
        lines.append("")
        lines.append(f"- Catalogue A5 seeds on ray: **{m['n_seeds']}**")
        lines.append(f"- Disc identically square (true Q-family): **{m['specialisations']['disc_true_identically_square']}**")
        lines.append(f"- Integer β-model in Z[s]: **{m['specialisations']['beta_in_Z_s']}**")
        lines.append(f"- α_true = `{m['specialisations']['alpha_true_Q']}`")
        lines.append(f"- β_true = `{m['specialisations']['beta_true_Q']}`")
        lines.append("")
        lines.append("| α | β | poly | on family | m |")
        lines.append("|--:|--:|------|:---------:|---|")
        for s in m["specialisations"]["seed_matches"]:
            lines.append(
                f"| {s.get('a')} | {s.get('b')} | `{s.get('poly','')}` | {s.get('on')} | {s.get('m')} |"
            )
        lines.append("")
        samples = m.get("family_A5_samples") or []
        if samples:
            lines.append("Extra A5 specialisations on the pure-even family:")
            for s in samples:
                lines.append(f"- m={s['m']}: `{s['poly']}`")
            lines.append("")

    lines += [
        "---",
        "",
        "## All multi-seed k-groups (A5)",
        "",
    ]
    for m in analysis["multi_seed_groups"]:
        pe = m["specialisations"]["multi_seed_pure_even"]
        lines.append(
            f"- **k={m['k']}**: {m['n_seeds']} seeds; pure-even multi: **{pe}** — "
            + ", ".join(f"({s['a']},{s['b']})" for s in m["seeds"])
        )
    if multi_not_pe:
        lines += [
            "",
            "Multi-seed groups failing pure-even family condition (disc□ seeds not on LSW param):",
            "",
        ]
        for m in multi_not_pe:
            lines.append(f"- k={m['k']}: matches={m['specialisations']['seed_matches']}")

    lines += [
        "",
        "---",
        "",
        "## Full A5 catalogue (unique)",
        "",
        "| α | β | k=β/α | poly | source |",
        "|--:|--:|-------|------|--------|",
    ]
    for s in sorted(a5_list, key=lambda r: (str(k_of(r["a"], r["b"]) or ""), r["a"], r["b"])):
        kk = k_of(s["a"], s["b"])
        lines.append(
            f"| {s['a']} | {s['b']} | {kk} | `{s.get('poly')}` | {s.get('source')} |"
        )

    lines += [
        "",
        "---",
        "",
        "## k-group size histogram",
        "",
        "```",
        str(analysis["groups_all_sizes"]),
        "```",
        "",
        "---",
        "",
        "## Conclusions",
        "",
        "1. Enlarged catalogue yields more A5 BJ lattice seeds to feed k-grouping.",
        "2. Every multi-seed pure-even slice is an LSW-type family on a fixed ray β=kα;",
        "   LSW itself is k=-4; flagship is the rational-k family k=-8/5.",
        "3. New multi-seed pure-even slices (if any above) are the natural fusion fuel:",
        "   one pure-even 1-parameter family carrying ≥2 HQCC A5 seeds.",
        "4. Still open for geometric fusion: a pure-even family joining flagship to a",
        "   seed with a *different* k (not on the same ray).",
        "",
        "_Generated by enlarge_seed_catalogue.py_",
    ]

    doc = "\n".join(lines)
    blob = {
        "elapsed_sec": elapsed,
        "lattice_size": len(lat),
        "scan": {
            "tested_pairs": scan["tested_pairs"],
            "n_even": scan["n_even"],
            "n_A5": len(a5_list),
            "n_D5": scan["n_D5"],
        },
        "A5_seeds": [
            {
                "a": s["a"],
                "b": s["b"],
                "k": str(k_of(s["a"], s["b"])),
                "poly": s.get("poly"),
                "galois": s.get("galois"),
                "status": s.get("status"),
                "source": s.get("source"),
                "tag": s.get("tag"),
            }
            for s in a5_list
        ],
        "analysis": {
            "n_groups": analysis["n_groups"],
            "n_multi_seed_groups": analysis["n_multi_seed_groups"],
            "groups_all_sizes": analysis["groups_all_sizes"],
            "multi_seed_pure_even_slices": analysis["multi_seed_pure_even_slices"],
            "multi_seed_groups": analysis["multi_seed_groups"],
        },
        "verdict": (
            f"A5={len(a5_list)}, k-groups={analysis['n_groups']}, "
            f"multi-seed={analysis['n_multi_seed_groups']}, "
            f"multi pure-even={len(analysis['multi_seed_pure_even_slices'])}"
        ),
    }

    write_md(OUT / "ENLARGED_SEED_CATALOGUE.md", doc)
    write_md(RESULTS / "ENLARGED_SEED_CATALOGUE.md", doc)
    write_md(ROOT / "ENLARGED_SEED_CATALOGUE.md", doc)
    write_json(OUT / "ENLARGED_SEED_CATALOGUE.json", blob)
    print(blob["verdict"], flush=True)
    print(f"Wrote ENLARGED_SEED_CATALOGUE.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
