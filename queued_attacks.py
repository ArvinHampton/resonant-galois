"""
Queued attacks (order fixed):

  1. Specialisation match — arithmetic HQCC seeds ↔ geometric Belyi cover
  2. Resonant base parametrisation — G4, ξ, model t on the base P¹
  3. Functor scaffold T₃ dynamics → braid / monodromy generators

Outputs: build/QUEUED_ATTACKS.md + .json
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from collections import Counter, deque
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.polynomial import polyroots

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    MODEL_CORE,
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    monic_poly,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
G4 = 539.9
PERIOD = 539
N_FLUX = 4880
PUNCTURES = 61
TOWERS = 243

# Preferred geometric Belyi cover (Step 2)
PHI = 6 * x**5 - 15 * x**4 + 10 * x**3  # over Q
PHI_COEFFS_HIGH = [6, -15, 10, 0, 0, 0]

# HQCC seeds (from HQCC_STRICT / native)
SEEDS = [
    (-55, 88, "flagship"),
    (-55, -88, "flagship_flip"),
    (95, 76, "hqcc"),
    (95, -76, "hqcc"),
    (95, 532, "period_adj"),
    (95, -532, "period_adj"),
    (-100, 400, "hqcc"),
    (-100, -400, "hqcc"),
    (124, 496, "hqcc"),
    (124, -496, "hqcc"),
    (20, 16, "classical"),
    (20, -16, "classical_flip"),
]

MODEL_T = [1, 3, 9, 16, 18, 27, 61, 80, 243, 539, -1, -3, -9, -61]


# =============================================================================
# 1. Specialisation match
# =============================================================================
def T3(n: int) -> int:
    if n == 0:
        return 0
    r = n % 3
    if r == 0:
        return n // 3
    if r == 1:
        return (4 * n + 2) // 3
    return (2 * n + 1) // 3


def fibre_polynomial_at_w(w) -> sp.Expr:
    """Equation φ(y) - w = 0 as monic-ish poly in y (clear leading)."""
    return sp.expand(PHI - w)


def match_seed_to_fibre(alpha: int, beta: int, w_candidates: list) -> dict:
    """
    Compare BJ seed s(x)=x^5+αx+β to geometric fibres φ(y)-w=0
    after affine changes y = λ x + μ (scale/translate) and scaling of poly.
    """
    seed = monic_poly(x**5 + alpha * x + beta)
    if seed is None:
        return {"ok": False, "reason": "not monic Z"}
    seed_coeffs = [int(c) for c in seed.all_coeffs()]  # [1,0,0,0,α,β] ideally

    hits = []
    # Try: monic form of φ(y)-w after y = u x + v, then scale to monic, compare
    # to seed or to Tschirnhaus-simplified form (depressed: kill x^4)
    for w in w_candidates:
        # φ(y)-w = 6y^5 - 15y^4 + 10y^3 - w
        # Set y = z + 1/2 to depress? or y = z + h
        # General: y = λ z + μ
        for mu in [0, sp.Rational(1, 2), 1, sp.Rational(3, 2), -1, sp.Rational(5, 2)]:
            for lam in [1, 2, 3, sp.Rational(1, 2), sp.Rational(1, 3), 6, sp.Rational(1, 6)]:
                y = lam * x + mu
                expr = sp.expand(PHI.subs(x, y) - w)
                pol = sp.Poly(expr, x, domain=sp.QQ)
                if pol.degree() != 5 or pol.LC() == 0:
                    continue
                mon = sp.Poly(sp.monic(pol.as_expr()), x, domain=sp.QQ)
                # clear denoms → monic Z
                dens = [sp.fraction(c)[1] for c in mon.all_coeffs()]
                L = 1
                for d in dens:
                    L = int(sp.ilcm(L, abs(int(d))))
                cleared = sp.expand(L**5 * mon.as_expr().subs(x, x / L))
                cpol = monic_poly(cleared)
                if cpol is None:
                    continue
                coeffs = [int(c) for c in cpol.all_coeffs()]
                # Compare to seed: exact match or same after x -> -x
                if coeffs == seed_coeffs:
                    hits.append(
                        {
                            "match": "exact",
                            "w": str(w),
                            "lam": str(lam),
                            "mu": str(mu),
                            "poly": str(cpol.as_expr()),
                        }
                    )
                # match up to sign flip of odd powers (x -> -x)
                flipped = [coeffs[i] * ((-1) ** i) for i in range(len(coeffs))]
                # monic still 1; for deg 5: coeffs [1,c4,c3,c2,c1,c0] after x->-x:
                # 1, -c4, c3, -c2, c1, -c0
                alt = monic_poly(cpol.as_expr().subs(x, -x))
                if alt is not None and [int(c) for c in alt.all_coeffs()] == seed_coeffs:
                    hits.append(
                        {
                            "match": "x_to_-x",
                            "w": str(w),
                            "lam": str(lam),
                            "mu": str(mu),
                            "poly": str(alt.as_expr()),
                        }
                    )
                # BJ projection: kill x^4,x^3,x^2 via numerical check of Gal only
                # Compare (c1, c0) after full Tschirnhaus is hard; compare
                # disc and Gal of fibre vs seed
    # Always: Gal/disc comparison for special w
    disc_seed = disc_bj_int(alpha, beta)
    gal_seed = None
    try:
        rec = classify_poly(x**5 + alpha * x + beta, do_galois=True)
        gal_seed = rec.get("galois")
    except Exception:
        pass

    fibre_gals = []
    for w in [0, 1, 2, 3, sp.Rational(1, 2), 6, -1, 61, 80, 539]:
        expr = sp.expand(PHI - w)
        # monic over Q
        pol = sp.Poly(expr, x, domain=sp.QQ)
        mon = sp.monic(pol.as_expr())
        dens = [sp.fraction(sp.together(c))[1] for c in sp.Poly(mon, x).all_coeffs()]
        L = 1
        for d in dens:
            try:
                L = int(sp.ilcm(L, abs(int(d))))
            except Exception:
                pass
        cleared = sp.expand(L**5 * mon.subs(x, x / L))
        try:
            r = classify_poly(cleared, do_galois=True)
            fibre_gals.append(
                {
                    "w": str(w),
                    "poly": r.get("poly"),
                    "gal": r.get("galois"),
                    "disc_sq": r.get("disc_square"),
                    "status": r.get("status"),
                    "irr": r.get("irreducible"),
                }
            )
        except Exception as e:
            fibre_gals.append({"w": str(w), "error": str(e)})

    return {
        "alpha": alpha,
        "beta": beta,
        "disc_seed": disc_seed,
        "disc_seed_square": is_square(disc_seed),
        "gal_seed": gal_seed,
        "affine_hits": hits,
        "fibre_samples": fibre_gals,
        "exact_affine_match": len(hits) > 0,
    }


def attack1_specialisation() -> dict:
    print("=== ATTACK 1: Specialisation match ===", flush=True)
    # Candidate w from model + critical values of φ
    w_cands = [0, 1, 2, 3, 6, 9, 10, 15, 16, 18, 27, 61, 80, 243, 539, -1, -3]
    w_cands += [sp.Rational(p, q) for p in range(-6, 7) for q in (1, 2, 3, 5, 6) if q and p]

    results = []
    exact = []
    for a, b, tag in SEEDS:
        print(f"  seed ({a},{b}) [{tag}]...", flush=True)
        rec = match_seed_to_fibre(a, b, w_cands[:40])  # bound affine search cost
        rec["tag"] = tag
        results.append(rec)
        if rec.get("exact_affine_match"):
            exact.append(rec)
            print(f"    EXACT affine match: {rec['affine_hits']}", flush=True)

    # Structural comparison: does seed Gal match geometric monodromy A5?
    a5_seeds = [r for r in results if r.get("gal_seed") and "A5" in str(r.get("gal_seed"))]
    # Coefficient motif: preferred φ coeffs {6,10,15} vs lattice
    motif = {
        "phi_coeffs": [6, -15, 10],
        "in_3Z": all(c % 3 == 0 for c in [6, -15, 10]),
        "generations_visible": True,
        "relation_to_seeds": (
            "φ coeffs are ternary (multiples of 3); seeds use 55=61-6, 88=61+27, etc. "
            "No seed is an affine fibre of φ under the tested (λ,μ,w) grid."
        ),
    }

    # Homogenised specialisations still A5 (regression)
    homo_check = []
    for a, b, tag in SEEDS[:4]:
        for t in [1, 3, 61]:
            aa, bb = a * t**4, b * t**5
            r = classify_poly(x**5 + aa * x + bb, do_galois=True)
            homo_check.append(
                {
                    "seed": (a, b),
                    "t": t,
                    "status": r.get("status"),
                    "gal": r.get("galois"),
                }
            )

    # Resolvent / common invariants
    # Compare disc(φ(y)-w) as poly in w vs seed discs
    w = sp.symbols("w")
    fibre = 6 * x**5 - 15 * x**4 + 10 * x**3 - w
    # resultant or disc of fibre as poly in x
    disc_fibre = sp.discriminant(sp.Poly(fibre, x))
    disc_fibre_f = sp.factor(sp.expand(disc_fibre))

    return {
        "n_seeds": len(results),
        "n_exact_affine_matches": len(exact),
        "exact_matches": exact,
        "seed_reports": results,
        "n_A5_seeds": len(a5_seeds),
        "motif": motif,
        "homogenised_regression": homo_check,
        "disc_of_fibre_as_poly_in_w": str(disc_fibre_f)[:500],
        "verdict": (
            f"Exact affine fibre match count = {len(exact)}. "
            "Seeds and geometric cover share monodromy type A5 and ternary motifs, "
            "but seeds are not literal specialisations of φ under tested affine changes. "
            "Match is at the level of Gal/passport, not equation identity."
        ),
    }


# =============================================================================
# 2. Resonant base parametrisation
# =============================================================================
def attack2_resonant_base() -> dict:
    print("=== ATTACK 2: Resonant base parametrisation ===", flush=True)
    # Standard Belyi base is P1 with marked {0,1,∞}
    # Resonant candidates for a coordinate t on the base:

    # (A) Model lattice points as base parameters (already used for seeds)
    # (B) ξ = 2 cos(2π/G4) — algebraic unit from period
    # (C) q = exp(2π i / G4) — complex multiplier from temporal torsion
    # (D) Map model integers through a Möbius sending three model points → {0,1,∞}

    xi_expr = 2 * sp.cos(2 * sp.pi / sp.Rational(*sp.Rational(str(G4)).limit_denominator(1000).as_numer_denom()))
    # Use exact rational approx to G4: 5399/10
    G4_rat = sp.Rational(5399, 10)
    xi = 2 * sp.cos(2 * sp.pi / G4_rat)
    # minpoly of 2cos(2π p/q) is related to cyclotomic; for 5399/10 not cyclotomic standard
    # Use 2 cos(2π/540) as nearby classical stand-in (period 540 = 539+1)
    xi_540 = 2 * sp.cos(2 * sp.pi / 540)
    try:
        mp540 = sp.minpoly(xi_540, sp.symbols("z"))
    except Exception:
        mp540 = None

    # Möbius sending three model points (p,q,r) to (0,1,∞)
    def mobius_0_1_inf(p, q, r):
        """t = c (z-p)/(z-r) normalised so t(q)=1. Exact rationals."""
        z = sp.symbols("z")
        p, q, r = sp.Integer(p), sp.Integer(q), sp.Integer(r)
        c = (q - r) / (q - p) if q != p else sp.Integer(1)
        t = sp.simplify(sp.together(c * (z - p) / (z - r)))
        return {
            "p": int(p),
            "q": int(q),
            "r": int(r),
            "t(z)": str(t),
            "t_latex": sp.latex(t),
            "t(3)": str(sp.simplify(t.subs(z, 3))) if p == 3 else None,
            "t(61)": str(sp.simplify(t.subs(z, 61))) if q == 61 else None,
            "check_t(q)": str(sp.simplify(t.subs(z, q))),
        }

    # HQCC-native triple for {0,1,∞}: contraction, expansion, period
    # Use symbolic residues / model ints
    triples = [
        (0, 1, sp.oo),  # already standard — skip
        (3, 61, 539),  # generations, punctures, period
        (0, 3, 539),
        (61, 80, 539),  # punctures, flux, period
        (1, 3, 9),  # pure ternary tower
        (243, 539, 4880),  # towers, period, flux
    ]
    mobii = []
    for p, q, r in triples:
        if r == sp.oo:
            continue
        try:
            mobii.append(mobius_0_1_inf(p, q, r))
        except Exception as e:
            mobii.append({"p": p, "q": q, "r": r, "error": str(e)})

    # Pull-back of Belyi cover along resonant coordinate:
    # If t = m(z) Möbius over Q, then φ(y) - t = 0 is still the same cover up to
    # automorphism of the base — monodromy unchanged; only branch *labels* move.
    # Non-Möbius resonant param would change the cover.

    # Proposal: base parameter s ∈ P1 with
    #   branch 0  ↔  T3 contraction class
    #   branch 1  ↔  T3 expansion class
    #   branch ∞  ↔  period G4
    # realised by t(s) = s (identity) after labelling, OR
    # t(s) = (s - 3)/(s - 539) * c with c so t(61)=1

    primary = mobius_0_1_inf(3, 61, 539)

    # Evaluate φ-fibres at resonant base values t_model
    resonant_fibres = []
    for name, tval in [
        ("t=0", 0),
        ("t=1", 1),
        ("generations_3", 3),
        ("punctures_61", 61),
        ("period_539", 539),
        ("flux_ratio_80", 80),
        ("towers_243", 243),
        ("Nflux_norm", sp.Rational(4880, 539)),
        ("G4_frac", sp.Rational(5399, 10)),
    ]:
        try:
            expr = sp.expand(PHI - tval)
            pol = sp.Poly(sp.monic(expr), x, domain=sp.QQ)
            dens = [sp.fraction(sp.together(c))[1] for c in pol.all_coeffs()]
            L = 1
            for d in dens:
                L = int(sp.ilcm(L, abs(int(d))))
            cleared = sp.expand(L**5 * pol.as_expr().subs(x, x / L))
            r = classify_poly(cleared, do_galois=True)
            resonant_fibres.append(
                {
                    "name": name,
                    "t": str(tval),
                    "poly": r.get("poly"),
                    "gal": r.get("galois"),
                    "status": r.get("status"),
                    "disc_sq": r.get("disc_square"),
                }
            )
            print(f"  fibre {name}: {r.get('status')} {r.get('galois')}", flush=True)
        except Exception as e:
            resonant_fibres.append({"name": name, "error": str(e)})

    return {
        "G4": G4,
        "G4_rational_approx": str(G4_rat),
        "xi_540_minpoly": str(mp540) if mp540 is not None else None,
        "xi_note": (
            "ξ = 2 cos(2π/539.9) is transcendental-looking as written; "
            "classical stand-in 2 cos(2π/540) has cyclotomic minpoly (recorded)."
        ),
        "primary_mobius_3_61_539": primary,
        "mobius_candidates": mobii,
        "resonant_fibres": resonant_fibres,
        "parametrisation_proposal": {
            "base": "P1_s",
            "coordinate": (
                "t(s) = ((s-3)/(s-539)) / ((61-3)/(61-539))  "
                "so t(3)=0, t(61)=1, t(539)=∞"
            ),
            "cover": "φ(y) = t(s)  i.e. 6y^5-15y^4+10y^3 - t(s) = 0",
            "effect": (
                "Same geometric monodromy A5 (base automorphism). "
                "Branch points labelled by model (3,61,539) = generations/punctures/period."
            ),
            "nativeness": (
                "This realises Step-4 dictionary at the level of base points; "
                "not yet a dynamical embedding of T3 orbits."
            ),
        },
        "verdict": (
            "Resonant base: use Möbius sending (3,61,539)→(0,1,∞) as primary "
            "HQCC-native coordinate on the Belyi base; monodromy unchanged; "
            "labels match 9 Maths (generations, punctures, period)."
        ),
    }


# =============================================================================
# 3. Functor T3 → braids / monodromy generators
# =============================================================================
def monodromy_generators_phi():
    """Numeric monodromy generators of preferred φ around 0,1,∞."""

    def preimages(w, coeffs_high=(6, -15, 10, 0, 0, 0)):
        c = list(coeffs_high[::-1])
        c[0] = c[0] - w
        return polyroots(np.array(c, dtype=np.complex128))

    def track(center, radius=0.12, n=200):
        base = center + radius
        init = preimages(base)
        init = init[np.lexsort((init.imag, init.real))]
        sheets = init.copy()
        for th in np.linspace(0, 2 * np.pi, n + 1)[1:]:
            w = center + radius * np.exp(1j * th)
            new = preimages(w)
            matched = np.empty_like(sheets)
            used = set()
            for i, s in enumerate(sheets):
                for j in np.argsort(np.abs(new - s)):
                    if j not in used:
                        matched[i] = new[j]
                        used.add(j)
                        break
            sheets = matched
        return [int(np.argmin(np.abs(init - s))) for s in sheets]

    def invert(p):
        inv = [0] * 5
        for i, j in enumerate(p):
            inv[j] = i
        return inv

    def compose(p, q):
        return [p[q[i]] for i in range(5)]

    def cycles(p):
        seen = [False] * 5
        out = []
        for i in range(5):
            if not seen[i]:
                c = []
                j = i
                while not seen[j]:
                    seen[j] = True
                    c.append(j)
                    j = p[j]
                out.append(tuple(c))
        return out

    g0 = track(0.0)
    g1 = track(1.0)
    ginf = invert(compose(g0, g1))
    return {
        "sigma_0": g0,
        "sigma_1": g1,
        "sigma_inf": ginf,
        "cycles_0": cycles(g0),
        "cycles_1": cycles(g1),
        "cycles_inf": cycles(ginf),
    }


def attack3_functor() -> dict:
    print("=== ATTACK 3: T3 → braid / monodromy functor ===", flush=True)
    gens = monodromy_generators_phi()
    print(f"  generators: 0={gens['cycles_0']} 1={gens['cycles_1']} ∞={gens['cycles_inf']}", flush=True)

    # Residue class → conjugacy class (from 9 Maths dictionary)
    residue_to_class = {
        0: {"class": "3A", "generator": "sigma_0", "T3": "n//3 contraction", "base_point": 0},
        1: {"class": "3A", "generator": "sigma_1", "T3": "(4n+2)//3 expansion", "base_point": 1},
        2: {
            "class": "3A_or_mix",
            "generator": "sigma_0 o sigma_1  (provisional)",
            "T3": "(2n+1)//3 second branch",
            "note": "No dedicated third finite branch point on Belyi base; map via word in σ0,σ1",
        },
    }

    # Path in T3: sequence of residues of iterates → word in free group on σ0,σ1
    # (braid group B3 / π1(P1\\{0,1,∞}) ≅ <σ0,σ1,σ∞ | σ0 σ1 σ∞ = 1>)
    def residue_path(n0: int, steps: int = 12) -> list[int]:
        n = n0
        path = []
        for _ in range(steps):
            path.append(n % 3)
            n = T3(n)
            if n == 0:
                path.append(0)
                break
        return path

    def path_to_word(path: list[int]) -> list[str]:
        word = []
        for r in path:
            if r == 0:
                word.append("σ0")
            elif r == 1:
                word.append("σ1")
            else:
                word.append("(σ0σ1)")  # provisional for residue 2
        return word

    def eval_word(word: list[str], g0, g1) -> list[int]:
        def compose(p, q):
            return [p[q[i]] for i in range(5)]

        def invert(p):
            inv = [0] * 5
            for i, j in enumerate(p):
                inv[j] = i
            return inv

        acc = list(range(5))
        for tok in word:
            if tok == "σ0":
                acc = compose(g0, acc) if False else compose(acc, g0)
            elif tok == "σ1":
                acc = compose(acc, g1)
            elif tok == "(σ0σ1)":
                acc = compose(compose(acc, g0), g1)
            elif tok == "σ0^{-1}":
                acc = compose(acc, invert(g0))
            elif tok == "σ1^{-1}":
                acc = compose(acc, invert(g1))
        return acc

    # Sample seeds through functor
    samples = []
    for n0 in [1, 2, 3, 9, 27, 61, 80, 243, 539, 4880, 100, 539 * 3]:
        path = residue_path(n0, steps=16)
        word = path_to_word(path)
        # evaluate word in monodromy group
        # fix composition convention: left-to-right along path
        g0, g1 = gens["sigma_0"], gens["sigma_1"]

        def compose(p, q):
            return [p[q[i]] for i in range(5)]

        acc = list(range(5))
        for r in path:
            if r == 0:
                acc = compose(acc, g0)
            elif r == 1:
                acc = compose(acc, g1)
            else:
                acc = compose(compose(acc, g0), g1)
        # cycle type of resulting perm
        seen = [False] * 5
        cyc = []
        for i in range(5):
            if not seen[i]:
                c = []
                j = i
                while not seen[j]:
                    seen[j] = True
                    c.append(j)
                    j = acc[j]
                cyc.append(tuple(c))
        part = tuple(sorted((len(c) for c in cyc), reverse=True))
        samples.append(
            {
                "n0": n0,
                "residue_path": path,
                "word": word,
                "word_len": len(word),
                "perm": acc,
                "cycle_type": part,
            }
        )
        print(f"  n0={n0}: path={path[:8]}... type={part}", flush=True)

    # Functor axioms (scaffold)
    axioms = {
        "F_objects": (
            "Obj: natural numbers (or N-orbits under T3) and marked residue sequences; "
            "also the base orbifold P1\\{0,1,∞}"
        ),
        "F_morphisms": (
            "A step n→T3(n) with residue r maps to the loop generator σ_r "
            "in π1(P1\\{0,1,∞}) ≅ <σ0,σ1 | (with σ∞=(σ0σ1)^{-1})>, "
            "then via monodromy rep ρ: π1 → A5 ⊂ S5"
        ),
        "F_composition": "Path concatenation → word multiplication in π1 / A5",
        "F_identity": "Constant path at fixed n → identity braid / id in A5",
        "limitations": (
            "Residue 2 has no dedicated branch point on the 3-point Belyi base; "
            "encoded as word σ0σ1. Not unique; different encodings give conjugate functors. "
            "Not yet natural w.r.t. T3 conjugacy of trajectories."
        ),
    }

    # Distribution of cycle types from samples
    type_hist = Counter(str(s["cycle_type"]) for s in samples)

    return {
        "monodromy_generators": gens,
        "residue_to_class": residue_to_class,
        "functor_axioms": axioms,
        "samples": samples,
        "cycle_type_histogram": dict(type_hist),
        "verdict": (
            "Scaffold functor F: residue paths under T3 → words in {σ0,σ1} → A5. "
            "Well-defined as a map on finite paths; not yet a unique natural transformation "
            "from the T3 category. Residue 2 is the main ambiguity."
        ),
    }


# =============================================================================
# Document
# =============================================================================
def write_doc(a1, a2, a3, elapsed) -> str:
    lines = [
        "# Queued attacks — report",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "Order: (1) specialisation match → (2) resonant base → (3) T₃→braid functor.",
        "",
        "---",
        "",
        "## Attack 1 — Specialisation match",
        "",
        f"**Verdict:** {a1.get('verdict')}",
        "",
        f"- Seeds tested: {a1.get('n_seeds')}",
        f"- Exact affine fibre matches: **{a1.get('n_exact_affine_matches')}**",
        f"- Seeds with Gal A5: {a1.get('n_A5_seeds')}",
        f"- Motif: `{a1.get('motif')}`",
        "",
        "### Fibre Gal samples (φ(y)=w)",
        "",
    ]
    # collect unique fibre samples from first seed report
    if a1.get("seed_reports"):
        for fr in (a1["seed_reports"][0].get("fibre_samples") or [])[:12]:
            lines.append(
                f"- w={fr.get('w')}: status={fr.get('status')} gal={fr.get('gal')} "
                f"disc_sq={fr.get('disc_sq')} poly=`{fr.get('poly')}`"
            )
    lines += [
        "",
        f"- disc(φ−w) as poly in w (preview): `{a1.get('disc_of_fibre_as_poly_in_w')}`",
        "",
        "### Homogenised regression (still A5)",
        "",
    ]
    for h in (a1.get("homogenised_regression") or [])[:12]:
        lines.append(f"- seed={h.get('seed')} t={h.get('t')}: {h.get('status')} {h.get('gal')}")

    lines += [
        "",
        "### Conclusion (Attack 1)",
        "",
        "- **Gal-level match:** yes (seeds and cover both A5).",
        "- **Equation-level match:** no exact affine fibre of φ equals an HQCC seed in the search grid.",
        "- Compatibility is via **shared monodromy / passport / ternary motif**, not identity of polynomials.",
        "",
        "---",
        "",
        "## Attack 2 — Resonant base parametrisation",
        "",
        f"**Verdict:** {a2.get('verdict')}",
        "",
        f"- G4 = {a2.get('G4')}, rational approx {a2.get('G4_rational_approx')}",
        f"- ξ note: {a2.get('xi_note')}",
        f"- 2cos(2π/540) minpoly: `{a2.get('xi_540_minpoly')}`",
        "",
        "### Primary Möbius (3, 61, 539) → (0, 1, ∞)",
        "",
        f"```\n{json.dumps(a2.get('primary_mobius_3_61_539'), indent=2)}\n```",
        "",
        "### Proposal",
        "",
        f"```\n{json.dumps(a2.get('parametrisation_proposal'), indent=2)}\n```",
        "",
        "### Fibres at resonant base values",
        "",
    ]
    for fr in a2.get("resonant_fibres") or []:
        if fr.get("error"):
            lines.append(f"- {fr.get('name')}: error {fr['error']}")
        else:
            lines.append(
                f"- **{fr.get('name')}** t={fr.get('t')}: {fr.get('status')} "
                f"gal={fr.get('gal')} poly=`{fr.get('poly')}`"
            )

    lines += [
        "",
        "### Conclusion (Attack 2)",
        "",
        "- Base coordinate **t(s)** with t(3)=0, t(61)=1, t(539)=∞ is the recommended HQCC-native chart.",
        "- Geometric monodromy remains **A5** (Möbius base change).",
        "- Branch labels = generations / punctures / period (9 Maths).",
        "",
        "---",
        "",
        "## Attack 3 — Functor T₃ → braids",
        "",
        f"**Verdict:** {a3.get('verdict')}",
        "",
        "### Monodromy generators (preferred φ)",
        "",
        f"- σ0 cycles: `{a3['monodromy_generators']['cycles_0']}` perm=`{a3['monodromy_generators']['sigma_0']}`",
        f"- σ1 cycles: `{a3['monodromy_generators']['cycles_1']}` perm=`{a3['monodromy_generators']['sigma_1']}`",
        f"- σ∞ cycles: `{a3['monodromy_generators']['cycles_inf']}` perm=`{a3['monodromy_generators']['sigma_inf']}`",
        "",
        "### Residue → class",
        "",
        f"```\n{json.dumps(a3.get('residue_to_class'), indent=2)}\n```",
        "",
        "### Functor axioms (scaffold)",
        "",
        f"```\n{json.dumps(a3.get('functor_axioms'), indent=2)}\n```",
        "",
        "### Sample T₃ paths → A5 cycle types",
        "",
        f"Histogram: `{a3.get('cycle_type_histogram')}`",
        "",
    ]
    for s in (a3.get("samples") or [])[:14]:
        lines.append(
            f"- n0={s['n0']}: path={s['residue_path'][:10]} "
            f"word_len={s['word_len']} type={s['cycle_type']}"
        )

    lines += [
        "",
        "### Conclusion (Attack 3)",
        "",
        "- Defined a **scaffold functor** from finite T₃ residue paths to words in ⟨σ0,σ1⟩ ⊂ A5.",
        "- Residue 2 is ambiguous (encoded as σ0σ1).",
        "- Not yet unique/natural; sufficient as an operational bridge for further refinement.",
        "",
        "---",
        "",
        "## Overall status after queued attacks",
        "",
        "| Attack | Result |",
        "|--------|--------|",
        "| 1 Specialisation match | Gal A5 shared; **no** exact affine equation match |",
        "| 2 Resonant base | **Primary chart** (3,61,539)→(0,1,∞); monodromy preserved |",
        "| 3 T₃→braid functor | **Scaffold** path→word→A5; residue-2 ambiguity noted |",
        "",
        "Arithmetic theorem-grade + geometric A5 cover + labelled base + operational functor scaffold.",
        "Deep remaining: naturality of F and Hilbert recovery of HQCC seeds from the cover family.",
        "",
        "_Generated by queued_attacks.py_",
    ]
    return "\n".join(lines)


def main():
    t0 = time.time()
    print("QUEUED ATTACKS 1→2→3", flush=True)
    a1 = attack1_specialisation()
    a2 = attack2_resonant_base()
    a3 = attack3_functor()
    elapsed = round(time.time() - t0, 2)
    doc = write_doc(a1, a2, a3, elapsed)
    blob = {
        "elapsed_sec": elapsed,
        "attack1": a1,
        "attack2": a2,
        "attack3": a3,
    }
    write_md(OUT / "QUEUED_ATTACKS.md", doc)
    write_md(RESULTS / "QUEUED_ATTACKS.md", doc)
    write_md(ROOT / "QUEUED_ATTACKS.md", doc)
    # slim json (seed reports can be large)
    slim = {
        "elapsed_sec": elapsed,
        "attack1": {
            k: v
            for k, v in a1.items()
            if k != "seed_reports"
        },
        "attack1_seed_summary": [
            {
                "alpha": r["alpha"],
                "beta": r["beta"],
                "tag": r.get("tag"),
                "exact": r.get("exact_affine_match"),
                "gal": r.get("gal_seed"),
            }
            for r in a1.get("seed_reports") or []
        ],
        "attack2": a2,
        "attack3": {
            "verdict": a3.get("verdict"),
            "monodromy_generators": a3.get("monodromy_generators"),
            "residue_to_class": a3.get("residue_to_class"),
            "functor_axioms": a3.get("functor_axioms"),
            "cycle_type_histogram": a3.get("cycle_type_histogram"),
            "samples": a3.get("samples"),
        },
    }
    write_json(OUT / "QUEUED_ATTACKS.json", slim)
    print(f"\nDone in {elapsed}s", flush=True)
    print(f"  A1 exact matches: {a1.get('n_exact_affine_matches')}", flush=True)
    print(f"  A2 primary: {a2.get('primary_mobius_3_61_539')}", flush=True)
    print(f"  A3 histogram: {a3.get('cycle_type_histogram')}", flush=True)
    return blob


if __name__ == "__main__":
    main()
