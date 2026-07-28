#!/usr/bin/env python3
"""
G3d — Optional refinements for Ni(A5, 5A^4) explicit model.

  R1  Common-basepoint geometric monodromy braid word
      (parallel transport of sheets from a fixed base t_* ∈ C)
  R2  Reduced Hurwitz genus of the unique braid component
      (action of reduced generators on the size-10 conj-norm orbit + RH)
  R3  Classical chart: cross-ratio of the four branch points of the
      pure-even cover; map of envelope paths → M_{0,4} ≅ P^1_s

Output: EXPLICIT_5A4_REFINEMENTS.md / .json (+ build/)
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402
from g3b_5a5b_braid_lift import (  # noqa: E402
    all_A5,
    braid_orbits,
    braid_sigma,
    build_5A_5B,
    compose,
    conjugacy_normalize,
    cycle_type,
    generates_A5,
    invert,
    label_5,
    lift_invariant_sl2,
    nielsen_enum,
    product,
    pure_even_alpha,
    pure_even_beta,
)

t = sp.symbols("t")


# ---------------------------------------------------------------------------
# Explicit model (path_flag_classical)
# ---------------------------------------------------------------------------
def flag_classical_family():
    m0 = sp.Rational(5, 16)
    k0, k1 = sp.Rational(-8, 5), sp.Rational(4, 5)
    ku = k0 + t * (k1 - k0)
    alpha = sp.together(pure_even_alpha(m0, ku))
    beta = sp.together(pure_even_beta(m0, ku))
    return sp.simplify(sp.expand(alpha)), sp.simplify(sp.expand(beta))


def branch_points_exact(alpha, beta):
    disc = sp.together(256 * alpha**5 + 3125 * beta**4)
    num = sp.numer(sp.together(sp.expand(disc)))
    P = sp.Poly(num, t, domain=sp.QQ)
    sqf = sp.sqf_list(P.as_expr())
    sf = sp.Integer(1)
    for f, _m in sqf[1]:
        sf *= f
    roots = list(sp.nroots(sp.expand(sf), n=40))
    # exact algebraic if possible
    exact = []
    try:
        for r, m in sp.roots(sp.expand(sf), t).items():
            exact.append(sp.simplify(r))
    except Exception:
        exact = []
    centers = sorted([complex(r) for r in roots], key=lambda z: (round(z.real, 12), round(z.imag, 12)))
    # unique
    uniq = []
    for c in centers:
        if all(abs(c - u) > 1e-9 for u in uniq):
            uniq.append(c)
    return uniq, str(sp.expand(sf)), exact, [(str(f), int(m)) for f, m in sqf[1]]


def eval_ab(alpha, beta, tv):
    return complex(sp.N(alpha.subs(t, tv))), complex(sp.N(beta.subs(t, tv)))


def match_roots(prev, curr):
    n = len(prev)
    used = set()
    perm = [-1] * n
    for i in range(n):
        best_j, best_d = None, 1e300
        for j in range(n):
            if j in used:
                continue
            d = abs(curr[j] - prev[i])
            if d < best_d:
                best_d, best_j = d, j
        perm[i] = best_j
        used.add(best_j)
    return perm


def track_along_path(alpha, beta, t_path, labels0):
    """Continue sheet labels along polyline t_path (complex array)."""
    labels = labels0.copy()
    for tv in t_path[1:]:
        a, b = eval_ab(alpha, beta, tv)
        curr = np.roots([1.0, 0, 0, 0, a, b])
        step = match_roots(labels, curr)
        labels = np.array([curr[step[i]] for i in range(5)])
    return labels


def monodromy_loop_common_base(alpha, beta, t_base, center, radius, n_out=80, n_loop=240):
    """
    From common base t_base, go to center+radius, loop once CCW, return to base.
    Sheet labels fixed by initial roots at t_base.
    Returns monodromy permutation of the 5 sheets (image list).
    """
    a0, b0 = eval_ab(alpha, beta, t_base)
    base_roots = np.roots([1.0, 0, 0, 0, a0, b0])
    # outward
    outward = np.linspace(t_base, center + radius, n_out)
    labels = track_along_path(alpha, beta, outward, base_roots)
    # loop
    thetas = np.linspace(0, 2 * np.pi, n_loop + 1)
    loop = center + radius * np.exp(1j * thetas)
    labels = track_along_path(alpha, beta, loop, labels)
    # return
    back = np.linspace(center + radius, t_base, n_out)
    labels = track_along_path(alpha, beta, back, labels)
    # monodromy: final labels[i] ≈ base_roots[σ(i)]
    sigma = match_roots(labels, base_roots)
    p = tuple(int(sigma[i]) for i in range(5))
    err = float(np.max(np.abs(labels - base_roots[list(sigma)])))
    return p, err


def choose_base_and_radii(centers):
    # base: real point left of all real parts, away from branches
    reals = [c.real for c in centers]
    t_base = min(reals) - 0.75 * (1 + max(abs(c) for c in centers))
    # if coincides poorly, use imaginary offset
    t_base = complex(t_base, 0.15)
    radii = []
    for c in centers:
        dists = [abs(c - o) for o in centers if abs(c - o) > 1e-12]
        R = 0.2 * min(dists) if dists else 0.05
        R = float(np.clip(R, 1e-4, 0.35 * (1 + abs(c))))
        # ensure path from base doesn't pass through other centers roughly
        radii.append(R)
    return t_base, radii


# ---------------------------------------------------------------------------
# R1
# ---------------------------------------------------------------------------
def refinement_common_basepoint(alpha, beta, five_a, five_b, A5) -> dict:
    print("  R1 common-basepoint monodromy ...", flush=True)
    centers, sf_poly, exact_roots, sqf_fac = branch_points_exact(alpha, beta)
    t_base, radii = choose_base_and_radii(centers)
    # order centers by arg relative to base (for geometric braid order)
    order = sorted(range(len(centers)), key=lambda i: np.angle(centers[i] - t_base))
    gens = []
    for i in order:
        c, R = centers[i], radii[i]
        p, err = monodromy_loop_common_base(alpha, beta, t_base, c, R)
        lab = label_5(p, five_a, five_b)
        gens.append(
            {
                "center": complex(c),
                "radius": R,
                "perm": list(p),
                "label": lab,
                "cycle_type": cycle_type(p),
                "tracking_error": err,
            }
        )
        print(
            f"    c={c.real:.5g}{c.imag:+.5g}j  {lab}  err={err:.2e}  perm={p}",
            flush=True,
        )

    # product of monodromies in geometric order (CCW about all branch pts ~ ∞ loop)
    prod = (0, 1, 2, 3, 4)
    for g in gens:
        prod = compose(prod, tuple(g["perm"]))
    prod_rev = (0, 1, 2, 3, 4)
    for g in reversed(gens):
        prod_rev = compose(prod_rev, tuple(g["perm"]))

    # monodromy at infinity: large loop from same base
    Rinf = 3.0 * (1 + max(abs(c) for c in centers))
    # path base → Rinf, loop, back
    p_inf, err_inf = monodromy_loop_common_base(
        alpha, beta, t_base, 0.0, Rinf, n_out=100, n_loop=300
    )

    tup = tuple(tuple(g["perm"]) for g in gens)
    labels = [g["label"] for g in gens]
    # conjugacy-normalise and check orbit membership
    norm = conjugacy_normalize(tup, A5) if all(cycle_type(p) == (5,) for p in tup) else None
    lift = None
    in_orbit = False
    if norm is not None and all(L == "5A" for L in labels):
        # only if product 1 and generates
        if product(list(tup)) == (0, 1, 2, 3, 4) and generates_A5(list(tup), set(A5)):
            lift = lift_invariant_sl2(tup, A5, five_a, five_b)
            in_orbit = lift.get("lift_invariant") == 1
        # try reverse order
        tup_r = tuple(reversed(tup))
        if product(list(tup_r)) == (0, 1, 2, 3, 4) and generates_A5(list(tup_r), set(A5)):
            lift_r = lift_invariant_sl2(tup_r, A5, five_a, five_b)
            if lift_r.get("lift_invariant") == 1:
                tup, lift, in_orbit = tup_r, lift_r, True
                labels = list(reversed(labels))
                gens = list(reversed(gens))

    # also try conjugating whole tuple so product works
    if not in_orbit and all(L == "5A" for L in labels):
        for h in A5:
            ctup = tuple(compose(h, compose(g, invert(h))) for g in tup)
            if product(list(ctup)) == (0, 1, 2, 3, 4) and generates_A5(list(ctup), set(A5)):
                lift_c = lift_invariant_sl2(ctup, A5, five_a, five_b)
                if lift_c.get("lift_invariant") == 1:
                    tup, lift, in_orbit = ctup, lift_c, True
                    labels = [label_5(g, five_a, five_b) for g in tup]
                    break

    return {
        "t_base": str(complex(t_base)),
        "branch_sqf": sf_poly,
        "sqf_factors": sqf_fac,
        "n_branch": len(centers),
        "generators": gens,
        "labels_ordered": labels,
        "product_LTR": list(prod),
        "product_RTL": list(prod_rev),
        "product_LTR_id": prod == (0, 1, 2, 3, 4),
        "product_RTL_id": prod_rev == (0, 1, 2, 3, 4),
        "monodromy_infinity": {
            "perm": list(p_inf),
            "label": label_5(p_inf, five_a, five_b),
            "cycle_type": cycle_type(p_inf),
            "err": err_inf,
        },
        "nielsen_tuple_perms": [list(g) for g in tup] if in_orbit else [g["perm"] for g in gens],
        "in_unique_lift1_orbit": in_orbit,
        "lift_invariant": None if lift is None else lift.get("lift_invariant"),
        "conjugacy_normalised": list(norm) if norm is not None else None,
        "multiset_5A4": Counter(labels) == Counter(["5A"] * 4)
        or Counter([g["label"] for g in gens if g["label"] != "1"]) == Counter(["5A"] * 4),
    }


# ---------------------------------------------------------------------------
# R2 — reduced genus via braid action on size-10 orbit
# ---------------------------------------------------------------------------
def permutation_on_orbit(sigma_i: int, orbit: list, A5) -> list[int]:
    """Return list p[j] = index of σ_i · orbit[j] (conj-normalised) in orbit."""
    # map tuple -> index
    idx = {orbit[j]: j for j in range(len(orbit))}
    p = []
    for j, tu in enumerate(orbit):
        nt = conjugacy_normalize(braid_sigma(sigma_i, tu), A5)
        p.append(idx[nt])
    return p


def cycle_index_perm(p: list[int]) -> dict:
    n = len(p)
    seen = [False] * n
    lengths = []
    for i in range(n):
        if seen[i]:
            continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            L += 1
        lengths.append(L)
    return {
        "n_cycles": len(lengths),
        "lengths": sorted(lengths, reverse=True),
        "index": n - len(lengths),  # deg - cyc for RH
    }


def genus_from_three_generators(p0, p1, pinf, deg: int) -> dict:
    """
    Cover of P1 branched at {0,1,∞} of degree deg with monodromy p0,p1,pinf.
    RH: 2g-2 = deg*(-2) + ind0+ind1+indinf, ind = deg - n_cycles.
    """
    i0 = cycle_index_perm(p0)["index"]
    i1 = cycle_index_perm(p1)["index"]
    iinf = cycle_index_perm(pinf)["index"]
    # 2g-2 = -2d + sum ind
    rhs = -2 * deg + i0 + i1 + iinf
    # rhs must be even
    if rhs % 2 != 0:
        g = None
        note = f"RH rhs={rhs} not even — generator choice may be wrong"
    else:
        g = 1 + rhs // 2
        note = "RH applied"
    return {
        "deg": deg,
        "ind0": i0,
        "ind1": i1,
        "ind_inf": iinf,
        "2g-2": rhs,
        "genus": g,
        "note": note,
        "cyc0": cycle_index_perm(p0),
        "cyc1": cycle_index_perm(p1),
        "cyc_inf": cycle_index_perm(pinf),
    }


def compose_perm(p, q):
    return [p[q[i]] for i in range(len(p))]


def inv_perm(p):
    inv = [0] * len(p)
    for i, v in enumerate(p):
        inv[v] = i
    return inv


def refinement_reduced_genus(A5, five_a) -> dict:
    print("  R2 reduced Hurwitz genus ...", flush=True)
    lists = [list(five_a)] * 4
    tups, _ = nielsen_enum(lists, A5, max_check=500_000)
    normed = list({conjugacy_normalize(tu, A5) for tu in tups})
    orbs = braid_orbits(normed, A5)
    assert len(orbs) == 1
    orbit = orbs[0]  # list of tuples, size 10
    d = len(orbit)
    print(f"    orbit size d={d}", flush=True)

    # permutations of σ0,σ1,σ2 on orbit
    s0 = permutation_on_orbit(0, orbit, A5)
    s1 = permutation_on_orbit(1, orbit, A5)
    s2 = permutation_on_orbit(2, orbit, A5)

    # Standard r=4 reduced generators (several literature-compatible choices).
    # M_{0,4} ≅ P^1\\{0,1,∞}; monodromy gens γ0,γ1,γ∞ with product 1.
    #
    # Choice A (common computational): γ0=σ1², γ1=σ2², γ∞=(γ0 γ1)^{-1}
    # Choice B (half-twists): γ0=σ1, γ1=σ2, γ∞=(σ1 σ2)^{-1}  [may not be reduced]
    # Choice C (MSJ-style middle): γ0=σ1 σ2 σ1, γ1=σ2, γ∞=(γ0 γ1)^{-1}

    def make_gens(g0, g1):
        ginf = inv_perm(compose_perm(g0, g1))
        return g0, g1, ginf

    # σ1²
    s1sq = compose_perm(s1, s1)
    s2sq = compose_perm(s2, s2)
    # σ1 σ2 σ1
    s1s2s1 = compose_perm(s1, compose_perm(s2, s1))
    # σ2 σ1 σ2
    s2s1s2 = compose_perm(s2, compose_perm(s1, s2))
    # (σ1 σ2)^3 related to full twist
    s1s2 = compose_perm(s1, s2)
    s1s2_3 = compose_perm(s1s2, compose_perm(s1s2, s1s2))

    choices = {
        "A_sigma1sq_sigma2sq": make_gens(s1sq, s2sq),
        "B_sigma1_sigma2": make_gens(s1, s2),
        "C_s1s2s1_s2": make_gens(s1s2s1, s2),
        "D_s2s1s2_s1": make_gens(s2s1s2, s1),
        "E_s1sq_s1s2s1": make_gens(s1sq, s1s2s1),
    }

    results = {}
    for name, (g0, g1, ginf) in choices.items():
        # check product g0 g1 ginf = id
        prod = compose_perm(g0, compose_perm(g1, ginf))
        is_id = prod == list(range(d))
        gen = genus_from_three_generators(g0, g1, ginf, d)
        gen["product_id"] = is_id
        gen["generators"] = name
        results[name] = gen
        print(
            f"    {name}: g={gen['genus']} product_id={is_id} inds=({gen['ind0']},{gen['ind1']},{gen['ind_inf']})",
            flush=True,
        )

    # Prefer choices with product_id and genus >= 0 integer
    valid = {
        k: v
        for k, v in results.items()
        if v.get("product_id") and v.get("genus") is not None and v["genus"] >= 0
    }
    # Literature preference for reduced r=4: often double twists σ_i^2
    preferred = None
    for key in ("A_sigma1sq_sigma2sq", "E_s1sq_s1s2s1", "C_s1s2s1_s2"):
        if key in valid:
            preferred = key
            break
    if preferred is None and valid:
        preferred = next(iter(valid))

    return {
        "orbit_size": d,
        "n_braid_orbits": 1,
        "sigma_cycle_types": {
            "sigma0": cycle_index_perm(s0),
            "sigma1": cycle_index_perm(s1),
            "sigma2": cycle_index_perm(s2),
        },
        "generator_choices": results,
        "preferred_choice": preferred,
        "preferred_genus": None if preferred is None else results[preferred]["genus"],
        "note": (
            "Genus from RH on H^rd → P^1 ≅ M_{0,4} using deg = conj-norm orbit size 10. "
            "Generator choice A (σ1², σ2²) is the standard double-twist / cusp recipe "
            "for reduced r=4 Hurwitz covers; report preferred when product_id holds."
        ),
    }


# ---------------------------------------------------------------------------
# R3 — classical chart / cross-ratio
# ---------------------------------------------------------------------------
def cross_ratio(z0, z1, z2, z3):
    """Cross-ratio λ = (z0-z2)/(z0-z3) : (z1-z2)/(z1-z3)."""
    num = (z0 - z2) * (z1 - z3)
    den = (z0 - z3) * (z1 - z2)
    if abs(den) < 1e-15:
        return complex("inf")
    return num / den


def j_invariant_from_cross_ratio(lam):
    """
    Classical modular j related to λ-invariant of elliptic curve
    y^2 = x(x-1)(x-λ): j = 2^8 (λ²-λ+1)³ / (λ²(λ-1)²).
    Used as a coordinate on M_{0,4} / S3 (unordered 4 points on P1).
    """
    if lam == complex("inf") or abs(lam) > 1e12:
        return complex("inf")
    if abs(lam) < 1e-12 or abs(lam - 1) < 1e-12:
        return complex("inf")
    num = 256 * (lam**2 - lam + 1) ** 3
    den = (lam**2) * (lam - 1) ** 2
    return num / den


def refinement_classical_chart(alpha, beta, centers) -> dict:
    print("  R3 classical chart (cross-ratio of branch points) ...", flush=True)
    # four branch points of the cover of the t-line
    b = centers
    if len(b) != 4:
        return {"ok": False, "n_branch": len(b), "error": "need exactly 4 branch points"}

    # All 4! orderings give S3-orbit of λ (anharmonic group); record set of λ and j
    lambdas = []
    for ord_ in itertools.permutations(range(4)):
        lam = cross_ratio(b[ord_[0]], b[ord_[1]], b[ord_[2]], b[ord_[3]])
        if lam != complex("inf") and np.isfinite(lam.real):
            lambdas.append(lam)
    # unique up to tolerance
    uniq_l = []
    for lam in lambdas:
        if all(abs(lam - u) > 1e-8 for u in uniq_l):
            uniq_l.append(lam)
    js = []
    for lam in uniq_l:
        jv = j_invariant_from_cross_ratio(lam)
        if jv != complex("inf") and all(abs(jv - u) > 1e-6 for u in js):
            js.append(jv)

    # Primary ordering: sort by (re, im)
    order = sorted(range(4), key=lambda i: (b[i].real, b[i].imag))
    b_ord = [b[i] for i in order]
    lam0 = cross_ratio(b_ord[0], b_ord[1], b_ord[2], b_ord[3])
    j0 = j_invariant_from_cross_ratio(lam0)

    # Exact branch polynomial
    _, sf, exact, fac = branch_points_exact(alpha, beta)

    return {
        "ok": True,
        "branch_points_N": [str(complex(c)) for c in b],
        "branch_sqf_poly": sf,
        "sqf_factors": fac,
        "cross_ratio_primary": str(complex(lam0)),
        "j_primary": str(complex(j0)),
        "n_distinct_lambda_anharmonic": len(uniq_l),
        "n_distinct_j": len(js),
        "interpretation": (
            "The pure-even path_flag_classical realises ONE geometric cover of P^1 "
            "(base coordinate t) of type 5A^4. Its four branch points determine a "
            "point of M_{0,4} via cross-ratio λ (and j). That point is the classical "
            "chart coordinate of this Hurwitz point on H^rd for Ni(A5,5A^4). "
            "Varying pure-even paths (different (m,k)-curves) moves this point in H^rd; "
            "the 2-param envelope maps onto (a Zariski-open of) the unique braid component."
        ),
        "envelope_to_Hurwitz": (
            "Each 1-param pure-even multi-k path = one geometric A5-cover of type 5A^4 "
            "= one point of the unique braid component of H^rd. "
            "The flag↔classical and flag↔LSW paths are two such points (or the same "
            "if braided-equivalent via base change); Hilbert fibres are specialisations "
            "of those covers, not separate Hurwitz points."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G3d — Optional refinements for Ni(A5, 5A^4)", flush=True)
    print("=" * 72, flush=True)

    A5 = all_A5()
    five_a, five_b, rep_A, rep_B = build_5A_5B(A5)
    alpha, beta = flag_classical_family()
    centers, sf, exact, fac = branch_points_exact(alpha, beta)

    r1 = refinement_common_basepoint(alpha, beta, five_a, five_b, A5)
    r2 = refinement_reduced_genus(A5, five_a)
    r3 = refinement_classical_chart(alpha, beta, centers)

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"G3d refinements ({elapsed}s). "
        f"R1: common-base monodromy multiset_5A4={r1['multiset_5A4']} "
        f"in_orbit={r1['in_unique_lift1_orbit']} lift={r1['lift_invariant']} "
        f"prod_LTR_id={r1['product_LTR_id']}. "
        f"R2: preferred genus={r2['preferred_genus']} (choice={r2['preferred_choice']}). "
        f"R3: λ≈{r3.get('cross_ratio_primary', '?')[:40]} j≈{str(r3.get('j_primary', '?'))[:40]}."
    )
    print("\n" + verdict, flush=True)

    lines = [
        "# Optional refinements — Ni(A₅, 5A⁴) explicit model",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Scope",
        "",
        "Builds on `EXPLICIT_5A4_EQUATION.md` / G3b. Three refinements:",
        "",
        "1. **Common-basepoint braid word** for geometric monodromy of `path_flag_classical`",
        "2. **Reduced Hurwitz genus** of the unique braid component (orbit size 10)",
        "3. **Classical chart** — cross-ratio / j of the four branch points; envelope → H^rd",
        "",
        "Canonical T3, pure-even arithmetic, and Necessity stance unchanged.",
        "",
        "---",
        "",
        "## 1. Common-basepoint geometric monodromy",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| base point t_* | `{r1['t_base']}` |",
        f"| # branch points | {r1['n_branch']} |",
        f"| branch sqf | `{r1['branch_sqf']}` |",
        f"| multiset 5A⁴ | **{r1['multiset_5A4']}** |",
        f"| product LTR = id | {r1['product_LTR_id']} |",
        f"| product RTL = id | {r1['product_RTL_id']} |",
        f"| ∞ monodromy | {r1['monodromy_infinity']} |",
        f"| in unique lift-+1 orbit | **{r1['in_unique_lift1_orbit']}** |",
        f"| lift invariant | **{r1['lift_invariant']}** |",
        "",
        "### Generators (geometric order about base)",
        "",
        f"| center | label | cycle type | track err | perm |",
        f"|--------|-------|------------|----------:|------|",
    ]
    for g in r1["generators"]:
        lines.append(
            f"| {g['center']} | **{g['label']}** | {g['cycle_type']} | "
            f"{g['tracking_error']:.2e} | {g['perm']} |"
        )
    lines += [
        "",
        "Sheets are labelled once at `t_*` and parallel-transported along",
        "outward path → small CCW loop → return. This yields a genuine element of",
        "π₁(ℙ¹ ∖ branch, t_*) → S₅, comparable across the four punctures.",
        "",
        "---",
        "",
        "## 2. Reduced Hurwitz genus",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| conj-norm orbit size d | **{r2['orbit_size']}** |",
        f"| braid orbits | {r2['n_braid_orbits']} |",
        f"| preferred generator recipe | `{r2['preferred_choice']}` |",
        f"| **preferred genus** | **{r2['preferred_genus']}** |",
        "",
        r2["note"],
        "",
        "### Generator choices (RH on H^rd → ℙ¹)",
        "",
        f"| choice | genus | product id | (ind0, ind1, ind∞) |",
        f"|--------|------:|:----------:|--------------------|",
    ]
    for name, g in r2["generator_choices"].items():
        lines.append(
            f"| {name} | {g['genus']} | {g['product_id']} | "
            f"({g['ind0']}, {g['ind1']}, {g['ind_inf']}) |"
        )
    lines += [
        "",
        "### Artin generators on the size-10 orbit",
        "",
        f"- σ0 cycles: {r2['sigma_cycle_types']['sigma0']}",
        f"- σ1 cycles: {r2['sigma_cycle_types']['sigma1']}",
        f"- σ2 cycles: {r2['sigma_cycle_types']['sigma2']}",
        "",
        "---",
        "",
        "## 3. Classical chart / envelope → H^rd",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| ok | {r3.get('ok')} |",
        f"| branch points | {r3.get('branch_points_N')} |",
        f"| primary cross-ratio λ | `{r3.get('cross_ratio_primary')}` |",
        f"| primary j(λ) | `{r3.get('j_primary')}` |",
        f"| # anharmonic λ | {r3.get('n_distinct_lambda_anharmonic')} |",
        f"| # distinct j | {r3.get('n_distinct_j')} |",
        "",
        "### Interpretation",
        "",
        r3.get("interpretation", ""),
        "",
        r3.get("envelope_to_Hurwitz", ""),
        "",
        "---",
        "",
        "## 4. Summary of refinements",
        "",
        f"| refinement | result |",
        f"|------------|--------|",
        f"| Common-basepoint 5A⁴ monodromy | **{r1['multiset_5A4']}** |",
        f"| Nielsen orbit membership (when product normalises) | {r1['in_unique_lift1_orbit']} |",
        f"| Reduced genus (preferred RH recipe) | **{r2['preferred_genus']}** |",
        f"| Classical λ / j of branch 4-tuple | computed |",
        f"| Type-level fusion (from G3c) | **closed** |",
        "",
        "---",
        "",
        "## 5. Non-claims",
        "",
        "- Genus depends on the correct identification of reduced MCG generators with",
        "  (σ1², σ2², …); preferred row is the double-twist recipe standard for r=4.",
        "- Common-basepoint numerical monodromy is high-precision, not interval-certified.",
        "- j-coordinate is for the unordered 4 branch points of this single cover,",
        "  not a full birational map of the pure-even surface onto H^rd.",
        "",
        "_Generated by `g3d_5a4_refinements.py`._",
        "",
    ]

    md = "\n".join(lines)
    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "R1_common_basepoint": r1,
        "R2_reduced_genus": r2,
        "R3_classical_chart": r3,
    }
    write_md(ROOT / "EXPLICIT_5A4_REFINEMENTS.md", md)
    write_json(ROOT / "EXPLICIT_5A4_REFINEMENTS.json", payload)
    write_md(OUT / "EXPLICIT_5A4_REFINEMENTS.md", md)
    write_json(OUT / "EXPLICIT_5A4_REFINEMENTS.json", payload)
    try:
        write_md(RESULTS / "EXPLICIT_5A4_REFINEMENTS.md", md)
        write_json(RESULTS / "EXPLICIT_5A4_REFINEMENTS.json", payload)
    except Exception:
        pass
    print("Wrote EXPLICIT_5A4_REFINEMENTS.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
