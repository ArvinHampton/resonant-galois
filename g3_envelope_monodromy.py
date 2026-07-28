#!/usr/bin/env python3
"""
G3 — monodromy identification of the pure-even envelope.

Arithmetic multi-k is finished: the pure-even envelope
  α = 256 m² − 3125 k⁴/256,  β = k α
has disc ≡ (256 α² m)² identically, and cross-k paths hit multiple catalogue
ratios. G3 asks: what is the *geometric* monodromy of these 1-parameter
families as covers of the t-line — can it be named as a Nielsen class
Ni(A5, C)?

Pipeline
--------
E0  Define envelope paths (flagship↔classical, flagship↔LSW, fixed-k rays, …)
E1  Symbolic disc D(t) = disc_x(f_t); square-free branch locus; factorisation
E2  Numerical local monodromy: loop t about each finite branch point + ∞;
    track 5 roots → cycle type in S5
E3  Assemble conjugacy multiset; match against r=4 A5 Nielsen shortlist
E4  Group generation / A5 check from sampled monodromy permutations
E5  Catalogue multi-k confirmation on the same paths (arithmetic control)
E6  2-parameter envelope structural notes (branch divisor in (m,k))

Output: G3_ENVELOPE_MONODROMY.md / .json (+ build/)
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

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj, disc_bj_int  # noqa: E402

t = sp.symbols("t")

# ---------------------------------------------------------------------------
# Catalogue / pure-even
# ---------------------------------------------------------------------------
CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("classical_m", 20, -16, Fraction(-4, 5)),
]
CAT_BY_AB = {(a, b): (tag, k) for tag, a, b, k in CATALOGUE}


def pure_even_alpha(m, k):
    return 256 * m**2 - sp.Rational(3125) * k**4 / 256


def pure_even_beta(m, k):
    return k * pure_even_alpha(m, k)


# ---------------------------------------------------------------------------
# E0 — families
# ---------------------------------------------------------------------------
def family_flag_classical():
    """Fixed m=5/16, k linear −8/5 → 4/5. Endpoints: flagship / classical scale."""
    m0 = sp.Rational(5, 16)
    k0, k1 = sp.Rational(-8, 5), sp.Rational(4, 5)
    ku = k0 + t * (k1 - k0)
    alpha = sp.together(pure_even_alpha(m0, ku))
    beta = sp.together(pure_even_beta(m0, ku))
    return {
        "id": "path_flag_classical",
        "m": str(m0),
        "k(t)": str(ku),
        "alpha": alpha,
        "beta": beta,
        "note": "same-m linear-k; arithmetic multi-k flagship↔classical",
        "expected_hits": [("flagship", 0), ("classical", 1)],
    }


def family_flag_lsw():
    m0 = sp.Rational(5, 16)
    # LSW sample (-100,400) has k=-4; m from pure-even: α = 256m² - 3125*256/256 = 256m²-3125
    # -100 = 256m² - 3125 ⇒ 256m² = 3025 ⇒ m² = 3025/256 ⇒ m = 55/16
    m1 = sp.Rational(55, 16)
    # path: (m,k) linear from (5/16,-8/5) to (55/16,-4)
    mu = m0 + t * (m1 - m0)
    ku = sp.Rational(-8, 5) + t * (sp.Rational(-4) - sp.Rational(-8, 5))
    alpha = sp.together(pure_even_alpha(mu, ku))
    beta = sp.together(pure_even_beta(mu, ku))
    return {
        "id": "path_flag_lsw",
        "m(t)": str(mu),
        "k(t)": str(ku),
        "alpha": alpha,
        "beta": beta,
        "note": "linear (m,k) flagship→LSW",
        "expected_hits": [("flagship", 0), ("lsw_m100", 1)],
    }


def family_lsw_ray():
    """Fixed k=-4, m = t (cleared LSW ray α=256t²-3125)."""
    ku = sp.Integer(-4)
    mu = t
    alpha = sp.together(pure_even_alpha(mu, ku))
    beta = sp.together(pure_even_beta(mu, ku))
    return {
        "id": "ray_lsw_k_m4",
        "k": "-4",
        "m(t)": "t",
        "alpha": alpha,
        "beta": beta,
        "note": "fixed-k pure-even ray (single k)",
        "expected_hits": [],
    }


def family_flagship_homog():
    """Homogenisation of flagship seed: α=-55 t^4, β=88 t^5."""
    alpha = -55 * t**4
    beta = 88 * t**5
    return {
        "id": "homog_flagship",
        "alpha": alpha,
        "beta": beta,
        "note": "homogenised flagship ray (Theorem 3); single seed k=-8/5",
        "expected_hits": [("flagship", 1)],  # t=1
    }


def family_classical_homog():
    alpha = 20 * t**4
    beta = 16 * t**5
    return {
        "id": "homog_classical",
        "alpha": alpha,
        "beta": beta,
        "note": "homogenised classical; k=4/5",
        "expected_hits": [("classical", 1)],
    }


# ---------------------------------------------------------------------------
# E1 — discriminant / branch locus
# ---------------------------------------------------------------------------
def disc_family(alpha, beta):
    """D(t) = 256 α^5 + 3125 β^4 as sympy expr in t."""
    return sp.together(256 * alpha**5 + 3125 * beta**4)


def analyze_branch_locus(alpha, beta) -> dict:
    D = disc_family(alpha, beta)
    D_exp = sp.expand(sp.together(D))
    # content / polynomial in t
    try:
        P = sp.Poly(sp.numer(sp.together(D_exp)), t, domain=sp.QQ)
    except Exception:
        # may be Laurent if denoms; clear
        D2 = sp.together(D_exp)
        num, den = sp.fraction(D2)
        P = sp.Poly(sp.expand(num), t, domain=sp.QQ)
    # square-free part
    try:
        sqf = sp.sqf_list(P.as_expr())
        # sqf = (content, [(factor, mult), ...])
        factors = [(str(f), int(m)) for f, m in sqf[1]]
        square_free = sp.Integer(1)
        for f, m in sqf[1]:
            square_free *= f
        sf_poly = sp.Poly(sp.expand(square_free), t, domain=sp.QQ)
    except Exception as e:
        factors = []
        sf_poly = P
        sqf = (1, [])

    # all multiplicities even? (consistent with pure-even identical square)
    mults = [m for _f, m in factors]
    all_even = all(m % 2 == 0 for m in mults) if mults else None

    # finite branch points: roots of square-free part
    roots = []
    try:
        for r, m in sp.roots(sf_poly.as_expr(), t).items():
            roots.append({"root": str(sp.simplify(r)), "root_N": complex(sp.N(r)), "mult_in_sf": int(m)})
    except Exception:
        # numerical
        try:
            for r in sp.nroots(sf_poly.as_expr(), n=20):
                roots.append({"root": str(r), "root_N": complex(r), "mult_in_sf": 1})
        except Exception:
            pass

    # Also zeros of α (often the pure-even branch source)
    try:
        Pa = sp.Poly(sp.numer(sp.together(sp.expand(alpha))), t, domain=sp.QQ)
        alpha_roots = [complex(sp.N(r)) for r in sp.nroots(Pa.as_expr(), n=15)]
    except Exception:
        alpha_roots = []

    return {
        "disc_degree": int(P.degree()) if P else None,
        "disc_leading": str(P.LC()) if P else None,
        "sqf_factors": factors,
        "all_mult_even": all_even,
        "n_finite_branch_sf": len(roots),
        "branch_roots": roots[:20],
        "alpha_zeros_N": alpha_roots[:12],
        "disc_expr_head": str(D_exp)[:200],
    }


# ---------------------------------------------------------------------------
# E2 — numerical monodromy
# ---------------------------------------------------------------------------
def roots_at(alpha_c, beta_c):
    """Roots of x^5 + α x + β."""
    # np.roots of [1,0,0,0,α,β]
    return np.roots([1.0, 0.0, 0.0, 0.0, float(np.real(alpha_c)), float(np.real(beta_c))])


def eval_ab(alpha_expr, beta_expr, tval):
    a = complex(sp.N(alpha_expr.subs(t, tval)))
    b = complex(sp.N(beta_expr.subs(t, tval)))
    return a, b


def match_roots(prev, curr):
    """Greedy nearest-neighbour matching prev→curr; return permutation of indices."""
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


def perm_compose(p, q):
    """(p∘q)(i) = p(q(i)) as lists of images."""
    return [p[q[i]] for i in range(len(p))]


def perm_to_cycles(perm):
    n = len(perm)
    seen = [False] * n
    cycles = []
    for i in range(n):
        if seen[i]:
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = perm[j]
        cycles.append(tuple(cyc))
    return cycles


def cycle_type_from_perm(perm):
    cycles = perm_to_cycles(perm)
    lengths = sorted([len(c) for c in cycles], reverse=True)
    return tuple(lengths)


def class_name_A5(ct: tuple) -> str:
    """Map cycle type to A5/S5 conjugacy label (as partitions of 5)."""
    # normalize padding
    parts = tuple(sorted([c for c in ct if c > 0], reverse=True))
    if parts == (5,):
        return "5A/5B"  # need lift to distinguish; report 5*
    if parts == (3, 2):
        return "odd:(3,2)"  # not in A5
    if parts == (3, 1, 1):
        return "3A"
    if parts == (2, 2, 1):
        return "2A"
    if parts == (2, 1, 1, 1):
        return "odd:transp"
    if parts == (4, 1):
        return "odd:4"
    if parts == (1, 1, 1, 1, 1):
        return "1"
    return f"type{parts}"


def local_monodromy(alpha_expr, beta_expr, center, radius=None, nsteps=240):
    """
    Walk t = center + R e^{iθ}, θ:0→2π; track roots; return cycle type + perm.
    """
    c = complex(center)
    if radius is None:
        radius = max(1e-3, 0.05 * (1.0 + abs(c)))
    # avoid radius hitting other branches roughly — user can refine
    thetas = np.linspace(0, 2 * np.pi, nsteps + 1)
    t0 = c + radius
    a0, b0 = eval_ab(alpha_expr, beta_expr, t0)
    # if coeffs huge, scale issue — still ok
    prev = roots_at(a0.real if abs(a0.imag) < 1e-10 else a0, b0.real if abs(b0.imag) < 1e-10 else b0)
    # Use complex coeffs properly
    prev = np.roots([1.0, 0, 0, 0, a0, b0])
    total_perm = list(range(5))  # identity on label space = initial order

    # We track: position list ordered by continuing labels 0..4
    labels = prev.copy()  # labels[i] = current complex root for label i

    for th in thetas[1:]:
        tv = c + radius * np.exp(1j * th)
        a, b = eval_ab(alpha_expr, beta_expr, tv)
        curr = np.roots([1.0, 0, 0, 0, a, b])
        # match labels[i] to curr
        perm_step = match_roots(labels, curr)
        new_labels = np.array([curr[perm_step[i]] for i in range(5)])
        # permutation of sheet indices: which old label went where
        # After step, label i is at curr[perm_step[i]]; reindex so
        # the monodromy sigma satisfies labels_after = curr reordered...
        # Standard: sigma(i) = j means root i continues to slot j in ordered curr after matching.
        # Accumulate: compose inverse matching
        # Here: root that was at labels[i] moves to new_labels[i] = curr[perm_step[i]]
        # The monodromy permutation on the set of sheets (initial root indices):
        # we keep fixed labels 0..4 attached to continuously tracked roots.
        # So the permutation of *positions in a fixed ordering of curr* is not needed —
        # at end of loop, compare final labels to start via matching.
        labels = new_labels

    # Final monodromy: match final labels to initial prev roots
    final_perm_match = match_roots(prev, labels)
    # final_perm_match[i] = j means initial root i is now at labels position that
    # matches prev[j]? match_roots(prev, labels)[i] = index in labels closest to prev[i]
    # After full loop, labels[k] is the continuation of initial labels[k]=prev[k].
    # So labels should be a permutation of prev: labels[i] ≈ prev[σ(i)]
    # Find σ such that labels[i] ≈ prev[σ(i)]
    sigma = match_roots(labels, prev)  # sigma[i] = j with prev[j]≈labels[i]
    # Wait: match_roots(A,B)[i] = j means B[j] closest to A[i]
    # We want σ with labels[i] ≈ prev[σ(i)], so σ(i) = match_roots(labels, prev)[i]
    sigma = match_roots(labels, prev)
    ct = cycle_type_from_perm(sigma)
    return {
        "center": str(center),
        "center_N": complex(c),
        "radius": radius,
        "cycle_type": ct,
        "class": class_name_A5(ct),
        "perm": sigma,
        "tracking_error": float(np.max([abs(labels[i] - prev[sigma[i]]) for i in range(5)])),
    }


def monodromy_at_infinity(alpha_expr, beta_expr, R=1e3, nsteps=300):
    return local_monodromy(alpha_expr, beta_expr, 0.0, radius=R, nsteps=nsteps)


# ---------------------------------------------------------------------------
# E3 — Nielsen matching
# ---------------------------------------------------------------------------
# Shortlist of r=4 types (from A5_HURWITZ_R4)
NIELSEN_SHORTLIST = {
    "3A,3A,3A,3A": {"classes": ["3A", "3A", "3A", "3A"], "orbit": 18, "g": 0},
    "2A,3A,3A,3A": {"classes": ["2A", "3A", "3A", "3A"], "orbit": 96, "g": 0},
    "2A,2A,3A,3A": {"classes": ["2A", "2A", "3A", "3A"], "orbit": 108, "g": 0},
    "2A,3A,3A,5*": {"classes": ["2A", "3A", "3A", "5A/5B"], "orbit": 240, "g": 0},
    "3A,3A,3A,5*": {"classes": ["3A", "3A", "3A", "5A/5B"], "orbit": "40-60", "g": 0},
    "2A,2A,3A,5*": {"classes": ["2A", "2A", "3A", "5A/5B"], "orbit": 180, "g": 0},
}


def match_nielsen(class_list: list[str]) -> list[dict]:
    """Compare multiset of local monodromy classes to shortlist (order-free)."""
    # filter out identity and failed
    cl = [c for c in class_list if c not in ("1",) and not c.startswith("type")]
    # map 5A/5B to 5*
    norm = []
    for c in cl:
        if c.startswith("5"):
            norm.append("5A/5B")
        elif c.startswith("odd"):
            norm.append(c)
        else:
            norm.append(c)
    cnt = Counter(norm)
    hits = []
    for name, info in NIELSEN_SHORTLIST.items():
        target = Counter(info["classes"])
        # allow 5* wildcards
        if cnt == target:
            hits.append({"nielsen": name, "match": "exact", **info})
        else:
            # partial: same support size
            if sum(cnt.values()) == sum(target.values()) and set(cnt) <= set(target) | {"5A/5B"}:
                # soft
                if sorted(cnt.elements()) == sorted(target.elements()):
                    hits.append({"nielsen": name, "match": "exact", **info})
    # also report sorted signature
    signature = "+".join(f"{k}×{v}" if v > 1 else k for k, v in sorted(cnt.items()))
    return {"signature": signature, "counter": dict(cnt), "matches": hits, "raw_classes": class_list}


# ---------------------------------------------------------------------------
# E4 — generate subgroup from perms
# ---------------------------------------------------------------------------
def generates_A5_from_perms(perms: list[list[int]]) -> dict:
    """BFS generate subgroup of S5; check order and alternateness."""
    if not perms:
        return {"order": 1, "is_A5": False, "is_S5": False}

    def compose(p, q):
        return [p[q[i]] for i in range(5)]

    def invert(p):
        inv = [0] * 5
        for i, v in enumerate(p):
            inv[v] = i
        return inv

    def sign(p):
        # sign of permutation
        seen = [False] * 5
        sig = 1
        for i in range(5):
            if seen[i]:
                continue
            # cycle length
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:  # even length cycle is odd perm
                sig *= -1
        return sig

    idp = list(range(5))
    seen = {tuple(idp)}
    queue = [idp]
    gens = list(perms) + [invert(p) for p in perms]
    while queue:
        g = queue.pop()
        for s in gens:
            h = compose(g, s)
            th = tuple(h)
            if th not in seen:
                seen.add(th)
                queue.append(h)
                if len(seen) > 120:
                    break
        if len(seen) > 120:
            break
    order = len(seen)
    all_even = all(sign(list(p)) == 1 for p in seen)
    return {
        "order": order,
        "is_A5": order == 60 and all_even,
        "is_S5": order == 120,
        "all_even": all_even,
        "n_generators": len(perms),
    }


# ---------------------------------------------------------------------------
# E5 — catalogue hits on path
# ---------------------------------------------------------------------------
def catalogue_hits_on_family(alpha, beta, t_vals) -> list[dict]:
    hits = []
    for tv in t_vals:
        try:
            a = sp.simplify(alpha.subs(t, tv))
            b = sp.simplify(beta.subs(t, tv))
            aR, bR = sp.Rational(a), sp.Rational(b)
            if aR.denominator != 1 or bR.denominator != 1:
                # try clear common — skip non-Z BJ
                continue
            aa, bb = int(aR), int(bR)
        except Exception:
            continue
        if (aa, bb) in CAT_BY_AB:
            tag, k = CAT_BY_AB[(aa, bb)]
            hits.append({"tag": tag, "k": str(k), "t": str(tv), "alpha": aa, "beta": bb})
    return hits


# ---------------------------------------------------------------------------
# Analyze one family end-to-end
# ---------------------------------------------------------------------------
def analyze_family(fam: dict) -> dict:
    print(f"\n=== {fam['id']} ===", flush=True)
    alpha, beta = fam["alpha"], fam["beta"]

    # E1
    print("  branch locus ...", flush=True)
    locus = analyze_branch_locus(alpha, beta)
    print(
        f"    deg D={locus['disc_degree']} factors={locus['sqf_factors']} "
        f"all_even_mult={locus['all_mult_even']} n_branch≈{locus['n_finite_branch_sf']}",
        flush=True,
    )

    # Collect branch centers for monodromy: finite roots + use alpha zeros
    centers = []
    for br in locus["branch_roots"]:
        centers.append(br["root_N"])
    for az in locus["alpha_zeros_N"]:
        # add if not duplicate
        if all(abs(az - c) > 1e-6 for c in centers):
            centers.append(az)

    # Cap number of centers for cost
    # Prefer finite small ones
    centers = sorted(centers, key=lambda z: abs(z))[:12]

    # E2 local monodromy
    print(f"  local monodromy at {len(centers)} centers + infinity ...", flush=True)
    locals_m = []
    for c in centers:
        if abs(c) > 1e6:
            continue
        # radius: half distance to nearest other center
        dists = [abs(c - o) for o in centers if abs(c - o) > 1e-12]
        R = 0.3 * min(dists) if dists else 0.05 * (1 + abs(c))
        R = max(R, 1e-4)
        R = min(R, 0.5 * (1 + abs(c)))
        try:
            lm = local_monodromy(alpha, beta, c, radius=R, nsteps=200)
            locals_m.append(lm)
            print(
                f"    c≈{c.real:.4g}{c.imag:+.4g}j  type={lm['cycle_type']} class={lm['class']} "
                f"err={lm['tracking_error']:.2e}",
                flush=True,
            )
        except Exception as e:
            locals_m.append({"center_N": complex(c), "error": str(e), "class": "fail"})

    try:
        minf = monodromy_at_infinity(alpha, beta, R=500.0, nsteps=280)
        print(f"    ∞  type={minf['cycle_type']} class={minf['class']}", flush=True)
    except Exception as e:
        minf = {"error": str(e), "class": "fail"}

    classes = [lm.get("class", "fail") for lm in locals_m if lm.get("class") not in ("fail", "1")]
    if minf.get("class") not in ("fail", "1", None):
        classes_with_inf = classes + [minf["class"]]
    else:
        classes_with_inf = classes

    # E3 match
    match_fin = match_nielsen(classes)
    match_all = match_nielsen(classes_with_inf)

    # E4 group from finite perms with low tracking error
    perms = [
        lm["perm"]
        for lm in locals_m
        if "perm" in lm and lm.get("tracking_error", 1) < 1e-4
    ]
    if "perm" in minf and minf.get("tracking_error", 1) < 1e-3:
        perms_inf = perms + [minf["perm"]]
    else:
        perms_inf = perms
    gen = generates_A5_from_perms(perms)
    gen_inf = generates_A5_from_perms(perms_inf)

    # E5 catalogue
    t_vals = [0, 1, sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(1, 4)]
    t_vals += list(range(-3, 4))
    hits = catalogue_hits_on_family(alpha, beta, t_vals)
    cat_k = sorted({h["k"] for h in hits})

    # Sample Gal at a regular t
    sample_gal = None
    for tv in [sp.Rational(1, 2), sp.Integer(2), sp.Rational(3, 2)]:
        try:
            a = sp.Rational(sp.simplify(alpha.subs(t, tv)))
            b = sp.Rational(sp.simplify(beta.subs(t, tv)))
            if a.denominator != 1 or b.denominator != 1:
                continue
            aa, bb = int(a), int(b)
            if aa == 0:
                continue
            d = disc_bj_int(aa, bb)
            if d <= 0 or not is_square(d):
                continue
            r = classify_poly(x**5 + aa * x + bb, do_galois=True)
            sample_gal = {"t": str(tv), "alpha": aa, "beta": bb, "status": r.get("status"), "gal": r.get("galois")}
            if (r.get("status") or "").startswith("HIT_A5"):
                break
        except Exception:
            continue

    return {
        "id": fam["id"],
        "note": fam.get("note"),
        "locus": locus,
        "local_monodromy": [
            {
                "center": str(lm.get("center_N", lm.get("center"))),
                "cycle_type": lm.get("cycle_type"),
                "class": lm.get("class"),
                "tracking_error": lm.get("tracking_error"),
                "error": lm.get("error"),
            }
            for lm in locals_m
        ],
        "monodromy_infinity": {
            "cycle_type": minf.get("cycle_type"),
            "class": minf.get("class"),
            "tracking_error": minf.get("tracking_error"),
            "error": minf.get("error"),
        },
        "class_signature_finite": match_fin,
        "class_signature_with_inf": match_all,
        "group_from_finite": gen,
        "group_from_finite_and_inf": gen_inf,
        "catalogue_hits": hits,
        "catalogue_k": cat_k,
        "multi_k": len(cat_k) >= 2,
        "sample_gal": sample_gal,
        "nielsen_id": (
            match_all["matches"][0]["nielsen"]
            if match_all["matches"]
            else (match_fin["matches"][0]["nielsen"] if match_fin["matches"] else None)
        ),
    }


# ---------------------------------------------------------------------------
# E6 — 2-param envelope structure
# ---------------------------------------------------------------------------
def envelope_structure() -> dict:
    m, k = sp.symbols("m k", nonzero=True)
    alpha = pure_even_alpha(m, k)
    beta = pure_even_beta(m, k)
    D = sp.expand(sp.together(256 * alpha**5 + 3125 * beta**4))
    ideal = sp.expand(sp.together((256 * alpha**2 * m) ** 2))
    id_ok = sp.expand(D - ideal) == 0
    # branch when α=0: 256 m² = 3125 k⁴/256
    branch_rel = sp.together(256 * m**2 - sp.Rational(3125) * k**4 / 256)
    return {
        "dimension": 2,
        "parameters": ["m", "k"],
        "disc_identical_square": id_ok,
        "branch_divisor": "α(m,k)=0 (and m=0 degenerate)",
        "branch_equation": str(branch_rel) + " = 0",
        "foliation": "fixed-k rays are pure-even 1-param; cross-k paths join multi-seed ratios",
        "note": (
            "The envelope is a rational surface of pure-even BJ polynomials. "
            "Geometric monodromy is well-defined on 1-param slices (paths/rays). "
            "A full Nielsen ID is an ID of those slices' monodromy types."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G3 — monodromy ID of the pure-even envelope", flush=True)
    print("=" * 72, flush=True)

    print("\n[E6] Envelope structure ...", flush=True)
    env = envelope_structure()
    print(f"  disc_id_square={env['disc_identical_square']} dim={env['dimension']}", flush=True)

    families = [
        family_flag_classical(),
        family_flag_lsw(),
        family_lsw_ray(),
        family_flagship_homog(),
        family_classical_homog(),
    ]

    results = []
    for fam in families:
        try:
            results.append(analyze_family(fam))
        except Exception as e:
            print(f"  FAIL {fam['id']}: {e}", flush=True)
            results.append({"id": fam["id"], "error": str(e)})

    elapsed = round(time.time() - t0, 2)

    # Global summary
    ids = []
    multi_paths = []
    for r in results:
        if r.get("nielsen_id"):
            ids.append((r["id"], r["nielsen_id"]))
        if r.get("multi_k"):
            multi_paths.append(r["id"])

    # Consistency of signatures across multi-k paths
    sigs = {
        r["id"]: r.get("class_signature_with_inf", {}).get("signature")
        for r in results
        if "class_signature_with_inf" in r
    }

    geometric_named = len(ids) > 0
    # Strong success: a multi-k path has a shortlist Nielsen name
    strong = any(
        r.get("multi_k") and r.get("nielsen_id")
        for r in results
    )

    verdict = (
        f"G3 envelope monodromy ({elapsed}s). "
        f"Envelope disc□ identity={env['disc_identical_square']}. "
        f"Multi-k paths: {multi_paths}. "
        f"Nielsen names found: {ids}. "
        f"Signatures: {sigs}. "
        f"Geometric multi-k named={strong}."
    )
    print("\n" + verdict, flush=True)

    # Report
    lines = [
        "# G3 — monodromy identification of the pure-even envelope",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Goal",
        "",
        "The pure-even envelope",
        "",
        "```text",
        "α = 256 m² − 3125 k⁴ / 256,   β = k α",
        "disc(x⁵+αx+β) = (256 α² m)²   (identical square)",
        "```",
        "",
        "already supplies **arithmetic multi-k** (cross-k paths hit flagship, classical, LSW, …).",
        "G3 asks for a **Nielsen name**: the geometric monodromy of natural 1-parameter",
        "slices as a permutation representation of π₁(ℙ¹ ∖ branch locus) → A₅,",
        "compared to the r=4 shortlist (3A⁴, 2A3A³, 2A²3A², …).",
        "",
        "---",
        "",
        "## 1. Envelope structure (2-parameter)",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| dimension | {env['dimension']} |",
        f"| disc identical square | **{env['disc_identical_square']}** |",
        f"| branch divisor | {env['branch_divisor']} |",
        f"| branch equation | `{env['branch_equation']}` |",
        f"| foliation | {env['foliation']} |",
        "",
        env["note"],
        "",
        "---",
        "",
        "## 2. One-parameter slices — monodromy",
        "",
    ]

    for r in results:
        lines.append(f"### `{r.get('id')}`")
        lines.append("")
        if r.get("error"):
            lines.append(f"**Error:** {r['error']}")
            lines.append("")
            continue
        lines.append(f"_{r.get('note')}_")
        lines.append("")
        loc = r.get("locus", {})
        lines.append(
            f"- disc degree: **{loc.get('disc_degree')}**; "
            f"sqf factors: `{loc.get('sqf_factors')}`; "
            f"all mult even: **{loc.get('all_mult_even')}**"
        )
        lines.append(
            f"- finite monodromy classes: "
            f"**{r.get('class_signature_finite', {}).get('signature')}**"
        )
        lines.append(
            f"- with ∞: **{r.get('class_signature_with_inf', {}).get('signature')}**"
        )
        lines.append(
            f"- group ⟨finite⟩: order={r.get('group_from_finite', {}).get('order')} "
            f"A5={r.get('group_from_finite', {}).get('is_A5')} "
            f"all_even={r.get('group_from_finite', {}).get('all_even')}"
        )
        lines.append(
            f"- group ⟨finite+∞⟩: order={r.get('group_from_finite_and_inf', {}).get('order')} "
            f"A5={r.get('group_from_finite_and_inf', {}).get('is_A5')}"
        )
        lines.append(
            f"- Nielsen match (with ∞): "
            f"**{r.get('nielsen_id') or r.get('class_signature_with_inf', {}).get('matches') or 'none'}**"
        )
        lines.append(
            f"- catalogue hits: {r.get('catalogue_hits')}; "
            f"k={r.get('catalogue_k')}; multi-k=**{r.get('multi_k')}**"
        )
        if r.get("sample_gal"):
            lines.append(f"- sample fibre Gal: {r['sample_gal']}")
        lines.append("")
        lines.append("Local monodromy (finite):")
        lines.append("")
        lines.append("| center | cycle type | class | track err |")
        lines.append("|--------|------------|-------|----------:|")
        for lm in r.get("local_monodromy", [])[:15]:
            lines.append(
                f"| {lm.get('center')} | {lm.get('cycle_type')} | {lm.get('class')} | "
                f"{lm.get('tracking_error')} |"
            )
        mi = r.get("monodromy_infinity", {})
        lines.append("")
        lines.append(
            f"Infinity: type={mi.get('cycle_type')} class=**{mi.get('class')}** "
            f"err={mi.get('tracking_error')}"
        )
        lines.append("")

    lines += [
        "---",
        "",
        "## 3. Nielsen shortlist comparison",
        "",
        "Target classes (from `A5_HURWITZ_R4.md`):",
        "",
        "| Nielsen type | class multiset | orbit | g |",
        "|--------------|----------------|------:|--:|",
        "| 3A⁴ | 3A×4 | 18 | 0 |",
        "| 2A 3A³ | 2A+3A×3 | 96 | 0 |",
        "| 2A² 3A² | 2A×2+3A×2 | 108 | 0 |",
        "| 2A 3A² 5* | 2A+3A×2+5* | 240 | 0 |",
        "| 3A³ 5* | 3A×3+5* | 40–60 | 0 |",
        "",
        "A **name** is reported when the multiset of local monodromy conjugacy classes",
        "(finite branch points, optionally including ∞) equals a shortlist multiset.",
        "",
        f"| family | multi-k? | signature (w/ ∞) | Nielsen ID |",
        f"|--------|:--------:|------------------|------------|",
    ]
    for r in results:
        if r.get("error"):
            lines.append(f"| {r['id']} | ? | error | — |")
            continue
        lines.append(
            f"| {r['id']} | {r.get('multi_k')} | "
            f"{r.get('class_signature_with_inf', {}).get('signature')} | "
            f"{r.get('nielsen_id') or '—'} |"
        )

    lines += [
        "",
        "---",
        "",
        "## 4. Multi-k / geometric conclusion",
        "",
        f"| test | result |",
        f"|------|--------|",
        f"| Pure-even envelope disc identity | **{env['disc_identical_square']}** |",
        f"| Multi-k arithmetic paths present | **{len(multi_paths) > 0}** ({multi_paths}) |",
        f"| Local monodromy computed | **True** |",
        f"| Nielsen shortlist name on any family | **{geometric_named}** ({ids}) |",
        f"| **Named geometric multi-k** (multi-k path + Nielsen ID) | **{strong}** |",
        "",
        (
            "**Success:** a multi-k pure-even path carries an explicit Nielsen label — "
            "arithmetic multi-k is upgraded to geometric multi-k for that slice."
            if strong
            else (
                "**Partial / open:** monodromy signatures are computed, but either they do not "
                "match the r=4 shortlist exactly, or multi-k paths remain unnamed. "
                "Possible causes: (i) more than 4 effective branch points (r≥5 Nielsen type); "
                "(ii) 5A vs 5B not resolved; (iii) numerical tracking merges/splits loops; "
                "(iv) the envelope slice is a pullback / not a primitive Hurwitz curve."
            )
        ),
        "",
        "### What this cut established",
        "",
        "1. Explicit disc/branch analysis of pure-even 1-param slices (square multiplicities).",
        "2. Numerical local monodromy cycle types at finite branch points and at ∞.",
        "3. Comparison to the G1/G2 Nielsen shortlist.",
        "4. Group generation check (A5 vs smaller even groups).",
        "5. Parallel arithmetic multi-k catalogue confirmation on the same paths.",
        "",
        "### Next if unnamed",
        "",
        "1. Refine radii / certified monodromy (interval tracking) for ambiguous centers.",
        "2. Distinguish 5A vs 5B via lift invariant / complex conjugation.",
        "3. Allow r≥5 signatures (match against A5_HURWITZ r=5 type list).",
        "4. Compute the monodromy of the *discriminant double cover* / resolvent of the",
        "   Galois closure as a cover of the t-line (degree 60 or 120 → quotient).",
        "5. Literature: is the pure-even BJ envelope a known Hurwitz family?",
        "",
        "---",
        "",
        "## 5. Non-claims",
        "",
        "- Numerical monodromy is not a certified braid factorization.",
        "- Matching class multisets is necessary but not sufficient for a full Hurwitz",
        "  component ID (braid orbit / lift invariant may refine further).",
        "- Does not reopen pure-even arithmetic, Canonical T3, or Necessity.",
        "",
        "_Generated by `g3_envelope_monodromy.py`._",
        "",
    ]

    md = "\n".join(lines)
    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "envelope": env,
        "families": results,
        "nielsen_ids": ids,
        "multi_k_paths": multi_paths,
        "signatures": sigs,
        "geometric_multi_k_named": strong,
    }

    write_md(ROOT / "G3_ENVELOPE_MONODROMY.md", md)
    write_json(ROOT / "G3_ENVELOPE_MONODROMY.json", payload)
    write_md(OUT / "G3_ENVELOPE_MONODROMY.md", md)
    write_json(OUT / "G3_ENVELOPE_MONODROMY.json", payload)
    try:
        write_md(RESULTS / "G3_ENVELOPE_MONODROMY.md", md)
        write_json(RESULTS / "G3_ENVELOPE_MONODROMY.json", payload)
    except Exception:
        pass

    print(f"\nWrote G3_ENVELOPE_MONODROMY.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
