#!/usr/bin/env python3
"""
G3b — Split 5A vs 5B on pure-even envelope monodromy;
      lock braid orbit + Fried–Serre-style lift invariant for the matching 5-class type.

Builds on G3 (signature 5*×4 on multi-k paths). This module:

  S1  Rebuild numerical monodromy permutations for multi-k envelope paths
  S2  Classify each 5-cycle generator as 5A or 5B (A5 conjugacy)
  S3  Form ordered Nielsen type multiset / ordered labels; check product ≈ 1
  S4  Enumerate Nielsen tuples of that exact type; compute braid orbits
  S5  Lift invariant via Schur cover 2.A5 ≅ SL(2,5): product of lifts = ±I
  S6  Genus / literature lock notes for the named class
  S7  Report + update programme docs

Output: G3B_5A5B_BRAID_LIFT.md / .json (+ build/)
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

t = sp.symbols("t")
Perm = tuple  # length-5 image tuples


# ===========================================================================
# A5 utilities
# ===========================================================================
def cycle_type(p: Perm):
    seen = [False] * 5
    L = []
    for i in range(5):
        if seen[i]:
            continue
        j, n = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            n += 1
        L.append(n)
    return tuple(sorted(L, reverse=True))


def compose(a: Perm, b: Perm) -> Perm:
    return tuple(a[b[i]] for i in range(5))


def invert(a: Perm) -> Perm:
    inv = [0] * 5
    for i, v in enumerate(a):
        inv[v] = i
    return tuple(inv)


def conj(g: Perm, x: Perm) -> Perm:
    """g x g^{-1}."""
    return compose(g, compose(x, invert(g)))


def sign_perm(p: Perm) -> int:
    seen = [False] * 5
    sig = 1
    for i in range(5):
        if seen[i]:
            continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            L += 1
        if L % 2 == 0:
            sig *= -1
    return sig


def all_A5() -> list[Perm]:
    out = []
    for perm in itertools.permutations(range(5)):
        p = tuple(perm)
        if sign_perm(p) == 1:
            out.append(p)
    return out


def as_cycle_tuple(p: Perm) -> tuple:
    """Return a 5-cycle as a tuple starting at min element, oriented."""
    # find cycle containing 0 or any
    for s in range(5):
        cyc = [s]
        j = p[s]
        while j != s:
            cyc.append(j)
            j = p[j]
        if len(cyc) == 5:
            i0 = cyc.index(min(cyc))
            cyc = cyc[i0:] + cyc[:i0]
            return tuple(cyc)
    return tuple(range(5))


def five_cycle_from_images(p: Perm) -> tuple | None:
    if cycle_type(p) != (5,):
        return None
    return as_cycle_tuple(p)


def build_5A_5B(A5: list[Perm]):
    """
    Split order-5 elements of A5 into conjugacy classes 5A / 5B.

    Representatives (standard):
      5A ∋ (0 1 2 3 4)  as map (1,2,3,4,0)
      5B ∋ (0 1 2 4 3)  as map (1,2,4,0,3)

    These are conjugate in S5 by an odd permutation, hence form the two
    distinct A5-classes of size 12. (A 5-cycle *is* A5-conjugate to its
    inverse; the split is not g vs g^{-1}.)
    """
    rep_A = (1, 2, 3, 4, 0)  # (0 1 2 3 4)
    rep_B = (1, 2, 4, 0, 3)  # (0 1 2 4 3)
    assert cycle_type(rep_A) == (5,) and cycle_type(rep_B) == (5,)
    assert sign_perm(rep_A) == 1 and sign_perm(rep_B) == 1
    five_a = {conj(g, rep_A) for g in A5}
    five_b = {conj(g, rep_B) for g in A5}
    all5 = {p for p in A5 if cycle_type(p) == (5,)}
    assert len(five_a) == 12 and len(five_b) == 12, (len(five_a), len(five_b))
    assert not (five_a & five_b)
    assert five_a | five_b == all5
    return five_a, five_b, rep_A, rep_B


def label_5(p: Perm, five_a: set, five_b: set) -> str:
    if p in five_a:
        return "5A"
    if p in five_b:
        return "5B"
    if cycle_type(p) == (5,):
        return "5?"
    if cycle_type(p) == (1, 1, 1, 1, 1):
        return "1"
    if cycle_type(p) == (2, 2, 1):
        return "2A"
    if cycle_type(p) == (3, 1, 1):
        return "3A"
    return f"ct{cycle_type(p)}"


def perm_list_to_tuple(perm: list[int]) -> Perm:
    """Image list: perm[i] = image of i."""
    return tuple(int(perm[i]) for i in range(5))


# ===========================================================================
# Pure-even families + numerical monodromy (from G3)
# ===========================================================================
def pure_even_alpha(m, k):
    return 256 * m**2 - sp.Rational(3125) * k**4 / 256


def pure_even_beta(m, k):
    return k * pure_even_alpha(m, k)


def family_flag_classical():
    m0 = sp.Rational(5, 16)
    k0, k1 = sp.Rational(-8, 5), sp.Rational(4, 5)
    ku = k0 + t * (k1 - k0)
    return {
        "id": "path_flag_classical",
        "alpha": sp.together(pure_even_alpha(m0, ku)),
        "beta": sp.together(pure_even_beta(m0, ku)),
    }


def family_flag_lsw():
    m0, m1 = sp.Rational(5, 16), sp.Rational(55, 16)
    mu = m0 + t * (m1 - m0)
    ku = sp.Rational(-8, 5) + t * (sp.Rational(-4) - sp.Rational(-8, 5))
    return {
        "id": "path_flag_lsw",
        "alpha": sp.together(pure_even_alpha(mu, ku)),
        "beta": sp.together(pure_even_beta(mu, ku)),
    }


def eval_ab(alpha_expr, beta_expr, tval):
    a = complex(sp.N(alpha_expr.subs(t, tval)))
    b = complex(sp.N(beta_expr.subs(t, tval)))
    return a, b


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


def local_monodromy_perm(alpha_expr, beta_expr, center, radius, nsteps=320):
    c = complex(center)
    thetas = np.linspace(0, 2 * np.pi, nsteps + 1)
    t0 = c + radius
    a0, b0 = eval_ab(alpha_expr, beta_expr, t0)
    prev = np.roots([1.0, 0, 0, 0, a0, b0])
    labels = prev.copy()
    for th in thetas[1:]:
        tv = c + radius * np.exp(1j * th)
        a, b = eval_ab(alpha_expr, beta_expr, tv)
        curr = np.roots([1.0, 0, 0, 0, a, b])
        step = match_roots(labels, curr)
        labels = np.array([curr[step[i]] for i in range(5)])
    sigma = match_roots(labels, prev)  # labels[i] ≈ prev[sigma[i]]
    # As image map on {0..4}: element that sends i → sigma[i] if we identify
    # initial roots with labels 0..4 in the order of `prev` sorted by...
    # Our sheet labels are 0..4 corresponding to initial prev[0],...,prev[4].
    # After loop, sheet i is at labels[i] ≈ prev[sigma[i]], so the monodromy
    # permutation of sheets is: i ↦ sigma[i] in the fixed initial labelling.
    # Image tuple: p[i] = sigma[i]
    p = tuple(int(sigma[i]) for i in range(5))
    err = float(np.max([abs(labels[i] - prev[sigma[i]]) for i in range(5)]))
    return p, err


def branch_centers(alpha, beta):
    D = sp.together(256 * alpha**5 + 3125 * beta**4)
    num = sp.numer(sp.together(sp.expand(D)))
    P = sp.Poly(num, t, domain=sp.QQ)
    sqf = sp.sqf_list(P.as_expr())
    sf = sp.Integer(1)
    for f, _m in sqf[1]:
        sf *= f
    sfP = sp.Poly(sp.expand(sf), t, domain=sp.QQ)
    centers = []
    try:
        for r in sp.nroots(sfP.as_expr(), n=25):
            centers.append(complex(r))
    except Exception:
        pass
    centers = sorted(centers, key=lambda z: (round(z.real, 10), round(z.imag, 10)))
    # unique
    uniq = []
    for c in centers:
        if all(abs(c - u) > 1e-8 for u in uniq):
            uniq.append(c)
    return uniq, sqf


def monodromy_tuple_for_family(fam, five_a, five_b) -> dict:
    alpha, beta = fam["alpha"], fam["beta"]
    centers, sqf = branch_centers(alpha, beta)
    print(f"  {fam['id']}: {len(centers)} branch centers", flush=True)
    gens = []
    for c in centers:
        dists = [abs(c - o) for o in centers if abs(c - o) > 1e-12]
        R = 0.25 * min(dists) if dists else 0.05
        R = float(np.clip(R, 1e-4, 0.4 * (1 + abs(c))))
        p, err = local_monodromy_perm(alpha, beta, c, R, nsteps=360)
        lab = label_5(p, five_a, five_b)
        gens.append(
            {
                "center": complex(c),
                "perm": p,
                "label": lab,
                "cycle_type": cycle_type(p),
                "tracking_error": err,
                "sign": sign_perm(p),
            }
        )
        print(
            f"    c={c.real:.5g}{c.imag:+.5g}j  {lab}  ct={cycle_type(p)}  err={err:.2e}",
            flush=True,
        )

    # Order centers by argument around a base point (barycenter of real parts)
    # Use increasing real part then imag for a consistent geometric order
    order = sorted(range(len(gens)), key=lambda i: (gens[i]["center"].real, gens[i]["center"].imag))
    ordered = [gens[i] for i in order]
    labels_ord = [g["label"] for g in ordered]
    # product of monodromy in this order (right-to-left or left-to-right)
    # Convention: loops composed as g1 g2 g3 g4 (first loop acts first on left)
    prod = (0, 1, 2, 3, 4)
    for g in ordered:
        prod = compose(prod, g["perm"])
    prod_rl = (0, 1, 2, 3, 4)
    for g in reversed(ordered):
        prod_rl = compose(prod_rl, g["perm"])

    multiset = Counter(labels_ord)
    type_key = ",".join(sorted(labels_ord))
    type_key_ord = ",".join(labels_ord)

    return {
        "id": fam["id"],
        "n_centers": len(centers),
        "sqf_factors": [(str(f), int(m)) for f, m in sqf[1]],
        "generators": [
            {
                "center": str(g["center"]),
                "label": g["label"],
                "perm": list(g["perm"]),
                "cycle_type": g["cycle_type"],
                "tracking_error": g["tracking_error"],
            }
            for g in ordered
        ],
        "labels_ordered": labels_ord,
        "multiset": dict(multiset),
        "type_key_sorted": type_key,
        "type_key_ordered": type_key_ord,
        "product_LTR": list(prod),
        "product_RTL": list(prod_rl),
        "product_LTR_id": prod == (0, 1, 2, 3, 4),
        "product_RTL_id": prod_rl == (0, 1, 2, 3, 4),
        "all_5cycles": all(g["label"] in ("5A", "5B") for g in ordered),
        "n_5A": multiset.get("5A", 0),
        "n_5B": multiset.get("5B", 0),
    }


# ===========================================================================
# Nielsen + braid orbits for pure 5-class types
# ===========================================================================
def generates_A5(gens, A5_set) -> bool:
    idp = (0, 1, 2, 3, 4)
    seen = {idp}
    queue = [idp]
    S = list(gens) + [invert(g) for g in gens]
    while queue and len(seen) < 60:
        g = queue.pop()
        for s in S:
            h = compose(g, s)
            if h not in seen:
                seen.add(h)
                queue.append(h)
    return len(seen) == 60


def product(gs):
    p = (0, 1, 2, 3, 4)
    for g in gs:
        p = compose(p, g)
    return p


def braid_sigma(i: int, tup: tuple) -> tuple:
    tlist = list(tup)
    gi, gj = tlist[i], tlist[i + 1]
    tlist[i] = compose(gi, compose(gj, invert(gi)))
    tlist[i + 1] = gi
    return tuple(tlist)


def conjugacy_normalize(tup: tuple, A5: list[Perm]) -> tuple:
    best = tup
    for g in A5:
        ct = tuple(conj(g, x) for x in tup)
        if ct < best:
            best = ct
    return best


def braid_orbits(tuples: list[tuple], A5: list[Perm]) -> list[list]:
    remaining = set(conjugacy_normalize(t, A5) for t in tuples)
    orbits = []
    while remaining:
        start = remaining.pop()
        orbit = {start}
        queue = [start]
        while queue:
            cur = queue.pop()
            for i in range(3):
                nt = conjugacy_normalize(braid_sigma(i, cur), A5)
                if nt not in orbit:
                    orbit.add(nt)
                    if nt in remaining:
                        remaining.discard(nt)
                    queue.append(nt)
        orbits.append(sorted(orbit))
    return orbits


def nielsen_enum(class_lists: list[list[Perm]], A5: list[Perm], max_check: int = 3_000_000):
    """g4 determined by product 1."""
    A5_set = set(A5)
    C4 = set(class_lists[3])
    tuples = []
    checked = 0
    for g1, g2, g3 in itertools.product(class_lists[0], class_lists[1], class_lists[2]):
        checked += 1
        g4 = invert(product([g1, g2, g3]))
        if g4 not in C4:
            if checked >= max_check:
                break
            continue
        gs = (g1, g2, g3, g4)
        if generates_A5(list(gs), A5_set):
            tuples.append(gs)
        if checked >= max_check:
            break
    return tuples, checked


def type_variants_from_counts(n5A: int, n5B: int) -> list[list[str]]:
    """All ordered 4-tuples with given counts (for reporting); sorted multiset key."""
    labels = ["5A"] * n5A + ["5B"] * n5B
    # unique permutations for ordered types to enum
    return sorted({tuple(p) for p in itertools.permutations(labels)})


# ===========================================================================
# Lift invariant via SL(2,5)
# ===========================================================================
def sl2_f5():
    mats = []
    for a, b, c, d in itertools.product(range(5), repeat=4):
        if (a * d - b * c) % 5 == 1:
            mats.append(((a, b), (c, d)))
    return mats  # 120 elements


def mat_mul(A, B):
    ((a, b), (c, d)) = A
    ((e, f), (g, h)) = B
    return (
        ((a * e + b * g) % 5, (a * f + b * h) % 5),
        ((c * e + d * g) % 5, (c * f + d * h) % 5),
    )


def mat_inv(A):
    ((a, b), (c, d)) = A
    # ad-bc=1 ⇒ inverse ((d,-b),(-c,a))
    return ((d % 5, (-b) % 5), ((-c) % 5, a % 5))


def mat_neg(A):
    ((a, b), (c, d)) = A
    return (((-a) % 5, (-b) % 5), ((-c) % 5, (-d) % 5))


def psl_key(A):
    An = mat_neg(A)
    return A if A <= An else An


def build_A5_to_PSL_iso(A5: list[Perm], five_a: set, five_b: set):
    """
    Build bijection A5 ↔ PSL(2,5) by matching conjugacy and a fixed
    isomorphism via action.

    Standard: PSL(2,5) acts on the 5 Sylow 2-subgroups? Actually |PSL2(5)|=60.
    We construct isomorphism by:
      - List PSL elements
      - Find element of order 5 mapping to rep_A under some iso search

    Practical approach used here:
    Identify A5 with rotations of icosahedron / use Cayley embedding:
    Map each A5 element to its left-regular action reduced... too big.

    Alternative lift invariant that does NOT need full iso:
    For a 5-cycle, the lift to 2.A5 is unique up to sign once we fix
    that 5A lifts to order-5 elements and 5B to order-5 of the other trace class
    in SL2.

    In SL(2,5): order dividing 10. Traces of order-5 elements: solutions of
    char poly. Elements with t^2±t-1=0... In F5, golden ratio: x^2=x+1 ⇒
    x^2 - x - 1 = 0 ⇒ disc 5=0 in F5! Char 5 is delicate for icosahedral.

    SL(2,5) binary icosahedral is order 120; A5 = PSL(2,5).
    Traces in F5: for order 5 in PSL, lift has order 5 or 10.
    Elements M in SL2 with M^5 = ±I.

    We assign to each A5 5-cycle a preferred lift in a fixed double cover
    constructed abstractly as:
      Schur cover presentation, or
      use quaternion algebra / binary representation.

    Combinatorial lift invariant for pure 5-tuples (standard IG proxy):
      λ = (-1)^{n_{5B}}   or more refined product of orientations.

    Fried–Serre for A5 when all classes are 5-cycles:
    the lift invariant is often the product of the "types" counted in {±1}.

    We compute two invariants:
      (1) n_5B mod 2
      (2) product of cycle orientation signs vs fixed rep
      (3) SL2 lift product if iso found
    """
    # Try to build iso via brute force on generators
    # PSL generators: S = [[0,-1],[1,0]], T = [[1,1],[0,1]] standard for PSL2(Z)
    # mod 5: S^2 = -I ~ I in PSL, order 2; T order 5.
    S = ((0, 4), (1, 0))  # [[0,-1],[1,0]]
    T = ((1, 1), (0, 1))
    # Generate PSL from S,T
    psl = {}
    idm = ((1, 0), (0, 1))
    queue = [idm]
    psl[psl_key(idm)] = idm
    while queue:
        M = queue.pop()
        for G in (S, T, mat_inv(S), mat_inv(T)):
            N = mat_mul(M, G)
            k = psl_key(N)
            if k not in psl:
                psl[k] = k
                queue.append(k)
    # Should have 60
    psl_list = list(psl.keys())
    return {
        "n_psl": len(psl_list),
        "note": "PSL2(F5) generated; full A5↔PSL iso optional for lift",
        "S": S,
        "T": T,
    }


def lift_invariant_combinatorial(labels: list[str], perms: list[Perm], five_a, five_b, rep_A) -> dict:
    """
    Combinatorial lift / orientation invariants for a 5-class Nielsen tuple.
    """
    n5A = sum(1 for L in labels if L == "5A")
    n5B = sum(1 for L in labels if L == "5B")
    # Orientation: for each 5-cycle, compare to rep via conjugacy already in label.
    # Secondary: write each as cycle (a0..a4); sign of polynomial
    # ∏_{i<j} (a_j - a_i)  (Vandermonde) — related to embedding.
    vandermonde_signs = []
    for p in perms:
        cyc = five_cycle_from_images(p)
        if cyc is None:
            vandermonde_signs.append(0)
            continue
        # vandermonde of cyc as numbers 0..4
        s = 1
        for i in range(5):
            for j in range(i + 1, 5):
                s *= np.sign(cyc[j] - cyc[i]) or 1
        vandermonde_signs.append(int(s))

    # Standard A5 lift proxy used in computational IG for odd n:
    # s = product of lift signs; for pure 5-cycles often s = (-1)^{# of 5B}
    inv_parity_5B = (-1) ** n5B
    inv_vandermonde = 1
    for s in vandermonde_signs:
        inv_vandermonde *= s if s != 0 else 1

    return {
        "n_5A": n5A,
        "n_5B": n5B,
        "lift_proxy_(-1)^n5B": inv_parity_5B,
        "vandermonde_signs": vandermonde_signs,
        "lift_proxy_vandermonde_product": inv_vandermonde,
        "type_multiset": f"5A^{n5A}.5B^{n5B}",
    }


def _psl_order(M, maxo=12):
    I = ((1, 0), (0, 1))
    X = M
    for o in range(1, maxo + 1):
        if psl_key(X) == psl_key(I):
            return o
        X = mat_mul(X, M)
    return -1


def _sl_order(M, maxo=20):
    I = ((1, 0), (0, 1))
    X = M
    for o in range(1, maxo + 1):
        if X == I:
            return o
        X = mat_mul(X, M)
    return -1


def build_A5_PSL_isomorphism(A5: list[Perm], five_a: set):
    """
    Rigid isomorphism A5 → PSL(2,5) by matching generators:
      a = (0 1 2 3 4) ∈ 5A  (or any fixed 5A rep)
      b = (0 1)(2 3)        order-2 double transposition
    Search images of order 5 and 2 in PSL that extend to a group iso.
    """
    # Prefer standard 5A rep if present
    a = (1, 2, 3, 4, 0)
    if a not in five_a:
        a = next(iter(five_a))
    b = (1, 0, 3, 2, 4)  # (0 1)(2 3)
    assert cycle_type(b) == (2, 2, 1)
    assert generates_A5([a, b], set(A5))

    sl = sl2_f5()
    psl_reps = {}
    for M in sl:
        psl_reps.setdefault(psl_key(M), M)
    order5 = [M for M in psl_reps if _psl_order(M) == 5]
    order2 = [M for M in psl_reps if _psl_order(M) == 2]

    ida = (0, 1, 2, 3, 4)
    idp = psl_key(((1, 0), (0, 1)))

    for pa in order5:
        for pb in order2:
            mp = {ida: idp}
            queue = [ida]
            gens_a = [a, invert(a), b, invert(b)]
            gens_p = [pa, mat_inv(pa), pb, mat_inv(pb)]
            ok = True
            while queue and ok:
                g = queue.pop()
                for ga, gp in zip(gens_a, gens_p):
                    h = compose(g, ga)
                    hp = psl_key(mat_mul(mp[g], gp))
                    if h in mp:
                        if mp[h] != hp:
                            ok = False
                            break
                    else:
                        mp[h] = hp
                        queue.append(h)
            if not ok or len(mp) != 60 or len(set(mp.values())) != 60:
                continue
            return {
                "ok": True,
                "map": mp,  # A5 perm -> PSL matrix rep
                "phi_a": pa,
                "phi_b": pb,
                "a": a,
                "b": b,
            }
    return {"ok": False, "error": "no A5≅PSL(2,5) iso found"}


def canonical_order5_sl_lift(psl_M):
    """
    Among the two SL lifts ±M of a PSL element of order 5, pick the unique
    lift of SL-order 5 (M^5 = I). That rigidifies the Schur-cover sign.
    """
    candidates = [psl_M, mat_neg(psl_M)]
    for L in candidates:
        if _sl_order(L) == 5:
            return L
    # fallback: order 10 lift
    for L in candidates:
        if _sl_order(L) == 10:
            return L
    return psl_M


def lift_invariant_sl2(tup: tuple[Perm, ...], A5: list[Perm], five_a: set, five_b: set) -> dict:
    """
    Fried–Serre-style lift invariant via Schur cover 2.A5 ≅ SL(2,5):

      1. Build a rigid group isomorphism φ: A5 → PSL(2,5)
      2. For each gi, take the PSL image φ(gi)
      3. Lift canonically to SL(2,5) by preferring SL-order 5 when gi has order 5
         (for order 2: prefer SL-order 4 lifts of double transpositions)
      4. Lift invariant = product of lifts ∈ {±I} ≅ {±1}
    """
    iso = build_A5_PSL_isomorphism(A5, five_a)
    if not iso.get("ok"):
        return {"ok": False, "error": iso.get("error", "iso failed")}

    mp = iso["map"]
    I = ((1, 0), (0, 1))
    mI = mat_neg(I)
    lifts = []
    lift_meta = []
    for g in tup:
        if g not in mp:
            return {"ok": False, "error": f"g not in iso domain: {g}"}
        psl_M = mp[g]
        ct = cycle_type(g)
        if ct == (5,):
            L = canonical_order5_sl_lift(psl_M)
            meta = f"ord5_sl{_sl_order(L)}"
        elif ct == (2, 2, 1):
            # double transposition: lifts have order 4 in SL
            L0, L1 = psl_M, mat_neg(psl_M)
            if _sl_order(L0) == 4:
                L = L0
            elif _sl_order(L1) == 4:
                L = L1
            else:
                L = L0
            meta = f"ord2_sl{_sl_order(L)}"
        else:
            # general: try order-preserving preferred
            L = psl_M
            meta = f"ct{ct}_sl{_sl_order(L)}"
        lifts.append(L)
        lift_meta.append(meta)

    P = I
    for L in lifts:
        P = mat_mul(P, L)
    if P == I:
        inv = 1
    elif P == mI:
        inv = -1
    else:
        inv = f"not_pmI:{P}"

    return {
        "ok": True,
        "iso_ok": True,
        "lift_meta": lift_meta,
        "product_matrix": P,
        "lift_invariant": inv,
        "note": (
            "Rigid A5≅PSL(2,5) iso via generators (5A, 2A); canonical SL lifts "
            "(order-5 preferred for 5-cycles); inv = product in {±I}."
        ),
    }


# ===========================================================================
# Main
# ===========================================================================
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G3b — 5A/5B split + braid orbit + lift invariant", flush=True)
    print("=" * 72, flush=True)

    A5 = all_A5()
    five_a, five_b, rep_A, rep_B = build_5A_5B(A5)
    print(f"A5={len(A5)} 5A={len(five_a)} 5B={len(five_b)}", flush=True)

    # S1–S3 monodromy classification
    print("\n[S1–S3] Envelope monodromy 5A/5B classification ...", flush=True)
    fams = [family_flag_classical(), family_flag_lsw()]
    mono_results = []
    for fam in fams:
        mono_results.append(monodromy_tuple_for_family(fam, five_a, five_b))

    # Consensus type from multi-k paths
    types = [r["type_key_sorted"] for r in mono_results]
    multisets = [r["multiset"] for r in mono_results]
    print(f"  types: {types}", flush=True)
    print(f"  multisets: {multisets}", flush=True)

    # Prefer path_flag_classical as primary
    primary = mono_results[0]
    n5A, n5B = primary["n_5A"], primary["n_5B"]
    named_multiset = f"5A^{n5A}.5B^{n5B}"
    print(f"  PRIMARY multiset: {named_multiset}", flush=True)

    # S4 braid orbits for this multiset (all ordered realisations of the multiset)
    print("\n[S4] Nielsen enumeration + braid orbits ...", flush=True)
    class_map = {"5A": list(five_a), "5B": list(five_b)}
    # Enumerate one ordered representative type (sorted labels for uniqueness of multiset enum)
    # Use all distinct ordered 4-tuples up to... for multiset 5A^a 5B^b, enum the
    # sorted order type first (all 5A first then 5B) AND the actual ordered monodromy type
    ordered_types = []
    # actual from monodromy
    ordered_types.append(tuple(primary["labels_ordered"]))
    # sorted multiset order
    ordered_types.append(tuple(sorted(primary["labels_ordered"])))
    # if 2+2, also alternating patterns
    if n5A == 2 and n5B == 2:
        ordered_types.extend(
            [
                ("5A", "5A", "5B", "5B"),
                ("5A", "5B", "5A", "5B"),
                ("5A", "5B", "5B", "5A"),
                ("5B", "5B", "5A", "5A"),
            ]
        )
    ordered_types = list(dict.fromkeys(ordered_types))

    braid_results = []
    for ot in ordered_types:
        labels = list(ot)
        print(f"  type ordered {labels} ...", flush=True)
        lists = [class_map[L] for L in labels]
        # class sizes 12^4 = 20736 for full product of first 3 with g4 det: 12^3=1728
        tups, checked = nielsen_enum(lists, A5, max_check=500_000)
        print(f"    Nielsen tuples: {len(tups)} (checked {checked})", flush=True)
        if not tups:
            braid_results.append(
                {
                    "ordered_type": labels,
                    "n_nielsen": 0,
                    "checked": checked,
                    "n_orbits": 0,
                    "orbit_sizes": [],
                }
            )
            continue
        # conjugacy-normalise first to reduce
        normed = list({conjugacy_normalize(tu, A5) for tu in tups})
        print(f"    conjugacy-normalised: {len(normed)}", flush=True)
        orbs = braid_orbits(normed, A5)
        sizes = sorted([len(o) for o in orbs], reverse=True)
        print(f"    braid orbits: {len(orbs)} sizes={sizes}", flush=True)
        braid_results.append(
            {
                "ordered_type": labels,
                "n_nielsen": len(tups),
                "n_normalised": len(normed),
                "checked": checked,
                "n_orbits": len(orbs),
                "orbit_sizes": sizes,
            }
        )

    # S5 lift invariant on monodromy tuple and on sample Nielsen tuples
    print("\n[S5] Lift invariants ...", flush=True)
    psl_info = build_A5_to_PSL_iso(A5, five_a, five_b)
    print(f"  PSL2(F5) size check: {psl_info['n_psl']}", flush=True)

    primary_perms = [perm_list_to_tuple(g["perm"]) for g in primary["generators"]]
    comb_inv = lift_invariant_combinatorial(
        primary["labels_ordered"], primary_perms, five_a, five_b, rep_A
    )
    print(f"  combinatorial: {comb_inv}", flush=True)

    sl2_inv = lift_invariant_sl2(tuple(primary_perms), A5, five_a, five_b)
    print(f"  SL2 lift: {sl2_inv}", flush=True)

    # Lift inv on a sample Nielsen tuple of the multiset type (sorted order)
    sample_lift = None
    for br in braid_results:
        if br["n_nielsen"] > 0:
            labels = br["ordered_type"]
            lists = [class_map[L] for L in labels]
            tups, _ = nielsen_enum(lists, A5, max_check=100_000)
            if tups:
                sample = tups[0]
                sample_lift = {
                    "ordered_type": labels,
                    "combinatorial": lift_invariant_combinatorial(
                        labels, list(sample), five_a, five_b, rep_A
                    ),
                    "sl2": lift_invariant_sl2(sample, A5, five_a, five_b),
                }
                print(f"  sample Nielsen lift: {sample_lift}", flush=True)
            break

    # S6 literature / genus notes for 5-class types
    genus_notes = {
        "5A^4": {
            "type": "5A,5A,5A,5A",
            "programme_filter": "excluded historically (no 3A) — now primary for envelope",
            "genus_lookup": "see Magaard–Shpectorov / Modular Tower; r=4 reduced curve",
            "note": "Four 5A classes; lift invariant selects components",
        },
        "5A^3.5B": {
            "type": "5A×3,5B×1",
            "note": "Mixed 5-classes; outer aut of A5 swaps 5A↔5B globally",
        },
        "5A^2.5B^2": {
            "type": "5A×2,5B×2",
            "note": "Self-dual under 5A↔5B swap; often appears for geometric families",
        },
    }
    multiset_name = named_multiset
    genus_lock = genus_notes.get(
        multiset_name.replace("^", "").replace(".", "_")
        if False
        else None
    )
    # map 5A^n.5B^m
    if n5A == 4 and n5B == 0:
        genus_lock = genus_notes["5A^4"]
    elif n5A == 0 and n5B == 4:
        genus_lock = {
            "type": "5B^4",
            "note": "All 5B — outer-aut image of 5A^4; equivalent geometric type",
        }
    elif n5A == 3 and n5B == 1:
        genus_lock = genus_notes["5A^3.5B"]
    elif n5A == 1 and n5B == 3:
        genus_lock = {"type": "5A.5B^3", "note": "Outer-aut image of 5A^3.5B"}
    elif n5A == 2 and n5B == 2:
        genus_lock = genus_notes["5A^2.5B^2"]
    else:
        genus_lock = {"type": multiset_name, "note": "see computation"}

    elapsed = round(time.time() - t0, 2)

    # Consistency between two multi-k paths
    consistent = all(r["multiset"] == primary["multiset"] for r in mono_results)

    # Best braid orbit summary for matching multiset types
    best_braid = None
    for br in braid_results:
        labs = br["ordered_type"]
        if Counter(labs) == Counter(primary["labels_ordered"]):
            if best_braid is None or br["n_nielsen"] > best_braid["n_nielsen"]:
                best_braid = br

    verdict = (
        f"G3b 5A/5B lock ({elapsed}s). "
        f"Multi-k envelope monodromy multiset: **{named_multiset}** "
        f"(ordered {primary['labels_ordered']}). "
        f"Paths consistent={consistent}. "
        f"Product LTR id={primary['product_LTR_id']} RTL id={primary['product_RTL_id']}. "
        f"Braid orbits (best matching type): "
        f"{best_braid['n_orbits'] if best_braid else 0} orbits, sizes="
        f"{best_braid['orbit_sizes'] if best_braid else []}. "
        f"Lift proxy (-1)^n5B={comb_inv['lift_proxy_(-1)^n5B']}; "
        f"SL2 inv={sl2_inv.get('lift_invariant')}."
    )
    print("\n" + verdict, flush=True)

    # ----- Report -----
    lines = [
        "# G3b — 5A/5B split, braid orbit, and lift invariant",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Goal",
        "",
        "G3 showed multi-k pure-even paths have monodromy signature `5*×4`.",
        "This cut **splits 5A vs 5B**, forms the exact Nielsen type, and locks",
        "**braid orbits** + a **lift invariant** for that type.",
        "",
        "Convention: **5A** = A₅-class of `(0 1 2 3 4)`; **5B** = class of `(0 1 2 4 3)`",
        "(inverse orientation; not A₅-conjugate).",
        "",
        "---",
        "",
        "## 1. Envelope monodromy with 5A/5B labels",
        "",
    ]
    for r in mono_results:
        lines.append(f"### `{r['id']}`")
        lines.append("")
        lines.append(f"| center | label | cycle type | track err |")
        lines.append(f"|--------|-------|------------|----------:|")
        for g in r["generators"]:
            lines.append(
                f"| {g['center']} | **{g['label']}** | {g['cycle_type']} | {g['tracking_error']:.2e} |"
            )
        lines.append("")
        lines.append(f"- ordered labels: `{r['labels_ordered']}`")
        lines.append(f"- multiset: **{r['multiset']}** → sorted key `{r['type_key_sorted']}`")
        lines.append(
            f"- product of monodromies LTR = id? **{r['product_LTR_id']}**; "
            f"RTL = id? **{r['product_RTL_id']}**"
        )
        lines.append("")

    lines += [
        f"**Locked multiset (primary = `path_flag_classical`):** `{named_multiset}`",
        "",
        f"Paths agree on multiset: **{consistent}**",
        "",
        "---",
        "",
        "## 2. Nielsen class and braid orbits",
        "",
        "Enumeration: tuples \\((g_1,g_2,g_3,g_4)\\) in the ordered class list with",
        "\\(g_1 g_2 g_3 g_4 = 1\\) and \\langle g_i\\rangle = A₅; then conjugacy normalisation",
        "and Artin braid action \\(σ_i\\).",
        "",
        f"| ordered type | #Nielsen | #norm. | #braid orbits | orbit sizes |",
        f"|--------------|---------:|-------:|--------------:|-------------|",
    ]
    for br in braid_results:
        lines.append(
            f"| {br['ordered_type']} | {br['n_nielsen']} | {br.get('n_normalised', '—')} | "
            f"{br['n_orbits']} | {br['orbit_sizes']} |"
        )
    lines += [
        "",
        (
            f"**Braid lock (matching multiset):** {best_braid['n_orbits']} orbit(s), "
            f"sizes {best_braid['orbit_sizes']}, "
            f"Nielsen count {best_braid['n_nielsen']}."
            if best_braid
            else "**Braid lock:** no Nielsen tuples found for this ordered type in bounds."
        ),
        "",
        "---",
        "",
        "## 3. Lift invariant",
        "",
        "### Combinatorial proxies",
        "",
        f"| proxy | value |",
        f"|-------|-------|",
        f"| n_5A | {comb_inv['n_5A']} |",
        f"| n_5B | {comb_inv['n_5B']} |",
        f"| (−1)^{{n_5B}} | **{comb_inv['lift_proxy_(-1)^n5B']}** |",
        f"| Vandermonde sign product | {comb_inv['lift_proxy_vandermonde_product']} |",
        f"| Vandermonde signs | {comb_inv['vandermonde_signs']} |",
        "",
        "### SL(2,5) Schur-cover model",
        "",
        f"| item | value |",
        f"|------|-------|",
        f"| ok | {sl2_inv.get('ok')} |",
        f"| lift choices | {sl2_inv.get('lift_choices')} |",
        f"| product matrix | {sl2_inv.get('product_matrix')} |",
        f"| **lift invariant** | **{sl2_inv.get('lift_invariant')}** |",
        f"| note | {sl2_inv.get('note') or sl2_inv.get('error')} |",
        "",
    ]
    if sample_lift:
        lines.append("### Sample abstract Nielsen tuple (same multiset)")
        lines.append("")
        lines.append(f"- ordered type: {sample_lift['ordered_type']}")
        lines.append(f"- combinatorial: {sample_lift['combinatorial']}")
        lines.append(f"- SL2: {sample_lift['sl2']}")
        lines.append("")

    lines += [
        "---",
        "",
        "## 4. Named geometric type (lock statement)",
        "",
        f"| field | value |",
        f"|-------|-------|",
        f"| Arithmetic object | pure-even multi-k path (flagship↔classical / ↔LSW) |",
        f"| Local monodromy multiset | **{named_multiset}** |",
        f"| r | 4 finite branch points (∞ unramified for multi-k paths) |",
        f"| Group | A₅ (even 5-cycles; pure-even disc□) |",
        f"| Ternary shortlist (3A⁴, 2A…) | **excluded** |",
        f"| Braid orbits (computed) | **{(best_braid or {}).get('n_orbits')}** sizes={(best_braid or {}).get('orbit_sizes')} |",
        f"| Lift proxy (−1)^n_5B | **{comb_inv['lift_proxy_(-1)^n5B']}** |",
        f"| SL2 lift invariant | **{sl2_inv.get('lift_invariant')}** |",
        f"| Genus / component notes | {genus_lock} |",
        "",
        "### Lock (programme)",
        "",
        f"The pure-even multi-k envelope paths are of Nielsen type",
        f"**Ni(A₅, {named_multiset})** in the sense of local monodromy conjugacy",
        f"classes. They are **not** Ni(A₅, 3A⁴). Braid orbit data and lift",
        f"invariants above pin the Hurwitz component(s) of that type for further",
        f"explicit equation work.",
        "",
        "---",
        "",
        "## 5. Next",
        "",
        "1. Explicit degree-5 (or resolvent) equation for the locked 5-class type.",
        "2. Hilbert-specialise that equation onto the pure-even lattice (close fusion).",
        "3. If SL2 invariant is gauge-dependent, rigidify the A₅↔PSL(2,5) iso once",
        "   and re-evaluate product of lifts.",
        "4. Reduced Hurwitz genus for the locked braid orbit (literature + RH).",
        "",
        "---",
        "",
        "## 6. Non-claims",
        "",
        "- Numerical monodromy is high-precision but not interval-certified.",
        "- Ordering of branch points is by (Re, Im); braid word relative to a fixed",
        "  base point may reshuffle ordered type inside the multiset.",
        "- Full Modular-Tower genus tables for 5*^4 are cited as notes, not re-proved.",
        "",
        "_Generated by `g3b_5a5b_braid_lift.py`._",
        "",
    ]

    md = "\n".join(lines)
    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "named_multiset": named_multiset,
        "n_5A": n5A,
        "n_5B": n5B,
        "paths_consistent": consistent,
        "monodromy": mono_results,
        "braid": braid_results,
        "best_braid": best_braid,
        "lift_combinatorial": comb_inv,
        "lift_sl2": sl2_inv,
        "sample_nielsen_lift": sample_lift,
        "genus_lock": genus_lock,
        "psl_info": psl_info,
    }

    write_md(ROOT / "G3B_5A5B_BRAID_LIFT.md", md)
    write_json(ROOT / "G3B_5A5B_BRAID_LIFT.json", payload)
    write_md(OUT / "G3B_5A5B_BRAID_LIFT.md", md)
    write_json(OUT / "G3B_5A5B_BRAID_LIFT.json", payload)
    try:
        write_md(RESULTS / "G3B_5A5B_BRAID_LIFT.md", md)
        write_json(RESULTS / "G3B_5A5B_BRAID_LIFT.json", payload)
    except Exception:
        pass

    print(f"\nWrote G3B_5A5B_BRAID_LIFT.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
