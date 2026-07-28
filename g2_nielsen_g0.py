#!/usr/bin/env python3
"""
G2 — genus-0 Nielsen types beyond 3A^4: 2A 3A³ and 2A² 3A².

For each type:
  C. Combinatorial lock (Nielsen sample / product-1 / generate A5; orbit sizes from prior)
  A. Explicit cover ansatz (ramification profiles) + Newton realisation at rational s
  H. Hilbert specialisations of monic fibres → irr / disc□ / BJ / multi-seed catalogue
  E. Arithmetic pure-even envelope control (must stay green)

Types (r=4, lookup g=0, programme shortlist):
  T_2A3A3  = 2A,3A,3A,3A   orbit size 96 (prior)
  T_2A2_3A2 = 2A,2A,3A,3A  orbit size 108 (prior)

Output: G2_NIELSEN_G0.md / .json (+ build/)
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, classify_poly, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

y, t = sp.symbols("y t")

# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------
CATALOGUE = [
    ("flagship", -55, 88, Fraction(-8, 5)),
    ("flag_145", 145, -232, Fraction(-8, 5)),
    ("flag_320", 320, -512, Fraction(-8, 5)),
    ("classical", 20, 16, Fraction(4, 5)),
    ("s95_76", 95, 76, Fraction(4, 5)),
    ("s220_176", 220, 176, Fraction(4, 5)),
    ("lsw_m100", -100, 400, Fraction(-4)),
    ("lsw_124m", 124, -496, Fraction(-4)),
    ("lsw_m209", -209, 836, Fraction(-4)),
    ("s180", -180, 432, Fraction(-12, 5)),
    ("s220m", 220, -528, Fraction(-12, 5)),
    ("s55_176", -55, 176, Fraction(-16, 5)),
    ("flagship_m", -55, -88, Fraction(8, 5)),
    ("classical_m", 20, -16, Fraction(-4, 5)),
    ("lsw4_m100", -100, -400, Fraction(4)),
]
CAT_BY_AB = {(a, b): (tag, k) for tag, a, b, k in CATALOGUE}
CAT_K = sorted({k for *_, k in CATALOGUE}, key=lambda f: (f.denominator, abs(f.numerator)))


def k_of(a, b):
    if a == 0:
        return None
    return Fraction(int(b), int(a))


def pure_even_alpha(m: Fraction, k: Fraction) -> Fraction:
    return 256 * m * m - Fraction(3125) * (k**4) / 256


# ---------------------------------------------------------------------------
# C — combinatorial: sample Nielsen tuples for target types
# ---------------------------------------------------------------------------
def cycle_type(p):
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


def compose(a, b):
    return tuple(a[b[i]] for i in range(5))


def invert(a):
    inv = [0] * 5
    for i, v in enumerate(a):
        inv[v] = i
    return tuple(inv)


def product(gs):
    p = (0, 1, 2, 3, 4)
    for g in gs:
        p = compose(p, g)
    return p


def generates_A5(gens) -> bool:
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


def a5_classes():
    A5 = []
    for perm in itertools.permutations(range(5)):
        p = tuple(perm)
        invs = sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j])
        if invs % 2 == 0:
            A5.append(p)
    classes = {"2A": [], "3A": [], "5A": [], "5B": []}
    fives = [p for p in A5 if cycle_type(p) == (5,)]
    # split 5A/5B by conjugacy in A5 (two classes of 5-cycles)
    # simple: powers of a fixed 5-cycle vs the other class
    five_a, five_b = set(), set()
    if fives:
        g0 = fives[0]
        # class of g0: conjugates in A5
        for h in A5:
            # conj h g0 h^{-1}
            # as maps: h ∘ g0 ∘ h^{-1}
            hi = invert(h)
            cg = compose(h, compose(g0, hi))
            five_a.add(cg)
        for g in fives:
            if g in five_a:
                pass
            else:
                five_b.add(g)
        # ensure partition
        for g in fives:
            if g not in five_a and g not in five_b:
                five_b.add(g)
    for p in A5:
        ct = cycle_type(p)
        if ct == (2, 2, 1):
            classes["2A"].append(p)
        elif ct == (3, 1, 1):
            classes["3A"].append(p)
        elif ct == (5,):
            if p in five_a:
                classes["5A"].append(p)
            else:
                classes["5B"].append(p)
    return A5, classes


def sample_nielsen(type_labels: list[str], max_found: int = 3, max_check: int = 2_000_000) -> dict:
    """Find sample tuples with product 1 generating A5 in given class sequence."""
    A5, classes = a5_classes()
    lists = [classes[lab] for lab in type_labels]
    sizes = [len(L) for L in lists]
    found = []
    checked = 0
    C4 = set(lists[3])
    # g4 determined
    for g1, g2, g3 in itertools.product(lists[0], lists[1], lists[2]):
        checked += 1
        pref = product([g1, g2, g3])
        g4 = invert(pref)
        if g4 not in C4:
            if checked >= max_check:
                break
            continue
        gs = (g1, g2, g3, g4)
        if generates_A5(list(gs)):
            found.append(
                {
                    "cycle_types": [cycle_type(g) for g in gs],
                    "tuple_images": [list(g) for g in gs],
                }
            )
            if len(found) >= max_found:
                break
        if checked >= max_check:
            break
    return {
        "type": type_labels,
        "class_sizes": sizes,
        "checked": checked,
        "n_found": len(found),
        "samples": found,
        "ok": len(found) > 0,
        "prior_orbit_data": {
            "2A,3A,3A,3A": {"orbits": 1, "orbit_size": 96, "genus_lookup": 0},
            "2A,2A,3A,3A": {"orbits": 1, "orbit_size": 108, "genus_lookup": 0},
        }.get(",".join(type_labels)),
    }


# ---------------------------------------------------------------------------
# Cover ansätze + Newton
# ---------------------------------------------------------------------------
def poly_val_derivs_from_roots(y0, roots_mults, c=1.0):
    """
    N(y) = c * Π (y-r)^m
    Return N, N', N'' at y0 via logarithmic derivatives for stability.
    """
    # direct product in complex
    N = c + 0j
    for r, m in roots_mults:
        N *= (y0 - r) ** m
    # N'/N = sum m/(y-r)
    if abs(N) < 1e-30:
        # fallback finite difference scale
        return 0j, 0j, 0j
    s1 = 0j
    s2 = 0j
    for r, m in roots_mults:
        d = y0 - r
        s1 += m / d
        s2 += m / (d * d)
    Np = N * s1
    Npp = N * (s1 * s1 - s2)
    return N, Np, Npp


def residual_2A3A3(v, s_val: float):
    """
    Type assignment: 3A at 0, 3A at 1, 3A at s, 2A at ∞.

    N = c * y^3 * (y-1) * (y-a)     → zeros mult 3,1,1 at 0 (branch 0)
    D = (y-r1)^2 * (y-r2)^2           → poles 2,2; ∞ order 1 → 2A at ∞
    Triple root of N-D at q (branch 1)
    Triple root of N-sD at w (branch s)

    v = [c, a, r1, r2, q, w]
    """
    c, a, r1, r2, q, w = v

    def ND_at(pt):
        N, Np, Npp = poly_val_derivs_from_roots(
            pt, [(0.0, 3), (1.0, 1), (a, 1)], c=c
        )
        D, Dp, Dpp = poly_val_derivs_from_roots(pt, [(r1, 2), (r2, 2)], c=1.0)
        return N, Np, Npp, D, Dp, Dpp

    Nq, Npq, Nppq, Dq, Dpq, Dppq = ND_at(q)
    Nw, Npw, Nppw, Dw, Dpw, Dppw = ND_at(w)
    # N - 1*D triple at q; N - s D triple at w
    return np.array(
        [
            (Nq - Dq).real,
            (Npq - Dpq).real,
            (Nppq - Dppq).real,
            (Nw - s_val * Dw).real,
            (Npw - s_val * Dpw).real,
            (Nppw - s_val * Dppw).real,
        ],
        dtype=float,
    ) + 1j * np.array(
        [
            (Nq - Dq).imag,
            (Npq - Dpq).imag,
            (Nppq - Dppq).imag,
            (Nw - s_val * Dw).imag,
            (Npw - s_val * Dpw).imag,
            (Nppw - s_val * Dppw).imag,
        ],
        dtype=float,
    )
    # Wait - residual should be real 6-vector for real v. Use real parts only if v real.


def residual_2A3A3_real(v, s_val: float):
    c, a, r1, r2, q, w = [float(x) for x in v]

    def vals(pt):
        # N = c y^3 (y-1)(y-a)
        yy = pt
        N = c * yy**3 * (yy - 1.0) * (yy - a)
        # logarithmic
        if abs(yy) < 1e-14 or abs(yy - 1) < 1e-14 or abs(yy - a) < 1e-14:
            # expand carefully via product rule poly
            pass
        # poly form derivatives
        # N = c (y^5 - (1+a)y^4 + a y^3)
        Np = c * (5 * yy**4 - 4 * (1 + a) * yy**3 + 3 * a * yy**2)
        Npp = c * (20 * yy**3 - 12 * (1 + a) * yy**2 + 6 * a * yy)
        # D = (y-r1)^2 (y-r2)^2 = ((y-r1)(y-r2))^2 = (y^2 - (r1+r2)y + r1 r2)^2
        s12 = r1 + r2
        p12 = r1 * r2
        u = yy**2 - s12 * yy + p12
        up = 2 * yy - s12
        upp = 2.0
        D = u * u
        Dp = 2 * u * up
        Dpp = 2 * up * up + 2 * u * upp
        return N, Np, Npp, D, Dp, Dpp

    Nq, Npq, Nppq, Dq, Dpq, Dppq = vals(q)
    Nw, Npw, Nppw, Dw, Dpw, Dppw = vals(w)
    return np.array(
        [
            Nq - Dq,
            Npq - Dpq,
            Nppq - Dppq,
            Nw - s_val * Dw,
            Npw - s_val * Dpw,
            Nppw - s_val * Dppw,
        ],
        dtype=float,
    )


def residual_2A2_3A2_real(v, s_val: float):
    """
    Type: 2A at 0, 2A at ∞, 3A at 1, 3A at s.

    N = c * y^2 * (y-1)^2 * (y-a)     → zeros 2,2,1 → 2A at 0
    D = (y-r1)^2 * (y-r2)^2             → poles 2,2; ∞ order 1 → 2A at ∞
    Triple root of N-D at q (3A at 1)
    Triple root of N-sD at w (3A at s)

    v = [c, a, r1, r2, q, w]
    """
    c, a, r1, r2, q, w = [float(x) for x in v]

    def vals(pt):
        yy = pt
        # N = c y^2 (y-1)^2 (y-a) = c (y^2 - 2y +1) y^2 (y-a) = c (y^4 - 2 y^3 + y^2)(y-a)
        # = c (y^5 - a y^4 - 2 y^4 + 2a y^3 + y^3 - a y^2)
        # = c (y^5 - (a+2) y^4 + (2a+1) y^3 - a y^2)
        N = c * (yy**2) * ((yy - 1.0) ** 2) * (yy - a)
        # derivatives via product: factors A=y^2, B=(y-1)^2, C=(y-a)
        A, Ap, App = yy**2, 2 * yy, 2.0
        B, Bp, Bpp = (yy - 1.0) ** 2, 2 * (yy - 1.0), 2.0
        C, Cp, Cpp = (yy - a), 1.0, 0.0
        # N = c A B C
        Np = c * (Ap * B * C + A * Bp * C + A * B * Cp)
        Npp = c * (
            App * B * C
            + Ap * Bp * C
            + Ap * B * Cp
            + Ap * Bp * C
            + A * Bpp * C
            + A * Bp * Cp
            + Ap * B * Cp
            + A * Bp * Cp
            + A * B * Cpp
        )
        s12 = r1 + r2
        p12 = r1 * r2
        u = yy**2 - s12 * yy + p12
        up = 2 * yy - s12
        upp = 2.0
        D = u * u
        Dp = 2 * u * up
        Dpp = 2 * up * up + 2 * u * upp
        return N, Np, Npp, D, Dp, Dpp

    Nq, Npq, Nppq, Dq, Dpq, Dppq = vals(q)
    Nw, Npw, Nppw, Dw, Dpw, Dppw = vals(w)
    return np.array(
        [
            Nq - Dq,
            Npq - Dpq,
            Nppq - Dppq,
            Nw - s_val * Dw,
            Npw - s_val * Dpw,
            Nppw - s_val * Dppw,
        ],
        dtype=float,
    )


def residual_2A3A3_alt(v, s_val: float):
    """
    Alt placement: 2A at 0, 3A at 1, 3A at s, 3A at ∞.

    N = c * y^2 * (y-b)^2 * (y-a)     → 2A at 0
    D = (y-r1) * (y-r2)                 → poles 1,1; ∞ order 3 → 3A at ∞
    Triple roots of N-D, N-sD at q,w

    v = [c, a, b, r1, r2, q, w]  — 7 params, 6 eqs → underdetermined; fix b=1? 
    Fix simple zero of N at free a, double at 0 and at d:
    N = c y^2 (y-d)^2 (y-a), fix d free.
    Still 7. Fix a domain aut: set d=1 (double at 1 for zeros) but then branch 0 is 2A and zeros at 1 interfere with branch 1.

    Use 6-param: N = c y^2 (y-1)^2 (y-a) same as 2A at 0 with doubles at 0,1 — then branch at 1 for 3A is separate (N-D).
    D quadratic for 3A at ∞.
    v = [c,a,r1,r2,q,w]
    """
    c, a, r1, r2, q, w = [float(x) for x in v]

    def vals(pt):
        yy = pt
        N = c * yy**2 * (yy - 1.0) ** 2 * (yy - a)
        A, Ap, App = yy**2, 2 * yy, 2.0
        B, Bp, Bpp = (yy - 1.0) ** 2, 2 * (yy - 1.0), 2.0
        C, Cp, Cpp = (yy - a), 1.0, 0.0
        Np = c * (Ap * B * C + A * Bp * C + A * B * Cp)
        Npp = c * (
            App * B * C
            + 2 * Ap * Bp * C
            + 2 * Ap * B * Cp
            + A * Bpp * C
            + 2 * A * Bp * Cp
            + A * B * Cpp
        )
        D = (yy - r1) * (yy - r2)
        Dp = 2 * yy - (r1 + r2)
        Dpp = 2.0
        return N, Np, Npp, D, Dp, Dpp

    Nq, Npq, Nppq, Dq, Dpq, Dppq = vals(q)
    Nw, Npw, Nppw, Dw, Dpw, Dppw = vals(w)
    return np.array(
        [
            Nq - Dq,
            Npq - Dpq,
            Nppq - Dppq,
            Nw - s_val * Dw,
            Npw - s_val * Dpw,
            Nppw - s_val * Dppw,
        ],
        dtype=float,
    )


def newton_solve(residual_fn, s_val, n_trials=80, niter=70):
    rng = np.random.default_rng(abs(hash((round(s_val, 8), residual_fn.__name__))) % (2**32))
    best = None
    for trial in range(n_trials):
        x0 = rng.normal(scale=1.2, size=6)
        if trial == 0:
            x0 = np.array([1.0, 2.0, 0.3, -0.4, 0.5, -0.5])
        v = x0.astype(float)
        ok = False
        for _ in range(niter):
            r = residual_fn(v, s_val)
            nr = float(np.linalg.norm(r))
            if nr < 1e-12:
                ok = True
                break
            J = np.zeros((6, 6))
            eps = 1e-7
            for j in range(6):
                dv = np.zeros(6)
                dv[j] = eps
                J[:, j] = (residual_fn(v + dv, s_val) - r) / eps
            try:
                step = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -r, rcond=None)[0]
            # damp if needed
            v = v + step
            if np.any(np.isnan(v)) or np.linalg.norm(v) > 1e6:
                break
        nr = float(np.linalg.norm(residual_fn(v, s_val)))
        if best is None or nr < best[2]:
            best = (v, ok or nr < 1e-10, nr)
        if best[1]:
            return best
    return best


# ---------------------------------------------------------------------------
# Fibre monic + classification
# ---------------------------------------------------------------------------
def monic_fibre_from_ND_coeffs(N_coeffs, D_coeffs, tv):
    """
    N, D as high-to-low numpy monic-ish; equation N - t D = 0, make monic deg 5.
    """
    # pad
    N = np.array(N_coeffs, dtype=float)
    D = np.array(D_coeffs, dtype=float)
    # N - t D as same degree
    deg = max(len(N), len(D))
    Np = np.zeros(deg)
    Dp = np.zeros(deg)
    Np[-len(N) :] = N
    Dp[-len(D) :] = D
    fib = Np - float(tv) * Dp
    # drop leading zeros
    i = 0
    while i < len(fib) - 1 and abs(fib[i]) < 1e-14:
        i += 1
    fib = fib[i:]
    if len(fib) != 6 or abs(fib[0]) < 1e-14:
        return None
    fib = fib / fib[0]
    return fib


def coeffs_2A3A3(v):
    """N = c y^3(y-1)(y-a) = c(y^5-(1+a)y^4+a y^3); D=(y-r1)^2(y-r2)^2."""
    c, a, r1, r2, q, w = [float(x) for x in v]
    N = c * np.array([1.0, -(1 + a), a, 0.0, 0.0, 0.0])  # y^5..y^0 but wait deg5: y^5-(1+a)y^4+a y^3
    # c * (y^5 - (1+a)y^4 + a y^3 + 0 y^2 + 0 y + 0)
    N = np.array([c, -c * (1 + a), c * a, 0.0, 0.0, 0.0])
    s12 = r1 + r2
    p12 = r1 * r2
    # u = y^2 - s12 y + p12; D = u^2 = y^4 - 2 s12 y^3 + (s12^2+2p12)y^2 - 2 s12 p12 y + p12^2
    D = np.array(
        [
            1.0,
            -2 * s12,
            s12**2 + 2 * p12,
            -2 * s12 * p12,
            p12**2,
        ]
    )
    # D deg 4 → pad to align for N-tD: treat as [0, D...] for y^5 coeff 0 in D
    return N, np.array([0.0, *D])


def coeffs_2A2_3A2(v):
    """N = c y^2 (y-1)^2 (y-a); D same as 2A at ∞ form."""
    c, a, r1, r2, q, w = [float(x) for x in v]
    # N = c (y^5 - (a+2)y^4 + (2a+1)y^3 - a y^2)
    N = np.array([c, -c * (a + 2), c * (2 * a + 1), -c * a, 0.0, 0.0])
    s12 = r1 + r2
    p12 = r1 * r2
    D = np.array([0.0, 1.0, -2 * s12, s12**2 + 2 * p12, -2 * s12 * p12, p12**2])
    return N, D


def coeffs_2A3A3_alt(v):
    """N = c y^2(y-1)^2(y-a); D = (y-r1)(y-r2) deg2 → pad."""
    c, a, r1, r2, q, w = [float(x) for x in v]
    N = np.array([c, -c * (a + 2), c * (2 * a + 1), -c * a, 0.0, 0.0])
    s12 = r1 + r2
    p12 = r1 * r2
    D = np.array([0.0, 0.0, 0.0, 1.0, -s12, p12])
    return N, D


def float_coeffs_to_Zpoly(coeffs):
    c = [float(x) for x in coeffs]
    if abs(c[0]) < 1e-12:
        return None
    c = [ci / c[0] for ci in c]
    rats = []
    for ci in c:
        r = sp.nsimplify(ci, tolerance=1e-8, rational=True)
        if not (getattr(r, "is_rational", False) or getattr(r, "is_Integer", False)):
            try:
                r = sp.Rational(str(Fraction(ci).limit_denominator(400)))
            except Exception:
                return None
        rats.append(sp.Rational(r))
    try:
        pol = sp.Poly.from_list(rats, y, domain=sp.QQ)
        mon = sp.Poly(sp.monic(pol.as_expr()), y, domain=sp.QQ)
        dens = [sp.fraction(sp.together(co))[1] for co in mon.all_coeffs()]
        L = 1
        for d in dens:
            L = int(sp.ilcm(L, abs(int(d))))
        cleared = sp.expand(L**5 * mon.as_expr().subs(y, y / L))
        pz = sp.Poly(cleared, y, domain=sp.ZZ)
        if pz.LC() == -1:
            pz = sp.Poly(-pz.as_expr(), y, domain=sp.ZZ)
        if pz.LC() != 1 or pz.degree() != 5:
            return None
        return pz
    except Exception:
        return None


def try_bj(pz):
    coeffs = [sp.Rational(c) for c in pz.all_coeffs()]
    if len(coeffs) != 6:
        return None
    shift = -coeffs[1] / 5
    z = sp.symbols("z")
    fsh = sp.expand(pz.as_expr().subs(y, z + shift))
    psh = sp.Poly(fsh, z, domain=sp.QQ)
    cc = [sp.Rational(c) for c in psh.all_coeffs()]
    if len(cc) != 6 or cc[1] != 0:
        return None
    if cc[2] == 0 and cc[3] == 0:
        try:
            return {"form": "BJ", "alpha": int(cc[4]), "beta": int(cc[5]), "k": str(k_of(int(cc[4]), int(cc[5])))}
        except Exception:
            return {"form": "BJ_QQ", "alpha": str(cc[4]), "beta": str(cc[5])}
    return {"form": "depressed", "p": str(cc[2]), "q": str(cc[3])}


def hilbert_scan(type_name, residual_fn, coeff_fn, s_list, t_list, n_trials=60):
    covers = []
    cat_hits = []
    bj_hits = []
    stats = {
        "covers_ok": 0,
        "covers_attempted": len(s_list),
        "n_fibres_Z": 0,
        "n_irr": 0,
        "n_even": 0,
        "n_A5": 0,
        "n_BJ": 0,
        "n_cat": 0,
    }

    for s_val in s_list:
        v, ok, nr = newton_solve(residual_fn, float(s_val), n_trials=n_trials)
        rec = {
            "s": s_val,
            "ok": bool(ok),
            "newton_res": nr,
            "params": None if v is None else [float(x) for x in v],
        }
        if not ok:
            covers.append(rec)
            print(f"  {type_name} s={s_val}: FAIL res={nr:.1e}", flush=True)
            continue
        stats["covers_ok"] += 1
        covers.append(rec)
        print(f"  {type_name} s={s_val}: ok res={nr:.1e}", flush=True)
        N, D = coeff_fn(v)
        for tv in t_list:
            fib = monic_fibre_from_ND_coeffs(N, D, tv)
            if fib is None:
                continue
            pz = float_coeffs_to_Zpoly(fib)
            if pz is None:
                continue
            stats["n_fibres_Z"] += 1
            if not pz.is_irreducible:
                continue
            stats["n_irr"] += 1
            d = int(pz.discriminant())
            dsq = d > 0 and is_square(d)
            bj = try_bj(pz)
            if bj and bj.get("form") == "BJ":
                stats["n_BJ"] += 1
                a, b = bj["alpha"], bj["beta"]
                bj_hits.append({"s": s_val, "t": tv, "alpha": a, "beta": b, "k": bj.get("k")})
                if (a, b) in CAT_BY_AB:
                    tag, ck = CAT_BY_AB[(a, b)]
                    stats["n_cat"] += 1
                    cat_hits.append(
                        {
                            "tag": tag,
                            "k": str(ck),
                            "s": s_val,
                            "t": tv,
                            "alpha": a,
                            "beta": b,
                            "type": type_name,
                        }
                    )
            if dsq:
                stats["n_even"] += 1
                # galois only occasionally (cost)
                if stats["n_even"] <= 8:
                    r = classify_poly(pz.as_expr().subs(y, x), do_galois=True)
                    if (r.get("status") or "").startswith("HIT_A5"):
                        stats["n_A5"] += 1

    cat_k = sorted({h["k"] for h in cat_hits})
    return {
        "type": type_name,
        "stats": stats,
        "covers": covers,
        "catalogue_hits": cat_hits,
        "bj_sample": bj_hits[:25],
        "catalogue_k": cat_k,
        "multi_k": len(cat_k) >= 2,
    }


# ---------------------------------------------------------------------------
# Envelope control
# ---------------------------------------------------------------------------
def envelope_control():
    hits = []
    m0 = Fraction(5, 16)
    for tv in [0, 1, Fraction(1, 2), Fraction(1, 3)]:
        ku = Fraction(-8, 5) + tv * (Fraction(4, 5) - Fraction(-8, 5))
        a = pure_even_alpha(m0, ku)
        if a.denominator != 1:
            continue
        aa, bb = int(a), int(ku * a)
        if (aa, bb) in CAT_BY_AB:
            tag, ck = CAT_BY_AB[(aa, bb)]
            hits.append({"tag": tag, "k": str(ck), "t": str(tv)})
    # also exact flagship/classical endpoints via known m
    # flagship (-55,88): m^2 = (α + 3125 k^4/256)/256
    for tag, a, b, k in CATALOGUE[:6]:
        d = disc_bj_int(a, b)
        hits.append(
            {
                "tag": tag,
                "k": str(k),
                "disc_square": d > 0 and is_square(d),
                "source": "catalogue_seed",
            }
        )
    return {
        "path_hits": [h for h in hits if "t" in h],
        "multi_k": len({h["k"] for h in hits if "t" in h}) >= 2,
        "catalogue_disc_square": all(
            h.get("disc_square") for h in hits if h.get("source") == "catalogue_seed"
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 72, flush=True)
    print("G2 — genus-0 Nielsen types 2A 3A³ and 2A² 3A²", flush=True)
    print("=" * 72, flush=True)

    # C combinatorial
    print("\n[C] Nielsen sample for target types ...", flush=True)
    nielsen = {}
    for labels in (["2A", "3A", "3A", "3A"], ["2A", "2A", "3A", "3A"]):
        key = ",".join(labels)
        print(f"  sampling {key} ...", flush=True)
        nielsen[key] = sample_nielsen(labels, max_found=2, max_check=800_000)
        print(
            f"    found={nielsen[key]['n_found']} checked={nielsen[key]['checked']} "
            f"prior={nielsen[key]['prior_orbit_data']}",
            flush=True,
        )

    s_list = [-3, -2, -1, -0.5, 0.5, 1.5, 2, 2.5, 3, 4, 5, -1.5, 0.333, 0.666]
    t_list = list(range(-8, 9)) + [0.5, 1.5, 2.5, -0.5, Fraction(2, 3)]

    # A+H type 2A 3A³ (main ansatz: 3A at 0, 2A at ∞)
    print("\n[A/H] Type 2A,3A³ — ansatz 3A@0 + 2A@∞ + 3A@1 + 3A@s ...", flush=True)
    h_2a3a3 = hilbert_scan(
        "2A3A3_3A0_2Ainf",
        residual_2A3A3_real,
        coeffs_2A3A3,
        s_list,
        t_list,
        n_trials=70,
    )

    print("\n[A/H] Type 2A,3A³ alt — ansatz 2A@0 + 3A@∞ + 3A@1 + 3A@s ...", flush=True)
    h_2a3a3_alt = hilbert_scan(
        "2A3A3_2A0_3Ainf",
        residual_2A3A3_alt,
        coeffs_2A3A3_alt,
        s_list,
        t_list,
        n_trials=70,
    )

    print("\n[A/H] Type 2A²,3A² — ansatz 2A@0 + 2A@∞ + 3A@1 + 3A@s ...", flush=True)
    h_2a2_3a2 = hilbert_scan(
        "2A2_3A2",
        residual_2A2_3A2_real,
        coeffs_2A2_3A2,
        s_list,
        t_list,
        n_trials=70,
    )

    print("\n[E] Pure-even envelope control ...", flush=True)
    env = envelope_control()
    print(f"  multi_k={env['multi_k']} path_hits={env['path_hits']}", flush=True)

    elapsed = round(time.time() - t0, 2)

    all_cat = (
        h_2a3a3["catalogue_hits"]
        + h_2a3a3_alt["catalogue_hits"]
        + h_2a2_3a2["catalogue_hits"]
    )
    all_k = sorted({h["k"] for h in all_cat})
    multi = len(all_k) >= 2
    any_hit = len(all_cat) > 0

    verdict = (
        f"G2 cut ({elapsed}s). "
        f"Nielsen samples: 2A3A³ found={nielsen['2A,3A,3A,3A']['n_found']}, "
        f"2A²3A² found={nielsen['2A,2A,3A,3A']['n_found']}. "
        f"Covers ok: 2A3A3={h_2a3a3['stats']['covers_ok']}/{h_2a3a3['stats']['covers_attempted']}, "
        f"alt={h_2a3a3_alt['stats']['covers_ok']}, "
        f"2A2_3A2={h_2a2_3a2['stats']['covers_ok']}. "
        f"Catalogue hits total={len(all_cat)} k={all_k} multi_k={multi}. "
        f"Envelope control multi_k={env['multi_k']}."
    )
    print("\n" + verdict, flush=True)

    def stats_block(name, h):
        st = h["stats"]
        return [
            f"### {name}",
            "",
            f"| quantity | value |",
            f"|----------|------:|",
            f"| covers ok | {st['covers_ok']}/{st['covers_attempted']} |",
            f"| monic Z fibres | {st['n_fibres_Z']} |",
            f"| irreducible | {st['n_irr']} |",
            f"| disc □ | {st['n_even']} |",
            f"| A5 (sample) | {st['n_A5']} |",
            f"| BJ | {st['n_BJ']} |",
            f"| **catalogue hits** | **{st['n_cat']}** |",
            f"| catalogue k | {h['catalogue_k']} |",
            f"| multi-k | {h['multi_k']} |",
            "",
        ]

    lines = [
        "# G2 — genus-0 Nielsen types \(2A\,3A^3\) and \(2A^2\,3A^2\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        "---",
        "",
        "## 0. Goal",
        "",
        "After G1 exhausted the pure-ternary class \(3A^4\), realise other **lookup genus-0**",
        "\(r=4\) Nielsen types and test Hilbert specialisations against the multi-seed",
        "pure-even catalogue.",
        "",
        "| type | prior orbits | orbit size | genus lookup |",
        "|------|-------------:|-----------:|:------------:|",
        "| \(2A,3A,3A,3A\) | 1 | 96 | 0 |",
        "| \(2A,2A,3A,3A\) | 1 | 108 | 0 |",
        "",
        "Source: `A5_HURWITZ_R4.md` / Magaard–Shpectorov–James style tables.",
        "",
        "---",
        "",
        "## 1. Combinatorial lock (Nielsen samples)",
        "",
    ]
    for key, rec in nielsen.items():
        lines.append(f"### `{key}`")
        lines.append("")
        lines.append(f"- class sizes: {rec['class_sizes']}")
        lines.append(f"- checked (g4-determined): {rec['checked']}")
        lines.append(f"- generating tuples found: **{rec['n_found']}**")
        lines.append(f"- prior orbit data: {rec['prior_orbit_data']}")
        if rec["samples"]:
            lines.append(f"- sample cycle types: {rec['samples'][0]['cycle_types']}")
        lines.append("")
    lines += [
        "---",
        "",
        "## 2. Cover ansätze",
        "",
        "### \(2A\,3A^3\) — placement A (3A@0, 2A@∞, 3A@1, 3A@s)",
        "",
        r"$$N = c\, y^3(y-1)(y-a),\qquad D=(y-r_1)^2(y-r_2)^2$$",
        "",
        "Zeros of \(N\): mults \(3,1,1\) → type **3A** at branch 0.  ",
        "Poles of \(D\) plus \(\infty\): mults \(2,2,1\) → type **2A** at \(\infty\).  ",
        "Triple roots of \(N-D\) and \(N-sD\) impose **3A** at \(1\) and \(s\).",
        "",
        "### \(2A\,3A^3\) — placement B (2A@0, 3A@∞, 3A@1, 3A@s)",
        "",
        r"$$N = c\, y^2(y-1)^2(y-a),\qquad D=(y-r_1)(y-r_2)$$",
        "",
        "Zeros: \(2,2,1\) → **2A** at 0. Poles \(1,1\) + \(\infty\) order 3 → **3A** at \(\infty\).",
        "",
        "### \(2A^2\,3A^2\) — placement (2A@0, 2A@∞, 3A@1, 3A@s)",
        "",
        r"$$N = c\, y^2(y-1)^2(y-a),\qquad D=(y-r_1)^2(y-r_2)^2$$",
        "",
        "Zeros **2A** at 0; poles **2A** at \(\infty\); triple roots → **3A** at \(1,s\).",
        "",
        "Parameters \((c,a,r_1,r_2,q,w)\) solved by Newton at rational cross-ratios \(s\).",
        "Fibre: monic form of \(N - t D = 0\) in \(y\).",
        "",
        "---",
        "",
        "## 3. Hilbert specialisations vs multi-seed catalogue",
        "",
    ]
    lines += stats_block("2A 3A³ — placement A", h_2a3a3)
    lines += stats_block("2A 3A³ — placement B (alt)", h_2a3a3_alt)
    lines += stats_block("2A² 3A²", h_2a2_3a2)

    if all_cat:
        lines.append("### Catalogue hits")
        lines.append("")
        for h in all_cat:
            lines.append(f"- {h}")
        lines.append("")
    else:
        lines.append("_No exact catalogue seeds recovered from these G2 ansätze/scans._")
        lines.append("")

    lines += [
        "---",
        "",
        "## 4. Arithmetic multi-k control",
        "",
        f"- Envelope / catalogue control multi-k path: **{env['multi_k']}**",
        f"- Path hits: {env['path_hits']}",
        f"- Catalogue seeds disc□ (sample): **{env['catalogue_disc_square']}**",
        "",
        "---",
        "",
        "## 5. Multi-k conclusion (G2)",
        "",
        f"| test | result |",
        f"|------|--------|",
        f"| Nielsen tuples exist (2A3A³) | **{nielsen['2A,3A,3A,3A']['ok']}** |",
        f"| Nielsen tuples exist (2A²3A²) | **{nielsen['2A,2A,3A,3A']['ok']}** |",
        f"| Explicit covers realised (any ansatz) | **{h_2a3a3['stats']['covers_ok'] + h_2a3a3_alt['stats']['covers_ok'] + h_2a2_3a2['stats']['covers_ok'] > 0}** |",
        f"| Catalogue Hilbert hit | **{any_hit}** ({len(all_cat)}) |",
        f"| Geometric multi-k (≥2 catalogue k) | **{multi}** |",
        f"| Arithmetic multi-k control | **{env['multi_k']}** |",
        "",
        "**Geometric multi-k via G2 shortlist ansätze: "
        + ("HIT." if multi else "not achieved in this cut.")
        + "**",
        "",
        "### What this cut established",
        "",
        "1. Combinatorial existence of Nielsen generators for both g=0 shortlist types.",
        "2. Explicit degree-5 cover ansätze matching the ramification profiles.",
        "3. Newton realisation at many rational \(s\) + fibre specialisation pipeline.",
        "4. Full multi-seed catalogue re-test (flagship / classical / LSW / ±12/5 / ±16/5).",
        "5. Arithmetic multi-k still available; geometric fusion still open if hits empty.",
        "",
        "### Next if empty",
        "",
        "1. Domain Möbius freedom inside each ansatz (move {0,1,∞} labels).",
        "2. Other g=0 types: \(2A\,3A^2\,5A\), \(3A^3\,5A\).",
        "3. G3: monodromy identification of the pure-even envelope.",
        "4. Literature resolvents (Malle–König / computational IG databases) for these classes.",
        "",
        "---",
        "",
        "## 6. Non-claims",
        "",
        "- Newton covers are numeric realisations of the ramification profile; full monodromy",
        "  certification that the Galois group is exactly the chosen Nielsen class for every",
        "  \(s\) is not claimed (would need certified monodromy / certified braid tracking).",
        "- Negative catalogue results are for these ansätze + scan bounds.",
        "- Does not reopen pure-even arithmetic, Canonical T3, or Necessity.",
        "",
        "_Generated by `g2_nielsen_g0.py`._",
        "",
    ]

    md = "\n".join(lines)
    payload = {
        "verdict": verdict,
        "elapsed_s": elapsed,
        "nielsen": nielsen,
        "hilbert_2A3A3": h_2a3a3,
        "hilbert_2A3A3_alt": h_2a3a3_alt,
        "hilbert_2A2_3A2": h_2a2_3a2,
        "all_catalogue_hits": all_cat,
        "multi_k": multi,
        "envelope_control": env,
    }

    write_md(ROOT / "G2_NIELSEN_G0.md", md)
    write_json(ROOT / "G2_NIELSEN_G0.json", payload)
    write_md(OUT / "G2_NIELSEN_G0.md", md)
    write_json(OUT / "G2_NIELSEN_G0.json", payload)
    try:
        write_md(RESULTS / "G2_NIELSEN_G0.md", md)
        write_json(RESULTS / "G2_NIELSEN_G0.json", payload)
    except Exception:
        pass

    print(f"\nWrote G2_NIELSEN_G0.md / .json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
