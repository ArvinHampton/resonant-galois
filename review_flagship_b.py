"""
Automated review checks for REVIEW_PACKAGE.md order:
  1. Flagship Mestre P_t
  2. B-embed identity + T-match + sample Gal
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import classify_poly, is_square, x  # noqa: E402

t, y, z = sp.symbols("t y z")
A = sp.symbols("A")


def check_flagship() -> None:
    P = x**5 - 55 * x + 88
    R = x**4 + 8 * x**3 - 32 * x**2 + 33
    # C1 Mestre
    W = sp.diff(P, x, 2) * R - 2 * sp.diff(P, x) * sp.diff(R, x)
    rem = sp.Poly(sp.expand(W), x).rem(sp.Poly(P, x))
    assert rem == 0, f"Mestre rem={rem}"
    # C2 family
    F = sp.expand(sp.resultant(P.subs(x, y), z - y - t * R.subs(x, y), y))
    pol = sp.Poly(F, z)
    assert pol.LC() == 1 and pol.degree() == 5
    # C3 t=0
    assert sp.expand(F.subs(t, 0) - (z**5 - 55 * z + 88)) == 0
    # C4 disc square
    D = sp.expand(pol.discriminant())
    cont, facs = sp.factor_list(D)
    assert all(m % 2 == 0 or not getattr(fi, "free_symbols", set()) for fi, m in facs)
    c = sp.Rational(cont)
    assert c > 0 and sp.integer_nthroot(int(c.p), 2)[1] and sp.integer_nthroot(int(c.q), 2)[1]
    # closed form coeffs: factored + expanded table
    c4, c3, c2, c1, c0 = [sp.expand(pol.coeff_monomial(z**k)) for k in (4, 3, 2, 1, 0)]
    assert c4 == -385 * t
    assert c3 == -167200 * t**2 - 1320 * t
    assert sp.expand(c3 - (-440 * t * (380 * t + 3))) == 0
    assert c2 == 63888000 * t**3 + 721600 * t**2 + 10560 * t
    assert sp.expand(c2 - 3520 * t * (18150 * t**2 + 205 * t + 3)) == 0
    assert c1 == 2509056000 * t**4 - 240064000 * t**3 - 1161600 * t**2 - 14080 * t - 55
    assert c0 == (
        24966656000 * t**5
        - 4894208000 * t**4
        + 232320000 * t**3
        + 774400 * t**2
        + 1815 * t
        + 88
    )
    # C5 Gal samples
    for tv in [0, 1, -1, 2, 3, 5, 9, 27, 61, 80]:
        chi = sp.expand(F.subs(t, tv)).subs(z, x)
        rec = classify_poly(chi, do_galois=True)
        assert rec.get("disc_square"), f"t={tv} disc not square"
        assert str(rec.get("status", "")).startswith("HIT_A5"), f"t={tv} {rec.get('status')}"
        print(f"  flagship t={tv}: {rec.get('status')}", flush=True)
    print("flagship P_t: PASS", flush=True)


def check_B() -> None:
    P = x**5 + 75 * x**3 + A * x**2 + 3 * A
    D = sp.expand(sp.Poly(P, x).discriminant())
    assert sp.expand(D - (18 * A * (A**2 + 84375)) ** 2) == 0
    b = sp.symbols("b", nonzero=True)
    chiT = x**5 - (-75) * x**3 - (-A) * x**2 + ((-A) * (-75) - b * (72 * A / b))
    assert sp.simplify(chiT - P) == 0
    for Av in [3, 9, 61, 80, 539, -3]:
        rec = classify_poly(sp.expand(P.subs(A, Av)), do_galois=True)
        assert rec.get("disc_square"), f"A={Av}"
        print(f"  B A={Av}: disc□={rec.get('disc_square')} {rec.get('status')}", flush=True)
    print("B-embed / avatars: PASS", flush=True)


def main():
    print("REVIEW PACKAGE CHECKS (flagship first)", flush=True)
    check_flagship()
    check_B()
    print("ALL REVIEW CHECKS PASSED", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
