# Flagship Mestre lift — review object (primary)

**Status:** Generative centre extension · **not** a necessity theorem  
**Review priority:** **1 — first**  
**Script check:** `python review_flagship_b.py` (or re-run block below)

---

## Theorem-shaped statement

Let
\[
P(x)=x^5-55x+88\in\mathbb{Z}[x]
\]
(the HQCC flagship seed: \(88=61+3^3\), \(\mathrm{Gal}(P/\mathbb{Q})=A_5\), \(\operatorname{disc}(P)\) a square).

Let
\[
R(x)=x^4+8x^3-32x^2+33.
\]

**Mestre condition (verified):**
\[
P''R-2P'R'\equiv 0\pmod{P}
\]
(the \(R\)-space is 1-dimensional over \(\mathbb{Q}\); this \(R\) is a generator).

**Family:**
\[
P_t(z)
=\operatorname{Res}_y\bigl(P(y),\; z-y-t\,R(y)\bigr)
\in\mathbb{Q}[t][z],
\]
monic of degree 5 in \(z\).

---

## Explicit coefficient table

Write
\[
P_t(z)=z^5+c_4(t)\,z^4+c_3(t)\,z^3+c_2(t)\,z^2+c_1(t)\,z+c_0(t).
\]

| Power | Coefficient \(c_i(t)\) |
|-------|------------------------|
| \(z^5\) | \(1\) |
| \(z^4\) | \(-385\,t\) |
| \(z^3\) | \(-440\,t\,(380\,t+3)=-167200\,t^2-1320\,t\) |
| \(z^2\) | \(3520\,t\,(18150\,t^2+205\,t+3)=63888000\,t^3+721600\,t^2+10560\,t\) |
| \(z^1\) | \(55\,(45619200\,t^4-4364800\,t^3-21120\,t^2-256\,t-1)\) \(=2509056000\,t^4-240064000\,t^3-1161600\,t^2-14080\,t-55\) |
| \(z^0\) | \(11\,(2269696000\,t^5-444928000\,t^4+21120000\,t^3+70400\,t^2+165\,t+8)\) \(=24966656000\,t^5-4894208000\,t^4+232320000\,t^3+774400\,t^2+1815\,t+88\) |

**Note on \(c_2\):** expanded leading term is \(63888000\,t^3\) (\(3520\times 18150\)), not \(63936000\).

### Compact factored form (locked)

\[
\begin{align*}
c_4&=-385\,t\\
c_3&=-440\,t\,(380\,t+3)\\
c_2&=3520\,t\,(18150\,t^2+205\,t+3)\\
c_1&=55\,(45619200\,t^4-4364800\,t^3-21120\,t^2-256\,t-1)\\
c_0&=11\,(2269696000\,t^5-444928000\,t^4+21120000\,t^3+70400\,t^2+165\,t+8)
\end{align*}
\]

### Sanity checks

| \(t\) | Leading form of \(P_t\) | Note |
|------:|-------------------------|------|
| \(0\) | \(z^5-55z+88\) | Recovers seed |
| \(1\) | substitute \(t=1\) in the table | Sample Gal \(A_5\) (verified) |

Disc identity (auto-check): \(\operatorname{disc}_z(P_t)=58564000000\cdot Q(t)^2\) in \(\mathbb{Q}[t]\) for an explicit \(Q\in\mathbb{Q}[t]\) of degree 10 (see below).

---

## Claims for review (checklist)

| # | Claim | Status |
|---|--------|:------:|
| C1 | \(P''R-2P'R'\equiv 0\pmod{P}\) | **Verified** (remainder \(0\)) |
| C2 | \(P_t\) monic deg 5, LC \(=1\) | **Verified** |
| C3 | \(P_0=P\) (recovery of seed) | **Verified** |
| C4 | \(\operatorname{disc}_z(P_t)\) is a square in \(\mathbb{Q}[t]\) | **Verified** |
| C5 | Sample specialisations: irr + disc□ + \(\mathrm{Gal}=A_5\) | **Verified** (table) |

### Discriminant identity (C4)

\[
\operatorname{disc}_z(P_t)
=
58564000000
\cdot
Q(t)^2,
\]
where
\[
\begin{aligned}
Q(t)
&=
20545741772554240000\,t^{10}
-2221801010626560000\,t^9
+38099706839040000\,t^8 \\
&\quad
+763228815360000\,t^7
-9746366080000\,t^6
-1882566400\,t^5 \\
&\quad
+1018662425\,t^4
-3953400\,t^3
+1.
\end{aligned}
\]

- Content \(58564000000=242000^2\) is a square.  
- No odd-multiplicity non-constant factors.  
- Hence \(\operatorname{disc}_z(P_t)\) is a square in \(\mathbb{Q}[t]\).

*(Note: \(58564000000=\operatorname{disc}(P)\), consistent with \(t=0\).)*

---

## Sample Galois table (C5)

Integer specialisations \(t\in\mathbb{Z}\) of the monic model above (already in \(\mathbb{Z}[z]\) for these \(t\)):

| \(t\) | disc □ | irreducible | Gal / status |
|------:|:------:|:-----------:|--------------|
| \(0\) | yes | yes | \(A_5\) (seed) |
| \(1\) | yes | yes | \(A_5\) |
| \(-1\) | yes | yes | \(A_5\) |
| \(2\) | yes | yes | \(A_5\) |
| \(3\) | yes | yes | \(A_5\) |
| \(5\) | yes | yes | \(A_5\) |
| \(7\) | yes | yes | \(A_5\) |
| \(9\) | yes | yes | \(A_5\) |
| \(27\) | yes | yes | \(A_5\) |
| \(61\) | yes | yes | \(A_5\) |
| \(80\) | yes | yes | \(A_5\) |

Lattice parameters \(t\in\{3,9,27,61,80\}\) appear in the table and remain \(A_5\).

---

## What this is / is not

| Is | Is not |
|----|--------|
| Explicit 1-param \(A_5\) family through the HQCC flagship seed | Forced alternating monodromy from HQCC axioms alone |
| Disc □ as a polynomial identity in \(t\) | A cut of template \(T\) forcing evenness |
| Generative enlargement of the pure-even centre | A necessity theorem |

**Dependence:** one even seed + classical Mestre deformation. HQCC enters by **choice of seed** \(P\), not by naming the differential condition.

**Contamination:** \(P\) and any specialisation \(t\in\mathbb{Z}\) (including lattice values such as \(61,80,539\)) are used only as **integers**. No physical period, detector claim, or 539-step dynamical law enters the resultant or disc identity. See `CONTAMINATION_BOUNDARY.md`.

---

## Minimal verification script

```python
import sympy as sp
from lib.common import classify_poly, x

t, y, z = sp.symbols("t y z")
P = x**5 - 55*x + 88
R = x**4 + 8*x**3 - 32*x**2 + 33
assert sp.Poly(sp.diff(P,x,2)*R - 2*sp.diff(P,x)*sp.diff(R,x), x).rem(sp.Poly(P,x)) == 0
F = sp.resultant(P.subs(x,y), z - y - t*R.subs(x,y), y)
pol = sp.Poly(sp.expand(F), z)
assert pol.LC() == 1 and pol.degree() == 5
assert sp.expand(F.subs(t,0) - (z**5 - 55*z + 88)) == 0
D = sp.expand(pol.discriminant())
cont, facs = sp.factor_list(D)
assert all(m % 2 == 0 or not fi.free_symbols for fi, m in facs)
assert sp.integer_nthroot(abs(int(sp.Rational(cont))), 2)[1]
for tv in [0, 1, -1, 2, 3, 61, 80]:
    chi = sp.expand(F.subs(t, tv)).subs(z, x)
    rec = classify_poly(chi, do_galois=True)
    assert rec["disc_square"] and str(rec.get("status","")).startswith("HIT_A5")
print("flagship P_t: PASS")
```

---

## Cross-links

- Pure-even centre: `PURE_EVEN_MULTI_K.md`  
- Mestre on full lattice: `L0_MESTRE_ORBIT.md`  
- Second review object (beyond BJ): `B_EMBED_LATTICE.md`, `EVENNESS_AVATAR.md`  
- Review order: `REVIEW_PACKAGE.md`

_Review-grade rewrite; identities re-verified 2026-07-24._
