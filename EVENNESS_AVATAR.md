# Evenness avatars — review companion (with B-embed)

**Status:** Built-in disc□ by parameterisation · **not** Crit-2 forcing  
**Review priority:** **2b — with `B_EMBED_LATTICE.md`**  
**Do not review before:** flagship Mestre lift (`MESTRE_FLAGSHIP_PT.md`)

---

## Design rule

Admit only matrix models whose \(\operatorname{disc}(\chi)\) is a **square as a polynomial** in free parameters — not sparse search on \(T\).

---

## Avatar PE — pure-even BJ matrix

\[
T(0,-\alpha,k,0,0,1),
\qquad
\alpha=256m^2-\frac{3125\,k^4}{256},
\quad
\beta=k\alpha.
\]

| Claim | Status |
|-------|:------:|
| \(\chi=x^5+\alpha x+\beta\) | **Identity** |
| \(\operatorname{disc}=(256\alpha^2 m)^2\) | **Identity** |
| Beyond BJ-embed? | No (\(d=0\)) |
| HQCC-native forcing? | **No** — classical pure-even; lattice specialises \((m,k)\) |

This is the matrix packaging of the **finished pure-even multi-\(k\) centre** (`PURE_EVEN_MULTI_K.md`).

---

## Avatar B — non-BJ degree-1 matrix

\[
T(-A,\,b,\,72A/b,\,-75,\,0,\,0)
\qquad(b\mid 72A).
\]

| Claim | Status |
|-------|:------:|
| \(\chi=x^5+75x^3+A x^2+3A\) | **Identity** |
| \(\operatorname{disc}=(18A(A^2+84375))^2\) | **Identity** |
| Beyond BJ-embed (\(d\neq 0\))? | **Yes** |
| HQCC-native forcing? | **No** — \(d=-75\), \(bc=72A\) are classical ansatz; lattice specialises \(A\) |

Full lattice tables: **`B_EMBED_LATTICE.md`**.

---

## Side-by-side

| | PE avatar | B avatar |
|--|-----------|----------|
| Free params | \(m,k\) | \(A\) (and divisor \(b\)) |
| Beyond BJ | No | **Yes** |
| Evenness | pure-even formula | \(P_A\) disc identity |
| Role vs Mestre flagship | Classical centre | Parallel generative track |

---

## What a review should **not** claim

- That either avatar is forced by unrestricted ternary matrix axioms.  
- That \(L_0\) specialisation upgrades these to a necessity theorem.  
- That PE and B are related by a canonical map \(\Phi\) (open; see `L0_PE_B_UNIFY.md`).

---

## Minimal check

```python
import sympy as sp
from lib.common import x
m, k, A, b = sp.symbols("m k A b", nonzero=True)
al = 256*m**2 - sp.Rational(3125)*k**4/256
be = k*al
# PE charpoly via BJ-embed
chi_pe = x**5 + al*x + be
Dpe = 256*al**5 + 3125*be**4
assert sp.expand(Dpe - (256*al**2*m)**2) == 0
PB = x**5 + 75*x**3 + A*x**2 + 3*A
assert sp.expand(sp.Poly(PB,x).discriminant() - (18*A*(A**2+84375))**2) == 0
print("avatars: PASS")
```

_Review companion to B_EMBED_LATTICE.md; re-verified 2026-07-24._
