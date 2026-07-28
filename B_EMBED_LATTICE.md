# B-embed lattice — review object (secondary)

**Status:** Generative centre extension · beyond BJ-embed · **not** necessity  
**Review priority:** **2 — after flagship \(P_t\)**  
**Companion:** `EVENNESS_AVATAR.md` (matrix packaging of the same identity)

---

## Theorem-shaped statement

### Polynomial family

For a parameter \(A\in\mathbb{Q}\setminus\{0\}\),
\[
P_A(x)=x^5+75x^3+A x^2+3A\in\mathbb{Q}[x].
\]

**Discriminant identity (verified):**
\[
\operatorname{disc}(P_A)
=
324\,A^2\,(A^2+84375)^2
=
\bigl(18\,A\,(A^2+84375)\bigr)^2.
\]
Hence \(\operatorname{disc}(P_A)\) is a square in \(\mathbb{Q}[A]\) for all \(A\).

### Matrix realisation (beyond BJ-embed)

The structural template
\[
\chi_T=x^5-d x^3-(a+ef)x^2-(bf+ce)x+(ad-bc)
\]
matches \(P_A\) under
\[
d=-75,\qquad e=f=0,\qquad a=-A,\qquad bc=72A.
\]
Explicitly, for any \(b\mid 72A\),
\[
T(-A,\,b,\,72A/b,\,-75,\,0,\,0)
\quad\Longrightarrow\quad
\chi_T=P_A.
\]
Here \(d=-75\neq 0\), so the realisation is **not** contained in the BJ-embed \(d=0\), \(a=-ef\).

---

## Claims for review (checklist)

| # | Claim | Status |
|---|--------|:------:|
| B1 | \(\operatorname{disc}(P_A)=(18A(A^2+84375))^2\) in \(\mathbb{Q}[A]\) | **Verified** |
| B2 | \(\chi_T=P_A\) under \(d=-75,e=f=0,a=-A,bc=72A\) | **Verified** |
| B3 | Beyond BJ-embed (\(d\neq 0\)) | **True** by construction |
| B4 | Lattice \(A\in L_0\) yields many irr + \(\mathrm{Gal}=A_5\) specialisations | **Sample table** |

---

## Lattice \(A\to A_5\) (sample)

Prefer \(b\in\{1,3,8,9,\ldots\}\) dividing \(72A\). Representative resonant \(A\):

| \(A\) | example \(T\) | disc □ | sample Gal |
|------:|---------------|:------:|------------|
| \(3\) | \(T(-3,8,27,-75,0,0)\) | yes | \(A_5\) |
| \(9\) | \(T(-9,3,216,-75,0,0)\) | yes | \(A_5\) |
| \(27\) | \(T(-27,3,648,-75,0,0)\) | yes | \(A_5\) |
| \(61\) | \(T(-61,3,1464,-75,0,0)\) | yes | \(A_5\) |
| \(80\) | \(T(-80,3,1920,-75,0,0)\) | yes | \(A_5\) |
| \(243\) | \(T(-243,3,5832,-75,0,0)\) | yes | \(A_5\) |
| \(539\) | \(T(-539,3,12936,-75,0,0)\) | yes | \(A_5\) |
| \(55\) | \(T(-55,3,1320,-75,0,0)\) | yes | \(A_5\) |
| \(88\) | \(T(-88,3,2112,-75,0,0)\) | yes | \(A_5\) |
| \(4880\) | \(T(-4880,3,117120,-75,0,0)\) | yes | \(A_5\) |

Sign flips \(A\mapsto -A\) (and corresponding \(b,c\)) behave similarly.  
Broader scan: **104** unique lattice \(A\), all disc□ by identity; **50** Gal-checked \(A_5\) in `l0_secondary_invariants` / prior B-embed run.

**Residual arithmetic:** irreducibility (and \(A_5\) vs smaller even groups) is not free — only evenness is identical.

**Contamination:** \(A\in\mathbb{Z}\) (e.g. \(61,80,539\)) is a specialisation parameter only. No physical period or experimental “hit” is used in the disc identity or embed. See `CONTAMINATION_BOUNDARY.md`.

---

## What this is / is not

| Is | Is not |
|----|--------|
| Second evenness avatar (parameter \(A\), \(d\neq 0\)) | Crit-2 forcing on unrestricted \(T\) |
| Systematic lattice specialisations | HQCC-axiom naming of \(d=-75\) or \(bc=72A\) |
| Parallel track to pure-even / Mestre | Necessity theorem |

---

## Minimal verification script

```python
import sympy as sp
from lib.common import classify_poly, x

A = sp.symbols("A")
P = x**5 + 75*x**3 + A*x**2 + 3*A
D = sp.expand(sp.Poly(P, x).discriminant())
assert sp.expand(D - (18*A*(A**2 + 84375))**2) == 0
b = sp.symbols("b", nonzero=True)
chiT = (
    x**5 - (-75)*x**3 - (-A)*x**2
    + ((-A)*(-75) - b*(72*A/b))
)
assert sp.simplify(chiT - P) == 0
for Av in [3, 61, 80, 539]:
    rec = classify_poly(P.subs(A, Av), do_galois=True)
    assert rec["disc_square"]
print("B-embed: PASS")
```

---

## Cross-links

- Primary review: `MESTRE_FLAGSHIP_PT.md`  
- Avatar packaging: `EVENNESS_AVATAR.md`  
- Review order: `REVIEW_PACKAGE.md`

_Review-grade rewrite; identities re-verified 2026-07-24._
