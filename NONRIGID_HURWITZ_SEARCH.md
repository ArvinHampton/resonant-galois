# Non-rigid pure-even multi-\(k\) search (abandon \(\varphi\))

_Elapsed: 119.19s_

**Verdict:** Arithmetic multi-k envelope: pure-even over Q(m,s), recovers 26/26 catalogue seeds on fixed-k slices — NOT Hurwitz-geometric. Cross-k rational curves (monomial/Bezier/quad) non-ray pure-even hits: 0. Mestre multi-seed pure-even: 0. Biparam poly multi-k: 0. Constructive envelope paths through cross-k seed pairs: 5/5 (arithmetic pure-even 1-param through two k's). Concrete positive-dim Hurwitz pure-even multi-k candidate: None. AMBITIOUS ROUTE: no geometric Hurwitz candidate yet; arithmetic envelope + paths give pure-even multi-k families without branch-cycle geometry.

---

## Goal

Find a **non-rigid geometric family** (positive-dimensional Hurwitz space) that is
**already pure-even** and **specialises onto several fixed-\(k\) arithmetic families**
(LSW \(k=-4\), flagship \(k=-8/5\), classical \(k=4/5\), …).

**Abandon** rigid \(\varphi\) (disc \(=5\cdot\square\) over \(\mathbb{Q}\)).

This is the most ambitious fusion route. Prior status: **no concrete Hurwitz candidate**.

---

## A. Arithmetic multi-\(k\) envelope (not Hurwitz)

Over \(\mathbb{Q}(m,s)\):

$$\alpha(m,s)=256 m^2-\frac{3125\, s^4}{256},\qquad \beta(m,s)=s\cdot\alpha(m,s)$$

$$\operatorname{disc}=(256\,\alpha^2 m)^2\quad\text{(identically square)}.$$

- Disc identity: **True**
- Catalogue seeds recovered: **26/26**
- Hurwitz-geometric? **False**
- This is the universal LSW-type parametrization of pure-even points on rays β=kα. It is pure-even and multi-k by construction, but it is an arithmetic parametrization of a subvariety of the BJ even surface — not a positive-dimensional Hurwitz space of branched covers with prescribed geometric monodromy data independent of the BJ embedding.

### Recovery by freezing \(s=k\)

- k=-4: 4/4 seeds on envelope
- k=4: 2/2 seeds on envelope
- k=-8/5: 4/4 seeds on envelope
- k=8/5: 2/2 seeds on envelope
- k=4/5: 3/3 seeds on envelope
- k=-4/5: 2/2 seeds on envelope
- k=-12/5: 3/3 seeds on envelope
- k=12/5: 2/2 seeds on envelope
- k=-16/5: 2/2 seeds on envelope
- k=16/5: 2/2 seeds on envelope

### Sample specialisations

- m=3 s=-4: α=-821 β=3284 → HIT_A5 S5TransitiveSubgroups.A5
- m=5 s=-8/5: α=6320 β=-10112 → HIT_A5 S5TransitiveSubgroups.A5
- m=5 s=4/5: α=6395 β=5116 → HIT_A5 S5TransitiveSubgroups.A5
- m=2 s=-12/5: α=619 β=-1485 → odd_monodromy S5TransitiveSubgroups.S5
- m=7 s=4: α=9419 β=37676 → HIT_A5 S5TransitiveSubgroups.A5

---

## B. Cross-\(k\) pure-even rational curves (ansatz search)

- Tested curves: **1680**
- Pure-even hits **not** contained in a single \(k\)-ray: **0**
- Priority seed pairs: 10

_No non-ray pure-even monomial/Bezier/quad bridge through two different-\(k\) catalogue seeds in the scanned ansätze. (Consistent with earlier even-surface scan.)_

---

## C. Hurwitz-adjacent / deformation ansätze

### C1. Mestre \(f-tr\)

- Tested: 396
- Multi-seed pure-even hits: **0**


### C2. Low-degree biparameter BJ

- Tested: 802
- Multi-\(k\) pure-even hits: **0**
- Envelope control disc□: True
- Low-degree biparam polynomial ansätze rarely have disc identically square except forms related to the k-ray envelope (or degenerate).


### C3. Survey of constructions

| Name | Pure-even? | Multi-\(k\)? | Hurwitz? | Status |
|------|:----------:|:------------:|:--------:|--------|
| LSW (Lavallee–Spearman–Williams) | True | False | arithmetic family with A5 specialisations; not a full Hurwitz moduli description | -4 |
| Fixed-k pure-even slices (enlarged catalogue) | True | False | False | one k each |
| Arithmetic multi-k envelope | True | True | False | arithmetic candidate only — not geometric Hurwitz |
| Homogenisation rays | True | False | False |  |
| Rigid φ (3A,3A,5A) | False | None | None | ABANDONED for fusion over Q (disc=5·□) |
| Mestre f-tr / Arala lines | rare in scans | False | None | no multi-seed pure-even hit in prior bounds |
| Positive-dim A5 Hurwitz space (abstract) | unknown | unknown | True | NO CONCRETE EQUATION candidate that is pure-even and hits several fixed-k families |

**Concrete Hurwitz pure-even multi-\(k\) candidate:** `None`

---

## D. Structure of the pure-even locus (why paths work arithmetically)

- Fixed-\(m\) varying-\(k\) curve disc□: **True**
- Catalogue seeds on m=5 curve: `[]`
- Catalogue seeds on m=55 curve: `[]`

Curves of fixed m and varying k are pure-even and cross all k-slices, but each meets C_k at only the single point with that m. Catalogue seeds sit at different m on different k, so fixed-m curves recover at most one seed per k and only when that seed's m matches. The 2-param envelope (m,k both free) recovers all; any 1-param curve through two catalogue seeds from different k with different m must be a non-constant path in (m,k)-space — i.e. a curve m(u), k(u) with disc still identically square (always true on the envelope).

Given seeds (αi,βi) on the envelope with parameters (mi,ki), i=1,2: any rational path (m(u),k(u)) with (m(0),k(0))=(m1,k1), (m(1),k(1))=(m2,k2) gives a pure-even 1-param family through both seeds via the envelope formulas. Example: linear path m=(1-u)m1+u m2, k=(1-u)k1+u k2.

### Explicit pure-even paths through cross-\(k\) seed pairs

#### flagship (\(k=-8/5\)) → lsw_m100 (\(k=-4\))
- m-path: `25*u/8 + 5/16`
- k-path: `-12*u/5 - 8/5`
- disc identically square: **True**
- endpoints match seeds: **True**
- midpoint sample: `None`
- Hurwitz-geometric? **False**
- Pure-even 1-param family through two different-k seeds via envelope path. Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family over P1 is A5-or-smaller for specialisations; no branch-cycle description.

#### flagship (\(k=-8/5\)) → classical (\(k=4/5\))
- m-path: `5/16`
- k-path: `12*u/5 - 8/5`
- disc identically square: **True**
- endpoints match seeds: **True**
- midpoint sample: `None`
- Hurwitz-geometric? **False**
- Pure-even 1-param family through two different-k seeds via envelope path. Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family over P1 is A5-or-smaller for specialisations; no branch-cycle description.

#### flagship (\(k=-8/5\)) → s180 (\(k=-12/5\))
- m-path: `5*u/8 + 5/16`
- k-path: `-4*u/5 - 8/5`
- disc identically square: **True**
- endpoints match seeds: **True**
- midpoint sample: `None`
- Hurwitz-geometric? **False**
- Pure-even 1-param family through two different-k seeds via envelope path. Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family over P1 is A5-or-smaller for specialisations; no branch-cycle description.

#### classical (\(k=4/5\)) → lsw_m100 (\(k=-4\))
- m-path: `25*u/8 + 5/16`
- k-path: `4/5 - 24*u/5`
- disc identically square: **True**
- endpoints match seeds: **True**
- midpoint sample: `{'u': '1/2', 'alpha': 820, 'beta': -1312, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}`
- Hurwitz-geometric? **False**
- Pure-even 1-param family through two different-k seeds via envelope path. Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family over P1 is A5-or-smaller for specialisations; no branch-cycle description.

#### flag_145 (\(k=-8/5\)) → lsw_124m (\(k=-4\))
- m-path: `21*u/8 + 15/16`
- k-path: `-12*u/5 - 8/5`
- disc identically square: **True**
- endpoints match seeds: **True**
- midpoint sample: `None`
- Hurwitz-geometric? **False**
- Pure-even 1-param family through two different-k seeds via envelope path. Arithmetic, not a Hurwitz space of covers. Geometric monodromy as a family over P1 is A5-or-smaller for specialisations; no branch-cycle description.

---

## E. Conclusions

1. **Abandoned \(\varphi\)** for this route: rigid, not pure-even over \(\mathbb{Q}\).

2. **Arithmetic multi-\(k\) envelope** over \(\mathbb{Q}(m,s)\) is pure-even and specialises onto **all** fixed-\(k\) slices (and essentially all catalogue seeds on them). It is **not** a Hurwitz space: no independent branch-cycle / geometric monodromy package.

3. **Constructive pure-even 1-param paths** through any two envelope seeds (including different \(k\)) exist by linear paths in \((m,k)\)-space. Example: flagship \(\leftrightarrow\) LSW, flagship \(\leftrightarrow\) classical. These solve the *arithmetic* multi-seed pure-even problem across \(k\).

4. **No concrete positive-dimensional Hurwitz candidate** was found that is already pure-even and maps onto several fixed-\(k\) families with geometric branch data. Cross-\(k\) polynomial curve ansätze on the even surface and Mestre/biparam searches produced **no** non-envelope geometric hit in bounds.

5. **Programme split:**
   - *Arithmetic fusion fuel:* envelope + paths (pure-even multi-\(k\) over \(\mathbb{Q}\)).
   - *Geometric fusion (ambitious):* still **open / no candidate** — would need a true family of covers with dim\(>0\) Hurwitz data whose BJ/Hilbert specialisations land on several \(C_k\).

### Recommended stance

- Treat the **\((m,s)\)-envelope** and **cross-\(k\) envelope paths** as the explicit pure-even multi-family arithmetic object.
- Do **not** claim Hurwitz geometry for them.
- Further geometric work must start from known positive-dim \(A_5\) Hurwitz strata (or other non-rigid constructions) and *test* pure-even + multi-\(k\) specialisation — not from \(\varphi\).

_Generated by nonrigid_hurwitz_search.py_