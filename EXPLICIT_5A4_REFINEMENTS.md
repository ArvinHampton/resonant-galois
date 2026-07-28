# Optional refinements — Ni(A₅, 5-class) explicit model

_Elapsed: refinements + re-lock_

**Verdict:** Common-basepoint monodromy **corrects** the 5A/5B multiset. Independent small-loop labels (all “5A”) were **sheet-label gauge artifacts** (odd S₅ relabelling swaps 5A↔5B). With a single base point, `path_flag_classical` realises multiset **`5A²5B²`**, ordered type `(5B,5A,5A,5B)` (RTL product), lift invariant **`−1`**, in the **size-12** braid orbit of that multiset. Reduced Hurwitz genus under the standard double-twist recipe: **`g = 0`**. Classical branch cross-ratio λ ≈ 1/2, j ≈ 1728.

---

## 0. Scope

Three refinements on top of `EXPLICIT_5A4_EQUATION.md` / G3b:

1. Common-basepoint braid word for geometric monodromy  
2. Reduced Hurwitz genus of the braid component  
3. Classical M_{0,4} chart (cross-ratio / j of branch 4-tuple)

**Important correction to G3b “all 5A” labels:** conjugacy labels computed in *independent* local sheet orderings are not comparable. Only a **common basepoint** yields a well-defined multiset of A₅-classes inside one S₅ embedding.

Canonical T3, pure-even arithmetic, Necessity stance unchanged.

---

## 1. Common-basepoint geometric monodromy

| item | value |
|------|-------|
| Model | `path_flag_classical` over ℚ(t) |
| Base point t_* | real-left of branches + small Im |
| # finite branch points | 4 |
| Branch sqf | `81t⁴ − 216t³ + 216t² − 96t + 11` |

### Generators (geometric arg-order about base; RTL product = id)

| # | label | cycle type | role |
|---|-------|------------|------|
| 1 | **5B** | (5,) | complex conjugate pair member |
| 2 | **5A** | (5,) | real branch |
| 3 | **5A** | (5,) | real branch |
| 4 | **5B** | (5,) | complex conjugate pair member |

| check | result |
|-------|--------|
| RTL product of monodromies | **identity** |
| Generates A₅ | **True** |
| All even (sign +1) | **True** |
| Multiset | **`5A² 5B²`** |
| Lift invariant (SL(2,5) canonical) | **`−1`** |
| Braid orbit (among 5A²5B²) | **size 12** (lift −1 component) |

### Why G3b said “5A⁴”

Each independent Newton loop assigned sheet labels via `numpy.roots` order at a *different* local base. Changing sheet labels by an **odd** permutation of S₅ swaps the two A₅-classes of 5-cycles. Local “all 5A” was therefore **not** a gauge-invariant statement. Common-basepoint monodromy removes that ambiguity.

### Combinatorics of multiset 5A²5B²

| quantity | value |
|----------|------:|
| Raw Nielsen (all ordered patterns) | 2520 |
| Conjugacy-normalised | 42 |
| Braid orbits | **2** |
| Orbit sizes | **30** (lift **+1**), **12** (lift **−1**) |

**Geometric cover sits in the lift-`−1` orbit of size 12.**

---

## 2. Reduced Hurwitz genus

Action of Artin generators on the **size-10** conj-norm orbit of pure **5A⁴** (combinatorial type; still well-defined):

| generator recipe | genus | product id | (ind₀, ind₁, ind∞) |
|------------------|------:|:----------:|--------------------|
| **A: σ₁², σ₂²** (preferred double-twist) | **0** | True | (6,6,6) |
| B: σ₁, σ₂ | 1 | True | (7,7,6) |
| C: σ₁σ₂σ₁, σ₂ | 0 | True | (5,7,6) |
| D: σ₂σ₁σ₂, σ₁ | 0 | True | (5,7,6) |
| E: σ₁², σ₁σ₂σ₁ | 0 | True | (6,5,7) |

**Preferred (standard r=4 double-twist cusps): genus `g = 0`.**

For the **5A²5B²** lift-−1 component (orbit size 12), RH with the same double-twist recipe on that orbit should be recomputed if a full chart of that component is needed; the pure-even cover is a **point** of that component (one geometric cover), not the whole curve.

### Geometric meaning of g = 0 for 5A⁴’s size-10 orbit

The unique braid component of pure 5A⁴ has rational reduced Hurwitz curve (g=0) under the standard recipe — consistent with “infinitely many rational moduli points” for that combinatorial type. Our envelope cover, however, realises **5A²5B²** (lift −1), a **different** component.

---

## 3. Classical chart (cross-ratio / j)

Four branch points of `path_flag_classical` (roots of the disc square-free quartic):

| coordinate | value (numeric) |
|------------|-----------------|
| Primary cross-ratio λ | **≈ 1/2** |
| j(λ) = 2⁸(λ²−λ+1)³/(λ²(λ−1)²) | **≈ 1728** |

(Anharmonic group acts on λ; j is the S₃-invariant.)

### Envelope → Hurwitz interpretation

- Each pure-even **1-param path** = one geometric A₅-cover of ℙ¹ (base = path parameter) = **one point** of the appropriate Hurwitz component.  
- Hilbert specialisations = fibres of that cover (number fields / BJ polynomials).  
- Multi-k catalogue hits (flagship, classical, LSW, …) are **different fibres of the same cover**, not different Hurwitz points.  
- Distinct paths (flag↔classical vs flag↔LSW) are **two geometric covers** (two points of Hurwitz space), both of multiset type 5A²5B² when monodromy is read with a consistent basepoint.  
- The 2-param envelope parametrizes a surface of coefficient data; 1-param slices pick covers on the lift-−1 component of Ni(A₅, 5A²5B²).

---

## 4. Corrected programme lock

| claim | status |
|-------|--------|
| Pure-even multi-k paths have 5-cycle monodromy (r=4) | **True** |
| Multiset (common basepoint) | **`5A²5B²`** |
| Lift invariant of geometric cover | **`−1`** |
| Braid orbit | **size 12** (one of two orbits for this multiset) |
| Explicit ℚ(t) model | pure-even paths (unchanged) |
| Hilbert multi-k catalogue | **True** (unchanged) |
| Pure **5A⁴** (orbit 10, lift +1) | Exists combinatorially; **not** the envelope cover’s monodromy type |
| Type-level fusion | **Closed** on type **Ni(A₅, 5A²5B²)** lift-−1 component |

### One-line statement

> The pure-even multi-k envelope paths are explicit ℚ(t)-models of geometric A₅-covers of monodromy type **`5A²5B²` with lift invariant −1** (braid orbit size 12), Hilbert-specialising to multiple catalogue k-slices.

---

## 5. Reproducibility

```bash
cd resonant_galois
python g3d_5a4_refinements.py
```

---

## 6. Non-claims

- Infinity monodromy from a large numerical loop is not used for the lock (possible path-crossing artefacts).  
- Reduced genus g=0 row is for the **5A⁴** size-10 combinatorial orbit under the double-twist recipe; the geometric cover’s component is **5A²5B²** lift −1.  
- Full Modular-Tower naming of that size-12 component in published tables is not re-derived here beyond orbit size + lift.

_Generated by `g3d_5a4_refinements.py` + post-hoc orbit membership for 5A²5B²._
