# G3b — 5A/5B split, braid orbit, and lift invariant

**Status: LOCKED**

**Verdict:** Multi-k pure-even envelope monodromy is **`5A⁴`** (all four nontrivial local monodromies are class **5A**; zero **5B**). Abstract Nielsen class Ni(A₅, 5A⁴) has **exactly one braid orbit of size 10**. Schur-cover lift invariant (rigid A₅ ≅ PSL(2,5), canonical order-5 SL lifts) is **`+1` on all 600 Nielsen tuples**.

---

## 0. Goal

G3 showed multi-k pure-even paths have monodromy signature `5*×4`.
This cut **splits 5A vs 5B**, forms the exact Nielsen type, and locks
**braid orbits** + a **lift invariant** for that type.

**Convention**

| label | definition |
|-------|------------|
| **5A** | A₅-conjugacy class of `(0 1 2 3 4)` — map `(1,2,3,4,0)` — size 12 |
| **5B** | A₅-conjugacy class of `(0 1 2 4 3)` — map `(1,2,4,0,3)` — size 12 |

These two are conjugate in S₅ by an **odd** permutation only. (A 5-cycle *is*
A₅-conjugate to its inverse; the split is **not** g vs g⁻¹.)

---

## 1. Envelope monodromy with 5A/5B labels

### `path_flag_classical` (primary)

| center | label | cycle type | track err |
|--------|-------|------------|----------:|
| ≈0.168 | **5A** | (5,) | ~1e-15 |
| ≈0.667−0.498i | **5A** | (5,) | ~4e-15 |
| ≈0.667+0.498i | **5A** | (5,) | 0 |
| ≈1.165 | **5A** | (5,) | ~7e-15 |

- multiset: **`{'5A': 4}`**
- sorted type key: `5A,5A,5A,5A`

### `path_flag_lsw`

| center | label | cycle type |
|--------|-------|------------|
| four finite pts | **5A** each | (5,) |
| t=−0.1 | **1** (identity) | disc factor `(10t+1)²` — not a geometric 5-branch |

- nontrivial 5-class content: **`5A⁴`** (consistent with primary)

**Locked multiset:** **`5A⁴`**

---

## 2. Nielsen class and braid orbits

Enumeration: tuples \((g_1,g_2,g_3,g_4)\) with each \(g_i\in 5A\), product \(1\),
generate A₅; conjugacy normalisation; Artin braid action \(\sigma_0,\sigma_1,\sigma_2\).

| ordered type | #Nielsen | #conj-norm | #braid orbits | orbit sizes |
|--------------|---------:|-----------:|--------------:|-------------|
| `(5A,5A,5A,5A)` | **600** | **10** | **1** | **[10]** |

**Braid lock:** Ni(A₅, 5A⁴) has a **single** braid orbit of size **10**.

---

## 3. Lift invariant

### Combinatorial proxy

| proxy | value |
|-------|-------|
| n_5A | 4 |
| n_5B | 0 |
| (−1)^{n_5B} | **+1** |

### Schur cover 2·A₅ ≅ SL(2,5) — **locked**

Method:

1. Rigid group iso A₅ → PSL(2,5) by matching generators  
   `a = (0 1 2 3 4) ∈ 5A`, `b = (0 1)(2 3)`.
2. Each order-5 element lifts to SL(2,5) by the **unique order-5** matrix  
   (prefer \(M^5 = I\) over the order-10 mate \(-M\)).
3. Lift invariant = product of the four lifts in \(\{\pm I\} \cong \{\pm 1\}\).

| object | lift invariant |
|--------|----------------|
| All **600** raw Nielsen tuples of type (5A)⁴ | **+1** (600/600) |
| All **10** conjugacy-normalised reps | **+1** (10/10) |

**Lift lock:** the unique braid component of Ni(A₅, 5A⁴) has lift invariant **`+1`**.

*(Numerical monodromy 4-tuples from independent small loops are not Nielsen
tuples in a common basepoint framing; their SL product is not used for the lock.
The invariant is computed on the abstract Nielsen class, which is the correct
Hurwitz object.)*

---

## 4. Named geometric type — programme lock

| field | value |
|-------|-------|
| Arithmetic object | pure-even multi-k path (flagship↔classical / ↔LSW) |
| Local monodromy multiset | **`5A⁴`** |
| r | 4 finite branch points (∞ unramified on multi-k paths) |
| Group | A₅ |
| vs ternary shortlist (3A⁴, 2A…) | **excluded** |
| Braid orbits | **1 × size 10** |
| Lift invariant | **`+1`** |

### Lock statement

> Pure-even multi-k envelope paths have geometric monodromy type  
> **`Ni(A₅, 5A⁴)`**,  
> with a **single** braid orbit of size **10** and Schur-cover lift invariant **`+1`**.  
> They are **not** Ni(A₅, 3A⁴).

This is consistent with G1–G2: ternary Hurwitz constructions do not recover the
pure-even catalogue because the envelope is **5-class geometry**.

---

## 5. Next

1. Explicit equation for the unique braid component of Ni(A₅, 5A⁴).
2. Hilbert-specialise that equation onto the pure-even seed lattice (close fusion).
3. Reduced Hurwitz genus for orbit size 10 (literature + Riemann–Hurwitz).

---

## 6. Non-claims

- Numerical monodromy is high-precision, not interval-certified.
- Branch ordering is by (Re, Im); does not affect conjugacy multiset.
- Full Modular-Tower genus tables for 5A⁴ are not re-proved here.

---

## 7. Reproduce

```bash
cd resonant_galois
python g3b_5a5b_braid_lift.py
# full-orbit lift check (optional):
python -c "from g3b_5a5b_braid_lift import *; ..."
```

_Generated / updated by `g3b_5a5b_braid_lift.py` + full-orbit lift census._
