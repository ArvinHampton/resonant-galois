# G3 — monodromy identification of the pure-even envelope

_Elapsed: 13.38s_

**Verdict:** G3 envelope monodromy (13.38s). Envelope disc□ identity=True. Multi-k paths `path_flag_classical` / `path_flag_lsw` have monodromy signature **5A/5B×4** (four 5-cycles; ∞ trivial) — **not** the ternary shortlist 3A⁴/2A*. Geometric multi-k named as 5-class r=4 multiset: **True** (5A vs 5B open). Ternary Nielsen ID: False.

---

## 0. Goal

The pure-even envelope

```text
α = 256 m² − 3125 k⁴ / 256,   β = k α
disc(x⁵+αx+β) = (256 α² m)²   (identical square)
```

already supplies **arithmetic multi-k** (cross-k paths hit flagship, classical, LSW, …).
G3 asks for a **Nielsen name**: the geometric monodromy of natural 1-parameter
slices as a permutation representation of π₁(ℙ¹ ∖ branch locus) → A₅,
compared to the r=4 shortlist (3A⁴, 2A3A³, 2A²3A², …).

---

## 1. Envelope structure (2-parameter)

| item | value |
|------|-------|
| dimension | 2 |
| disc identical square | **True** |
| branch divisor | α(m,k)=0 (and m=0 degenerate) |
| branch equation | `(-3125*k**4 + 65536*m**2)/256 = 0` |
| foliation | fixed-k rays are pure-even 1-param; cross-k paths join multi-seed ratios |

The envelope is a rational surface of pure-even BJ polynomials. Geometric monodromy is well-defined on 1-param slices (paths/rays). A full Nielsen ID is an ID of those slices' monodromy types.

---

## 2. One-parameter slices — monodromy

### `path_flag_classical`

_same-m linear-k; arithmetic multi-k flagship↔classical_

- disc degree: **16**; sqf factors: `[('81*t**4 - 216*t**3 + 216*t**2 - 96*t + 11', 4)]`; all mult even: **True**
- finite monodromy classes: **5A/5B×4**
- with ∞: **5A/5B×4**
- group ⟨finite⟩: order=5 A5=False all_even=True
- group ⟨finite+∞⟩: order=5 A5=False
- Nielsen match (with ∞): **none**
- catalogue hits: [{'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}, {'tag': 'classical_m', 'k': '-4/5', 't': '1/3', 'alpha': 20, 'beta': -16}, {'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}]; k=['-4/5', '-8/5', '4/5']; multi-k=**True**
- sample fibre Gal: {'t': '2', 'alpha': -1255, 'beta': -4016, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}

Local monodromy (finite):

| center | cycle type | class | track err |
|--------|------------|-------|----------:|
| (0.16821707292625981+0j) | (5,) | 5A/5B | 2.886579864025407e-15 |
| (0.6666666666666666-0.49844959374040687j) | (5,) | 5A/5B | 5.7902148899581444e-15 |
| (0.6666666666666666+0.49844959374040687j) | (5,) | 5A/5B | 0.0 |
| (1.1651162604070735+0j) | (5,) | 5A/5B | 4.441239029644353e-15 |

Infinity: type=(1, 1, 1, 1, 1) class=**1** err=1.7563101230182143e-07

### `path_flag_lsw`

_linear (m,k) flagship→LSW_

- disc degree: **18**; sqf factors: `[('10*t + 1', 2), ('81*t**4 + 216*t**3 - 284*t**2 - 4*t + 11', 4)]`; all mult even: **True**
- finite monodromy classes: **5A/5B×4**
- with ∞: **5A/5B×4**
- group ⟨finite⟩: order=60 A5=True all_even=True
- group ⟨finite+∞⟩: order=60 A5=True
- Nielsen match (with ∞): **none**
- catalogue hits: [{'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'lsw_m100', 'k': '-4', 't': '1', 'alpha': -100, 'beta': 400}, {'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'lsw_m100', 'k': '-4', 't': '1', 'alpha': -100, 'beta': 400}]; k=['-4', '-8/5']; multi-k=**True**
- sample fibre Gal: {'t': '2', 'alpha': -9455, 'beta': 60512, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}

Local monodromy (finite):

| center | cycle type | class | track err |
|--------|------------|-------|----------:|
| (-0.1+0j) | (1, 1, 1, 1, 1) | 1 | 5.617333549722722e-15 |
| (-0.19104873969093802+0j) | (5,) | 5A/5B | 3.953391794989626e-15 |
| (0.20774114394235957+0j) | (5,) | 5A/5B | 6.358389842764979e-15 |
| (0.9434454977240734+0j) | (5,) | 5A/5B | 6.52013858905055e-15 |
| (-3.6268045686421617+0j) | (5,) | 5A/5B | 1.9560100690747752e-14 |

Infinity: type=(1, 1, 1, 1, 1) class=**1** err=2.6236815868147128e-11

### `ray_lsw_k_m4`

_fixed-k pure-even ray (single k)_

- disc degree: **10**; sqf factors: `[('t', 2), ('256*t**2 - 3125', 4)]`; all mult even: **True**
- finite monodromy classes: **5A/5B×2**
- with ∞: **2A+5A/5B×2**
- group ⟨finite⟩: order=5 A5=False all_even=True
- group ⟨finite+∞⟩: order=60 A5=True
- Nielsen match (with ∞): **none**
- catalogue hits: []; k=[]; multi-k=**False**
- sample fibre Gal: {'t': '1/2', 'alpha': -3061, 'beta': 12244, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}

Local monodromy (finite):

| center | cycle type | class | track err |
|--------|------------|-------|----------:|
| 0j | (1, 1, 1, 1, 1) | 1 | 7.628290615314396e-15 |
| (-3.4938562148434213+0j) | (5,) | 5A/5B | 1.3877787807814457e-14 |
| (3.4938562148434213+0j) | (5,) | 5A/5B | 7.944109290391274e-15 |

Infinity: type=(2, 2, 1) class=**2A** err=1.5583418786108402e-13

### `homog_flagship`

_homogenised flagship ray (Theorem 3); single seed k=-8/5_

- disc degree: **20**; sqf factors: `[('t', 20)]`; all mult even: **True**
- finite monodromy classes: ****
- with ∞: ****
- group ⟨finite⟩: order=1 A5=False all_even=True
- group ⟨finite+∞⟩: order=1 A5=False
- Nielsen match (with ∞): **none**
- catalogue hits: [{'tag': 'flagship', 'k': '-8/5', 't': '1', 'alpha': -55, 'beta': 88}, {'tag': 'flagship', 'k': '-8/5', 't': '1', 'alpha': -55, 'beta': 88}]; k=['-8/5']; multi-k=**False**
- sample fibre Gal: {'t': '2', 'alpha': -880, 'beta': 2816, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}

Local monodromy (finite):

| center | cycle type | class | track err |
|--------|------------|-------|----------:|
| 0j | (1, 1, 1, 1, 1) | 1 | 1.153297605929378e-16 |

Infinity: type=(1, 1, 1, 1, 1) class=**1** err=3.0926163410039053e-12

### `homog_classical`

_homogenised classical; k=4/5_

- disc degree: **20**; sqf factors: `[('t', 20)]`; all mult even: **True**
- finite monodromy classes: ****
- with ∞: ****
- group ⟨finite⟩: order=1 A5=False all_even=True
- group ⟨finite+∞⟩: order=1 A5=False
- Nielsen match (with ∞): **none**
- catalogue hits: [{'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}, {'tag': 'classical_m', 'k': '-4/5', 't': '-1', 'alpha': 20, 'beta': -16}, {'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}]; k=['-4/5', '4/5']; multi-k=**True**
- sample fibre Gal: {'t': '2', 'alpha': 320, 'beta': 512, 'status': 'HIT_A5', 'gal': 'S5TransitiveSubgroups.A5'}

Local monodromy (finite):

| center | cycle type | class | track err |
|--------|------------|-------|----------:|
| 0j | (1, 1, 1, 1, 1) | 1 | 1.8670653977139892e-16 |

Infinity: type=(1, 1, 1, 1, 1) class=**1** err=2.7754494148515964e-12

---

## 3. Nielsen shortlist comparison

Target classes (from `A5_HURWITZ_R4.md`):

| Nielsen type | class multiset | orbit | g |
|--------------|----------------|------:|--:|
| 3A⁴ | 3A×4 | 18 | 0 |
| 2A 3A³ | 2A+3A×3 | 96 | 0 |
| 2A² 3A² | 2A×2+3A×2 | 108 | 0 |
| 2A 3A² 5* | 2A+3A×2+5* | 240 | 0 |
| 3A³ 5* | 3A×3+5* | 40–60 | 0 |

A **name** is reported when the multiset of local monodromy conjugacy classes
(finite branch points, optionally including ∞) equals a shortlist multiset.

| family | multi-k? | signature (w/ ∞) | Nielsen ID |
|--------|:--------:|------------------|------------|
| path_flag_classical | True | 5A/5B×4 | — |
| path_flag_lsw | True | 5A/5B×4 | — |
| ray_lsw_k_m4 | False | 2A+5A/5B×2 | — |
| homog_flagship | False |  | — |
| homog_classical | True |  | — |

---

## 4. Multi-k / geometric conclusion

| test | result |
|------|--------|
| Pure-even envelope disc identity | **True** |
| Multi-k arithmetic paths present | **True** (`path_flag_classical`, `path_flag_lsw`, …) |
| Local monodromy computed | **True** |
| Match to **ternary** shortlist (3A⁴, 2A3A³, 2A²3A²) | **False** |
| **Geometric multiset ID of multi-k paths** | **`5A/5B × 4`** (four 5-cycles; ∞ trivial) |

### Main positive result

The multi-k pure-even paths have numerical monodromy signature

```text
5* × 4     (four finite branch points, each a 5-cycle; monodromy at ∞ = id)
```

That is **not** the pure-ternary / mixed shortlist used in G1–G2. It is the class
multiset of an **r=4 type built from A₅'s 5-classes** (5A⁴ / 5A³5B / 5A²5B² up to
outer automorphism swapping 5A↔5B — not resolved in this cut).

| family | multi-k | signature | reading |
|--------|:-------:|-----------|---------|
| `path_flag_classical` | True | 5*×4 | **5-class r=4 family** |
| `path_flag_lsw` | True | 5*×4 | **5-class r=4 family** |
| `ray_lsw_k_m4` | False | 2A+5*×2 | single-k ray; ∞ = 2A |
| homog rays | single k | degenerate at 0 | not cross-k monodromy |

**Consequence for G1/G2:** negative catalogue hits on 3A⁴ / 2A* covers are
**consistent** with this ID — the pure-even multi-k envelope is not trying to be
those Hurwitz components.

**Named geometric multi-k (ternary shortlist):** False.  
**Named geometric multi-k (5-class multiset `5*^4`):** **True** for the multi-k
envelope paths, pending 5A vs 5B refinement and braid-orbit / lift-invariant lock.

### What this cut established

1. Disc/branch analysis of pure-even 1-param slices (square multiplicities).
2. Numerical local monodromy: multi-k paths → **four 5-cycles**.
3. Ternary shortlist comparison: **no match** (by design of the geometry).
4. Arithmetic multi-k catalogue confirmation on the same paths.
5. Clear redirect: further geometric fusion work on the envelope should target
   **5-class Nielsen types**, not 3A⁴.

### Next

1. Distinguish **5A vs 5B** (lift invariant / complex conjugation / fixed 5-cycle class).
2. Lock braid orbit + reduced Hurwitz genus for the matching 5*^4 type.
3. Produce an explicit Hurwitz equation for that type and re-test lattice specialisations.
4. Literature: pure-even BJ pencil vs known icosahedral / 5-cycle A₅ families.

---

## 5. Non-claims

- Numerical monodromy is not a certified braid factorization.
- Matching class multisets is necessary but not sufficient for a full Hurwitz
  component ID (braid orbit / lift invariant may refine further).
- Does not reopen pure-even arithmetic, Canonical T3, or Necessity.

_Generated by `g3_envelope_monodromy.py`._
