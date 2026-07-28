# Can Nielsen labels take values in conjugacy classes of \(\mathrm{Gal}(f/R_n)\) with split Frobenius as dictionary?

_Elapsed: 1.88s_

**Answer:** YES as a type-dictionary after monodromy identification; NO as a literal equality Nielsen class = Frob class without geometry.

---

## 0. Three objects

### Nielsen label

An r-tuple C=(C_1,...,C_r) of conjugacy classes in a finite group G ≤ S_d (or abstract G with a permutation representation), such that there exist g_i ∈ C_i with g_1...g_r=1 generating G (Nielsen / Hurwitz).

- Nature: **geometric — monodromy of a cover of P¹\{branch points}**
- Takes values in: conjugacy classes of G_geom (geometric monodromy group)

### \(\mathrm{Gal}(f/R_n)\)

For monic separable f ∈ R_n[x], Gal(f/R_n) is Gal(K/R_n) where K is the splitting field of f over R_n, as a permutation group on the roots.

- Nature: **arithmetic — Aut of a Galois extension of number fields**
- Conjugacy classes: conjugacy classes of G_arith = Gal(f/R_n) ≤ S_deg

### Split Frobenius

A rational prime p that splits completely in R_n/Q, unramified in K/Q (or K/R_n). Then O_{R_n}/P ≅ F_p for P|p, and Frob_P ∈ Gal(K/R_n) is well-defined up to conjugacy; its cycle type on roots equals the factorisation type of f mod P ≅ f mod p (after reducing a Z-model of f).

- Role: Chebotarev dictionary over R_n, computable via F_p factorisation when p splits in R_n and f has a model reducing well.

---

## 1. The dictionary (what “yes” means)

```
  Nielsen class C_i  ⊂  G_geom
         │ specialise cover over R_n
         ▼
  conjugacy class in G_arith = Gal(f/R_n)  ≥  G_geom
         │ Chebotarev
         ▼
  Frob_P  (P | p,  p split in R_n/Q)
         │ reduce
         ▼
  factorisation type of f mod p   ←→  cycle type of class
```

For a cover φ: X→P¹ defined over R_n with geometric monodromy G_geom, and a fibre f = fibre polynomial at a non-branch R_n-point, G_geom ↪ G_arith = Gal(f/R_n) ⊆ S_d (after identifying sheets with roots). Often equality after Hilbert specialisation (arithmetic monodromy = geometric).

Nielsen classes label branch monodromy generators (loops around branch points). Frob classes label arithmetic primes. They are NOT the same loops: the dictionary identifies both as conjugacy classes in a common group G after specialisation, with cycle type as the coarse invariant.

Coarse dictionary: conjugacy class ↦ cycle type in S_d. Nielsen C_i has cycle type τ_i; Frob_p has cycle type = factorisation type of f mod p. Same G ⇒ same possible types; Chebotarev predicts densities of Frob types from class sizes in G_arith.

### Affirmative list

- Nielsen labels are conjugacy classes in G_geom.
- After specialising a cover over R_n to f ∈ R_n[x], G_geom ≤ Gal(f/R_n).
- Both Nielsen classes and Frob classes are then conjugacy classes in (a group identified with a subgroup of) Gal(f/R_n).
- Split primes p in R_n give a Chebotarev dictionary: Frob cycle types on f are readable by factoring over F_p, labelling classes in Gal(f/R_n).
- Thus Nielsen and split-Frob share a common value space: conjugacy classes / cycle types in that monodromy group — a legitimate dictionary of types.

### Negative list (do not over-claim)

- Nielsen labels are not defined as Frob classes; they come from π_1 of the punctured base, not from Spec O_K.
- Without a cover and specialisation, Gal(f/R_n) alone has no Nielsen data.
- Split Frob does not invent branch cycle types; it samples arithmetic classes.
- For f ∈ Q[x], Gal(f/R_n) may be smaller than Gal(f/Q) if R_n meets the splitting field — dictionary must track base change carefully.
- Cannot assign Nielsen class C_3^4 to an arbitrary BJ seed just from split-Frob histograms without geometric construction.

---

## 2. Coarse vs fine labels

| level | invariant | Nielsen | split Frob | match? |
|-------|-----------|---------|------------|--------|
| **Coarse** | cycle type in \(S_d\) | class \(\mapsto\) type | factorisation type mod \(p\) | **Yes** |
| **Fine** | class in \(G\) (e.g. 5A vs 5B) | full class | Artin symbol | only with extra structure |

**5A/5B:** A5 has two classes of 5-cycles, fused in S5. Cycle type (5) is shared. Nielsen labels for (3A,3A,5A) vs (3A,3A,5B) require the fine distinction. Split-Frob factorisation type cannot choose 5A vs 5B by itself.

- Coarse: YES — standard, computable, used throughout the programme
- Fine: POSSIBLE in principle via class field / resolvents, but split Frob cycle type alone does not separate same-cycle-type classes. Extra structure (e.g. action on resolvent rings, or Frobenius in extensions containing R_n) needed.

---

## 3. Base change caution

For f ∈ Q[x] with Gal(f/Q)=A5, if the splitting field K is linearly disjoint from R_n over Q, then Gal(f/R_n) ≅ Gal(f/Q) ≅ A5. If K ∩ R_n ≠ Q, Gal(f/R_n) is a proper quotient / subgroup situation (actually [K R_n : R_n] = [K:K∩R_n], so Gal(f/R_n) ≅ Gal(K/K∩R_n) may be smaller than A5).

For typical A5 fields (disc not divisible by special primes of R_n only), K ∩ R_n = Q for small n, so Gal(f/R_n)=A5 still. Always check compositum.

If f ∈ R_n[x] \ Q[x], Gal(f/R_n) is the native group; split Frob still labels its conjugacy classes via reductions. Nielsen labels require a cover over R_n whose fibre is f.

---

## 4. Demo — flagship \(A_5\) seed + split Frob in \(R_n\)

Fibre \(f=x^5-55x+88\), \(\mathrm{Gal}(f/\mathbb{Q})=A_5\). Nielsen-relevant \(A_5\) cycle types vs observed split-Frob types.

### Abstract \(A_5\) dictionary (cycle type ↔ Nielsen role)

| cycle type | \(A_5\) class | Nielsen role |
|------------|---------------|--------------|
| \((1^5)\) | id | not a branch class |
| \((2,2,1)\) | 2A | e.g. in \((2A,3A,5A)\) |
| \((3,1,1)\) | 3A | ternary / \(C_3\) / \(3A^4\) |
| \((5)\) | 5A or 5B | period class; fine label open |

### \(n=5\) — split primes: 20 (sample [11, 19, 29, 31, 41, 59, 61, 71])

| Frob pattern (split \(p\)) | count | \(A_5\) / Nielsen |
|---------------------------|------:|------------------|
| `(5,)` | 9 | 5A/5B (5-cycles; two classes, same cycle type) — period / 5-cycle branch class; 5A vs 5B not split by cycle type |
| `(1, 1, 3)` | 6 | 3A (3-cycle) — ternary branch class C_3 / 3A — central to 3A^4, (3A,3A,5A) |
| `(1, 2, 2)` | 4 | 2A (double transposition) — possible branch class (e.g. in (2A,3A,5A)) |

Control (non-split) patterns: `{'(5,)': 10, '(1, 1, 3)': 9, '(1, 2, 2)': 4, '(1, 1, 1, 1, 1)': 1}`

### \(n=7\) — split primes: 13 (sample [13, 29, 41, 43, 71, 83, 97, 113])

| Frob pattern (split \(p\)) | count | \(A_5\) / Nielsen |
|---------------------------|------:|------------------|
| `(1, 1, 3)` | 5 | 3A (3-cycle) — ternary branch class C_3 / 3A — central to 3A^4, (3A,3A,5A) |
| `(5,)` | 5 | 5A/5B (5-cycles; two classes, same cycle type) — period / 5-cycle branch class; 5A vs 5B not split by cycle type |
| `(1, 2, 2)` | 2 | 2A (double transposition) — possible branch class (e.g. in (2A,3A,5A)) |
| `(1, 1, 1, 1, 1)` | 1 | id — not a branch class (trivial) |

Control (non-split) patterns: `{'(5,)': 14, '(1, 1, 3)': 10, '(1, 2, 2)': 6}`

### \(n=11\) — split primes: 8 (sample [23, 43, 67, 89, 109, 131, 197, 199])

| Frob pattern (split \(p\)) | count | \(A_5\) / Nielsen |
|---------------------------|------:|------------------|
| `(1, 2, 2)` | 3 | 2A (double transposition) — possible branch class (e.g. in (2A,3A,5A)) |
| `(1, 1, 3)` | 3 | 3A (3-cycle) — ternary branch class C_3 / 3A — central to 3A^4, (3A,3A,5A) |
| `(5,)` | 2 | 5A/5B (5-cycles; two classes, same cycle type) — period / 5-cycle branch class; 5A vs 5B not split by cycle type |

Control (non-split) patterns: `{'(5,)': 15, '(1, 1, 3)': 10, '(1, 2, 2)': 4, '(1, 1, 1, 1, 1)': 1}`

### \(n=15\) — split primes: 9 (sample [29, 31, 59, 61, 89, 149, 151, 179])

| Frob pattern (split \(p\)) | count | \(A_5\) / Nielsen |
|---------------------------|------:|------------------|
| `(5,)` | 5 | 5A/5B (5-cycles; two classes, same cycle type) — period / 5-cycle branch class; 5A vs 5B not split by cycle type |
| `(1, 1, 3)` | 3 | 3A (3-cycle) — ternary branch class C_3 / 3A — central to 3A^4, (3A,3A,5A) |
| `(1, 2, 2)` | 1 | 2A (double transposition) — possible branch class (e.g. in (2A,3A,5A)) |

Control (non-split) patterns: `{'(5,)': 13, '(1, 1, 3)': 9, '(1, 2, 2)': 7, '(1, 1, 1, 1, 1)': 1}`

**Types seen at split primes (union over n):** `['(1, 1, 1, 1, 1)', '(1, 1, 3)', '(1, 2, 2)', '(5,)']`

Cycle types of split Frob on the flagship match A5 class cycle types (3,1,1), (5), (2,2,1), id — the same types that index Nielsen classes for A5 covers. Fine Nielsen labels (5A vs 5B) are NOT determined by cycle type alone; need roots of unity / class field refinements.

---

## 5. Programme consequences

| use case | allowed? |
|----------|:--------:|
| Label Frob types of pure-even \(A_5\) fibres by Nielsen-relevant cycle types (3A, 5A/B, 2A) | **Yes** (coarse) |
| Use only \(p\) split in \(R_n\) as the Chebotarev sample for “resonant” arithmetic | **Yes** |
| Assert a seed is a \(3A^4\) fibre from split-Frob histogram alone | **No** |
| Identify \(\mathrm{Gal}(f/R_n)\) classes with Nielsen classes after a cover specialisation over \(R_n\) | **Yes** (standard) |
| Separate 5A from 5B with split-Frob cycle type only | **No** |
| Replace geometric Nielsen theory by split Frob | **No** |

---

## 6. Locked answer

> Can Nielsen labels take values in conjugacy classes of \(\mathrm{Gal}(f/R_n)\)
> with split Frobenius as dictionary?

**Yes, as a coarse type-dictionary after monodromy identification; no, as a literal substitution.**

1. Nielsen labels \(\in\) conjugacy classes of \(G_{\mathrm{geom}}\).
2. Specialisation over \(R_n\) embeds \(G_{\mathrm{geom}}\le\mathrm{Gal}(f/R_n)\).
3. Split Frob supplies Chebotarev sampling of conjugacy classes of
   \(\mathrm{Gal}(f/R_n)\) via factorisation over \(\mathbb{F}_p\).
4. Shared invariant: **cycle type** (and, with more work, fine class in \(G\)).
5. Split Frob does not create Nielsen data without a cover; it reads arithmetic
   classes in the same group-theoretic currency geometric monodromy uses.

```bash
python nielsen_gal_rn_dictionary.py
```

_Generated by nielsen_gal_rn_dictionary.py_