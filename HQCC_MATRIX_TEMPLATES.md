# HQCC matrix templates (locked exploration)

**Status:** templates verified; necessity theorem still open; pure-even multi-k remains the finished arithmetic centre.

---

## 1. Base structural matrix M

```
M = [[0, 1, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [3, 0, 0, 80, 0],
     [0, 0, 0, 0, 1],
     [61, 0, 0, -3, 0]]
```

| Item | Value |
|------|--------|
| Characteristic polynomial | χ_M = x⁵ + 3x³ − 3x² − 4889 |
| Discriminant | 1781436218361244736 — **not a square** |
| Galois group | **S₅** (3-cycles present, odd monodromy) |
| Flux fingerprint | 4889 = 4880 + 3²; 80 = 4880/61 |

**Block reading**

| Block | Role |
|-------|------|
| UL 3×3 | Companion-like ternary block (entry 3) |
| Couplings | 80 and 61 (flux / puncture integers) |
| LR 2×2 | [[0,1], [−3,0]] |

This is the concrete **evenness obstruction**: full ternary + flux structure, disc odd ⇒ S₅.

---

## 2. Structural template T(a,b,c,d,e,f)

```
T(a,b,c,d,e,f) = [[0, 1, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [a, 0, 0, b, e],
                  [0, 0, 0, 0, 1],
                  [c, f, 0, d, 0]]
```

Base M is the specialisation (a,b,c,d,e,f) = (3, 80, 61, −3, 0, 0).

**Characteristic polynomial (exact)**

```
χ_T = x⁵ − d x³ − (a + ef) x² − (bf + ce) x + (ad − bc)
```

Parameters from the resonant / model lattice {3, 9, 27, 61, 80, 243, 539, …} and short combinations.

**Deformation result**  
Lattice scan + square-disc gate produced explicit A₅ realisations. Example:

```
T(3, 0, 0, −3, 1, 3)  →  χ = x⁵ + 3x³ − 6x² − 9
disc = 3470769 = 1863²,  Gal = A₅
```

Same template shape; parameters chosen so disc is square ⇒ A₅. Structure alone does not force evenness.

---

## 3. BJ-embed subclass

Impose

```
d = 0,  a = −ef
```

Then

```
χ = x⁵ − (bf + ce) x − bc
```

i.e. Bring–Jerrard with α = −(bf+ce), β = −bc.

On this thin subclass the pure-even theory applies: fix k = β/α and run the classical envelope. Even monodromy becomes an identity.

| Property | Status |
|----------|--------|
| Disc identically square on pure-even rays | Yes (classical BJ) |
| Forced by full unrestricted T | **No** |
| Native HQCC labelling of (b,c,e,f) | Ansatz on top of the template |

---

## 4. Degree-6 enlargement T₆

Base:

```
χ = x⁶ + 3x⁴ − 3x² − 4889
Gal ≈ S₄ × C₂  (disc not square)
```

Sparse ternary specialisations produced 4 verified A₆ polynomials, e.g.

```
x⁶ − 3x⁴ + 9x² ± 18x + 9   (disc = 75816²)
```

Same lesson: ternary structure + disc gate works; structure alone does not force even monodromy.

---

## 5. Summary

| Claim | Verdict |
|-------|--------|
| Templates encode order-3 / flux data | Yes |
| Templates produce A₅ / A₆ after disc gate | Yes |
| Templates force disc □ | **No** (base M, base T₆) |
| BJ-embed recovers pure-even arithmetic | Yes, on a thin subclass |
| Templates alone yield a necessity theorem | **No** |

**Criterion 2 implication:** any axiom list claiming “resonant matrix shape ⇒ alternating monodromy” must be strictly stronger than membership in T(a,…,f). The known strengthening (BJ-embed + pure-even condition) is exactly the finished pure-even multi-k theory — not a new HQCC-native necessity theorem.

See also: `NECESSITY_THEOREM.md`, `EVENNESS_OBSTRUCTION.md`, `ARITHMETIC_MULTI_K.md`.
