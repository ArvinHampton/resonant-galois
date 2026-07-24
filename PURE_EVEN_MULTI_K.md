# Pure-even multi-k centre (theorem-grade core)

**Citable arithmetic centre of resonant-galois.**  
HQCC / ternary lattice appears only as motivation and specialisation source — not as a proved necessity claim.

---

## 1. Classical identities

**BJ discriminant.** For monic Bring–Jerrard quintics
```
f = x⁵ + α x + β ∈ F[x]  (char F ≠ 2,5)
disc(f) = 256 α⁵ + 3125 β⁴
```

**Operational A₅ criterion.**  
Irreducible + disc □ + Frobenius type (3,1,1) ⇒ Gal = A₅.

**Homogenisation lemma.**  
For seed s(x) = x⁵ + α x + β with disc(s) □,
```
f_t = x⁵ + α t⁴ x + β t⁵  ⇒  disc(f_t) = t²⁰ disc(s)
```
hence disc(f_t) is a square for all t ≠ 0.

---

## 2. Pure-even families

**Fixed-k slice.** For any k ≠ 0,
```
α(m) = 256 m² − (3125 k⁴)/256
β(m) = k · α(m)
```
Then disc is an identical square in m. Irreducible specialisations have even monodromy; sampling yields a positive density of Gal = A₅.

**Two-parameter envelope.**
```
α(m,s) = 256 m² − (3125 s⁴)/256
β(m,s) = s · α(m,s)
```
Pure-even over F(m,s). Freezing s = k recovers every fixed-k slice.

**Cross-k paths.** Any two envelope seeds (distinct k allowed) lie on an explicit pure-even 1-parameter family via a rational path in (m,k)-space, e.g.
```
m(u) = (1−u)m₁ + u m₂
k(u) = (1−u)k₁ + u k₂
```

---

## 3. Lattice specialisations

Coefficients drawn from the ternary / resonant lattice
```
L₀ ⊂ {3, 9, 27, 61, 80, 243, 539, …} and short combinations
```
produce many A₅ seeds. Grouping by k = β/α yields multi-seed pure-even slices (e.g. k = −4 LSW class, k = −8/5 flagship, k = 4/5 classical).

Stage D (density / height): disc-height asymptotics theorem; irreducibility density conjecture + computational evidence; Chebotarev cycle-type proxy consistent with A₅. See `STAGE_D_DENSITY.md`.

---

## 4. Explicit non-claims

| Claim | Status |
|-------|--------|
| Necessity theorem (HQCC axioms force Gal = A_n) | **Not claimed** — open research |
| Geometric Nielsen multi-k fusion | **Not claimed** — open |
| Physical S²-11DM²ET-X / G₄ / 539-step dynamics | **Not claimed** in this package |
| Structural template T forces disc □ | **False** (base M = S₅) |

---

## 5. One-sentence summary

Pure-even multi-k is a finished, theorem-grade arithmetic theory of even A₅ families over ℚ, organised by the ratio invariant k = β/α; the ternary lattice is a rich specialisation source, not a proved source of necessity.
