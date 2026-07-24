# Rigid fibre t=3 — negative control (locked)

## Polynomial

```
φ(y) − 3 = 6y⁵ − 15y⁴ + 10y³ − 3
```

Monic ℚ model and primitive ℤ model:

```
monic_Q:  y⁵ − (5/2)y⁴ + (5/3)y³ − 1/2
monic_Z:  z⁵ − 15z⁴ + 60z³ − 3888   (z = 6y)
```

## Checks

| Check | Result |
|-------|--------|
| Disc form | disc = 5 · (square) — **not** a square in ℚ |
| Irreducibility | True (monic_Q and monic_Z) |
| Galois group | **S₅** (odd monodromy; Frobenius types include (1,4), (2,3), (1,1,1,2)) |

## Family identity

For the preferred rigid Belyi cover φ(y) = 6y⁵ − 15y⁴ + 10y³ over ℚ with passport (3,1,1)(3,1,1)(5):

```
disc(monic(φ − t)) = 5 · (square in t)   for rational t ∉ {0,1}
```

Consequence: every irreducible specialisation over ℚ has **odd** arithmetic monodromy (S₅).  
At t ∈ {0,1}, disc = 0 (critical / reducible).

## Contrast (locked)

| Side | Example | Disc | Parity | Gal |
|------|---------|------|--------|-----|
| Pure-even resonant | x⁵ − 55x + 88 (flagship) | □ | even | A₅ |
| Rigid cover fibre | monic(φ − 3) | 5·□ not □ | odd | S₅ |

## Programme reading

Further surgery on the same rigid φ/ℚ cannot produce even irreducible fibres.  
Geometric multi-k fusion must use other routes (non-rigid Hurwitz, number-field base change, or different geometry).

Status: **locked negative control**.
