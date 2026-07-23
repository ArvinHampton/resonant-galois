# Arithmetic Multi-k Pure-Even Theory

## Fixed-k pure-even families

For any fixed rational k ≠ 0 the Bring–Jerrard family

```
α(m) = 256 m² − (3125 k⁴)/256
β(m) = k · α(m)
f_m(x) = x⁵ + α(m) x + β(m)
```

has discriminant that is an identical square in m. Consequently every irreducible specialisation has even monodromy. Computational sampling shows that a positive-density set of rational m yields Gal = A₅.

## Multi-seed slices (enlarged catalogue)

An exhaustive lattice scan produced ≈60 distinct A₅ seeds in Bring–Jerrard form. Grouping by the ratio k = β/α yields 16 distinct ratios, of which 10 are multi-seed:

| k | # A₅ seeds | Notable members |
|---|------------|-----------------|
| −4 (LSW class) | 11 | (−100, 400), (124, −496), … |
| +4 | 11 | sign-flips of LSW |
| −8/5 (flagship) | 4 | (−55, 88), (145, −232), (320, −512), … |
| +8/5 | 4 | flagship flips |
| 4/5 (classical) | 4 | (20, 16), (95, 76), … |
| −4/5 | 4 | classical flips |
| ±12/5 | 5 each | |
| ±16/5 | 3 each | |

Single-seed ratios also appear (e.g. ±28/5, ±12, ±8).

## Two-parameter envelope

```
α(m,s) = 256 m² − (3125 s⁴)/256
β(m,s) = s · α(m,s)
```

- Pure-even over ℚ(m,s)
- Freezing s = k recovers every fixed-k slice
- All listed multi-seed catalogue points are recovered

## Cross-k paths

Any two envelope seeds (even with distinct k) lie on an explicit pure-even 1-parameter family obtained by a rational path in (m,k)-space, for example the linear path

```
m(u) = (1−u) m₁ + u m₂
k(u) = (1−u) k₁ + u k₂
```

Verified examples include flagship (k = −8/5) ↔ classical (k = 4/5) and flagship ↔ LSW (k = −4). Discriminants remain squares along these paths.

## Relation to earlier generative experiment

The matrix-template catalogues (36 A₅ + 4 A₆) remain valid and are retained for historical and comparative purposes. The pure-even BJ theory is strictly stronger: it supplies infinite families rather than finite lists, and organises them by a transparent arithmetic invariant.

## What is still open

A pure geometric (Nielsen-labelled) A₅ family whose Hilbert specialisations recover seeds belonging to more than one fixed-k class. The arithmetic multi-k theory is complete; the geometric multi-k identification is not.
