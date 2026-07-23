# resonant-galois

**Constructive arithmetic generation of pure-even A₅ families from a ternary / resonant lattice**

This repository documents an experimental but systematic contribution to constructive Inverse Galois theory.

## Core result

An explicit arithmetic theory of pure-even Bring–Jerrard quintics

```
x⁵ + α x + β ∈ ℚ[x]
```

organised by the ratio invariant

```
k = β / α.
```

For every fixed rational k ≠ 0 the family

```
α(m) = 256 m² - (3125 k⁴)/256 ,   β(m) = k · α(m)
```

has identically square discriminant. Specialisations at rational m produce infinitely many A₅ polynomials (and occasionally D₅). Several distinct values of k each carry multiple catalogue seeds, yielding multi-seed pure-even slices. A two-parameter envelope and explicit cross-k paths join seeds belonging to different ratio classes while remaining pure-even over ℚ(t).

The numerical lattice that supplies the seeds is the same ternary / resonant lattice that appears in the HQCC and “9 Maths of Unification” arithmetic (order-3 structure, three-generation data, Ad_SO(3) branching).

## What this is

- A finished generative arithmetic construction of pure-even A₅ families over ℚ.
- An organisation of those families by a simple rational invariant k.
- Explicit multi-seed and multi-k pure-even paths.
- A carefully documented experimental pipeline that also produced A₆ polynomials from 6×6 template enlargements.

## What this is not

- Not a claim of priority for realising A₅ or A₆ (both groups have many known realisations).
- Not a general solution of the Inverse Galois Problem.
- Not a conceptual theorem that the ternary structure *must* produce alternating monodromy.
- Not a verification of any physical claims associated with the source arithmetic or with the larger S²-11DM²ET-X framework.

## Status (July 2026)

| Layer | Status |
|-------|--------|
| Pure-even fixed-k slices | Finished |
| Multi-seed pure-even families | Finished (10 distinct k) |
| 2-parameter envelope + cross-k paths | Finished |
| Early matrix-template catalogues (A₅ / A₆) | Finished |
| Geometric multi-k (Nielsen-labelled Hurwitz family) | Open |
| Structural necessity theorem | Open |

## Key files

- `ARITHMETIC_MULTI_K.md` — fixed-k slices, envelope, cross-k paths, and the enlarged seed catalogue
- `CATALOGUE.md` — early matrix-template A₅ (36) and A₆ (4) lists
- `FUSION_GAP.md` — precise statement of the remaining geometric fusion problem
- `RESOLUTION.md` — three conceptual criteria that would turn the experiment into a theorem
- `EVENNESS_OBSTRUCTION.md` — explicit counter-examples showing structural axioms alone do not force square discriminant
- `RESONANT_ALGEBRA.md` — implications for resonant algebra and the generative use of the ternary lattice

## Pipeline (arithmetic)

```
ternary / resonant lattice seeds
        ↓
Bring–Jerrard form x⁵ + α x + β
        ↓
ratio k = β/α
        ↓
pure-even 1-parameter family for each k
        ↓
envelope / cross-k paths → multi-k pure-even families over ℚ(t)
```

## License

MIT
