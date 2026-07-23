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

## Verification status

**Stages A + B** (identities, evenness, multi-k recovery) — computationally confirmed.  
**Stage D** (density & height) — completed:

- D3 discriminant height asymptotics: **theorem**
- D1 irreducibility density on the slices: conjecture + strong computational evidence
- D2 Chebotarev cycle-type proxy: empirical histograms consistent with A₅

See `STAGE_D_DENSITY.md` and `data/STAGE_D_SUMMARY.md`.

## What this is

- A finished generative arithmetic construction of pure-even A₅ families over ℚ.
- An organisation of those families by a simple rational invariant k.
- Explicit multi-seed and multi-k pure-even paths.
- Computational verification of core identities, sample rates, and height asymptotics.
- A carefully documented experimental pipeline that also produced A₆ polynomials.

## What this is not

- Not a claim of priority for realising A₅ or A₆.
- Not a general solution of the Inverse Galois Problem.
- Not a conceptual theorem that the ternary structure *must* produce alternating monodromy.
- Not a verification of any physical claims associated with the source arithmetic or with the larger S²-11DM²ET-X framework.

## Status

| Layer | Status |
|-------|--------|
| Pure-even fixed-k slices | Finished + verified |
| Multi-seed pure-even families | Finished + verified |
| 2-parameter envelope + cross-k paths | Finished + verified |
| Stage A+B arithmetic verification | Recorded |
| Stage D density & height | Recorded (D3 theorem; D1–D2 conjecture + evidence) |
| Geometric multi-k (Nielsen-labelled) | Open |
| Structural necessity theorem | Open |

## Key files

- `ARITHMETIC_MULTI_K.md` — fixed-k slices, envelope, cross-k paths
- `STAGE_D_DENSITY.md` — density evidence and height theorem
- `RESONANT_NUMBER_THEORY.md` — Stage A+B verification report
- `data/` — outsider-checkable summaries
- `FUSION_GAP.md` — remaining geometric fusion problem
- `RESOLUTION.md` — three conceptual criteria for a theorem
- `EVENNESS_OBSTRUCTION.md` — explicit counter-examples
- `RESONANT_ALGEBRA.md` — implications for resonant algebra
- `CATALOGUE.md` — early matrix-template lists

## License

MIT
