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
α(m) = 256 m² − (3125 k⁴)/256 ,   β(m) = k · α(m)
```

has identically square discriminant. Specialisations at rational m produce infinitely many A₅ polynomials (and occasionally D₅). Several distinct values of k each carry multiple catalogue seeds, yielding multi-seed pure-even slices. A two-parameter envelope and explicit cross-k paths join seeds belonging to different ratio classes while remaining pure-even over ℚ(t).

The numerical lattice that supplies the seeds is the same ternary / resonant lattice that appears in the HQCC and “9 Maths of Unification” arithmetic (order-3 structure, three-generation data).

## Verification status (locked 2026-07-24)

| Layer | Status |
|-------|--------|
| Pure-even fixed-k slices | **Finished + verified** |
| Multi-seed pure-even families | **Finished + verified** |
| 2-parameter envelope + cross-k paths | **Finished + verified** |
| Stage A+B arithmetic verification | **Recorded** |
| Stage D density & height | **Recorded** (D3 theorem; D1–D2 conjecture + evidence) |
| Rigid φ/ℚ even fibres (t=3 control) | **Ruled out** (disc = 5·□ → odd → S₅) |
| Pure-even specialisations (positive control) | **PASS** (flagship, LSW, classical, cross-k paths) |
| Non-classical Rₙ track | Scaffolded (field-agnostic pure-even; cosine orbits ≠ multi-k paths) |
| Explicit 3A⁴ resolvent model | **In hand** (H^(rd) ≅ ℙ¹; eliminant P(q,w); fibre s=−1 over ℚ(√5)) |
| Geometric multi-k (Nielsen-labelled) | **Open** |
| Structural necessity theorem | **Open** |

**Stages A + B** — computationally confirmed.  
**Stage D** — completed: D3 disc-height asymptotics theorem; D1 irreducibility density conjecture + evidence; D2 Chebotarev proxy histograms consistent with A₅.

See `STAGE_D_DENSITY.md`, `RIGID_FIBRE_T3.md`, `PURE_EVEN_SPECIALISATIONS.md`, `EXPLICIT_3A4_EQUATION.md`, `NONCLASSICAL_RESONANT_FIELD.md`.

## What this is

- A finished generative arithmetic construction of pure-even A₅ families over ℚ.
- Organisation of those families by the rational invariant k = β/α.
- Explicit multi-seed and multi-k pure-even paths.
- Computational verification of core identities, sample rates, and height asymptotics.
- A documented experimental pipeline that also produced A₆ polynomials.
- Explicit geometric 3A⁴ cover equations (resolvent still multi-sheeted in s alone).

## What this is not

- Not a claim of priority for realising A₅ or A₆.
- Not a general solution of the Inverse Galois Problem.
- Not a conceptual theorem that ternary structure *must* produce alternating monodromy.
- Not a verification of any physical claims of the S²-11DM²ET-X framework.

## Key files

| File | Role |
|------|------|
| `ARITHMETIC_MULTI_K.md` | Fixed-k slices, envelope, cross-k paths |
| `STAGE_D_DENSITY.md` | Density evidence and height theorem |
| `RIGID_FIBRE_T3.md` | Negative control: monic(φ−3) → S₅ |
| `PURE_EVEN_SPECIALISATIONS.md` | Positive control: pure-even A₅ samples |
| `EXPLICIT_3A4_EQUATION.md` | 3A⁴ cover / resolvent equations |
| `NONCLASSICAL_RESONANT_FIELD.md` | Rₙ track (non-classical enrichment) |
| `FUSION_GAP.md` | Geometric multi-k still open |
| `RESONANT_NUMBER_THEORY.md` | Stage A+B report |
| `RESOLUTION.md` | Three conceptual criteria for a theorem |
| `data/` | Outsider-checkable summaries |

## License

MIT
