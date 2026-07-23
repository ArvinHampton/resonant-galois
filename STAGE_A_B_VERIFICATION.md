# Stage A + B Arithmetic Verification

**Date:** 23 July 2026  
**Status:** Confirmed (computational)

This document records independent verification of the pure-even multi-k arithmetic core of Resonant Algebra.

## Summary

| Mode | Stage | Result |
|------|-------|--------|
| DIG | A1 core identities | PASS |
| GROW | A3 beyond A5 | PASS |
| BUILD | B outsider-checkable tables | PASS |
| Arithmetic grounding of pure-even multi-k | A + B | CONFIRMED |

## DIG (Stage A1) — Core identities

Verified:

- Bring–Jerrard discriminant identity holds numerically.
- Homogenisation lemma: `disc(f_t) = t^{20} · disc(seed)` holds.
- General fixed-k slice has identically square discriminant.
- Catalogue seeds: 21/21 pure-even and Gal = A5 (in tested set).
- Fixed-k slices: zero evenness failures across sampled points.
- Cross-k paths (flagship ↔ classical, flagship ↔ LSW, classical ↔ LSW): pure-even and recover multiple catalogue ratios.

Sample Hilbert rates on fixed-k slices (excerpt):

| k | Sample points | Irreducible | A5 (checked) | Even failures |
|---|---------------|-------------|--------------|---------------|
| −4 | 60 | 60 | ~32 | 0 |
| 4 | 60 | 60 | ~32 | 0 |
| −8/5 | 12 | 12 | 12 | 0 |
| 4/5 | 12 | 12 | 12 | 0 |
| −12/5 | 12 | 12 | 12 | 0 |

## GROW (Stage A3) — Beyond A5

- Explicit even irreducible degree-6 polynomials obtained.
- Two concrete examples with Gal = A6:  
  `x^6 − 6x^2 ± 6x + 2`
- Envelope construction continues to harvest A5 seeds across multiple ratios.

## BUILD (Stage B) — Outsider-checkable material

- Irreducibility rates along fixed-k slices (tables).
- Galois histograms on cross-k paths (tables).
- Re-confirmation of the rigid φ obstruction: disc = 5 · (square).
- Thin even A6 examples.

Machine-readable companion data: `build/RNT_STAGE_B_DATA.json` (local package).

## Interpretation (strict)

1. The pure-even multi-k arithmetic theory is computationally verified and strengthened.
2. The method is not limited to A5; concrete A6 examples exist.
3. Rate tables and identities are available for independent checking.
4. Geometric multi-k fusion and structural necessity theorems remain open.
5. No physical or experimental claims are made or supported by these data.

## Programme stance

The citable mathematical centre remains the pure-even multi-k theory.  
This verification constitutes Stage A + Stage B arithmetic grounding of that core.  
It does not close the geometric fusion gap or establish a broader empirical (physical) foundation for Resonant Number Theory.
