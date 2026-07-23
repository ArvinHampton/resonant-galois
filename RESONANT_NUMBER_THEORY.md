# Resonant Number Theory — Stage A + B Arithmetic Verification

**Date:** 23 July 2026  
**Status:** Arithmetic core verified; structural and geometric gaps remain open.

## Scope of this document

This file records computational and symbolic verification of the pure-even multi-k theory that currently constitutes the mathematical centre of Resonant Algebra.  
“Empirical” here means arithmetically / computationally verified. It does **not** mean experimental or physical confirmation.

## DIG (Stage A1) — Mathematical core

| Check | Result |
|-------|--------|
| BJ discriminant identity (numeric) | match |
| Homogenisation lemma | proved |
| General k-slice disc identity | True |
| Catalogue seeds disc square | 21/21 |
| Catalogue HIT_A5 | 21/21 |
| k-slice even_fail | 0 on all sampled slices |
| Homogenisation t^{20} factor | True |
| Cross-k path disc identity | True (×3 paths) |
| Cross-k multi-catalogue k | True (×3 paths) |

### Sample Hilbert rates on fixed-k slices

| k | #Z pts | irr | A5 (checked) | even fail |
|---|--------|-----|--------------|-----------|
| −4 | 60 | 60 | ~32 | 0 |
| +4 | 60 | 60 | ~32 | 0 |
| −8/5 | 12 | 12 | 12 | 0 |
| +4/5 | 12 | 12 | 12 | 0 |
| −12/5 | 12 | 12 | 12 | 0 |
| … | … | … | … | 0 |

Cross-k paths tested: flagship ↔ classical, flagship ↔ LSW, classical ↔ LSW. All remain pure-even and recover multiple catalogue ratios.

## GROW (Stage A3) — Reach beyond A₅

| Item | Result |
|------|--------|
| Degree-6 even irreducible polynomials | many |
| Explicit Gal = A₆ | 2 examples: x⁶ − 6x² ± 6x + 2 |
| Envelope A5 harvest | 50 seeds across 10 ratios k |

The method is not confined to degree 5.

## BUILD (Stage B) — Outsider-checkable statements

| ID | Content | Status |
|----|---------|--------|
| B1 | Irreducibility rates along k-slices | empirical tables |
| B2 | Galois histograms on cross-k paths | empirical tables |
| B3 | Rigid φ disc = 5 · (square) | proved (re-checked) |
| B4 | Thin even A₆ examples | empirical |

Machine-readable companion: `data/RNT_STAGE_B_DATA.json` (to be added when serialised).

## Interpretation

1. The pure-even multi-k arithmetic theory is computationally grounded: identities hold, catalogue seeds are even + A₅, fixed-k slices show no evenness failures in sample, cross-k paths are multi-k.
2. The generative method produces alternating groups beyond A₅ (explicit A₆ examples with square discriminant).
3. Stage B supplies rates, histograms and a proved obstruction that independent workers can re-verify.
4. Geometric multi-k fusion and the three structural resolution criteria remain open and are **not** claimed by this verification.

## Programme stance (unchanged)

- Citable mathematical centre = pure-even multi-k theory + verification data.
- Geometric fusion = open research problem.
- Physical / experimental claims = outside the scope of this repository and of the present verification.

See also: `ARITHMETIC_MULTI_K.md`, `FUSION_GAP.md`, `RESOLUTION.md`, `RESONANT_ALGEBRA.md`.
