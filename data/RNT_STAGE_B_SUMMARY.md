# Stage B Data Summary

Outsider-checkable arithmetic statements produced by the Stage A+B verification run (23 July 2026).

## Fixed-k slice rates (sampled)

- k = −4: 60 integer points tested → 60 irreducible, 0 evenness failures, substantial A5 yield.
- k = +4: same pattern.
- k = −8/5 (flagship class): 12 points, 12 irreducible, 12 A5, 0 evenness failures.
- k = +4/5 (classical class): 12 points, 12 irreducible, 12 A5, 0 evenness failures.
- Additional multi-seed ratios (±12/5, ±16/5, …) likewise show zero evenness failures in sample.

## Cross-k paths

Three explicit paths (flagship–classical, flagship–LSW, classical–LSW) remain pure-even over the parameter and recover seeds belonging to more than one catalogue ratio k.

## Degree 6

At least two monic sextics with square discriminant and Gal = A6 were obtained:

- x⁶ − 6x² + 6x + 2
- x⁶ − 6x² − 6x + 2

## Rigid obstruction (re-confirmed)

For the preferred rigid cover φ of passport (3A,3A,5A),

disc(monic(φ(y) − t)) = 5 · (square)

for rational t ∉ {0,1}. Consequently there are no even irreducible specialisations over Q.

## Note on machine-readable data

A full JSON dump of the sampled points, discriminants and Galois labels can be added in a subsequent commit as `data/RNT_STAGE_B_DATA.json` when the serialisation is finalised.
