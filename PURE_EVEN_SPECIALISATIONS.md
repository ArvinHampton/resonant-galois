# Pure-even specialisations — positive control (locked)

## Two-sided contrast

| Side | Example | Disc | Parity | Gal |
|------|---------|------|--------|-----|
| Pure-even | x⁵ − 55x + 88 | □ | even | A₅ |
| Rigid t=3 | monic(φ − 3) | 5·□ | odd | S₅ |

## Fresh even-side data (representative)

| Block | Result |
|-------|--------|
| k-slices | Dense Z-coeff points; even_fail = 0 on tested lattice; A₅ among Gal checks |
| Flagship k = −8/5 | Samples e.g. (−55,88), (320,−512), (820,−1312), … |
| LSW k = −4 | Dense A₅ hits; recovers (−100,400), (124,−496), … |
| Classical k = 4/5 | Recovers (20,16), (95,76), … |
| Cross-k paths | disc identity True; multi-catalogue k on flag↔classical, flag↔LSW, classical↔LSW |
| Homogenisation | Seed families sample-even; classical lemma proved |

## Sample path specialisations

- flag ↔ classical: (−55,88) → (20,16) (and k = −4/5 midpoint)
- flag ↔ LSW: (−55,88) → (−100,400)
- classical ↔ LSW: (20,16) → (−100,400)

## Independent checks (external)

- Flagship x⁵ − 55x + 88: Gal = A₅, disc = 242000² (square).
- x⁵ − 100x + 400: Gal = A₅.
- s=2 fibre x⁵ + 20x − 32: D₅ (even monodromy includes dihedral; consistent).

## Status

Pure-even arithmetic multi-k: **PASS / locked**.  
Geometric Nielsen multi-k: still open (see `FUSION_GAP.md`, `EXPLICIT_3A4_EQUATION.md`).
