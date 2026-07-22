# resonant-galois

**Constructive generation of explicit polynomials with Galois groups A₅ and A₆ from ternary-matrix templates**

This repository documents an experimental arithmetic construction in Inverse Galois theory.

## What this is

A systematic method that produces monic irreducible polynomials over ℚ whose Galois groups are the alternating groups A₅ and A₆. The method uses matrix templates that embed an order-3 / ternary block structure, followed by classical filters (irreducibility, square discriminant, cycle-type census) and computational group identification.

**Catalogues**
- 36 unique polynomials with Gal = A₅
- 4 unique polynomials with Gal = A₆

The templates and numerical lattice are motivated by ternary arithmetic arising in the "9 Maths of Unification" / resonant algebra (HQCC qutrit structure, three-generation data, Ad_SO(3) branching).

## What this is not

- Not a claim of priority for realizing A₅ or A₆ (both groups have many known realizations).
- Not a general solution of the Inverse Galois Problem.
- Not a verification of any physical claims associated with the source arithmetic.
- Not a conceptual theorem that the ternary structure *must* produce alternating groups; it is a successful generative experiment whose evenness obstruction is documented.

## Pipeline

```
ternary / companion matrix template
        ↓
characteristic polynomial
        ↓
fast filters: irreducible?  disc square?  cycle types?
        ↓
survivors only → full Galois identification
        ↓
cross-check
```

## Key files

- `RESOLUTION.md` — status of the three conceptual criteria that would turn the experiment into a theorem
- `CATALOGUE.md` — summary of the A₅ and A₆ polynomials
- `EVENNESS_OBSTRUCTION.md` — explicit counter-examples (base matrices) showing structural axioms alone do not force square discriminant

## Status

| Layer                        | Status                          |
|-----------------------------|---------------------------------|
| Generative experiment       | Finished and documented         |
| Catalogues                  | 36 A₅ + 4 A₆                    |
| Evenness obstruction        | Explicit                        |
| Conceptual theorem          | Open                            |

## License

MIT
