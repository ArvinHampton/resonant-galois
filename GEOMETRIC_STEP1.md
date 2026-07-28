# Geometric cover — Step 1: candidate rigid tuples in \(A_5\)

_Elapsed: 2.64s_

## Goal

List conjugacy-class tuples in \(A_5\) that:
1. admit product-one representatives generating \(A_5\);
2. are (heuristically) absolutely rigid;
3. can be **labelled** by HQCC ternary / flux / period motifs.

Cover construction (Step 2) is **not** performed here.

---

## Conjugacy classes of \(A_5\)

| Class | Size | Cycle type | HQCC label |
|-------|-----:|------------|------------|
| 1A | 1 | 1+1+1+1+1 | unramified / identity |
| 2A | 15 | 2+2+1 | T-complementarity / flux involution (type 2+2+1) |
| 3A | 20 | 3+1+1 | ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1) |
| 5A | 12 | 5 | period / pentagonal sector — 5-cycle class A (use sparingly) |
| 5B | 12 | 5 (inverse class) | period / pentagonal sector — 5-cycle class B (inverse class) |

---

## Candidate signatures (ranked)

| Signature | prod-1 | gen \(A_5\) | conj. orbits (sample) | abs. rigid? | HQCC score |
|-----------|-------:|----------:|----------------------:|:-----------:|-----------:|
| `(3A,3A,5A)` | 60 | 60 | 1 | True | 18 |
| `(3A,3A,5B)` | 60 | 60 | 1 | True | 18 |
| `(2A,3A,5A)` | 60 | 60 | 1 | True | 11 |
| `(2A,3A,5B)` | 60 | 60 | 1 | True | 11 |
| `(3A,3A,3A,3A)` | 2940 | 1080 | 18 | False | 40 |
| `(3A,3A,3A,2A)` | 1920 | 1440 | 22 | False | 33 |
| `(2A,3A,3A,3A)` | 1920 | 1440 | 22 | False | 33 |
| `(2A,2A,3A,3A)` | 1560 | 1080 | 18 | False | 26 |
| `(3A,3A,2A,5A)` | 1200 | 1200 | 20 | False | 21 |
| `(3A,3A,5A,5B)` | 1020 | 1020 | 17 | False | 16 |
| `(2A,3A,5A,5B)` | 720 | 720 | 12 | False | 9 |
| `(3A,3A,3A)` | 140 | 0 | None | None | 30 |
| `(2A,2A,3A)` | 60 | 0 | None | None | 16 |

### HQCC labelling of top candidates

#### `(3A,3A,5A)`
- HQCC score: **18**
- gen \(A_5\): 60, rigid≈True
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - period / pentagonal sector — 5-cycle class A (use sparingly)
- sample cycle types: `[[(3, 1, 1), (3, 1, 1), (5,)], [(3, 1, 1), (3, 1, 1), (5,)]]`

#### `(3A,3A,5B)`
- HQCC score: **18**
- gen \(A_5\): 60, rigid≈True
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - period / pentagonal sector — 5-cycle class B (inverse class)
- sample cycle types: `[[(3, 1, 1), (3, 1, 1), (5,)], [(3, 1, 1), (3, 1, 1), (5,)]]`

#### `(2A,3A,5A)`
- HQCC score: **11**
- gen \(A_5\): 60, rigid≈True
  - T-complementarity / flux involution (type 2+2+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - period / pentagonal sector — 5-cycle class A (use sparingly)
- sample cycle types: `[[(2, 2, 1), (3, 1, 1), (5,)], [(2, 2, 1), (3, 1, 1), (5,)]]`

#### `(2A,3A,5B)`
- HQCC score: **11**
- gen \(A_5\): 60, rigid≈True
  - T-complementarity / flux involution (type 2+2+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - period / pentagonal sector — 5-cycle class B (inverse class)
- sample cycle types: `[[(2, 2, 1), (3, 1, 1), (5,)], [(2, 2, 1), (3, 1, 1), (5,)]]`

#### `(3A,3A,3A,3A)`
- HQCC score: **40**
- gen \(A_5\): 1080, rigid≈False
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
- sample cycle types: `[[(3, 1, 1), (3, 1, 1), (3, 1, 1), (3, 1, 1)], [(3, 1, 1), (3, 1, 1), (3, 1, 1), (3, 1, 1)]]`

#### `(3A,3A,3A,2A)`
- HQCC score: **33**
- gen \(A_5\): 1440, rigid≈False
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - T-complementarity / flux involution (type 2+2+1)
- sample cycle types: `[[(3, 1, 1), (3, 1, 1), (3, 1, 1), (2, 2, 1)], [(3, 1, 1), (3, 1, 1), (3, 1, 1), (2, 2, 1)]]`

#### `(2A,3A,3A,3A)`
- HQCC score: **33**
- gen \(A_5\): 1440, rigid≈False
  - T-complementarity / flux involution (type 2+2+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
- sample cycle types: `[[(2, 2, 1), (3, 1, 1), (3, 1, 1), (3, 1, 1)], [(2, 2, 1), (3, 1, 1), (3, 1, 1), (3, 1, 1)]]`

#### `(2A,2A,3A,3A)`
- HQCC score: **26**
- gen \(A_5\): 1080, rigid≈False
  - T-complementarity / flux involution (type 2+2+1)
  - T-complementarity / flux involution (type 2+2+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
  - ternary branch / Z/3 / n↦n/3 or Ad residue (type 3+1+1)
- sample cycle types: `[[(2, 2, 1), (2, 2, 1), (3, 1, 1), (3, 1, 1)], [(2, 2, 1), (2, 2, 1), (3, 1, 1), (3, 1, 1)]]`

---

## Interpretation for HQCC-nativeness

### Preferred (ternary-heavy)

- Signatures with **multiple 3A** factors match the ternary / HQCC branch motif.
- `(3A,3A,3A,2A)` and `(3A,3A,3A,3A)` maximise ternary content; rigidity must be checked carefully
  (four branch points → Hurwitz space may be positive-dimensional → families, not a single cover).

### Classical rigid (known)

- `(2A,3A,5A)` / `(2A,3A,5B)`: The signature (2A,3A,5A) is the classical icosahedral triple for A5; it is absolutely rigid in the literature. Our conjugacy-orbit heuristic should find a single orbit when sampling is complete.
- HQCC labelling is **mixed**: one ternary (3A), one flux involution (2A), one period-like 5-cycle.
- Geometrically standard (icosahedral); **nativeness** requires Step 4 (express branch data via \(\{n/3,3n\pm1\}\)),
  not only abstract cycle types.

### Tension with arithmetic programme

- Pure 5-cycle emphasis risks \(D_5\)-type specialisations (G4 Heavy lesson).
- Arithmetic HQCC **seeds** already force 3-cycles via the operational criterion;
  geometry should **not** drop the ternary class.

---

## Recommended tuples for Step 2

1. **`(3A,3A,5A)`** — absolutely rigid (heuristic); HQCC score 18.
2. **`(3A,3A,5B)`** — absolutely rigid (heuristic); HQCC score 18.
3. **`(2A,3A,5A)`** — absolutely rigid (heuristic); HQCC score 11.
4. **`(2A,3A,5B)`** — absolutely rigid (heuristic); HQCC score 11.

### Immediate Step 2 options

1. **Classical path:** realise `(2A,3A,5A)` as the icosahedral cover / Belyi map over \(\mathbb{Q}(t)\)
   (known existence); then attempt Step 4 HQCC labelling of branch points.
2. **Ternary path:** if a rigid (or 0-dim braid orbit) multi-`3A` signature exists, construct its
   Hurwitz cover preferentially — better HQCC narrative.
3. **Computational:** Magma/GAP `HurwitzClassNumber` / braid orbit routines for signatures
   marked rigid=False but gen>0 (may still have rational points).

---

## What Step 1 does *not* claim

- No geometric monodromy proof yet (needs cover in Step 2).
- No automatic compatibility with HQCC seeds \(x^5-55x+88\), etc. (Step 3).
- No closed form in \(\xi=2\cos(2\pi/539.9)\).

---

## Raw results

```json
[
  {
    "signature": [
      "3A",
      "3A",
      "5A"
    ],
    "class_sizes": [
      20,
      20,
      12
    ],
    "cartesian_total": 4800,
    "checked": 400,
    "n_product_one": 60,
    "n_generating_A5": 60,
    "n_conjugacy_orbits_sample": 1,
    "absolutely_rigid_heuristic": true,
    "hqcc_score": 18,
    "hqcc_labels": [
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "period / pentagonal sector \u2014 5-cycle class A (use sparingly)"
    ],
    "sample_tuple_cycle_types": [
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ]
    ]
  },
  {
    "signature": [
      "3A",
      "3A",
      "5B"
    ],
    "class_sizes": [
      20,
      20,
      12
    ],
    "cartesian_total": 4800,
    "checked": 400,
    "n_product_one": 60,
    "n_generating_A5": 60,
    "n_conjugacy_orbits_sample": 1,
    "absolutely_rigid_heuristic": true,
    "hqcc_score": 18,
    "hqcc_labels": [
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "period / pentagonal sector \u2014 5-cycle class B (inverse class)"
    ],
    "sample_tuple_cycle_types": [
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          3,
          1,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ]
    ]
  },
  {
    "signature": [
      "2A",
      "3A",
      "5A"
    ],
    "class_sizes": [
      15,
      20,
      12
    ],
    "cartesian_total": 3600,
    "checked": 300,
    "n_product_one": 60,
    "n_generating_A5": 60,
    "n_conjugacy_orbits_sample": 1,
    "absolutely_rigid_heuristic": true,
    "hqcc_score": 11,
    "hqcc_labels": [
      "T-complementarity / flux involution (type 2+2+1)",
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "period / pentagonal sector \u2014 5-cycle class A (use sparingly)"
    ],
    "sample_tuple_cycle_types": [
      [
        [
          2,
          2,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          2,
          2,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ],
      [
        [
          2,
          2,
          1
        ],
        [
          3,
          1,
          1
        ],
        [
          5
        ]
      ]
    ]
  },
  {
    "signature": [
      "2A",
      "3A",
      "5B"
    ],
    "class_sizes": [
      15,
      20,
      12
    ],
    "cartesian_total": 3600,
    "checked": 300,
    "n_product_one": 60,
    "n_generating_A5": 60,
    "n_conjugacy_orbits_sample": 1,
    "absolutely_rigid_heuristic": true,
    "hqcc_score": 11,
    "hqcc_labels": [
      "T-complementarity / flux involution (type 2+2+1)",
      "ternary branch / Z/3 / n\u21a6n/3 or Ad residue (type 3+1+1)",
      "period / pentagonal sector \u2014 5-cycle class B (inverse class)"
    ],
    "sample_tuple_cycle_types": [
      [
        [
          2,
          2,
          1
        ],
        [
          3,
          1
```

_Generated by geometric_step1.py_