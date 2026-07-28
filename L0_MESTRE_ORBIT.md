# Direction 4 — Mestre orbit of the lattice \(L_0\)

_Elapsed: 4.55s_

**Verdict:** L0 Mestre orbit (4.55s). seeds=24, with R-space=24, families disc□ in Q(t)=24, sample A5 specs=35. Even lattice-t pairs=480/480. Graph nodes=528 edges=504. Second hops=3. Lattice t remains a stable parameter set under Mestre deformation of L0 seeds.

Context: `TERNARY_LATTICE_DIRECTIONS.md`. Necessity paused.

---

## Setup

For each seed \(P\) (PE multi-seed representatives + B-embed \(P_A\), \(A\in L_0\)):
1. Solve \(P''R-2P'R'\equiv 0\pmod{P}\), \(\deg R<\deg P\).
2. Form \(P_t=\operatorname{Res}_y(P(y),z-y-t R(y))\).
3. Specialise at lattice \(t\in L_0\); record disc□ / Gal.
4. Optional second Mestre hop on an \(A_5\) specialisation.

- Seeds: **24**
- With nontrivial \(R\): **24**
- Families with disc□ in \(\mathbb{Q}(t)\): **24**
- Sample \(A_5\) specialisations: **35**
- Graph: **528** nodes, **504** edges

### By kind

| kind | seeds | with R | fam disc□ | A5 specs |
|------|------:|-------:|----------:|---------:|
| PE | 9 | 9 | 9 | 35 |
| B | 15 | 15 | 15 | 0 |

---

## Orbit table (per seed)

| seed | kind | dim R | fam □ | t=0 ok | A5 | even specs | R |
|------|------|------:|:-----:|:------:|---:|-----------:|---|
| flagship_explicit | PE | 1 | True | True | 10 | 20 | `x**4 + 8*x**3 - 32*x**2 + 33` |
| classical_20_16 | PE | 1 | True | True | 10 | 20 | `x**4 - 4*x**3 - 8*x**2 - 12` |
| classical_95_76 | PE | 1 | True | True | 10 | 20 | `x**4 - 4*x**3 - 8*x**2 - 57` |
| lsw_100_400 | PE | 1 | True | True | 5 | 20 | `x**4 + 20*x**3 - 200*x**2 + 60` |
| flagship | PE | 1 | True | True | 0 | 20 | `x**4 + 8*x**3 - 32*x**2 - 192` |
| lsw | PE | 1 | True | True | 0 | 20 | `5*x**4 + 100*x**3 - 1000*x**2 + 9363` |
| lsw_m | PE | 1 | True | True | 0 | 20 | `5*x**4 - 100*x**3 - 1000*x**2 + 9363` |
| s12 | PE | 1 | True | True | 0 | 20 | `x**4 + 12*x**3 - 72*x**2 + 183` |
| s16 | PE | 1 | True | True | 0 | 20 | `x**4 + 16*x**3 - 128*x**2 + 708` |
| B_A=3 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=9 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=27 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=61 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=80 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=243 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=539 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=55 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=88 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=95 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=-3 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=-61 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=-80 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=18 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |
| B_A=54 | B | 1 | True | True | 0 | 20 | `x**3 + 15*x` |

---

## Lattice-\(t\) specialisations (sample rows with Gal)

| seed | \(t\) | disc□ | Gal | lattice-like coeffs |
|------|----:|:-----:|-----|:-------------------:|
| flagship_explicit | 0 | True | A5 | True |
| flagship_explicit | 1 | True | A5 | False |
| flagship_explicit | -1 | True | A5 | False |
| flagship_explicit | 2 | True | A5 | False |
| flagship_explicit | 3 | True | A5 | False |
| flagship_explicit | 5 | True | A5 | False |
| flagship_explicit | 9 | True | A5 | False |
| flagship_explicit | 27 | True | A5 | False |
| flagship_explicit | 61 | True | A5 | False |
| flagship_explicit | 80 | True | A5 | False |
| classical_20_16 | 0 | True | A5 | True |
| classical_20_16 | 1 | True | A5 | False |
| classical_20_16 | -1 | True | A5 | False |
| classical_20_16 | 2 | True | A5 | False |
| classical_20_16 | 3 | True | A5 | False |
| classical_20_16 | 5 | True | A5 | False |
| classical_20_16 | 9 | True | A5 | False |
| classical_20_16 | 27 | True | A5 | False |
| classical_20_16 | 61 | True | A5 | False |
| classical_20_16 | 80 | True | A5 | False |
| classical_95_76 | 0 | True | A5 | False |
| classical_95_76 | 1 | True | A5 | False |
| classical_95_76 | -1 | True | A5 | False |
| classical_95_76 | 2 | True | A5 | False |
| classical_95_76 | 3 | True | A5 | False |
| classical_95_76 | 5 | True | A5 | False |
| classical_95_76 | 9 | True | A5 | False |
| classical_95_76 | 27 | True | A5 | False |
| classical_95_76 | 61 | True | A5 | False |
| classical_95_76 | 80 | True | A5 | False |
| lsw_100_400 | 0 | True | A5 | True |
| lsw_100_400 | 1 | True | A5 | False |
| lsw_100_400 | -1 | True | A5 | False |
| lsw_100_400 | 2 | True | A5 | False |
| lsw_100_400 | 3 | True | A5 | False |

---

## Graph summary

- Nodes: **528** (seed / family / specialisation)
- Edges: **504**
- Gal labels on specialisation edges: `{'A5': 35, 'even_unchecked': 445}`

Structure:
```
seed P  --Mestre-->  family P_t  --t∈L0-->  specialisation P_t0
```

## Second Mestre hops

| parent | t | dim R₂ | disc□ at t₂=1 | status |
|--------|--:|-------:|:-------------:|--------|
| flagship_explicit | 1 | 1 | True | None |
| flagship_explicit | 2 | 1 | True | None |
| flagship_explicit | 3 | 1 | True | None |

---

## Lattice stability

- Specialisation pairs \((\mathrm{seed},t)\): **480**
- With disc□: **480**
- Parameter \(t\) drawn from \(L_0\) **stays a valid lattice parameter** for the Mestre family whenever the seed is even (disc□ identity in \(t\)).
- Coefficient vectors of \(P_t\) for \(t\\neq 0\) are **typically not** elementwise in \(L_0\) (Mestre mixes degrees) — stability is at the **parameter** level, not the raw coefficient monoid.

## Conclusions (Dir 4)

1. Every tested even PE / B seed on \(L_0\) has \(\dim R\ge 1\) Mestre space.
2. `shift_y_tR` families inherit disc□ in \(\mathbb{Q}(t)\) from the seed.
3. Lattice \(t\) produces a systematic cloud of even (often \(A_5\)) specialisations.
4. The **Mestre orbit graph** is a well-defined generative structure on \(L_0\): seeds → 1-param families → lattice fibres.
5. Still **not** a necessity theorem: evenness enters via seed disc□ / B-identity.

## Sequence status

| Dir | Status |
|----:|--------|
| 1 Secondary invariants | Done |
| 2 Resonant monoid | Done |
| 3 PE ↔ B unify | First cut |
| 4 Mestre orbit | **Done** (this file) |
| 5 Necessity avatar | Paused |

```bash
python l0_mestre_orbit.py
```

_Generated by l0_mestre_orbit.py_