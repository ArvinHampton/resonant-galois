# Direction 4 — Mestre orbit of L₀ (locked)

Runner: `l0_mestre_orbit.py` (~4.5s)

---

## Results

| Metric | Value |
|--------|--------|
| Seeds (PE multi-k + B-embed A ∈ L₀) | 24 |
| Nontrivial Mestre R-space | 24/24 |
| Families with disc □ in ℚ(t) | 24/24 |
| Even pairs (seed, t ∈ L₀) | 480/480 |
| Sample Gal A₅ | 35 |
| Graph | 528 nodes · 504 edges |
| Second Mestre hops | 3 (all disc □ at t₂ = 1) |

---

## Structure

```
seed P  --Mestre-->  family P_t  --t ∈ L₀-->  specialisation
```

- **Lattice stability** is at the parameter t ∈ L₀ (always even when the seed is even).
- Raw coefficients of P_t for t ≠ 0 are generally **not** in the monoid — coefficient mixing is expected.

---

## Reading

Every tested lattice seed admits a Mestre deformation with identically square disc; every lattice parameter t keeps the specialised fibre even. The orbit graph is a generative object on L₀, not a necessity theorem.

See also: `MESTRE_FLAGSHIP_PT.md`, `NEW_ALGEBRAIC_IDEAS.md`, `TERNARY_LATTICE_DIRECTIONS.md`.
