# Direction 2 — Resonant monoid / saturation (locked)

Runner: `l0_monoid_pe_b.py` (shared with Dir 3).

---

## Objects (bound ≤ 30k)

| Object | Size | Meaning |
|--------|------|--------|
| M₀ = ⟨3, 61, 80, 243, 539⟩_× | 32 | Pure multiplicative resonant monoid |
| Resonant-smooth products | 32 | Matches M₀ |
| M₂ (+ short additive closure) | ~4504 | Catalogue-style enrichment |

---

## Findings

1. **B-embed on monoid sample:** irr ≈ 0.99 (disc □ free by identity).
2. **PE cleared |α| ∈ M₂:** 19/42 — monoid organises parameters (k, A, t) more tightly than raw (α, β).
3. **Residue support mod 3ᵐ:** proper thin subset (density proxy < 1).
4. **No empty irr classes mod 9** for A ∈ [−200, 200] ∖ {0}.

---

## Reading

L₀ is not an arbitrary search pool: the multiplicative monoid M₀ is tiny and exact; additive enrichment M₂ recovers catalogue scale. Residue thinness is a first saturation-type constraint.

Not a necessity claim. Structure theory of the lattice as a monoid.

See also: `L0_SECONDARY_INVARIANTS.md`, `TERNARY_LATTICE_DIRECTIONS.md`.
