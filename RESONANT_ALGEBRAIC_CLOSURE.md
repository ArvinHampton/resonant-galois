# Resonant algebraic closure — locked

**Status (2026-07-24)**

Four candidates for “resonant algebraic closure.” Only A is the operational closure behind the finished pure-even multi-k theory. C is the design mirror of classical binary. B and D do not force monodromy necessity.

---

## Four candidates

| Candidate | Content | What it closes | Status |
|-----------|---------|----------------|--------|
| **A. Lattice + pure-even envelope** | L₀ → combos → k → envelope + cross-k paths | (α,β) under pure-even ops | **In use; multi-k finished** |
| **B. Cyclotomic ℛ** | Adjoin 2cos(2π/539) | Coefficient field | Defined; no monodromy necessity |
| **C. Mod-2 ↔ mod-3 design mirror** | Same branched-contraction scheme, p=3 | Design analogy | Clear; **functor missing** |
| **D. Homotopy / cobordism 539** | HQCC narrative | — | Not outsider-checkable |

Evenness never needs B–D: it is BJ geometry on the output of A.

---

## Definitions worth adopting

**Resonant lattice closure (A).**  
Smallest set of (α,β) obtained from the resonant / ternary lattice L₀ = {3, 9, 27, 61, 80, 243, 539, …} and stable under:
1. integer combinations,
2. formation of k = β/α,
3. pure-even envelope α(m) = 256m² − (3125 k⁴)/256, β = k·α,
4. cross-k paths in (m,k)-space.

Closed, explicit, and identical with the finished multi-k machine.

**Resonant design mirror (C).**  
A mod-3 branched map with contraction on 0 mod 3, expansion on 1,2 mod 3, and coefficient lattice generated from order-3 data. HQCC/T₃ fits this definition; classical Collatz is the mod-2 twin. This is a design analogy, **not** a constructed functor.

**Non-definition (D).**  
“Topological resonant algebraic closure” remains undefined until a space and boundary maps exist that let an outsider recompute 539.

---

## Relation to necessity

To force alternating monodromy one still needs Criterion 1–3:
- a canonical HQCC object with proved alternating monodromy,
- structural matrix axioms ⇒ disc □ + 3-cycles always, or
- an HQCC-invariant sign character.

Candidate A does not supply necessity. B–D do not either (yet).

---

## Bottom line

| Claim | Verdict |
|-------|--------|
| Resonant algebraic closure (data + finished theory) | **Candidate A** |
| Mirror of classical binary | **Candidate C** (design only) |
| Generative success = forced Aₙ from HQCC axioms | **False** |

**Next depth (optional):** write an explicit functor for C — binary height/valuation data → ternary lattice element or template T(a,…,f) — then test even monodromy on the image.

See also: `NECESSITY_THEOREM.md`, `ARITHMETIC_MULTI_K.md`, `HQCC_MATRIX_TEMPLATES.md`.
