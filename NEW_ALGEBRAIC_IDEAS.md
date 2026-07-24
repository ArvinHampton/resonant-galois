# New algebraic ideas A–F — executed (locked)

**Status (2026-07-24).** Generative machine enlarged. Criterion 2 forcing on T remains closed. Necessity remains paused.

Runner: `new_algebraic_ideas.py` (~112s).

---

## Primary — Idea A (Mestre) HIT

For each even HQCC BJ seed, solve
```
P″R − 2P′R′ ≡ 0 (mod P),  deg R < deg P.
```
All 10 seeds give dim R-space = 1.

**Example (flagship):**
```
R = x⁴ + 8x³ − 32x² + 33
```

**Family construction that works:**
```
P_t(z) = Res_y( P(y), z − y − t R(y) )   # shift_y_tR
```
- disc identically □ in ℚ(t)
- specialisations t = ±1, 2, 3, 5, 7 → Gal A₅

Other resultants (uPp_*) also disc□ but need clearing to monic ℤ[x].

**Why new:** even families first, then ask for matrix realisation — not cuts of T.

---

## Primary companion — Idea F (embed) HIT

| Source | Embed in T? |
|--------|-------------|
| BJ seeds | Classical BJ-embed |
| Depressed Mestre specialisations | Yes (90 specs checked) |
| Non-BJ family B | Yes, with d ≠ 0 (below) |

---

## Secondary — Idea B (non-BJ deg-1) HIT

```
P_A = x⁵ + 75 x³ + A x² + 3A
disc(P_A) = 324 A² (A² + 84375)²   # identically a square
```

All 16 resonant lattice values of A tested → **HIT_A5**.

**Sparse T-realisation beyond BJ-embed:**
```
d = −75,  e = f = 0,  a = −A,  b c = 72 A
```
Examples: A = 3 ⇒ c = 216; A = 61 ⇒ c = 4392.

---

## If A+F had failed — Ideas C / D / E

| Idea | Outcome |
|------|--------|
| C — ternary companion / transfer graph / 3-cycle block | No identical-square-by-construction avatar; old T stays closed negative for Crit-2 forcing |
| D — icosahedral scan | disc□ = 0 in probe |
| E — T₃ orbit polys | Not systematic A₅ sources |

---

## Do not retry (locked)

| Approach | Why |
|----------|-----|
| More linear cuts of same T | Exhausted |
| F → T hoping disc□ → 1 | Rate 0 |
| Rigid φ/ℚ surgery | Factor 5 permanent |
| Collatz ⇒ evenness without pure-even | Composite only |

---

## Synthesis

| Claim | Status |
|-------|--------|
| Arithmetic centre (pure-even multi-k) | **Finished** |
| Crit 2 forcing on this T | **Closed** (negative) |
| Generative enlargement | **Yes** — Mestre lifts + non-BJ P_A |
| Necessity (HQCC axioms force A_n) | **Still paused** |

Embeds still use classical evenness (Mestre / disc identity), not HQCC-axiom naming alone. These enlarge the generative machine; they do not prove necessity.

---

## Natural next (optional)

1. Closed-form flagship Mestre family P_t
2. Systematic HQCC-lattice points on the B-embed b c = 72 A
3. Optional new matrix avatar only with a built-in evenness identity

One-line: **Change the equation, not recut the same variety — Mestre + non-BJ family succeed as generative tools; Crit-2 forcing and necessity stay where they were.**
