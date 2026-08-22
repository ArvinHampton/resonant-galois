# Research roadmap — post-publish stance

_Do not reopen settled layers. Citable centre = pure-even multi-k + four-face organising principle. Necessity paused. Canonical T3 production-locked._

---

## Status baseline (do not reopen)

| Layer | Status |
|-------|--------|
| Pure-even multi-k arithmetic | **Finished + published core** (`PURE_EVEN_MULTI_K.md`) |
| Four-face organising principle | **Locked** (`TERNARY_ORGANIZING_PRINCIPLE.md`) |
| Stage D density / height | **Locked** (`STAGE_D_DENSITY.md`) |
| Matrix templates + evenness obstruction | **Documented** (`HQCC_MATRIX_TEMPLATES.md`) |
| Tier 1.1 identical-square subclass | **Deepened / locked** (`TIER11_DEEPEN.md`) — no Crit-2 fragment |
| Tier 1.2 Candidate C functor | **First cut done** (`CANDIDATE_C_FUNCTOR.md`) |
| Tier 1.3 / blowup genus of P | **g=1** (`GENUS_P_BLOWUP.md`) |
| Criterion 3 sign character | **Deepened** — no forcing χ (`CRITERION3_DEEPEN.md`) |
| Arboreal T₃ vs catalogue | **Probed** (`ARBOREAL_T3.md`) |
| **Artin conductor / ramification census** | **Done** (`ARTIN_CONDUCTOR_RAMIFICATION.md`) — support + Frob proxy; Swan open |
| **Canonical T3** | **Production lock** |
| **Necessity theorem** | **Paused** (`NECESSITY_THEOREM.md`) |
| **Geometric multi-k fusion** | **Principal open problem** |
| Physical S²-11DM²ET-X claims | **Separate** (not Galois inputs) |

**Slogan:** generative success ≠ forced alternating monodromy from HQCC axioms.

---

## Citable package

1. **`PURE_EVEN_MULTI_K.md`** — theorem-grade BJ / pure-even / envelope / paths; HQCC as lattice only.
2. **`TERNARY_ORGANIZING_PRINCIPLE.md`** — four faces; structural reading.
3. **`THEOREMS.md`** — ledger.
4. **`ARTIN_CONDUCTOR_RAMIFICATION.md`** — disc support, persistent vs moving primes, Chebotarev proxy.
5. Verification: `lib/lemmas.py`, `pure_even_specialisations.py`, `stage_d_density.py`.
6. Negative control: `RIGID_FIBRE_T3.md`.

---

## Artin conductor census (locked 2026-08-22)

| Family | Persistent primes | Moving primes |
|--------|-------------------|---------------|
| PE / flagship seed | {2,5} or {2,5,11} | seed coeffs |
| Mestre P_t | {2,5,11} | Q(t) in disc = C·Q(t)² |
| B-avatar | {2,3} | A and A²+84375 |

- All census discs perfect squares.
- Frob samples: only even types; consistent with A₅.
- Full Swan/Artin f_p exponents: **not** computed.
- M₀ organises parameters, not full conductor support.

Optional follow-on: Magma/Sage conductor(); inertia at 2,3,5; radical(disc) vs height.

---

## Active priorities

### Tier G — Geometric multi-k fusion (principal open)

See `GEOMETRIC_MULTI_K_FUSION.md`, G1–G3 docs. Necessity stays orthogonal.

### Tier A — Strengthen published centre

| # | Task | Status |
|---|------|--------|
| A1 | Pure-even multi-k core | **Done** |
| A2 | L₀ secondary invariants | **Done** |
| A3 | Outsider verification package | Partial |
| A4 | Artin conductor support census | **Done** |

### Tier L — Ternary lattice

| Dir | Status |
|----:|--------|
| 1–4 | **Done** (invariants, monoid, PE↔B first cut, Mestre orbit) |
| 5 | Necessity avatar — paused |

### Tier A+ — Generative enlargements

Mestre flagship P_t, B-embed lattice, evenness avatars PE+B — **executed and review-checked**.

### Tier C — Necessity

**Paused.** Resume only with HQCC-native naming without classical evenness ansätze.

---

## What not to do

| Do not | Why |
|--------|-----|
| Re-litigate pure-even needing HQCC for disc identity | Settled: no |
| Surgery on rigid φ/Q for even fibres | Ruled out |
| Treat necessity as centre homework | Paused |
| Claim full Artin f(ρ) from this census | Support only |
| Re-run same linear cuts of T | Exhausted |

---

## Priority table (now)

| Track | Priority |
|-------|:--------:|
| Cite pure-even + organising principle | **Centre** |
| Geometric multi-k fusion | Principal open |
| Outsider verification package | High (support) |
| Magma full conductors (optional) | Low |
| Necessity Crit 1–3 | **Paused** |
