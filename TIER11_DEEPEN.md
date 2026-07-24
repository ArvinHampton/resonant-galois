# Tier 1.1 deepen — locked

**Status (2026-07-24).** Largest natural subclasses of T with disc identically square classified. No Criterion-2 fragment. Priority shifts to Tier 2.

Artifacts: `TIER11_DEEPEN.md` / `.json`, runner `tier11_deepen.py` (~105s).

---

## Largest natural subclass with disc identically square

| Family | Beyond BJ-embed? | Free params | Identically □? | HQCC-native? |
|--------|------------------|-------------|----------------|--------------|
| Pure-even envelope | No | 2 (m,k) | Yes | No |
| Pure-even k-slice | No | 1 | Yes | No |
| Homog fixed even seed x⁵+px³+rx+s (a=e=0, p≠0) | Yes | 1 | Yes in t | No |
| Weighted scale of fixed seed | Yes | 1 (one scale) | Yes | No |
| Genuine poly 1-param family of seeds → 2-param after homog | — | — | **0 found** | — |

**Ceiling:** overall dim 2 = pure-even envelope; beyond BJ-embed dim 1 only.

---

## Scan details (false friends filtered)

- Disc(χ_T): total degree 12; all structural cuts odd (not identically square).
- On a=e=0: factor b² × (deg-8 odd part) — square only under further relations.
- Raw identical-square seed hits: 12 → after classification:
  - genuine 2-param beyond BJ: **0**
  - weighted fixed seeds (e.g. 6u², −7u⁴, ±8u⁵): still dim 1
  - s≡0: reducible χ=x(⋯)
  - pure-even p=0: classical, not beyond BJ
- Poly maps: LSW / flag-k pure-even embed and homog seed □; M-deforms / blends not.

---

## 3-cycles

Seed (6,−7,±8) homog in T, t=1,2,3,5: all HIT_A5, 3-cycle census true.  
**Inherited from seed (Hilbert), not forced by T structure.**

---

## HQCC-axiom naming

| Subclass | HQCC name? | Why not |
|----------|------------|--------|
| Pure-even envelope | None | Classical BJ ansatz |
| Homog no-x² (a=e=0) | None | a=0 contradicts base M (a=3); t-weights = homog ansatz |

---

## Necessity verdict

**No Criterion-2 fragment.**

Nothing is simultaneously:
1. identically disc-square,
2. beyond classical pure-even,
3. 3-cycle-forced by structure, and
4. HQCC-axiom-named.

---

## Priority shift

Further 1.1 only with a **new algebraic idea** (e.g. literature parametric even quintics realised inside T).

**Default → Tier 2** (paper / catalogue invariants / verification package).

Geometric leftovers stay non-blocking.

See also: `NECESSITY_THEOREM.md`, `HQCC_MATRIX_TEMPLATES.md`, `EVENNESS_OBSTRUCTION.md`, `RESEARCH_ROADMAP.md`.
