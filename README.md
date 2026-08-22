# Resonant Galois Programme

Structural constructions linking 9 Maths / HQCC ternary–flux arithmetic
to explicit polynomials with Galois groups A₅, A₆, and a finished
**pure-even multi-k** arithmetic centre.

## Citable centre (start here)

| Doc | Role |
|-----|------|
| **`PURE_EVEN_MULTI_K.md`** | Theorem-grade pure-even multi-k core (HQCC as lattice only) |
| **`TERNARY_ORGANIZING_PRINCIPLE.md`** | Four-face structural reading of generative success |
| **`THEOREMS.md`** | Full theorem ledger |
| **`ARTIN_CONDUCTOR_RAMIFICATION.md`** | Disc factorisations, Frob types, persistent vs moving ramification |
| **`NECESSITY_THEOREM.md`** | Forced Aₙ from axioms — **open / paused** |
| **`REPORT.md`** | Status report |
| **`RESEARCH_ROADMAP.md`** | Priorities |

## Review package (generative extensions)

**Order:** flagship Mestre first, then B-avatar. See **`REVIEW_PACKAGE.md`**.

| Order | Doc | Check |
|------:|-----|-------|
| 1 | `MESTRE_FLAGSHIP_PT.md` | `python review_flagship_b.py` |
| 2 | `B_EMBED_LATTICE.md`, `EVENNESS_AVATAR.md` | same |

**Contamination boundary:** **`CONTAMINATION_BOUNDARY.md`** — lattice integers may be model-motivated; proofs use only arithmetic properties. No G₄/539.9 s, GW/Belle II, or 539-step dynamics as Gal inputs.

```bash
cd resonant_galois
python -c "from lib.lemmas import verify_disc_formulas; print(verify_disc_formulas())"
python pure_even_specialisations.py
python stage_d_density.py
```

## Quick full rebuild

```bash
python build_all.py
```

Outputs land in `build/` (and optional mirror `../a5_brute_results/`).

## Artin conductors and ramification (2026-08-22)

Doc: **`ARTIN_CONDUCTOR_RAMIFICATION.md`**

- All census fibres: disc perfect square (even image)
- Flagship / PE persistent support: {2,5,11}
- B-avatar persistent: {2,3}; moving primes from A and A²+84375
- Mestre P_t: persistent {2,5,11}; moving from Q(t)
- Frobenius samples: only even types (5), (3,1,1), (2,2,1), id — consistent with A₅
- Full Swan/Artin exponents: not computed (support + Chebotarev proxy only)
- Monoid M₀ organises parameters, not full conductor support

## Necessity criteria (open research — not the finished centre)

| # | Criterion | Module | Doc |
|---|-----------|--------|-----|
| 1 | Canonical HQCC monodromy object | `criterion1_hqcc.py` | `build/CRITERION1_HQCC.md` |
| 2 | Axioms ⇒ disc² + 3-cycles | `criterion2_axioms.py` / `tier11_deepen.py` | `TIER11_DEEPEN.md` |
| 3 | Sign character / ternary invariant | `criterion3_deepen.py` | `CRITERION3_DEEPEN.md` |

Master synthesis: **`RESOLUTION.md`**.  
Theorem ledger: **`THEOREMS.md`**.  
HQCC seed definition: **`HQCC_SEED.md`**.  
Geometric: **`GEOMETRIC_COVER.md`**, genus of P: **`GENUS_P_BLOWUP.md`**.  
Arboreal probe: **`ARBOREAL_T3.md`**.

## Experimental engine (already validated)

```
matrix template (ternary) → char poly → irr / disc² / cycles → galois_group
```

Catalogues: `build/CATALOGUE.md` (A5/A6 freeze).

## Layout

```
resonant_galois/
  build_all.py
  lib/common.py
  criterion1_hqcc.py
  criterion2_axioms.py
  criterion3_sign.py
  build/           # generated
  RESOLUTION.md
  README.md
```

Prior heavy scans remain in `../a5_brute_results/` (DEFORM_M, A6_T6, fingerprints, …).

## License

MIT
