# Resonant Galois Programme

Structural constructions linking 9 Maths / HQCC ternary–flux arithmetic
to explicit polynomials with Galois groups \(A_5\), \(A_6\), and a finished
**pure-even multi-\(k\)** arithmetic centre.

## Citable centre (start here)

| Doc | Role |
|-----|------|
| **`PURE_EVEN_MULTI_K.md`** | Theorem-grade pure-even multi-\(k\) core (HQCC as lattice only) |
| **`TERNARY_ORGANIZING_PRINCIPLE.md`** | Four-face structural reading of generative success |
| **`THEOREMS.md`** | Full theorem ledger |
| **`GEOMETRIC_MULTI_K_FUSION.md`** | **Principal open problem** — Nielsen \(A_5\) family + lattice Hilbert recovery |
| **`G1_3A4_TRIPLE_ROOT.md`** | G1 first cut: 3A⁴ triple-root elim + multi-seed Hilbert (no geometric multi-\(k\) yet) |
| **`G1_SEED_TSCHIRNHAUS.md`** | G1 seed-first + quadratic Tschirnhaus (still no catalogue locus hits) |
| **`G1_PARAM_FIELD_RESOLVENT.md`** | G1 param-field resolvent: \(f\in K(s)[x]\), Norm, multi-sheet Hilbert re-test |
| **`G2_NIELSEN_G0.md`** | G2: \(2A\,3A^3\) / \(2A^2\,3A^2\) covers + multi-seed Hilbert (no cat hits) |
| **`G3_ENVELOPE_MONODROMY.md`** | G3: pure-even envelope monodromy = 5-cycle type on multi-k paths |
| **`G3B_5A5B_BRAID_LIFT.md`** | G3b: **5A⁴** lock; braid orbit 1×10; lift inv +1 |
| **`EXPLICIT_5A4_EQUATION.md`** | Explicit \(\mathbb{Q}(t)\) pure-even model + catalogue Hilbert |
| **`EXPLICIT_5A4_REFINEMENTS.md`** | Common-base monodromy **5A²5B²** lift−1; genus; λ/j chart |
| **`NECESSITY_THEOREM.md`** | Forced \(A_n\) from axioms — **open / paused** |
| **`REPORT.md`** | Status report |
| **`RESEARCH_ROADMAP.md`** | Priorities after publish (Tier G first) |

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

## Necessity criteria (open research — not the finished centre)

| # | Criterion | Module | Doc |
|---|-----------|--------|-----|
| 1 | Canonical HQCC monodromy object | `criterion1_hqcc.py` | `build/CRITERION1_HQCC.md` |
| 2 | Axioms ⇒ disc² + 3-cycles | `criterion2_axioms.py` / `tier11_deepen.py` | `TIER11_DEEPEN.md` |
| 3 | Sign character / ternary invariant | `criterion3_deepen.py` | `CRITERION3_DEEPEN.md` |

Master synthesis: **`RESOLUTION.md`**.  
Theorem ledger: **`THEOREMS.md`**.  
HQCC seed definition: **`HQCC_SEED.md`**.  
Geometric: **`GEOMETRIC_COVER.md`**, genus of \(P\): **`GENUS_P_BLOWUP.md`**.  
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
