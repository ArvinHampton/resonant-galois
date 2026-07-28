# Resonant Galois Programme — Status Report

_Package: `Desktop/resonant_galois/` · Resolution path locked_

---

## 1. Bottom line (locked)

### Citable centre (publish)

1. **Pure-even multi-\(k\)** — **`PURE_EVEN_MULTI_K.md`**: theorem-grade BJ disc, \(k\)-slices, envelope, cross-\(k\) paths, homogenisation; HQCC as **lattice motivation only**.  
2. **Four-face organising principle** — **`TERNARY_ORGANIZING_PRINCIPLE.md`**: dynamics / lattice / matrices / Galois ternary branching + pure-even sign lock. Structural reading of generative success.  
3. **Necessity** (forced \(A_n\) from HQCC axioms) — **paused** as citable claim; remains open research (`NECESSITY_THEOREM.md`).  
4. **Principal next problem** — **geometric multi-\(k\) fusion** (`GEOMETRIC_MULTI_K_FUSION.md`): one Nielsen-labelled pure \(A_5\) family recovering the pure-even / HQCC seed lattice as Hilbert specialisations.  
5. **Canonical T3** — production dynamical baseline locked; coefficient choice not reopened.

**Slogan:** generative success + organising principle ≠ forced alternating monodromy.

### Probes completed this cycle

| Probe | Doc | Result |
|-------|-----|--------|
| Crit 3 sign character | `CRITERION3_DEEPEN.md` | No rate-1 ternary/HQCC \(\chi\) on unrestricted \(T\); pure-even control 40/40 |
| Genus of \(P\) (blowup) | `GENUS_P_BLOWUP.md` | **\(g=1\)** under ordinary sings at \((1,1)\) + infinity; not genus 0 |
| Arboreal T₃ vs catalogue | `ARBOREAL_T3.md` | Evenness from pure-even identity; Frob TV≈0.04; encoding-dependent path monodromy |
| Tier 1.1 identical-square | `TIER11_DEEPEN.md` | Beyond BJ: dim 1 homog only; HQCC naming fails |

### Earlier Tier 1 (locked)

- **1.2** Candidate C functor — `CANDIDATE_C_FUNCTOR.md` / `TIER12_SHARP_NEXT.md`  
- **1.3** genus estimate — `GENUS_P_QW.md` (superseded for exact \(g\) by blowup doc)

### New algebraic ideas A–F — **`NEW_ALGEBRAIC_IDEAS.md`**

| Idea | Result |
|------|--------|
| **A** Mestre on HQCC seeds | **HIT:** \(P''R-2P'R'\equiv0\pmod P\); `shift_y_tR` 1-param \(A_5\) families |
| **F** Embed \(\to T\) | Seeds + Mestre specs embed; non-BJ hosts |
| **B** \(x^5+75x^3+Ax^2+3A\) | **HIT:** disc \(324A^2(A^2+84375)^2\); \(d=-75\) in \(T\) |
| **C–E** | No extra identical-square-by-search avatar |

### Follow-through — **`NEXT_MESTRE_B_AVATAR.md`**

| Deliverable | Doc | Result |
|-------------|-----|--------|
| Closed-form flagship \(P_t\) | `MESTRE_FLAGSHIP_PT.md` | Explicit; disc□; 10/10 sample \(A_5\) incl. model \(t=61,80\) |
| B-embed lattice \(bc=72A\) | `B_EMBED_LATTICE.md` | 104 unique \(A\), all disc□; 50 \(A_5\) checked |
| Evenness avatars | `EVENNESS_AVATAR.md` | PE + B identities proved; B beyond BJ-embed |

### Ternary lattice \(L_0\) — structure (not necessity)

**Roles:** specialisation source · Mestre parameter · template coordinates.  
**Directions:** `TERNARY_LATTICE_DIRECTIONS.md`.

**Direction 1 done** — `L0_SECONDARY_INVARIANTS.md`:
- PE: 148 fibres; checked Gal A5=36, **D5=4** (refinement signal)
- B: 316 \(A\)-points; A5=45 among checks; disc□ by identity
- Mestre \(t\in L_0\): 18 A5 checked
- Strict PE↔B coordinate overlap: 14 (no map \(\Phi\) yet)

**Directions 2–3** — `L0_MONOID_SATURATION.md`, `L0_PE_B_UNIFY.md`:
- Monoid \(M_0=\langle 3,61,80,243,539\rangle_\times\) size 32; additive \(M_2\) size ~4504
- B irr ~0.99 on monoid sample; PE raw coeffs rarely in monoid
- PE↔B: no rate-1 \(\Phi\); best elementary \(A=k_n k_d\), \(A=\pm\alpha\)

**Direction 4** — `L0_MESTRE_ORBIT.md`: Mestre orbit on 24 PE/B seeds; all have \(R\)-space and disc□ families; 480/480 lattice-\(t\) pairs even; graph 528 nodes; second hops disc□.

**Lattice programme (Dirs 1–4):** invariants · monoid · PE↔B cut · Mestre orbit **done**. Necessity **paused**. Optional Dir 5 only with new axiom-named avatar.

### Review package (generative extensions)

**Order:** (1) flagship Mestre \(P_t\) → (2) B-avatar. Index: **`REVIEW_PACKAGE.md`**.

| Priority | Object | Doc | Auto-check |
|:--------:|--------|-----|------------|
| 1 | Flagship Mestre lift | `MESTRE_FLAGSHIP_PT.md` | `python review_flagship_b.py` |
| 2 | Non-BJ B-embed + avatars | `B_EMBED_LATTICE.md`, `EVENNESS_AVATAR.md` | same script |

Re-verified: Mestre rem=0; disc□ in \(\mathbb{Q}[t]\); \(t=0\) recovers seed; Gal \(A_5\) on lattice \(t\); B disc identity + \(T\)-match.  
Integrity: `MATH_INTEGRITY_REVIEW.md` (30/30). **Contamination boundary:** `CONTAMINATION_BOUNDARY.md` — lattice integers only as arithmetic data; no G₄/539.9 s, GW/Belle II, or 539-step dynamics inside disc proofs.

| Layer | Status |
|-------|--------|
| **Citable centre (Resonant Algebra)** | **Pure-even multi-\(k\) arithmetic theory** |
| **Arithmetic foundations** | **Theorem-grade** |
| **Necessity theorem** (alternating monodromy **forced** by HQCC axioms) | **Open** — `NECESSITY_THEOREM.md` |
| **Geometry (\(\varphi\), monodromy \(A_5\))** | **Done** (construction, not necessity) |
| **Surgery on same rigid \(\varphi/\mathbb{Q}\)** | **Ruled out** (\(\mathrm{disc}=5\cdot\square\)) |
| **BJ pencils / U4** | **Done** |
| **Fused pure \(A_5\) + Hilbert seed lattice** | **Open — structural** (not required for A1 success) |

**Necessity vs generative.** Catalogues and pure-even slices show the lattice *can* feed the classical even locus. They do **not** prove that every (or any canonical) HQCC-built object *must* have alternating monodromy. That is Criteria 1–3 — all open.

**HQCC matrix templates** — **`HQCC_MATRIX_TEMPLATES.md`**: base \(M\) is \(S_5\) (odd disc); \(T(a,\ldots,f)\) produces \(A_5\) only after disc gate; BJ-embed recovers pure-even on a thin subclass — **not** Criterion-2 necessity.

**Resonant algebraic closure** — **`RESONANT_ALGEBRAIC_CLOSURE.md`**: Candidate **A** (lattice + pure-even envelope) = operational closure of finished multi-\(k\); **B** cyclotomic \(\mathcal{R}\); **C** mod-2↔mod-3 design mirror (functor missing); **D** homotopy 539 (not checkable). A does not imply necessity.

**Resolution path:** **`RESOLUTION_PATH.md`**  
**Empirical grounding (A+B):** **`RESONANT_NUMBER_THEORY.md`** — **CONFIRMED**

| Stage | Focus | Empirical |
|-------|--------|-----------|
| **A near-term** | Multi-\(k\) core; structural criteria; generative reach | DIG+GROW **PASS** |
| **B medium-term** | Checkable predictions + data tables | BUILD **PASS** |
| **D medium** | Density / asymptotics (D1–D3) | **COMPLETE** (`STAGE_D_DENSITY.md`) |

**Established (centre)**

1. BJ disc identity; operational \(A_5\); homogenisation lemma.  
2. **Pure-even \(k\)-slices + 2-param envelope + cross-\(k\) paths** (disc identically square).  
3. Catalogue: multi-seed slices (flagship \(k=-8/5\), LSW \(k=-4\), classical \(k=4/5\), …).  
4. Rigid \(\varphi/\mathbb{Q}\) even fibres ruled out; base change \(\mathbb{Q}(\sqrt5)\) side-only.  

**Still open (structural / geometric)**

> Nielsen-labelled geometric multi-\(k\); full Criterion-1 fusion object.

**Do not confuse** arithmetic multi-\(k\) (done) with geometric multi-\(k\) (open).

**Negative control (rigid odd fibre \(t=3\))** — **LOCKED** (`RIGID_FIBRE_T3.md`):
monic(\(\varphi-3\)) irr, disc \(=5\cdot\square\) not □, Gal \(=S_5\). Contrast with pure-even resonant seeds (disc □, Gal \(A_5\)) independently verified on both sides.

**Pure-even specialisations** — **PASS** (`PURE_EVEN_SPECIALISATIONS.md`): 1728 Z-coeff slice pts, irr=1728, even_fail=0, A5 among checks=152; cross-\(k\) paths multi-catalogue; homogenised seeds all sample-even.

**Explicit \(3A^4\) equation** — **`EXPLICIT_3A4_EQUATION.md`**: \(H^{\mathrm{rd}}\cong\mathbb{P}^1\); eliminant \(P(q,w)=0\); deg-5 resolvent \(N-tD\); exact fibre \(s=-1\).

**\(k(s)\) plot** — **`K_OF_S.md`**, `build/k_of_s.png`: 577 BJ fibres, 58 \(s\)-values; median \(k(s)\) does **not** hit ≥2 multi-seed catalogue ratios (strict geometric multi-\(k\) still open).

**Non-classical track** — scaffold **LOCKED** (`NONCLASSICAL_RESONANT_FIELD.md`): polys over \(\mathcal{R}= \mathbb{Q}(2\cos 2\pi/N)\) (proxies \(n=5,7,11,15\); full \(N=539\) deg 210); split-prime Frob labels; field-agnostic pure-even + cosine \(k\)-trials. Enrichment only — does not replace Z/Q centre.

---

## 2. Arithmetic (theorem-grade)

| Result | Notes |
|--------|--------|
| \(\operatorname{disc}(x^5+ax+b)=256a^5+3125b^4\) | Symbolic identity |
| BJ even \(\iff\) that integer is square; +\((3,1,1)\Rightarrow A_5\) | Operational |
| Homogenisation lemma | \(f_t=x^5+\alpha t^4 x+\beta t^5\), disc \(=t^{20}\operatorname{disc}(s)\) |
| **HQCC seeds** | e.g. \(x^5-55x+88\) (\(88=61+3^3\)), \((95,\pm76)\), \((95,\pm532)\), … |
| Empirical | Classical ray 33/33 \(A_5\); flagship ray 9/9 model \(t\) \(A_5\) |
| **Enlarged BJ A5 catalogue** | **60** unique; **10** multi-seed pure-even \(k=\beta/\alpha\) slices (`ENLARGED_SEED_CATALOGUE.md`) |
| Flagship multi-seed | \(k=-8/5\): flagship + \((145,-232)\), \((320,-512)\), \((1145,-1832)\) |
| LSW multi-seed | \(k=-4\): **11** seeds (was 2) |

Catalogues: **36 \(A_5\)** (mixed shapes) + **60 BJ \(A_5\)** seeds; **4 \(A_6\)**.

---

## 3. Geometry

| Step | Result |
|------|--------|
| Rigid triples | Preferred **(3A,3A,5A)**; fallback **(2A,3A,5A)** |
| Preferred Belyi cover | \(\varphi(y)=6y^5-15y^4+10y^3\in\mathbb{Q}[y]\) |
| Branch locus | \(\{0,1,\infty\}\), types \((3,1,1),(3,1,1),(5)\) |
| Geometric monodromy | **\(A_5\)** |
| 9 Maths labels | T₃ / T₃ / \(G_4=539.9\) |
| Resonant base chart | \((3,61,539)\to(0,1,\infty)\) |

---

## 4. Fusion probes

| Probe | Outcome |
|-------|---------|
| Affine fibres of \(\varphi\) vs seeds | **0** equation matches |
| Mild \(\varphi\) twists | **Blocked** (cannot stay BJ) |
| Linear BJ pencils | Seeds at \(t=0,1\); disc **not** □ in \(\mathbb{Q}(t)\) |
| 38 BJ families scan | Pure even = **homogenisation rays only** |
| U4 residue \(2\mapsto\sigma_\infty\) | Canonical on 3-point cover |
| U5 4-point cover | \(A_5\)-generating, not abs. rigid |

---

## 5. Gap A — even surface search (**complete**, ~99 s)

Script: `gap_a_even_surface.py` → **`GAP_A_EVEN_SURFACE.md`**  
Surface: \(256\alpha^5 + 3125\beta^4 = \gamma^2\)

| Method | Tested | Pure-even multi-seed hits |
|--------|-------:|--------------------------:|
| Homogenisation rays (1 seed) | 12 | 0 (1 seed each by design) |
| Linear pencils in \((\alpha,\beta)\) | 66 | **0** |
| Lines through \((\alpha,\beta,\gamma)\) lifts | 264 | **0** |
| Monomial bridges | 45 | **0** |
| Quadratic Bezier midpoints | 1350 | **0** |
| Rational quadratic (sample+algebra) | 44145 | **0** |
| Plane sections (3 seed lifts) | 20 | **0** (all \(F\) irreducible deg 5) |

**Conclusion:** No rational curve of the scanned types carries **two or more** distinct HQCC seeds while keeping disc a square in \(\mathbb{Q}(u)\). Pure-even curves remain **single-seed homogenisation rays**.

---

## 6. Remaining viable routes only

1. **Non-rigid \(A_5\) families** (positive-dimensional Hurwitz spaces).  
2. **Arithmetic over a number field \(K\)** in which the permanent factor 5 becomes a square.  
3. **Entirely different geometric constructions** not arising from the rigid triples already examined.

**Do not:** further surgery on the same rigid \(\varphi/\mathbb{Q}\).

---

## 7. Key files

| File | Role |
|------|------|
| `THEOREMS.md` | Master ledger |
| `FUSION_GAP.md` | Principal open fusion problem |
| `HQCC_SEED.md` | Seed definition |
| `GEOMETRIC_STEP2.md` | Cover \(\varphi\) |
| `GAP_A_BJ_FAMILIES.md` | 38-family disc/specialisation scan |
| `FUSION_NEXT.md` | Pencils + U4/U5 |
| `QUEUED_ATTACKS.md` | Spec match / base / scaffold functor |
| `gap_a_even_surface.py` | Surface curve search (re-run to complete) |

---

## 8. One-paragraph summary

The programme has **theorem-grade arithmetic** (HQCC seeds + homogenisation) and an **explicit \(\mathbb{Q}\)-Belyi cover with monodromy \(A_5\)**. Fusion progress includes equation-level BJ pencils and U4’s canonical braid assignment. Systematic scans show that **pure-even BJ families are essentially homogenisation rays through a single seed**; linear and low-degree multi-seed bridges leave the even surface. A **single pure geometric \(A_5\) family that Hilbert-recovers the full HQCC seed lattice** remains the principal open problem.

_End of report._
