# Realistic resolution path — Resonant Algebra / Resonant Galois

_Locked programme stance: mathematical integrity + scientific method._  
_Stages **A–B** executed and empirically grounded; stages **C–Z** complete the programme ledger._

---

## Bottom line

| Layer | Status |
|-------|--------|
| **Published / citable centre** | **Pure-even multi-\(k\) arithmetic theory** |
| **Necessity theorem** (forced \(A_n\) from HQCC axioms) | **Open** — `NECESSITY_THEOREM.md`; **not** Stage A centre |
| **Research roadmap** | **`RESEARCH_ROADMAP.md`** — Tier 1–3; do not reopen finished layers |
| **Structural fusion** (geometry + lattice as one object) | **Open** — Stages F, M, N, R; not required for A success |
| **Generative reach** beyond \(A_5\) | **Stage A near-term** (scaffold) → **G, S** |
| **Independent arithmetic predictions** | **Stage B** (empirical tables) → **D, O, V** |
| **Programme completion (C–Z)** | **Roadmap locked** + inventory run (`RESONANT_NUMBER_THEORY_CZ.md`) |

The programme does **not** wait on a Nielsen-labelled geometric multi-\(k\) resolvent to have a finished positive result. That geometric problem remains open and high-value; the **citable core** is the arithmetic theory already proved and catalogued.

---

## Stage ladder (A–Z at a glance)

| Stage | Horizon | Theme | Centre? |
|-------|---------|-------|:-------:|
| **A** | near | Secure mathematical core (multi-\(k\), criteria attacks, generative scaffold) | **Yes** |
| **B** | medium | Independent checkable arithmetic predictions | **Yes** |
| **C** | near–med | Structural **criteria** campaign (Crit 1–3 / necessity — `NECESSITY_THEOREM.md`) | open work |
| **D** | medium | **Density** / asymptotic arithmetic — **COMPLETE** | B upgrade |
| **E** | near | **External** reproducibility package | A/B support |
| **F** | ongoing | **Fusion** research track (non-blocking) | open |
| **G** | medium | **Generative** \(A_n\) theorem-grade (\(n\ge 6\)) | A3 upgrade |
| **H** | long | **Hilbert** modular / icosahedral \(\mathbb{Q}(\sqrt5)\) | enrichment |
| **I** | near | **Integration** pipeline (`build_all`, mirrors) | infra |
| **J** | medium | **Journal**-ready citable artefacts | publication |
| **K** | medium | **\(k\)-class** taxonomy & characterisation | arithmetic |
| **L** | near–med | **Lattice** / HQCC dictionary formalisation | arithmetic |
| **M** | long | **Multi-\(k\) geometric** (Nielsen-labelled) | open |
| **N** | long | **Nielsen / Hurwitz** \(r=4,5+\) systematic | open |
| **O** | medium | **Obstruction** theory completeness | arithmetic+geo |
| **P** | near | **Portable** pure-even library package | A1 product |
| **Q** | med–long | **Quadratic** & number-field base change | side route |
| **R** | long | **Resolvent** closed form \(f_s\in\mathbb{Q}(s)[x]\) | open |
| **S** | long | **Simple groups** beyond alternating | generative |
| **T** | medium | **Theorem** stack hardening / machine-checked ids | core |
| **U** | med–long | **U4 / U5** naturality & functoriality | structural |
| **V** | near | **Verification** regression suite (all PASS gates) | infra |
| **W** | ongoing | **Writing** / exposition / monographs | publication |
| **X** | medium | e**X**otic poly shapes & other disc identities; **non-classical \(\mathcal{R}\)** (`NONCLASSICAL_RESONANT_FIELD.md`) | generative |
| **Y** | near | **Yield** KPIs / programme metrics dashboard | meta |
| **Z** | terminal | **Zenith** synthesis + parked-research registry | ledger |

---

## Stage A — Secure the mathematical core (near-term)

### A1. Keep pure-even multi-\(k\) theory as the published centre

**Object.** Bring–Jerrard quintics \(x^5+\alpha x+\beta\in\mathbb{Z}[x]\) with
\[
\operatorname{disc}=256\alpha^5+3125\beta^4
\]
a perfect square (even monodromy among transitive subgroups of \(S_5\)), and
\[
k=\beta/\alpha\in\mathbb{Q}
\]
the **ratio class**.

**Theorems (citable).**

1. **BJ disc identity.** \(\operatorname{disc}(x^5+\alpha x+\beta)=256\alpha^5+3125\beta^4\).
2. **Operational \(A_5\).** Irr + disc square + Frobenius type \((3,1,1)\) \(\Rightarrow\mathrm{Gal}=A_5\).
3. **Homogenisation lemma.** If disc(seed) is square, then
   \(f_t=x^5+\alpha t^4 x+\beta t^5\) has disc \(=t^{20}\operatorname{disc}(\mathrm{seed})\) square for all \(t\neq0\).
4. **Pure-even \(k\)-slice.** For any fixed \(k\in\mathbb{Q}\setminus\{0\}\),
   \[
   \alpha(m)=256 m^2-\frac{3125\,k^4}{256},\qquad
   \beta(m)=k\cdot\alpha(m)
   \]
   satisfies
   \[
   \operatorname{disc}=(256\,\alpha(m)^2\, m)^2
   \]
   identically in \(\mathbb{Q}(m)\). Hence every specialisation with \(m\neq0\) and \(\alpha(m)\neq0\) is even (when irreducible: Gal \(\le A_5\); with \((3,1,1)\): \(A_5\)).
5. **2-parameter envelope.** The same formulae with free \((m,s)\) and \(k=s\) give a pure-even family over \(\mathbb{Q}(m,s)\) that recovers every fixed-\(k\) slice by freezing \(s=k\).
6. **Cross-\(k\) paths.** Any two envelope points (including different \(k\)) lie on a pure-even 1-param path via a rational path \((m(u),k(u))\) in parameter space (e.g. linear). Specialisations recover multiple catalogue ratio classes (flagship \(k=-8/5\), classical \(k=4/5\), LSW \(k=-4\), …).

**Catalogue (empirical, reproducible).**  
~60 BJ \(A_5\) seeds; 10 multi-seed pure-even \(k\)-slices; explicit paths flagship↔classical, flagship↔LSW, etc.  
Docs: `ENLARGED_SEED_CATALOGUE.md`, `NONRIGID_MULTI_SEED.md`, `NONRIGID_HURWITZ_SEARCH.md`, `REALISE_3A4_SPECIALISE.md`.  
Code: `lib/lemmas.py`, `enlarge_seed_catalogue.py`, `nonrigid_multi_seed.py`, `realise_3a4_specialise.py`.

**What this is not.** It is **not** a claim that these families are Nielsen-labelled Hurwitz components (e.g. \(\mathrm{Ni}(A_5,C_3^4)\)). Geometric multi-\(k\) remains open.

**Publication stance.** Present theorems 1–6 + catalogue as **Resonant Algebra / arithmetic core**. Mark geometric fusion as open problem, not as a missing proof of the arithmetic claims.

---

### A2. Continue attacks on structural criteria

Keep working, without forcing premature closure:

| Criterion / gap | Near-term attack | Status |
|-----------------|------------------|--------|
| **Canonical geometric object** | \(3A^4\) resolvent \(f_s\in\mathbb{Q}(s)[x]\); other g=0 shortlist types | Open (`TRIPLE_ROOT_ELIMINATE.md`, `EXPLICIT_3A4_RESOLVENT.md`) |
| **Sign character / evenness** | Disc square gate; base change \(K=\mathbb{Q}(\sqrt5)\); permanent factor 5 on \(\varphi\) | Partial: arithmetic gate closed; geometric evenness over \(\mathbb{Q}\) blocked for rigid \(\varphi\) |
| **U4 / residue→braid naturality** | 4-point covers; functoriality | Partial (`FUSION_NEXT.md`) |
| **Hilbert recovery of HQCC lattice** | Specialise geometric families onto catalogue \(k\)-slices | Arithmetic paths yes; Nielsen-labelled no |

**Scientific method.** Separate:

- *proved arithmetic statements* (A1),
- *structural open problems* (this table),
- *high-effort speculative avenues* (Hilbert modular / icosahedral for \(\sqrt5\), higher-rank Hurwitz, envelope Nielsen ID — `HILBERT_MODULAR_A5.md`, `AVENUE_RANK_EXECUTE.md`).

Do not blur arithmetic multi-\(k\) with geometric multi-\(k\).

---

### A3. Enlarge generative reach (not a one-off for \(A_5\))

Make the method visibly portable:

| Direction | Near-term concrete work |
|-----------|-------------------------|
| **Higher alternating groups \(A_n\)** | Depressed / sparse forms; disc-square searches; operational criteria (e.g. \(A_n\) via even + transitive + suitable cycle types); catalogues for \(n=6,7\) where feasible |
| **Other simple groups** | Thin subclasses (e.g. \(A_6\) already in catalogue); document disc / Frobenius gates analogously |
| **Homogenisation / envelope pattern** | Abstract: “seed with square disc + weighted 1-param family with disc \(= t^{N}\) disc(seed)” as a reusable lemma pattern |
| **Ratio classes** | For other poly shapes, identify analogous “thin pure-even loci” (not only BJ \(k=\beta/\alpha\)) |

**Success metric for A3.** At least one theorem-grade statement of the form:  
*for a stated polynomial class, disc-square is an identity on an explicit parametric family, and specialisations yield \(\mathrm{Gal}=A_n\) (or other \(G\)) under an operational criterion* — for some \(n\neq5\) or some \(G\neq A_5\).

Scaffold: `generative_reach.py` / `GENERATIVE_REACH.md` (A6 catalogue already exists; extend systematically).

---

## Stage B — Independent arithmetic predictions (medium-term)

Produce statements **other mathematicians can check** without accepting the full Resonant narrative:

### B1. New families

- Further multi-seed pure-even slices (more \(k\), denser lattice).
- Explicit cross-\(k\) paths with proved disc identity and listed Hilbert \(A_5\) specialisations.
- Optional: pure-even families outside BJ (icosahedral, \(pqr\) forms) with closed disc formulae.

### B2. Density / asymptotic statements

Candidates (to be proved or carefully conjectured with numerical support):

- Density of integer \(m\) for which \(x^5+\alpha(m)x+\beta(m)\) is irreducible (fixed \(k\)-slice).
- Chebotarev: distribution of Frobenius types along a pure-even family.
- Growth of \(|\operatorname{disc}|\) along \(k\)-slices / envelope.

### B3. Constraints on Galois groups and discriminants

- **Necessary** conditions: square disc for even monodromy on BJ; permanent factor 5 for monic(\(\varphi-t\)) over \(\mathbb{Q}\).
- **Catalogue constraints:** which rational \(k\) appear among HQCC-lattice \(A_5\) seeds (empirical → conjectural characterisation).
- **Obstruction theorems** already usable by others: no even irr rational fibre of preferred \(\varphi/\mathbb{Q}\).

### B4. Reproducibility

- Public lemmas in `lib/lemmas.py` with symbolic checks.
- JSON catalogues + scripts that regenerate `ENLARGED_SEED_CATALOGUE`, multi-seed paths, disc identities.
- Clear “theorem vs empirical vs open” labels in `THEOREMS.md` and this file.

---

## What is deliberately not Stage A/B centre

| Item | Role |
|------|------|
| Nielsen-labelled geometric multi-\(k\) | Open research; high value; not required for citable arithmetic core |
| Hilbert modular forms for \(\mathbb{Q}(\sqrt5)\) | Enrichment of base-change avenue; high effort |
| Full Criterion-1 fusion object | Principal open problem; structural attacks continue under A2 / F / M / R |
| Claiming \(\varphi\)-surgery or Q-even rigid fibres | **Ruled out** — do not reopen without new geometry |

---

## Immediate action list (Stage A) — executed

1. **Freeze A1** as the citable centre: theorems 1–6 + catalogue pointers in `THEOREMS.md` / this file.  
2. **A2:** one active structural attack at a time (prefer 3A⁴ closed form *or* sign-character theorem sketch — not both at full intensity).  
3. **A3:** ship one portable lemma + one non-\(A_5\) catalogue slice (e.g. \(A_6\) or deg-6 pure-even family probe).  
4. **Stage B prep:** list 3 checkable predictions (density / new family / obstruction) with scripts that emit verifiable data.

### Empirical grounding run (executed)

Script: `resonant_number_theory_ab.py` → **`RESONANT_NUMBER_THEORY.md`**, `RNT_STAGE_B_DATA.json`.

| block | result |
|-------|--------|
| DIG (A1) | **PASS** — identities, catalogue 21/21 A5, paths multi-\(k\), even_fail=0 |
| GROW (A3) | **PASS** — \(A_6\) examples; envelope harvest 50 A5 |
| BUILD (B) | **PASS** — B1–B4 tables; \(\varphi\) obstruction identity |
| **RNT empirical grounding** | **CONFIRMED** |

---

# Stages C–Z — Full programme ledger

_These stages complete the resolution path beyond A/B. Status tags:_

| Tag | Meaning |
|-----|---------|
| **LOCKED** | Deliverable frozen; treat as programme constant |
| **ACTIVE** | In progress / next work queue |
| **SCAFFOLD** | Artefacts exist; theorem-grade not yet claimed |
| **OPEN** | Research open; high value; non-blocking for centre |
| **PARKED** | Deliberately deferred or ruled out under current geometry |
| **INFRA** | Engineering / reproducibility support |

Inventory runner: `python resonant_number_theory_cz.py` → **`RESONANT_NUMBER_THEORY_CZ.md`**.

---

## Stage C — Structural criteria campaign

**Horizon:** near–medium · **Tag:** ACTIVE  
**Extends:** A2 · **Feeds:** F, M, R, U

### Goal
Turn the Criterion-1/2/3 gap table into a **single campaign plan** with one primary attack at a time, explicit stop rules, and no confusion with arithmetic multi-\(k\).

### Deliverables
1. Ranked structural attack queue (current: prefer \(3A^4\) closed form **or** sign-character theorem sketch).
2. Written separation: arithmetic evenness (done) vs geometric monodromy object (open).
3. Status table updated after each attack cycle (`FUSION_GAP.md`, `AVENUE_RANK_EXECUTE.md`).

### Success criterion
At least one structural criterion advances from “partial” to “proved” **or** is rigorously **ruled out** with a published obstruction (as \(\varphi/\mathbb{Q}\) already is).

### Existing artefacts
`criterion1_hqcc.py`, `criterion2_axioms.py`, `criterion3_sign.py`, `FUSION_GAP.md`, `THEOREM_ATTACK.md`, `AVENUE_RANK_EXECUTE.md`.

---

## Stage D — Density and asymptotic arithmetic

**Horizon:** medium · **Tag:** **COMPLETE** (executed)  
**Extends:** B2 · **Feeds:** J, Y

### Goal
Upgrade empirical rates (irr along \(k\)-slices; Gal histograms on paths) into **stated theorems or precise conjectures** with numerical support and clear hypotheses.

### Deliverables
1. Irreducibility density statements for fixed multi-seed \(k\)-slices (conjecture + large-\(m\) tables).
2. Chebotarev / Frobenius-type histograms along pure-even families.
3. Disc-height growth bounds along envelope rays.

### Success criterion
≥1 density statement published as **theorem** or as **conjecture with machine-checkable evidence table** outsiders can regenerate.

### Executed

Script: `stage_d_density.py` → **`STAGE_D_DENSITY.md`**, `STAGE_D_DATA.json`.

| block | status | result |
|-------|--------|--------|
| **D1** | conjecture + evidence | irr rate **1.0** on all multi-seed slices sampled; even_fail=0; A5 among Gal checks ≈28–30/30 |
| **D2** | Chebotarev proxy | 100 irr fibres; 4500 prime factorisations; types compatible with \(A_5\) classes |
| **D3** | **proved** | \(\log\|\mathrm{disc}\|\sim 10\log\|m\|+48\log 2\); \(\sqrt{\mathrm{disc}}=\|256\alpha^2 m\|\) match 100% |
| **Stage D** | **COMPLETE** | |

### Existing artefacts
`STAGE_D_DENSITY.md`, `build/STAGE_D_DATA.json`, prior B1/B2 in `RESONANT_NUMBER_THEORY.md`.

---

## Stage E — External reproducibility package

**Horizon:** near · **Tag:** ACTIVE / INFRA  
**Extends:** B4 · **Feeds:** J, V, W

### Goal
A third party can re-run the arithmetic core and Stage B tables with only standard Python + sympy, without accepting Resonant narrative.

### Deliverables
1. Minimal dependency list + pinned run scripts.
2. “Outsider checklist”: identities → catalogue → paths → \(\varphi\) obstruction.
3. JSON schemas for catalogues and Stage B data.

### Success criterion
Clean re-run of `resonant_number_theory_ab.py` + key lemma checks from a cold clone; docs point to exact commands.

### Existing artefacts
`lib/lemmas.py`, `lib/common.py`, `resonant_number_theory_ab.py`, `build/*.json`, `README.md`.

---

## Stage F — Fusion research track (non-blocking)

**Horizon:** ongoing · **Tag:** OPEN  
**Extends:** A2 geometric half · **Feeds:** M, N, R, Z

### Goal
Continue the **principal open fusion problem** without blocking A/B success: a single geometric object that is pure \(A_5\) and recovers the HQCC lattice as Hilbert specialisations.

### Deliverables
1. Living status in `FUSION_GAP.md` (locked bottom line).
2. Viable-route only list: non-rigid Hurwitz; base-change side; other geometry.
3. Explicit **do not claim** list (rigid \(\varphi\) surgery ruled out).

### Success criterion (optional bonus)
Fusion object found — **bonus**, not definition of programme success under this path.

### Existing artefacts
`FUSION_GAP.md`, `FUSION_NEXT.md`, `GEOMETRIC_*.md`, `NONRIGID_*.md`, `REALISE_3A4_SPECIALISE.md`.

---

## Stage G — Generative \(A_n\) theory (\(n\ge 6\))

**Horizon:** medium · **Tag:** SCAFFOLD  
**Extends:** A3 · **Feeds:** S, X, J

### Goal
Portable pure-even / homogenisation method with **theorem-grade** statements for at least one \(n\neq 5\).

### Deliverables
1. Operational criterion template for \(A_n\) (even + transitive + cycle types).
2. Explicit pure-even or disc-square parametric families in deg 6+ where feasible.
3. Catalogue slices with Gal \(A_6\) (and \(A_7\) probes if tractable).

### Success criterion
One theorem-grade identity + operational Gal gate for some \(n\neq 5\), with ≥1 verified family.

### Existing artefacts
`generative_reach.py`, `GENERATIVE_REACH.md`, RNT GROW \(A_6\) examples, `CATALOGUE.md` A6 freeze.

---

## Stage H — Hilbert modular / icosahedral \(\mathbb{Q}(\sqrt5)\)

**Horizon:** long · **Tag:** OPEN / enrichment  
**Extends:** Avenue 5 · **Feeds:** Q, F

### Goal
Use the classical Hirzebruch–Klein / Hilbert modular surface link for \(\mathbb{Q}(\sqrt5)\) as **enrichment** of the permanent-factor-5 obstruction and 3A⁴ base-change, not as a short-cut past arithmetic multi-\(k\).

### Deliverables
1. Narrative + references map (`HILBERT_MODULAR_A5.md`).
2. Explicit dictionary: disc factor 5 ↔ units / CM / icosahedral invariants (where rigorous).
3. Optional: special values recovering catalogue \(k\) over \(K=\mathbb{Q}(\sqrt5)\).

### Success criterion
Published enrichment note usable by specialists; multi-\(k\) over \(\mathbb{Q}\) not claimed from modular forms alone.

### Existing artefacts
`HILBERT_MODULAR_A5.md`, `K_SQRT5_EVEN.md`, `hilbert_modular_a5.py`, `k_sqrt5_even.py`.

---

## Stage I — Integration pipeline

**Horizon:** near · **Tag:** INFRA  
**Extends:** package engineering · **Feeds:** E, V, Y

### Goal
Single entry that rebuilds core docs and records pass/fail (`build_all.py` + C–Z inventory).

### Deliverables
1. `build_all.py` covers theorem / geometric / fusion modules.
2. C–Z inventory script mirrors outputs to `build/`.
3. Stable paths: `ROOT`, `OUT`, optional `RESULTS` mirror.

### Success criterion
One-command rebuild of critical md/json without silent identity failures.

### Existing artefacts
`build_all.py`, `lib/common.py`, `resonant_number_theory_ab.py`, `resonant_number_theory_cz.py`.

---

## Stage J — Journal-ready citable artefacts

**Horizon:** medium · **Tag:** ACTIVE  
**Extends:** A1 + B · **Feeds:** W, Z

### Goal
Artefacts suitable for citation: theorem statements, proofs/sketches, catalogues, obstruction theorems — cleanly separated from open geometric fusion.

### Deliverables
1. Paper-shaped outline: **Resonant Algebra — pure-even multi-\(k\) families**.
2. Lemma file + theorem ledger export.
3. Data appendix (JSON + hash or version stamp).

### Success criterion
A mathematician outside the programme can cite theorems 1–6 + \(\varphi\) obstruction without reading fusion notes.

### Existing artefacts
`THEOREMS.md`, `RESOLUTION_PATH.md`, `RESONANT_NUMBER_THEORY.md`, `lib/lemmas.py`.

---

## Stage K — \(k\)-class taxonomy

**Horizon:** medium · **Tag:** SCAFFOLD  
**Extends:** catalogue · **Feeds:** D, L, J

### Goal
Characterise which rational \(k=\beta/\alpha\) appear among HQCC-lattice \(A_5\) seeds; organise multi-seed slices.

### Deliverables
1. Complete multi-seed slice table (current: 10 slices, 16 ratios).
2. Conjectural constraints on \(k\) (denominator, sign-flip symmetry, lattice origin).
3. Cross-\(k\) path inventory linking named rays (flagship, classical, LSW, …).

### Success criterion
Taxonomy document with conjectures labelled as such; all multi-seed \(k\) regenerable from scripts.

### Existing artefacts
`ENLARGED_SEED_CATALOGUE.md`, `enlarge_seed_catalogue.py`, path docs in RNT DIG.

---

## Stage L — Lattice / HQCC dictionary

**Horizon:** near–medium · **Tag:** SCAFFOLD  
**Extends:** HQCC seed def · **Feeds:** J, F, M

### Goal
Formal dictionary: model integers \(\{3,9,61,80,243,539,\ldots\}\) ↔ seeds ↔ homogenised families ↔ (optional) branch labels.

### Deliverables
1. Frozen vocabulary: **HQCC seed** (`HQCC_SEED.md`).
2. Explicit flagship arithmetic: \(88=61+3^3\), etc.
3. Separation of lattice origin from geometric cover claims.

### Success criterion
Seed definition usable without monodromy; geometric labels optional enrichment only.

### Existing artefacts
`HQCC_SEED.md`, `HQCC_NATIVE.md`, `HQCC_STRICT.md`, model core in `lib/common.py`.

---

## Stage M — Geometric multi-\(k\) (Nielsen-labelled)

**Horizon:** long · **Tag:** OPEN  
**Extends:** F · **Feeds:** R, N, Z

### Goal
A Nielsen-labelled geometric family whose specialisations realise **multiple** catalogue \(k\) (flagship + classical + LSW, …).

### Deliverables
1. Candidate Hurwitz type + reduced space model.
2. Proof or strong numerical evidence of multi-\(k\) Hilbert specialisations.
3. Explicit comparison to arithmetic envelope (control, not substitute).

### Success criterion
Documented geometric multi-\(k\) with Nielsen type; **or** obstruction ruling out a listed type.

### Existing artefacts
`A5_HURWITZ_R4.md`, `REALISE_3A4_SPECIALISE.md`, `AVENUE_RANK_EXECUTE.md` (geometric multi-\(k\) still open).

---

## Stage N — Nielsen / Hurwitz systematic (\(r=4,5+\))

**Horizon:** long · **Tag:** OPEN  
**Extends:** M · **Feeds:** R, F

### Goal
Systematic treatment of filter-pass Nielsen types for \(A_5\): genus, braid orbits, equation models where feasible.

### Deliverables
1. \(r=4\) type census (done scaffold: 19 filter-pass; g=0 shortlist).
2. Explicit models for shortlist beyond \(3A^4\) (\(2A3A^3\), \(2A^2 3A^2\), …).
3. \(r\ge 5\) effort estimate + only pursue if dim / class-space justifies.

### Success criterion
≥1 additional g=0 type with explicit cover equations over a number field (ideally \(\mathbb{Q}\)).

### Existing artefacts
`a5_hurwitz_r4.py`, `A5_HURWITZ_R4.md`, `GENUS_3A4_LOCK.md`, Avenue 2/6 notes.

---

## Stage O — Obstruction theory completeness

**Horizon:** medium · **Tag:** ACTIVE  
**Extends:** B3 · **Feeds:** C, F, J

### Goal
Catalogue **proved** obstructions (what cannot work) alongside constructive theorems.

### Deliverables
1. \(\varphi/\mathbb{Q}\): disc \(=5\cdot\square\) — no even irr rational fibre (**proved**).
2. Mild surgery / further surgery on same rigid cover — **ruled out**.
3. Base-change: even-over-\(K=\mathbb{Q}(\sqrt5)\) without descent to \(\mathbb{Q}\) (**proved side fact**).

### Success criterion
Obstruction ledger section in `THEOREMS.md` / this file; each entry has status proved/partial/open.

### Existing artefacts
`GEOMETRIC_RIGID_DEFORM.md`, `K_SQRT5_EVEN.md`, RNT B3, `FUSION_GAP.md`.

---

## Stage P — Portable pure-even library package

**Horizon:** near · **Tag:** ACTIVE / INFRA  
**Extends:** A1 · **Feeds:** E, G, X

### Goal
`lib/lemmas.py` (+ thin API) as the reusable pure-even engine: disc, homogenisation, \(k\)-slice, envelope, path disc checks.

### Deliverables
1. Documented public functions with identity self-tests.
2. No Float pollution: exact `Rational` parameters only.
3. Import path usable by all programme scripts.

### Success criterion
All Stage A1 identities re-verified via library calls alone in the C–Z inventory.

### Existing artefacts
`lib/lemmas.py`, `lib/common.py`, imports across catalogue scripts.

---

## Stage Q — Quadratic and number-field base change

**Horizon:** med–long · **Tag:** SCAFFOLD / side  
**Extends:** Avenue 5 · **Feeds:** H, O

### Goal
Theory of evenness after base change: when permanent disc factors become squares; what descends.

### Deliverables
1. \(K=\mathbb{Q}(\sqrt5)\) case complete for preferred \(\varphi\) (even-over-\(K\), no descent).
2. Template for other square-free factors.
3. Explicit non-claims: no automatic HQCC \(\mathbb{Z}\)-lattice recovery.

### Success criterion
Base-change note with proved identities; lattice recovery attempts documented fail-closed.

### Existing artefacts
`K_SQRT5_EVEN.md`, `k_sqrt5_even.py`, Avenue 5 execute block.

---

## Stage R — Resolvent closed form \(f_s\in\mathbb{Q}(s)[x]\)

**Horizon:** long · **Tag:** OPEN  
**Extends:** M, N · **Feeds:** F, Z

### Goal
Trusted closed-form BJ (or equivalent) model for \(\mathrm{Ni}(A_5,C_3^4)\) / \(3A^4\) reduced Hurwitz \(\cong\mathbb{P}^1_s\).

### Deliverables
1. Genus 0 lock + orbit size (done: Bailey–Fried; orbit 18).
2. Explicit cover form / numeric fibres (partial).
3. Closed form over \(\mathbb{Q}(s)\) — **still open** (triple-root elim collapsed / untrusted fits).

### Success criterion
Algebraically verified \(f_s\in\mathbb{Q}(s)[x]\) with correct monodromy type; specialisations match Nielsen data.

### Existing artefacts
`EXPLICIT_3A4_RESOLVENT.md`, `build_3a4_resolvent.py`, `TRIPLE_ROOT_ELIMINATE.md`, `GENUS_3A4_LOCK.md`.

---

## Stage S — Simple groups beyond alternating

**Horizon:** long · **Tag:** OPEN  
**Extends:** G · **Feeds:** X, J

### Goal
Where thin subclasses exist, document disc / Frobenius gates for other simple \(G\le S_n\).

### Deliverables
1. Scope note: which \(G\) are realistic for sparse forms.
2. At least exploratory catalogue hits (if any) with honest status.
3. Do not force \(G\) claims without operational criteria.

### Success criterion
Either a verified non-\(A_n\) simple example in programme scope, or a documented negative survey in bounds.

### Existing artefacts
Catalogue machinery in `lib/common.py`; A5/A6 focus to date.

---

## Stage T — Theorem stack hardening

**Horizon:** medium · **Tag:** ACTIVE  
**Extends:** A1 · **Feeds:** J, E, V

### Goal
Every citable arithmetic claim has: statement, proof sketch or symbolic check, code path, status label.

### Deliverables
1. `THEOREMS.md` ledger complete for theorems 1–9 (multi-\(k\) core).
2. Symbolic random + identity checks in `lib/lemmas.py`.
3. No theorem/empirical/open label collisions.

### Success criterion
Third-party audit can mark each ledger line as theorem / empirical / open without ambiguity.

### Existing artefacts
`THEOREMS.md`, `theorem_attack.py`, `lib/lemmas.py`, RNT DIG identities.

---

## Stage U — U4 / U5 naturality

**Horizon:** med–long · **Tag:** PARTIAL / OPEN  
**Extends:** structural criteria · **Feeds:** F, C

### Goal
Canonical residue→braid assignment (U4) and 4-point cover naturality (U5) as functorial data, not ad hoc labels.

### Deliverables
1. U4 canonical on present 3-point cover (done scaffold).
2. U5: \(A_5\)-generating 4-point covers documented.
3. Naturality statements under base change / Hurwitz braid action.

### Success criterion
Published naturality lemma for U4; U5 either naturalised or scoped as partial.

### Existing artefacts
`FUSION_NEXT.md`, `fusion_next.py`, geometric step docs.

---

## Stage V — Verification regression suite

**Horizon:** near · **Tag:** INFRA / ACTIVE  
**Extends:** E · **Feeds:** I, Y, Z

### Goal
One script (or suite) that re-checks all **PASS** gates: disc identities, catalogue A5, path multi-\(k\), even_fail=0, \(\varphi\) obstruction, A6≥1.

### Deliverables
1. RNT A+B empirical grounding (done).
2. C–Z inventory against artefact map (this stage set).
3. Optional: `build_all` exit code policy for identity failures.

### Success criterion
Regression green ⇔ programme arithmetic centre still sound.

### Existing artefacts
`resonant_number_theory_ab.py`, `resonant_number_theory_cz.py`, `build_all.py`.

---

## Stage W — Writing and exposition

**Horizon:** ongoing · **Tag:** ACTIVE  
**Extends:** J · **Feeds:** external audience

### Goal
Clear prose for outsiders: what is proved, what is open, what is ruled out.

### Deliverables
1. `REPORT.md` bottom line locked.
2. `RESOLUTION_PATH.md` (this file) as method + ladder.
3. Per-probe md files with verdicts.

### Success criterion
Reader of `REPORT.md` + this file understands centre vs open fusion in one sitting.

### Existing artefacts
All top-level `*.md` reports; `IMPLICATIONS.md`.

---

## Stage X — Exotic polynomial shapes

**Horizon:** medium · **Tag:** SCAFFOLD  
**Extends:** G · **Feeds:** S, P

### Goal
Disc identities and pure-even loci outside BJ quintics (icosahedral, \(pqr\), sparse sextics, …).

### Deliverables
1. List of shapes with closed disc formulae.
2. Thin pure-even loci where found.
3. Negative scans documented with bounds.

### Success criterion
≥1 non-BJ shape with identity-level even family **or** exhaustive-in-bounds negative report.

### Existing artefacts
Deg-6 probes in generative reach; Mestre/quad scans noted empty in fusion notes; T5 embed subclass.

---

## Stage Y — Yield KPIs / metrics dashboard

**Horizon:** near · **Tag:** INFRA / meta  
**Extends:** all · **Feeds:** Z

### Goal
Quantitative programme health: seed counts, multi-seed slices, identity pass rates, open problem count, empirical grounding flags.

### Deliverables
1. Scorecard tables in RNT A+B and C–Z docs.
2. Catalogue counts (60 BJ A5, 10 multi-seed slices, …).
3. Stage status histogram (LOCKED / ACTIVE / OPEN / …).

### Success criterion
C–Z inventory emits KPI table regenerated from artefacts + light checks.

### Existing artefacts
Catalogue md/json, RNT scorecards, this ladder.

---

## Stage Z — Zenith synthesis and parked-research registry

**Horizon:** terminal ledger · **Tag:** LOCKED (as process)  
**Extends:** entire programme · **Feeds:** future work only

### Goal
Single place that states: (1) what success means under the realistic path; (2) what is parked; (3) what would reopen a ruled-out route.

### Deliverables
1. **Success definition (realistic path):**
   - Pure-even multi-\(k\) published as citable centre (**A1**).
   - Structural opens listed without confusion (**A2 / C / O**).
   - Generative reach beyond \(A_5\) at least scaffolded (**A3 / G**).
   - Independent predictions + data (**B / D**).
   - Empirical grounding **CONFIRMED** (DIG/GROW/BUILD).
   - Fusion geometric multi-\(k\) remains **optional bonus**.
2. **Parked / ruled-out registry:**
   - \(\varphi\)-surgery for Q-even irr fibres — **ruled out**.
   - Claiming arithmetic multi-\(k\) = Nielsen multi-\(k\) — **forbidden**.
   - Hilbert modular as short-cut past envelope — **parked as enrichment only**.
3. **Reopen rules:** new geometric construction not reducible to rigid \(\varphi/\mathbb{Q}\); or closed-form resolvent with multi-\(k\) specialisations.

### Success criterion
This section + inventory scorecard are current; A+B CONFIRMED; C–Z stages all defined with tags.

### Existing artefacts
This file; `FUSION_GAP.md`; `REPORT.md`; `RESONANT_NUMBER_THEORY.md`; `RESONANT_NUMBER_THEORY_CZ.md`.

---

## Success criteria (full ladder)

| Stage | Done when |
|-------|-----------|
| **A** | Pure-even multi-\(k\) cleanly published as core; structural opens listed; ≥1 generative result beyond \(A_5\) |
| **B** | Independent predictions exist that outsiders can verify/refute with standard NT + released formulae |
| **C** | Structural campaign has ranked queue + stop rules; ≥1 criterion advanced or ruled out |
| **D** | ≥1 density theorem or precise conjecture + regenerable evidence |
| **E** | Cold-clone outsider re-run path documented |
| **F** | Fusion status living; non-blocking; viable routes only |
| **G** | Theorem-grade \(A_n\) (\(n\neq5\)) identity + gate, or documented scaffold gap |
| **H** | Enrichment note shipped; no false multi-\(k\) claims |
| **I** | One-command rebuild healthy |
| **J** | Citeable paper-shaped core artefact |
| **K** | \(k\)-taxonomy + regenerable multi-seed table |
| **L** | HQCC seed dictionary frozen |
| **M** | Geometric multi-\(k\) found **or** type-wise obstruction |
| **N** | Further Hurwitz models beyond \(3A^4\) scaffold |
| **O** | Obstruction ledger complete for known routes |
| **P** | Library identities re-verify Stage A1 alone |
| **Q** | Base-change theory written fail-closed on descent |
| **R** | Trusted \(f_s\in\mathbb{Q}(s)[x]\) **or** recorded open with attack log |
| **S** | Non-\(A_n\) simple survey or example |
| **T** | Theorem ledger audit-clean |
| **U** | U4 naturality lemma; U5 scoped |
| **V** | Regression suite green |
| **W** | Outsider-readable report stack |
| **X** | Non-BJ shape identity or bounded negative |
| **Y** | KPI dashboard regenerated |
| **Z** | Success + parked registry current; A–Z ladder complete |
| **Fusion (optional)** | Geometric object recovers multi-\(k\) / HQCC lattice — **bonus** |

---

## Document map

| Doc | Role |
|-----|------|
| **`RESOLUTION_PATH.md`** | This file — Stage A–Z lock |
| **`RESONANT_NUMBER_THEORY.md`** | A+B empirical grounding (CONFIRMED) |
| **`STAGE_D_DENSITY.md`** | Stage D density/asymptotics (COMPLETE) |
| **`RESONANT_NUMBER_THEORY_CZ.md`** | C–Z inventory scorecard |
| `THEOREMS.md` | Theorem-grade ledger + open problems |
| `ENLARGED_SEED_CATALOGUE.md` | Multi-seed pure-even catalogue |
| `NONRIGID_HURWITZ_SEARCH.md` / `REALISE_3A4_SPECIALISE.md` | Envelope + paths |
| `FUSION_GAP.md` | Geometric fusion still open |
| `HILBERT_MODULAR_A5.md` | Speculative √5 enrichment |
| `GENERATIVE_REACH.md` | Stage A3 / G scaffold |
| `REPORT.md` | Programme status report |

### Scripts

| Script | Stages |
|--------|--------|
| `resonant_number_theory_ab.py` | A + B (DIG / GROW / BUILD) |
| `resonant_number_theory_cz.py` | C–Z inventory + KPI |
| `build_all.py` | Integration rebuild |
| `generative_reach.py` | A3 / G |
| `enlarge_seed_catalogue.py` | A1 / K |
| `realise_3a4_specialise.py` | A1 paths / M control |

---

_Generated as programme resolution path — Resonant Algebra / Resonant Galois — Stages A–Z complete as ledger._
