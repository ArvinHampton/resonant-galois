# Research roadmap — post-publish stance

_Do not reopen settled layers. Citable centre = pure-even multi-\(k\) + four-face organising principle. Necessity paused. Canonical T3 production-locked._

---

## Status baseline (do not reopen)

| Layer | Status |
|-------|--------|
| Pure-even multi-\(k\) arithmetic | **Finished + published core** (`PURE_EVEN_MULTI_K.md`) |
| Four-face organising principle | **Locked** (`TERNARY_ORGANIZING_PRINCIPLE.md`) |
| Stage D density / height | **Locked** (`STAGE_D_DENSITY.md`) |
| Matrix templates + evenness obstruction | **Documented** (`HQCC_MATRIX_TEMPLATES.md`) |
| Tier 1.1 identical-square subclass | **Deepened / locked** (`TIER11_DEEPEN.md`) — no Crit-2 fragment |
| Tier 1.2 Candidate C functor | **First cut done** (`CANDIDATE_C_FUNCTOR.md`) — secondary if fusion deferred |
| Tier 1.3 / blowup genus of \(P\) | **\(g=1\)** under ordinary accounting (`GENUS_P_BLOWUP.md`) |
| Criterion 3 sign character | **Deepened** — no forcing \(\chi\) (`CRITERION3_DEEPEN.md`) |
| Arboreal T₃ vs catalogue | **Probed** (`ARBOREAL_T3.md`) — design-consistent, not necessity |
| **Canonical T3** (residue-1: \((4n+2)//3\)) | **Production lock** — do not reopen coefficient choice |
| **Necessity theorem** | **Paused** as citable claim (`NECESSITY_THEOREM.md`) |
| **Geometric multi-\(k\) fusion** | **Principal open problem** (`GEOMETRIC_MULTI_K_FUSION.md`) |
| Physical S²-11DM²ET-X claims | **Separate** (not Galois inputs) |

**Slogan:** generative success ≠ forced alternating monodromy from HQCC axioms.

---

## Citable package (what to point outsiders at)

1. **`PURE_EVEN_MULTI_K.md`** — theorem-grade BJ / pure-even / envelope / paths; HQCC as lattice only.  
2. **`TERNARY_ORGANIZING_PRINCIPLE.md`** — four faces; structural reading.  
3. **`THEOREMS.md`** — ledger.  
4. Verification: `lib/lemmas.py`, `pure_even_specialisations.py`, `stage_d_density.py`.  
5. Negative control: `RIGID_FIBRE_T3.md`.

---

## Active priorities (after pause)

### Tier G — Geometric multi-\(k\) fusion (**principal next problem**)

Problem card: **`GEOMETRIC_MULTI_K_FUSION.md`**.  
Ledger: `FUSION_GAP.md`. Avenues 1–7: `AVENUE_RANK_EXECUTE.md`.

| # | Route | Status | Goal |
|---|-------|--------|------|
| G1 | \(3A^4\) geometric model + Hilbert to catalogue | **Cuts done:** triple-root, Tschirnhaus, **param-field resolvent** (`G1_PARAM_FIELD_RESOLVENT.md`) — Norm deg 10 at \(s=-1\); 0/15 seeds divide Norm; multi-sheet cat=0 | Domain Möbius; cubic τ on deg-10; G2/G3 |
| G2 | Other g=0 Nielsen types (\(2A3A^3\), \(2A^2 3A^2\)) | **First cut done** (`G2_NIELSEN_G0.md`) — covers 14/14×3 ansätze; catalogue hits **0**; multi-\(k\) False | Domain Möbius; \(2A3A^2 5A\); G3 |
| G3 | Geometric monodromy + explicit pure-even model | **Corrected lock:** multiset **`5A²5B²`**, lift **−1**, braid orbit **size 12** (common basepoint); pure-even path = \(\mathbb{Q}(t)\)-model; multi-k Hilbert **True** (`EXPLICIT_5A4_EQUATION.md`, `EXPLICIT_5A4_REFINEMENTS.md`) | Optional: full chart of size-12 component |
| G– | Rigid \(\varphi\) surgery / \(\mathbb{Q}(\sqrt{5})\) descent | **Ruled out / side-only** | Do not retry as main line |

**Success:** construction **or** sharp obstruction. Necessity stays open and orthogonal.

### Tier A — Strengthen the published centre

| # | Task | Status |
|---|------|--------|
| A1 | Pure-even multi-\(k\) paper / manuscript | **Done** in-repo (`PURE_EVEN_MULTI_K.md`) |
| A2 | Catalogue secondary invariants on \(L_0\) | **Done** (`L0_SECONDARY_INVARIANTS.md`) |
| A3 | Outsider verification package (curated scripts + README) | Partial (scripts exist) |

### Tier L — Ternary lattice as theorem-bearing object

Frame: **`TERNARY_LATTICE_DIRECTIONS.md`**. Three roles of \(L_0\): specialisation / Mestre \(t\) / template coords.

| Dir | Goal | Status |
|----:|------|--------|
| **1** | Secondary invariants \(v_3\), disc primes, height, Gal | **Done** |
| **2** | Resonant monoid / saturation | **Done** (`L0_MONOID_SATURATION.md`) |
| **3** | Unify PE ↔ B avatars on \(L_0\) | **First cut** (`L0_PE_B_UNIFY.md`) — no canonical \(\Phi\) |
| **4** | Mestre orbit graph | **Done** (`L0_MESTRE_ORBIT.md`) — 24/24 seeds R-space; 480/480 even lattice-\(t\) |
| **5** | Necessity-facing avatar | Optional / paused |

### Tier A+ — New algebraic ideas (**executed**)

Doc: **`NEW_ALGEBRAIC_IDEAS.md`** / `new_algebraic_ideas.py`.

| Idea | Rank | Result |
|------|------|--------|
| **A** Mestre on HQCC seeds | Primary | **HIT:** dim-1 \(R\)-space; `shift_y_tR` families disc□ in \(\mathbb{Q}(t)\); sample \(A_5\) |
| **F** Embed into \(T\) | Primary | **HIT:** seeds + depressed Mestre specs embed; BJ-embed classical; non-BJ via B |
| **B** Non-BJ deg-1 \(x^5+75x^3+Ax^2+3A\) | Secondary | **HIT:** disc \(=324A^2(A^2+84375)^2\); lattice \(A\Rightarrow A_5\); \(d=-75\) in \(T\) |
| **C** New matrix avatar | If A+F fail | No identical-square-by-construction avatar yet |
| **D** Icosahedral params | Low | Sporadic disc□=0 in scan |
| **E** T₃-native polys | Low | Orbit polys not systematic \(A_5\) |

**Do not retry:** more cuts of same \(T\); \(F\to T\) disc□ hope; rigid \(\varphi\) surgery; Collatz⇒evenness without pure-even.

**Follow-through (executed):** `next_mestre_b_avatar.py` → **`MESTRE_FLAGSHIP_PT.md`**, **`B_EMBED_LATTICE.md`**, **`EVENNESS_AVATAR.md`**.

| Item | Result |
|------|--------|
| Closed-form flagship \(P_t\) | Explicit coeffs; disc□ in \(\mathbb{Q}(t)\); \(t=0\) = seed; \(t\in\{0,\pm1,2,3,5,7,9,61,80\}\) all sample \(A_5\) |
| B-embed lattice \(bc=72A\) | **104** unique lattice \(A\), all disc□; **50** checked \(A_5\) (model \(3,9,27,61,80,243,539,\ldots\)) |
| Evenness avatars | PE: \(T(0,-\alpha,k,0,0,1)\) disc\((256\alpha^2 m)^2\); B: \(T(-A,b,72A/b,-75,0,0)\) disc\((18A(A^2+84375))^2\); beyond BJ |

**Still open on this track:** HQCC-axiom *naming* of the embed relations (necessity); not required for generative use.

### Tier S — Secondary (if geometric fusion deferred)

| # | Task | Status |
|---|------|--------|
| S1 | Enlarge generative machine to \(A_7+\) (pure-even / lattice) | Open, secondary |
| S2 | Explicit design-mirror functor (Candidate C) binary→ternary | First cut done; write/test if deferred |
| S3 | Sharpen effective HIT / density inside pure-even families | Stage D baseline locked; effective constants open |

### Tier B — Optional depth (not blocking)

| # | Task | Status |
|---|------|--------|
| B1 | Criterion 3 further characters | Paused with necessity |
| B2 | Genus of \(P\) | \(g=1\) ordinary accounting (`GENUS_P_BLOWUP.md`) — known block on direct eliminant |
| B3 | Arboreal T₃ | Probe (`ARBOREAL_T3.md`) |
| B4 | Non-classical \(\mathcal{R}\)-coefficients | Low priority |

### Tier C — Necessity (explicitly open research)

Criteria 1–3 remain open. Mestre + non-BJ families **enlarge the generative machine** inside/near \(T\) but do **not** by themselves prove every HQCC axiom object has Gal \(A_n\). Resume necessity only if embed relations are **named by HQCC axioms** without classical evenness ansätze.

---

## Probe results (this package cycle)

| Probe | Result |
|-------|--------|
| Crit 2 / Tier 1.1 | Beyond BJ: dim 1 homog only; HQCC naming fails |
| Crit 3 deepen | No rate-1 ternary/HQCC character on unrestricted \(T\) |
| Genus \(P\) blowup | \(g=1\) if ordinary sings at \((1,1)\) and \([1:0:0],[0:1:0]\) |
| Arboreal T₃ | Evenness from pure-even identity; path monodromy encoding-dependent; Frob TV(B,C)≈0.04 |

---

## What not to do

| Do not | Why |
|--------|-----|
| Re-litigate pure-even “needing” HQCC for disc identity | Settled: **no** |
| Surgery on rigid \(\varphi/\mathbb{Q}\) for even fibres | Ruled out |
| Treat necessity as unfinished homework for the centre | **Paused** — open research |
| Expand catalogues without a new invariant | Noise |
| Re-run same linear cuts of \(T\) for identical-square | Exhausted |

---

## Priority table (now)

| Track | Priority |
|-------|:--------:|
| Publish / cite pure-even + organising principle | **Centre** (finished) |
| **Geometric multi-\(k\) fusion** (G1–G3) | **Principal next** |
| Catalogue invariants + verification package | High (support) |
| Secondary S1–S3 if fusion deferred | Medium |
| Arboreal deeper models | Low |
| Necessity Crit 1–3 | **Paused** |
