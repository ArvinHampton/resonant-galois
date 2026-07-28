# Theorem-grade ledger — Resonant Galois Programme

## Programme position (locked)

> **Arithmetic foundations are theorem-grade.**  
> **Citable centre of Resonant Algebra: pure-even multi-\(k\) theory** (theorems 1–6 below + catalogue).  
> **Explicit rigid cover \(\varphi/\mathbb{Q}\) with monodromy \(A_5\) is done.**  
> **Mild / further surgery on the same rigid \(\varphi/\mathbb{Q}\) is ruled out**  
> (disc \(=5\cdot\square\) on rational fibres ⇒ always odd when irreducible).  
> **BJ pencils** give equation-level inclusion of HQCC seeds.  
> **U4** gives a canonical residue-to-braid map on the present cover.  
>
> **Principal open fusion problem (structural, not required for A1 success):** a **single object** that is simultaneously a  
> **pure geometric \(A_5\) family** and **recovers the HQCC seed lattice as Hilbert specialisations**.
>
> **Citable centre:** pure-even multi-\(k\) (**`PURE_EVEN_MULTI_K.md`**) + four-face organising principle (**`TERNARY_ORGANIZING_PRINCIPLE.md`**).
>
> **Necessity theorem (open / paused as citable claim):** alternating monodromy as a **forced** consequence of HQCC / resonant axioms (Criteria 1–3). See **`NECESSITY_THEOREM.md`**. Catalogues and pure-even slices are **not** necessity.
>
> **Organizing principle (structural reading):** ternary branching in dynamics / lattice / matrices / Galois aims at \(A_n=\langle\text{3-cycles}\rangle\); pure-even locks the sign.

### Realistic resolution path

See **`RESOLUTION_PATH.md`**. Necessity is **not** Stage A centre; pure-even multi-\(k\) is.

| Stage | Goal |
|-------|------|
| **A (near-term)** | Secure mathematical core: multi-\(k\) as published centre; structural criteria attacks; generative reach beyond \(A_5\) |
| **B (medium-term)** | Independent arithmetic predictions checkable by outsiders |
| **D (medium)** | Density / asymptotics — **COMPLETE** (`STAGE_D_DENSITY.md`): D3 proved height; D1–D2 conjectures + tables |

### Remaining viable routes (geometric fusion only)

1. **Non-rigid \(A_5\) families** (Hurwitz / resolvents) — arithmetic multi-\(k\) **done**; Nielsen multi-\(k\) open.  
2. **Base change** \(K=\mathbb{Q}(\sqrt5)\) (+ Hilbert modular enrichment) — side route.  
3. **Other geometry** — open / high effort.

**Route 1 — pure-even multi-\(k\) (citable core)** — LSW, envelope, paths:

\[
f_t^{\mathrm{LSW}}=x^5+(t^2-3125)x-4(t^2-3125),\qquad
\operatorname{disc}=(16(t^2-3125)^2 t)^2.
\]

General \(k\)-slice / envelope: \(\alpha=256m^2-3125k^4/256\), \(\beta=k\alpha\), disc \((256\alpha^2 m)^2\).  
Cross-\(k\) paths join flagship (\(k=-8/5\)), classical (\(k=4/5\)), LSW (\(k=-4\)), ….

| Layer | Status |
|-------|--------|
| **Arithmetic / multi-\(k\)** | **Theorem-grade — programme centre** |
| **Geometry (\(\varphi\), monodromy \(A_5\))** | **Done** |
| **Surgery on same rigid \(\varphi/\mathbb{Q}\)** | **Ruled out** |
| **BJ pencils (equation-level seeds)** | **Done** (not pure even generic) |
| **U4 residue→braid** | **Canonical** |
| **Fused pure \(A_5\) + Hilbert seed lattice** | **Open — structural** |
| **Necessity theorem (Crit 1–3)** | **Open / paused** (`NECESSITY_THEOREM.md`) |
| **Publishable pure-even core** | **`PURE_EVEN_MULTI_K.md`** |
| **Genus of \(P(q,w)\)** | **\(g=1\)** ordinary accounting (`GENUS_P_BLOWUP.md`) |
| **Crit 3 deepen** | No forcing character (`CRITERION3_DEEPEN.md`) |
| **Arboreal T₃** | Probe only (`ARBOREAL_T3.md`) |
| **HQCC templates \(M\), \(T\), BJ-embed** | Explored (`HQCC_MATRIX_TEMPLATES.md`): structure ⇏ disc □ |
| **Resonant algebraic closure** | Candidate A = lattice+pure-even (`RESONANT_ALGEBRAIC_CLOSURE.md`); not necessity |
| **Generative reach beyond \(A_5\)** | **Stage A3** (`GENERATIVE_REACH.md`) |
| **Catalogues / pipeline** | Finished |

See **`RESOLUTION_PATH.md`**, **`FUSION_GAP.md`**, **`GEOMETRIC_RIGID_DEFORM.md`**.

### Vocabulary: HQCC seed

An **HQCC seed** is a BJ polynomial \(s(x)=x^5+\alpha x+\beta\in\mathbb{Z}[x]\) with:

1. **Lattice** — \(\alpha,\beta\) from the HQCC / resonant lattice (or short integer combinations of \(\{3,9,27,61,80,243,539,\ldots\}\));
2. **Arithmetic** — \(\operatorname{disc}(s)\) square and \(\mathrm{Gal}(s/\mathbb{Q})=A_5\);
3. **Homogenisation** — \(f_t=x^5+\alpha t^4 x+\beta t^5\) then has \(\operatorname{disc}(f_t)=t^{20}\operatorname{disc}(s)\) (even for all \(t\neq0\)).

**Not required:** geometric cover / branch-cycle origin; formula in \(\xi=2\cos(2\pi/539.9)\); full \(T_5\) outside BJ embed.

Full definition: **`HQCC_SEED.md`**. Flagship: \(s(x)=x^5-55x+88\) (\(88=61+3^3\)).

---

## What is now proved (arithmetic)

| # | Result | Criteria |
|---|--------|----------|
| **1** | \(\operatorname{disc}(x^5 + a x + b) = 256a^5 + 3125b^4\) (symbolic identity) | 2 + 3 |
| **2** | On the Bring–Jerrard class: even monodromy \(\iff\) that integer is a square; + Frobenius type \((3,1,1)\) \(\Rightarrow A_5\) | 2 + 3 |
| **3** | Homogenised classical family \(f_t = x^5 + 20 t^4 x + 16 t^5\): \(\operatorname{disc}(f_t) = t^{20}\cdot\operatorname{disc}(x^5+20x+16)\) is always a square for \(t\neq 0\) | 1 + 2 + 3 |
| **4** | Operational group-theoretic criterion for \(A_5\) (irr + disc² + type \((3,1,1)\)) | pipeline |
| **5** | **General homogenisation lemma.** For any integers \(\alpha,\beta\) with \(\operatorname{disc}(x^5+\alpha x+\beta)\) a square, the family \(f_t=x^5+\alpha t^4 x+\beta t^5\) has \(\operatorname{disc}(f_t)=t^{20}\operatorname{disc}(\mathrm{seed})\) a square for all \(t\neq 0\). Hence Gal even whenever \(f_t\) is irreducible. | 1 + 2 + 3 |
| **6** | **HQCC seeds exist** (def. `HQCC_SEED.md`). E.g. \((\alpha,\beta)=(-55,88)\), \((95,76)\), \((95,532)\), …. Flagship seed \(x^5-55x+88\) (\(88=61+3^3\)) → family \(x^5-55 t^4 x+88 t^5\). | 1 |
| **7** | **T5 embed.** Every such homogenised family realises inside the structural \(T_5\) template on the thin subclass \(d=0\), \(a=-ef\) (BJ embed), so Crit 2 on that subclass reduces to theorems 1–2. | 2 |
| **8** | **Pure-even \(k\)-slice.** For \(k\in\mathbb{Q}\setminus\{0\}\), \(\alpha(m)=256m^2-3125k^4/256\), \(\beta=k\alpha\) has disc \((256\alpha^2 m)^2\) identically in \(\mathbb{Q}(m)\). | multi-\(k\) core |
| **9** | **2-param envelope + cross-\(k\) paths.** Free \((m,s)\) with \(k=s\) is pure-even over \(\mathbb{Q}(m,s)\); rational paths \((m(u),k(u))\) join any two envelope seeds (different \(k\) allowed) with disc identically square. | multi-\(k\) core |
| **10** | **Disc height (Stage D3).** On pure-even \(k\)-slice, \(\operatorname{disc}=(256\alpha(m)^2 m)^2\); as \(\|m\|\to\infty\), \(\log\|\operatorname{disc}\|=10\log\|m\|+48\log 2+o(1)\). | Stage D |
| **11** | **Rigid fibre \(t=3\) negative control.** \(\mathrm{monic}(\varphi-3)=y^5-\frac52 y^4+\frac53 y^3-\frac12\); \(\operatorname{disc}=3125/36=5\cdot(25/6)^2\) not a square in \(\mathbb{Q}\); irr; \(\mathrm{Gal}=S_5\). Locks odd side of pure-even \(\leftrightarrow\) rigid contrast. | `RIGID_FIBRE_T3.md` |
| **12** | **Field-agnostic pure-even.** Over any field \(F\) with \(\mathrm{char}\notin\{2,5\}\), \(\alpha=256m^2-3125k^4/256\), \(\beta=k\alpha\) give \(\operatorname{disc}=(256\alpha^2 m)^2\) in \(F(m,k)\). Non-classical track may take \(m,k\in\mathcal{R}\). | `NONCLASSICAL_RESONANT_FIELD.md` |

### Conjectures with evidence (Stage D)

| # | Statement | Evidence |
|---|-----------|----------|
| **D1** | Positive irreducibility density on each multi-seed pure-even \(k\)-slice lattice \(L_k\) | irr rate 1.0 on large sampled \(Z\)-coefficient sets; even_fail=0 (`STAGE_D_DATA.json`) |
| **D2** | Chebotarev equidistribution in \(A_5\) on Gal=\(A_5\) fibres of pure-even slices | cycle-type histograms on 100 fibres / 4500 primes; class densities \(1/3,2/5,1/4\) predicted |

Proof sketches and code: `lib/lemmas.py`, `theorem_attack.py`, `hqcc_native.py`, `hqcc_strict_analysis.py`, `t5_disc_ideal.py`, `nonrigid_multi_seed.py`, `enlarge_seed_catalogue.py`, `realise_3a4_specialise.py`, `stage_d_density.py`.

### Empirical confirmation (arithmetic)

- Classical \(f_t\): **33/33** tested specialisations → \(\mathrm{Gal}=A_5\) (incl. model \(t\in\{3,9,61,80,243,539\}\)).
- HQCC flagship \(g_t=x^5-55 t^4 x+88 t^5\): **9/9** model \(t\) → \(\mathrm{Gal}=A_5\).

Model integers enter as **Hilbert specialisations of proved-even families**, not free coefficient search. On these classes the square-discriminant gate is an **identity**.

---

## Principal open problem — geometric cover

### What would count as a solution

A finite cover \(X\to Y\) (or a family \(F_t\in\mathbb{Q}(t)[x]\) with geometric monodromy)
built from **HQCC branch data** — Möbius maps \(\{n/3,\,3n\pm1,\,(4n+2)/3,\,(2n+1)/3\}\),
ternary / \(\mathrm{Ad}_{\mathrm{SO}(3)}\) structure, T-complementarity — such that:

1. the **geometric** monodromy group is \(A_n\) (or contains \(A_n\)) by a proof, not by
   specialisation census alone; and
2. model integers arise as specialisations (Hilbert), recovering catalogue hits where appropriate.

### What does *not* yet solve it

| Object | Why it is arithmetic, not geometric |
|--------|-------------------------------------|
| BJ / homogenised \(f_t\) with HQCC \((\alpha,\beta)\) | Coefficients in the lattice; no branch-cycle type prescribed by HQCC |
| Dense \(T_5/T_6\) lattice search | Generator + filters; base \(M\), \(T_6\) odd |
| Cubic resultants \(y^3-3sy-t\), \(x=y+m/y\) | Prefer \(A_3\); not yet an \(A_5\) cover |
| Möbius block char polys | Typically reducible / odd |

### Geometric Steps 1–2 (done)

- **Step 1:** Rigid tuples — **`GEOMETRIC_STEP1.md`**. Preferred `(3A,3A,5A)`; fallback `(2A,3A,5A)`.
- **Step 2:** Explicit Belyi covers — **`GEOMETRIC_STEP2.md`**.
  - Preferred over **\(\mathbb{Q}\):** \(\varphi(y)=6y^5-15y^4+10y^3\), monodromy **\(A_5\)**, passport (3,1,1)(3,1,1)(5).
  - Fallback over \(\mathbb{Q}(2^{1/5},3^{1/5})\): radical quintic Belyi, monodromy **\(A_5\)**, passport (3,1,1)(2,2,1)(5).
  - Branch points labelled by 9 Maths / T₃ / G₄ (dictionary level).
- **Queued attacks 1–3 done** (`QUEUED_ATTACKS.md`): Gal-level seed↔cover match only (no equation identity); resonant base \((3,61,539)\to(0,1,\infty)\); \(T_3\to A_5\) path functor scaffold.
- **From geometry** (`GAP_A_FROM_GEOMETRY.md`): deformations/pull-backs of \(\varphi\) + Tschirnhaus; **no** HQCC BJ seed recovered; lattice fibres of monic(\(\varphi-t\)) odd in scan.
- **Still open (deep):** pure geometric \(A_5\) family with Hilbert multi-seed recovery; naturality of \(T_3\) functor; \(\xi=2\cos(2\pi/539.9)\) modular identification.

### Secondary open items (subordinate to geometry)

- Full unrestricted \(T_5/T_6\): no axiom list forces disc² (evenness obstruction remains).
- Crit 3 invariant on unrestricted structural \(M\) (solved only on BJ / homogenised classes).
- Non-BJ thin class with \(\operatorname{disc}(\chi_{T_5})\equiv S^2\) identically.

---

## Resolution criteria — status map

| Criterion | Arithmetic status | Geometric status |
|-----------|-------------------|------------------|
| **1** Canonical HQCC object | Infinite HQCC-lattice families with proved even monodromy + empirical \(A_5\) | **Open** — proved monodromy cover from branch data |
| **2** Axioms ⇒ disc² + 3-cycles | Proved on BJ class + BJ-embed in \(T_5\); homogenisation lemma | Open for full structural templates |
| **3** Sign / ternary invariant | Proved on BJ + homogenised (disc formula) | Open for unrestricted monodromy representations |

---

## Paths

| Path | Content |
|------|---------|
| `THEOREMS.md` | This ledger (programme position) |
| `HQCC_SEED.md` | **Definition** of HQCC seed |
| `THEOREM_ATTACK.md` | Criteria 1–3 attack report |
| `HQCC_NATIVE.md` / `HQCC_STRICT.md` | HQCC-lattice families + provenance |
| `T5_DISC_IDEAL.md` | Template disc + BJ embed |
| `IMPLICATIONS.md` | Claims / non-claims |
| `RESOLUTION.md` | Build status + criteria gaps |

_Arithmetic closed at theorem grade for the infinite HQCC-lattice family class.
Geometry is the remaining resolution frontier._
