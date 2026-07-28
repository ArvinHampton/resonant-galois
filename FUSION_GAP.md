# Fused Criterion-1 gap

## Bottom line (locked)

1. **Mild surgery on the existing Belyi map remains insufficient** (already known).  
2. **BJ pencils** give **equation-level inclusion** of the HQCC seeds in a geometric family.  
3. **U4** supplies a **canonical residue-to-braid assignment** on the present cover.  
4. **Further surgery on the same rigid \(\varphi/\mathbb{Q}\)** is **ruled out**:  
   \(\operatorname{disc}(\mathrm{monic}(\varphi-t))=5\cdot(\mathrm{square})\) for \(t\in\mathbb{Q}\setminus\{0,1\}\)  
   ⇒ no even irreducible arithmetic specialisations (`GEOMETRIC_RIGID_DEFORM.md`).

> A **single object** that is simultaneously a **pure geometric \(A_5\) family** and  
> **recovers the HQCC seed lattice as Hilbert specialisations** is **still missing**.  
> That remains the **principal open fusion problem** (structural).  
>
> **Problem card (active):** `GEOMETRIC_MULTI_K_FUSION.md` — geometric multi-\(k\) fusion  
> is now the **principal next problem** after the finished pure-even multi-\(k\) centre.  
> Necessity stays open/paused; **Canonical T3** is the production dynamical baseline  
> (do not reopen residue-1 coefficients).  
>
> Under the **realistic resolution path** (`RESOLUTION_PATH.md`), programme success  
> for **Stage A** is defined by the **pure-even multi-\(k\) arithmetic core** (citable).  
> Closing or sharply obstructing geometric fusion is the strongest IG-adjacent advance  
> that does not require Necessity.
>
> **Negative control locked** (`RIGID_FIBRE_T3.md`): fibre \(t=3\) of \(\varphi\) has
> monic disc \(=5\cdot\square\) (odd), irr, \(\mathrm{Gal}=S_5\). Pure-even resonant
> slices remain the even/\(A_5\) side of the arithmetic distinction.

## Remaining viable routes only

1. **Non-rigid \(A_5\) families** (positive-dimensional Hurwitz spaces).  
   — **Arithmetic multi-\(k\) solved** (`NONRIGID_HURWITZ_SEARCH.md`): envelope
   \(\alpha=256m^2-3125s^4/256\), \(\beta=s\alpha\) pure-even over \(\mathbb{Q}(m,s)\);
   recovers all fixed-\(k\) slices; linear paths in \((m,k)\) give pure-even 1-param
   families through flagship\(\leftrightarrow\)LSW, flagship\(\leftrightarrow\)classical, etc.
   — **\(r=4\) Hurwitz strata** (`A5_HURWITZ_R4.md`): 19 filter-pass Nielsen types;
   g=0 shortlist includes \(3A^4\) (1 orbit), etc.
   — **All 4 steps executed** (`REALISE_3A4_SPECIALISE.md`): Nielsen \(3A^4\) ok;
   **Step 4 SUCCESS** for envelope paths (flagship\(\leftrightarrow\)classical/LSW/…).
   — **Genus lock** (`GENUS_3A4_LOCK.md`): \(3A^4\) reduced Hurwitz curve **genus 0**
   over \(\mathbb{Q}\), infinitely many rational points (Bailey–Fried; orbit size 18).
   — **Explicit model** (`EXPLICIT_3A4_RESOLVENT.md`): \(H^{\mathrm{rd}}\cong\mathbb{P}^1_s\);
   covers + exact \(s=-1\) over \(\mathbb{Q}(\sqrt{5})\); closed form \(f_s\in\mathbb{Q}(s)[x]\) open.
   — **Avenues 1–7 executed** (`AVENUE_RANK_EXECUTE.md`): geometric multi-\(k\) still open;
   arithmetic multi-\(k\) solid (envelope). Rigid/base-change routes blocked/side-only.
   \(\varphi\) (\(r=3\)) abandoned for \(\mathbb{Q}\)-fusion.
2. **Arithmetic over a number field \(K\)** in which the permanent factor 5 becomes a square.  
   — **Probed** (`K_SQRT5_EVEN.md`): over \(K=\mathbb{Q}(\sqrt{5})\) disc of monic(\(\varphi-t\)) is identically square; **153/153** rational \(t\) even-over-\(K\); **no descent** to even-over-\(\mathbb{Q}\); no lattice recovery of HQCC seeds. Side route only.  
   — **Hilbert modular / icosahedral enrichment** (`HILBERT_MODULAR_A5.md`): classical
   Hirzebruch–Klein link \(\mathbb{Q}(\sqrt{5})\) Hilbert modular surface \(\leftrightarrow A_5\);
   Klein \(A,B,C,D\) (weights 2,6,10,15). Natural home for the quadratic obstruction
   (\(\varphi\)-disc factor 5; 3A⁴ cover at \(s=-1\)). High-effort / speculative for multi-\(k\)
   over \(\mathbb{Q}\); **not** a short-cut past arithmetic multi-\(k\) or the 3A⁴ resolvent.
3. **Entirely different geometric constructions** not arising from the rigid triples already examined.

### Non-rigid route — multi-seed pure-even slices (`ENLARGED_SEED_CATALOGUE.md`)

Enlarged HQCC lattice BJ scan + grouping by \(k=\beta/\alpha\):

| quantity | value |
|----------|------:|
| Unique A5 BJ seeds | **60** |
| Distinct ratios \(k=\beta/\alpha\) | **16** |
| **Multi-seed pure-even slices** | **10** |

General pure-even family on the ray \(\beta=k\alpha\):

\[
\alpha(m)=256 m^2 - \tfrac{3125\,k^4}{256},\qquad \beta(m)=k\cdot\alpha(m),\qquad
\operatorname{disc}=(256\,\alpha^2 m)^2.
\]

| \(k=\beta/\alpha\) | # A5 seeds | Notable members |
|-------------------|----------:|-----------------|
| \(-4\) (LSW) | **11** | \((-100,400)\), \((124,-496)\), … |
| \(4\) | **11** | sign-flips of LSW ray |
| \(-12/5\) | **5** | \((-180,432)\), \((220,-528)\), … |
| \(12/5\) | **5** | sign-flips |
| **\(-8/5\)** (flagship) | **4** | **flagship** \((-55,88)\), \((145,-232)\), \((320,-512)\), \((1145,-1832)\) |
| \(8/5\) | **4** | flagship flip ray |
| \(4/5\) (classical) | **4** | classical \((20,16)\), \((95,76)\), … |
| \(-4/5\) | **4** | classical flip ray |
| \(\pm16/5\) | **3** each | includes \((-55,\pm176)\) |

**Resolved (arithmetic multi-seed):** flagship now sits on a pure-even family with **three further A5 seeds** (same \(k=-8/5\)). Classical and \((95,76)\) cohabit \(k=4/5\). LSW expanded from 2 to 11 seeds on \(k=-4\).

**Still open:** pure-even family joining seeds of **different** \(k\) (e.g. flagship \(k=-8/5\) with classical \(k=4/5\)); geometric fusion with rigid \(\varphi\).

Earlier notes (`NONRIGID_A5_FAMILIES.md`, `NONRIGID_MULTI_SEED.md`): LSW + flagship slice discovery; Mestre/quad scans empty in bounds.

The principal open fusion problem remains open; further surgery on the same rigid \(\varphi/\mathbb{Q}\) is ruled out.

---

## Status detail

Arithmetic and geometric **existence** pieces are in hand. What is missing is a
**single fused Criterion-1 object** that is both:

- geometric (HQCC-labelled branch data / monodromy \(A_5\)), and  
- arithmetic (specialises to the known HQCC seeds).

## Either / or (both still open at full strength)

### Gap A — Hilbert recovery of seeds

Find a geometric family \(F_u\in\mathbb{Q}(u)[x]\) whose **generic** fibre has monodromy
\(A_5\) (pure even geometric monodromy) and whose Hilbert specialisations recover

\[
x^5 - 55x + 88,\quad
x^5 + 95x \pm 76,\quad \ldots
\]

| Progress | Limitation |
|----------|------------|
| Affine fibres of \(\varphi\) ≠ seeds | Mild \(\varphi\)-surgery blocked (BJ obstruction) |
| Linear BJ pencils through seeds | Seeds at \(t=0,1\); disc **not** □ in \(\mathbb{Q}(t)\) |
| Homogenised rays through one seed | Disc \(=t^{20}\times\)square in \(\mathbb{Q}(t)\); A5 along ray — **only pure-even class found** |
| Systematic scan of 38 BJ families (`GAP_A_BJ_FAMILIES.md`) | No pure-even multi-seed family; exploratory forms fail disc □ test |
| Even-surface curves (`GAP_A_EVEN_SURFACE.md`, complete) | Rays 12; linear 66/0; space 264/0; monomial 45/0; Bezier 1350/0; ratquad 44145/0; planes 20/0 low-deg factors. **No multi-seed pure-even rational curve found.** |
| **From geometry** (`GAP_A_FROM_GEOMETRY.md`) | Start from \(\varphi\): pull-backs/Tschirnhaus; lattice fibres odd; **0** seed hits. |
| **Other rigid triples + mild A5 moves** (`GEOMETRIC_RIGID_DEFORM.md`) | (3A,3A,5*), (2A,3A,5*); pull-backs; dense even search. **Theorem:** monic(\(\varphi-t\)) has \(\mathrm{disc}=5\cdot(\mathrm{square})\) over \(\mathbb{Q}(t)\) ⇒ **no** even irr. rational specialisation. Absolute rigidity ⇒ no continuous A5-preserving deformation of one triple. |

### Gap B — natural \(T_3\to\) braid functor

| Progress | Limitation |
|----------|------------|
| **U4:** residue \(2\mapsto\sigma_\infty=(\sigma_0\sigma_1)^{-1}\) | Canonical on 3-point \(\pi_1\); formal naturality still open |
| **U5:** 4-point cover dictionary | Generates \(A_5\); not absolutely rigid; no explicit \(\mathbb{Q}\)-equation yet |

## What is *not* missing

- BJ disc formula and homogenisation lemma  
- Existence of HQCC seeds and proved-even infinite families  
- Explicit \(\varphi/\mathbb{Q}\) with geometric monodromy \(A_5\)  
- Resonant base chart \((3,61,539)\to(0,1,\infty)\)  
- BJ pencils with equation-level seed inclusion  
- U4 canonical residue-2 assignment  
- Catalogues and operational pipeline  

## Paths

| Doc | Role |
|-----|------|
| `THEOREMS.md` | Master ledger |
| `FUSION_NEXT.md` | BJ pencil + U4/U5 probes |
| `FUSION_DEPTH.md` | Mild \(\varphi\)-twist obstruction |
| `QUEUED_ATTACKS.md` | Specialisation / base / scaffold functor |
| `GEOMETRIC_STEP2.md` | Cover \(\varphi\) |
| `GAP_A_FROM_GEOMETRY.md` | Start from \(\varphi\); Tschirnhaus; seed test |
| `GAP_A_EVEN_SURFACE.md` | Classical BJ even-surface curve search (stopped as main line) |
| `HQCC_SEED.md` | Seed definition |

_Locked: principal open problem = pure geometric \(A_5\) family with Hilbert recovery of the HQCC seed lattice.
Diophantine search on the classical BJ even surface is **not** the active attack path; geometry-first from \(\varphi\) is._
