# Geometric multi-\(k\) fusion — principal open problem

**Resonant Galois Programme · problem card**  
**Status:** Open (active priority).  
**Depends on:** finished pure-even multi-\(k\) centre (`PURE_EVEN_MULTI_K.md`).  
**Does not depend on:** Necessity theorem (remains open / paused).  
**Dynamical baseline:** **Canonical T3** is locked as the production iterative map  
(\((4n+2)//3\) on residue 1). Do **not** reopen the coefficient choice (T4121 stays experimental / FPGA only).

---

## 1. Statement

Construct a **single geometric object** that is at once:

1. a **pure Nielsen-labelled \(A_5\) family**  
   (discriminant a square in the parameter ring / pure-even geometric monodromy), and  
2. **recovers the HQCC / pure-even seed lattice as Hilbert specialisations**.

In symbols: find \(F_u\in\mathbb{Q}(u)[x]\) (or a mild extension of \(\mathbb{Q}\)) such that

- the geometric monodromy of the Galois closure over \(\overline{\mathbb{Q}}(u)\) is \(A_5\)
  (equivalently: a Nielsen-realised class with even arithmetic specialisations), and  
- there exist specialisations \(u\mapsto u_i\in\mathbb{Q}\) recovering the catalogue seeds, e.g.
  \[
  x^5-55x+88,\quad
  x^5+20x+16,\quad
  x^5+95x\pm76,\quad
  \ldots
  \]
  and, ideally, fibres on **two or more** fixed-\(k\) pure-even rays
  (flagship \(k=-8/5\), classical \(k=4/5\), LSW \(k=-4\), …).

**Success criterion (strong):** multi-\(k\) Hilbert hits from an explicit Nielsen-class equation.  
**Success criterion (sharp obstruction):** prove that no genus-0 / rational Nielsen model in the
programme shortlist can recover two distinct catalogue \(k\)-slices — a clean negative theorem
is publishable value equal to a construction.

Closing or **sharply obstructing** this gap is the strongest Inverse-Galois-adjacent advance
available that does **not** require solving the full Necessity theorem.

---

## 2. Why this is the natural next problem

| Fact | Doc |
|------|-----|
| Arithmetic multi-\(k\) is **finished and citable** | `PURE_EVEN_MULTI_K.md` |
| Geometric multi-\(k\) is the **explicit fusion gap** | `FUSION_GAP.md`, `A5_HURWITZ_R4.md` |
| Avenues 1–7 executed: geometric Nielsen hit = **False** | `AVENUE_RANK_EXECUTE.md` |
| Rigid \(\varphi/\mathbb{Q}\) surgery ruled out (disc \(=5\cdot\square\)) | `GEOMETRIC_RIGID_DEFORM.md`, `RIGID_FIBRE_T3.md` |
| Linear BJ pencils include seeds equation-level but disc **not** □ in \(\mathbb{Q}(t)\) | `FUSION_NEXT.md` |
| Even-surface multi-seed pure-even rational curves: **none found** | `GAP_A_EVEN_SURFACE.md` |

Arithmetic and geometric **existence** pieces sit side by side; the missing object is their
**fusion** into one Criterion-1-grade geometric family with lattice Hilbert recovery.

---

## 3. What is already locked (do not reopen)

| Layer | Status |
|-------|--------|
| Pure-even multi-\(k\) (disc identity, \(k\)-slices, envelope, cross-\(k\) paths) | **Finished** |
| Four-face organising principle | **Locked** |
| Stage D density / height | **Locked** |
| Necessity Criteria 1–3 | **Open / paused** — not required for this problem |
| Canonical T3 coefficient on residue 1: \((4n+2)//3\) | **Production lock** |
| Mild surgery on rigid \(\varphi/\mathbb{Q}\) | **Ruled out** |
| Diophantine search on classical BJ even surface as main line | **Stopped** |

---

## 4. Concrete attack routes (bounded)

The most direct geometric route already examined for the eliminant \(P(q,w)\) is
**blocked by positive genus** (`GENUS_P_BLOWUP.md`: \(g=1\) under ordinary accounting).
Remaining options are concrete:

### Route G1 — Closed form for \(3A^4\) (highest leverage)

| Item | Status |
|------|--------|
| Nielsen class \(\mathrm{Ni}(A_5,C_3^4)\) | Locked; 1 braid orbit, size 18 |
| Reduced Hurwitz curve | **Genus 0** over \(\mathbb{Q}\); infinitely many \(\mathbb{Q}\)-points (`GENUS_3A4_LOCK.md`) |
| Model \(H^{\mathrm{rd}}\cong\mathbb{P}^1_s\) | Locked (`EXPLICIT_3A4_RESOLVENT.md`) |
| Exact fibre at \(s=-1\) | Over \(\mathbb{Q}(\sqrt{5})\); descent issues |
| Closed form \(f_s\in\mathbb{Q}(s)[x]\) | **Not achieved** (this normal form) |
| Multi-\(k\) catalogue Hilbert hits from Nielsen equation | **None** in G1 cut |

**G1 first cut executed** (`G1_3A4_TRIPLE_ROOT.md` / `g1_3a4_triple_root.py`):

| Finding | Result |
|---------|--------|
| Rational pts on eliminant \(P(q,w)\) (height ≤ 24) | 4, all **degenerate**; **0** non-deg |
| Physical points on \(P\) | \(q,w=\pm1/\sqrt5\) on \(P\) (s=−1 fibre) |
| Catalogue seeds reverse-compatible with \(N-tD\) | **16/16** priority seeds (2 sols each) |
| Reverse sols on triple-root locus | **0/16** |
| Match to s=−1 params \((p_2,\sigma,\pi)=(-1,0,-1/25)\) | **0** |
| Norm of s=−1 fibre → monic deg-5 over \(\mathbb{Q}\) | **Never** (deg-10 irr. or quadratic factors) |
| Geometric multi-\(k\) | **False** |
| Arithmetic multi-\(k\) control | **True** |

**G1 Tschirnhaus cut** (`G1_SEED_TSCHIRNHAUS.md` / `g1_seed_tschirnhaus.py`):

| Track | Result |
|-------|--------|
| T1: s=−1 fibre + quadratic τ | **0/10** seeds |
| T2: locus + τ reverse (NumPy LM) | **0/10** seeds (best residual ≃ 0.24 for flagship, not 0) |
| T3: forward covers + τ → BJ | 21 BJ fibres; **0** catalogue hits |
| Translation-only reverse (control) | **10/10** still |

**Conclusion:** quadratic fibre Tschirnhaus \(\tau=c_0+c_1 y+c_2 y^2\) is **not enough** to place pure-even catalogue seeds on this \(3A^4\) normal-form locus. Obstruction is deeper than the G1 translation reverse.

**G1 parameter-field resolvent** (`G1_PARAM_FIELD_RESOLVENT.md` / `g1_param_field_resolvent.py`):

| Item | Result |
|------|--------|
| Model at \(s=-1\) | \(f\in\mathbb{Q}(\sqrt5)(t)[y]\), \([K:\mathbb{Q}]=2\) |
| Norm over \(\mathbb{Q}(t)\) | deg **10** \(=5\cdot[K:\mathbb{Q}]\) |
| R1: catalogue seed divides Norm | **0/15** |
| R2: forward Norm factors | degs \(\{10,2,1\}\) only — **no** deg-5; cat=0 |
| R3: multi-sheet product norms | 48 Z-products; cat=0; multi-\(k\)=False |
| Sheet count proxy | max ≈ **2** (matches quadratic field) |

**Conclusion:** accepting \(f\in K(s)[x]\) and norming to \(\mathbb{Q}\) does **not** recover the pure-even catalogue for this \(3A^4\) normal form. The Q-model is genuinely deg-10 at \(s=-1\) and still misses the lattice seeds.

**Next work on G1/G2/G3:** (i) domain Möbius / different normal form; (ii) cubic Tschirnhaus on the deg-10 Q-model; (iii) G2 other g=0 Nielsen types; (iv) G3 envelope monodromy ID.

### Route G2 — Other Nielsen types aimed at genus-0 charts

Shortlist from `A5_HURWITZ_R4.md` / Avenue 2:

- \(2A\,3A^3\) (orbit size 96, g=0 lookup)
- \(2A^2\,3A^2\) (orbit size 108, g=0 lookup)
- \(2A\,3A^2\,5A\), \(3A^3\,5A\), … (not yet cut)

**G2 first cut executed** (`G2_NIELSEN_G0.md` / `g2_nielsen_g0.py`):

| Item | Result |
|------|--------|
| Nielsen samples \(2A3A^3\), \(2A^2 3A^2\) | Found generating tuples (prior orbits 96 / 108) |
| Cover ansätze Newton at rational \(s\) | **14/14** each (3 placements) |
| Catalogue Hilbert hits | **0** |
| Geometric multi-\(k\) | **False** |
| Arithmetic multi-\(k\) control | **True** |

**Conclusion:** profile-correct degree-5 covers for both shortlist types realise over many rational \(s\), but their fibres do not recover the pure-even multi-seed catalogue in this scan.

### Route G3 — Geometric identification of the pure-even envelope

**G3 + G3b + G3c executed** (`G3_ENVELOPE_MONODROMY.md`, `G3B_5A5B_BRAID_LIFT.md`, `EXPLICIT_5A4_EQUATION.md`):

| Item | Result |
|------|--------|
| Envelope disc identity | **True** (theorem-grade) |
| Multi-k paths | Arithmetic multi-\(k\) **True** |
| Local monodromy | four 5-cycles; ∞ = id |
| **5A vs 5B split** | **all four are 5A** → multiset **`5A⁴`** |
| vs ternary shortlist | **excluded** |
| **Braid orbits on Ni(A₅, 5A⁴)** | **1 orbit**, size **10** |
| **Lift invariant (SL(2,5))** | **`+1`** on all 600 Nielsen tuples |
| **Explicit \(\mathbb{Q}(t)\)-model** | pure-even `path_flag_classical` (and siblings) |
| Monodromy re-check on model | **5A⁴** |
| Hilbert catalogue multi-\(k\) | **True** (flagship, classical, LSW, …) |

**G3d refinements** (`EXPLICIT_5A4_REFINEMENTS.md`): common-basepoint monodromy
**corrects** the multiset. Independent-loop “all 5A” labels were sheet-label gauge
artefacts. With a fixed base point:

| Item | Value |
|------|-------|
| Multiset | **`5A²5B²`** |
| Lift invariant | **`−1`** |
| Braid orbit | **size 12** (of two orbits 30/+1 and 12/−1) |
| Reduced genus (5A⁴ size-10, double-twist RH) | **g = 0** (combinatorial type) |
| Branch cross-ratio / j | λ ≈ 1/2, j ≈ 1728 |

**Type-level fusion (corrected):** pure-even multi-\(k\) paths are explicit
\(\mathbb{Q}(t)\)-models of **`Ni(A₅, 5A²5B²)` lift −1**, with multi-seed Hilbert hits.

### Ruled out / side-only

| Route | Why |
|-------|-----|
| Rigid \(\varphi\) fibre surgery | Permanent factor 5; no even irr. \(\mathbb{Q}\)-specialisations |
| Base change \(\mathbb{Q}(\sqrt{5})\) alone | Even over \(K\); no descent; no lattice recovery |
| Linear BJ pencils alone | Seeds at endpoints; disc not identically square |
| Same linear cuts of matrix template \(T\) | Exhausted |

---

## 5. Secondary options (if geometric fusion is deferred)

Use only if G1–G3 stall after a bounded effort, or in parallel at low cost:

| # | Option | Notes |
|---|--------|-------|
| S1 | Enlarge the generative machine to \(A_7+\) by pure-even / lattice methods | Same disc-evenness + Frobenius strategy; new envelopes |
| S2 | Design-mirror functor (Candidate C) — binary data → ternary lattice / matrix templates | First cut exists (`CANDIDATE_C_FUNCTOR.md`); write/test explicit maps |
| S3 | Sharpen effective Hilbert-irreducibility / density inside pure-even families | Strengthen Stage D; effective constants; Chebotarev rates |

These do **not** close the fusion gap; they enlarge or polish the arithmetic centre.

---

## 6. Explicit non-goals

- Do **not** reopen the pure-even disc identity or claim it “needs” HQCC / T3.
- Do **not** treat Necessity as unfinished homework for the centre.
- Do **not** change production dynamics off Canonical T3 for Galois work that needs an iterative map.
- Do **not** expand catalogues without a new geometric invariant or Nielsen equation.

---

## 7. Document map

| Doc | Role |
|-----|------|
| **This file** | Locked problem card + attack list |
| `FUSION_GAP.md` | Historical gap ledger (Gaps A/B) |
| `PURE_EVEN_MULTI_K.md` | Finished arithmetic centre |
| `A5_HURWITZ_R4.md` | Nielsen types / braid orbits |
| `GENUS_3A4_LOCK.md` | \(3A^4\) genus 0 lock |
| `EXPLICIT_3A4_RESOLVENT.md` / `EXPLICIT_3A4_EQUATION.md` | Explicit model status |
| `AVENUE_RANK_EXECUTE.md` | Avenues 1–7 scorecard |
| `NECESSITY_THEOREM.md` | Open / paused (orthogonal) |
| `RESEARCH_ROADMAP.md` | Priority table |

---

## 8. One-line characterisation

**Geometric multi-\(k\) fusion** = one Nielsen-labelled pure \(A_5\) family whose Hilbert specialisations recover the pure-even / HQCC seed lattice across more than one ratio \(k\) — the arithmetic centre already knows how to live on many \(k\); geometry has not yet named that surface with a braid orbit.

_Adopted as principal open problem after pure-even multi-\(k\) finish; Canonical T3 production lock; Necessity remains open._
