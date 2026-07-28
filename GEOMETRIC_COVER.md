# Geometric (branch-cycle) lift — Criterion 1 open problem

## Programme position

Arithmetic Crit 1–3 is theorem-grade on **HQCC seeds** and their homogenisations
(`HQCC_SEED.md`, `THEOREMS.md`). This document is the attack plan for the
**principal open problem**: a cover whose **branch cycle types** are HQCC-native
and whose **geometric monodromy** is alternating.

---

## Target

Construct a cover (or rigid geometric family) whose branch cycle types are
dictated by HQCC data

\[
\{n/3,\ 3n+1,\ 3n-1,\ \ldots\}
\]

and **prove** that its geometric monodromy group is alternating (\(A_5\) or larger \(A_n\)).

Arithmetic specialisations of BJ / homogenised families are already theorem-grade;
the geometric origin is not.

---

## Attack plan

### Step 1 — Rigid cycle-type data labelled by HQCC

Select a rigid tuple of conjugacy classes in \(A_5\) (or \(S_5\)) whose cycle types
can be naturally labelled by the ternary branches of the HQCC map.

Minimal working example for \(A_5\):

- one or more **3-cycles** (the ternary signature),
- a **5-cycle** or double transposition of type \(2{+}2{+}1\),
- chosen so the tuple is **rigid** (Hurwitz dimension zero) and product-one +
  generation conditions hold.

**Status: complete (computation).** See `GEOMETRIC_STEP1.md`.

| Signature | gen \(A_5\) | abs. rigid? | HQCC score | Notes |
|-----------|------------:|:-----------:|-----------:|-------|
| **`(3A,3A,5A)`** / `(3A,3A,5B)` | 60 | **Yes** | **18** | **Preferred:** two ternary classes + one 5-cycle |
| `(2A,3A,5A)` / `(2A,3A,5B)` | 60 | **Yes** | 11 | Classical icosahedral triple |
| `(3A,3A,3A,2A)` etc. (r=4) | >0 | No | high | Families (positive-dim Hurwitz), not a single rigid cover |
| `(3A,3A,3A)`, `(2A,2A,3A)` | 0 | — | — | Do not generate \(A_5\) |

**Step 2 recommendation:** construct the cover for **`(3A,3A,5A)`** first (ternary-heavier rigid triple); keep classical `(2A,3A,5A)` as fallback.

### Step 2 — Geometric realisation — **DONE**

See **`GEOMETRIC_STEP2.md`**.

**Preferred (3A,3A,5A)** over \(\mathbb{Q}\):
\[
\varphi(y)=6y^5-15y^4+10y^3=y^3(6y^2-15y+10)
\]
Branch locus \(\{0,1,\infty\}\), types \((3,1,1),(3,1,1),(5)\), **geometric monodromy \(A_5\)**
(numeric group order 60 + rigidity). HQCC labels: T₃ contraction / T₃ expansion / G₄ period.

**Fallback (3A,2A,5A)** over \(\mathbb{Q}(2^{1/5},3^{1/5})\):
\(\varphi(x)=x^5+a x^4+b x^3\) with radical \(a,b\); types \((3,1,1),(2,2,1),(5)\); monodromy \(A_5\).

### Step 3 — Compatibility with arithmetic

Check specialisations against known HQCC seeds
\(x^5-55x+88\), \(x^5+95x\pm76\), … and homogenisations.
Desirable, not required for the geometric theorem.

### Step 4 — Native labelling

Branch points / local monodromy generators expressed via HQCC operations
\(\{n/3,\,3n\pm1\}\) — not merely “some 3-cycles”. This is HQCC-nativeness
beyond a generic rigid \(A_5\)-cover.

---

## Method constraints

- Classical Hurwitz theory and rigidity methods.
- **Does not** depend on free coefficient search.
- Distinct from HQCC **seed** definition (arithmetic BJ lattice points).

---

## HQCC labelling dictionary (provisional)

| Cycle type in \(A_5\) | HQCC motif |
|----------------------|------------|
| 3-cycle \((3,1,1)\) | Ternary / generations / \(\mathbb{Z}/3\) / branch \(n\mapsto n/3\) or residue of \(\mathrm{Ad}\) |
| Double transp. \((2,2,1)\) | T-complementarity involution / flux pairing |
| 5-cycle \((5)\) | Period / pentagonal sector (use sparingly — \(D_5\) tension) |
| Identity | Unramified |

**Design preference:** maximise **3-cycle** content (ternary signature); treat 5-cycles
as secondary (G4 Heavy \(D_5\) lesson: 5-fold without 3-cycles is the wrong monodromy).

---

## Paths

| Path | Role |
|------|------|
| `GEOMETRIC_COVER.md` | This plan |
| `geometric_step1.py` | Step 1 computation |
| `build/GEOMETRIC_STEP1.md` | Candidate rigid tuples |
| `HQCC_SEED.md` | Arithmetic seed def (separate) |
| `THEOREMS.md` | Programme position |
