# Resonant algebraic closure — detailed exploration

_Turn the phrase into candidate mathematical definitions, see what each closes, and match them to objects already in hand._

---

## Status (one line)

| Claim | Status |
|-------|--------|
| **Candidate A** (lattice + pure-even envelope) | **In use** — effective arithmetic closure of the finished multi-\(k\) theory |
| **Candidate B** (cyclotomic \(\mathcal{R}\)) | **Defined** — coefficient enrichment; no monodromy necessity |
| **Candidate C** (mod-2 ↔ mod-3 design mirror) | **Analogy clear** — functor not written |
| **Candidate D** (homotopy / cobordism 539) | **Narrative** — not outsider-checkable as a closure |
| **Necessity theorem** | Still **open** (`NECESSITY_THEOREM.md`) |

**Bottom line.** In the only sense that currently touches the data and the finished theory, resonant algebraic closure is **Candidate A**. The design mirror of classical binary is **Candidate C**; it explains lattice + 3-cycles fit, not forced alternating monodromy.

**Sharpened picture:** `TERNARY_ORGANIZING_PRINCIPLE.md` — why ternary (not random density) feeds \(A_n\): four faces of the same branching + pure-even sign lock.

---

## 1. Working vocabulary

| Term | Working meaning |
|------|-----------------|
| **Binary side** | Mod-2 (or Collatz-type) branching: contraction on even residue, expansion on odd |
| **Ternary / HQCC side** | Mod-3 branching: \(n/3\), expanding residues via T₃ / balanced-ternary style rules |
| **Resonant data** | Period/lattice integers \(\{3,9,27,61,80,243,539,\ldots\}\), ring \(\mathcal{R}=\mathbb{Z}[2\cos(2\pi/539)]\) |
| **Closure** | An operation that starts from a class of objects and returns a class stable under stated rules |

“Resonant algebraic closure” is **not** a standard term. Below are four precise candidates.

---

## 2. Candidate A — Lattice closure (arithmetic, already effective)

**Objects.** Pairs \((\alpha,\beta)\in\mathbb{Z}^2\).

**Generators.** Resonant / ternary integers
\[
L_0=\{3,9,27,61,80,243,539,\ldots\}.
\]

**Operations.**

1. Short \(\mathbb{Z}\)-combinations (sums/products of bounded length).  
2. Ratio \(k=\beta/\alpha\) when \(\alpha\neq 0\).  
3. Pure-even envelope: for fixed \(k\),
   \[
   \alpha(m)=256m^2-\frac{3125k^4}{256},\qquad \beta=k\cdot\alpha.
   \]
4. Cross-\(k\) paths: linear (or rational) paths in \((m,k)\)-space.

**Closure \(\overline{L}\).** Smallest set containing \(L_0\)-coefficients and stable under (1)–(4).

| Property | Status |
|----------|--------|
| Well-defined | **Yes** |
| Produces pure-even multi-\(k\) | **Yes** (by construction) |
| Uses HQCC lattice | **Yes** |
| Forces \(\mathrm{Gal}=A_5\) | **No** (only evenness; irr + 3-cycles still needed) |
| Mirror of binary | Only by choice of \(L_0\) |

This is the **operational closure the programme already uses**. It is algebraic and resonant; it is **not** a duality with binary Collatz.

**Artefacts:** `ENLARGED_SEED_CATALOGUE.md`, pure-even theory in `RESOLUTION_PATH.md` / `THEOREMS.md`, `PURE_EVEN_SPECIALISATIONS.md`, `lib/lemmas.py`.

---

## 3. Candidate B — Cyclotomic / cosine closure (field-theoretic)

**Base field.** \(\mathbb{Q}\).

**Adjoin.** \(\xi=2\cos(2\pi/539)\), get \(\mathcal{R}=\mathbb{Q}(\xi)\) (degree \(\varphi(539)/2=210\)).

**Close under.** Field operations; optionally polynomials \(f\in\mathcal{R}[x]\) and \(\mathrm{Gal}(f/\mathcal{R})\).

| Property | Status |
|----------|--------|
| Standard algebraic step | **Yes** (finite extension; **not** \(\overline{\mathbb{Q}}\)) |
| Pure-even identity extends to \(\mathcal{R}(m,k)\) | **Yes** (field-agnostic disc identity) |
| Catalogue rational \(k\) lie in \(\mathcal{R}\) | **Yes** (\(\mathbb{Q}\subset\mathcal{R}\)) |
| Cosine \(k=2\cos(2\pi p/n)\) as geometric candidates | Pure-even over \(R_n\); do **not** match multi-seed rational \(k\) |
| Forces alternating monodromy | **No** |

**Limitation (locked).** \(\mathrm{Gal}(R_n/\mathbb{Q})\) does **not** turn multi-\(k\) paths into unions of cosine orbits (orbits finite, paths infinite) — `GAL_RN_MULTI_K_ORBITS.md`. Cosine closure enriches coefficients; it does not reorganise the pure-even envelope into a Galois-orbit picture.

**Artefacts:** `NONCLASSICAL_RESONANT_FIELD.md`, `PURE_EVEN_K_CROSS_RATIO.md`, `R539_A5_HQCC_LATTICE.md`.

---

## 4. Candidate C — Monodromy / branching closure (design mirror)

**Idea.** Binary and ternary are both “modulus-\(p\) branched contractions.” Closure = pass from \(p=2\) data to \(p=3\) data by a fixed rule.

| Binary prototype | Ternary image |
|------------------|---------------|
| Residue mod 2 | Residue mod 3 |
| Contract on \(0\) mod 2 | Contract on \(0\) mod 3 |
| Expand on \(1\) mod 2 | Expand on \(1,2\) mod 3 |
| Height \(\sim\log_2\) | Height \(\sim\log_3\) |
| Cycle-type bias (2-power structure) | Cycle-type bias (3-cycles) |

**Possible formalisation.**

- Replace 2-adic valuation by 3-adic valuation in height functions.  
- Replace \(S_n\) generators built from transpositions by generators built from 3-cycles (as in matrix templates \(T\)).  
- Replace binary partition of \(\mathbb{N}\) by ternary partition (T₃).

| Property | Status |
|----------|--------|
| Clear design analogy | **Yes** |
| Matches Galois pattern (3-cycles ↔ ternary) | **Yes** |
| Explicit functor (binary objects) → (ternary objects) | **Yes** — F1–F3 (`CANDIDATE_C_FUNCTOR.md`) |
| Proves even monodromy | **No** — F1/F3 import pure-even; F2 does not enrich disc□ |

This is the strongest reading of “HQCC mirrors classical binary”: same scheme, modulus 3 instead of 2. Functors are now **constructed and tested**; they do **not** yield necessity.

**Artefacts:** `CANDIDATE_C_FUNCTOR.md`, `candidate_c_functor.py`.

---

## 5. Candidate D — Topological / homotopy closure (HQCC-native narrative)

**Claim in HQCC texts.** 539 arises as cardinality of homotopy classes in a cobordism / flux ensemble (4880 units, 243 towers, three generations).

**Closure idea.** Paths in a physical / resonant configuration space are closed under concatenation up to homotopy; the resonant clock discretises each class to one step.

| Property | Status |
|----------|--------|
| Internal to HQCC narrative | **Yes** |
| Independent mathematical definition of the space and cobordism class | **Not supplied** in the arithmetic repo |
| Connects to pure-even multi-\(k\) identities | **No direct map** |
| Connects to \(\mathrm{Gal}(f/\mathbb{Q})=A_5\) | **No** |

Until the topological space and boundary maps are defined so that an outsider can recompute “539 classes,” this remains **programme language**, not algebraic closure in the sense of Candidates A–C.

---

## 6. How the candidates sit relative to pure-even multi-\(k\)

```
Binary Collatz-type design
         │  (analogy: mod 2 → mod 3)
         ▼
Ternary / HQCC branching  ──►  resonant lattice L₀
         │
         ▼
Candidate A: lattice + pure-even envelope  ──►  multi-k theory (FINISHED)
         │
         ├── Candidate B: base change to ℛ ──►  enrichment, no necessity
         │
         ├── Candidate C: functor (missing) ──►  would make “mirror” theorem-grade
         │
         └── Candidate D: topology (missing) ──►  would make “539 from homotopy” checkable
```

**Evenness itself never needs B–D:** it is BJ geometry on the output of A.

---

## 7. Minimal definitions worth adopting

### Definition (resonant lattice closure)

The **resonant algebraic closure** of \(L_0\) in the pure-even sense is the set of all \((\alpha,\beta)\) obtained from \(L_0\) by

1. integer combinations,  
2. formation of \(k=\beta/\alpha\), and  
3. the pure-even envelope and cross-\(k\) paths.

This object is **closed**, **explicit**, and already generates the finished multi-\(k\) theory.  
**(= Candidate A.)**

### Definition (resonant design mirror)

A **resonant design mirror** of a binary branched map is a mod-3 branched map whose contracting residue is \(0\bmod 3\) and whose expanding residues are the nonzero classes mod 3, together with a coefficient lattice generated by order-3 data.

HQCC/T₃ fits this definition; classical Collatz fits the mod-2 analogue.  
**(= Candidate C as design, not as functor.)**

### Non-definition (for now)

“Topological resonant algebraic closure of binary into HQCC” — not yet a single operation with input/output types and a proved invariance (e.g. even monodromy).  
**(= Candidate D until formalised.)**

---

## 8. What would make closure force necessity

To upgrade the mirror from design analogy to **necessity theorem** (`NECESSITY_THEOREM.md`), one needs at least one of:

1. **Functor** \(F:\{\text{binary branched data}\}\to\{\text{ternary / HQCC data}\}\) such that Gal or monodromy of \(F(X)\) is alternating whenever \(X\) satisfies a stated binary hypothesis; or  
2. **Invariant** \(\chi\) built only from HQCC axioms with \(\mathrm{sgn}\circ\rho=\chi\) (Criterion 3); or  
3. **Axiom list** on resonant matrices implying disc □ + 3-cycles identically (Criterion 2), beyond BJ-embed (`HQCC_MATRIX_TEMPLATES.md`, `T_SUBCLASS_IDENTICAL_SQUARE.md`).

**Candidate A does not do this.** Candidates B–D do not yet either.

---

## 9. Bottom line

| Candidate | Content | Status |
|-----------|---------|--------|
| **A.** Lattice + pure-even envelope | Effective arithmetic closure | **In use; finished multi-\(k\)** |
| **B.** Cyclotomic \(\mathcal{R}\) | Coefficient field enrichment | Defined; no monodromy necessity |
| **C.** Mod-2 ↔ mod-3 design mirror | Parallel branching schemes | Clear analogy; **functor missing** |
| **D.** Homotopy / cobordism 539 | HQCC narrative | Not yet outsider-checkable closure |

**Resonant algebraic closure**, in the only sense that currently touches the data and the finished theory, is **Candidate A**: close the ternary lattice under pure-even operations.

The **mirror of classical binary** is **Candidate C**: same branched-contraction design with modulus 3. That mirror explains why the lattice and the 3-cycles fit together; it does **not**, by itself, prove that HQCC axioms force alternating monodromy.

**Slogan (aligned with necessity lock):**  
generative success (A) ≠ forced alternating monodromy from HQCC axioms.

---

## 10. Next depth (if pursued)

Write an **explicit functor for Candidate C**:

| | |
|--|--|
| **Input** | Binary height / valuation data (e.g. 2-adic height, Collatz-type itinerary, or binary branched cover data) |
| **Output** | Ternary lattice element, HQCC seed coefficients, or matrix template \(T(a,\ldots,f)\) |
| **Test** | Whether even monodromy is preserved or forced on the image |

Until that functor exists with checkable input/output types, Candidate C remains design language.

---

## Document map

| Doc | Role |
|-----|------|
| **`RESONANT_ALGEBRAIC_CLOSURE.md`** | This file |
| `NECESSITY_THEOREM.md` | Forced \(A_n\) from axioms — open |
| `RESOLUTION_PATH.md` / `THEOREMS.md` | Pure-even multi-\(k\) finished centre |
| `NONCLASSICAL_RESONANT_FIELD.md` | Candidate B |
| `HQCC_MATRIX_TEMPLATES.md` | Templates / Crit 2 |
| `GAL_RN_MULTI_K_ORBITS.md` | Cosine orbits ≠ multi-\(k\) paths |

_Generated as programme lock — Resonant algebraic closure candidates A–D._
