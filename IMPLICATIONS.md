# Implications of the resonant Galois findings

This note extracts **what follows** from the finished generative experiment
(36 unique \(A_5\), 4 unique \(A_6\), explicit evenness obstruction, three
resolution criteria). It is not a theorem paper; it is the map from evidence
to claims, non-claims, and the next mathematical moves.

Companion docs: `RESOLUTION.md`, `build/CATALOGUE.md`,
`build/CRITERION{1,2,3}_*.md`, and the archives in `../a5_brute_results/`.

---

## 0. Programme position (locked)

> The **arithmetic** side of the resolution criteria now rests on **proved
> identities** for an infinite family whose coefficients are visibly
> **HQCC-lattice**. The **geometric cover** remains the **principal open problem**.

One-sentence technical summary: the ternary/flux lattice both (i) supplies
**theorem-grade** infinite BJ/homogenised families with model specialisations
\(A_5\), and (ii) still fails to force even monodromy for unrestricted structural
\(T_n\) — while a **proved geometric monodromy cover** from HQCC branch data is
not yet constructed.

---

## 1. Implications for specific polynomials

### 1.1 Eliminated / corrected candidates

| Object | Status | Implication |
|--------|--------|-------------|
| \(g(x)=x^5-80x^4+4880x^3-539x^2+61x-3\) | irr, disc **not** square → odd Gal | Model coefficients alone ⇏ \(A_5\) |
| G4 Heavy \(p_1,p_2\) | irr + disc² but Gal \(=D_5\) | disc² alone ⇏ \(A_5\); need a 3-cycle |
| Base \(M\): \(x^5+3x^3-3x^2-4889\) | Gal \(S_5\) | Companion + ternary + model entries can be **odd** |
| Base \(T_6\): \(x^6+3x^4-3x^2-4889\) | Gal \(S_4\times C_2\) | Same obstruction at degree 6 |

**Hard lesson.** Presence of 3, 61, 80, 539, 4880 in coefficients is **not**
evidence of alternating monodromy. Those constants tag the lattice; they do
not gate parity.

### 1.2 Confirmed realisations (frozen catalogues)

| Group | Count | Representative sources |
|-------|------:|------------------------|
| \(A_5\) | **36** | Lattice search, \(M\)-deformations, model seeds (61/243/539), near-rigid BJ specialisations |
| \(A_6\) | **4** | \(T_6\) enlargement; two sign-twins at disc \(5748065856\), two at \(2300833089\) |
| \(D_5\) | sample **7** | Square-disc even groups **without** 3-cycles |

**Implication.** Alternating groups at model-adjacent degrees are **not rare
accidents** of pure coefficient noise; they recur under controlled
deformations of a fixed structural template. The generator is reproducible
(`build_all.py` + archived scans).

### 1.3 Fingerprint implication

Across ~48 catalogue polys, coefficient tags are dominated by
**ternary / \(3^k\)** (generations, near \(\pm3\), \(3^2\), towers).
Discriminants that are squares factor in varied ways — no single “magic”
square root unifies all hits.

**Implication.** The model leaves a **visible arithmetic fingerprint** on
successful polys, but disc² is **not** explained by a single model integer.
Evenness remains an independent condition.

---

## 2. Implications for Galois theory / methodology

### 2.1 Operational criterion (degree 5) is the correct experimental law

> For irreducible monic \(f\in\mathbb{Z}[x]\) of degree 5:  
> \(\operatorname{disc}(f)\) square **and** some unramified \(p\) of type \((3,1,1)\)  
> \(\Rightarrow \mathrm{Gal}(f/\mathbb{Q})=A_5\).

This is classical group theory (transitive subgroups of \(A_5\): only \(A_5\)
contains a 3-cycle). Empirically:

- All catalogue \(A_5\) hits that were cycle-censused show type \((3,1,1)\).
- All catalogue \(D_5\) hits lack 3-cycles (only 5-cycles and \((2,2,1)\)).

**Implication.** The search pipeline’s design is correct:

```
template → χ → irr → disc² → (3,1,1) census → Gal ID
```

Detection algorithms confirm; they do not invent the group. Structure should
**engineer** the 3-cycle; disc² selects evenness.

### 2.2 Degree 6 is strictly harder

Square disc + irr \(\Rightarrow G\le A_6\) **or** another even transitive
subgroup of \(S_6\). There is **no** single cheap cycle type that forces
\(A_6\) the way \((3,1,1)\) forces \(A_5\).

Empirically: **340** square-disc \(T_6\) survivors vs only **4** true \(A_6\).
Even groups include \(A_4\), \(S_4\)-type embeddings, etc.

**Implication.** Inductive enlargement \(T_5\to T_6\to T_n\) works as a
**generator**, but identification cost and false-positive even groups rise
fast. For \(n\gtrsim 12\) and for sporadics, pure lattice search is the wrong
tool; geometric monodromy / rigidity is required (see inductive roadmap).

### 2.3 Structural ternary injects 3-cycles *a priori* — at a price

Constructions S1–S5 (equivariant matrices, HQCC blocks, \(\omega\)-norms,
orbit partitions, BJ-shaped deformations) make \(\mathbb{Z}/3\) structure
part of the definition. When they hit irr + disc² + \((3,1,1)\), Gal is
\(A_5\) by the operational theorem.

**Tension.** Strong equivariance often forces **reducibility** (block
factorisation). The productive middle ground is:

- residual 3-structure or deformations that **break** reducibility while
  keeping a 3-cycle in reduction (S4/S5 style; \(M\)-deformations);
- not pure block-diagonal \(\mathrm{Ad}_{\mathrm{SO}(3)}\) char polys alone.

**Implication.** “Force 3-cycles from ternary data” is **partially solved
constructively**; “force them **and** keep irr **and** disc²” remains open as
a theorem.

---

## 3. Implications for the 9 Maths / HQCC model

### 3.1 What the model is good for (supported)

1. **Lattice design.** Constants \(\{3,9,18,61,80,243,520,539,4880\}\) define
   a finite search space that is **rich enough** to hit \(A_5\) and \(A_6\)
   repeatedly.
2. **Template geometry.** Companion chains + last-block couplings
   (ternary/flux slots) are a workable encoding of “structural matrix”
   for monodromy experiments.
3. **Narrative alignment.** Ternary (generations / qutrit / \(\mathbb{Z}/3\))
   correctly points at **3-cycles** as the group-theoretic feature that
   separates \(A_5\) from \(D_5\). The G4 Heavy \(D_5\) episode is exactly
   the pentagonal vs ternary distinction the model cares about.
4. **Inductive hint.** Same template class produces \(A_5\) and \(A_6\);
   the method is not degree-5-only.

### 3.2 What the model does **not** yet buy (unsupported)

1. **Evenness.** Ternary weight, \(\det\pm1\), trace-zero, and companion-chain
   axioms do **not** force disc². Measured rates among irreducibles in the
   tiny ternary pool: disc² \(\sim 0.7\%\); ternary weight \(\ge 3\) gave rate
   **0** in the Criterion-3 sample.
2. **Canonical monodromy object.** Möbius/HQCC blocks and cubic resultants
   recover classical small groups (\(A_3\), one near-rigid \(A_5\)) but do not
   yet yield a **proved** geometric monodromy \(A_n\) for a single HQCC cover.
3. **Sign character = ternary invariant.** Criterion 3 finds **no** simple
   ternary statistic implying \(\operatorname{sgn}\circ\rho=1\).

### 3.3 Conceptual reframing

| Old informal claim | Refined claim supported by data |
|--------------------|----------------------------------|
| “Model numbers make \(A_5\)” | Model lattice **samples** many \(A_5/A_6\); not a forcing law |
| “Ternary forces alternating” | Ternary structure **targets 3-cycles**; parity is separate |
| “Base structural \(M\) is the monodromy” | Base \(M\) and \(T_6\) are **odd**; monodromy lives in a **subclass** or a **different object** |
| “Experiment = resolution” | Experiment = **evidence + obstruction**; resolution = Criteria 1–3 theorems |

**Implication for the 9 Maths.** The productive interface with Galois theory is
not “every resonant matrix has Gal \(A_n\)” (false). It is:

> There should exist a **canonical** resonant object (cover / representation /
> moduli specialisation) for which 3-cycles and even monodromy are
> **theorems**, with the experimental catalogues as Hilbert specialisations
> and regression tests.

That is exactly Criteria 1–3.

---

## 4. Implications for the three resolution criteria

### Criterion 1 — Canonical HQCC object

**Status:** scaffolds exist; no monodromy proof.

**Implication of findings:**

- Free matrix deformation is **too loose** to be the theorem object
  (base \(M\) is \(S_5\)).
- Rigid / near-rigid families and resultants from cubic (ternary) data are
  closer to a **geometric** statement.
- A theorem-shaped target:

  > ∃ finite cover \(X\to Y\) (or \(f_t\in\mathbb{Q}(t)[x]\)) built only from
  > HQCC branch data with geometric monodromy \(A_n\); model integers arise
  > by specialisation (Hilbert).

**Best next attack:** fix one rigid \(A_5\) branch-cycle type compatible with
HQCC ternary branching; prove monodromy; specialise into the model lattice
and match catalogue hits.

### Criterion 2 — Axioms ⇒ disc² + 3-cycles

**Status:** obstruction explicit; no forcing axiom list.

**Implication of findings:**

- Weak axioms A1–A3 (companion, ternary entry, model entries) are **refuted**
  as a disc² theorem by base \(M\) and base \(T_6\).
- Stronger candidates (B1–B5 in `CRITERION2_AXIOMS.md`) point away from
  “integer matrix with 3’s” toward **algebraic groups / volume forms /
  norms / rigid branch types**.
- Shrinking the class until disc² becomes **provable** — even if much smaller
  than experimental \(T_n\) — is higher value than more brute force.

**Best next attack:** prove disc² for a **thin** subclass (e.g. self-adjoint
w.r.t. a fixed model quadratic form, or char polys that are norms from a
controlled cubic étale algebra), then check whether catalogue hits lie in
that subclass.

### Criterion 3 — Sign character / ternary invariant

**Status:** correlations quantified; no invariant yet.

**Implication of findings:**

- Sign/evenness **is** the disc² gate; inventing a fancy name does not remove
  the obstruction.
- \(\det(M)=\pm1\) is **not** a useful proxy for disc² in samples.
- Candidate mechanisms (volume form on root space, T-complementarity
  involution, quadratic character of the flux lattice, \(\omega\)-norms)
  should be tested as **definitions of subclasses**, not as tags on the full
  lattice.

**Best next attack:** pick one candidate invariant, define the subclass of
matrices/polys that satisfy it, recompute disc² rate; aim for rate \(1\)
(or a proof), not another \(0.5\%\).

---

## 5. Implications for research priority (what to do / not do)

### Do next (theorem track)

| Priority | Action | Why |
|----------|--------|-----|
| P1 | Shrink class → **prove** disc² (Crit 2) | Evenness is the dominant obstruction |
| P2 | One rigid HQCC-linked \(A_5\) family with proved monodromy (Crit 1) | Turns generator into geometry |
| P3 | Define a sign invariant with rate → 1 (Crit 3) | Unifies Crit 2/3 language |
| P4 | Keep catalogues as **regression** | Any proposed theorem must recover known hits |

### Do not invest heavily (diminishing returns)

| Avoid | Why |
|-------|-----|
| Larger pure lattice scans at deg 5–6 | Generator already works; more hits do not close the gap |
| Claiming \(A_n\) from model coefficients alone | Refuted by \(g\), base \(M\), \(T_6\) |
| Monster / sporadic coefficient search | Wrong scale; needs representation / geometric monodromy |
| Treating disc² alone as \(A_5\) | \(D_5\) counterexamples |

### Optional experimental side-quests (only if they feed a theorem)

- Parametric families through best 61/243/539 seeds (Hilbert-style families).
- Controlled off-diagonal mixing in S1 to reduce reducibility while keeping
  residual 3-cycles.
- \(A_7/A_8\) probes only as **ceiling tests**, not as mainline.

---

## 6. Claims vs non-claims (locked)

### Claimed

1. The ternary/flux **method systematically finds** \(A_5\) and \(A_6\).
2. The **obstruction to a theorem is measurable** (base \(M\), base \(T_6\),
   subclass rates, sign correlations).
3. Degree-5 identification law **irr + disc² + \((3,1,1)\) ⇒ \(A_5\)** is the
   correct operational theorem for this programme.
4. The package is a **complete handoff surface** for audit or continuation
   (`resonant_galois/`, catalogues, criteria builds).

### Not claimed

1. Structural axioms alone force disc².
2. Ternary weight or \(\det\pm1\) forces even monodromy.
3. Criterion 1 already supplies a proved HQCC monodromy \(A_n\).
4. Inductive lattice search will reach sporadics.

---

## 7. Theorem promotion — results of the full attack

Executed via `theorem_attack.py` (all three criteria). Full write-up:
`THEOREM_ATTACK.md`.

### Proved / lemma-grade

1. **BJ disc formula.** \(\operatorname{disc}(x^5+ax+b)=256a^5+3125b^4\)
   (symbolic identity verified).
2. **Evenness on the BJ class.** For irreducible \(x^5+ax+b\in\mathbb{Z}[x]\),
   Gal is even iff that integer is a square; with type \((3,1,1)\), Gal \(=A_5\).
3. **Homogenised A5 family (theorem).**
   \[
   f_t = x^5 + 20 t^4 x + 16 t^5
   \]
   satisfies
   \(\operatorname{disc}(f_t)=t^{20}\cdot\operatorname{disc}(x^5+20x+16)\),
   which is a square for all \(t\in\mathbb{Z}\setminus\{0\}\).
   Hence \(\mathrm{Gal}(f_t/\mathbb{Q})\le A_5\) whenever \(f_t\) is
   irreducible. Empirically **all 33** tested specialisations (including
   model \(t\in\{3,9,61,80,243,539\}\)) are irreducible with Gal \(=A_5\).
4. **Operational A5 criterion** (group theory) remains the filter that
   upgrades even+3-cycle to \(A_5\).

### Still open

| Gap | Why it remains |
|-----|----------------|
| Full \(T_5/T_6\) axioms ⇒ disc² | Rates under det±1 / self-adjoint / ternary still \(\sim 1\%\) |
| HQCC-native *geometric cover* | Resultant cubics prefer \(A_3\); BJ seeds are arithmetic, not branch-cycle covers |
| Sign invariant on unrestricted structural \(M\) | Best empirical rates ≪ 1; only thin classes solved |

### HQCC-native update (executed)

Strict HQCC lattice search found **12 non-classical** BJ seeds with Gal \(A_5\)
(beyond classical \(x^5\pm20x\pm16\)). Flagship examples:

| Seed | Note |
|------|------|
| \(x^5-55x\pm88\) | \(88=61+27\) (punctures \(+\) \(3^3\)) |
| \(x^5+95x\pm532\) | \(532=539-7\) (period-adjacent) |
| \(x^5+95x\pm76\), \(x^5-100x\pm400\), \(x^5+124x\pm496\) | strict lattice |

**Theorem (same proof as classical homogenisation):** for seed \((\alpha,\beta)\) with
square disc,
\(f_t=x^5+\alpha t^4 x+\beta t^5\) has square disc for all \(t\neq0\).
Primary family \(x^5-55 t^4 x+88 t^5\): **proved even**, **9/9** model \(t\) specialisations Gal \(A_5\).

Docs: `HQCC_NATIVE.md`, `HQCC_STRICT.md`. Run: `python hqcc_native.py`.

### Recommended next cut

1. Canonicalise the shortest seed in \(\{3,61,80,243,539,4880\}\)-arithmetic.
2. Geometric cover from HQCC branch cycle types (still open).
3. Gröbner form of `disc(χ_T5)` on HQCC template slots.

---

## 8. Paths

| Path | Role |
|------|------|
| `resonant_galois/IMPLICATIONS.md` | This document |
| `resonant_galois/RESOLUTION.md` | Resolution criteria + build status |
| `resonant_galois/HQCC_NATIVE.md` | HQCC-native attack report |
| `resonant_galois/HQCC_STRICT.md` | Strict lattice provenance |
| `resonant_galois/build/` | Regenerable criterion reports + catalogue |
| `a5_brute_results/` | Heavy scan archives (DEFORM_M, A6_T6, fingerprints, …) |
| `python build_all.py` | Full rebuild |
| `python hqcc_native.py` | HQCC-native pass |

_Implications updated after HQCC-native theorem class (strict A5 seeds + homogenised families)._
