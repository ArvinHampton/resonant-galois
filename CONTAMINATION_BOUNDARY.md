# Contamination boundary — locked maths package

**Rule:** Lattice integers may be *motivated* by the model; proofs may only use their **arithmetic properties**.

Necessity remains **paused**. This file protects the citable / review-grade layer from non-mathematical inputs.

---

## Contamination risks (from wider project history)

| Risk | Status in locked maths package | How to stay clean |
|------|--------------------------------|-------------------|
| Using **G₄ / 539.9 s** as a *proved physical period* inside Gal proofs | **Avoided** if \(539\) is only a **lattice integer** | Write \(539\in\mathbb{Z}\); never “period 539.9 s implies Gal \(=A_5\)” |
| Citing unverifiable **GW / Belle II “hits”** as support for \(A_5\) families | **Out of scope** — must stay out | No experimental citations in pure-even / Mestre / B-avatar theorems |
| Treating HQCC **539-step dynamical claims** as input to disc identities | **Not required** — good | Disc identities live in \(\mathbb{Q}(m,k)\), \(\mathbb{Q}[t]\), \(\mathbb{Q}[A]\) with no dynamical hypothesis |
| Blurring **“motivation”** with **“hypothesis of a theorem”** | **Watch language** in any external write-up | Use “seed chosen from lattice” not “HQCC forces evenness” |

---

## Allowed vs forbidden language

| Allowed (motivation / data) | Forbidden (as theorem input) |
|-----------------------------|------------------------------|
| “Flagship seed uses \(88=61+3^3\) from the resonant lattice” | “Because 539 is the G₄ period, disc is square” |
| “Specialise \(t\in\{3,9,61,80,539,\ldots\}\)” | “Physical 539.9 s enters the resultant” |
| “HQCC motivates searching ternary-rich integers” | “HQCC axioms imply \(\mathrm{Gal}=A_5\)” (necessity; paused) |
| “Model integers \(\{3,9,27,61,80,243,539,4880\}\)” as a set \(L_0\subset\mathbb{Z}\) | “Belle II / GW confirmation of the family” |
| Four-face organising principle as **structural reading** | Organising principle as a **proved** monodromy theorem |

---

## What proofs may use

For objects in the **locked maths package** (below), a correct proof may use only:

1. Polynomial algebra over \(\mathbb{Q}\) (resultants, discriminants, factorisation).  
2. Group theory of \(S_5/A_5\) (operational criterion: irr + disc □ + type \((3,1,1)\)).  
3. Integer / rational specialisations and Hilbert-type sampling.  
4. Named parameters that are **elements of \(\mathbb{Z}\) or \(\mathbb{Q}\)** (including \(539\), \(61\), \(80\), …) **without** physical units or dynamical laws.

A proof **must not** use:

- SI periods, detector claims, entanglement narratives, or 539-step maps as *hypotheses*.  
- “9 Maths” labels as *logical premises* (labels on branch points in geometric docs are dictionary-level only).

---

## Locked maths package (contamination-clean by design)

| Doc | Role | Physical/GW/Belle? |
|-----|------|:------------------:|
| `PURE_EVEN_MULTI_K.md` | Theorem-grade pure-even centre | **No** |
| `MESTRE_FLAGSHIP_PT.md` | Flagship Mestre lift (review #1) | **No** — \(539\) only if as lattice \(t\) |
| `B_EMBED_LATTICE.md` | B-family + embed (review #2) | **No** |
| `EVENNESS_AVATAR.md` | Matrix packaging of identities | **No** |
| `REVIEW_PACKAGE.md` | Review order + pass criteria | **No** |
| `MATH_INTEGRITY_REVIEW.md` | Machine-checked identities | **No** |
| `lib/lemmas.py` | BJ / pure-even algebra | **No** |
| `review_flagship_b.py` | Auto-check | **No** |

Supporting lattice research (`L0_*`, `TERNARY_LATTICE_DIRECTIONS.md`) stays generative: monoid, invariants, Mestre orbit — still **arithmetic only**.

---

## Docs that may contain model labels (not theorem premises)

These files may mention G₄, 539.9, T₃, flux, etc. as **labels, motivation, or open fusion language**. They are **not** inputs to disc identities in the locked package:

| Area | Examples |
|------|----------|
| Geometric / Hurwitz scaffolding | `GEOMETRIC_STEP1.md`, `GEOMETRIC_STEP2.md`, `GEOMETRIC_COVER.md` |
| Criteria / fusion narrative | `NECESSITY_THEOREM.md`, `FUSION_*.md`, `HQCC_SEED.md` (seed *definition* may mention ξ; not used in pure-even proof) |
| Historical catalogues | Older build logs, avenue notes |

**External write-up rule:** If a sentence could be read as “physics ⇒ Galois,” rewrite it as “lattice integer chosen for model reasons; the proof uses only …”

---

## Integer 539 specifically

| Use | OK? |
|-----|:---:|
| \(t=539\) or \(A=539\) as a specialisation parameter in \(\mathbb{Z}\) | **Yes** |
| Generator of monoid \(M_0=\langle 3,61,80,243,539\rangle\) | **Yes** |
| “G₄ = 539.9 seconds is a proved period used in the disc proof” | **No** |
| \(\xi=2\cos(2\pi/539.9)\) as coefficient field for the flagship theorem | **No** (non-classical track only; separate) |

---

## Checklist for any new theorem write-up

1. Can every hypothesis be stated in pure algebra / number theory?  
2. Does any step mention detectors, SI units, or dynamical period laws? → **remove or demote to motivation footnote**.  
3. Is HQCC used only to *select* seeds/parameters, not to *prove* disc □?  
4. Is necessity still labelled open/paused if not proved?

---

## One-line rule (repeat)

> **Lattice integers may be motivated by the model; proofs may only use their arithmetic properties.**

_Locked boundary for the Resonant Galois maths package. Physical / HQCC dynamical narratives are a separate track._
