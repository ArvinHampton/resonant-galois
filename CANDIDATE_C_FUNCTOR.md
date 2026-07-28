# Candidate C functor — Tier 1.2

_Elapsed: 96.29s_

**Verdict:** Candidate C functors F1–F3 (96.29s). F1 BJ disc□=1.00 A5=1.00; F2 BJ disc□=0.00 A5=0.00; F3 BJ disc□=1.00 A5=1.00; random control disc□=0.00. Functor explicit; even monodromy on image is imported pure-even or partial — not a necessity theorem.

Turns the mod-2 ↔ mod-3 design mirror into **explicit maps** with typed input/output. Does **not** prove necessity.

---

## Type signature

### Input: `BinaryData`

| field | type | meaning |
|-------|------|---------|
| `n` | \(\mathbb{N}_{>0}\) | binary seed |
| `v2` | \(\mathbb{N}_0\) | \(v_2(n)\) |
| `odd_part` | odd \(\mathbb{N}\) | \(n/2^{v_2}\) |
| `collatz_itinerary` | \(\{0,1\}^*\) | Collatz word: 0=even→\(n/2\), 1=odd→\(3n+1\) |
| `collatz_length` | \(\mathbb{N}_0\) | truncation length |

### Output: `TernaryOutput`

| field | type | meaning |
|-------|------|---------|
| `lattice` | `TernaryLatticeElement` | integer from ternary digits / model mix |
| `template` | `TemplateParams` | \((a,b,c,d,e,f)\) for \(T\) |
| `bj` | `BJPair` | \((\alpha,\beta)\) for \(x^5+\alpha x+\beta\) |
| `functor_id` | string | F1 / F2 / F3 |

$$F:\ \mathbf{BinaryData}\ \longrightarrow\ \mathbf{TernaryOutput}.$$

---

## Three concrete functors

### F1 — itinerary evaluation (design-faithful)

1. Map Collatz word bit \(0\mapsto 0\), \(1\mapsto 1\) or \(2\) (alternating).  
2. Evaluate as base-3 integer; add model generator indexed by \(v_2\).  
3. Template: deform base \(M\) by \(f=\ell \bmod 243\).  
4. BJ: LSW pure-even ray \(k=-4\) with \(m=3^{\min(v_2,8)}+(\ell\bmod 16)\).

- BJ disc□ rate: **1.000** (pure-even by construction)  
- BJ A5 rate: **1.000**  
- Template disc□ rate: **0.000**  
- Template A5 rate: **0.000**

### F2 — valuation transport

1. \(\ell = 3^{v_2}\cdot(\mathrm{odd}\bmod 3^5)+61\).  
2. BJ-embed template: \(d=0\), \(a=-ef\), \(e=3\), \(f=v_2+1\), \(c=61\), \(b\in\{0,80\}\).  
3. BJ pair from embed formulae \(\alpha=-(bf+ce)\), \(\beta=-bc\).

- BJ disc□ rate: **0.000**  
- BJ A5 rate: **0.000**  
- Template disc□ rate: **0.000**  
- Template A5 rate: **0.000**

### F3 — catalogue \(k\) + pure-even envelope

1. Choose multi-seed \(k\) by popcount(itinerary) mod 8.  
2. Pure-even \((\alpha,\beta)\) with integer coefficients when possible.  
3. Lattice label from itinerary base-3.

- BJ disc□ rate: **1.000** (by pure-even design)  
- BJ A5 rate: **1.000**  
- Template disc□ rate: **0.000**

---

## Control and comparison

Random BJ control (\(|lpha|,|eta|\leq H\)): disc□ rate **0.000**, A5 rate **0.000**.

| map | BJ disc□ | BJ A5 | vs control disc□ |
|-----|--------:|------:|:----------------:|
| F1 | 1.000 | 1.000 | enriched (pure-even) |
| F2 | 0.000 | 0.000 | partial |
| F3 | 1.000 | 1.000 | enriched (pure-even) |
| random BJ | 0.000 | 0.000 | baseline |

---

## Sample images

### F1 samples

| \(n\) | lattice | BJ status | disc□ | tmpl status |
|----:|--------:|-----------|:-----:|-------------|
| 1 | 4 | HIT_A5 | True | odd_monodromy |
| 2 | 9 | HIT_A5 | True | odd_monodromy |
| 3 | 22 | HIT_A5 | True | odd_monodromy |
| 4 | 18 | HIT_A5 | True | odd_monodromy |
| 5 | 4 | HIT_A5 | True | odd_monodromy |
| 6 | 66 | HIT_A5 | True | odd_monodromy |
| 7 | 181624 | HIT_A5 | True | odd_monodromy |
| 8 | 61 | HIT_A5 | True | odd_monodromy |
| 9 | 9629662 | HIT_A5 | True | odd_monodromy |
| 10 | 12 | HIT_A5 | True | odd_monodromy |
| 11 | 39631 | HIT_A5 | True | odd_monodromy |
| 12 | 189 | HIT_A5 | True | odd_monodromy |

### F2 samples

| \(n\) | lattice | BJ status | disc□ | tmpl status |
|----:|--------:|-----------|:-----:|-------------|
| 1 | 62 | odd_monodromy | False | odd_monodromy |
| 2 | 64 | odd_monodromy | False | odd_monodromy |
| 3 | 64 | odd_monodromy | False | odd_monodromy |
| 4 | 70 | odd_monodromy | False | odd_monodromy |
| 5 | 66 | odd_monodromy | False | odd_monodromy |
| 6 | 70 | odd_monodromy | False | odd_monodromy |
| 7 | 68 | odd_monodromy | False | odd_monodromy |
| 8 | 88 | odd_monodromy | False | odd_monodromy |
| 9 | 70 | odd_monodromy | False | odd_monodromy |
| 10 | 76 | odd_monodromy | False | odd_monodromy |
| 11 | 72 | odd_monodromy | False | odd_monodromy |
| 12 | 88 | odd_monodromy | False | odd_monodromy |

### F3 samples

| \(n\) | lattice | BJ status | disc□ | tmpl status |
|----:|--------:|-----------|:-----:|-------------|
| 1 | 1 | HIT_A5 | True | odd_monodromy |
| 2 | 0 | HIT_A5 | True | odd_monodromy |
| 3 | 10 | HIT_A5 | True | odd_monodromy |
| 4 | 0 | HIT_A5 | True | odd_monodromy |
| 5 | 1 | HIT_A5 | True | odd_monodromy |
| 6 | 30 | HIT_A5 | True | odd_monodromy |
| 7 | 179425 | HIT_A5 | True | odd_monodromy |
| 8 | 0 | HIT_A5 | True | odd_monodromy |
| 9 | 4844476 | HIT_A5 | True | odd_monodromy |
| 10 | 3 | HIT_A5 | True | odd_monodromy |
| 11 | 19936 | HIT_A5 | True | odd_monodromy |
| 12 | 90 | HIT_A5 | True | odd_monodromy |

---

## Interpretation

Design-faithful (itinerary mirror). BJ image uses pure-even LSW so disc□ rate=1.00 (by construction on LSW ray). Template M-deform disc□ rate=0.00 (not forced even — same obstruction as base M).

Valuation transport. BJ-embed image has disc□ only when (α,β) hit even locus; observed BJ disc□=0.00, A5=0.00. Template often non-BJ or odd.

Strongest even monodromy by design: pure-even catalogue k. BJ disc□ rate=1.00, A5=1.00. Evenness is pure-even theory, not a new force from the binary input.

**Preservation vs enrichment:** Binary Collatz data has no native Gal monodromy to 'preserve'. Tests are enrichment rates on F(image) vs random BJ control. F3 enriches disc□ to ~1 by invoking pure-even; that is not necessity from binary hypotheses alone — pure-even is inserted in the codomain construction.

**Necessity:** None of F1–F3 proves that HQCC axioms force A_n. They make Candidate C checkable: explicit maps exist; even monodromy on the image is either imported (F1/F3 pure-even) or partial (F2).

---

## Relation to the four-face principle

F1–F3 make the **dynamics → lattice → matrices → (optional) Galois** pipeline into total functions on binary data. The 3-cycle face is encouraged by ternary digits and template couplings; the sign face is supplied by pure-even when we choose F1/F3 BJ outputs — i.e. the organising principle is **implemented**, not newly forced.

---

## What to do next (functor track)

**Executed** in `TIER12_SHARP_NEXT.md` / `tier12_sharp_next.py`:

1. Binary \(k(n)\) + pure-even under \(\mathcal{H}\) — disc□ rate 1 when Z model exists (composite lemma).
2. \(F\to T\) only disc□→1 — **not achieved** (best honest rate 0).
3. Tier 1.1 return — no new identical-square subclass.

Further functor work only if a new T-only construction appears.

```bash
python candidate_c_functor.py
```

_Generated by candidate_c_functor.py — Tier 1.2_