# Candidate C functor — locked

**Status (2026-07-24).** Design mirror is no longer only an analogy: F1–F3 are total, typed, reproducible maps. Even monodromy on the image is imported (pure-even / catalogue) or not enriched (design-faithful path).

---

## Type signature

```
F: BinaryData ⟶ TernaryOutput
```

**Input BinaryData:** n, v₂(n), odd part, Collatz itinerary {0,1}*  
(0 = even → n/2, 1 = odd → 3n+1).

**Output TernaryOutput:** ternary lattice element + optional T(a,…,f) + optional BJ (α,β).

---

## Three functors (92 seeds, ~96 s)

| Map | Construction | BJ disc□ | BJ A₅ | Template disc□ |
|-----|--------------|----------|-------|----------------|
| **F1** | Itinerary → ternary base-3 + model; M-deform; LSW pure-even k=−4 | 1.00 | 1.00 | 0.00 |
| **F2** | ℓ = 3^(v₂)·(odd mod 3⁵)+61; BJ-embed only (no pure-even) | 0.00 | 0.00 | 0.00 |
| **F3** | Popcount → catalogue k + pure-even envelope | 1.00 | 1.00 | 0.00 |
| Random BJ control | comparable coeff size | ~0 | ~0 | — |

---

## Interpretation

1. **Candidate C is concrete.** F1–F3 are total functions with fixed types — the design mirror is implemented, not only described.

2. **Even monodromy is not forced by binary data alone.**
   - F1 / F3 achieve disc□ and A₅ rate 1.00 by *calling* the pure-even / LSW / catalogue-k machine on the codomain.
   - F2 is the design-faithful path (ternary height → BJ-embed, no pure-even insert) and yields disc□ rate ≈ 0.

3. **No native binary Galois property is preserved.** Binary Collatz data has no Gal monodromy to carry forward. Tests measure *enrichment of the image*, not preservation of a binary Galois invariant.

4. **Not a necessity theorem.** High A₅ rates on F1/F3 are explained by the pure-even codomain construction already finished in the arithmetic centre — not by a hypothesis on n alone.

---

## Relation to the four-face principle

F1–F3 implement

```
dynamics → lattice → matrices → (optional) Galois
```

as functions of binary data. The 3-cycle face is encouraged by ternary digits; the sign face is supplied only when the output is forced into a pure-even BJ ray.

---

## Relation to necessity (Criteria 1–3)

| Route | Effect of F1–F3 |
|-------|-----------------|
| Crit 1 (canonical object) | Functors exist; none is an HQCC-native cover with proved alternating monodromy |
| Crit 2 (structural axioms ⇒ disc□) | F2 shows design-faithful BJ-embed from binary data does **not** give disc□ rate → 1 |
| Crit 3 (sign character) | Not addressed |

---

## Next on this track (optional)

1. Merge F2 height with a pure-even k *determined by binary data* (not hard-coded catalogue k).
2. Prove lemmas of the form: n ∈ ℋ ⇒ disc(F(n)) □ for a stated binary class ℋ.
3. Seek an F that outputs only T(…) (no pure-even insert) with disc□ rate → 1 — that would be Criterion 2 news.

Until (3), Candidate C remains: **implemented design mirror; evenness still classical pure-even on the image.**

Script: `candidate_c_functor.py`

See also: `RESONANT_ALGEBRAIC_CLOSURE.md`, `NECESSITY_THEOREM.md`, `HQCC_MATRIX_TEMPLATES.md`.
