# Sharp next options — executed (locked)

**Status (2026-07-24).** Options A–C on the Candidate-C / Tier-1.1 frontier are closed out. Necessity fragment still open. Pure-even multi-k centre unchanged.

---

## A. Binary chooses k, then pure-even under ℋ

Map:
```
k(n) = (−1)^(pop) · (2·odd+1) / (2^(min(v₂,4)) + 1 + (len mod 4))
```
then pure-even (α,β)∈ℚ and monic ℤ model z⁵+Az+B.

| ℋ | Z-coeff rate | disc□ rate | ℋ ⇒ disc□? |
|---|--------------|------------|------------|
| H_all | 1.00 | 1.00 | Yes (on scan) |
| H_core (v₂≤6, len≤20, odd≤64) | 1.00 | 1.00 | Yes (on scan) |
| H_small_v2 / short Collatz / bounded odd | 1.00 | 1.00 | Yes |

**Reading:** Once a pure-even ℤ model exists for k(n), disc□ is classical. So
```
ℋ ⇒ disc□  ⇔  ℋ ⇒ “ℤ pure-even model for k(n)”
```
That is a **composite lemma about F**, not HQCC necessity (pure-even is still in the codomain).

---

## B. F → T(…) only — disc□ → 1?

Honest variants (M-deform, embed shape, full mix, binary homog seed, …):

| Best disc□ rate | Crit-2 signal (>0.5)? |
|-----------------|------------------------|
| **0.00** | **False** |

No template-only functor hits disc□ at high rate. Evenness obstruction on T stands: without pure-even / pre-chosen even seed, **no Criterion 2 signal**.

---

## C. Tier 1.1 return — identical-square subclass of T

| Search | Result |
|--------|--------|
| New multi-param identical-square cuts | **None** |
| Homogenisation no-x² (known) | 1-param, beyond BJ-embed when p≠0 |
| Pure-even envelope | 2-param, inside BJ-embed |
| Bilinear cuts | all fail |

**No new Criterion 2 fragment.**

---

## Locked conclusions

| Option | Outcome |
|--------|--------|
| A. Binary k + pure-even under ℋ | Works as composite; **not necessity** |
| B. F → T only, disc□ → 1 | **Failed** |
| C. New identical-square subclass of T | **None** beyond known families |

| Claim | Status |
|-------|--------|
| Necessity fragment | **Still open** |
| Pure-even multi-k centre | **Finished / unchanged** |
| Design mirror (Candidate C) | Implemented; does not force evenness alone |
| Evenness obstruction on unrestricted T | **Reinforced** |

Script: `tier12_sharp_next.py`

See also: `CANDIDATE_C_FUNCTOR.md`, `NECESSITY_THEOREM.md`, `HQCC_MATRIX_TEMPLATES.md`, `EVENNESS_OBSTRUCTION.md`.
