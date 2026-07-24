# Tier 1.1 deepen — locked

**Status (2026-07-24)**

No Criterion-2 fragment. Nothing is simultaneously:
- identically disc-square,
- beyond classical pure-even,
- 3-cycle-forced by structure,
- and HQCC-axiom-named.

**Priority shift:** further 1.1 only with a new algebraic idea. Default \(\to\) **Tier 2** (paper / catalogue invariants / verification).

---

## Largest natural subclass with disc identically square

| Family | Beyond BJ-embed? | Free params | Identically \(\square\)? | HQCC-native? |
|--------|------------------|-------------|------------------------|--------------|
| Pure-even envelope | No | 2 (m,k) | Yes | No |
| Pure-even k-slice | No | 1 | Yes | No |
| Homog fixed even seed \(x^5+px^3+rx+s\) (a=e=0, p\(\neq\)0) | Yes | 1 | Yes in t | No |
| Weighted scale of fixed seed | Yes | 1 | Yes | No |
| Genuine poly 1-param family of seeds \(\to\) 2-param after homog | — | — | 0 found | — |

**Ceiling:** overall dim 2 = pure-even envelope; beyond BJ-embed dim 1 only.

---

## Scan summary

- Disc(\(\chi_T\)): total degree 12; structural cuts odd (not identically square).
- On a=e=0: factor \(b^2\times\)(deg-8 odd part) — square only under further relations.
- Genuine 2-param beyond BJ: **0**.
- 3-cycles on homog seeds: HIT_A5 by Hilbert inheritance from seed, not forced by T.
- HQCC naming: pure-even = classical BJ ansatz; homog no-x\(^2\) contradicts base M (a=3).

---

## Necessity verdict

| Criterion 2 demand | Met? |
|--------------------|------|
| Identically disc \(\square\) | Only on classical pure-even / homog rays |
| Beyond BJ-embed at dim \(\ge 2\) | No |
| 3-cycles forced by T structure | No |
| HQCC-axiom named | No |

**Explicit failed probe** for Criterion 2 on the current template class.

Script: `tier11_deepen.py`
