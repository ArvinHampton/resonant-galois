# G1 — parameter-field resolvent

_Elapsed: 22.72s_

**Verdict:** G1 param-field resolvent (22.72s). s=-1: f∈Q(√5)(t)[y], Norm deg 10 over Q(t). R1 seed|norm exact hits=0/15 (common-root signals=0). R2 forward cat=0 BJ=0 even5=0. R3 multi-sheet cat=0 multi_k=False. max_sheets≈2. geometric_catalogue_hit=False.

---

## 0. Goal

Accept the cover fibre as

```text
f ∈ K(s)[x] with [K : Q(s)] > 1
```

form the **norm / multi-sheet resolvent** with coefficients in Q(s)
(or in Q after specialising s), then re-test Hilbert specialisations against
the multi-seed pure-even catalogue (flagship −8/5, classical 4/5, LSW −4, …).

This is the G1 path that does **not** require a single-valued f_s ∈ Q(s)[x].

---

## 1. R0 — exact parameter field at s = −1

| item | value |
|------|-------|
| K | Q(sqrt(5)) |
| [K:Q] | 2 |
| s | -1 |
| params | p2=-1, c=-sqrt(5), σ=0, π=-1/25 |
| f over K(t) | `sqrt(5)*t*y**2/5 - sqrt(5)*t/125 + y**5 - y**3` |
| Norm over Q(t) | `-t**2*y**4 + 2*t**2*y**2/25 - t**2/625 + 5*y**10 - 10*y**8 + 5*y**6` |
| deg_y(Norm) | **10** = 5 × [K:Q] |

So the natural Q-model of this geometric fibre is **degree 10**, not 5.

---

## 2. R1 — reverse: does a catalogue seed divide the s=−1 norm?

For S = y⁵ + α y + β, reduce Norm(y,t) modulo S and require the remainder
(coeffs in Q(t)) to vanish identically in y — solve for t.

| seed | k | exact t hit? | #t | gcd(coeffs) deg | common-root candidates |
|------|---|:------------:|---:|----------------:|-------------------------|
| flagship | -8/5 | False | 0 | 0 | [] |
| flag_145 | -8/5 | False | 0 | 0 | [] |
| flag_320 | -8/5 | False | 0 | 0 | [] |
| classical | 4/5 | False | 0 | 0 | [] |
| s95_76 | 4/5 | False | 0 | 0 | [] |
| s220_176 | 4/5 | False | 0 | 0 | [] |
| lsw_m100 | -4 | False | 0 | 0 | [] |
| lsw_124m | -4 | False | 0 | 0 | [] |
| lsw_m209 | -4 | False | 0 | 0 | [] |
| s180 | -12/5 | False | 0 | 0 | [] |
| s220m | -12/5 | False | 0 | 0 | [] |
| s55_176 | -16/5 | False | 0 | 0 | [] |
| flagship_m | 8/5 | False | 0 | 0 | [] |
| classical_m | -4/5 | False | 0 | 0 | [] |
| lsw4_m100 | 4 | False | 0 | 0 | [] |

**R1 exact hits: 0/15** (seeds with common-root signal on coeff gcd: 0).

---

## 3. R2 — forward specialisation of the s=−1 norm

| quantity | value |
|----------|------:|
| factor degree histogram | {10: 28, 2: 6, 1: 3} |
| even disc deg-5 factors | 0 |
| BJ deg-5 | 0 |
| A5 (among classified even) | 0 |
| **exact catalogue hits** | **0** |
| BJ on a catalogue k-ray (not exact seed) | 0 |

_No catalogue seed among deg-5 factors of the s=−1 norm specialisations._

---

## 4. R4 — sheet counts (degree proxy)

| s | # Newton sheets |
|--:|----------------:|
| -3 | 1 |
| -2 | 1 |
| -1 | 2 |
| -0.5 | 2 |
| 0.5 | 2 |
| 1.5 | 1 |
| 2 | 2 |
| 3 | 2 |
| 4 | 1 |
| 5 | 2 |
| -1.5 | 2 |
| 2.5 | 1 |

- max sheets observed: **2**
- exact [Q(√5):Q] at s=−1: **2**
- note: Sheet count is a geometric upper bound on the number of covers of this normal form over C at fixed s; parameter field degree may be smaller after identifying Galois orbits. Exact [Q(√5):Q]=2 at s=-1.

---

## 5. R3 — multi-sheet product norms at rational s

At each rational s, collect distinct Newton solutions (sheets), form monic
fibres F_i(y;t), and take the product ∏ F_i as a numerical stand-in for
Norm_{K/Q}(f) when the sheets form a Galois orbit. Recognise rational monic
polys, factor, match catalogue.

| quantity | value |
|----------|------:|
| sheet counts | {'-2': 1, '-1': 2, '0.5': 2, '2': 1, '3': 2, '-0.5': 1, '1.5': 1} |
| product recognised over Z | 48 |
| deg-5 factors from products | 0 |
| BJ found | 0 |
| **catalogue hits** | **0** |
| catalogue k | [] |
| multi-k | False |

_No catalogue hits from multi-sheet norms in this scan._

---

## 6. R5 — arithmetic multi-k control

| quantity | value |
|----------|------:|
| pure-even Z samples | 40 |
| disc □ | 16 |
| exact catalogue among samples | 0 |

Arithmetic multi-k remains available; geometric fusion is the open gap.

---

## 7. Multi-k conclusion

| test | result |
|------|--------|
| Parameter-field model f∈K(s)[x] constructed (s=−1) | **True** (K=Q(√5)) |
| Norm over Q of deg 5·[K:Q] | **True** (deg 10) |
| R1 seed divides Norm | **False** (0/15) |
| R2 forward catalogue hit | **False** |
| R3 multi-sheet catalogue hit | **False** |
| Geometric multi-k | **False** |
| Arithmetic multi-k control | **True** |

**Geometric multi-k via parameter-field resolvent of this 3A⁴ normal form: not achieved in this cut.**

### What this cut established

1. **Explicit f ∈ K(t)[y]** at the known geometric fibre s=−1 with [K:Q]=2.
2. **Canonical Norm** is deg 10 over Q(t) — the correct Q-model degree for this fibre.
3. **Reverse division test** of all multi-seed catalogue BJ seeds against that Norm.
4. **Forward factorisation** of Norm specialisations (deg histogram, even/BJ/A5).
5. **Multi-sheet product** proxy for Norm at other rational s + catalogue re-test.

### If still empty — meaning

The pure-even catalogue is not among the Hilbert fibres of the normed 3A⁴
normal-form covers sampled here. Remaining geometric routes:

1. Domain Möbius / different normal form (change {0,1,∞} labels) before norming.
2. Cubic Tschirnhaus on the deg-10 Q-model (not only on deg-5 K-fibres).
3. G2: other genus-0 Nielsen types with possibly rational parameter fields.
4. G3: monodromy identification of the pure-even envelope (arithmetic multi-k
   already has multi-k; give it a Nielsen name).

---

## 8. Non-claims

- Not a proof that no Nielsen realisation of the pure-even lattice exists.
- Negative for this normal form’s parameter-field norm + scan bounds.
- Does not reopen pure-even arithmetic, Canonical T3, or Necessity.

_Generated by `g1_param_field_resolvent.py`._
