# G1 — 3A⁴ triple-root elimination + multi-seed Hilbert test

_Elapsed: 27.73s_

**Verdict:** G1 cut (27.73s). P(q,w) rational pts=4 (non-degenerate=0). Seed reverse=16/16; on triple-root locus=0/16. Known s=-1 param matches=0. Norm fibres: even=0 BJ=0 cat=0. Single-valued f_s in Q(s)[x]=False. Geometric multi-k=False.

---

## 0. Goal

Push the triple-root model for Ni(A₅, C₃⁴) toward a single-valued family fₛ ∈ ℚ(s)[x], then test Hilbert specialisations against the multi-seed pure-even catalogue (flagship −8/5, classical 4/5, LSW −4, and other multi-seed ratios).

Locks: pure-even multi-k finished; Canonical T3 production; Necessity paused. See `GEOMETRIC_MULTI_K_FUSION.md`.

---

## 1. Eliminant chart P(q,w)

Physical component of the triple-root eliminant (from `EXPLICIT_3A4_EQUATION.md`):

```
20*q**3*w**3 - 40*q**3*w**2 + 27*q**3*w - 6*q**3 - 40*q**2*w**3 + 73*q**2*w**2 - 45*q**2*w + 9*q**2 + 27*q*w**3 - 45*q*w**2 + 26*q*w - 5*q - 6*w**3 + 9*w**2 - 5*w + 1
```

- Rational points height ≤ 24: **4**
- Non-degenerate (q,w ∉ {0,1}, q≠w): **0**

| q | w | height | note |
|---|---|-------:|------|
| 1 | 1 | 1 | singular (1,1) |
| 0 | 1/2 | 2 | degenerate |
| 1/2 | 0 | 2 | degenerate |
| 1/2 | 1/2 | 2 | degenerate |

### Quadratic probes

| label | on P? | q | w |
|-------|:-----:|---|---|
| s_m1_physical | True | sqrt(5)/5 | -sqrt(5)/5 |
| s_m1_swap | True | -sqrt(5)/5 | sqrt(5)/5 |
| diag_pos | False | sqrt(5)/5 | sqrt(5)/5 |
| diag_neg | False | -sqrt(5)/5 | -sqrt(5)/5 |

**Obstruction.** Up to height 24 the only rational points are degenerate (zero/pole collisions or the singular point (1,1)). The known physical fibre s=−1 sits on **Q(√5)** points (±1/√5). This blocks a stream of Q-covers from the (q,w)-chart and explains the failure of polyfits for (c,p₂,rᵢ)(s)∈ℚ(s).

---

## 2. Single-valued fₛ ∈ ℚ(s)[x]

- **Achieved:** `False`
- Reason: Eliminant chart P(q,w)=0 has no non-degenerate rational points up to height 24 (only degenerate (0,1/2),(1/2,0),(1,1),(1/2,1/2)). Physical covers live on quadratic points (e.g. s=-1 over Q(sqrt(5))). Normal-form parameters are multi-sheeted over Q(s); closed form f_s in Q(s)[x] not obtained.

H^rd ≅ P¹_s still guarantees some rational moduli coordinate, but **this normal form** is multi-sheeted over ℚ(s). The eliminant chart is not a rational parameter source (g=1 after ordinary blowup; sparse rational points). Exact model at s=−1: monic over ℚ(√5); norm to ℚ(t) is degree 10:

```
5*(y**5 - y**3)**2 - t**2*(y**2 - 1/25)**2 = 0
```

---

## 3. Seed-first attack (reverse BJ → triple-root locus)

Every priority catalogue seed is **compatible** with the N−tD normal form after the y⁴-kill shift (2 reverse solutions each). Imposing the triple-root chart equations (σ,π) = chart(p₂,q) and P(q,w)=0 is the geometric filter.

| seed | k | reverse sols | on triple-root locus? | #locus hits |
|------|---|-------------:|:---------------------:|------------:|
| flagship | -8/5 | 2 | False | 0 |
| flag_145 | -8/5 | 2 | False | 0 |
| flagship_m | 8/5 | 2 | False | 0 |
| classical | 4/5 | 2 | False | 0 |
| s95_76 | 4/5 | 2 | False | 0 |
| classical_m | -4/5 | 2 | False | 0 |
| lsw_m100 | -4 | 2 | False | 0 |
| lsw_124m | -4 | 2 | False | 0 |
| s180 | -12/5 | 2 | False | 0 |
| s55_176 | -16/5 | 2 | False | 0 |
| flag_320 | -8/5 | 2 | False | 0 |
| flag_1145 | -8/5 | 2 | False | 0 |
| s220_176 | 4/5 | 2 | False | 0 |
| s395_316 | 4/5 | 2 | False | 0 |
| s95_m76 | -4/5 | 2 | False | 0 |
| lsw_m209 | -4 | 2 | False | 0 |

**Summary:** reverse 16/16; on locus **0/16**.

_No catalogue seed’s reverse parameters lie on the triple-root locus in this normal form._ Compatibility with N−tD is necessary but not sufficient for a geometric 3A⁴ specialisation.

### Known s=−1 fibre parameter match

Reverse sols with (p₂,σ,π)=(−1, 0, −1/25): **0** (would place the seed on the known geometric fibre for some t ∈ ℚ(√5)).

_None of the catalogue seeds match the exact s=−1 cover parameters._

---

## 4. Hilbert specialisations (s=−1 norm over ℚ)

| quantity | value |
|----------|------:|
| deg-5 factors tested | 0 |
| irreducible | 0 |
| disc square | 0 |
| BJ after y⁴-shift | 0 |
| exact catalogue hits | 0 |

_No catalogue seed recovered from deg-5 factors of the s=−1 norm._

**Structural note (norm factorisation).** For rational \(t\), the cleared norm
\(625\cdot\bigl(5(y^5-y^3)^2-t^2(y^2-1/25)^2\bigr)\) is degree 10. Sample
factor types over \(\mathbb{Q}\): irreducible deg-10 (e.g. \(t=2,3,5\)), or
products of quadratics (e.g. \(t=\pm1\)) — **never a monic deg-5 factor**.
So the geometric fibre over \(\mathbb{Q}(\sqrt5)\) does **not** descend to a
degree-5 model over \(\mathbb{Q}(t)\) by taking norms; a \(\mathbb{Q}\)-model of
this fibre would need a different construction (resolvent / different chart).

---

## 5. Arithmetic multi-k control (must stay green)

| quantity | value |
|----------|------:|
| pure-even Z samples | 80 |
| disc □ | 40 |
| exact catalogue among samples | 2 |

Pure-even multi-k arithmetic continues to supply disc-square BJ fibres on all catalogue ratios. Geometric fusion is the open gap — not arithmetic evenness.

---

## 6. Multi-k conclusion

| test | result |
|------|--------|
| Single-valued f_s ∈ Q(s)[x] | **False** |
| Seeds on triple-root locus | **0/16** |
| Catalogue hit via s=−1 norm | **False** |
| Geometric multi-k | **False** |
| Arithmetic multi-k (control) | **True** |

**Geometric multi-k via this 3A⁴ normal-form cut: not achieved.**

### What this cut established

1. **Sparse rational geometry of P:** no useful Q-points → no Q-cover flood from (q,w).
2. **Normal-form compatibility of all priority seeds** (reverse always solvable).
3. **Triple-root locus is a real filter:** reverse sols generally miss the locus.
4. **s=−1 geometric fibre** does not carry catalogue seeds in its (p₂,σ,π) slot.
5. Arithmetic centre remains healthy (control).

### Next steps (ordered)

1. **Seed-first residual equations:** for reverse (p₂,σ,π), solve the overdetermined
   triple-root system allowing a **mild Tschirnhaus / coordinate change** on the
   domain (not pure y) so that catalogue seeds can sit on a Nielsen fibre.
2. **Parameter-field resolvent:** build f ∈ K(s)[x] with [K:Q(s)]>1 from the
   multi-sheeted normal form; norm to a higher-degree model over Q(s); re-test.
3. **G2:** explicit equations for other g=0 types (2A 3A³, 2A² 3A²).
4. **G3:** monodromy ID of the pure-even envelope.
5. Do **not** reopen pure-even arithmetic, Canonical T3, or Necessity.

---

## 7. Non-claims

- Not a proof that geometric multi-k is impossible.
- Not a change to the pure-even multi-k theorem.
- Negative only for this normal form + height-24 rational search + seed-first filter.

_Generated by `g1_3a4_triple_root.py`._
