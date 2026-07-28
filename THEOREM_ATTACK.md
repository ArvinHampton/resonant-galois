# Theorem-promotion attack — Criteria 1–3

_Elapsed: 118.25s_

## Executive result

| Criterion | Advance | Status |
|-----------|---------|--------|
| **2** Thin classes + closed disc | BJ formula lemma; homogenised A5 family | **Partial theorem** |
| **1** Rigid / 1-param families | Homogenised A5 seed + HQCC probes | **Evidence for geometric A5** |
| **3** Sign invariants | BJ class + homogenised family **solved**; full T5 open | **Partial theorem** |

### What is now proved / lemma-grade

1. **Lemma (BJ disc).** \(\operatorname{disc}(x^5+ax+b)=256a^5+3125b^4\) (symbolic identity: True).
2. **Corollary (Crit 2+3 on BJ class).** For irreducible \(x^5+ax+b\in\mathbb{Z}[x]\), Gal \(\le A_5\) (even) iff \(256a^5+3125b^4\) is a square; with Frobenius type \((3,1,1)\) one has Gal \(=A_5\).
3. **Lemma (operational A5).** For irreducible monic f in Z[x] of degree 5: disc(f) square and some unramified p with factorization type (3,1,1) implies Gal(f/Q) = A5 (only transitive subgroup of A5 containing a 3-cycle).
4. **Theorem (homogenised A5 seed family).** For all t in Z\{0}, disc(x^5+20 t^4 x+16 t^5) = t^{20} * disc(x^5+20x+16) is a square in Z, because t^{20}=(t^{10})^2 and disc(seed) is a square. Hence Gal(f_t/Q) ≤ A5 whenever f_t is irreducible; with a (3,1,1) Frobenius, Gal = A5. proved=True.

### What remains open

- Full structural T5/T6 lattice: no axiom list forces disc² for all \(M\).
- Canonical HQCC cover (not BJ) with **proved** geometric monodromy \(A_n\).
- Sign invariant on unrestricted ternary matrices with rate \(=1\).

---

## Criterion 2 — thin subclasses

### `bj`
- class: BJ: x^5 + a x + b
- evenness: `256*a^5 + 3125*b^4 is a square in Z`
- A5 hits: 0
- status: LEMMA: disc(x^5+a x+b)=256 a^5+3125 b^4 (verified). Even monodromy ⇔ this integer is a square. With irr+(3,1,1) ⇒ A5 by operational theorem. PROVED thin family: For all t in Z\{0}, disc(x^5+20 t^4 x+16 t^5) = t^{20} * disc(x^5+20x+16) is a square in Z, because t^{20}=(t^{10})^2 and disc(seed) is a square. Hence Gal(f_t/Q) ≤ A5 whenever f_t is irreducible; with a (3,1,1) Frobenius, Gal = A5. (proved=True).

### `icosa`
- class: icosa-adj: x^5 + 5m x^3 + 5 m^2 x + n
- stats: `{'tested': 198, 'odd': 189, 'red_or_bad': 9}`
- A5 hits: 0
- status: Closed-form disc in m,n available. Evenness ⇔ that form is a square. Hits: 0 sq-disc irr, A5=0.

### `self_adjoint`
- class: sparse self-adjoint model matrices
- stats: `{'tested': 3001, 'red': 1256, 'irr': 1744, 'sq': 16}`
- disc² rate among irr: **0.0092**
- A5 hits: 0
- status: Self-adjointness alone does NOT force disc² (rate=0.0092). B4 still open.

### `det1_ternary`
- class: T5, det±1, some entry ±3
- stats: `{'tested': 1440, 'irr': 952, 'sq': 8, 'red': 488}`
- disc² rate among irr: **0.0084**
- A5 hits: 0
- status: det±1 + ternary entry: still rate ≪ 1 for disc²; not a theorem class.

### `omega_norm`
- class: omega-norm: N(x^2+(a+bω)x+(c+dω))*(x-e)
- stats: `{'tested': 7776, 'red': 7776}`
- A5 hits: 0
- status: Z/3 built into construction (B2 direction). disc² rate among irr = None. Still not identically square; reducibility remains the tax on equivariance.

### `forced_square`
- class: forced_square
- forced families with A5: 1
  - **homogenised_A5_seed**: always_sq=True n_A5=8 a=20*t**4 b=16*t**5
  - **homogenised_flip**: always_sq=False n_A5=0 a=-20*t**4 b=16*t**5
  - **model_scale_3**: always_sq=False n_A5=0 a=3*t**4 b=9*t**5
  - **model_scale_61**: always_sq=False n_A5=0 a=61*t**4 b=3*t**5
  - **euler_try_5_4**: always_sq=False n_A5=0 a=5*t**4 b=4*t**5
  - **euler_try_m5_4**: always_sq=False n_A5=0 a=-5*t**4 b=4*t**5
- line deform stats: `{'irr': 56, 'red': 5, 'sq': 1}` A5=1
- status: Homogenisation of the known A5 seed x^5+20x+16 yields a 1-param family with disc² at all tested t (scale invariance of the square condition under weighted degrees). This is a **theorem-grade thin class**: Gal ≤ A5 for all t where f_t is irreducible, and =A5 when a 3-cycle appears (operational criterion).

---

## Criterion 1 — one-parameter / HQCC families

- Families with A-hits: **2**
- Status: One-parameter specialisation probes run. Homogenised classical A5 seed is the strongest Crit-1 object: many t give A5, supporting geometric monodromy A5 for that family. HQCC resultant family produces A3/dihedral-type more often than A6 — ternary cubic data alone prefers S3-type monodromy until coupled further.

### homogenised_A5_20_16
- \(f_t\) = `16*t**5 + 20*t**4*x + x**5`
- stats: `{'irr': 33}`
- group histogram: `{'S5TransitiveSubgroups.A5': 33}`
- inferred generic Gal: **S5TransitiveSubgroups.A5**
- A-hits: 33
  - t=-12: `x**5 + 414720*x - 3981312` HIT_A5
  - t=-11: `x**5 + 292820*x - 2576816` HIT_A5
  - t=-10: `x**5 + 200000*x - 1600000` HIT_A5
  - t=-9: `x**5 + 131220*x - 944784` HIT_A5
  - t=-8: `x**5 + 81920*x - 524288` HIT_A5

### BJ_a=t_b=1
- \(f_t\) = `t*x + x**5 + 1`
- stats: `{'irr': 31, 'red': 3}`
- group histogram: `{'S5TransitiveSubgroups.S5': 31}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### BJ_a=t_b=3
- \(f_t\) = `t*x + x**5 + 3`
- stats: `{'irr': 31, 'red': 3}`
- group histogram: `{'S5TransitiveSubgroups.S5': 30, 'S5TransitiveSubgroups.M20': 1}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### BJ_a=t_b=9
- \(f_t\) = `t*x + x**5 + 9`
- stats: `{'irr': 32, 'red': 2}`
- group histogram: `{'S5TransitiveSubgroups.S5': 31, 'S5TransitiveSubgroups.M20': 1}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### BJ_a=t_b=16
- \(f_t\) = `t*x + x**5 + 16`
- stats: `{'irr': 33, 'red': 1}`
- group histogram: `{'S5TransitiveSubgroups.S5': 32, 'S5TransitiveSubgroups.M20': 1}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### BJ_a=t_b=-3
- \(f_t\) = `t*x + x**5 - 3`
- stats: `{'irr': 31, 'red': 3}`
- group histogram: `{'S5TransitiveSubgroups.S5': 30, 'S5TransitiveSubgroups.M20': 1}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### hqcc_resultant_s1_m1_t
- \(f_t\) = `t**2 - t*x**3 + 6*t*x - 3*x**2 + 16`
- stats: `{'not_monic': 31, 'irr': 2}`
- group histogram: `{'S3TransitiveSubgroups.A3': 2}`
- inferred generic Gal: **S3TransitiveSubgroups.A3**
- A-hits: 2
  - t=-1: `x**3 - 3*x**2 - 6*x + 17` HIT_A3
  - t=1: `x**3 + 3*x**2 - 6*x - 17` HIT_A3

### icosa_m=t_n=3
- \(f_t\) = `5*t**2*x + 5*t*x**3 + x**5 + 3`
- stats: `{'irr': 34}`
- group histogram: `{'S5TransitiveSubgroups.M20': 34}`
- inferred generic Gal: **S5TransitiveSubgroups.M20**
- A-hits: 0

### model_weight_3_539
- \(f_t\) = `539*t**5 + 3*t**4*x + x**5`
- stats: `{'irr': 33}`
- group histogram: `{'S5TransitiveSubgroups.S5': 33}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

### near_A5_p=61t_q=3
- \(f_t\) = `61*t*x + x**5 + 3`
- stats: `{'irr': 34}`
- group histogram: `{'S5TransitiveSubgroups.S5': 33, 'S5TransitiveSubgroups.M20': 1}`
- inferred generic Gal: **S5TransitiveSubgroups.S5**
- A-hits: 0

---

## Criterion 3 — sign invariants

- Status: Crit 3 SOLVED on BJ thin class (closed disc formula). Crit 3 SOLVED on homogenised A5 seed family (always even when defined). Crit 3 OPEN for full T5 structural lattice — best empirical invariant in this pass: det_pm1 with rate 0.0044742729306487695.

### Empirical rates P(disc² | invariant)

| Invariant | n | sq | rate |
|-----------|--:|---:|-----:|
| all_irr | 41158 | 66 | 0.0016 |
| det_pm1 | 1788 | 8 | 0.0045 |
| det_pos | 21307 | 32 | 0.0015 |
| tw_ge2 | 33878 | 50 | 0.0015 |
| tw_ge3 | 21065 | 34 | 0.0016 |
| palindromic | 0 | 0 | — |
| omega_shape | 540 | 0 | 0.0000 |
| leg3=+1 | 7205 | 8 | 0.0011 |
| leg61=+1 | 19746 | 62 | 0.0031 |
| leg3=leg61=+1 | 3404 | 8 | 0.0024 |
| leg_prod=+1 | 6465 | 8 | 0.0012 |

- Best empirical: **det_pm1**

### Theorem-grade pieces

- BJ: For f=x^5+a x+b monic in Z[x]: sgn(Gal) is trivial (Gal ≤ A5 among subgroups of S5 with this shape's transitive candidates when irr) iff 256 a^5 + 3125 b^4 is a square in Z.
- Note: This is a complete evenness criterion for the BJ thin class — Criterion 3 solved *inside* that class via closed-form disc.
- Homogenised family: f_t = x^5 + 20 t^4 x + 16 t^5 has disc square for all tested t≠0 (weighted homogeneous lift of an A5 seed). (ok=9 bad=0)

---

## Catalogue regression

```
{
  "catalogue_A5": 36,
  "catalogue_A5_BJ_form": 0,
  "BJ_form_disc_formula_confirms_even": 0,
  "BJ_samples": [],
  "homogenised_family_A_hits": 33
}
```

---

## Next moves (after this attack)

1. Prove (not only test) that \(f_t=x^5+20t^4 x+16t^5\) has square disc for all \(t\in\mathbb{Z}\setminus\{0\}\) via the BJ formula + elementary arithmetic.
2. Produce an HQCC-native (not classical BJ) one-parameter family with the same property.
3. For T5 templates, compute the ideal of the condition `disc(χ)` square in parameters and identify generators (Gröbner) — algebraic form of Crit 2.
4. Keep catalogues as regression: every new theorem class must recover BJ-shaped catalogue hits.

_Generated by theorem_attack.py_