# Criterion 3 deepen — ternary / HQCC sign character search

_Elapsed: 3.08s_

**Verdict:** Criterion 3 deepen (3.08s). irr=1046, disc□=3, baseline P(□)=0.0029. Best bucket: {'character': 'sf_kernel_sign', 'value': '1', 'n': 413, 'P_sq': 0.007263922518159807, 'lift_over_baseline': 0.004395853684507799}. Rate-1 buckets (n≥10): 0. Pure-even control disc□=40/40. Base M disc□=False. No ternary/HQCC character forces disc□ at rate 1 on unrestricted T. Crit-3 necessity fragment: not obtained.

---

## Goal

Find a quadratic character \(\chi\) built from ternary / HQCC data such that
\(\operatorname{sgn}\circ\rho=\chi\) (or \(=1\)) on monodromy of \(T\)-specialisations,
i.e. disc square is forced or equal to a computable HQCC invariant.

Prior result (`CRITERION3_SIGN.md`): ternary weight and \(\det\) sign do **not** force disc□.

---

## 1. Sample

- Irreducible monic deg-5 \(\chi_T\): **1046**
- Disc square: **3**
- Baseline \(P(\square)\): **0.002868**

---

## 2. Candidate characters (conditional \(P(\square\mid\chi=v)\))

### `leg_det_3`

- Best \(n\ge 20\): value=`-1`, n=233, P(□)=0.0043, lift=0.001423776659910224
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 233 | 0.0043 |
| 0 | 581 | 0.0034 |
| 1 | 232 | 0.0000 |

### `leg_det_5`

- Best \(n\ge 20\): value=`-1`, n=494, P(□)=0.0040, lift=0.0011805141622994093
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 494 | 0.0040 |
| 1 | 552 | 0.0018 |

### `leg_det_61`

- Best \(n\ge 20\): value=`1`, n=1046, P(□)=0.0029, lift=0.0
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 1 | 1046 | 0.0029 |

### `leg_prod_3`

- Best \(n\ge 20\): value=`0`, n=966, P(□)=0.0031, lift=0.00023752122845979327
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 37 | 0.0000 |
| 0 | 966 | 0.0031 |
| 1 | 43 | 0.0000 |

### `leg_prod_5`

- Best \(n\ge 20\): value=`1`, n=472, P(□)=0.0042, lift=0.0013692193019412126
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 474 | 0.0021 |
| 0 | 100 | 0.0000 |
| 1 | 472 | 0.0042 |

### `leg_const_3`

- Best \(n\ge 20\): value=`1`, n=233, P(□)=0.0043, lift=0.001423776659910224
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 232 | 0.0000 |
| 0 | 581 | 0.0034 |
| 1 | 233 | 0.0043 |

### `leg_const_5`

- Best \(n\ge 20\): value=`-1`, n=494, P(□)=0.0040, lift=0.0011805141622994093
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 494 | 0.0040 |
| 1 | 552 | 0.0018 |

### `leg_x2_3`

- Best \(n\ge 20\): value=`0`, n=894, P(□)=0.0034, lift=0.00048763586433456935
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 76 | 0.0000 |
| 0 | 894 | 0.0034 |
| 1 | 76 | 0.0000 |

### `sf_kernel_sign`

- Best \(n\ge 20\): value=`1`, n=413, P(□)=0.0073, lift=0.004395853684507799
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 633 | 0.0000 |
| 1 | 413 | 0.0073 |

### `v3_prod_parity`

- Best \(n\ge 20\): value=`0`, n=522, P(□)=0.0038, lift=0.0009633487908690648
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 0 | 522 | 0.0038 |
| 1 | 524 | 0.0019 |

### `det_sign`

- Best \(n\ge 20\): value=`-1`, n=481, P(□)=0.0042, lift=0.0012899353243521505
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 481 | 0.0042 |
| 1 | 565 | 0.0018 |

### `ternary_weight_mod2`

- Best \(n\ge 20\): value=`0`, n=517, P(□)=0.0039, lift=0.0010004031199263288
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 0 | 517 | 0.0039 |
| 1 | 529 | 0.0019 |

### `a_mod3`

- Best \(n\ge 20\): value=`0`, n=1046, P(□)=0.0029, lift=0.0
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 0 | 1046 | 0.0029 |

### `has_model_61`

- Best \(n\ge 20\): value=`0`, n=847, P(□)=0.0035, lift=0.000673843799169716
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 0 | 847 | 0.0035 |
| 1 | 199 | 0.0000 |

### `has_model_80`

- Best \(n\ge 20\): value=`0`, n=946, P(□)=0.0032, lift=0.0003031785236418611
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| 0 | 946 | 0.0032 |
| 1 | 100 | 0.0000 |

### `chi_flux`

- Best \(n\ge 20\): value=`-1`, n=233, P(□)=0.0043, lift=0.001423776659910224
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 233 | 0.0043 |
| 0 | 581 | 0.0034 |
| 1 | 232 | 0.0000 |

### `chi_ternary_det`

- Best \(n\ge 20\): value=`1`, n=524, P(□)=0.0038, lift=0.0009487250594777631
- Rate-1 buckets: **none**

| value | n | P(□) |
|------:|--:|-----:|
| -1 | 522 | 0.0019 |
| 1 | 524 | 0.0038 |

### Combo \(\chi_{\mathrm{flux}}=1\) and ternary-weight even

- n=113, P(□)=0.0000

---

## 3. Controls

- **Pure-even LSW control:** disc□ on **40/40** sampled \(m\) (identity).
- **Base \(M\):** disc□=**False** (odd monodromy despite ternary/flux).

---

## 4. Conclusion

1. **No candidate character** achieves \(P(\square\mid\chi=v)=1\) on a nontrivial unrestricted \(T\) bucket.
2. Best lifts remain small (order of baseline ~0.5%); no HQCC sign character found.
3. Pure-even continues to force disc□ by **classical identity**, not by a ternary character of \(T\).
4. Criterion 3 remains **open**; under programme stance it is **paused** as a necessity route.

```bash
python criterion3_deepen.py
```

_Generated by criterion3_deepen.py_