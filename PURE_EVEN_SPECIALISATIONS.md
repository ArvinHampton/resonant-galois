# Pure-even specialisations

_Elapsed: 86.85s_

**Verdict:** Pure-even specialisations (86.85s). Z-pts=1728, irr=1728, even_fail=0, A5_checked=152; paths multi-k id=True; homogenisation all_even=True; contrast locked. PASS.

More data on the **even** side after rigid \(t=3\) locked the **odd** negative control.

Formula (theorem): for fixed \(k\in\mathbb{Q}\setminus\{0\}\),
$$\alpha(m)=256m^2-\frac{3125\,k^4}{256},\qquad\beta(m)=k\cdot\alpha(m),\qquad\operatorname{disc}=(256\,\alpha(m)^2 m)^2.$$

---

## Contrast (both sides locked)

| side | example | disc | parity | Gal |
|------|---------|------|--------|-----|
| **Pure-even** | \(x^5-55x+88\) | `58564000000` = □ | even | \(A_5\) |
| **Rigid \(t=3\)** | monic(\(\varphi-3\)) | `3125/36` = \(5\cdot\square\) | odd | \(S_5\) |

*pure-even resonant slices (disc identically square, Gal≤A5 / A5) ↔ rigid φ fibres (disc=5·□, odd, typically S5)*

---

## 1. \(k\)-slice specialisations

| \(k\) | name | #Z | irr | red | even fail | A5 / checked | sample hits |
|------|------|---:|----:|----:|----------:|-------------:|-------------|
| -4 | LSW | 540 | 540 | 0 | 0 | 20/20 | m=-1/8→(-3121,12484); m=1/8→(-3121,12484) |
| 4 | LSW_flip | 540 | 540 | 0 | 0 | 20/20 | m=-1/8→(-3121,-12484); m=1/8→(-3121,-12484) |
| -8/5 | flagship | 108 | 108 | 0 | 0 | 18/20 | m=-5/4→(320,-512); m=5/4→(320,-512) |
| 8/5 | flagship_flip | 108 | 108 | 0 | 0 | 18/20 | m=-5/4→(320,512); m=5/4→(320,512) |
| 4/5 | classical | 108 | 108 | 0 | 0 | 20/20 | m=-5/8→(95,76); m=5/8→(95,76) |
| -4/5 | classical_flip | 108 | 108 | 0 | 0 | 20/20 | m=-5/8→(95,-76); m=5/8→(95,-76) |
| -12/5 | s12 | 108 | 108 | 0 | 0 | 18/20 | m=-5/8→(-305,732); m=5/8→(-305,732) |
| 12/5 | s12_flip | 108 | 108 | 0 | 0 | 18/20 | m=-5/8→(-305,-732); m=5/8→(-305,-732) |

**Totals:** Z=1728, irr=1728, even_fail=**0**, A5 among Gal checks=**152**

### Flagship \(k=-8/5\) A5 hits (sample)

- m=`-5/4`: `x^5 + (320)x + (-512)` — **HIT_A5** disc=1073741824000000
- m=`5/4`: `x^5 + (320)x + (-512)` — **HIT_A5** disc=1073741824000000
- m=`-15/8`: `x^5 + (820)x + (-1312)` — **HIT_A5** disc=104168853504000000
- m=`15/8`: `x^5 + (820)x + (-1312)` — **HIT_A5** disc=104168853504000000
- m=`-5/2`: `x^5 + (1520)x + (-2432)` — **HIT_A5** disc=2186423566336000000
- m=`5/2`: `x^5 + (1520)x + (-2432)` — **HIT_A5** disc=2186423566336000000
- m=`-25/8`: `x^5 + (2420)x + (-3872)` — **HIT_A5** disc=21950349414400000000
- m=`25/8`: `x^5 + (2420)x + (-3872)` — **HIT_A5** disc=21950349414400000000

### LSW \(k=-4\) A5 hits (sample)

- m=`-1/8`: `x^5 + (-3121)x + (12484)` — **HIT_A5**
- m=`1/8`: `x^5 + (-3121)x + (12484)` — **HIT_A5**
- m=`-1/4`: `x^5 + (-3109)x + (12436)` — **HIT_A5**
- m=`1/4`: `x^5 + (-3109)x + (12436)` — **HIT_A5**
- m=`-3/8`: `x^5 + (-3089)x + (12356)` — **HIT_A5**
- m=`3/8`: `x^5 + (-3089)x + (12356)` — **HIT_A5**

### Classical \(k=4/5\) A5 hits (sample)

- m=`-5/8`: `x^5 + (95)x + (76)` — **HIT_A5**
- m=`5/8`: `x^5 + (95)x + (76)` — **HIT_A5**
- m=`-5/4`: `x^5 + (395)x + (316)` — **HIT_A5**
- m=`5/4`: `x^5 + (395)x + (316)` — **HIT_A5**
- m=`-15/8`: `x^5 + (895)x + (716)` — **HIT_A5**
- m=`15/8`: `x^5 + (895)x + (716)` — **HIT_A5**

---

## 2. Cross-\(k\) pure-even paths

| path | disc id | multi catalogue \(k\) | hist | A5 samples |
|------|:-------:|:---------------------:|------|------------|
| flag_classical | **True** | **True** ['-8/5', '-4/5', '4/5'] | `{'even_Z': 4, 'irr': 3, 'HIT_A5': 3, 'non_Z': 27, 'reducible': 1}` | 3 |
| flag_lsw | **True** | **True** ['-4', '-8/5'] | `{'even_Z': 2, 'irr': 2, 'HIT_A5': 2, 'non_Z': 29}` | 2 |
| classical_lsw | **True** | **True** ['-4', '-8/5', '4/5'] | `{'even_Z': 3, 'irr': 3, 'HIT_A5': 3, 'non_Z': 28}` | 3 |

### Path specialisations (A5)

- **flag_classical** u=`0` k=`-8/5`: x^5+(-55)x+(88) — HIT_A5
- **flag_classical** u=`1/3` k=`-4/5`: x^5+(20)x+(-16) — HIT_A5
- **flag_classical** u=`1` k=`4/5`: x^5+(20)x+(16) — HIT_A5
- **flag_lsw** u=`0` k=`-8/5`: x^5+(-55)x+(88) — HIT_A5
- **flag_lsw** u=`1` k=`-4`: x^5+(-100)x+(400) — HIT_A5
- **classical_lsw** u=`0` k=`4/5`: x^5+(20)x+(16) — HIT_A5
- **classical_lsw** u=`1/2` k=`-8/5`: x^5+(820)x+(-1312) — HIT_A5
- **classical_lsw** u=`1` k=`-4`: x^5+(-100)x+(400) — HIT_A5

---

## 3. Homogenised seed families

Classical lemma proved: **True** — `For all t in Z\{0}, disc(x^5+20 t^4 x+16 t^5) = t^{20} * disc(x^5+20x+16) is a square in Z, because t^{20}=(t^{10})^2 an…`

General: if disc(seed) is square, then \(f_t=x^5+\alpha t^4 x+\beta t^5\) has disc \(=t^{20}\operatorname{disc}(\mathrm{seed})\) square.

| seed | \(k\) | family | all sample even | specs |
|------|------|--------|:---------------:|-------|
| (-100, 400) | -4 | `x^5 + (-100) t^4 x + (400) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (124, -496) | -4 | `x^5 + (124) t^4 x + (-496) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (-100, -400) | 4 | `x^5 + (-100) t^4 x + (-400) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (124, 496) | 4 | `x^5 + (124) t^4 x + (496) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (-55, 88) | -8/5 | `x^5 + (-55) t^4 x + (88) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (145, -232) | -8/5 | `x^5 + (145) t^4 x + (-232) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (-55, -88) | 8/5 | `x^5 + (-55) t^4 x + (-88) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (145, 232) | 8/5 | `x^5 + (145) t^4 x + (232) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (20, 16) | 4/5 | `x^5 + (20) t^4 x + (16) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (95, 76) | 4/5 | `x^5 + (95) t^4 x + (76) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (20, -16) | -4/5 | `x^5 + (20) t^4 x + (-16) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |
| (95, -76) | -4/5 | `x^5 + (95) t^4 x + (-76) t^5` | **True** | t=2:HIT_A5, t=3:HIT_A5, t=5:HIT_A5, t=9:even |

---

## Scorecard

| check | pass |
|-------|:----:|
| even_fail = 0 on all slices | **True** |
| A5 harvest ≥ 20 | **True** (152) |
| paths disc id + multi-\(k\) | **True** |
| homogenisation samples even | **True** |
| classical homo lemma | **True** |
| contrast vs rigid \(t=3\) | **True** |
| **Pure-even data PASS** | **True** |

```bash
python pure_even_specialisations.py
```

_Generated by pure_even_specialisations.py_