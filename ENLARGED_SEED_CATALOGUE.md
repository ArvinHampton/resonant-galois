# Enlarged HQCC A₅ seed catalogue — grouped by \(k=\beta/\alpha\)

_Elapsed: 15.97s_

**Verdict:** Unique A5 BJ seeds: **60**. Distinct ratios \(k=\beta/\alpha\): **16**. Multi-seed k-groups: **10**. Multi-seed pure-even slices: **10**.

---

## Method

1. Build enlarged HQCC lattice (model core, powers of 3, short combos, dense small integers).
2. Scan Bring–Jerrard \(x^5+\alpha x+\beta\): keep disc perfect square; classify Gal.
3. Retain irreducible \(\mathrm{Gal}=A_5\) seeds.
4. Group by reduced rational \(k=\beta/\alpha\).
5. For each multi-seed \(k\), test the pure-even family condition
   \(\alpha + 3125 k^4/256 = \square_{\mathbb{Q}}\) and record the LSW-type family.

General pure-even family on the ray \(\beta=k\alpha\):

$$\alpha(m)=256 m^2 - \frac{3125\,k^4}{256},\qquad \beta(m)=k\cdot\alpha(m)$$

$$\operatorname{disc}=(256\,\alpha(m)^2\, m)^2\quad\text{(identically square in }\mathbb{Q}(m)).$$

---

## Scan stats

| quantity | value |
|----------|------:|
| lattice size | 9318 |
| pairs tested | 4300002 |
| disc□ (all sources) | 72 |
| A5 irreducible | 60 |
| D5 (scan label) | 12 |
| other even irr | 0 |
| k-groups | 16 |
| multi-seed k-groups | 10 |
| multi-seed pure-even slices | 10 |

---

## Multi-seed pure-even slices (primary output)

### \(k = -4\)

- Catalogue A5 seeds on ray: **11**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `256*s**2 - 3125`
- β_true = `4*(3125 - 256*s**2)`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -625 | 2500 | `x**5 - 625*x + 2500` | True | 50 |
| -524 | 2096 | `x**5 - 524*x + 2096` | True | 51 |
| -421 | 1684 | `x**5 - 421*x + 1684` | True | 52 |
| -316 | 1264 | `x**5 - 316*x + 1264` | True | 53 |
| -209 | 836 | `x**5 - 209*x + 836` | True | 54 |
| -100 | 400 | `x**5 - 100*x + 400` | True | 55 |
| 124 | -496 | `x**5 + 124*x - 496` | True | 57 |
| 239 | -956 | `x**5 + 239*x - 956` | True | 58 |
| 356 | -1424 | `x**5 + 356*x - 1424` | True | 59 |
| 475 | -1900 | `x**5 + 475*x - 1900` | True | 60 |
| 596 | -2384 | `x**5 + 596*x - 2384` | True | 61 |

Extra A5 specialisations on the pure-even family:
- m=1: `x**5 - 2869*x + 11476`
- m=1: `x**5 - 2869*x + 11476`
- m=1/16: `x**5 - 3124*x + 12496`
- m=1/16: `x**5 - 3124*x + 12496`
- m=5/16: `x**5 - 3100*x + 12400`
- m=2: `x**5 - 2101*x + 8404`

### \(k = 4\)

- Catalogue A5 seeds on ray: **11**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `256*s**2 - 3125`
- β_true = `4*(256*s**2 - 3125)`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -625 | -2500 | `x**5 - 625*x - 2500` | True | 50 |
| -524 | -2096 | `x**5 - 524*x - 2096` | True | 51 |
| -421 | -1684 | `x**5 - 421*x - 1684` | True | 52 |
| -316 | -1264 | `x**5 - 316*x - 1264` | True | 53 |
| -209 | -836 | `x**5 - 209*x - 836` | True | 54 |
| -100 | -400 | `x**5 - 100*x - 400` | True | 55 |
| 124 | 496 | `x**5 + 124*x + 496` | True | 57 |
| 239 | 956 | `x**5 + 239*x + 956` | True | 58 |
| 356 | 1424 | `x**5 + 356*x + 1424` | True | 59 |
| 475 | 1900 | `x**5 + 475*x + 1900` | True | 60 |
| 596 | 2384 | `x**5 + 596*x + 2384` | True | 61 |

Extra A5 specialisations on the pure-even family:
- m=1: `x**5 - 2869*x - 11476`
- m=1: `x**5 - 2869*x - 11476`
- m=1/16: `x**5 - 3124*x - 12496`
- m=1/16: `x**5 - 3124*x - 12496`
- m=5/16: `x**5 - 3100*x - 12400`
- m=2: `x**5 - 2101*x - 8404`

### \(k = -12/5\)

- Catalogue A5 seeds on ray: **5**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `5*(32000*s**2 - 81)`
- β_true = `972 - 384000*s**2`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -380 | 912 | `x**5 - 380*x + 912` | True | 5 |
| -305 | 732 | `x**5 - 305*x + 732` | True | 10 |
| -180 | 432 | `x**5 - 180*x + 432` | True | 15 |
| 220 | -528 | `x**5 + 220*x - 528` | True | 25 |
| 820 | -1968 | `x**5 + 820*x - 1968` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 220*x - 528`
- m=5/16: `x**5 - 380*x + 912`
- m=25/8: `x**5 + 2095*x - 5028`
- m=5/8: `x**5 - 305*x + 732`
- m=75/16: `x**5 + 5220*x - 12528`
- m=15/16: `x**5 - 180*x + 432`

### \(k = 12/5\)

- Catalogue A5 seeds on ray: **5**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `5*(32000*s**2 - 81)`
- β_true = `384000*s**2 - 972`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -380 | -912 | `x**5 - 380*x - 912` | True | 5 |
| -305 | -732 | `x**5 - 305*x - 732` | True | 10 |
| -180 | -432 | `x**5 - 180*x - 432` | True | 15 |
| 220 | 528 | `x**5 + 220*x + 528` | True | 25 |
| 820 | 1968 | `x**5 + 820*x + 1968` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 220*x + 528`
- m=5/16: `x**5 - 380*x - 912`
- m=25/8: `x**5 + 2095*x + 5028`
- m=5/8: `x**5 - 305*x - 732`
- m=75/16: `x**5 + 5220*x + 12528`
- m=15/16: `x**5 - 180*x - 432`

### \(k = -4/5\)

- Catalogue A5 seeds on ray: **4**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `5*(32000*s**2 - 1)`
- β_true = `4 - 128000*s**2`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| 20 | -16 | `x**5 + 20*x - 16` | True | 5 |
| 95 | -76 | `x**5 + 95*x - 76` | True | 10 |
| 220 | -176 | `x**5 + 220*x - 176` | True | 15 |
| 1220 | -976 | `x**5 + 1220*x - 976` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 620*x - 496`
- m=5/16: `x**5 + 20*x - 16`
- m=25/8: `x**5 + 2495*x - 1996`
- m=5/8: `x**5 + 95*x - 76`
- m=75/16: `x**5 + 5620*x - 4496`
- m=15/16: `x**5 + 220*x - 176`

### \(k = -8/5\)

- Catalogue A5 seeds on ray: **4**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `80*(2000*s**2 - 1)`
- β_true = `128 - 256000*s**2`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -55 | 88 | `x**5 - 55*x + 88` | True | 5 |
| 145 | -232 | `x**5 + 145*x - 232` | True | 15 |
| 320 | -512 | `x**5 + 320*x - 512` | True | 20 |
| 1145 | -1832 | `x**5 + 1145*x - 1832` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 545*x - 872`
- m=5/16: `x**5 - 55*x + 88`
- m=25/8: `x**5 + 2420*x - 3872`
- m=75/16: `x**5 + 5545*x - 8872`
- m=15/16: `x**5 + 145*x - 232`
- m=25/4: `x**5 + 9920*x - 15872`

### \(k = 4/5\)

- Catalogue A5 seeds on ray: **4**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `5*(32000*s**2 - 1)`
- β_true = `128000*s**2 - 4`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| 20 | 16 | `x**5 + 20*x + 16` | True | 5 |
| 95 | 76 | `x**5 + 95*x + 76` | True | 10 |
| 220 | 176 | `x**5 + 220*x + 176` | True | 15 |
| 1220 | 976 | `x**5 + 1220*x + 976` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 620*x + 496`
- m=5/16: `x**5 + 20*x + 16`
- m=25/8: `x**5 + 2495*x + 1996`
- m=5/8: `x**5 + 95*x + 76`
- m=75/16: `x**5 + 5620*x + 4496`
- m=15/16: `x**5 + 220*x + 176`

### \(k = 8/5\)

- Catalogue A5 seeds on ray: **4**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `80*(2000*s**2 - 1)`
- β_true = `256000*s**2 - 128`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -55 | -88 | `x**5 - 55*x - 88` | True | 5 |
| 145 | 232 | `x**5 + 145*x + 232` | True | 15 |
| 320 | 512 | `x**5 + 320*x + 512` | True | 20 |
| 1145 | 1832 | `x**5 + 1145*x + 1832` | True | 35 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 + 545*x + 872`
- m=5/16: `x**5 - 55*x - 88`
- m=25/8: `x**5 + 2420*x + 3872`
- m=75/16: `x**5 + 5545*x + 8872`
- m=15/16: `x**5 + 145*x + 232`
- m=25/4: `x**5 + 9920*x + 15872`

### \(k = -16/5\)

- Catalogue A5 seeds on ray: **3**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `1280*(125*s**2 - 1)`
- β_true = `4096 - 512000*s**2`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -655 | 2096 | `x**5 - 655*x + 2096` | True | 25 |
| -55 | 176 | `x**5 - 55*x + 176` | True | 35 |
| 745 | -2384 | `x**5 + 745*x - 2384` | True | 45 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 - 655*x + 2096`
- m=5/16: `x**5 - 1255*x + 4016`
- m=25/8: `x**5 + 1220*x - 3904`
- m=5/8: `x**5 - 1180*x + 3776`
- m=75/16: `x**5 + 4345*x - 13904`
- m=15/16: `x**5 - 1055*x + 3376`

### \(k = 16/5\)

- Catalogue A5 seeds on ray: **3**
- Disc identically square (true Q-family): **True**
- Integer β-model in Z[s]: **True**
- α_true = `1280*(125*s**2 - 1)`
- β_true = `512000*s**2 - 4096`

| α | β | poly | on family | m |
|--:|--:|------|:---------:|---|
| -655 | -2096 | `x**5 - 655*x - 2096` | True | 25 |
| -55 | -176 | `x**5 - 55*x - 176` | True | 35 |
| 745 | 2384 | `x**5 + 745*x + 2384` | True | 45 |

Extra A5 specialisations on the pure-even family:
- m=25/16: `x**5 - 655*x - 2096`
- m=5/16: `x**5 - 1255*x - 4016`
- m=25/8: `x**5 + 1220*x + 3904`
- m=5/8: `x**5 - 1180*x - 3776`
- m=75/16: `x**5 + 4345*x + 13904`
- m=15/16: `x**5 - 1055*x - 3376`

---

## All multi-seed k-groups (A5)

- **k=-4**: 11 seeds; pure-even multi: **True** — (-625,2500), (-524,2096), (-421,1684), (-316,1264), (-209,836), (-100,400), (124,-496), (239,-956), (356,-1424), (475,-1900), (596,-2384)
- **k=4**: 11 seeds; pure-even multi: **True** — (-625,-2500), (-524,-2096), (-421,-1684), (-316,-1264), (-209,-836), (-100,-400), (124,496), (239,956), (356,1424), (475,1900), (596,2384)
- **k=-12/5**: 5 seeds; pure-even multi: **True** — (-380,912), (-305,732), (-180,432), (220,-528), (820,-1968)
- **k=12/5**: 5 seeds; pure-even multi: **True** — (-380,-912), (-305,-732), (-180,-432), (220,528), (820,1968)
- **k=-4/5**: 4 seeds; pure-even multi: **True** — (20,-16), (95,-76), (220,-176), (1220,-976)
- **k=-8/5**: 4 seeds; pure-even multi: **True** — (-55,88), (145,-232), (320,-512), (1145,-1832)
- **k=4/5**: 4 seeds; pure-even multi: **True** — (20,16), (95,76), (220,176), (1220,976)
- **k=8/5**: 4 seeds; pure-even multi: **True** — (-55,-88), (145,232), (320,512), (1145,1832)
- **k=-16/5**: 3 seeds; pure-even multi: **True** — (-655,2096), (-55,176), (745,-2384)
- **k=16/5**: 3 seeds; pure-even multi: **True** — (-655,-2096), (-55,-176), (745,2384)

---

## Full A5 catalogue (unique)

| α | β | k=β/α | poly | source |
|--:|--:|-------|------|--------|
| -116 | 1392 | -12 | `x**5 - 116*x + 1392` | scan |
| -380 | 912 | -12/5 | `x**5 - 380*x + 912` | scan |
| -305 | 732 | -12/5 | `x**5 - 305*x + 732` | scan |
| -180 | 432 | -12/5 | `x**5 - 180*x + 432` | known:s180 |
| 220 | -528 | -12/5 | `x**5 + 220*x - 528` | known:s220_m |
| 820 | -1968 | -12/5 | `x**5 + 820*x - 1968` | scan_pass2 |
| -655 | 2096 | -16/5 | `x**5 - 655*x + 2096` | scan_pass2 |
| -55 | 176 | -16/5 | `x**5 - 55*x + 176` | scan |
| 745 | -2384 | -16/5 | `x**5 + 745*x - 2384` | scan_pass2 |
| 95 | -532 | -28/5 | `x**5 + 95*x - 532` | known:s95_m532 |
| -625 | 2500 | -4 | `x**5 - 625*x + 2500` | scan_pass2 |
| -524 | 2096 | -4 | `x**5 - 524*x + 2096` | scan_pass2 |
| -421 | 1684 | -4 | `x**5 - 421*x + 1684` | scan_pass2 |
| -316 | 1264 | -4 | `x**5 - 316*x + 1264` | scan |
| -209 | 836 | -4 | `x**5 - 209*x + 836` | scan |
| -100 | 400 | -4 | `x**5 - 100*x + 400` | known:s100 |
| 124 | -496 | -4 | `x**5 + 124*x - 496` | known:s124_m |
| 239 | -956 | -4 | `x**5 + 239*x - 956` | scan |
| 356 | -1424 | -4 | `x**5 + 356*x - 1424` | scan |
| 475 | -1900 | -4 | `x**5 + 475*x - 1900` | scan_pass2 |
| 596 | -2384 | -4 | `x**5 + 596*x - 2384` | scan_pass2 |
| 20 | -16 | -4/5 | `x**5 + 20*x - 16` | known:classical_m |
| 95 | -76 | -4/5 | `x**5 + 95*x - 76` | known:s95_m76 |
| 220 | -176 | -4/5 | `x**5 + 220*x - 176` | scan |
| 1220 | -976 | -4/5 | `x**5 + 1220*x - 976` | scan_pass2 |
| -271 | 2168 | -8 | `x**5 - 271*x + 2168` | scan_pass2 |
| -55 | 88 | -8/5 | `x**5 - 55*x + 88` | known:flagship |
| 145 | -232 | -8/5 | `x**5 + 145*x - 232` | scan |
| 320 | -512 | -8/5 | `x**5 + 320*x - 512` | scan |
| 1145 | -1832 | -8/5 | `x**5 + 1145*x - 1832` | scan_pass2 |
| -116 | -1392 | 12 | `x**5 - 116*x - 1392` | scan |
| -380 | -912 | 12/5 | `x**5 - 380*x - 912` | scan |
| -305 | -732 | 12/5 | `x**5 - 305*x - 732` | scan |
| -180 | -432 | 12/5 | `x**5 - 180*x - 432` | known:s180_m |
| 220 | 528 | 12/5 | `x**5 + 220*x + 528` | known:s220 |
| 820 | 1968 | 12/5 | `x**5 + 820*x + 1968` | scan_pass2 |
| -655 | -2096 | 16/5 | `x**5 - 655*x - 2096` | scan_pass2 |
| -55 | -176 | 16/5 | `x**5 - 55*x - 176` | scan |
| 745 | 2384 | 16/5 | `x**5 + 745*x + 2384` | scan_pass2 |
| 95 | 532 | 28/5 | `x**5 + 95*x + 532` | known:s95_532 |
| -625 | -2500 | 4 | `x**5 - 625*x - 2500` | scan_pass2 |
| -524 | -2096 | 4 | `x**5 - 524*x - 2096` | scan_pass2 |
| -421 | -1684 | 4 | `x**5 - 421*x - 1684` | scan_pass2 |
| -316 | -1264 | 4 | `x**5 - 316*x - 1264` | scan |
| -209 | -836 | 4 | `x**5 - 209*x - 836` | scan |
| -100 | -400 | 4 | `x**5 - 100*x - 400` | known:s100_m |
| 124 | 496 | 4 | `x**5 + 124*x + 496` | known:s124 |
| 239 | 956 | 4 | `x**5 + 239*x + 956` | scan |
| 356 | 1424 | 4 | `x**5 + 356*x + 1424` | scan |
| 475 | 1900 | 4 | `x**5 + 475*x + 1900` | scan_pass2 |
| 596 | 2384 | 4 | `x**5 + 596*x + 2384` | scan_pass2 |
| 20 | 16 | 4/5 | `x**5 + 20*x + 16` | known:classical |
| 95 | 76 | 4/5 | `x**5 + 95*x + 76` | known:s95_76 |
| 220 | 176 | 4/5 | `x**5 + 220*x + 176` | scan |
| 1220 | 976 | 4/5 | `x**5 + 1220*x + 976` | scan_pass2 |
| -271 | -2168 | 8 | `x**5 - 271*x - 2168` | scan_pass2 |
| -55 | -88 | 8/5 | `x**5 - 55*x - 88` | known:flagship_m |
| 145 | 232 | 8/5 | `x**5 + 145*x + 232` | scan |
| 320 | 512 | 8/5 | `x**5 + 320*x + 512` | scan |
| 1145 | 1832 | 8/5 | `x**5 + 1145*x + 1832` | scan_pass2 |

---

## k-group size histogram

```
{'-4': 11, '4': 11, '-12/5': 5, '12/5': 5, '-8/5': 4, '8/5': 4, '4/5': 4, '-4/5': 4, '16/5': 3, '-16/5': 3, '28/5': 1, '-28/5': 1, '12': 1, '-12': 1, '8': 1, '-8': 1}
```

---

## Conclusions

1. Enlarged catalogue yields more A5 BJ lattice seeds to feed k-grouping.
2. Every multi-seed pure-even slice is an LSW-type family on a fixed ray β=kα;
   LSW itself is k=-4; flagship is the rational-k family k=-8/5.
3. New multi-seed pure-even slices (if any above) are the natural fusion fuel:
   one pure-even 1-parameter family carrying ≥2 HQCC A5 seeds.
4. Still open for geometric fusion: a pure-even family joining flagship to a
   seed with a *different* k (not on the same ray).

_Generated by enlarge_seed_catalogue.py_