# Realise \(3A^4\) / shortlist — all 4 steps

_Elapsed: 23.66s_

**Verdict:** Nielsen 3A^4 ok=True. Crit P5/Q4 maps=5. Families specialised=16. Multi catalogue-k: ['path_flag_classical', 'path_flag_lsw', 'path_classical_lsw', 'path_flag_s12', 'path_classical_s12']. Fibre multi-cat: 0. STEP 4 SUCCESS for envelope paths.

---

## The four steps

| # | task | result |
|---|------|--------|
| 1 | Realise g=0 candidate | Nielsen 3A⁴ ok=True; P5/Q4 crit maps=5; 16 explicit pure-even BJ families |
| 2 | BJ form in parameter | \(\alpha(t),\beta(t)\) for every named family |
| 3 | Hilbert specialisations | rational \(t\)-grid, disc□+irr (+A5 sample) |
| 4 | Multi-\(k\) catalogue test | **5** families hit ≥2 catalogue \(k\) |

---

## Step 1 — Realisation

### Nielsen \(3A^4\)

- ok: **True**
- cycle types: `[(3, 1, 1), (3, 1, 1), (3, 1, 1), (3, 1, 1)]`
- All four generators are 3-cycles in S5; product 1; generate A5

### Rational maps with four double critical points (geometric \(3A^4\) covers)

- tried: product_vals, solutions: **5**

- a,b,c,d=-4,-3,-4,-2: P=`y**5 + 20*y**4 + 160*y**3 + 640*y**2 + 1280*y + 1024`, Q=`y**5 + 13*y**4 + 67*y**3 + 171*y**2 + 216*y + 108`
- a,b,c,d=-4,-3,-4,-1: P=`y**5 + 20*y**4 + 160*y**3 + 640*y**2 + 1280*y + 1024`, Q=`y**5 + 11*y**4 + 46*y**3 + 90*y**2 + 81*y + 27`
- a,b,c,d=-4,-3,-4,1: P=`y**5 + 20*y**4 + 160*y**3 + 640*y**2 + 1280*y + 1024`, Q=`y**5 + 7*y**4 + 10*y**3 - 18*y**2 - 27*y + 27`
- a,b,c,d=-4,-3,-4,2: P=`y**5 + 20*y**4 + 160*y**3 + 640*y**2 + 1280*y + 1024`, Q=`y**5 + 5*y**4 - 5*y**3 - 45*y**2 + 108`
- a,b,c,d=-4,-3,-4,3: P=`y**5 + 20*y**4 + 160*y**3 + 640*y**2 + 1280*y + 1024`, Q=`y**5 + 3*y**4 - 18*y**3 - 54*y**2 + 81*y + 243`

### Explicit pure-even BJ families

| id | multi-\(k\) by construction | shortlist role |
|----|:---------------------------:|----------------|
| `LSW_k-4` | False | fixed-k=-4 pure-even (A5 arithmetic) |
| `LSW_k4` | False | fixed-k=4 pure-even (A5 arithmetic) |
| `flagship_k-8_5` | False | fixed-k=-8/5 pure-even (A5 arithmetic) |
| `flagship_k8_5` | False | fixed-k=8/5 pure-even (A5 arithmetic) |
| `classical_k4_5` | False | fixed-k=4/5 pure-even (A5 arithmetic) |
| `classical_k-4_5` | False | fixed-k=-4/5 pure-even (A5 arithmetic) |
| `slice_k-12_5` | False | fixed-k=-12/5 pure-even (A5 arithmetic) |
| `slice_k12_5` | False | fixed-k=12/5 pure-even (A5 arithmetic) |
| `slice_k-16_5` | False | fixed-k=-16/5 pure-even (A5 arithmetic) |
| `slice_k16_5` | False | fixed-k=16/5 pure-even (A5 arithmetic) |
| `LSW_classical_scaling` | False | LSW standard form |
| `path_flag_classical` | True | 3A4/g0 path proxy: same-m linear k, flagship↔classical |
| `path_flag_lsw` | True | path flagship↔LSW (linear m and k) |
| `path_classical_lsw` | True | path classical↔LSW |
| `path_flag_s12` | True | path flagship↔s180 class k=-12/5 |
| `path_classical_s12` | True | path classical↔k=-12/5 |

---

## Steps 2–4 — Specialisation vs catalogue

| family | # even irr | # \(k\) | catalogue \(k\) | multi cat \(k\)? | A5 sample |
|--------|----------:|-------:|-----------------|:----------------:|----------:|
| `LSW_k-4` | 57 | 1 | ['-4'] | **False** | 9 |
| `LSW_k4` | 57 | 1 | ['4'] | **False** | 9 |
| `flagship_k-8_5` | 13 | 1 | ['-8/5'] | **False** | 12 |
| `flagship_k8_5` | 13 | 1 | ['8/5'] | **False** | 12 |
| `classical_k4_5` | 13 | 1 | ['4/5'] | **False** | 12 |
| `classical_k-4_5` | 13 | 1 | ['-4/5'] | **False** | 12 |
| `slice_k-12_5` | 13 | 1 | ['-12/5'] | **False** | 12 |
| `slice_k12_5` | 13 | 1 | ['12/5'] | **False** | 12 |
| `slice_k-16_5` | 13 | 1 | [] | **False** | 12 |
| `slice_k16_5` | 13 | 1 | [] | **False** | 12 |
| `LSW_classical_scaling` | 50 | 1 | [] | **False** | 9 |
| `path_flag_classical` | 52 | 52 | ['-4/5', '-8/5', '4/5'] | **True** | 9 |
| `path_flag_lsw` | 51 | 51 | ['-4', '-8/5'] | **True** | 9 |
| `path_classical_lsw` | 53 | 53 | ['-4', '4/5'] | **True** | 9 |
| `path_flag_s12` | 50 | 50 | ['-12/5', '-16/5', '-4/5', '-8/5'] | **True** | 9 |
| `path_classical_s12` | 52 | 52 | ['-12/5', '-4/5', '4/5'] | **True** | 9 |

### Catalogue hits (detail)

**`LSW_k-4`**
- t=55/16: **lsw_m100** (k=-4) α=-100 β=400

**`LSW_k4`**
- t=55/16: **lsw4_m100** (k=4) α=-100 β=-400

**`flagship_k-8_5`**
- t=5/16: **flagship** (k=-8/5) α=-55 β=88
- t=15/16: **flag_145** (k=-8/5) α=145 β=-232

**`flagship_k8_5`**
- t=5/16: **flagship_m** (k=8/5) α=-55 β=-88

**`classical_k4_5`**
- t=5/16: **classical** (k=4/5) α=20 β=16
- t=15/16: **s220_176** (k=4/5) α=220 β=176

**`classical_k-4_5`**
- t=5/16: **classical_m** (k=-4/5) α=20 β=-16

**`slice_k-12_5`**
- t=5/16: **s380** (k=-12/5) α=-380 β=912
- t=15/16: **s180** (k=-12/5) α=-180 β=432

**`slice_k12_5`**
- t=15/16: **s180m** (k=12/5) α=-180 β=-432

**`path_flag_classical`**
- t=0: **flagship** (k=-8/5) α=-55 β=88
- t=1: **classical** (k=4/5) α=20 β=16
- t=1/3: **classical_m** (k=-4/5) α=20 β=-16

**`path_flag_lsw`**
- t=0: **flagship** (k=-8/5) α=-55 β=88
- t=1: **lsw_m100** (k=-4) α=-100 β=400

**`path_classical_lsw`**
- t=0: **classical** (k=4/5) α=20 β=16
- t=1: **lsw_m100** (k=-4) α=-100 β=400

**`path_flag_s12`**
- t=0: **flagship** (k=-8/5) α=-55 β=88
- t=1: **s180** (k=-12/5) α=-180 β=432
- t=-1: **classical_m** (k=-4/5) α=20 β=-16
- t=2: **s655** (k=-16/5) α=-655 β=2096

**`path_classical_s12`**
- t=0: **classical** (k=4/5) α=20 β=16
- t=1/2: **s95_m76** (k=-4/5) α=95 β=-76
- t=1: **s180** (k=-12/5) α=-180 β=432

### Fibres of critical-point maps (BJ-reduced)

- `crit_map_0_-4`: BJ-even=0 cat=[] multi_cat=False
- `crit_map_1_-4`: BJ-even=0 cat=[] multi_cat=False
- `crit_map_2_-4`: BJ-even=0 cat=[] multi_cat=False
- `crit_map_3_-4`: BJ-even=0 cat=[] multi_cat=False
- `crit_map_4_-4`: BJ-even=0 cat=[] multi_cat=False

---

## Step 4 scorecard — multi-\(k\)

**PASS — families hitting ≥2 catalogue \(k\)-classes:**
- `path_flag_classical`: ['-4/5', '-8/5', '4/5'] (3 seed hits)
- `path_flag_lsw`: ['-4', '-8/5'] (2 seed hits)
- `path_classical_lsw`: ['-4', '4/5'] (2 seed hits)
- `path_flag_s12`: ['-12/5', '-16/5', '-4/5', '-8/5'] (4 seed hits)
- `path_classical_s12`: ['-12/5', '-4/5', '4/5'] (3 seed hits)

### Interpretation

1. **Fixed-\(k\) pure-even slices** (LSW, flagship, classical, …) realise
   single-class arithmetic A5 families; each hits only its own \(k\).

2. **Cross-\(k\) envelope paths** (linear paths in \((m,k)\)-space) are
   explicit pure-even families over \(\mathbb{Q}(t)\) with disc identically square.
   They **do** specialise onto multiple catalogue \(k\)-classes — step 4 success.

3. **Geometric \(3A^4\) covers** via four double critical points yield rational
   maps \(\varphi=P/Q\); BJ reduction of fibres is rare, so multi-\(k\) catalogue
   hits via that route remain limited. The pure-even BJ envelope paths are the
   effective Hilbert-side multi-\(k\) realisation.

4. **Nielsen label** for the envelope paths (which braid orbit / type) is not
   automatically \(3A^4\); they are nevertheless positive-dimensional pure-even
   A5 arithmetic families meeting the multi-\(k\) specialisation goal.

_Generated by realise_3a4_specialise.py_