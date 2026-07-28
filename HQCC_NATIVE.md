# HQCC-native monodromy — Criterion 1 attack

_Elapsed: 66.49s_

## Goal

Build a one-parameter family \(f_t\in\mathbb{Q}(t)[x]\) whose coefficients
are **generated only from HQCC / 9 Maths data** (ternary branches, model core,
flux/period integers), such that:

1. \(\operatorname{disc}(f_t)\) is a square for all \(t\neq 0\) (even monodromy), and
2. geometric / generic specialisation has Gal \(A_5\) (or \(A_n\)).

Classical reference (not HQCC-native): \(x^5+20t^4 x+16t^5\).

---

## HQCC data used

### Branches (Möbius \(n\mapsto (An+B)/C\))

| Name | (A,B,C) | Map |
|------|---------|-----|
| div3 | (1,0,3) | \(n/3\) |
| Ad_plus | (3,1,1) | \(3n+1\) |
| Ad_minus | (3,-1,1) | \(3n-1\) |
| Syr1 | (4,2,3) | \((4n+2)/3\) |
| Syr2 | (2,1,3) | \((2n+1)/3\) |

### Model core / seeds
`[1, 3, 9, 18, 20, 21, 27, 61, 80, 223, 243, 520, 539, 4880]`

### Lattice size: **2012** integers

---

## 1. Diophantine BJ search on HQCC lattice

- tested pairs: 40001
- square-disc irr: **8**
- A5 seeds: **2**
- D5 seeds: 6
  - a=20 b=-16: `x**5 + 20*x - 16` S5TransitiveSubgroups.A5
  - a=20 b=16: `x**5 + 20*x + 16` S5TransitiveSubgroups.A5

### Parametric even-seed search (α,β HQCC ⇒ disc square)
- tested: 260640
- even seeds (α,β): **24**
- A5 seeds: **16**

- **A5 seed** a=20 b=-16 disc=1024000000: `x**5 + 20*x - 16`
- **A5 seed** a=20 b=16 disc=1024000000: `x**5 + 20*x + 16`
- **A5 seed** a=-55 b=-88 disc=58564000000: `x**5 - 55*x - 88`
- **A5 seed** a=-55 b=88 disc=58564000000: `x**5 - 55*x + 88`
- **A5 seed** a=95 b=-76 disc=2085136000000: `x**5 + 95*x - 76`
- **A5 seed** a=95 b=76 disc=2085136000000: `x**5 + 95*x + 76`
- **A5 seed** a=95 b=-532 disc=252301456000000: `x**5 + 95*x - 532`
- **A5 seed** a=95 b=532 disc=252301456000000: `x**5 + 95*x + 532`
- **A5 seed** a=-100 b=-400 disc=77440000000000: `x**5 - 100*x - 400`
- **A5 seed** a=-100 b=400 disc=77440000000000: `x**5 - 100*x + 400`
- **A5 seed** a=124 b=-496 disc=196642060959744: `x**5 + 124*x - 496`
- **A5 seed** a=124 b=496 disc=196642060959744: `x**5 + 124*x + 496`
- **A5 seed** a=-180 b=-432 disc=60466176000000: `x**5 - 180*x - 432`
- **A5 seed** a=-180 b=432 disc=60466176000000: `x**5 - 180*x + 432`
- **A5 seed** a=220 b=-528 disc=374809600000000: `x**5 + 220*x - 528`
- **A5 seed** a=220 b=528 disc=374809600000000: `x**5 + 220*x + 528`
- even (not A5 yet) a=-5 b=-12 disc=64000000
- even (not A5 yet) a=-5 b=12 disc=64000000
- even (not A5 yet) a=11 b=-44 disc=11754029056
- even (not A5 yet) a=11 b=44 disc=11754029056
- even (not A5 yet) a=20 b=-32 disc=4096000000
- even (not A5 yet) a=20 b=32 disc=4096000000
- even (not A5 yet) a=-80 b=-384 disc=67108864000000
- even (not A5 yet) a=-80 b=384 disc=67108864000000

## 2. Homogenised HQCC-native families

### seed (20, -16)
- family: `x**5 + (20)*t**4*x + (-16)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 1024000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 25920*x + 124416` HIT_A5
  - t=-5: `x**5 + 12500*x + 50000` HIT_A5
  - t=-4: `x**5 + 5120*x + 16384` HIT_A5
  - t=-3: `x**5 + 1620*x + 3888` HIT_A5
  - t=-2: `x**5 + 320*x + 512` HIT_A5

### seed (20, 16)
- family: `x**5 + (20)*t**4*x + (16)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 1024000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 25920*x - 124416` HIT_A5
  - t=-5: `x**5 + 12500*x - 50000` HIT_A5
  - t=-4: `x**5 + 5120*x - 16384` HIT_A5
  - t=-3: `x**5 + 1620*x - 3888` HIT_A5
  - t=-2: `x**5 + 320*x - 512` HIT_A5

### seed (-55, -88)
- family: `x**5 + (-55)*t**4*x + (-88)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 58564000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 71280*x + 684288` HIT_A5
  - t=-5: `x**5 - 34375*x + 275000` HIT_A5
  - t=-4: `x**5 - 14080*x + 90112` HIT_A5
  - t=-3: `x**5 - 4455*x + 21384` HIT_A5
  - t=-2: `x**5 - 880*x + 2816` HIT_A5

### seed (-55, 88)
- family: `x**5 + (-55)*t**4*x + (88)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 58564000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 71280*x - 684288` HIT_A5
  - t=-5: `x**5 - 34375*x - 275000` HIT_A5
  - t=-4: `x**5 - 14080*x - 90112` HIT_A5
  - t=-3: `x**5 - 4455*x - 21384` HIT_A5
  - t=-2: `x**5 - 880*x - 2816` HIT_A5

### seed (95, -76)
- family: `x**5 + (95)*t**4*x + (-76)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 2085136000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 123120*x + 590976` HIT_A5
  - t=-5: `x**5 + 59375*x + 237500` HIT_A5
  - t=-4: `x**5 + 24320*x + 77824` HIT_A5
  - t=-3: `x**5 + 7695*x + 18468` HIT_A5
  - t=-2: `x**5 + 1520*x + 2432` HIT_A5

### seed (95, 76)
- family: `x**5 + (95)*t**4*x + (76)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 2085136000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 123120*x - 590976` HIT_A5
  - t=-5: `x**5 + 59375*x - 237500` HIT_A5
  - t=-4: `x**5 + 24320*x - 77824` HIT_A5
  - t=-3: `x**5 + 7695*x - 18468` HIT_A5
  - t=-2: `x**5 + 1520*x - 2432` HIT_A5

### seed (95, -532)
- family: `x**5 + (95)*t**4*x + (-532)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 252301456000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 123120*x + 4136832` HIT_A5
  - t=-5: `x**5 + 59375*x + 1662500` HIT_A5
  - t=-4: `x**5 + 24320*x + 544768` HIT_A5
  - t=-3: `x**5 + 7695*x + 129276` HIT_A5
  - t=-2: `x**5 + 1520*x + 17024` HIT_A5

### seed (95, 532)
- family: `x**5 + (95)*t**4*x + (532)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 252301456000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 123120*x - 4136832` HIT_A5
  - t=-5: `x**5 + 59375*x - 1662500` HIT_A5
  - t=-4: `x**5 + 24320*x - 544768` HIT_A5
  - t=-3: `x**5 + 7695*x - 129276` HIT_A5
  - t=-2: `x**5 + 1520*x - 17024` HIT_A5

### seed (-100, -400)
- family: `x**5 + (-100)*t**4*x + (-400)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 77440000000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 129600*x + 3110400` HIT_A5
  - t=-5: `x**5 - 62500*x + 1250000` HIT_A5
  - t=-4: `x**5 - 25600*x + 409600` HIT_A5
  - t=-3: `x**5 - 8100*x + 97200` HIT_A5
  - t=-2: `x**5 - 1600*x + 12800` HIT_A5

### seed (-100, 400)
- family: `x**5 + (-100)*t**4*x + (400)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 77440000000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 129600*x - 3110400` HIT_A5
  - t=-5: `x**5 - 62500*x - 1250000` HIT_A5
  - t=-4: `x**5 - 25600*x - 409600` HIT_A5
  - t=-3: `x**5 - 8100*x - 97200` HIT_A5
  - t=-2: `x**5 - 1600*x - 12800` HIT_A5

### seed (124, -496)
- family: `x**5 + (124)*t**4*x + (-496)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 196642060959744
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 160704*x + 3856896` HIT_A5
  - t=-5: `x**5 + 77500*x + 1550000` HIT_A5
  - t=-4: `x**5 + 31744*x + 507904` HIT_A5
  - t=-3: `x**5 + 10044*x + 120528` HIT_A5
  - t=-2: `x**5 + 1984*x + 15872` HIT_A5

### seed (124, 496)
- family: `x**5 + (124)*t**4*x + (496)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 196642060959744
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 160704*x - 3856896` HIT_A5
  - t=-5: `x**5 + 77500*x - 1550000` HIT_A5
  - t=-4: `x**5 + 31744*x - 507904` HIT_A5
  - t=-3: `x**5 + 10044*x - 120528` HIT_A5
  - t=-2: `x**5 + 1984*x - 15872` HIT_A5

### seed (-180, -432)
- family: `x**5 + (-180)*t**4*x + (-432)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 60466176000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 233280*x + 3359232` HIT_A5
  - t=-5: `x**5 - 112500*x + 1350000` HIT_A5
  - t=-4: `x**5 - 46080*x + 442368` HIT_A5
  - t=-3: `x**5 - 14580*x + 104976` HIT_A5
  - t=-2: `x**5 - 2880*x + 13824` HIT_A5

### seed (-180, 432)
- family: `x**5 + (-180)*t**4*x + (432)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 60466176000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 - 233280*x - 3359232` HIT_A5
  - t=-5: `x**5 - 112500*x - 1350000` HIT_A5
  - t=-4: `x**5 - 46080*x - 442368` HIT_A5
  - t=-3: `x**5 - 14580*x - 104976` HIT_A5
  - t=-2: `x**5 - 2880*x - 13824` HIT_A5

### seed (220, -528)
- family: `x**5 + (220)*t**4*x + (-528)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 374809600000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 285120*x + 4105728` HIT_A5
  - t=-5: `x**5 + 137500*x + 1650000` HIT_A5
  - t=-4: `x**5 + 56320*x + 540672` HIT_A5
  - t=-3: `x**5 + 17820*x + 128304` HIT_A5
  - t=-2: `x**5 + 3520*x + 16896` HIT_A5

### seed (220, 528)
- family: `x**5 + (220)*t**4*x + (528)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 374809600000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **21**
- groups: `{'S5TransitiveSubgroups.A5': 21}`
  - t=-6: `x**5 + 285120*x - 4105728` HIT_A5
  - t=-5: `x**5 + 137500*x - 1650000` HIT_A5
  - t=-4: `x**5 + 56320*x - 540672` HIT_A5
  - t=-3: `x**5 + 17820*x - 128304` HIT_A5
  - t=-2: `x**5 + 3520*x - 16896` HIT_A5

### seed (-5, -12)
- family: `x**5 + (-5)*t**4*x + (-12)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 64000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

### seed (-5, 12)
- family: `x**5 + (-5)*t**4*x + (12)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 64000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

### seed (11, -44)
- family: `x**5 + (11)*t**4*x + (-44)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 11754029056
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

### seed (11, 44)
- family: `x**5 + (11)*t**4*x + (44)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 11754029056
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

### seed (20, -32)
- family: `x**5 + (20)*t**4*x + (-32)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 4096000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

### seed (20, 32)
- family: `x**5 + (20)*t**4*x + (32)*t**5`
- disc identity: disc(f_t)=t**20 * disc(seed) = t**20 * 4096000000
- identity_ok=True proved_even=True
- specialisations: 21  A5: **0**
- groups: `{'S5TransitiveSubgroups.D5': 21}`

## 3. HQCC one-parameter constructions (non-BJ)

### Cubic resultant (Z/3 + x=y+m/y)
- **res_s1_m1**: stats=`{'tested': 21, 'irr': 19, 'red': 2}` groups=`{'S3TransitiveSubgroups.S3': 17, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=-1: `x**3 - 3*x**2 - 6*x + 17` HIT_A3
  - t=1: `x**3 + 3*x**2 - 6*x - 17` HIT_A3
- **res_s1_m3**: stats=`{'tested': 21, 'irr': 19, 'red': 2}` groups=`{'S3TransitiveSubgroups.S3': 17, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=-1: `x**3 - 9*x**2 - 12*x + 109` HIT_A3
  - t=1: `x**3 + 9*x**2 - 12*x - 109` HIT_A3
- **res_s3_m1**: stats=`{'tested': 21, 'irr': 20, 'red': 1}` groups=`{'S3TransitiveSubgroups.S3': 18, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=9: `x**3 + 9*x**2 - 972*x - 14661` HIT_A3
  - t=-9: `x**3 - 9*x**2 - 972*x + 14661` HIT_A3
- **res_s3_m3**: stats=`{'tested': 21, 'irr': 20, 'red': 1}` groups=`{'S3TransitiveSubgroups.S3': 18, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=9: `x**3 + 3*x**2 - 18*x - 57` HIT_A3
  - t=-9: `x**3 - 3*x**2 - 18*x + 57` HIT_A3
- **res_s9_m1**: stats=`{'tested': 21, 'irr': 21}` groups=`{'S3TransitiveSubgroups.S3': 20, 'S3TransitiveSubgroups.A3': 1}` A-hits=1
  - t=27: `x**3 + 27*x**2 - 21870*x - 1102977` HIT_A3
- **res_s61_m3**: stats=`{'tested': 21, 'irr': 21}` groups=`{'S3TransitiveSubgroups.S3': 19, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=61: `x**3 + 549*x**2 - 714432*x - 400040989` HIT_A3
  - t=-61: `x**3 - 549*x**2 - 714432*x + 400040989` HIT_A3
- **res_s3_m9**: stats=`{'tested': 21, 'irr': 20, 'red': 1}` groups=`{'S3TransitiveSubgroups.S3': 18, 'S3TransitiveSubgroups.A3': 2}` A-hits=2
  - t=9: `x**3 + 9*x**2 - 36*x - 333` HIT_A3
  - t=-9: `x**3 - 9*x**2 - 36*x + 333` HIT_A3

### Branch-transfer 5×5 matrices
- stats: `{'tested': 380, 'red': 380}`
- even: 0  A5: **0**

### T5 HQCC template lines
- lines with A5: **0**

### Ω-norm lines
- stats: `{'tested': 4096, 'red': 4096}` A5: **0**

### Branch-value deformations
- stats: `{'tested': 56, 'red': 7, 'irr': 49, 'odd': 49}` A5: **0**

## 4. Catalogue overlap (HQCC coefficients)
- catalogue A5: 36
- HQCC-coeff-native: **36**

- `x**5 + x**3 + 3*x**2 - 3` src=LATTICE
- `x**5 + x**3 - 3*x**2 + 3` src=LATTICE
- `x**5 + 9*x**3 - 9*x**2 - 54*x - 81` src=DEFORM_M
- `x**5 + 3*x**3 - 6*x**2 - 9` src=DEFORM_M
- `x**5 - 3*x**3 - 3*x**2 + 9*x + 18` src=DEFORM_M
- `x**5 + 3*x**3 - 12*x**2 - 18*x - 18` src=DEFORM_M
- `x**5 + 3*x**3 - 3*x**2 + 27*x + 72` src=DEFORM_M
- `x**5 + 3*x**3 + 6*x**2 + 9` src=DEFORM_M
- `x**5 - 3*x**3 + 3*x**2 + 9*x - 18` src=DEFORM_M
- `x**5 + 3*x**3 + 3*x**2 + 27*x - 72` src=DEFORM_M
- `x**5 - 9*x**3 - 6*x**2 + 81` src=DEFORM_M
- `x**5 - 9*x**3 + 6*x**2 - 81` src=DEFORM_M
- `x**5 - 3*x**3 + x**2 + x + 3` src=LATTICE
- `x**5 - 3*x**3 - x**2 + x - 3` src=LATTICE
- `x**5 - 3*x**4 - 3*x**3 - x**2 - 3` src=LATTICE
- `x**5 + 3*x**4 - 3*x**3 + x**2 + 3` src=LATTICE
- `x**5 - 3*x**4 + 61*x**3 + 3*x**2 - x + 1` src=LATTICE
- `x**5 - x**4 + 3*x**3 + 61*x**2 - 3*x + 1` src=LATTICE
- `x**5 + x**4 + 3*x**3 - 61*x**2 - 3*x - 1` src=LATTICE
- `x**5 + 3*x**4 + 61*x**3 - 3*x**2 - x - 1` src=LATTICE

---

## 5. Classical reference (not HQCC-native)

```
{
  "seed_poly": "x**5 + 20*x + 16",
  "seed_disc": 1024000000,
  "seed_disc_is_square": true,
  "family": "x**5 + 20*t**4*x + 16*t**5",
  "disc_factored": "1024000000*t**20",
  "disc_over_t20": "1024000000",
  "disc_over_t20_equals_seed": true,
  "theorem": "For all t in Z\\{0}, disc(x^5+20 t^4 x+16 t^5) = t^{20} * disc(x^5+20x+16) is a square in Z, because t^{20}=(t^{10})^2 and disc(seed) is a square. Hence Gal(f_t/Q) \u2264 A5 whenever f_t is irreducible; with a (3,1,1) Frobenius, Gal = A5.",
  "proved": true
}
```

---

## 6. Status / theorem claim

**SUCCESS (partial resolution of Crit 1):** found 18 HQCC-lattice A5 seed(s); 16 homogenised family(ies) with proved even monodromy and A5 specialisations. This is the HQCC-native analogue of the classical (20,16) theorem.

## 7. Next steps

1. If an HQCC (α,β) A5 seed exists: homogenise and promote to the same theorem grade as (20,16).
2. If not: the obstruction is Diophantine — HQCC lattice may not meet the BJ square-disc locus; then the native object must be **non-BJ** (T5 template line or geometric cover).
3. Gröbner form of disc(χ_T5) restricted to HQCC slots.
4. Rigid branch-cycle covers with HQCC-labelled conjugacy classes (geometric monodromy).

_Generated by hqcc_native.py_