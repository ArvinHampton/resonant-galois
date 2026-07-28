# Gap A — BJ families with A5 passport intent

_Elapsed: 107.2s_

## Task

1. Low-dimensional BJ families \(x^5+\alpha(t)x+\beta(t)\).
2. Disc square in \(\mathbb{Q}(t)\), or square × fixed square-free \(c\in\mathbb{Z}\).
3. Specialise at HQCC lattice; match known seeds; count \(A_5\) + 3-cycles.

**Verdict:** Families scanned: 38. Disc square in Q(t): 9 ['H_flagship', 'H_classical', 'H_hqcc_95_76', 'H_hqcc_95_532', 'H_hqcc_100_400', 'H_hqcc_124_496', 'P_flag_plain_t4_t5', 'P_class_plain_t4_t5', 'P_eulerish_m_t4_t5']. Even after const square-free twist: 1. Recover known seeds under lattice specs: 21. A5-rich (≥3 lattice A5 specs): ['H_flagship', 'H_classical', 'H_hqcc_95_76', 'H_hqcc_95_532', 'H_hqcc_100_400', 'H_hqcc_124_496', 'L_class_95', 'P_flag_plain_t4_t5', 'P_class_plain_t4_t5']. Homogenised seeds remain the only systematically pure-even families (disc = t^20 * square). Linear pencils recover seeds at endpoints but are not pure even. No new pure-even BJ family beyond homogenisation was found that both has disc square in Q(t) and hits multiple lattice seeds.

---

## Summary counts

| Metric | Count |
|--------|------:|
| Families | 38 |
| Disc □ in Q(t) | 9 |
| Even after const twist | 1 |
| Recover known seeds | 21 |
| A5-rich (≥3) | 9 |

Pure even IDs: `['H_flagship', 'H_classical', 'H_hqcc_95_76', 'H_hqcc_95_532', 'H_hqcc_100_400', 'H_hqcc_124_496', 'P_flag_plain_t4_t5', 'P_class_plain_t4_t5', 'P_eulerish_m_t4_t5']`
Twist list: `[('P_eulerish_t4_t5', 10)]`
Seed-recovering IDs: `['H_flagship', 'H_classical', 'H_hqcc_95_76', 'H_hqcc_95_532', 'H_hqcc_100_400', 'H_hqcc_124_496', 'L_flag_class', 'L_flag_95_76', 'L_flag_period', 'L_class_95', 'L_flag_100', 'L_95_period', 'P_flag_plain_t4_t5', 'P_class_plain_t4_t5', 'R_quad_flag', 'R_quad_class', 'R_beta_cubic_flag', 'M_flag_to_class_reparam', 'M_flag_steep', 'M_class_to_flag', 'M_punct_flux_line']`

---

## Ranked families

### `H_classical` (homogenised_seed)
- α=`20*t**4`, β=`16*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 4, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=20 β=-16
  - SEED t=1: α=20 β=16
  - SEED t=-2: α=320 β=-512
  - SEED t=2: α=320 β=512
  - A5 t=-1: `x**5 + (20)*x + (-16)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (320)*x + (-512)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (320)*x + (512)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (1620)*x + (-3888)` 3-cyc=True HIT_A5

### `P_class_plain_t4_t5` (power_balance_4_5)
- α=`20*t**4`, β=`16*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 4, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=20 β=-16
  - SEED t=1: α=20 β=16
  - SEED t=-2: α=320 β=-512
  - SEED t=2: α=320 β=512
  - A5 t=-1: `x**5 + (20)*x + (-16)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (320)*x + (-512)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (320)*x + (512)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (1620)*x + (-3888)` 3-cyc=True HIT_A5

### `H_flagship` (homogenised_seed)
- α=`-55*t**4`, β=`88*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=-55 β=-88
  - SEED t=1: α=-55 β=88
  - A5 t=-1: `x**5 + (-55)*x + (-88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (-880)*x + (-2816)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (-880)*x + (2816)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (-4455)*x + (-21384)` 3-cyc=True HIT_A5

### `H_hqcc_95_76` (homogenised_seed)
- α=`95*t**4`, β=`76*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=95 β=-76
  - SEED t=1: α=95 β=76
  - A5 t=-1: `x**5 + (95)*x + (-76)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (76)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (1520)*x + (-2432)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (1520)*x + (2432)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (7695)*x + (-18468)` 3-cyc=True HIT_A5

### `H_hqcc_95_532` (homogenised_seed)
- α=`95*t**4`, β=`532*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=95 β=-532
  - SEED t=1: α=95 β=532
  - A5 t=-1: `x**5 + (95)*x + (-532)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (532)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (1520)*x + (-17024)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (1520)*x + (17024)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (7695)*x + (-129276)` 3-cyc=True HIT_A5

### `H_hqcc_100_400` (homogenised_seed)
- α=`-100*t**4`, β=`400*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=-100 β=-400
  - SEED t=1: α=-100 β=400
  - A5 t=-1: `x**5 + (-100)*x + (-400)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (-100)*x + (400)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (-1600)*x + (-12800)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (-1600)*x + (12800)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (-8100)*x + (-97200)` 3-cyc=True HIT_A5

### `H_hqcc_124_496` (homogenised_seed)
- α=`124*t**4`, β=`496*t**5`
- note: disc = t^20 * disc(seed); square for t≠0 if seed even
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=124 β=-496
  - SEED t=1: α=124 β=496
  - A5 t=-1: `x**5 + (124)*x + (-496)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (124)*x + (496)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (1984)*x + (-15872)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (1984)*x + (15872)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (10044)*x + (-120528)` 3-cyc=True HIT_A5

### `P_flag_plain_t4_t5` (power_balance_4_5)
- α=`-55*t**4`, β=`88*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **True** (deg=20, odd_factors=[])
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'sq': 48, 'A5': 48}`
- known seeds hit: 2, A5: 48, A5+3-cycle: 48
  - SEED t=-1: α=-55 β=-88
  - SEED t=1: α=-55 β=88
  - A5 t=-1: `x**5 + (-55)*x + (-88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=-2: `x**5 + (-880)*x + (-2816)` 3-cyc=True HIT_A5
  - A5 t=2: `x**5 + (-880)*x + (2816)` 3-cyc=True HIT_A5
  - A5 t=-3: `x**5 + (-4455)*x + (-21384)` 3-cyc=True HIT_A5

### `P_eulerish_m_t4_t5` (power_balance_4_5)
- α=`-5*t**4`, β=`4*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **True** (deg=-1, odd_factors=None)
- const twist: found=True c=1 
- lattice specs: `{'tested': 48, 'nonpos_disc': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_eulerish_t4_t5` (power_balance_4_5)
- α=`5*t**4`, β=`4*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=True c=10 
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `L_class_95` (linear_pencil)
- α=`75*t + 20`, β=`60*t + 16`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=5, odd_factors=[('3*t + 1', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 6, 'A5': 6, 'nonpos_disc': 15, 'odd': 28}`
- known seeds hit: 2, A5: 6, A5+3-cycle: 6
  - SEED t=0: α=20 β=16
  - SEED t=1: α=95 β=76
  - A5 t=0: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (76)` 3-cyc=True HIT_A5
  - A5 t=5: `x**5 + (395)*x + (316)` 3-cyc=True HIT_A5
  - A5 t=8: `x**5 + (620)*x + (496)` 3-cyc=True HIT_A5
  - A5 t=16: `x**5 + (1220)*x + (976)` 3-cyc=True HIT_A5

### `L_flag_class` (linear_pencil)
- α=`75*t - 55`, β=`88 - 72*t`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=5, odd_factors=[('759375*t**5 - 2679399*t**4 + 3570534*t**3 - 2053854*t**2 + 331419*t + 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'nonpos_disc': 15, 'odd': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=20 β=16
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5

### `L_flag_95_76` (linear_pencil)
- α=`150*t - 55`, β=`88 - 12*t`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=5, odd_factors=[('24300000*t**5 - 44549919*t**4 + 32667624*t**3 - 11952864*t**2 + 2068374*t + 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'nonpos_disc': 15, 'odd': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=95 β=76
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (76)` 3-cyc=True HIT_A5

### `L_flag_period` (linear_pencil)
- α=`150*t - 55`, β=`444*t + 88`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=5, odd_factors=[('24300000*t**5 + 107257041*t**4 + 153021528*t**3 + 23801184*t**2 + 6923862*t + 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'nonpos_disc': 15, 'odd': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=95 β=532
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (532)` 3-cyc=True HIT_A5

### `L_flag_100` (linear_pencil)
- α=`-45*t - 55`, β=`312*t + 88`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=5, odd_factors=[('59049*t**5 - 36654201*t**4 - 40878486*t**3 - 16589826*t**2 - 2663331*t - 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'odd': 45, 'nonpos_disc': 2}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=-100 β=400
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (-100)*x + (400)` 3-cyc=True HIT_A5

### `L_95_period` (linear_pencil)
- α=`95`, β=`456*t + 76`
- note: contains endpoints at t=0,1
- disc □ in Q(t): **False** (deg=4, odd_factors=[('324*t**4 + 216*t**3 + 54*t**2 + 6*t + 5', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'odd': 47}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=95 β=76
  - SEED t=1: α=95 β=532
  - A5 t=0: `x**5 + (95)*x + (76)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (532)` 3-cyc=True HIT_A5

### `M_flag_to_class_reparam` (linear_general)
- α=`75*t - 55`, β=`88 - 72*t`
- note: general linear BJ
- disc □ in Q(t): **False** (deg=5, odd_factors=[('759375*t**5 - 2679399*t**4 + 3570534*t**3 - 2053854*t**2 + 331419*t + 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'nonpos_disc': 15, 'odd': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=20 β=16
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5

### `M_flag_steep` (linear_general)
- α=`150*t - 55`, β=`88 - 12*t`
- note: general linear BJ
- disc □ in Q(t): **False** (deg=5, odd_factors=[('24300000*t**5 - 44549919*t**4 + 32667624*t**3 - 11952864*t**2 + 2068374*t + 73205', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'nonpos_disc': 15, 'odd': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=-55 β=88
  - SEED t=1: α=95 β=76
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (95)*x + (76)` 3-cyc=True HIT_A5

### `M_class_to_flag` (linear_general)
- α=`20 - 75*t`, β=`72*t + 16`
- note: general linear BJ
- disc □ in Q(t): **False** (deg=5, odd_factors=[('759375*t**5 - 1117476*t**4 + 446688*t**3 - 175104*t**2 + 14592*t - 1280', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 2, 'A5': 2, 'odd': 15, 'nonpos_disc': 32}`
- known seeds hit: 2, A5: 2, A5+3-cycle: 2
  - SEED t=0: α=20 β=16
  - SEED t=1: α=-55 β=88
  - A5 t=0: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5
  - A5 t=1: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5

### `R_quad_flag` (low_deg_poly)
- α=`t**2 - 55`, β=`t**3 + 88`
- note: exploratory low-degree BJ pencil
- disc □ in Q(t): **False** (deg=12, odd_factors=[('3125*t**12 + 256*t**10 + 1100000*t**9 - 70400*t**8 + 152944000*t**6 - 425920000*t**4 + 8518400000*t**3 + 11712800000*t**2 + 58564000000', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 1, 'A5': 1, 'odd': 47, 'nonpos_disc': 1}`
- known seeds hit: 1, A5: 1, A5+3-cycle: 1
  - SEED t=0: α=-55 β=88
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5

### `R_quad_class` (low_deg_poly)
- α=`t**2 + 20`, β=`t**3 + 16`
- note: exploratory low-degree BJ pencil
- disc □ in Q(t): **False** (deg=12, odd_factors=[('3125*t**12 + 256*t**10 + 200000*t**9 + 25600*t**8 + 5824000*t**6 + 20480000*t**4 + 51200000*t**3 + 204800000*t**2 + 1024000000', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 1, 'A5': 1, 'odd': 48}`
- known seeds hit: 1, A5: 1, A5+3-cycle: 1
  - SEED t=0: α=20 β=16
  - A5 t=0: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5

### `R_beta_cubic_flag` (low_deg_poly)
- α=`-55`, β=`t**3 + 88`
- note: exploratory low-degree BJ pencil
- disc □ in Q(t): **False** (deg=12, odd_factors=[('t**12 + 352*t**9 + 46464*t**6 + 2725888*t**3 + 18740480', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'sq': 1, 'A5': 1, 'odd': 46, 'nonpos_disc': 2}`
- known seeds hit: 1, A5: 1, A5+3-cycle: 1
  - SEED t=0: α=-55 β=88
  - A5 t=0: `x**5 + (-55)*x + (88)` 3-cyc=True HIT_A5

### `M_punct_flux_line` (linear_general)
- α=`61 - 41*t`, β=`80 - 64*t`
- note: general linear BJ
- disc □ in Q(t): **False** (deg=5, odd_factors=[('115856201*t**5 - 1066657105*t**4 + 3588550410*t**3 - 5735550610*t**2 + 4438397405*t - 1344596301', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'odd': 19, 'sq': 1, 'A5': 1, 'nonpos_disc': 29}`
- known seeds hit: 1, A5: 1, A5+3-cycle: 1
  - SEED t=1: α=20 β=16
  - A5 t=1: `x**5 + (20)*x + (16)` 3-cyc=True HIT_A5

### `P_ternary_t4_t5` (power_balance_4_5)
- α=`3*t**4`, β=`9*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_punct_3_t4_t5` (power_balance_4_5)
- α=`61*t**4`, β=`3*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_3_punct_t4_t5` (power_balance_4_5)
- α=`3*t**4`, β=`61*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_unit_t4_t5` (power_balance_4_5)
- α=`t**4`, β=`t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_punct_flux_t4_t5` (power_balance_4_5)
- α=`61*t**4`, β=`80*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_flux_punct_t4_t5` (power_balance_4_5)
- α=`80*t**4`, β=`61*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_tower_period_t4_t5` (power_balance_4_5)
- α=`243*t**4`, β=`539*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_3_period_t4_t5` (power_balance_4_5)
- α=`3*t**4`, β=`539*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `P_3tower_t4_t5` (power_balance_4_5)
- α=`9*t**4`, β=`27*t**5`
- note: weighted deg for disc terms both t^{20}
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `R_ternary_qd` (low_deg_poly)
- α=`3*t**2`, β=`9*t**3`
- note: exploratory low-degree BJ pencil
- disc □ in Q(t): **False** (deg=12, odd_factors=[('84375*t**2 + 256', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `R_alpha_quad_beta_const` (low_deg_poly)
- α=`t**2`, β=`88`
- note: exploratory low-degree BJ pencil
- disc □ in Q(t): **False** (deg=10, odd_factors=[('t**10 + 732050000', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'odd': 49}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `C_icosa_adj` (classical_shape)
- α=`5*t**4`, β=`12*t**5`
- note: classical coefficient shape scan
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `C_phi_motif` (classical_shape)
- α=`10*t**4`, β=`6*t**5`
- note: coeffs from preferred Belyi φ motif
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `C_phi_motif2` (classical_shape)
- α=`15*t**4`, β=`10*t**5`
- note: φ motif 15,10
- disc □ in Q(t): **False** (deg=20, odd_factors=[])
- const twist: found=False c=None no small square-free c
- lattice specs: `{'tested': 48, 'odd': 48}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

### `M_ternary_to_punctish` (linear_general)
- α=`58*t + 3`, β=`79*t + 9`
- note: general linear BJ
- disc □ in Q(t): **False** (deg=5, odd_factors=[('168027332608*t**5 + 165174347765*t**4 + 59962267980*t**3 + 9711038430*t**2 + 725900940*t + 20565333', 1)])
- const twist: found=False c=None non-constant factors of odd multiplicity — need different family shape
- lattice specs: `{'tested': 49, 'odd': 34, 'nonpos_disc': 15}`
- known seeds hit: 0, A5: 0, A5+3-cycle: 0

---

## Interpretation for fusion

| Family class | Disc □ in Q(t)? | HQCC seeds | Role |
|--------------|:---------------:|:-----------|-----|
| Homogenised seed \(H_*\) | **Yes** (\(t^{20}\times\mathrm{const}\)) | ray through one seed | Theorem-grade even; arithmetic fusion |
| Linear pencil \(L_*\) | No | endpoints \(t=0,1\) | Equation-level inclusion; not pure A5 cover |
| Power \(P_*\) lattice coeffs | Yes iff seed disc □ | if \((p,q)\) is seed | Same as homogenisation |
| Exploratory \(R_*,M_*,C_*\) | Rarely | occasional | No pure-even multi-seed family found |

### Passport (3A,3A,5A)

For BJ fibres, geometric monodromy of the *family* is not the Belyi passport of \(\varphi\);
the operational proxy is: **even + 3-cycle Frobenius ⇒ A5**. Homogenised families
inherit A5 from the seed along the ray (empirically and by Hilbert for most t).

### Principal open (unchanged)

A pure geometric \(A_5\) family (disc □ identically in the parameter, or square-free
twist only) that Hilbert-recovers **multiple** distinct HQCC seeds — not merely
the homogenised ray through one seed — remains missing.

_Generated by gap_a_bj_families.py_