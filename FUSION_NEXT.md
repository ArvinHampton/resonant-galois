# Fusion next — BJ pencil (Gap A) & natural functor (Gap B)

_Elapsed: 2.19s_

Context: mild surgery on \(\varphi\) cannot produce BJ seeds (`FUSION_DEPTH.md`).
This module attacks the two productive directions.

---

## Gap A — BJ geometric pencil

**Verdict:** Gap A partial success: the linear pencil through any two HQCC seeds is a geometric family in Q(t)[x] of BJ type containing those seeds by construction. It is NOT a Belyi pull-back of φ; monodromy of the full t-cover is not pure A5 (gal hist sample: {'odd_or_neg': 45, 'S5TransitiveSubgroups.A5': 2}). No linear pencil has disc identically a square polynomial (identical-square list: []). Fusion at equation level for seeds: YES via pencils; fusion as pure geometric A5 cover specialising only to even fibres: still open.

- Structural note: Linear BJ pencils through two seeds always recover those two seeds at t=0,1 (equation-level inclusion). Geometric monodromy of the family over Q(t) is A5 only if the generic fibre has Gal A5 (often true when disc is square on a Zariski-open set of the A5 locus — here disc is NOT identically square, so the cover of the t-line has geometric monodromy in S5 with A5 on the even locus).
- Disc identically square on a linear pencil: **none**
- Flagship–classical fibre Gal histogram: `{'odd_or_neg': 45, 'S5TransitiveSubgroups.A5': 2}`

### Pencils

#### `flagship__classical` (linear_coeff_pencil)
- α(t)=`75*t - 55`, β(t)=`88 - 72*t`
- sq-disc specialisations: 2, A5: 2
- seed hits: [{'t': 0, 'alpha': -55, 'beta': 88, 'disc': 58564000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 55*x + 88', 'seed_tag': 'flagship'}, {'t': 1, 'alpha': 20, 'beta': 16, 'disc': 1024000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 20*x + 16', 'seed_tag': 'classical'}]
- endpoints: `{'t=0': (-55, 88), 't=1': (20, 16)}`
  - t=0: α=-55 β=88 HIT_A5
  - t=1: α=20 β=16 HIT_A5

#### `flagship__hqcc_95_76` (linear_coeff_pencil)
- α(t)=`150*t - 55`, β(t)=`88 - 12*t`
- sq-disc specialisations: 2, A5: 2
- seed hits: [{'t': 0, 'alpha': -55, 'beta': 88, 'disc': 58564000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 55*x + 88', 'seed_tag': 'flagship'}, {'t': 1, 'alpha': 95, 'beta': 76, 'disc': 2085136000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 95*x + 76', 'seed_tag': 'hqcc'}]
- endpoints: `{'t=0': (-55, 88), 't=1': (95, 76)}`
  - t=0: α=-55 β=88 HIT_A5
  - t=1: α=95 β=76 HIT_A5

#### `flagship__period` (linear_coeff_pencil)
- α(t)=`150*t - 55`, β(t)=`444*t + 88`
- sq-disc specialisations: 2, A5: 2
- seed hits: [{'t': 0, 'alpha': -55, 'beta': 88, 'disc': 58564000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 55*x + 88', 'seed_tag': 'flagship'}, {'t': 1, 'alpha': 95, 'beta': 532, 'disc': 252301456000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 95*x + 532', 'seed_tag': 'period'}]
- endpoints: `{'t=0': (-55, 88), 't=1': (95, 532)}`
  - t=0: α=-55 β=88 HIT_A5
  - t=1: α=95 β=532 HIT_A5

#### `classical__hqcc_95_76` (linear_coeff_pencil)
- α(t)=`75*t + 20`, β(t)=`60*t + 16`
- sq-disc specialisations: 4, A5: 4
- seed hits: [{'t': 0, 'alpha': 20, 'beta': 16, 'disc': 1024000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 20*x + 16', 'seed_tag': 'classical'}, {'t': 1, 'alpha': 95, 'beta': 76, 'disc': 2085136000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 95*x + 76', 'seed_tag': 'hqcc'}]
- endpoints: `{'t=0': (20, 16), 't=1': (95, 76)}`
  - t=0: α=20 β=16 HIT_A5
  - t=1: α=95 β=76 HIT_A5
  - t=5: α=395 β=316 HIT_A5
  - t=16: α=1220 β=976 HIT_A5

#### `flagship__hqcc_m100` (linear_coeff_pencil)
- α(t)=`-45*t - 55`, β(t)=`312*t + 88`
- sq-disc specialisations: 2, A5: 2
- seed hits: [{'t': 0, 'alpha': -55, 'beta': 88, 'disc': 58564000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 55*x + 88', 'seed_tag': 'flagship'}, {'t': 1, 'alpha': -100, 'beta': 400, 'disc': 77440000000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 100*x + 400', 'seed_tag': 'hqcc'}]
- endpoints: `{'t=0': (-55, 88), 't=1': (-100, 400)}`
  - t=0: α=-55 β=88 HIT_A5
  - t=1: α=-100 β=400 HIT_A5

#### `flagship__hqcc_124` (linear_coeff_pencil)
- α(t)=`179*t - 55`, β(t)=`408*t + 88`
- sq-disc specialisations: 2, A5: 2
- seed hits: [{'t': 0, 'alpha': -55, 'beta': 88, 'disc': 58564000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 - 55*x + 88', 'seed_tag': 'flagship'}, {'t': 1, 'alpha': 124, 'beta': 496, 'disc': 196642060959744, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 124*x + 496', 'seed_tag': 'hqcc'}]
- endpoints: `{'t=0': (-55, 88), 't=1': (124, 496)}`
  - t=0: α=-55 β=88 HIT_A5
  - t=1: α=124 β=496 HIT_A5

#### `W_flagship__classical` (weighted_homo_pencil)
- α(t)=`20*t**4 - 55`, β(t)=`16*t**5 + 88`
- sq-disc specialisations: 0, A5: 0
- seed hits: [{'t': 'homo_flagship_1', 'alpha': -55, 'beta': 88, 'seed_tag': 'homo'}, {'t': 'homo_classical_1', 'alpha': 20, 'beta': 16, 'seed_tag': 'homo'}]

#### `W_flagship__hqcc_95_76` (weighted_homo_pencil)
- α(t)=`95*t**4 - 55`, β(t)=`76*t**5 + 88`
- sq-disc specialisations: 0, A5: 0
- seed hits: [{'t': 'homo_flagship_1', 'alpha': -55, 'beta': 88, 'seed_tag': 'homo'}, {'t': 'homo_hqcc_95_76_1', 'alpha': 95, 'beta': 76, 'seed_tag': 'homo'}]

#### `W_classical__zero` (weighted_homo_pencil)
- α(t)=`20`, β(t)=`t**5 + 16`
- sq-disc specialisations: 1, A5: 1
- seed hits: [{'t': -2, 'alpha': 20, 'beta': -16, 'disc': 1024000000, 'disc_sq': True, 'gal': 'S5TransitiveSubgroups.A5', 'status': 'HIT_A5', 'irr': True, 'poly': 'x**5 + 20*x - 16', 'seed_tag': 'classical'}, {'t': 'homo_classical_1', 'alpha': 20, 'beta': 16, 'seed_tag': 'homo'}, {'t': 'homo_zero_1', 'alpha': 0, 'beta': 1, 'seed_tag': 'homo'}]
  - t=-2: α=20 β=-16 HIT_A5

### Gap A — what this gives for fusion

1. **Equation-level inclusion of seeds:** any linear pencil through two seeds contains them at t=0,1.
2. **Geometric family:** \(f_t\in\mathbb{Q}(t)[x]\) of BJ type (not a twist of \(\varphi\)).
3. **Not a pure A5 Belyi specialisation:** disc not identically square; many fibres odd (\(S_5\)).
4. Homogenised rays through a single seed remain the **proved-even** theorem-grade families.

---

## Gap B — natural T₃ → braid

**Verdict:** Best uniqueness rule without leaving the 3-point cover: assign residue 2 ↦ σ∞=(σ0 σ1)^{-1} (Rule U4). It is unique given the standard π1 presentation and matches the period / third-sector role in 9 Maths. Path histograms differ from σ0σ1 (encoding dependence of the old scaffold). Fully natural 1–1 residue→branch requires a 4-point cover (U5); Step 1 shows such signatures generate A5 but are not absolutely rigid.

**Recommended rule:** `U4_residue2_to_sigma_inf`

### Uniqueness rules

```
{
  "U1_shortest_3A_not_pure": {
    "candidates": [
      "00",
      "11"
    ],
    "unique": false,
    "note": "May still have several length-min words"
  },
  "U2_length2_3A": {
    "candidates": [
      "0i",
      "1i",
      "00",
      "11"
    ],
    "unique": false
  },
  "U3_class_only": {
    "assignment": "residue 2 \u21a6 conjugacy class 3A (not a specific word)",
    "natural": true,
    "loses": "path-ordering / actual braid lift",
    "note": "Natural as a map to conjugacy classes in A5; not a functor to braids"
  },
  "U4_residue2_to_sigma_inf": {
    "assignment": "2 \u21a6 \u03c3\u221e = (\u03c30 \u03c31)^{-1}",
    "cycle_type": [
      5
    ],
    "motivation": "Third generator of \u03c01; 9 Maths period sector / second T3 branch",
    "unique": true,
    "natural_wrt_presentation": true,
    "note": "Uses the relation \u03c30 \u03c31 \u03c3\u221e=1; no auxiliary word choice among many 3A elements"
  },
  "U5_four_point_cover": {
    "signature": "(3A,3A,3A,2A) or (3A,3A,3A,3A)",
    "from_step1": "generates A5, not absolutely rigid (positive-dim Hurwitz)",
    "residue_map": "0\u2192C1, 1\u2192C2, 2\u2192C3, period\u2192C4",
    "natural": true,
    "cost": "Lose single rigid Belyi \u03c6; gain dedicated branch for each T3 residue",
    "status": "existence of Q-cover not constructed here; combinatorial monodromy available"
  }
}
```

### Encoding comparison (path histograms)

```
{
  "sigma0sigma1": {
    "(5,)": 8,
    "(3, 1, 1)": 6
  },
  "sigma_inf_U4": {
    "(3, 1, 1)": 10,
    "(5,)": 3,
    "(2, 2, 1)": 1
  },
  "sigma1sigma0": {
    "(5,)": 10,
    "(3, 1, 1)": 4
  }
}
```

### Four-point cover route

```
{
  "signatures": [
    "(3A,3A,3A,2A)",
    "(3A,3A,3A,3A)",
    "(2A,3A,3A,3A)"
  ],
  "step1_result": "generate A5; conjugacy orbits \u226b 1 (not absolutely rigid)",
  "implication": "A 4-point cover with dedicated branch per T3 residue can exist as a positive-dimensional Hurwitz family over Q-bar; rational points may give Q-covers. Not a single rigid Belyi map like \u03c6.",
  "construction_status": "combinatorial existence (Step 1); explicit equation not built"
}
```

### Gap B — what this gives for fusion

1. **U4 (recommended on 3-point cover):** residue \(2\mapsto\sigma_\infty=(\sigma_0\sigma_1)^{-1}\)
   is unique given the standard \(\pi_1\) presentation — removes arbitrary word search.
2. **U5 (4-point):** fully natural residue↔branch dictionary; costs absolute rigidity of \(\varphi\).
3. Class-only maps (U3) are natural but land in conjugacy classes, not braids.

---

## Combined status

| Gap | Progress | Still open |
|-----|----------|------------|
| A BJ pencil | Pencils through seeds give equation-level family containing seeds | Pure geometric A5 family with only even fibres / Hilbert to seeds only |
| B Natural functor | **U4** uniqueness rule on 3-point cover; U5 4-point route stated | Prove naturality formally; build explicit 4-point Q-cover |

Arithmetic foundations + rigid \(\varphi\) remain solid. Fusion advances on both tracks without mild φ-surgery.

_Generated by fusion_next.py_