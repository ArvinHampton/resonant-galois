# Stage D — Density and asymptotic arithmetic

_Elapsed: 67.5s_

**Verdict:** Stage D density/asymptotics (67.5s). D1=PASS: min_irr_rate=1.0, even_fail_total=0. D2=PASS: fibres=100, prime_specs=4500. D3=PASS: disc-height theorem + formula match. Stage D: COMPLETE.

Extends Stage B empirical rates into **stated conjectures (D1–D2)** and a **proved height theorem (D3)**, with machine-checkable tables.

Script: `stage_d_density.py` · Data: `build/STAGE_D_DATA.json`

---

## Scorecard

| block | pass | note |
|-------|:----:|------|
| D1 irreducibility density | **True** | conjecture + tables; min irr rate 1.0 |
| D2 Chebotarev proxy | **True** | Frobenius cycle-type histograms |
| D3 disc-height growth | **True** | **proved** asymptotic from pure-even identity |
| **Stage D** | **True** | |

---

## D1 — Irreducibility density (conjecture + evidence)

**Status:** `conjecture_with_evidence`

**Statement.** Conjecture D1 (irreducibility density on pure-even k-slices). For each fixed multi-seed ratio k=β/α in the HQCC pure-even catalogue, let L_k be the set of rational m such that α(m)=256m²−3125k⁴/256 and β(m)=k·α(m) lie in Z, with α(m)≠0. Order L_k by height H(m)=max(|num|,den) in lowest terms (or by |m|). Then the natural density of m∈L_k for which f_m=x⁵+α(m)x+β(m) is irreducible over Q is positive. Moreover, among those irreducible fibres, Gal=A5 for a positive-density subset (equivalently: Frobenius type (3,1,1) appears).

**Proved support.** Evenness is not conjectural: disc(f_m)=(256 α(m)² m)² identically, so every irr fibre has Gal ≤ A5. Irreducibility and A5 density remain analytic/number-theoretic (Hilbert irreducibility applies to the 2-param envelope; fixed-k slices are 1-param specialisations).

**Hilbert remark.** By Hilbert irreducibility, the 2-parameter pure-even envelope over Q(m,s) has a Zariski-dense set of specialisations with Gal=A5 (when geometric monodromy is A5 on a fibre, or when operational criteria hold). D1 is the thinner 1-param statement along fixed-k rays.

### Evidence tables (integer $Z$-coefficient fibres)

| $k$ | name | #Z pts | irr | red | even fail | irr rate | A5 / checked |
|------|------|-------:|----:|----:|----------:|---------:|-------------:|
| -4 | LSW | 840 | 840 | 0 | 0 | 1.0 | 29/30 |
| 4 | LSW_flip | 840 | 840 | 0 | 0 | 1.0 | 29/30 |
| -8/5 | flagship | 168 | 168 | 0 | 0 | 1.0 | 30/30 |
| 8/5 | flagship_flip | 168 | 168 | 0 | 0 | 1.0 | 30/30 |
| 4/5 | classical | 168 | 168 | 0 | 0 | 1.0 | 30/30 |
| -4/5 | classical_flip | 168 | 168 | 0 | 0 | 1.0 | 30/30 |
| -12/5 | s12 | 168 | 168 | 0 | 0 | 1.0 | 28/30 |
| 12/5 | s12_flip | 168 | 168 | 0 | 0 | 1.0 | 28/30 |
| -16/5 | s16 | 168 | 168 | 0 | 0 | 1.0 | 28/30 |
| 16/5 | s16_flip | 168 | 168 | 0 | 0 | 1.0 | 28/30 |

### Density profile by $|m|$ cap (selected)

**k=-4** (`LSW`): `{'10': {'n_int': 80, 'n_irr': 80, 'irr_rate': 1.0}, '20': {'n_int': 200, 'n_irr': 200, 'irr_rate': 1.0}, '40': {'n_int': 480, 'n_irr': 480, 'irr_rate': 1.0}, '80': {'n_int': 840, 'n_irr': 840, 'irr_rate': 1.0}}`

**k=-8/5** (`flagship`): `{'10': {'n_int': 16, 'n_irr': 16, 'irr_rate': 1.0}, '20': {'n_int': 40, 'n_irr': 40, 'irr_rate': 1.0}, '40': {'n_int': 96, 'n_irr': 96, 'irr_rate': 1.0}, '80': {'n_int': 168, 'n_irr': 168, 'irr_rate': 1.0}}`

**k=4/5** (`classical`): `{'10': {'n_int': 16, 'n_irr': 16, 'irr_rate': 1.0}, '20': {'n_int': 40, 'n_irr': 40, 'irr_rate': 1.0}, '40': {'n_int': 96, 'n_irr': 96, 'irr_rate': 1.0}, '80': {'n_int': 168, 'n_irr': 168, 'irr_rate': 1.0}}`

### Sample $A_5$ fibres

- k=-4: `[{'m': '-1/4', 'a': -3109, 'b': 12436, 'status': 'HIT_A5'}, {'m': '1/4', 'a': -3109, 'b': 12436, 'status': 'HIT_A5'}, {'m': '-1/2', 'a': -3061, 'b': 12244, 'status': 'HIT_A5'}]`
- k=-8/5: `[{'m': '-5/4', 'a': 320, 'b': -512, 'status': 'HIT_A5'}, {'m': '5/4', 'a': 320, 'b': -512, 'status': 'HIT_A5'}, {'m': '-5/2', 'a': 1520, 'b': -2432, 'status': 'HIT_A5'}]`
- k=4/5: `[{'m': '-5/4', 'a': 395, 'b': 316, 'status': 'HIT_A5'}, {'m': '5/4', 'a': 395, 'b': 316, 'status': 'HIT_A5'}, {'m': '-5/2', 'a': 1595, 'b': 1276, 'status': 'HIT_A5'}]`

---

## D2 — Chebotarev / Frobenius types (conjecture + histograms)

**Status:** `empirical_chebotarev_proxy`

**Statement.** Conjecture D2 (Chebotarev along pure-even k-slices). Let f_m be an irreducible fibre on a multi-seed pure-even k-slice with Gal(f_m/Q)=A5. For unramified primes p, the Frobenius conjugacy class is equidistributed among conjugacy classes of A5 with natural densities equal to class sizes / |A5|. In particular, factorization type (3,1,1) occurs with density 20/60=1/3, type (5) with density 24/60=2/5, and type (2,2,1) with density 15/60=1/4.

**Proved support.** Conditional on Gal=A5, Chebotarev density theorem supplies the class densities. The programme already uses type (3,1,1) as the operational A5 witness. D2 asserts equidistribution empirically along the pure-even families (not a new group-theory theorem).

### Predicted $A_5$ class frequencies (Chebotarev)

`{'(1,1,1,1,1)': 0.016666666666666666, '(1,2,2)': 0.25, '(1,1,3)': 0.3333333333333333, '(5,)': 0.4}`

### Global observed factorization-type frequencies

(Aggregate unramified prime factorisations across sampled irr fibres.)

| pattern | count | freq |
|---------|------:|-----:|
| `(5,)` | 1916 | 0.425778 |
| `(1, 1, 3)` | 1376 | 0.305778 |
| `(1, 2, 2)` | 1148 | 0.255111 |
| `(1, 1, 1, 1, 1)` | 60 | 0.013333 |

- Fibres sampled: **100**
- A5 among Gal subsample: **52**
- Total prime specialisations: **4500**

### Per-$k$ top patterns (abbrev.)

- **k=-4** fibres=10: `[('(5,)', {'count': 190, 'freq': 0.422222}), ('(1, 1, 3)', {'count': 162, 'freq': 0.36}), ('(1, 2, 2)', {'count': 96, 'freq': 0.213333}), ('(1, 1, 1, 1, 1)', {'count': 2, 'freq': 0.004444})]`
- **k=4** fibres=10: `[('(5,)', {'count': 190, 'freq': 0.422222}), ('(1, 1, 3)', {'count': 162, 'freq': 0.36}), ('(1, 2, 2)', {'count': 96, 'freq': 0.213333}), ('(1, 1, 1, 1, 1)', {'count': 2, 'freq': 0.004444})]`
- **k=-8/5** fibres=10: `[('(5,)', {'count': 190, 'freq': 0.422222}), ('(1, 1, 3)', {'count': 134, 'freq': 0.297778}), ('(1, 2, 2)', {'count': 122, 'freq': 0.271111}), ('(1, 1, 1, 1, 1)', {'count': 4, 'freq': 0.008889})]`
- **k=8/5** fibres=10: `[('(5,)', {'count': 190, 'freq': 0.422222}), ('(1, 1, 3)', {'count': 134, 'freq': 0.297778}), ('(1, 2, 2)', {'count': 122, 'freq': 0.271111}), ('(1, 1, 1, 1, 1)', {'count': 4, 'freq': 0.008889})]`
- **k=4/5** fibres=10: `[('(5,)', {'count': 186, 'freq': 0.413333}), ('(1, 1, 3)', {'count': 152, 'freq': 0.337778}), ('(1, 2, 2)', {'count': 108, 'freq': 0.24}), ('(1, 1, 1, 1, 1)', {'count': 4, 'freq': 0.008889})]`
- **k=-4/5** fibres=10: `[('(5,)', {'count': 186, 'freq': 0.413333}), ('(1, 1, 3)', {'count': 152, 'freq': 0.337778}), ('(1, 2, 2)', {'count': 108, 'freq': 0.24}), ('(1, 1, 1, 1, 1)', {'count': 4, 'freq': 0.008889})]`

**Interpretation.** Observed masses concentrate on types compatible with subgroups of $A_5$ (e.g. $(1,1,3)$, $(5,)$, $(1,2,2)$). Exact match to class proportions $1/3,2/5,1/4$ is asymptotic in the prime; finite samples are a **proxy**, not a proof of equidistribution.

---

## D3 — Disc-height growth (**theorem**)

**Status:** `proved`

**Statement.** Theorem D3 (disc height on pure-even k-slices). For k∈Q\{0} and m∈Q\{0} with α(m)=256m²−3125k⁴/256 ≠ 0, the Bring–Jerrard fibre f_m=x⁵+α(m)x+β(m), β=kα, satisfies disc(f_m)=(256 α(m)² m)². Consequently log|disc(f_m)| = 2 log|256 α(m)² m|. As |m|→∞ with k fixed, α(m)∼256 m², hence √|disc| ∼ 2^{24} |m|^5 and |disc| ∼ 2^{48} |m|^{10}, i.e. log|disc(f_m)| = 10 log|m| + 48 log 2 + o(1).

- Leading $\sqrt{|\mathrm{disc}|}$ monomial in $m$: `16777216*m**5`
- $\deg_m |\mathrm{disc}|$ (leading): **10**
- $\log|\mathrm{disc}| \sim 10\,\log|m| + 48\log 2$

### Numerical check of $\sqrt{\mathrm{disc}}=|256\alpha^2 m|$

| $k$ | samples | formula match | mean residual $\log|\mathrm{disc}|-\mathrm{asymp}$ |
|------|--------:|:-------------:|----------------------------------------------------------:|
| -4 | 20 | 20/20 | 14.8749 |
| 4 | 20 | 20/20 | 14.8749 |
| -8/5 | 20 | 20/20 | -0.831584 |
| 8/5 | 20 | 20/20 | -0.831584 |
| 4/5 | 20 | 20/20 | -0.031554 |
| -4/5 | 20 | 20/20 | -0.031554 |
| -12/5 | 20 | 20/20 | -1.878777 |
| 12/5 | 20 | 20/20 | -1.878777 |

### Sample rows (\(k=-8/5\) flagship)

`[{'m': '-5/8', 'a': 20, 'b': -32, 'disc': 4096000000, 'log_disc': 22.133277, 'log_m': -0.470004, 'asymp_10logm_48log2': 28.571028, 'residual_log': -6.437752, 'sqrt_formula_match': True}, {'m': '5/8', 'a': 20, 'b': -32, 'disc': 4096000000, 'log_disc': 22.133277, 'log_m': -0.470004, 'asymp_10logm_48log2': 28.571028, 'residual_log': -6.437752, 'sqrt_formula_match': True}, {'m': '-5/4', 'a': 320, 'b': -512, 'disc': 1073741824000000, 'log_disc': 34.609926, 'log_m': 0.223144, 'asymp_10logm_48log2': 35.5025, 'residual_log': -0.892574, 'sqrt_formula_match': True}, {'m': '5/4', 'a': 320, 'b': -512, 'disc': 1073741824000000, 'log_disc': 34.609926, 'log_m': 0.223144, 'asymp_10logm_48log2': 35.5025, 'residual_log': -0.892574, 'sqrt_formula_match': True}, {'m': '-15/8', 'a': 820, 'b': -1312, 'disc': 104168853504000000, 'log_disc': 39.18479, 'log_m': 0.628609, 'asymp_10logm_48log2': 39.557151, 'residual_log': -0.372362, 'sqrt_formula_match': True}]`

---

## Outsider checklist (regenerate)

```bash
cd resonant_galois
python stage_d_density.py
```

Inspect `build/STAGE_D_DATA.json` for tables. No Resonant narrative required: only BJ disc formula, pure-even parametrisation, and standard irreducibility / Frobenius factorisation.

## Relation to Stages A–C, E+

| Stage | Link |
|-------|------|
| A1 | Pure-even identity is the **proved engine** for D3 and evenness in D1 |
| B1/B2 | Empirical rates upgraded to D1/D2 conjectures + larger tables |
| C | Structural criteria unchanged; density is arithmetic, not fusion |
| E/V | JSON tables are the reproducibility surface |
| J | D3 is citable; D1–D2 citable as conjectures with data |

## Success criterion (Stage D)

≥1 density statement as **theorem** or **conjecture with machine-checkable evidence** — **met**: D3 theorem + D1/D2 conjectures with regenerable tables.

**Stage D complete:** **True**

_Generated by stage_d_density.py — Resonant Number Theory Stage D._