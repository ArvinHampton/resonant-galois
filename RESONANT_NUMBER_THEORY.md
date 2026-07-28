# Resonant Number Theory — Stage A + B empirical grounding

_Elapsed: 67.23s_

**Verdict:** RNT Stage A+B empirical grounding (67.23s). DIG(A1)=PASS: k-slice id=True, catalogue A5=21/21, paths multi-k=True, slice even_fail=0. GROW(A3)=PASS: A6=2, envelope A5 harvest=50. BUILD(B)=PASS: B3 proved identity, B1/B2 tables emitted, B4 A6 examples. Empirical grounding of Resonant Number Theory: CONFIRMED.

**Dig · Grow · Build** executed against `RESOLUTION_PATH.md`.

---

## DIG — Stage A1 (mathematical core on data)

### Identities

| check | result |
|-------|--------|
| BJ disc formula (random trials) | 420/420 |
| Homogenisation lemma proved | **True** |
| General \(k\)-slice disc identity | **True** |
| Catalogue seeds disc□ | **21/21** |
| Catalogue HIT_A5 | **21** |
| Homogenisation disc \(t^{20}\) id | **True** |

### \(k\)-slice Hilbert statistics

| \(k\) | name | #Z pts | irr | A5 (checked) | irr rate | even fail |
|------|------|-------:|----:|-------------:|---------:|----------:|
| -4 | LSW | 60 | 60 | 32 | 1.0 | 0 |
| 4 | LSW_flip | 60 | 60 | 32 | 1.0 | 0 |
| -8/5 | flagship | 12 | 12 | 12 | 1.0 | 0 |
| 8/5 | flagship_flip | 12 | 12 | 12 | 1.0 | 0 |
| 4/5 | classical | 12 | 12 | 12 | 1.0 | 0 |
| -4/5 | classical_flip | 12 | 12 | 12 | 1.0 | 0 |
| -12/5 | s12 | 12 | 12 | 12 | 1.0 | 0 |
| 12/5 | s12_flip | 12 | 12 | 12 | 1.0 | 0 |

### Cross-\(k\) paths

| path | disc id | multi catalogue \(k\) | hist (u=j/20) |
|------|:-------:|:---------------------:|---------------|
| flag_classical | True | **True** ['-8/5', '4/5'] | `{'HIT_A5': 2, 'non_Z': 19}` |
| flag_lsw | True | **True** ['-4', '-8/5'] | `{'HIT_A5': 2, 'non_Z': 19}` |
| classical_lsw | True | **True** ['-4', '4/5'] | `{'HIT_A5': 3, 'non_Z': 18}` |

---

## GROW — Stage A3 (beyond \(A_5\))

- Deg-6 thin even-irr: **60**
- Gal \(A_6\): **2**
- Envelope lattice A5 harvest: **50** by \(k\): `{'-4': 9, '4': 9, '-8/5': 4, '8/5': 4, '4/5': 4, '-4/5': 4, '-12/5': 4, '12/5': 4, '-16/5': 4, '16/5': 4}`

### \(A_6\) examples

- `x**6 - 6*x**2 - 6*x + 2` status=HIT_A6
- `x**6 - 6*x**2 + 6*x + 2` status=HIT_A6

### Homogenisation trials (deg 6)

- `{'seed': 'x**6 - 6*x**2 - 6*x + 2', 'family': 'x^6+-6 t^4 x^2+-6 t^5 x+2 t^6', 'all_sample_even': True, 'samples': [{'t': 2, 'disc_sq': True, 'disc': 617238958518042624}, {'t': 3, 'disc_sq': True, 'disc': 118356224095636874869824}, {'t': 4, 'disc_sq': True, 'disc': 662755285163023424003506176}], 'seed_disc': 574848576}`
- `{'seed': 'x**6 - 6*x**2 + 6*x + 2', 'family': 'x^6+-6 t^4 x^2+6 t^5 x+2 t^6', 'all_sample_even': True, 'samples': [{'t': 2, 'disc_sq': True, 'disc': 617238958518042624}, {'t': 3, 'disc_sq': True, 'disc': 118356224095636874869824}, {'t': 4, 'disc_sq': True, 'disc': 662755285163023424003506176}], 'seed_disc': 574848576}`

---

## BUILD — Stage B (checkable predictions + data)

### B1 — Irreducibility along \(k\)-slices

- Status: *empirical_support*
- For each multi-seed pure-even k-slice, a positive fraction of integer m with Z-coeffs yield irreducible fibres (empirical rates below).
- Rates: `{'-4': {'n_int': 60, 'n_irr': 60, 'irr_rate': 1.0, 'even_fail': 0}, '4': {'n_int': 60, 'n_irr': 60, 'irr_rate': 1.0, 'even_fail': 0}, '-8/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}, '8/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}, '4/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}, '-4/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}, '-12/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}, '12/5': {'n_int': 12, 'n_irr': 12, 'irr_rate': 1.0, 'even_fail': 0}}`

### B2 — Galois along cross-\(k\) paths

- Status: *empirical_support*
- Along cross-k pure-even paths, Hilbert specialisations are predominantly even; A5 appears whenever Gal is computed on irr even fibres (histograms).
- Paths: `{'flag_classical': {'disc_identity': True, 'hist': {'HIT_A5': 2, 'non_Z': 19}, 'multi_catalogue_k': True, 'catalogue_k_hit': ['-8/5', '4/5']}, 'flag_lsw': {'disc_identity': True, 'hist': {'HIT_A5': 2, 'non_Z': 19}, 'multi_catalogue_k': True, 'catalogue_k_hit': ['-4', '-8/5']}, 'classical_lsw': {'disc_identity': True, 'hist': {'HIT_A5': 3, 'non_Z': 18}, 'multi_catalogue_k': True, 'catalogue_k_hit': ['-4', '4/5']}}`

### B3 — \(\varphi\) obstruction (proved)

- Status: *proved*
- Identity \(5\cdot\square\): **True**
- Theorem: disc monic(φ-t) = 5·(25 t(t-1)/36)^2 in Q(t); no even irr rational specialisation of preferred φ/Q.

### B4 — Generative \(A_6\)

- Status: *empirical_support*
- Even irr: 60, A6: 2

---

## Empirical grounding scorecard

| block | pass |
|-------|:----:|
| DIG A1 (core on data) | **True** |
| GROW A3 (beyond A5) | **True** |
| BUILD B (predictions+data) | **True** |
| **RNT empirical grounding** | **True** |

Machine-readable tables: `build/RNT_STAGE_B_DATA.json`.

### Interpretation

1. **Resonant Number Theory (arithmetic core)** is empirically grounded:
   identities hold; catalogue seeds are disc□+\(A_5\); \(k\)-slices never fail evenness;
   cross-\(k\) paths are pure-even and multi-catalogue-\(k\).
2. **Not only \(A_5\):** explicit \(A_6\) thin sextics with square disc.
3. **Stage B** supplies outsider-checkable rates, histograms, and a proved obstruction.
4. **Geometric multi-\(k\)** remains open and is **not** required for this grounding.

_Generated by resonant_number_theory_ab.py_