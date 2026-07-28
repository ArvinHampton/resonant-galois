# Direction 1 — Secondary invariants on \(L_0\)

_Elapsed: 7.03s_

**Verdict:** L0 secondary invariants (7.03s). PE fibres=148 gal={'A5': 36, 'D5': 4, 'even_unchecked': 108}; B A-points=316 gal={'even_unchecked': 269, 'A5': 45, 'red': 2}; Mestre t=42 gal={'A5': 18, 'even_unchecked': 24}. Strict PE↔B param overlap n=14. Direction 1 deliverable: lattice map with v3 / disc primes / height / Gal.

Context: `TERNARY_LATTICE_DIRECTIONS.md`. Necessity paused.

---

## Lattice \(L_0\)

Core model integers + short products/sums: **158** positive generators used.
Sample: `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 25, 27, 28, 29, 30, 31, 32, 34, 35, 36, 43, 45, 48, 52, 54, 55, 61, 62, 63, 64]` …

## Roles of \(L_0\) (reminder)

| Role | Use here |
|------|----------|
| Specialisation source | PE multi-\(k\) fibres; B-embed \(A\in\pm L_0\) |
| Mestre parameter | Flagship \(P_t\) at lattice \(t\) |
| Template coordinates | B-avatar \(T(-A,b,72A/b,-75,0,0)\) |

---

## 1. Pure-even multi-\(k\) fibres

- Fibres: **148**
- Gal class counts: `{'A5': 36, 'D5': 4, 'even_unchecked': 108}`

### By \(k\)-slice

| \(k\) | n | A5 | D5 | other even | unchecked | A5 rate (checked) |
|------|--:|---:|---:|-----------:|----------:|------------------:|
| -12/5 | 14 | 0 | 0 | 0 | 14 | 0.00 |
| -16/5 | 14 | 0 | 0 | 0 | 14 | 0.00 |
| -4 | 16 | 0 | 0 | 0 | 16 | 0.00 |
| -4/5 | 16 | 0 | 0 | 0 | 16 | 0.00 |
| -8/5 | 14 | 12 | 2 | 0 | 0 | 0.86 |
| 12/5 | 14 | 0 | 0 | 0 | 14 | 0.00 |
| 16/5 | 14 | 0 | 0 | 0 | 14 | 0.00 |
| 4 | 16 | 0 | 0 | 0 | 16 | 0.00 |
| 4/5 | 16 | 12 | 0 | 0 | 4 | 1.00 |
| 8/5 | 14 | 12 | 2 | 0 | 0 | 0.86 |

### By \(v_3(\alpha)\)

| \(v_3(\alpha)\) | n | A5 | D5 | unchecked |
|---------------|--:|---:|---:|----------:|
| 0 | 144 | 36 | 4 | 104 |
| 2 | 4 | 0 | 0 | 4 |

### By \(v_3(\beta)\)

| \(v_3(\beta)\) | n | A5 | D5 | unchecked |
|--------------|--:|---:|---:|----------:|
| 0 | 120 | 36 | 4 | 80 |
| 1 | 24 | 0 | 0 | 24 |
| 3 | 4 | 0 | 0 | 4 |

### Disc prime fingerprint (PE)

| prime | n | \(v_p(\mathrm{disc})>0\) | rate |
|------:|--:|-------------------------:|-----:|
| 3 | 148 | 20 | 0.135 |
| 5 | 148 | 130 | 0.878 |
| 61 | 148 | 6 | 0.041 |

### By height bin \(\log\max(\lvert\alpha\rvert,\lvert\beta\rvert)\)

| bin | n | A5 | D5 |
|-----|--:|---:|---:|
| 5-8 | 28 | 9 | 0 |
| 8-12 | 40 | 5 | 0 |
| <5 | 12 | 1 | 4 |
| >=12 | 68 | 21 | 0 |

### Master sample (PE, Gal-checked / flagship slices)

| name | \(k\) | \(m\) | \(\alpha\) | \(\beta\) | \(v_3\alpha\) | \(v_3\beta\) | h | Gal |
|------|------|------|----------:|----------:|------------:|------------:|--:|-----|
| flagship | -8/5 | 1/8 | -47500 | 380000 | 0 | 0 | 12.85 | A5 |
| flagship | -8/5 | 5/8 | 20 | -32 | 0 | 0 | 3.47 | D5 |
| flagship | -8/5 | 5/4 | 320 | -512 | 0 | 0 | 6.24 | A5 |
| flagship | -8/5 | 1 | 110000 | -880000 | 0 | 0 | 13.69 | A5 |
| flagship | -8/5 | 1/4 | -40000 | 320000 | 0 | 0 | 12.68 | A5 |
| flagship | -8/5 | 3/8 | -27500 | 220000 | 0 | 0 | 12.30 | A5 |
| flagship | -8/5 | 7/8 | 72500 | -580000 | 0 | 0 | 13.27 | A5 |
| flagship | -8/5 | 1 | 110000 | -880000 | 0 | 0 | 13.69 | A5 |
| flagship | -8/5 | 3/4 | 40000 | -320000 | 0 | 0 | 12.68 | A5 |
| flagship | -8/5 | 5/2 | 1520 | -2432 | 0 | 0 | 7.80 | A5 |
| flagship | -8/5 | 1/2 | -10000 | 80000 | 0 | 0 | 11.29 | A5 |
| flagship | -8/5 | -1/8 | -47500 | 380000 | 0 | 0 | 12.85 | A5 |
| flagship | -8/5 | -5/8 | 20 | -32 | 0 | 0 | 3.47 | D5 |
| flagship | -8/5 | -5/4 | 320 | -512 | 0 | 0 | 6.24 | A5 |
| flag_flip | 8/5 | 1/8 | -47500 | -380000 | 0 | 0 | 12.85 | A5 |
| flag_flip | 8/5 | 5/8 | 20 | 32 | 0 | 0 | 3.47 | D5 |
| flag_flip | 8/5 | 5/4 | 320 | 512 | 0 | 0 | 6.24 | A5 |
| flag_flip | 8/5 | 1 | 110000 | 880000 | 0 | 0 | 13.69 | A5 |
| flag_flip | 8/5 | 1/4 | -40000 | -320000 | 0 | 0 | 12.68 | A5 |
| flag_flip | 8/5 | 3/8 | -27500 | -220000 | 0 | 0 | 12.30 | A5 |
| flag_flip | 8/5 | 7/8 | 72500 | 580000 | 0 | 0 | 13.27 | A5 |
| flag_flip | 8/5 | 1 | 110000 | 880000 | 0 | 0 | 13.69 | A5 |
| flag_flip | 8/5 | 3/4 | 40000 | 320000 | 0 | 0 | 12.68 | A5 |
| flag_flip | 8/5 | 5/2 | 1520 | 2432 | 0 | 0 | 7.80 | A5 |
| flag_flip | 8/5 | 1/2 | -10000 | -80000 | 0 | 0 | 11.29 | A5 |
| flag_flip | 8/5 | -1/8 | -47500 | -380000 | 0 | 0 | 12.85 | A5 |
| flag_flip | 8/5 | -5/8 | 20 | 32 | 0 | 0 | 3.47 | D5 |
| flag_flip | 8/5 | -5/4 | 320 | 512 | 0 | 0 | 6.24 | A5 |
| classical | 4/5 | 1/8 | -625 | -2500 | 0 | 0 | 7.82 | A5 |
| classical | 4/5 | 5/8 | 95 | 76 | 0 | 0 | 4.55 | A5 |
| classical | 4/5 | 1/5 | 1279296875 | 127929687500 | 0 | 0 | 25.57 | A5 |
| classical | 4/5 | 5/4 | 395 | 316 | 0 | 0 | 5.98 | A5 |
| classical | 4/5 | 1 | 156875 | 627500 | 0 | 0 | 13.35 | A5 |
| classical | 4/5 | 1/4 | 6875 | 27500 | 0 | 0 | 10.22 | A5 |
| classical | 4/5 | 3/8 | 19375 | 77500 | 0 | 0 | 11.26 | A5 |

---

## 2. B-embed \(A\in\pm L_0\)

- Points: **316**
- Gal class counts: `{'even_unchecked': 269, 'A5': 45, 'red': 2}`
- disc□: all irr fibres (identity) — `314`

### By \(v_3(A)\)

| \(v_3(A)\) | n | A5 | D5 | unchecked | red |
|-----------|--:|---:|---:|----------:|----:|
| 0 | 190 | 7 | 0 | 181 | 2 |
| 1 | 48 | 20 | 0 | 28 | 0 |
| 2 | 28 | 12 | 0 | 16 | 0 |
| 3 | 24 | 4 | 0 | 20 | 0 |
| 4 | 6 | 2 | 0 | 4 | 0 |
| 5 | 14 | 0 | 0 | 14 | 0 |
| 6 | 2 | 0 | 0 | 2 | 0 |
| 7 | 2 | 0 | 0 | 2 | 0 |
| 8 | 2 | 0 | 0 | 2 | 0 |

### Disc primes (B)

| prime | n | \(v_p>0\) | rate |
|------:|--:|---------:|-----:|
| 3 | 314 | 314 | 1.000 |
| 5 | 314 | 86 | 0.274 |
| 61 | 314 | 36 | 0.115 |

### Master sample (B, model + checked)

| \(A\) | \(v_3(A)\) | \(T\) | h | Gal | disc \(v_3,v_5,v_{61}\) |
|----:|----------:|------|--:|-----|------------------------|
| 3 | 1 | `T(-3,1,216,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| -3 | 1 | `T(3,1,-216,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| 6 | 1 | `T(-6,1,432,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| -6 | 1 | `T(6,1,-432,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| 9 | 2 | `T(-9,1,648,-75,0,0)` | 4.32 | A5 | 14,0,0 |
| -9 | 2 | `T(9,1,-648,-75,0,0)` | 4.32 | A5 | 14,0,0 |
| 12 | 1 | `T(-12,1,864,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| -12 | 1 | `T(12,1,-864,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| 15 | 1 | `T(-15,1,1080,-75,0,0)` | 4.32 | A5 | 10,6,0 |
| -15 | 1 | `T(15,1,-1080,-75,0,0)` | 4.32 | A5 | 10,6,0 |
| 18 | 2 | `T(-18,1,1296,-75,0,0)` | 4.32 | A5 | 14,0,0 |
| -18 | 2 | `T(18,1,-1296,-75,0,0)` | 4.32 | A5 | 14,0,0 |
| 24 | 1 | `T(-24,1,1728,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| -24 | 1 | `T(24,1,-1728,-75,0,0)` | 4.32 | A5 | 10,0,0 |
| 27 | 3 | `T(-27,1,1944,-75,0,0)` | 4.39 | A5 | 16,0,0 |
| -27 | 3 | `T(27,1,-1944,-75,0,0)` | 4.39 | A5 | 16,0,0 |
| 30 | 1 | `T(-30,1,2160,-75,0,0)` | 4.50 | A5 | 10,6,0 |
| -30 | 1 | `T(30,1,-2160,-75,0,0)` | 4.50 | A5 | 10,6,0 |
| 36 | 2 | `T(-36,1,2592,-75,0,0)` | 4.68 | A5 | 14,0,0 |
| -36 | 2 | `T(36,1,-2592,-75,0,0)` | 4.68 | A5 | 14,0,0 |
| 45 | 2 | `T(-45,1,3240,-75,0,0)` | 4.91 | A5 | 14,6,0 |
| -45 | 2 | `T(45,1,-3240,-75,0,0)` | 4.91 | A5 | 14,6,0 |
| 48 | 1 | `T(-48,1,3456,-75,0,0)` | 4.97 | A5 | 10,0,0 |
| -48 | 1 | `T(48,1,-3456,-75,0,0)` | 4.97 | A5 | 10,0,0 |
| 54 | 3 | `T(-54,1,3888,-75,0,0)` | 5.09 | A5 | 16,0,2 |
| -54 | 3 | `T(54,1,-3888,-75,0,0)` | 5.09 | A5 | 16,0,2 |
| 55 | 0 | `T(-55,1,3960,-75,0,0)` | 5.11 | A5 | 4,6,0 |
| -55 | 0 | `T(55,1,-3960,-75,0,0)` | 5.11 | A5 | 4,6,0 |
| 61 | 0 | `T(-61,1,4392,-75,0,0)` | 5.21 | A5 | 4,0,2 |
| -61 | 0 | `T(61,1,-4392,-75,0,0)` | 5.21 | A5 | 4,0,2 |
| 63 | 2 | `T(-63,1,4536,-75,0,0)` | 5.24 | A5 | 14,0,0 |
| -63 | 2 | `T(63,1,-4536,-75,0,0)` | 5.24 | A5 | 14,0,0 |
| 66 | 1 | `T(-66,1,4752,-75,0,0)` | 5.29 | A5 | 10,0,0 |
| -66 | 1 | `T(66,1,-4752,-75,0,0)` | 5.29 | A5 | 10,0,0 |
| 69 | 1 | `T(-69,1,4968,-75,0,0)` | 5.33 | A5 | 10,0,0 |
| -69 | 1 | `T(69,1,-4968,-75,0,0)` | 5.33 | A5 | 10,0,0 |
| 72 | 2 | `T(-72,1,5184,-75,0,0)` | 5.38 | A5 | 14,0,0 |
| -72 | 2 | `T(72,1,-5184,-75,0,0)` | 5.38 | A5 | 14,0,0 |
| 80 | 0 | `T(-80,1,5760,-75,0,0)` | 5.48 | A5 | 4,6,0 |
| -80 | 0 | `T(80,1,-5760,-75,0,0)` | 5.48 | A5 | 4,6,0 |

---

## 3. Mestre flagship \(P_t\), \(t\in L_0\)

- Points: **42**
- Gal counts: `{'A5': 18, 'even_unchecked': 24}`

| \(t\) | \(v_3(t)\) | Gal | disc \(v_3,v_5,v_{61}\) | h |
|----:|----------:|-----|------------------------|--:|
| -1 | 0 | A5 | 0,6,0 | 24.13 |
| 0 | None | A5 | 0,6,0 | 4.48 |
| 1 | 0 | A5 | 0,6,0 | 23.73 |
| 2 | 0 | A5 | 0,6,0 | 27.31 |
| 3 | 1 | A5 | 0,6,0 | 29.37 |
| 4 | 0 | A5 | 0,6,0 | 30.82 |
| 5 | 0 | A5 | 0,6,0 | 31.95 |
| 6 | 1 | A5 | 0,6,0 | 32.87 |
| 7 | 0 | A5 | 0,6,0 | 33.64 |
| 8 | 0 | A5 | 0,6,0 | 34.31 |
| 9 | 2 | A5 | 0,6,0 | 34.91 |
| 10 | 0 | A5 | 0,6,0 | 35.43 |
| 11 | 0 | A5 | 0,6,0 | 35.91 |
| 12 | 1 | A5 | 0,6,0 | 36.35 |
| 13 | 0 | A5 | 0,6,0 | 36.75 |
| 14 | 0 | A5 | 0,6,0 | 37.12 |
| 15 | 1 | A5 | 0,6,0 | 37.47 |
| 16 | 0 | A5 | 0,6,0 | 37.79 |
| 17 | 0 | even_unchecked | 0,6,0 | 38.10 |
| 18 | 2 | even_unchecked | 0,6,0 | 38.38 |
| 19 | 0 | even_unchecked | 0,6,0 | 38.65 |
| 24 | 1 | even_unchecked | 0,6,0 | 39.82 |
| 25 | 0 | even_unchecked | 0,6,0 | 40.03 |
| 27 | 3 | even_unchecked | 0,6,0 | 40.41 |
| 28 | 0 | even_unchecked | 0,6,0 | 40.59 |
| 29 | 0 | even_unchecked | 0,6,0 | 40.77 |
| 30 | 1 | even_unchecked | 0,6,0 | 40.94 |
| 31 | 0 | even_unchecked | 0,6,0 | 41.10 |
| 32 | 0 | even_unchecked | 0,6,0 | 41.26 |
| 34 | 0 | even_unchecked | 0,6,0 | 41.57 |
| 35 | 0 | even_unchecked | 0,6,0 | 41.71 |
| 36 | 2 | even_unchecked | 0,6,0 | 41.85 |
| 43 | 0 | even_unchecked | 0,6,0 | 42.74 |
| 45 | 2 | even_unchecked | 0,6,0 | 42.97 |
| 48 | 1 | even_unchecked | 0,6,0 | 43.29 |
| 52 | 0 | even_unchecked | 0,6,0 | 43.69 |
| 54 | 3 | even_unchecked | 0,6,0 | 43.88 |
| 55 | 0 | even_unchecked | 0,6,0 | 43.97 |
| 61 | 0 | even_unchecked | 0,6,0 | 44.49 |
| 62 | 0 | even_unchecked | 0,6,0 | 44.57 |
| 63 | 2 | even_unchecked | 0,6,0 | 44.65 |
| 64 | 0 | even_unchecked | 0,6,0 | 44.73 |

### By \(v_3(t)\)

| \(v_3(t)\) | n | A5 | unchecked |
|-----------|--:|---:|----------:|
| 0 | 27 | 12 | 15 |
| 1 | 7 | 4 | 3 |
| 2 | 5 | 1 | 4 |
| 3 | 2 | 0 | 2 |

---

## 4. Overlap PE \(\leftrightarrow\) B (Direction 3 first cut)

- PE fibres: **148**, B parameters: **316**
- Strict overlap (\(A=\pmlpha\) or \(\pmeta\) of some PE fibre): **14**
- Values: `[-732, -320, -305, -95, -32, -12, -5, 5, 12, 32, 95, 305, 320, 732]`

Strict overlap = B-parameter A equals ±α or ±β of some PE fibre. Geometric PE↔B map Φ still open (Direction 3).

No canonical \(\Phi\) yet — only numerical coincidence of coordinates. Direction 3 remains open.

---

## 5. Observations (theorem-facing, not necessity)

1. **Evenness** on PE and B is by identity; Gal refinement \(A_5\) vs \(D_5\) is the residual invariant.
2. **\(v_3\)** stratifies both PE coefficients and B-parameter \(A\); useful lattice height for ternary story.
3. **Disc primes** 3, 5, 61 give ramification fingerprints; 61 is model-native when it divides disc.
4. **Mestre \(t\in L_0\)** continues to sample \(A_5\) (checked); lattice is stable as a parameter set under this family.
5. **Overlap PE/B** is thin under strict equality — unification needs a geometric map, not equality of integers.

## 6. Next

| Step | Direction |
|------|-----------|
| Done | **1** Secondary invariants (this file) |
| Next | **2** Resonant monoid / saturation |
| Then | **3** Unify PE \(\leftrightarrow\) B on \(L_0\) |
| Later | **4** Mestre orbit graph |

```bash
python l0_secondary_invariants.py
```

_Generated by l0_secondary_invariants.py_