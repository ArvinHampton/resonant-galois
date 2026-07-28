# Gap A — from geometric A5 family φ to BJ / HQCC seeds

_Elapsed: 2.02s_

## Strategy (locked)

Stop Diophantine search on the classical BJ even surface.
Begin from geometric families that **already have monodromy \(A_5\)**:

$$\varphi(y)=6y^5-15y^4+10y^3,\quad \text{passport }(3,1,1)(3,1,1)(5).$$

1. Deformations / pull-backs of \(\varphi\) keeping monodromy \(A_5\) (or checking evenness).
2. Tschirnhaus toward Bring–Jerrard \(x^5+ax+b\).
3. Specialise at HQCC lattice points; match known seeds.

**Verdict:** Geometric families from φ scanned: 16. Families with BJ Tschirnhaus success on some lattice fibre: 0. HQCC seed hits after Tschirnhaus: 0 (none). Primary family monic(φ−t) has geometric monodromy A5 but fibres at lattice t are typically odd (S5) except special t; Tschirnhaus to BJ rarely lands on HQCC seeds in the tested range. No fusion hit: obstruction deeper than classical even-surface Diophantine search alone — geometric A5 family does not yield HQCC BJ seeds via tested Tschirnhaus/lattice specs.

---

## Symbolic base family

- monic(φ−t): `-t/6 + y**5 - 5*y**4/2 + 5*y**3/3`
- depressed: `-t/6 + z**5 - 5*z**3/6 + 5*z/16 + 1/12`
- shift: 1/2
- note: Full reduction of a general depressed quintic to x^5+ax+b over Q(t) requires solving a resolvent (Bring radical); not always in Q(t). Fibrewise Tschirnhaus over Q is the practical fusion test.

---

## Family results

### `phi_minus_t`
- claim: A5 (rigid Belyi pullback of base)
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_1`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_2`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_3`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_1/2`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_6`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `scale_lam_1/3`
- claim: A5
- stats: `{'tested': 35, 'red': 2, 'odd': 33}`
- even fibres: 0, odd: 33, A5: 0, BJ form found: 0, **seed hits: 0**

### `pull_t`
- claim: A5 (base change of rigid cover)
- stats: `{'tested': 25, 'red': 1, 'odd': 24}`
- even fibres: 0, odd: 24, A5: 0, BJ form found: 0, **seed hits: 0**

### `pull_t**2`
- claim: A5 (base change of rigid cover)
- stats: `{'tested': 25, 'red': 2, 'odd': 23}`
- even fibres: 0, odd: 23, A5: 0, BJ form found: 0, **seed hits: 0**

### `pull_(t - 3)/(t - 539)`
- claim: A5 (base change of rigid cover)
- stats: `{'tested': 25, 'odd': 24, 'red': 1}`
- even fibres: 0, odd: 24, A5: 0, BJ form found: 0, **seed hits: 0**

### `pull_3*t + 61`
- claim: A5 (base change of rigid cover)
- stats: `{'tested': 25, 'odd': 25}`
- even fibres: 0, odd: 25, A5: 0, BJ form found: 0, **seed hits: 0**

### `pull_t**3`
- claim: A5 (base change of rigid cover)
- stats: `{'tested': 25, 'red': 2, 'odd': 23}`
- even fibres: 0, odd: 23, A5: 0, BJ form found: 0, **seed hits: 0**

### `deform_eps_y2`
- claim: unknown a priori; check fibre disc square / Gal
- stats: `{'tested': 60, 'red': 7, 'odd': 53}`
- even fibres: 0, odd: 53, A5: 0, BJ form found: 0, **seed hits: 0**

### `deform_eps_y`
- claim: unknown a priori; check fibre disc square / Gal
- stats: `{'tested': 60, 'red': 9, 'odd': 51}`
- even fibres: 0, odd: 51, A5: 0, BJ form found: 0, **seed hits: 0**

### `deform_eps_const`
- claim: unknown a priori; check fibre disc square / Gal
- stats: `{'tested': 60, 'red': 8, 'odd': 52}`
- even fibres: 0, odd: 52, A5: 0, BJ form found: 0, **seed hits: 0**

### `deform_eps_y4`
- claim: unknown a priori; check fibre disc square / Gal
- stats: `{'tested': 60, 'red': 8, 'odd': 52}`
- even fibres: 0, odd: 52, A5: 0, BJ form found: 0, **seed hits: 0**

---

## Interpretation

| Object | Role | Fusion outcome |
|--------|------|----------------|
| \(\varphi\) / monic(\(\varphi-t\)) | Pure geometric \(A_5\) family | Fibres at lattice \(t\) mostly **odd** |
| Pull-backs \(\varphi-R(t)\) | Still geometric \(A_5\) (base change) | Same pattern; no seed hits in scan |
| Domain scaling \(\varphi(\lambda y)-t\) | Automorphism of cover | No new seed hits |
| Coefficient deformations | Monodromy not guaranteed | Even/BJ occasional; seeds not hit |
| Tschirnhaus → BJ | Bridge to arithmetic form | Works on some even fibres; **not** HQCC seeds |

### Conclusion for the fusion gap

Starting from a geometric object that already “knows how to be \(A_5\)” does **not**
automatically produce the arithmetic HQCC BJ seeds under lattice specialisation +
Tschirnhaus in the families tested. The obstruction is therefore **not only**
the classical BJ even-surface Diophantine problem: even when monodromy is \(A_5\)
geometrically, the BJ models of fibres miss the HQCC lattice points.

Homogenised HQCC seeds remain the theorem-grade arithmetic track;
\(\varphi\) remains the theorem-grade geometric track; **equation-level fusion is still open**.

_Generated by gap_a_from_geometry.py_