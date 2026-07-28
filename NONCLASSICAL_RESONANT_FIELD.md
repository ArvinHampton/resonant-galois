# Non-classical possibilities — polynomials over \(\mathcal{R}\)

_Elapsed: 25.19s_

**Verdict:** Non-classical resonant field probes (25.19s). R_539 deg=210; proxies n=5,7,11,15. N1 samples ok; N2 split-Frob labels on flagship; N3 pure-even identity over generic fields=PASS. Scaffold LOCKED — research open, not classical closure.

Stays inside algebraic number theory but **leaves** the classical
setting \(f\in\mathbb{Z}[x]\), \(\mathrm{Gal}(f/\mathbb{Q})\).

## Working definition of \(\mathcal{R}\)

- **Period:** \(N=539\) (model G4≈539.9); \(539={'7': 2, '11': 1}\)
- **Full field:** \(\mathcal{R} = \mathbb{Q}(2\cos(2\pi/539))\), degree **210** (= \(\varphi(539)/2\))
- **Proxies (computable):** real subfields \(R_n=\mathbb{Q}(2\cos(2\pi/n))\) for \(n\in\{3,5,7,11,15,\ldots\}\) along the divisor tower

| \(n\) | \(\deg R_n\) | minpoly of \(2\cos(2\pi/n)\) | HQCC tag |
|------|------------:|-------------------------------|----------|
| 3 | 1 | `x + 1` | ternary/generations |
| 5 | 2 | `x**2 + x - 1` | — |
| 7 | 3 | `x**3 + x**2 - 2*x - 1` | period_factor |
| 11 | 5 | `x**5 + x**4 - 4*x**3 - 3*x**2 + 3*x + 1` | period_factor |
| 15 | 4 | `x**4 - x**3 - 4*x**2 + 4*x + 1` | — |
| 33 | 10 | `x**10 - x**9 - 10*x**8 + 10*x**7 + 34*x**6 - 34*x**5 - 43...` | — |
| 35 | 12 | `x**12 - x**11 - 12*x**10 + 11*x**9 + 54*x**8 - 43*x**7 - ...` | — |
| 55 | 20 | `_(deg too large / skip)_` | — |
| 539 | 210 | `_(deg too large / skip)_` | period |

**Stance.** Classical HQCC seeds remain in Z[x]. Non-classical work studies f∈R[x] and Gal(f/R), with Z-seeds as specialisations under embeddings R→R or traces/norms R→Q.

---

## N1 — Polynomials with coefficients in \(\mathcal{R}\)

Generate \(f\in R_n[x]\) (BJ shape \(x^5+\alpha x+\beta\), \(\alpha,\beta\in\mathbb{Z}[\xi]\))
and study Galois action **over the real subfield** \(R_n\), not only over \(\mathbb{Q}\).

- When \(\alpha,\beta\in\mathbb{Q}\), recover classical \(\mathrm{Gal}(f/\mathbb{Q})\).
- When \(\alpha,\beta\notin\mathbb{Q}\), \(\mathrm{Gal}(f/R_n)\) is the correct object;
  Weil restriction / norm forms give a polynomial over \(\mathbb{Q}\) of degree
  \(5\cdot[R_n:\mathbb{Q}]\) encoding the same arithmetic.

### \(n=5\) — \(\deg=2\)`

Minpoly: `x**2 + x - 1`

| tag | \(\alpha\) | \(\beta\) | in \(\mathbb{Q}\)? | classical Gal |
|-----|------|------|:----------------:|---------------|
| flagship_descended | `-55` | `88` | True | HIT_A5 |
| flag_plus_ternary_xi | `3*xi - 55` | `88 - xi` | False | — |
| classical_plus_xi | `xi + 20` | `16 - 3*xi` | False | — |
| LSW_plus_powers3_xi | `9*xi - 100` | `400 - 27*xi` | False | — |
| model_core_linear | `61*xi` | `80*xi + 3` | False | — |
| chebyshev_adjacent | `-xi - 4` | `2*xi` | False | — |

*Gal(f/R) for f∈R[x] is the decomposition group picture relative to R; when coeffs ∈ Q this recovers classical Gal(f/Q). Non-rational samples need relative resolvents / norm forms for a full Gal computation.*

### \(n=7\) — \(\deg=3\)`

Minpoly: `x**3 + x**2 - 2*x - 1`

| tag | \(\alpha\) | \(\beta\) | in \(\mathbb{Q}\)? | classical Gal |
|-----|------|------|:----------------:|---------------|
| flagship_descended | `-55` | `88` | True | HIT_A5 |
| flag_plus_ternary_xi | `3*xi - 55` | `88 - xi` | False | — |
| classical_plus_xi | `xi + 20` | `16 - 3*xi` | False | — |
| LSW_plus_powers3_xi | `9*xi - 100` | `400 - 27*xi` | False | — |
| model_core_linear | `61*xi` | `80*xi + 3` | False | — |
| chebyshev_adjacent | `xi**2 - 5` | `2*xi` | False | — |

*Gal(f/R) for f∈R[x] is the decomposition group picture relative to R; when coeffs ∈ Q this recovers classical Gal(f/Q). Non-rational samples need relative resolvents / norm forms for a full Gal computation.*

### \(n=11\) — \(\deg=5\)`

Minpoly: `x**5 + x**4 - 4*x**3 - 3*x**2 + 3*x + 1`

| tag | \(\alpha\) | \(\beta\) | in \(\mathbb{Q}\)? | classical Gal |
|-----|------|------|:----------------:|---------------|
| flagship_descended | `-55` | `88` | True | HIT_A5 |
| flag_plus_ternary_xi | `3*xi - 55` | `88 - xi` | False | — |
| classical_plus_xi | `xi + 20` | `16 - 3*xi` | False | — |
| LSW_plus_powers3_xi | `9*xi - 100` | `400 - 27*xi` | False | — |
| model_core_linear | `61*xi` | `80*xi + 3` | False | — |
| chebyshev_adjacent | `xi**2 - 5` | `2*xi` | False | — |

*Gal(f/R) for f∈R[x] is the decomposition group picture relative to R; when coeffs ∈ Q this recovers classical Gal(f/Q). Non-rational samples need relative resolvents / norm forms for a full Gal computation.*

### \(n=15\) — \(\deg=4\)`

Minpoly: `x**4 - x**3 - 4*x**2 + 4*x + 1`

| tag | \(\alpha\) | \(\beta\) | in \(\mathbb{Q}\)? | classical Gal |
|-----|------|------|:----------------:|---------------|
| flagship_descended | `-55` | `88` | True | HIT_A5 |
| flag_plus_ternary_xi | `3*xi - 55` | `88 - xi` | False | — |
| classical_plus_xi | `xi + 20` | `16 - 3*xi` | False | — |
| LSW_plus_powers3_xi | `9*xi - 100` | `400 - 27*xi` | False | — |
| model_core_linear | `61*xi` | `80*xi + 3` | False | — |
| chebyshev_adjacent | `xi**2 - 5` | `2*xi` | False | — |

*Gal(f/R) for f∈R[x] is the decomposition group picture relative to R; when coeffs ∈ Q this recovers classical Gal(f/Q). Non-rational samples need relative resolvents / norm forms for a full Gal computation.*

---

## N2 — Frobenius at primes split in the resonant field

Extra labelling of conjugacy classes: only use \(p\) that **split completely**
in \(R_n\) (minpoly of \(\xi_n\) splits into linears mod \(p\)).

This refines Stage D2 Chebotarev histograms: classes become pairs
\((\mathrm{cycle\ type},\ \mathrm{split\ in\ }R_n)\).

### \(n=5\) — split primes (p < 180): **17** found; sample `[11, 19, 29, 31, 41, 59, 61, 71, 79, 89, 101, 109]`
- Flagship at **split** primes: `{'primes_used': 16, 'patterns': {'(5,)': 7, '(1, 1, 3)': 5, '(1, 2, 2)': 4}}`
- Flagship at **control** primes: `{'primes_used': 20, 'patterns': {'(5,)': 8, '(1, 1, 3)': 7, '(1, 2, 2)': 4, '(1, 1, 1, 1, 1)': 1}}`

### \(n=7\) — split primes (p < 180): **11** found; sample `[13, 29, 41, 43, 71, 83, 97, 113, 127, 139, 167]`
- Flagship at **split** primes: `{'primes_used': 11, 'patterns': {'(1, 1, 3)': 4, '(5,)': 4, '(1, 2, 2)': 2, '(1, 1, 1, 1, 1)': 1}}`
- Flagship at **control** primes: `{'primes_used': 20, 'patterns': {'(5,)': 10, '(1, 2, 2)': 5, '(1, 1, 3)': 5}}`

### \(n=11\) — split primes (p < 180): **6** found; sample `[23, 43, 67, 89, 109, 131]`
- Flagship at **split** primes: `{'primes_used': 6, 'patterns': {'(1, 2, 2)': 3, '(1, 1, 3)': 2, '(5,)': 1}}`
- Flagship at **control** primes: `{'primes_used': 20, 'patterns': {'(5,)': 12, '(1, 1, 3)': 5, '(1, 2, 2)': 2, '(1, 1, 1, 1, 1)': 1}}`

### \(n=15\) — split primes (p < 180): **8** found; sample `[29, 31, 59, 61, 89, 149, 151, 179]`
- Flagship at **split** primes: `{'primes_used': 8, 'patterns': {'(5,)': 4, '(1, 1, 3)': 3, '(1, 2, 2)': 1}}`
- Flagship at **control** primes: `{'primes_used': 20, 'patterns': {'(5,)': 10, '(1, 2, 2)': 5, '(1, 1, 3)': 4, '(1, 1, 1, 1, 1)': 1}}`

**Programme use.** Conjugacy classes in Gal(f/Q) acquire an extra label: Frob_p for p split in R, vs inert/ramified. Classes that appear only at split primes are 'R-visible'; this refines Chebotarev histograms (Stage D2) without changing the group.

---

## N3 — Pure-even families over \(\mathcal{R}\); cosine constraints

**Field-agnostic pure-even identity:** **True**

Theorem (field-agnostic pure-even). Over any field F with char∉{2,5}, the formulae α=256m²−3125k⁴/256, β=kα give disc(x⁵+αx+β)=(256 α² m)² in F(m,k). Hence pure-even is not special to Q.

### Cosine-special \(k=2\cos(2\pi p/n)\)

Evenness holds automatically; the question is **which** such \(k\) arise
geometrically (monodromy constrained by multi-angle / Chebyshev relations).

| \(n\) | \(p\) | \(m\) | \(\alpha\) (abbrev) | in \(\mathbb{Q}\)? |
|------|------|------|---------------------|-------------------|
| 5 | 1 | 1 | `9375*sqrt(5)/512 + 109197/512` | False |
| 5 | 1 | 5/16 | `-9075/512 + 9375*sqrt(5)/512` | False |
| 5 | 1 | 5/4 | `9375*sqrt(5)/512 + 182925/512` | False |
| 5 | 2 | 1 | `109197/512 - 9375*sqrt(5)/512` | False |
| 5 | 2 | 5/16 | `-9375*sqrt(5)/512 - 9075/512` | False |
| 5 | 2 | 5/4 | `182925/512 - 9375*sqrt(5)/512` | False |
| 15 | 1 | 1 | `256 - 3125*(1 + sqrt(5) + sqrt(6)*sqrt(5 - sqrt...` | — |
| 15 | 1 | 5/16 | `25 - 3125*(1 + sqrt(5) + sqrt(6)*sqrt(5 - sqrt(...` | — |

### Model over \(R_5=\mathbb{Q}(\sqrt5)\)

- \(k=2cos(2π/5)=(-1+√5)/2\), \(m=5/16\)
- \(\alpha=-9075/512 + 9375*sqrt(5)/512\)
- \(\beta=27975/512 - 9225*sqrt(5)/512\)
- Model over Q(√5); recovers evenness over R_5. Links to K=Q(√5) base-change side route for φ (K_SQRT5_EVEN) but does not force descent to even-over-Q.

### Chebyshev / monodromy

If branch coordinates or ratio k are required to lie among {2 cos(2π a/n)}, Chebyshev recurrences T_r constrain which k-slices can appear as geometric specialisations of a cover defined over R_n. Pure-even arithmetic still holds for any k∈R; cosine selection is a geometric filter, not an evenness condition.

### Open questions (research)

- Q1. Which pure-even k∈R_n arise as cross-ratios of a cover defined over R_n?
- Q2. Does Gal(R_n/Q) act on the parameter space so that multi-k paths are unions of Galois orbits of cosine-special k?
- Q3. Can Nielsen labels be valued in conjugacy classes of Gal(f/R_n) with Frob at split primes as the arithmetic side of the dictionary?
- Q4. Is there f∈R_539[x] with monodromy A5 whose reduction at split primes recovers the HQCC Z-lattice seeds?

---

## What this is / is not

| Claim | Status |
|-------|--------|
| Pure-even identity over any field (char ≠2,5) | **Proved** (same algebra) |
| Classical Z-seed catalogues / Gal over Q | **Unchanged centre** |
| Rigid t=3 odd control over Q | **Locked** (`RIGID_FIBRE_T3.md`) |
| Full Gal(f/R_539) for non-rational f | **Open** (deg 210) |
| Cosine-forced multi-k geometric monodromy | **Open question** |
| Split-prime Frob labelling | **Operational** on proxies |
| Replacement of arithmetic multi-k over Q | **No** — enrichment only |

## Scorecard

| probe | pass |
|-------|:----:|
| R tower defined | **True** |
| N1 samples over R_n | **True** |
| N2 split primes + labels | **True** |
| N3 pure-even over generic + cosine trials | **True** |
| **Non-classical scaffold** | **True** |

```bash
python nonclassical_resonant_field.py
```

_Generated by nonclassical_resonant_field.py — Resonant Galois non-classical track._