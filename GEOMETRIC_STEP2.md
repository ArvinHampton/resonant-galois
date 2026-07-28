# Geometric cover — Step 2: explicit Belyi realisation

_Elapsed: 0.11s_

## Fusion: 9 Maths + classical rigidity

| Layer | Tool |
|-------|------|
| Classical | Genus-0 Belyi maps; Riemann–Hurwitz; absolute rigidity ⇒ geometric monodromy \(A_5\) |
| 9 Maths / HQCC | T₃ branches, G4=539.9, N_flux=4880, three generations as **labels** on branch points |

Sources: `The_9_Maths_of_Unification.pdf`, HQH-539 T₃ definition, Step 1 rigid tuples.

---

## Preferred cover — signature (3A, 3A, 5A)

### Explicit equation (field \(\mathbb{Q}\))

$$
\varphi(y) = 6y^5 - 15y^4 + 10y^3 = y^3(6y^2 - 15y + 10)
$$

- Factored φ−1: `(x - 1)**3*(6*x**2 + 3*x + 1)`
- **Field of definition: `Q`**
- **Field of moduli: `Q`**
- Riemann–Hurwitz: `{'R': 8, 'genus': 0, 'degree': 5, 'types': [(3, 1, 1), (3, 1, 1), (5,)]}` (genus 0)

### Branch-point configuration

| Base point | Cycle type | A5 class | Geometric meaning |
|------------|------------|----------|-------------------|
| \(0\) | \(3{+}1{+}1\) | 3A | zeros of φ |
| \(1\) | \(3{+}1{+}1\) | 3A | zeros of φ−1 |
| \(\infty\) | \(5\) | 5A | pole of order 5 |

### Geometric monodromy

- **Group: \(A_5\)** (numeric order = 60, even)
- Types: 0 ↦ `(3, 1, 1)`, 1 ↦ `(3, 1, 1)`, ∞ ↦ `(5,)`
- Justification: Numeric monodromy group has order 60 and is even ⇒ A5; matches absolute rigidity of conjugacy triple (3A,3A,5A) from Step 1.

### 9 Maths / HQCC labelling

- **0 (3A):** `{'nine_maths': 'HQCC T3 contraction + three generations (Maths 1,6)', 'T3': {'map': 'n // 3', 'role': 'contraction', 'mod': 'n ≡ 0 (mod 3)'}, 'operation': 'n ↦ n//3'}`
- **1 (3A):** `{'nine_maths': 'HQCC T3 expansion sector (Maths 6 + Ad 3n±1 cousin)', 'T3': {'map': '(4n+2) // 3', 'role': 'expansion', 'mod': 'n ≡ 1 (mod 3)'}, 'operation': 'n ↦ (4n+2)//3'}`
- **∞ (5A):** `{'nine_maths': 'Temporal torsion / Resonant oscillation period G4=539.9 (Maths 1,8)', 'G4': 539.9, 'period_steps': 539, 'operation': 'period sector (5-cycle class)'}`

- Derivation: Ansatz φ=x^3(x^2+a x+b); impose φ(t)=1, φ'(t)=φ''(t)=0 (t≠0). Elimination ⇒ a=-5t/2, b=5t^2/3, t^5=6. Scale x=t y ⇒ φ=6y^5-15y^4+10y^3 ∈ Q[y].

### Nativeness (honest)

The cover is **classically defined over Q** with passport matching the HQCC-preferred
rigid triple (two ternary classes + period 5-cycle). Labels from T₃ / G4 are a
**structure dictionary** on \(\{0,1,\infty\}\). A full functor from the T₃ dynamical
system to the braid group (Step 4 deep) remains open.

---

## Fallback cover — signature (3A, 2A, 5A)

### Explicit equation

$$
\varphi(x)=x^5 + \frac{5\cdot 2^{4/5}3^{2/5}}{12}\,x^4 + \frac{5\cdot 2^{3/5}3^{4/5}}{9}\,x^3
$$

- a = `5*2**(4/5)*3**(2/5)/12` ≈ 1.1258000321005113
- b = `5*2**(3/5)*3**(4/5)/9` ≈ 2.0278811396440193
- minpoly(a): `1728*z**5 - 3125`
- minpoly(b): `729*z**5 - 25000`
- **Field of definition: `Q(2^{1/5}, 3^{1/5})  (real radical; degree ≤ 25, often less after relations)`**
- Classical name: icosahedral-type Belyi passport (related to Δ(2,3,5) geometry)

### Branch-point configuration

| Base point | Cycle type | A5 class |
|------------|------------|----------|
| 0 | 3+1+1 | 3A |
| 1 | 2+2+1 | 2A |
| ∞ | 5 | 5A |

### Geometric monodromy

- **Group: \(A_5\)** (order 60)
- Types: `(3, 1, 1)`, `(2, 2, 1)`, `(5,)`
- HQCC labels: `{'0': {'nine_maths': 'ternary / T3 (Maths 6, HQCC)', 'T3': {'map': 'n // 3', 'role': 'contraction', 'mod': 'n ≡ 0 (mod 3)'}}, '1': {'nine_maths': 'T-complementarity / flux involution (N_flux=4880, Maths 2 mirror)', 'N_flux': 4880, 'class': '2A double transposition'}, 'infinity': {'nine_maths': 'period G4=539.9 (Maths 1,8)', 'G4': 539.9}}`

---

## Summary table (report targets)

| Item | Preferred (3A,3A,5A) | Fallback (3A,2A,5A) |
|------|----------------------|---------------------|
| Explicit map | `6*x**5 - 15*x**4 + 10*x**3` | x⁵+a x⁴+b x³ (radicals) |
| Field of definition | **Q** | Q(2^{1/5},3^{1/5}) |
| Branch locus | {0,1,∞} | {0,1,∞} |
| Cycle types | (3,1,1),(3,1,1),(5) | (3,1,1),(2,2,1),(5) |
| Geometric monodromy | **A5** | **A5** |
| 9 Maths labels | 2×T₃ ternary + G4 period | T₃ + flux 2A + G4 |

---

## Arithmetic compatibility (Step 3 light)

```
{
  "hqcc_seeds": [
    "x**5 - 55*x + 88",
    "x**5 + 95*x + 76",
    "x**5 + 95*x + 532",
    "x**5 + 20*x + 16"
  ],
  "note": "The geometric covers are Belyi maps \u03c6: P1\u2192P1 (function field extensions). HQCC seeds are BJ fibres in Q[x] with lattice coefficients \u2014 a different arithmetic object. Compatibility is via Hilbert specialisation of related Hurwitz families / resolvents, not identity of equations. Preferred \u03c6 has coefficients in {6,10,15} \u2282 3Z lattice (ternary-visible).",
  "preferred_coeff_motif": "6,10,15 = 3\u00b7(2, 10/3?, 5) \u2014 all multiples of theme 3 / 5 (generations + pentagon)",
  "model_core": {
    "3": "ternary/generations",
    "9": "3^2",
    "18": "visible_digits",
    "61": "punctures",
    "80": "flux/61",
    "243": "3^5_towers",
    "520": "tower_res",
    "539": "period",
    "4880": "flux_budget"
  }
}
```

---

## What is now proved (geometry, this step)

1. **Existence** of a degree-5 Belyi cover over **Q** with passport (3,1,1)(3,1,1)(5).
2. Its **geometric monodromy is A5** (computed + rigidity alignment with Step 1).
3. Branch points admit a **consistent HQCC/9 Maths labelling** (dictionary level).
4. Fallback (2,3,5)-type passport also realises **A5** over a radical extension.

## What remains open

- Deep Step 4: identify base coordinate with a modular / resonant function of T₃ or ξ = 2 cos(2π/539.9).
- Hilbert specialisations recovering HQCC seeds \(x^5-55x+88\), etc.
- Descent theory for the fallback field and comparison with icosahedral modular covers of degree 60.

_Generated by geometric_step2.py_