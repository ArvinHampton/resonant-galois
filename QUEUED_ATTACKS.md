# Queued attacks — report

_Elapsed: 142.91s_

Order: (1) specialisation match → (2) resonant base → (3) T₃→braid functor.

---

## Attack 1 — Specialisation match

**Verdict:** Exact affine fibre match count = 0. Seeds and geometric cover share monodromy type A5 and ternary motifs, but seeds are not literal specialisations of φ under tested affine changes. Match is at the level of Gal/passport, not equation identity.

- Seeds tested: 12
- Exact affine fibre matches: **0**
- Seeds with Gal A5: 12
- Motif: `{'phi_coeffs': [6, -15, 10], 'in_3Z': False, 'generations_visible': True, 'relation_to_seeds': 'φ coeffs are ternary (multiples of 3); seeds use 55=61-6, 88=61+27, etc. No seed is an affine fibre of φ under the tested (λ,μ,w) grid.'}`

### Fibre Gal samples (φ(y)=w)

- w=0: status=reducible gal=None disc_sq=None poly=`x**5 - 15*x**4 + 60*x**3`
- w=1: status=reducible gal=None disc_sq=None poly=`x**5 - 15*x**4 + 60*x**3 - 1296`
- w=2: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 2592`
- w=3: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 3888`
- w=1/2: status=reducible gal=None disc_sq=None poly=`x**5 - 30*x**4 + 240*x**3 - 20736`
- w=6: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 7776`
- w=-1: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 + 1296`
- w=61: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 79056`
- w=80: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 103680`
- w=539: status=odd_monodromy gal=S5TransitiveSubgroups.S5 disc_sq=False poly=`x**5 - 15*x**4 + 60*x**3 - 698544`

- disc(φ−w) as poly in w (preview): `4050000*w**2*(w - 1)**2`

### Homogenised regression (still A5)

- seed=(-55, 88) t=1: HIT_A5 S5TransitiveSubgroups.A5
- seed=(-55, 88) t=3: HIT_A5 S5TransitiveSubgroups.A5
- seed=(-55, 88) t=61: HIT_A5 S5TransitiveSubgroups.A5
- seed=(-55, -88) t=1: HIT_A5 S5TransitiveSubgroups.A5
- seed=(-55, -88) t=3: HIT_A5 S5TransitiveSubgroups.A5
- seed=(-55, -88) t=61: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, 76) t=1: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, 76) t=3: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, 76) t=61: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, -76) t=1: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, -76) t=3: HIT_A5 S5TransitiveSubgroups.A5
- seed=(95, -76) t=61: HIT_A5 S5TransitiveSubgroups.A5

### Conclusion (Attack 1)

- **Gal-level match:** yes (seeds and cover both A5).
- **Equation-level match:** no exact affine fibre of φ equals an HQCC seed in the search grid.
- Compatibility is via **shared monodromy / passport / ternary motif**, not identity of polynomials.

---

## Attack 2 — Resonant base parametrisation

**Verdict:** Resonant base: use Möbius sending (3,61,539)→(0,1,∞) as primary HQCC-native coordinate on the Belyi base; monodromy unchanged; labels match 9 Maths (generations, punctures, period).

- G4 = 539.9, rational approx 5399/10
- ξ note: ξ = 2 cos(2π/539.9) is transcendental-looking as written; classical stand-in 2 cos(2π/540) has cyclotomic minpoly (recorded).
- 2cos(2π/540) minpoly: `z**72 - 72*z**70 + 2484*z**68 - 54672*z**66 + 862290*z**64 - 10378368*z**62 + 99118656*z**60 - 771164928*z**58 + 4979436039*z**56 - 27048788359*z**54 + 124860697182*z**52 - 493488179055*z**50 + 1679230591920*z**48 - 4939170892236*z**46 + 12591234674136*z**44 - 27862699620789*z**42 + 53548580107737*z**40 - 89342899774245*z**38 + 129234068042378*z**36 - 161701290086427*z**34 + 174446876713788*z**32 - 161571652971864*z**30 + 127776527051184*z**28 - 85697684743788*z**26 + 48336595129890*z**24 - 22692664633986*z**22 + 8754890596452*z**20 - 2731802488955*z**18 + 675633064770*z**16 - 129028404735*z**14 + 18376162194*z**12 - 1859896818*z**10 + 124688565*z**8 - 4955895*z**6 + 96174*z**4 - 648*z**2 + 1`

### Primary Möbius (3, 61, 539) → (0, 1, ∞)

```
{
  "p": 3,
  "q": 61,
  "r": 539,
  "t(z)": "((-239/29)(z-3))/(z - 539)",
  "t_latex": "\\frac{24.7241379310345 - 8.24137931034483 z}{z - 539}"
}
```

### Proposal

```
{
  "base": "P1_s",
  "coordinate": "t(s) = ((s-3)/(s-539)) / ((61-3)/(61-539))  so t(3)=0, t(61)=1, t(539)=\u221e",
  "cover": "\u03c6(y) = t(s)  i.e. 6y^5-15y^4+10y^3 - t(s) = 0",
  "effect": "Same geometric monodromy A5 (base automorphism). Branch points labelled by model (3,61,539) = generations/punctures/period.",
  "nativeness": "This realises Step-4 dictionary at the level of base points; not yet a dynamical embedding of T3 orbits."
}
```

### Fibres at resonant base values

- **t=0** t=0: reducible gal=None poly=`x**5 - 15*x**4 + 60*x**3`
- **t=1** t=1: reducible gal=None poly=`x**5 - 15*x**4 + 60*x**3 - 1296`
- **generations_3** t=3: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 15*x**4 + 60*x**3 - 3888`
- **punctures_61** t=61: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 15*x**4 + 60*x**3 - 79056`
- **period_539** t=539: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 15*x**4 + 60*x**3 - 698544`
- **flux_ratio_80** t=80: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 15*x**4 + 60*x**3 - 103680`
- **towers_243** t=243: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 15*x**4 + 60*x**3 - 314928`
- **Nflux_norm** t=4880/539: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 8085*x**4 + 17431260*x**3 - 533801616089575680`
- **G4_frac** t=5399/10: odd_monodromy gal=S5TransitiveSubgroups.S5 poly=`x**5 - 150*x**4 + 6000*x**3 - 69971040000`

### Conclusion (Attack 2)

- Base coordinate **t(s)** with t(3)=0, t(61)=1, t(539)=∞ is the recommended HQCC-native chart.
- Geometric monodromy remains **A5** (Möbius base change).
- Branch labels = generations / punctures / period (9 Maths).

---

## Attack 3 — Functor T₃ → braids

**Verdict:** Scaffold functor F: residue paths under T3 → words in {σ0,σ1} → A5. Well-defined as a map on finite paths; not yet a unique natural transformation from the T3 category. Residue 2 is the main ambiguity.

### Monodromy generators (preferred φ)

- σ0 cycles: `[(0, 1, 2), (3,), (4,)]` perm=`[1, 2, 0, 3, 4]`
- σ1 cycles: `[(0,), (1,), (2, 4, 3)]` perm=`[0, 1, 4, 2, 3]`
- σ∞ cycles: `[(0, 3, 4, 2, 1)]` perm=`[3, 0, 1, 4, 2]`

### Residue → class

```
{
  "0": {
    "class": "3A",
    "generator": "sigma_0",
    "T3": "n//3 contraction",
    "base_point": 0
  },
  "1": {
    "class": "3A",
    "generator": "sigma_1",
    "T3": "(4n+2)//3 expansion",
    "base_point": 1
  },
  "2": {
    "class": "3A_or_mix",
    "generator": "sigma_0 o sigma_1  (provisional)",
    "T3": "(2n+1)//3 second branch",
    "note": "No dedicated third finite branch point on Belyi base; map via word in \u03c30,\u03c31"
  }
}
```

### Functor axioms (scaffold)

```
{
  "F_objects": "Obj: natural numbers (or N-orbits under T3) and marked residue sequences; also the base orbifold P1\\{0,1,\u221e}",
  "F_morphisms": "A step n\u2192T3(n) with residue r maps to the loop generator \u03c3_r in \u03c01(P1\\{0,1,\u221e}) \u2245 <\u03c30,\u03c31 | (with \u03c3\u221e=(\u03c30\u03c31)^{-1})>, then via monodromy rep \u03c1: \u03c01 \u2192 A5 \u2282 S5",
  "F_composition": "Path concatenation \u2192 word multiplication in \u03c01 / A5",
  "F_identity": "Constant path at fixed n \u2192 identity braid / id in A5",
  "limitations": "Residue 2 has no dedicated branch point on the 3-point Belyi base; encoded as word \u03c30\u03c31. Not unique; different encodings give conjugate functors. Not yet natural w.r.t. T3 conjugacy of trajectories."
}
```

### Sample T₃ paths → A5 cycle types

Histogram: `{'(5,)': 7, '(3, 1, 1)': 5}`

- n0=1: path=[1, 2, 1, 2, 1, 2, 1, 2, 1, 2] word_len=16 type=(5,)
- n0=2: path=[2, 1, 2, 1, 2, 1, 2, 1, 2, 1] word_len=16 type=(5,)
- n0=3: path=[0, 1, 2, 1, 2, 1, 2, 1, 2, 1] word_len=16 type=(5,)
- n0=9: path=[0, 0, 1, 2, 1, 2, 1, 2, 1, 2] word_len=16 type=(3, 1, 1)
- n0=27: path=[0, 0, 0, 1, 2, 1, 2, 1, 2, 1] word_len=16 type=(3, 1, 1)
- n0=61: path=[1, 1, 2, 1, 2, 2, 1, 1, 0, 2] word_len=16 type=(5,)
- n0=80: path=[2, 2, 2, 2, 0, 2, 0, 1, 2, 1] word_len=16 type=(3, 1, 1)
- n0=243: path=[0, 0, 0, 0, 0, 1, 2, 1, 2, 1] word_len=16 type=(5,)
- n0=539: path=[2, 2, 2, 0, 2, 2, 2, 0, 2, 0] word_len=16 type=(5,)
- n0=4880: path=[2, 1, 0, 0, 2, 0, 2, 2, 2, 1] word_len=16 type=(3, 1, 1)
- n0=100: path=[1, 2, 2, 2, 0, 1, 0, 0, 2, 1] word_len=16 type=(3, 1, 1)
- n0=1617: path=[0, 2, 2, 2, 0, 2, 2, 2, 0, 2] word_len=16 type=(5,)

### Conclusion (Attack 3)

- Defined a **scaffold functor** from finite T₃ residue paths to words in ⟨σ0,σ1⟩ ⊂ A5.
- Residue 2 is ambiguous (encoded as σ0σ1).
- Not yet unique/natural; sufficient as an operational bridge for further refinement.

---

## Overall status after queued attacks

| Attack | Result |
|--------|--------|
| 1 Specialisation match | Gal A5 shared; **no** exact affine equation match |
| 2 Resonant base | **Primary chart** (3,61,539)→(0,1,∞); monodromy preserved |
| 3 T₃→braid functor | **Scaffold** path→word→A5; residue-2 ambiguity noted |

Arithmetic theorem-grade + geometric A5 cover + labelled base + operational functor scaffold.
Deep remaining: naturality of F and Hilbert recovery of HQCC seeds from the cover family.

_Generated by queued_attacks.py_