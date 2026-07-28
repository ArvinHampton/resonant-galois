# New algebraic ideas A–F

_Elapsed: 112.23s_

**Verdict:** New algebraic ideas (112.23s). **Primary path A+F HIT.** Mestre: all 10 even HQCC seeds have 1-dim \(R\)-space; `shift_y_tR` families have disc identically □ in \(\mathbb{Q}(t)\) with sample Gal \(A_5\). Embed into \(T\): BJ seeds + depressed Mestre specialisations embed (90 specs). **Secondary B HIT:** \(x^5+75x^3+A x^2+3A\) has disc \(=324 A^2(A^2+84375)^2\) identically □; all 16 lattice \(A\) tested give \(A_5\); matches \(\chi_T\) with \(d=-75\). C/D/E: no new identical-square avatar by construction. Old \(T\) cuts remain closed negative experiment.

---

## Stance

- Arithmetic centre (pure-even multi-\(k\)) **finished** — `PURE_EVEN_MULTI_K.md`.
- Criterion 2 on template \(T(a,\ldots,f)\) **closed** negative — `TIER11_DEEPEN.md`.
- This document: **new equations / new matrices**, not more cuts of \(T\).

### Do not retry

| Approach | Reason |
|----------|--------|
| More linear/bilinear cuts of same \(T\) | Exhausted |
| \(F\to T\) only hoping disc□→1 | Rate 0 |
| Surgery on rigid \(\varphi/\mathbb{Q}\) | Permanent factor 5 |
| Binary Collatz forces evenness without pure-even | Composite only |

---

## Idea A — Mestre deformation (primary)

For irreducible seed \(P\) with disc □, solve
$$P''R-2P'R'\equiv 0\pmod{P},\quad \deg R<\deg P$$
(\(\deg R\le n-1\); for \(n=5\) one finds \(\dim=1\) on pure-even BJ seeds).
Then build families via resultants:
- `shift_y_tR`: \(\operatorname{Res}_y(P(y),\,z-y-t R(y))\)
- `uPp_tR`: \(\operatorname{Res}_y(P(y),\,u P'(y)-t R(y))\)
- `uPp_R_t`: \(\operatorname{Res}_y(P(y),\,u P'(y)-R(y)-t)\)

- Seeds processed: **10**
- Seeds with nontrivial \(R\)-space: **10**
- Family constructions with disc identically □ in \(\mathbb{Q}(t)\): **30**

### Seed `flagship` — \(P=x^5+(-55)x+(88)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 + 8*x**3 - 32*x**2 + 33']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 + 8*x**3 - 32*x**2 + 33` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 + 8*x**3 - 32*x**2 + 33` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 + 8*x**3 - 32*x**2 + 33` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `flagship_m` — \(P=x^5+(-55)x+(-88)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 - 8*x**3 - 32*x**2 + 33']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 - 8*x**3 - 32*x**2 + 33` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 - 8*x**3 - 32*x**2 + 33` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 - 8*x**3 - 32*x**2 + 33` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `classical_95_76` — \(P=x^5+(95)x+(76)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 - 4*x**3 - 8*x**2 - 57']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 - 4*x**3 - 8*x**2 - 57` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 - 4*x**3 - 8*x**2 - 57` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 - 4*x**3 - 8*x**2 - 57` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `classical_95_m76` — \(P=x^5+(95)x+(-76)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 + 4*x**3 - 8*x**2 - 57']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 + 4*x**3 - 8*x**2 - 57` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 + 4*x**3 - 8*x**2 - 57` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 + 4*x**3 - 8*x**2 - 57` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `classical_20_16` — \(P=x^5+(20)x+(16)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 - 4*x**3 - 8*x**2 - 12']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 - 4*x**3 - 8*x**2 - 12` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 - 4*x**3 - 8*x**2 - 12` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 - 4*x**3 - 8*x**2 - 12` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `classical_20_m16` — \(P=x^5+(20)x+(-16)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 + 4*x**3 - 8*x**2 - 12']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 + 4*x**3 - 8*x**2 - 12` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 + 4*x**3 - 8*x**2 - 12` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 + 4*x**3 - 8*x**2 - 12` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `lsw_slice_100` — \(P=x^5+(-100)x+(400)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 + 20*x**3 - 200*x**2 + 60']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 + 20*x**3 - 200*x**2 + 60` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 + 20*x**3 - 200*x**2 + 60` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 + 20*x**3 - 200*x**2 + 60` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `lsw_slice_124` — \(P=x^5+(124)x+(-496)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['5*x**4 + 100*x**3 - 1000*x**2 - 372']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `5*x**4 + 100*x**3 - 1000*x**2 - 372` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `5*x**4 + 100*x**3 - 1000*x**2 - 372` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `5*x**4 + 100*x**3 - 1000*x**2 - 372` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `flag_m_cleared` — \(P=x^5+(320)x+(-512)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['x**4 + 8*x**3 - 32*x**2 - 192']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `x**4 + 8*x**3 - 32*x**2 - 192` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `x**4 + 8*x**3 - 32*x**2 - 192` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `x**4 + 8*x**3 - 32*x**2 - 192` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

### Seed `lsw_m` — \(P=x^5+(-3121)x+(12484)\)

- seed disc□: **True**, null_dim(R): **1**
- R basis: `['5*x**4 + 100*x**3 - 1000*x**2 + 9363']`

| construction | R | disc□ in Q(t)? | samples |
|--------------|---|:--------------:|---------|
| shift_y_tR | `5*x**4 + 100*x**3 - 1000*x**2 + 9363` | **True** | t=1:HIT_A5,t=2:HIT_A5,t=3:HIT_A5,t=-1:HIT_A5 |
| uPp_tR | `5*x**4 + 100*x**3 - 1000*x**2 + 9363` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |
| uPp_R_t | `5*x**4 + 100*x**3 - 1000*x**2 + 9363` | **True** | t=1:not_monic_Z,t=2:not_monic_Z,t=3:not_monic_Z,t=-1:not_monic_Z |

---

## Idea F — embed even families into \(T\) (primary companion)

Solve \(\chi_T(a,b,c,d,e,f)=P(x)\) (after killing \(x^4\) by shift if needed).

- Static embeddable (seeds / non-BJ points): **10**
- Mestre specialisation embeddable: **90**

- **seed:flagship**: embeddable=**True** ({'a': '-e*f', 'b': '(55/2 - sqrt(352*e*f + 3025)/2)/f', 'c': 'sqrt(11)*sqrt(32*e*f + 275)/(2*e) + 55/(2*e)'})
- **seed:flagship_m**: embeddable=**True** ({'a': '-e*f', 'b': '(55/2 - sqrt(-352*e*f + 3025)/2)/f', 'c': 'sqrt(11)*sqrt(-32*e*f + 275)/(2*e) + 55/(2*e)'})
- **seed:classical_95_76**: embeddable=**True** ({'a': '-e*f', 'b': '(-sqrt(304*e*f + 9025)/2 - 95/2)/f', 'c': 'sqrt(19)*sqrt(16*e*f + 475)/(2*e) - 95/(2*e)'})
- **seed:classical_95_m76**: embeddable=**True** ({'a': '-e*f', 'b': '(-sqrt(-304*e*f + 9025)/2 - 95/2)/f', 'c': 'sqrt(19)*sqrt(-16*e*f + 475)/(2*e) - 95/(2*e)'})
- **seed:classical_20_16**: embeddable=**True** ({'a': '-e*f', 'b': '2*(-sqrt(4*e*f + 25) - 5)/f', 'c': '2*sqrt(4*e*f + 25)/e - 10/e'})
- **seed:classical_20_m16**: embeddable=**True** ({'a': '-e*f', 'b': '2*(-sqrt(-4*e*f + 25) - 5)/f', 'c': '2*sqrt(-4*e*f + 25)/e - 10/e'})
- **nonBJ_A=3**: embeddable=**True** ({'a': '-e*f - 3', 'b': 'sqrt(3)*e*sqrt(f*(-25*e*f - 72)/e)/f', 'c': '-sqrt(3)*sqrt(-f*(25*e*f + 72)/e)'})
- **nonBJ_A=9**: embeddable=**True** ({'a': '-e*f - 9', 'b': 'sqrt(3)*e*sqrt(f*(-25*e*f - 216)/e)/f', 'c': '-sqrt(3)*sqrt(-f*(25*e*f + 216)/e)'})
- **nonBJ_A=27**: embeddable=**True** ({'a': '-e*f - 27', 'b': 'sqrt(3)*e*sqrt(f*(-25*e*f - 648)/e)/f', 'c': '-sqrt(3)*sqrt(-f*(25*e*f + 648)/e)'})
- **nonBJ_A=61**: embeddable=**True** ({'a': '-e*f - 61', 'b': 'sqrt(3)*e*sqrt(f*(-25*e*f - 1464)/e)/f', 'c': '-sqrt(3)*sqrt(-f*(25*e*f + 1464)/e)'})
- **mestre_hit:flagship:shift_y_tR**: specs=[(1, True), (2, True), (3, True)]
- **mestre_hit:flagship:uPp_tR**: specs=[(1, True), (2, True), (3, True)]

**Note:** BJ seeds embed via classical BJ-embed (\(d=0,a=-ef\)). That recovers pure-even, not a new HQCC-native necessity fragment. Non-BJ embeds (when they exist) give \(d\neq 0\) realisations inside \(T\).

---

## Idea B — non-BJ degree-1 \(A_5\) family (secondary)

Family: `$x^5+75*x^3+A*x^2+3*A$`

- Disc identically square in parameter \(A\)? **True**
- Factored disc (preview): `324*A**2*(A**2 + 84375)**2`
- Lattice samples disc□: **16**, A5: **16**
- Generic \(T\)-match solutions: **2**
- Sparse match sample: `[{'e': '0', 'a': '-A', 'b': '72*A/c', 'f': '0', 'd': '-75'}, {'f': '0', 'a': '-A', 'b': '72*A/c', 'e': '0', 'd': '-75'}, {'a': '0', 'b': '-sqrt(3)*I*A/f', 'c': '-sqrt(3)*I*f', 'e': '-A/f', 'd': '-75'}]`

| \(A\) | disc□ | status |
|----:|:-----:|--------|
| 1 | True | HIT_A5 |
| 3 | True | HIT_A5 |
| 9 | True | HIT_A5 |
| 27 | True | HIT_A5 |
| 61 | True | HIT_A5 |
| 80 | True | HIT_A5 |
| 243 | True | HIT_A5 |
| 539 | True | HIT_A5 |
| -3 | True | HIT_A5 |
| -9 | True | HIT_A5 |
| 18 | True | HIT_A5 |
| 54 | True | HIT_A5 |
| 4880 | True | HIT_A5 |
| 55 | True | HIT_A5 |
| 88 | True | HIT_A5 |
| 95 | True | HIT_A5 |

---

## Idea C — change the matrix avatar (if A+F stall)

- C1 ternary-coeff companion scan: tested=401, irr=332, disc□=0, A5=0
- C2 transfer graph: disc□=False, status=None
- C3 3-cycle block deform: `[{'u': 0, 'irr': False}, {'u': 1, 'irr': False}, {'u': 3, 'irr': False}, {'u': 9, 'irr': False}, {'u': 61, 'irr': False}, {'u': 80, 'irr': False}, {'u': 243, 'irr': False}]`
- **Identically square by construction?** **False**

Old \(T\) remains a **closed negative experiment** for Crit 2. New avatars need a built-in evenness identity (as pure-even has), not another sparse scan.

---

## Idea D — icosahedral / invariant parameters

- Family `x^5+5m x^3+5m^2 x+n`
- irr samples=50, disc□=0
- Even hits (sample): `[]`
- Icosahedral-adjacent scans; disc□ sporadic not identical in (m,n)

---

## Idea E — HQCC-native polynomial from T₃

- Orbit polys: `[{'n0': 1, 'orbit': [1, 2], 'deg': 2, 'irr': False}, {'n0': 2, 'orbit': [2, 1], 'deg': 2, 'irr': False}, {'n0': 3, 'orbit': [3, 1, 2], 'deg': 3, 'irr': False}, {'n0': 61, 'orbit': [61, 82, 110, 73, 98, 65, 43, 58, 78, 26, 17, 11], 'deg': 12, 'irr': False}, {'n0': 80, 'orbit': [80, 53, 35, 23, 15, 5, 3, 1, 2], 'deg': 9, 'irr': False}, {'n0': 539, 'orbit': [539, 359, 239, 159, 53, 35, 23, 15, 5, 3, 1, 2], 'deg': 12, 'irr': False}]`
- deg5 irr=0, disc□=0
- Orbit polys from T3 are design probes; not systematic A5 machines

---

## Recommended resolution path (updated)

| Priority | Idea | Outcome this run | Next |
|:--------:|------|------------------|------|
| 1 | **A Mestre + F embed** | **HIT:** \(R\)-space dim 1 on all 10 seeds; `shift_y_tR` disc□ + sample \(A_5\); embeds into \(T\) after depress | Extract closed form \(P_t\) for flagship; ask which embed coords are HQCC-native (beyond BJ-embed) |
| 2 | **B non-BJ deg-1** | **HIT:** disc \(=324 A^2(A^2+84375)^2\); lattice \(\Rightarrow A_5\); \(\chi_T\) with \(d=-75\) | Sparse embed on resonant lattice |
| 3 | C new matrix | no identical-square avatar yet | only with evenness identity by design |
| 4 | D, E | disc□=0 in scans | low priority |

### Concrete wins (citable additions to the centre)

1. **Mestre lift of HQCC seeds:** every pure-even BJ seed tested admits a 1-param \(A_5\) family over \(\mathbb{Q}(t)\) via `shift_y_tR` (build evenness first; not a cut of \(T\)).
2. **Non-BJ even family in \(T\):** \(x^5+75x^3+A x^2+3A\) is identically disc-square and realisable with \(d=-75\neq 0\) — **beyond pure-even BJ-embed**.
3. Template \(T\) still does not *force* disc□ (Crit 2 closed negative), but it **hosts** these new even families.

### Sparse \(T\)-realisation of Idea B

\[
d=-75,\quad e=f=0,\quad a=-A,\quad bc=72A
\quad\Rightarrow\quad
\chi_T=x^5+75x^3+A x^2+3A.
\]
Example: \(A=3\), \(b=1\), \(c=216\); or \(A=61\), \(b=1\), \(c=4392\).

## One-line synthesis

The Model’s arithmetic centre is finished; Criterion 2 on the present template is closed as a *forcing* claim. **New progress:** Mestre lifts of HQCC seeds and a non-BJ degree-1 family with identically square disc embedding in \(T\) at \(d\neq 0\) — by changing the equation, not by cutting the same variety.

```bash
python new_algebraic_ideas.py
```

_Generated by new_algebraic_ideas.py_