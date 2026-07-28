# Which pure-even \(k\in R_n\) arise as cross-ratios of a cover over \(R_n\)?

_Elapsed: 13.98s_

**Question (N3 / Q1).** Among pure-even ratio classes \(k=\beta/\alpha\in R_n\),
which arise from covers defined over the real subfield
\(R_n=\mathbb{Q}(2\cos 2\pi/n)\)?

---

## 0. Two different quantities (do not conflate)

| symbol | meaning | lives in |
|--------|---------|----------|
| \(s\) | **Cross-ratio** of four branch points (Hurwitz / \(M_{0,4}\)) | \(\mathbb{P}^1\) |
| \(k\) | **Pure-even ratio** \(\beta/\alpha\) of a BJ fibre \(x^5+\alpha x+\beta\) | field of coeffs |

A cover over \(R_n\) has branch cross-ratio \(s\in\mathbb{P}^1(R_n)\) when the
branch locus is \(R_n\)-rational. Separately, a BJ model of a fibre may have
ratio \(k\in R_n\). **In general \(k\neq s\)**; at best \(k=\kappa(s)\) for an
unknown algebraic map \(\kappa\) attached to the cover type.

---

## 1. Arithmetic answer (complete)

**Identity:** `True`

For every k ∈ R_n (any n) and m ∈ R_n \ {0} with α(m,k)≠0, the BJ fibre x⁵+αx+β over R_n is pure-even: disc is a square in R_n.

$$\alpha(m)=256m^2-\frac{3125k^4}{256},\quad\beta=k\alpha,\quad\operatorname{disc}=(256\alpha^2 m)^2\quad\text{in }R_n(m).$$

**Corollary.** *Every* \(k\in R_n\setminus\{0\}\) arises as a pure-even parameter
over \(R_n\). There is no arithmetic restriction beyond \(k\in R_n\).

---

## 2. Branch cross-ratio \(s\) over \(R_n\) (complete, almost trivial)

**Theorem.** For any \(s\in R_n\setminus\{0,1\}\), the ordered 4-tuple
\((0,1,\infty,s)\) consists of \(R_n\)-rational points of \(\mathbb{P}^1\), and \(s\)
is their cross-ratio. Hence **every** such \(s\) arises as the branch
cross-ratio of a 4-point configuration over \(R_n\).

So if the question is read as “which \(s\in R_n\) are cross-ratios of covers
over \(R_n\)?”, the answer is: **all \(s\in R_n\setminus\{0,1,\infty\}\)**
(subject only to the cover existing for that Nielsen type at that \(s\)).

Existence of an \(A_5\) cover of type e.g. \(3A^4\) at a given \(s\in R_n\) is a
**Hurwitz** question (genus-0 reduced space \(\cong\mathbb{P}^1_s\) over \(\mathbb{Q}\)
already known; base change to \(R_n\) is free).

---

## 3. Non-trivial reading: pure-even \(k\) from geometric fibres

**Intended meaning:** which \(k\in R_n\) appear as \(\beta/\alpha\) for a BJ model
of a fibre of a cover \(X\to\mathbb{P}^1\) defined over \(R_n\)?

Every k ∈ R_n is pure-even over R_n (arithmetic). Every s ∈ R_n\{0,1} is a branch cross-ratio over R_n (trivial). Which pure-even k arise as BJ ratios of fibres of covers over R_n is open in general; known constraints and cases are listed below.

### Classification

| class | status | meaning |
|-------|--------|---------|
| `all_of_Rn_as_pure_even_params` | **yes** | k ∈ R_n ⇒ pure-even family over R_n exists |
| `all_of_Rn_as_branch_cross_ratios_s` | **yes** | s ∈ R_n\{0,1} ⇒ 4-branch chart over R_n |
| `catalogue_Q_k_inside_Rn` | **yes** | all multi-seed k ∈ Q ⊂ R_n |
| `cosine_k_as_geometric_BJ_ratio` | **open** | k=2cos(2πp/n) pure-even over R_n; geometric origin unknown |
| `k_from_3A4_fibres_over_R5` | **open** | s=-1 over R_5 known; map fibre→BJ k not closed-formed |
| `k_equals_s_identification` | **not_forced** | Identifying pure-even k with the Hurwitz cross-ratio s is a coordinate choice, not a theorem. In general k = κ(s) for some rational function / algebraic function κ of the cover moduli. |

### Necessary conditions

- NEC1 (arithmetic): k ∈ R_n is necessary and sufficient for the pure-even formulae to be defined over R_n with α,β ∈ R_n(m) when m ∈ R_n.
- NEC2 (branch locus): For a 4-point cover over R_n with R_n-rational branch locus, the cross-ratio s lies in P¹(R_n). This constrains s, not k.
- NEC3 (BJ link): A pure-even k arises geometrically only if some fibre of a cover over R_n is R_n-birational to a BJ quintic with β/α = k (after coordinate change over R_n).
- NEC4 (descent): If the cover is defined over R_n but the fibre field is a proper extension, k may lie in that extension, not in R_n.
- NEC5 (cosine geometry): If branch points are constrained to cyclotomic real loci (cosine relations), then s (and possibly k) lie in a thin subset of R_n — typically multi-angle / Chebyshev values — not all of R_n.

---

## 4. Catalogue \(k\in\mathbb{Q}\) inside every \(R_n\)

All multi-seed pure-even ratios from the HQCC catalogue lie in \(\mathbb{Q}\subset R_n\):

`['-4', '4', '-8/5', '8/5', '4/5', '-4/5', '-12/5', '12/5', '-16/5', '16/5']`

They **do** arise as pure-even parameters over \(R_n\), but their known origin is
**arithmetic** (envelope over \(\mathbb{Q}\)), not as Nielsen-labelled cross-ratios.
Whether a cover over \(R_n\) specialises to these \(k\) remains open
(geometric multi-\(k\) problem).

---

## 5. Cosine candidates \(k=2\cos(2\pi p/d)\in R_n\)

These are the \(k\) “visibly constrained by cosine relations.”
They are pure-even over \(R_n\) whenever \(d\mid n\) (so \(k\in R_d\subset R_n\)
in the cyclotomic plus tower). Matching a catalogue rational is accidental.

### \(n=5\) (deg \(R_n=2\))

| form | deg | numeric | = catalogue \(k\)? |
|------|----:|--------:|:------------------:|
| `2cos(2π·1/5)` | 2 | 0.618034 | False |
| `2cos(2π·2/5)` | 2 | -1.618034 | False |

### \(n=7\) (deg \(R_n=3\))

| form | deg | numeric | = catalogue \(k\)? |
|------|----:|--------:|:------------------:|
| `2cos(2π·1/7)` | 3 | 1.246980 | False |
| `2cos(2π·2/7)` | 3 | -0.445042 | False |
| `2cos(2π·3/7)` | 3 | -1.801938 | False |

### \(n=11\) (deg \(R_n=5\))

| form | deg | numeric | = catalogue \(k\)? |
|------|----:|--------:|:------------------:|
| `2cos(2π·1/11)` | 5 | 1.682507 | False |
| `2cos(2π·2/11)` | 5 | 0.830830 | False |
| `2cos(2π·3/11)` | 5 | -0.284630 | False |
| `2cos(2π·4/11)` | 5 | -1.309721 | False |
| `2cos(2π·5/11)` | 5 | -1.918986 | False |

### \(n=15\) (deg \(R_n=4\))

| form | deg | numeric | = catalogue \(k\)? |
|------|----:|--------:|:------------------:|
| `2cos(2π·1/3)` | 1 | -1.000000 | False |
| `2cos(2π·1/5)` | 2 | 0.618034 | False |
| `2cos(2π·2/5)` | 2 | -1.618034 | False |
| `2cos(2π·1/15)` | 4 | 1.827091 | False |
| `2cos(2π·2/15)` | 4 | 1.338261 | False |
| `2cos(2π·4/15)` | 4 | -0.209057 | False |
| `2cos(2π·7/15)` | 4 | -1.956295 | False |

**No** cosine value \(2\cos(2\pi p/d)\) in the scanned range equals a catalogue rational \(k\) (as expected: catalogue \(k\) are rational; non-rational cosines are irrational).

---

## 6. Known geometric cases in the programme

### Rigid \(\varphi/\mathbb{Q}\) (r=3)

- Fibres odd over Q (disc=5·□); not pure-even source

### \(3A^4\) at \(s=-1\) over \(R_5=\mathbb{Q}(\sqrt5)\)

- Branch cross-ratio \(s=-1\) ∈ \(R_5\): **True**
- Cover params: `{'c': '-√5', 'p2': -1, 'r1': '1/5', 'r2': '-1/5'}`
- Closed form \(f_s\to\mathrm{BJ}\): **False**
- s=-1 ∈ Q ⊂ R_5 is a cross-ratio over R_5. Without f_s ∈ R_5(s)[x] → BJ, no proven list of pure-even k from this cover. Numeric fibres gave even=0 over Q.

### Arithmetic envelope (not geometric)

- Pure-even \(k\): `['-4', '4', '-8/5', '8/5', '4/5', '-4/5', '-12/5', '12/5', '-16/5', '16/5']`
- All catalogue k arise arithmetically; not as Nielsen cross-ratios

### Cosine \(k\) as chart value

- \(k=2cos(2π/5)=(-1+√5)/2\) ∈ \(R_5\), pure-even **True**
- As branch \(s=k\) in chart \(\{0,1,\infty,s\}\): trivial yes
- As BJ ratio from a known cover: **unknown**

---

## 7. Partial answer (lock)

| reading of the question | answer |
|-------------------------|--------|
| Which \(k\in R_n\) admit pure-even families over \(R_n\)? | **All** \(k\in R_n\setminus\{0\}\) |
| Which \(s\in R_n\) are branch cross-ratios over \(R_n\)? | **All** \(s\in R_n\setminus\{0,1\}\) |
| Which catalogue \(k\in\mathbb{Q}\) lie in \(R_n\)? | **All** of them (\(\mathbb{Q}\subset R_n\)) |
| Which \(k\in R_n\) are BJ ratios of fibres of a cover over \(R_n\)? | **Open** — no closed \(s\mapsto k\) |
| Are cosine values distinguished? | **Yes as geometric candidates**; not forced by evenness |
| Is \(k=s\)? | **Not in general** |

### What would finish the non-trivial reading

1. Closed form \(f_s\in R_n(s)[x]\) (or over \(\mathbb{Q}(s)\)) for a Nielsen type.
2. BJ reduction: after Möbius in \(x\), read \(k(s)=\beta(s)/\alpha(s)\in R_n(s)\).
3. Image of \(k: H^{\mathrm{rd}}(R_n)\to\mathbb{P}^1\) — that image **is** the answer.

**Until then:** arithmetic multi-\(k\) over \(\mathbb{Q}\) remains the citable source of
explicit pure-even \(k\); over \(R_n\), pure-even is free in \(k\), while geometric
origin of specific \(k\) is constrained by Hurwitz data not yet converted to BJ.

---

## 8. Direct answers to the one-line question

> Which pure-even \(k\in R_n\) arise as cross-ratios of a cover over \(R_n\)?

**If “cross-ratio” means the Hurwitz parameter \(s\):**  
every \(s\in R_n\setminus\{0,1\}\) (and the cover type must exist at that \(s\)).

**If “cross-ratio” is used loosely for pure-even \(k\) coming from such a cover:**  
unknown list; must equal the image of the moduli map \(s\mapsto k(s)\) once a
BJ model exists. Currently known rigorously:

- all \(k\in R_n\) work **arithmetically** as pure-even parameters;
- all catalogue \(k\in\mathbb{Q}\) embed in every \(R_n\);
- \(s=-1\in R_5\) is a genuine geometric cross-ratio of a \(3A^4\) cover over \(R_5\);
- no theorem yet names a non-rational \(k\in R_n\setminus\mathbb{Q}\) as a BJ ratio
  of an \(R_n\)-cover fibre.

```bash
python pure_even_k_cross_ratio.py
```

_Generated by pure_even_k_cross_ratio.py_