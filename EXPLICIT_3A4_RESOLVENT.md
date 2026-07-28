# Explicit model of \(\mathrm{Ni}(A_5,C_3^4)\) and specialisation test

_Elapsed: 29.51s_

**Verdict:** Reduced Hurwitz curve modelled as P^1_s (g=0). Explicit (3,1,1)^4 covers constructed numerically for 20 rational s (Newton on triple-root equations); exact model at s=-1 over Q(√5). Fibres: irr=620, even=0, A5≈0, BJ-hits=0, catalogue hits=0, catalogue k=[], multi_cat=False. No geometric multi-k catalogue hit yet; multi-k remains arithmetic (envelope).

---

## 1. Explicit model of the reduced Hurwitz curve

By the genus lock (`GENUS_3A4_LOCK.md`), the reduced Hurwitz curve is
**isomorphic to \(\mathbb{P}^1\) over \(\mathbb{Q}\)**. An explicit model is therefore

$$H^{\mathrm{rd}}\ \cong\ \mathbb{P}^1_s,\qquad s=\text{cross-ratio of the four branch points.}$$

Infinitely many rational points: all \(s\in\mathbb{Q}\setminus\{0,1\}\) (and \(\infty\)).

The map \(H^{\mathrm{rd}}\to M_{0,4}\cong\mathbb{P}^1\) may be taken as the identity in this coordinate
(branch points placed at \(0,1,\infty,s\)).

---

## 2. Degree-5 resolvent form

For each \(s\), a genus-0 cover of type \((3,1,1)^4\) with branch values \(\{0,1,\infty,s\}\)
is realised in normal form

$$\varphi_s(y)=\frac{c\, y^3(y-1)(y-p_2)}{(y-r_1)(y-r_2)},$$

where \((c,p_2,r_1,r_2)\) are determined by the conditions that \(\varphi_s-1\) and
\(\varphi_s-s\) each have a **triple root** (remaining two preimages simple).

The fibre polynomial (resolvent in the base coordinate \(t=\varphi_s(y)\)) is

$$N_s(y)-t\,D_s(y)=0,\qquad N_s=c y^3(y-1)(y-p_2),\quad D_s=(y-r_1)(y-r_2),$$

cleared to a monic polynomial \(f_{s,t}(y)\in\mathbb{C}[y]\) of degree 5.

**Form:** `f_{s,t}(y) = clear( c(s) y^3 (y-1)(y-p2(s)) - t (y-r1(s))(y-r2(s)) ), monic in y, where (c,p2,r1,r2)(s) solve the triple-root conditions for branch values 1 and s.`

### Exact point \(s=-1\) (over \(\mathbb{Q}(\sqrt{5})\))

- \(c=-\sqrt{5}\), \(p_2=-1\), \(r_1=1/5\), \(r_2=-1/5\)
- \(N=-sqrt(5)*y**5 + sqrt(5)*y**3\)
- \(D=y**2 - 1/25\)
- Triple-root conditions verified symbolically.
- Descent of this single fibre to a deg-5 model over \(\mathbb{Q}(t)\) fails in general;
  the norm to \(\mathbb{Q}(t)\) yields a degree-10 equation
  \(5(y^5-y^3)^2-t^2(y^2-1/25)^2=0\).

### Numerical covers at rational \(s\)

- Requested: 21, solved: **20**

| \(s\) | \(p_2\) | \(c\) | \(r_1\) | \(r_2\) | Newton res |
|------|--------|------|--------|--------|------------|
| -3 | -1.627702 | -1.851043 | -0.330475 | 0.203032 | 3.7e-15 |
| -2 | -1.361936 | -1.983405 | -0.274090 | 0.201250 | 1.1e-15 |
| -1 | -1.000000 | -2.236068 | -0.200000 | 0.200000 | 8.5e-15 |
| -1/2 | -0.734249 | -2.505254 | -0.147768 | 0.201250 | 2.0e-15 |
| 1/2 | 2.618034 | -0.263932 | 1.047214 | 0.400000 | 1.2e-15 |
| 3/2 | 2.394613 | -0.506086 | 1.021345 | 0.426518 | 1.1e-15 |
| 2 | 2.618034 | -0.527864 | 1.047214 | 0.400000 | 3.8e-15 |
| 5/2 | 2.794757 | -0.546346 | 1.072301 | 0.383683 | 2.0e-15 |
| 3 | 2.944332 | -0.561919 | 1.095975 | 0.372232 | 1.8e-15 |
| 4 | 3.193144 | -0.586834 | 1.139178 | 0.356758 | 5.2e-15 |
| 5 | 3.399123 | -0.606168 | 1.177744 | 0.346485 | 4.0e-15 |
| -5/2 | -0.665424 | 6.478719 | 0.202147 | -0.134513 | 2.0e-15 |
| 5/3 | — | — | — | — | FAIL 0.016970877187570496 |
| 7/2 | 3.075501 | -0.575246 | 0.363592 | 1.118226 | 8.8e-15 |
| 1/3 | 0.339636 | -4.780934 | 0.126423 | 0.372232 | 3.2e-15 |
| 2/3 | 0.417604 | -4.632751 | 0.178115 | 0.426518 | 3.2e-16 |
| 4/3 | 2.301192 | -0.498551 | 1.013009 | 0.440211 | 5.2e-15 |
| 5/4 | 2.248124 | -0.495080 | 1.009062 | 0.448846 | 9.8e-16 |
| -3/2 | -0.834343 | 3.589567 | -0.167231 | 0.200435 | 1.8e-15 |
| 7/3 | 2.739486 | -0.540544 | 0.388426 | 1.064087 | 2.6e-15 |
| 8/3 | 2.847082 | -0.551825 | 0.379460 | 1.080355 | 2.0e-15 |

### Interpolation of parameters (numeric)

```{
  "attempted": true,
  "p2": {
    "poly_coeffs_high_to_low": [
      0.007130658257197161,
      -0.059103760009810664,
      -0.01676554481145891,
      1.228049480897495,
      0.538228538482208
    ],
    "max_abs_fit_error": 1.4769143998819534,
    "deg": 4
  },
  "c": {
    "poly_coeffs_high_to_low": [
      -0.02784504361866591,
      0.08018039759798518,
      0.4875922074273065,
      -0.8902073975855007,
      -1.7207586913393134
    ],
    "max_abs_fit_error": 5.267023338829051,
    "deg": 4
  },
  "r1": {
    "poly_coeffs_high_to_low": [
      0.0058042305800894535,
      -0.028592862546019594,
      -0.04858676432598846,
      0.39874447341736385,
      0.35967896686349454
    ],
    "max_abs_fit_error": 0.5035204264162757,
    "deg": 4
  },
  "r2": {
    "poly_coeffs_high_to_low": [
      -0.0012318025368440527,
      -0.004173910603053401,
      0.024221113480424352,
      0.14043161052682407,
      0.28807676997523585
    ],
    "max_abs_fit_error": 0.40600270800692373,
    "deg": 4
  }
}```

**Closed form** \(c(s),p_2(s),r_1(s),r_2(s)\in\mathbb{Q}(s)\) was **not** obtained in this run;
coefficients involve algebraic functions of \(s\) (often in quadratic extensions).
A global model over \(\mathbb{Q}(s)\) may require a different coordinate on \(H^{\mathrm{rd}}\)
or a resolvent of degree \(>5\) over \(\mathbb{Q}(s)\).

---

## 3. Specialisation at rational \((s,t)\)

| quantity | value |
|----------|------:|
| irreducible fibres | 620 |
| even disc | 0 |
| A5 (among even, galois run) | 0 |
| BJ reductions | 0 |
| catalogue seed hits | 0 |
| distinct catalogue \(k\) | [] |
| **multi catalogue \(k\)** | **False** |

### Catalogue hits

_None._ Geometric fibres of this normal form did not land on catalogue BJ seeds in the scan.

### Sample A5 fibres (even)

_No A5 even fibres recorded (or galois not triggered)._

### Sample BJ reductions

_No pure BJ (x^5+αx+β) fibres after y^4-shift in the scan._

---

## 4. Multi-\(k\) conclusion

**Geometric multi-\(k\) catalogue hit: False**

- The **reduced Hurwitz curve** has an explicit rational model \(\mathbb{P}^1_s\).
- **Degree-5 covers** of type \((3,1,1)^4\) are constructed for many rational \(s\)
  (numeric parameters; exact at \(s=-1\) over \(\mathbb{Q}(\sqrt{5})\)).
- Specialisations tested against the pure-even fixed-\(k\) catalogue:
  multi catalogue \(k\) = **False**.

Until a closed-form \(f_s\in\mathbb{Q}(s)[x]\) produces a multi-\(k\) catalogue hit,
**multi-\(k\) success remains arithmetic** (envelope/paths in
`REALISE_3A4_SPECIALISE.md` / `NONRIGID_HURWITZ_SEARCH.md`), **not geometric**
(Nielsen-labelled).

### Next to close the geometric gap

1. Determine \(c,p_2,r_1,r_2\) as algebraic functions of \(s\) in closed form
   (resultants on the triple-root ideal), or find a better coordinate on \(H^{\mathrm{rd}}\).
2. Produce \(f_s\in\mathbb{Q}(s)[x]\) (possibly after a resolvent of the parameter field).
3. Re-run the catalogue specialisation test on that closed form.

_Generated by build_3a4_resolvent.py_