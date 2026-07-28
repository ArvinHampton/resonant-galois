# Explicit equation — genus-0 \(3A^4\) reduced Hurwitz data

_Elapsed: 2.18s_

**Verdict:** Explicit 3A^4 model (2.18s). H^rd ≅ P^1 (g=0). Cover params via plane curve P(q,w)=0 with rational formulae for (c,σ,π,s) and quadratic p2(q,w). Known fibre s=-1 verified (True). Resolvent: monic_y(N-tD), deg 5. Single-valued f_s∈Q(s)[y] still multi-valued algebraic over s-line.

---

## 1. Reduced Hurwitz curve (moduli)

For \(\mathrm{Ni}(A_5,C_3^4)\) the **reduced Hurwitz curve** \(H^{\mathrm{rd}}\) is
irreducible of **genus 0** over \(\mathbb{Q}\), with infinitely many rational points
(Bailey–Fried; programme orbit size 18, single braid orbit — `GENUS_3A4_LOCK.md`).

As a curve,
$$H^{\mathrm{rd}}\ \cong\ \mathbb{P}^1$$
over \(\mathbb{Q}\). The map to moduli of 4 branch points is the cross-ratio
$$s:H^{\mathrm{rd}}\longrightarrow M_{0,4}\cong\mathbb{P}^1_s.$$
Thus the moduli space itself needs no higher-degree plane model: it **is**
\(\mathbb{P}^1\). What requires equations is the **cover / resolvent** over this base.

---

## 2. Cover normal form

Place branch values at \(\{0,1,\infty,s\}\). After automorphisms of the domain,
$$\varphi(y)=\frac{c\, y^3(y-1)(y-p_2)}{(y-r_1)(y-r_2)}
= \frac{N(y)}{D(y)},$$
with \(N=c\,y^3(y-1)(y-p_2)\) and \(D=y^2-\sigma y+\pi\)
(\(\sigma=r_1+r_2\), \(\pi=r_1 r_2\)).

Type \((3,1,1)\) at \(0\) and \(\infty\) is built in; type \((3,1,1)\) at \(1\) and \(s\)
is imposed by requiring \(\varphi-1\) and \(\varphi-s\) each to have a **triple root**
(denoted \(q\) and \(w\) respectively):
$$(\varphi-1)(q)=(\varphi-1)'(q)=(\varphi-1)''(q)=0,$$
$$(\varphi-s)(w)=(\varphi-s)'(w)=(\varphi-s)''(w)=0.$$

---

## 3. Explicit parameter formulae

Eliminating the linear variable \(c\) from the second-derivative equations yields
(for the \(q\)-chart; \(w\)-chart analogous):

$$c = -\frac{1}{q\bigl(6 p_2 q - 3 p_2 - 10 q^2 + 6 q\bigr)},$$

$$\sigma = \frac{q\bigl(8 p_2 q - 3 p_2 - 15 q^2 + 8 q\bigr)}{6 p_2 q - 3 p_2 - 10 q^2 + 6 q},$$

$$\pi = \frac{q^2\bigl(3 p_2 q - p_2 - 6 q^2 + 3 q\bigr)}{6 p_2 q - 3 p_2 - 10 q^2 + 6 q}.$$

Equating \(\sigma\) (resp. \(\pi\)) from the \(q\)- and \(w\)-charts produces two
polynomials \(F_1(p_2,q,w)=0\), \(F_2(p_2,q,w)=0\), quadratic in \(p_2\).
Their resultant in \(p_2\) factors as
$$-4 q^2 w^2 (q-w)^2\,(10qw-5q-5w+3)\,P(q,w),$$
with **physical component** the plane curve

### The eliminant curve \(P(q,w)=0\)

$$P(q,w)=20q^3w^3 - 40q^3w^2 + 27q^3w - 6q^3 - 40q^2w^3 + 73q^2w^2 - 45q^2w + 9q^2 + 27qw^3 - 45qw^2 + 26qw - 5q - 6w^3 + 9w^2 - 5w + 1$$

(equivalently in code: `P = 20*q**3*w**3 - 40*q**3*w**2 + 27*q**3*w - 6*q**3 - 40*q**2*w**3 + 73*q**2*w**2 - 45*q**2*w + 9*q**2 + 27*q*w**3 - 45*q*w**2 + 26*q*w - 5*q - 6*w**3 + 9*w**2 - 5*w + 1`).

Singular at \((q,w)=(1,1)\) (degenerate chart: `{'P': 0, 'Pq': 0, 'Pw': 0}`).

### Cross-ratio \(s\)

$$s = \frac{6 p_2 w^2 - 3 p_2 w - 10 w^3 + 6 w^2}{6 p_2 q^2 - 3 p_2 q - 10 q^3 + 6 q^2}.$$

### \(p_2\) on the physical branch

Solve \(F_1(p_2,q,w)=0\) (quadratic). The physical root is the branch with
\(p_2(1/\sqrt5,-1/\sqrt5)=-1\):

`p2 = (30*q**2*w - 15*q**2 + 30*q*w**2 - 37*q*w + 8*q - 15*w**2 + 8*w - sqrt(900*q**4*w**2 - 900*q**4*w + 225*q**4 - 1400*q**3*w**3 + 400*q**3*w**2 + 630*q**3*w - 240*q**3 + 900*q**2*w**4 + 400*q**2*w**3 - ...`

(`F1` poly: `16*p2**2*q*w - 8*p2**2*q - 8*p2**2*w + 3*p2**2 - 30*p2*q**2*w + 15*p2*q**2 - 30*p2*q*w**2 + 37*p2*q*w - 8*p2*q + 15*p2*w...`)

---

## 4. Degree-5 resolvent

For parameters \((c,p_2,\sigma,\pi)\) as above and fibre coordinate \(t\),
$$N(y)-t\,D(y)=0,\qquad N=c\,y^3(y-1)(y-p_2),\quad D=y^2-\sigma y+\pi.$$

Monic in \(y\) over the coefficient field:
$$f(y)=\frac{1}{c}\bigl(N(y)-t D(y)\bigr)= y^5-(1+p_2)y^4+p_2 y^3-\frac{t}{c}(y^2-\sigma y+\pi).$$

This is the **degree-5 resolvent** of the cover in the affine coordinate \(y\),
with coefficients in the function field of the \((q,w)\)-model
(and the free fibre parameter \(t\)).

- N: `c*y**3*(-p2 + y)*(y - 1)`
- D: `pi - sigma*y + y**2`
- monic/c: `-p2*y**4 + p2*y**3 + y**5 - y**4 - pi*t/c + sigma*t*y/c - t*y**2/c`

**Relation to \(s\) only.** Because \(H^{\mathrm{rd}}\to\mathbb{P}^1_s\) has degree \(>1\)
in this normal form (multiple covers / sheets for one \(s\)), a single-valued
\(f_s\in\mathbb{Q}(s)[y]\) is not expected without choosing a rational section of
\(H^{\mathrm{rd}}\to\mathbb{P}^1_s\). The **explicit** model is:
$$P(q,w)=0,\quad p_2=p_2(q,w),\quad s=s(p_2,q,w),\quadf_{q,w,t}(y)=\mathrm{monic}(N-tD).$$

---

## 5. Exact fibre at \(s=-1\) (over \(\mathbb{Q}(\sqrt5)\))

| param | value |
|-------|-------|
| \(s\) | \(-1\) |
| \(c\) | \(-\sqrt5\) |
| \(p_2\) | \(-1\) |
| \(r_1,r_2\) | \(\pm1/5\) |
| \(q,w\) | \(\pm1/\sqrt5\) |
| \(\sigma,\pi\) | \(0,\ -1/25\) |

$$N=-\sqrt5\, y^3(y^2-1),\qquad D=y^2-\frac1{25}.$$

Monic form over \(\mathbb{Q}(\sqrt5)(t)\):
$$y^5-y^3+\frac{t}{\sqrt5}\left(y^2-\frac1{25}\right)=0.$$

Norm to \(\mathbb{Q}(t)\) (eliminate \(\sqrt5\)):
$$5(y^5-y^3)^2-t^2\left(y^2-\frac1{25}\right)^2=0$$
(degree 10 over \(\mathbb{Q}(t)\), as expected from \([ \mathbb{Q}(\sqrt5):\mathbb{Q}]=2\)).

### Verification of formulae at this point

| check | pass |
|-------|:----:|
| P_at_known | **True** |
| p2_at_known | **True** |
| s_at_known | **True** |
| c_at_known | **True** |
| sigma_at_known | **True** |
| pi_at_known | **True** |

---

## 6. How to use the model

1. Pick \((q,w)\) on \(P(q,w)=0\) with \(q\neq w\), \(q,w\notin\{0,1\}\).
2. Set \(p_2\) to the physical root of \(F_1(p_2,q,w)=0\).
3. Compute \(c,\sigma,\pi\) from the \(q\)-chart formulae; \(s\) from the ratio above.
4. For each fibre parameter \(t\), form \(f(y)=\mathrm{monic}(N-tD)\).
5. Specialise to number fields; test Galois / BJ reduction / catalogue \(k\).

```bash
python explicit_3a4_equation.py
```

---

## 7. What is / is not closed

| item | status |
|------|--------|
| \(H^{\mathrm{rd}}\) genus 0 / \(\cong\mathbb{P}^1\) | **Locked** |
| Cover normal form \(\varphi=N/D\) | **Locked** |
| Eliminant \(P(q,w)=0\) | **Explicit polynomial** |
| \(c,\sigma,\pi,s\) as rational functions of \((p_2,q,w)\) | **Explicit** |
| \(p_2(q,w)\) | **Explicit quadratic formula** |
| Deg-5 resolvent \(N-tD\) | **Explicit** |
| Exact fibre \(s=-1\) over \(\mathbb{Q}(\sqrt5)\) | **Explicit** |
| Single-valued \(f_s\in\mathbb{Q}(s)[y]\) | **Blocked** for param of \(P(q,w)\): genus \(>0\) (`GENUS_P_QW.md`) |
| Geometric multi-\(k\) catalogue hit | **Open** (strict multi-seed false on \(P\)-point samples) |

_Generated by explicit_3a4_equation.py_