# Tier 1.3 — Genus / parameterisation of \(P(q,w)=0\)

_Elapsed: 25.32s_

**Verdict:** Tier 1.3 P(q,w) (25.32s). deg=6, pa=10, irr over Q=True. Ordinary triple at (1,1) δ=3; singularities at infinity; genus estimate g≈1 (not 0). No rational param from (1,1)-pencil (residual deg 3). Q-points=4, quad/extra samples; k-samples=120; catalogue tags=['class_flip', 'classical']; strict multi-seed multi-k=False.

---

## 1. Curve model

Physical eliminant for the \(3A^4\) cover (from `EXPLICIT_3A4_EQUATION.md`):

$$P(q,w)=0\quad (\deg 6,\ \deg_q=\deg_w=3).$$

| quantity | value |
|----------|-------|
| Total degree | **6** |
| Arithmetic genus \(p_a=(d-1)(d-2)/2\) | **10** |
| Irreducible over \(\mathbb{Q}\) | **True** |
| Affine singularity | \((1,1)\), multiplicity **3** |
| Lowest form at \((1,1)\) | `(Q + W)*(Q**2 + Q*W + W**2)` |
| \(\delta\) if ordinary triple | **3** |
| Infinity locus \(z=0\) | `20*q**3*w**3` |
| Projection residual deg from \((1,1)\) | **3** |
| Genus estimate (affine+inf lower) | **1** |
| Rational parameterisation found | **False** |

### Genus conclusion

g > 0 likely (estimate ~1 or higher if inf less severe; not g=0 from ordinary triple alone: pa-3=7). Projection from (1,1) has residual degree 3 → not a birational param by lines through (1,1).

**Not genus 0** under the ordinary-singularity estimate: \(p_a-3=7\) from the affine triple alone; infinity adds further \(\delta\). A global rational parameterisation is **not** expected, and the \((1,1)\)-line pencil does **not** birationally parameterise the curve (residual degree 3).

Groebner of singular ideal: `['q**2 - 2*q - w**2 + 2*w', 'q*w - q + 5*w**2/4 - 7*w/2 + 9/4', 'w**3 - 3*w**2 + 3*w - 1']`

---

## 2. Points

### Rational points (den ≤ 10): **4** verified

| \(q\) | \(w\) |
|------|------|
| 0 | 1/2 |
| 1/2 | 0 |
| 1/2 | 1/2 |
| 1 | 1 |

### Quadratic / known special: **5**

- q=`sqrt(5)/5`, w=`-sqrt(5)/5` (Q)
- q=`-sqrt(5)/5`, w=`sqrt(5)/5` (Q)
- q=`1/2`, w=`1/2` (Q)
- q=`0`, w=`1/2` (Q)
- q=`1/2`, w=`0` (Q)

Numeric real samples used for \(k\): **40** (subsampled).

---

## 3. Resolvent path and \(k\) vs catalogue

For each smooth point \((q,w)\): physical \(p_2\) from \(F_1\), then \(c,\sigma,\pi,s\), fibres \(N-tD\), numeric BJ \(k=\beta/\alpha\).

| quantity | value |
|----------|------:|
| \(k\) samples (real) | 120 |
| Catalogue near-hits (tol 0.08) | 40 |
| Catalogue tags | ['class_flip', 'classical'] |
| Strict multi-seed multi-\(k\) (≥2 families) | **False** |

### Sample \(k\) rows

| \(q\) | \(w\) | \(s\) | \(t\) | \(k\) |
|------|------|------|----:|-----:|
| sqrt(5)/5 | -sqrt(5)/5 | -0.999999999999997 | -2 | 0.5876 |
| sqrt(5)/5 | -sqrt(5)/5 | -0.999999999999997 | 2 | -0.5876 |
| sqrt(5)/5 | -sqrt(5)/5 | -0.999999999999997 | 3 | -0.7422 |
| sqrt(5)/5 | -sqrt(5)/5 | -0.999999999999997 | 4 | -0.8052 |
| -sqrt(5)/5 | sqrt(5)/5 | -0.999999999999997 | -2 | -0.5876 |
| -sqrt(5)/5 | sqrt(5)/5 | -0.999999999999997 | 2 | 0.5876 |
| -2.5000 | 0.4667 | -0.01898947305668241 | -2 | -0.8104 |
| -2.5000 | 0.4667 | -0.01898947305668241 | -1 | -0.6826 |
| -2.5000 | 0.4667 | -0.01898947305668241 | 4 | 0.7998 |
| -2.3980 | 0.4659 | -0.02111271392785044 | -2 | -0.8109 |
| -2.3980 | 0.4659 | -0.02111271392785044 | -1 | -0.6829 |
| -2.3980 | 0.4659 | -0.02111271392785044 | 4 | 0.8010 |
| -2.2959 | 0.4649 | -0.023567661983789207 | -2 | -0.8113 |
| -2.2959 | 0.4649 | -0.023567661983789207 | -1 | -0.6831 |
| -2.2959 | 0.4649 | -0.023567661983789207 | 4 | 0.8023 |
| -2.1939 | 0.4640 | -0.02642244915602465 | -2 | -0.8117 |
| -2.1939 | 0.4640 | -0.02642244915602465 | -1 | -0.6833 |
| -2.1939 | 0.4640 | -0.02642244915602465 | 4 | 0.8037 |
| -1.9898 | 0.4619 | -0.033698447078494016 | -2 | -0.8123 |
| -1.9898 | 0.4619 | -0.033698447078494016 | -1 | -0.6834 |

### Catalogue near-hits

- q=sqrt(5)/5, s=-0.999999999999997, t=3: k=-0.7422 ≈ class_flip
- q=sqrt(5)/5, s=-0.999999999999997, t=4: k=-0.8052 ≈ class_flip
- q=-2.5000, s=-0.01898947305668241, t=-2: k=-0.8104 ≈ class_flip
- q=-2.5000, s=-0.01898947305668241, t=4: k=0.7998 ≈ classical
- q=-2.3980, s=-0.02111271392785044, t=-2: k=-0.8109 ≈ class_flip
- q=-2.3980, s=-0.02111271392785044, t=4: k=0.8010 ≈ classical
- q=-2.2959, s=-0.023567661983789207, t=-2: k=-0.8113 ≈ class_flip
- q=-2.2959, s=-0.023567661983789207, t=4: k=0.8023 ≈ classical
- q=-2.1939, s=-0.02642244915602465, t=-2: k=-0.8117 ≈ class_flip
- q=-2.1939, s=-0.02642244915602465, t=4: k=0.8037 ≈ classical
- q=-1.9898, s=-0.033698447078494016, t=-2: k=-0.8123 ≈ class_flip
- q=-1.9898, s=-0.033698447078494016, t=4: k=0.8068 ≈ classical
- q=-1.8878, s=-0.038369671889085355, t=-2: k=-0.8125 ≈ class_flip
- q=-1.8878, s=-0.038369671889085355, t=4: k=0.8086 ≈ classical
- q=-1.7857, s=-0.04395995954225305, t=-2: k=-0.8126 ≈ class_flip

---

## 4. Locked outcome (Tier 1.3)

| question | answer |
|----------|--------|
| Genus 0? | **No** (estimate \(g>0\); not a rational curve from this analysis) |
| Global rational param? | **Not found**; (1,1)-pencil residual deg 3 |
| Single-valued \(f_s\in\mathbb{Q}(s)[y]\) via param of \(P\)? | **Blocked** by \(g>0\) |
| Point enumeration + BJ \(k\) | **Done** on Q / Q(\(\sqrt5\)) / real samples |
| Geometric multi-\(k\) catalogue (≥2 multi-seed families) | **False** so far |

**Geometric options left** (non-blocking; arithmetic multi-\(k\) remains citable centre):

1. **Full blowup** → exact geometric genus of \(P\).  
2. **Degree-3 function field** over \(\mathbb{P}^1\) (accept coefficients algebraic in \(s\)).  
3. **Other Nielsen types** (\(2A\,3A^3\), \(2A^2\,3A^2\)) aimed at a genus-0 chart.  
4. Leave **Nielsen fusion** as open research.

**Priority note:** With Tier 1.2 and 1.3 done, roadmap priority returns to **Tier 1.1** (identical-square subclass of \(T\)) for necessity — see `RESEARCH_ROADMAP.md`.

```bash
python genus_p_qw.py
```

_Generated by genus_p_qw.py — Tier 1.3_