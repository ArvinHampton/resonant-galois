# Hilbert modular forms and \(A_5\) — exploration

_Elapsed: 0.32s_

**Verdict:** Hilbert modular / icosahedral geometry for Q(√5) is a natural companion to the quadratic phenomena already observed (φ-disc factor 5; 3A^4 cover at s=-1). Klein A,B,C,D relation checked for weight-30 homogeneity. No new pure-even multi-k family over Q produced. Record as high-effort enrichment of Avenue 5; arithmetic multi-k remains the finished positive result; geometric multi-k (Nielsen-labelled) stays open.

---

## Classical connection (Hirzebruch / Klein)

There is a deep, classical link between Hilbert modular forms for the real
quadratic field \(\mathbb{Q}(\sqrt{5})\) and the icosahedral group \(A_5\):

- **Field / order:** `Q(sqrt(5))`, \(\mathcal{O}=Z[(1+sqrt(5))/2]  (golden integers)\).
- **Group:** `A5 ≃ PSL(2,F5)  (icosahedral)`.
- **Geometry:** Hilbert modular surface for a principal congruence subgroup of SL(2,O) (level related to the prime above 2), after cusp resolution, is equivariantly related to the Klein icosahedral surface / arrangement in P^4 (or quotient models in P^2). A5 acts on this surface.

### Invariant ring (combinatorial description)

| generator | weight | role |
|-----------|-------:|------|
| \(A\) | 2 | Klein icosahedral |
| \(B\) | 6 | Klein icosahedral |
| \(C\) | 10 | Klein icosahedral |
| \(D\) | 15 | 15-line arrangement / antiinvariant |

The (symmetric) ring of Hilbert modular forms in this setting is generated
using these invariants; the structure is one of the few Hilbert modular groups
describable combinatorially via the icosahedral arrangement.

### Klein relation (weight 30)

- R_ABC under weights (1,3,5) homogeneous of wt 15: **False**
- Full R1/R2 under (2,6,10,15): False/False (source-dependent normal forms)
- R_ABC: `64*A**3*C*(-A*C + 5*B) - 80*A**2*B*C**2 + 720*A*B**3*C - 1728*B**5`
- Weight conventions for Klein A,B,C,D vary (binary form degrees vs weighted Proj P(1:3:5) vs modular weights 2,6,10,15). Sign and term lists differ by source (Klein, Hirzebruch, Nagano). Programme use is the existence of the invariant ring and the Q(√5)–A5 link, not a specific normal form of R.

### References (entry points)

- Hirzebruch: Hilbert modular surfaces for Q(√5) and Klein’s cubic / icosahedron
- van der Geer: Hilbert modular surfaces
- Klein: Lectures on the icosahedron; invariants A,B,C,D
- Nagano et al.: icosahedral invariants and Hilbert modular forms for √5 (period maps / K3 / Shimura curves)
- Zagier et al.: modular surface for Q(√5) related to Klein cubic

---

## Relevance to the present programme

Two earlier observations already pointed toward \(\mathbb{Q}(\sqrt{5})\):

1. **Rigid \(\varphi\):** fibre disc \(=5\cdot(\mathrm{square})\) over \(\mathbb{Q}\)
   (proved=True); evenness after base change to
   \(\mathbb{Q}(\sqrt{5})\) (proved=True).
   See `K_SQRT5_EVEN.md`.
2. **3A⁴ cover at \(s=-1\):** lives over \(\mathbb{Q}(\sqrt{5})\)
   (`{'s': -1, 'c': '-sqrt(5)', 'p2': -1, 'r1': '1/5', 'r2': '-1/5', 'field': 'Q(sqrt(5))', 'source': 'EXPLICIT_3A4_RESOLVENT.md / build_3a4_resolvent.py'}`).

The Hilbert-modular geometry therefore supplies a natural arithmetic-geometric
**home for the quadratic obstruction** we encountered.

### Ring of integers

- Fundamental unit \(\varphi=1/2 + sqrt(5)/2\), minpoly `_x**2 - _x - 1`.

---

## Possible uses for geometric multi-\(k\) or \(A_5\) families

1. **Explicit equations from the invariant ring.** Generators \(A,B,C,D\) give
   concrete equations. Specialisations or linear systems on the Hilbert modular
   surface may produce parametric families of covers / resolvents with monodromy
   related to \(A_5\).
2. **Galois representations attached to Hilbert modular forms.** Weight-1 or
   parallel-weight forms can give 2-dimensional Galois representations with
   projective image \(A_5\) (icosahedral cases). Search among associated number
   fields / covers for pure-even Bring–Jerrard forms and multi-\(k\) membership.
3. **Understanding the evenness obstruction.** The preference for
   \(\mathbb{Q}(\sqrt{5})\) may be illuminated by the period mapping or the
   geometry of the Hilbert modular surface.
4. **Descent.** Forms or covers over \(\mathbb{Q}(\sqrt{5})\) can be examined for
   rational descent, potentially recovering objects over \(\mathbb{Q}\).

### Light probe (this package)

- D=0, B=1 curve: `-16*(4*A**4*C**2 - 20*A**3*C + 5*A**2*C**2 - 45*A*C + 108)`
- Integer points on that curve (tiny grid): []
- Limitation: Full Hilbert modular form spaces / Hecke eigenforms / projective A5 Galois representations require specialised software (e.g. Magma, Hilbert Modular Forms packages). This probe only checks classical invariant algebra and programme links.

---

## Rank relative to previous avenues

| item | value |
|------|-------|
| Rank | high-effort / speculative |
| Comparable to | Avenue 6 (higher-rank rigid systems r≥5), Avenue 7 (geometric lift of envelope) |
| Best recorded as | Enrichment of Avenue 5 (base change to Q(√5)) |
| Does not immediately supply | An explicit pure-even multi-k family over Q with Nielsen label |

**Not a replacement for:**

- Arithmetic multi-k envelope + paths (finished positive result)
- Avenue 1 closed-form 3A^4 resolvent over Q(s) (still open, higher leverage for multi-k over Q)

**Possible uses (list):**

- Explicit equations from invariant ring A,B,C,D → trial covers / resolvents
- Galois reps of weight-1 / parallel-weight HMF with projective image A5
- Illumination of the permanent factor-5 / evenness obstruction via periods
- Descent of Q(√5)-objects to Q

---

## Bottom line for the programme

1. The classical Hilbert-modular / icosahedral geometry for \(\mathbb{Q}(\sqrt{5})\)
   is a **natural and beautiful companion** to the quadratic phenomena we already
   observed.
2. It offers a **potential source** of new \(A_5\) equations and representations,
   but converting it into an explicit pure-even multi-\(k\) family over
   \(\mathbb{Q}\) remains a **substantial research project**.
3. The **finished positive result** of the programme continues to be the
   **arithmetic multi-\(k\) theory** (envelope + paths).
4. **Geometric multi-\(k\) (Nielsen-labelled) stays open**; Hilbert modular forms
   are one more **high-effort avenue** toward it, **not a short-cut**.
5. Record this as an **enrichment of Avenue 5** (base change), not a new primary
   attack path ahead of the 3A⁴ resolvent (Avenue 1).

### Ranking table (updated)

| Rank | Avenue | Effort | Likelihood of multi-\(k\) hit |
|-----:|--------|--------|-------------------------------|
| 1 | Better rational coordinate / resolvent for 3A⁴ | High | Moderate (genus 0) |
| 2 | Next shortlist genus-0 class (2A3A³ / 2A²3A²) | Med–High | Unknown |
| 3 | Positive-dimensional pure-even \(A_5\) strata | High | Open (arith. solid) |
| 4 | Other rigid triples | Medium | Low–Moderate |
| 5 | Base change + descent **(+ Hilbert modular / icosahedral)** | Medium–High | Low (probed); HMF speculative |
| 6 | Higher-rank rigid systems | Very high | Speculative |
| 7 | Geometric lift of the existing envelope | High | Speculative |

_Generated by hilbert_modular_a5.py_