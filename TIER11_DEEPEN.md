# Tier 1.1 deepen — identical-square subclasses of \(T\)

_Elapsed: 105.34s_

**Verdict:** Tier 1.1 deepen (105.34s). Disc not square unrestricted (deg 12; all structural cuts odd). Seed scan hits: 12 identical-square (incl. pure-even p=0, weighted fixed seeds, reducible s≡0). Genuine non-quasihomog 2-param beyond BJ: **0**. Poly maps with identical square: 4 (expected: homog fixed even seed ± pure-even). Largest beyond BJ-embed: **dim 1** (`homogenised_no_x2_fixed_even_seed`). 3-cycles on homog A5 seeds: inherited (not forced by T). HQCC-axiom naming: **fails**. Necessity fragment: **not obtained**.

---

## 0. Scope

Largest natural subclass of \(T(a,b,c,d,e,f)\) on which
\(\operatorname{disc}(\chi_T)\) is a **square in the polynomial ring** of free
parameters — **beyond** pure-even envelope (BJ-embed) and the known
1-param homogenisation of fixed even seeds \(x^5+px^3+rx+s\).

Then: forced 3-cycles? HQCC-axiom naming?

---

## 1. Discriminant structure

- Total degree: **12**
- Unrestricted already square? **False** (cuts_not_all_square)
- Factor note: factored on a=0,e=0 (no-x2 slice)
- Number of irreducible factors on cut: **2**

### Square tests on structural cuts

| cut | identical square? | info |
|-----|:-----------------:|------|
| `unrestricted_sample_a` | **False** | content_neg:-4 |
| `a_e_0` | **False** | content_neg:-1 |
| `BJ_raw` | **False** | content_neg:-1 |
| `e_f_0` | **False** | content_neg:-1 |

| factor (preview) | mult | deg |
|-------------------|-----:|----:|
| `b` | 2 | 1 |
| `256*b**3*f**5 - 3125*b**2*c**4 + 2000*b**2*c**2*d*f**2 + 128*b**2*d**2*f**4 + 900*b*c**2*d**3*f + 16*b*d**4*f**3 + 108*c` | 1 | 8 |

**Note:** “disc is a square number” is **not** a polynomial equation on
\((a,\ldots,f)\). Identical-square families must make every odd-multiplicity
factor vanish or pair up under the constraint ideal.

Univariate analysis: `{'deg_in_a': 6, 'gcd_with_derivative_deg': 0, 'gcd_expr_preview': '1', 'note': 'Nontrivial gcd(Disc, ∂_a Disc) generates an ideal of singular locus in a, not the full even-monodromy locus (which is not cut out by polynomials alone).'}`

---

## 2. Known baseline families

### pure_even_envelope_BJ

`{'name': 'pure_even_envelope_BJ', 'beyond_BJ_embed': False, 'free_params': ['m', 'k'], 'dim': 2, 'identical_square': True, 'note': 'Largest known; classical; not HQCC-native necessity'}`

### homogenised_no_x2_in_T

`{'name': 'homogenised_no_x2_in_T', 'beyond_BJ_embed': True, 'when': 'p≠0', 'free_params': ['t'], 'seed_params': ['p', 'r', 's'], 'dim_continuous': 1, 'identical_square_in_t': True, 'note': 'Square in t iff disc(seed) constant square'}`

---

## 3. Search: 1-parameter even seeds \(\Rightarrow\) 2-param after homog

If \(p(u),r(u),s(u)\) are polynomials with
\(\operatorname{disc}(x^5+p x^3+r x+s)\) identically square in \(u\), and \(p\not\equiv 0\),
then the \(T\)-family
$$a=e=0,\ d=-p(u)t^2,\ b=1,\ f=-r(u)t^4,\ c=-s(u)t^5$$
has \(\operatorname{disc}=t^{20}\operatorname{disc}(\mathrm{seed}(u))\) identically square in
\((t,u)\) — a **2-parameter** family **beyond BJ-embed**.

**Raw identical-square seed hits (any class):** **12**
**Genuine 2-param beyond BJ (non-quasihomog, irreducible candidate):** **0**
**Weighted fixed-seed (still dim 1 after scale):** **2**
**Reducible \(s\equiv 0\):** **4**
**Classical pure-even \(p=0\):** **4**

### Classification of scan hits

| class | count | meaning |
|-------|------:|---------|
| genuine 2-param beyond BJ | **0** | would give new dim-2 in \(T\) with \(d\neq 0\) |
| weighted fixed seed | 2 | \(p=p_0 u^2,r=r_0 u^4,s=s_0 u^5\) = scale of constant seed |
| reducible \(s\equiv 0\) | 4 | \(\chi=x(\cdots)\); disc may square but not \(A_5\) source |
| pure-even BJ seed \(p=0\) | 4 | classical; not beyond BJ-embed |
| constant even seeds | 2 | 1-param after homog (known) |

**No genuine 2-param even-seed polynomial family** found beyond quasihomogeneous scales of fixed seeds / pure-even / reducible cases.

Illustrative false friends:

| \(p\) | \(r\) | \(s\) | class |
|------|------|------|-------|
| 6*u**2 | -7*u**4 | 8*u**5 | weighted_fixed_seed_dim1 |
| 6*u**2 | -7*u**4 | -8*u**5 | weighted_fixed_seed_dim1 |
| 1 | u**2 | 0 | reducible_s0 |
| 2 | u**2 | 0 | reducible_s0 |
| 3 | u**2 | 0 | reducible_s0 |
| 6 | u**2 | 0 | reducible_s0 |

---

## 4. Polynomial maps \(\mathbb{A}^r\to T\)-parameters

| map | identical square? | beyond BJ? | info |
|-----|:-----------------:|:----------:|------|
| scale_M | **False** | True | content_neg:-4 |
| scale_all_ternary | **False** | True | odd=[('u', 7), ('324*u**2 - 947*u + 1152', 1)] |
| homog_weights | **False** | True | odd=[('256*u**12 + 4997*u**8 - 884*u**4 + 108', 1)] |
| LSW_embed_path | **True** | False | square |
| flag_k_embed_cleared | **True** | False | square |
| homog_seed_6_7_8 | **True** | True | square |
| M_deform_a_only | **False** | True | odd=[('9*u + 4880', 1), ('6561*u**5 + 39366*u**3 + 266814000*u**2 + 30140105', 1)] |
| envelope_LSW_flag_path | **False** | None | odd=[('u**2 - 3125', 5), ('3125*u**6*v**4 + 50000*u**6*v**3 + 300000*u**6*v**', 1)] |
| homog_two_seeds_blend | **False** | None | odd=[('432*v**5 - 343*v**4 + 25200*v**3 - 19208*v**2 + 39', 1)] |
| true_envelope_cleared | **False** | None | odd=[('256*u**2 - 3125*v**4', 5), ('52428800000*u**6*v**4 - 1920000000000*u**4*v**8 + ', 1)] |
| homog_two_param_blend_fixed_seeds | **False** | None | odd=[('50000*v**4 - 100000*v**3 + 249288*v**2 - 199288*v ', 1)] |
| pure_even_envelope_symbolic_mk | **True** | None | square |

---

## 5. 3-cycles on beyond-BJ homogenisation

Seed \(x^5+6x^3-7x\pm 8\) (Gal \(A_5\)) homogenised in \(T\):

| seed \(s\) | \(t\) | disc□ | status | has 3-cycle census |
|----------:|----:|:-----:|--------|:------------------:|
| -8 | 1 | True | HIT_A5 | True |
| -8 | 2 | True | HIT_A5 | True |
| -8 | 3 | True | HIT_A5 | True |
| -8 | 5 | True | HIT_A5 | True |
| 8 | 1 | True | HIT_A5 | True |
| 8 | 2 | True | HIT_A5 | True |
| 8 | 3 | True | HIT_A5 | True |
| 8 | 5 | True | HIT_A5 | True |

**Not forced by \(T\):** 3-cycles are inherited from the seed (Hilbert), same as lattice search.

---

## 6. HQCC-axiom naming

| subclass | dim | HQCC name? | reason |
|----------|----:|:----------:|--------|
| pure-even envelope (BJ-embed) | 2 | None | Classical BJ pure-even; not forced by ternary/flux axioms |
| homogenised no-x2 (a=e=0) | 1 | None | a=0 contradicts base M (a=3); t-weights are homogenisation ansatz |
| homogenised_no_x2_fixed_even_seed | 1 | None | No HQCC axiom forces a=e=0 and even seed |

---

## 7. Largest subclass — locked comparison

| family | beyond BJ-embed? | free continuous params | identical square? | HQCC-native? |
|--------|:----------------:|-----------------------:|:-----------------:|:------------:|
| Pure-even envelope | No | **2** | Yes | No |
| Pure-even \(k\)-slice | No | 1 | Yes | No |
| Homog fixed even seed \(x^5+px^3+rx+s\) | **Yes** (if \(p\neq0\)) | 1 | Yes in \(t\) | No |
| Homog of 1-param even-seed family | **Yes** *if* genuine family exists | **2** | Yes in \((t,u)\) | No |
| Weighted scale of fixed seed | **Yes** (\(p\neq 0\)) | **1** (one scale) | Yes | No |

**Genuine 2-param beyond BJ found in scan:** **0**.

**Largest beyond BJ-embed established constructively:** `homogenised_no_x2_fixed_even_seed` (dim **1**).

**Largest overall (including BJ-embed):** `pure_even_envelope_BJ_embed` (dim **2**).

---

## 8. Conclusion (Tier 1.1 deepen)

1. **No Crit-2 necessity fragment:** nothing found that is simultaneously identically disc-square, beyond classical pure-even, forced 3-cycles, and HQCC-axiom-named.
2. **Dimension ceiling (locked):**
   - overall identical-square: **dim 2** = pure-even envelope (BJ-embed);
   - **beyond** BJ-embed: **dim 1** = homogenisation of a *fixed* even seed \(x^5+px^3+rx+s\) with \(p\neq 0\) (and its quasihomogeneous reparametrisations).
   - Scan found **no** genuine polynomial 1-param *family of seeds* with identically square disc and \(p\not\equiv 0\), \(s\not\equiv 0\) that would lift to a new dim-2 beyond BJ.
3. **3-cycles** remain operational (seed + Hilbert), not structural from \(T\).
4. **HQCC naming fails** for every identical-square subclass examined: each requires pure-even and/or \(a=e=0\) homogenisation ansätze foreign to base \(M\).
5. **Priority:** further 1.1 only with a new algebraic idea (literature parametric even quintics realised in \(T\)). Otherwise shift to Tier 2 (paper / catalogue invariants) or optional geometric leftovers.

```bash
python tier11_deepen.py
```

_Generated by tier11_deepen.py_