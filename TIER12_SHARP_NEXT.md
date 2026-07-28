# Sharp next options — Tier 1.2 follow-ups + 1.1 return

_Elapsed: 0.46s_

**Verdict:** Sharp next A/B/C (0.46s). A: under H_core disc□ rate=1.000 (Z-rate=1.000); H⇒disc□ via pure-even when Z exists. B: best T-only disc□ rate=0.000 (crit2_signal=False). C: no new identical-square subclass beyond known homogenisation/envelope. Necessity fragment: not obtained.

---

## A. Binary data chooses \(k\), then pure-even

**Map:** k = (-1)^{pop}(2*odd+1) / (2^{min(v2,4)}+1+(len mod 4)); then pure-even search for Z α,β

**Lemma shape:** If pure-even Z-coeffs exist for k=k(n), then disc□ always (classical). H ⇒ disc□ reduces to H ⇒ Z-coeff existence for that k.

| hypothesis \(\mathcal{H}\) | \(n\) | Z-coeff rate | disc□ rate | \(A_5\) rate | \(\mathcal{H}\Rightarrow\) disc□? |
|---------------------------|----:|-------------:|-----------:|------------:|:----------------------------------:|
| `H_all` | 70 | 1.000 | 1.000 | 0.000 | **True** |
| `H_small_v2` | 68 | 1.000 | 1.000 | 0.000 | **True** |
| `H_short_collatz` | 46 | 1.000 | 1.000 | 0.000 | **True** |
| `H_bounded_odd` | 66 | 1.000 | 1.000 | 0.000 | **True** |
| `H_core` | 44 | 1.000 | 1.000 | 0.000 | **True** |

**Reading.** When Z-coefficients exist on the pure-even slice for \(k=k(n)\), disc□ holds **identically** (classical). So \(\mathcal{H}\Rightarrow\mathrm{disc}\square\) is equivalent to \(\mathcal{H}\Rightarrow\) “pure-even Z model exists for \(k(n)\)”. On the scanned seeds, Z-finding succeeds often under `H_core` / `H_all` (see rates). This is a **composite lemma about \(F\)**, not necessity from HQCC axioms alone — pure-even is still inserted in the codomain.

### Sample images

| \(n\) | \(k(n)\) | Z? | disc□ | status |
|----:|----------|:--:|:-----:|--------|
| 1 | 3/2 | True | True | disc_sq_skip_gal |
| 2 | 3/4 | True | True | disc_sq_skip_gal |
| 3 | 7/5 | True | True | disc_sq_skip_gal |
| 4 | 3/7 | True | True | disc_sq_skip_gal |
| 5 | -11/3 | True | True | disc_sq_skip_gal |
| 6 | 7/3 | True | True | disc_sq_skip_gal |
| 7 | -15/2 | True | True | disc_sq_skip_gal |
| 8 | 1/4 | True | True | disc_sq_skip_gal |
| 9 | 19/5 | True | True | disc_sq_skip_gal |
| 10 | -11/5 | True | True | disc_sq_skip_gal |
| 11 | 23/4 | True | True | disc_sq_skip_gal |
| 12 | 7/6 | True | True | disc_sq_skip_gal |

---

## B. \(F\to T(\ldots)\) only — disc□ rate hunt (Crit-2 signal)

**Binary→T(a..f) only, no pure-even envelope insert**

| variant | disc□ rate | \(A_5\) rate | Crit-2 signal (>0.5)? |
|---------|----------:|-------------:|:---------------------:|
| `B1_M_deform_f` | 0.000 | 0.000 | **False** |
| `B2_a_power3` | 0.000 | 0.000 | **False** |
| `B3_embed_shape` | 0.000 | 0.000 | **False** |
| `B4_full_mix` | 0.000 | 0.000 | **False** |
| `B5_homog_binary_seed` | 0.000 | 0.000 | **False** |
| `B6_ef0` | 0.000 | 0.000 | **False** |
| `B7_embed_ell` | 0.000 | 0.000 | **False** |

**Best disc□ rate:** 0.000
**Any Crit-2 signal:** **False**

### Square samples (if any)

_No disc□ hits in the T-only variants on this seed set._

**Reading.** No template-only functor achieved disc□ rate → 1. Variants that look like base \(M\) or free mix stay odd (evenness obstruction). Snapshots on the non-BJ homogenisation family hit disc□ only when the **seed** was chosen even — again an external gate, not Crit-2 from \(T\) shape.

---

## C. Tier 1.1 return — identically square subclasses of \(T\)

No new multi-parameter identical-square subclass beyond the known homogenisation (1-param) and pure-even envelope (2-param, BJ-embed). Bilinear cuts tested still fail. Crit-2 necessity fragment not found.

### Findings

- **homogenised_no_x2 (known):** `{'name': 'homogenised_no_x2 (known)', 'beyond_BJ_embed': True, 'identical_square_in': 't (iff disc(seed) constant square)', 'disc_identity': True, 'free_continuous_params': 1}`
- **a=-ef, ce=-bf, d=b^2/e^2:** `{'name': 'a=-ef, ce=-bf, d=b^2/e^2', 'beyond_BJ_embed': True, 'identical_square': True, 'info': '0', 'note': 'reduces to -108 d^5 e^2 type factor — not square in free params'}`
- **full weighted homogenisation in T:** `{'name': 'full weighted homogenisation in T', 'beyond_BJ_embed': 'when p or q nonzero', 'identical_square_in': 't iff disc(seed) square', 'realizable': 'Yes when parameters solve the 4 coupling equations as polys in t; no-x2 (q=0,a=e=0) and BJ (p=q=0,d=0,a=-ef) are the clean cases. Generic (p,q,r,s) may need non-polynomial (rational) T-params.', 'free_continuous_params': 1}`
- **pure-even envelope (m,k):** `{'name': 'pure-even envelope (m,k)', 'beyond_BJ_embed': False, 'identical_square': True, 'free_continuous_params': 2, 'note': 'Largest known identical-square family overall; not beyond BJ-embed'}`

### Bilinear / sparse cuts

| cut | identically square? |
|-----|:-------------------:|
| a=3,d=-3,f=1 free b,c,e | **False** |
| a=-ef,f=1,d=e free b,c,e | **False** |
| a=-ef,f=1,d=e**2 free b,c,e | **False** |
| b=c,e=f,d=0 free a,b,e | **False** |
| e=f,a=-e**2,d=0 free b,c,e | **False** |

---

## Locked conclusions

| option | outcome |
|--------|---------|
| A. Binary \(k\) + pure-even under \(\mathcal{H}\) | Works as composite; disc□ when Z exists; **not** HQCC necessity |
| B. \(F\to T\) only, disc□→1 | **Not achieved** on tested variants |
| C. New identical-square subclass of \(T\) | **None** beyond known homogenisation / pure-even envelope |

**Necessity fragment:** still open.  
**Finished centre:** pure-even multi-\(k\) untouched.  
**Organising principle:** still explains generative efficiency, not force.

```bash
python tier12_sharp_next.py
```

_Generated by tier12_sharp_next.py_