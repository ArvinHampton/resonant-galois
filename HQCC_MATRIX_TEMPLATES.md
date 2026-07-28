# HQCC matrix templates — explicit exploration

_Elapsed: 0.93s_

**Verdict:** HQCC matrix templates (0.93s). chi_T formula identity=True. Base M: disc not square, Gal=S5. Example T(3,0,0,-3,1,3): A5. BJ-embed → classical pure-even (thin subclass). Templates alone do NOT force disc□ / necessity. Verification=PASS.

Links: necessity target `NECESSITY_THEOREM.md` (Criterion 2); pure-even theory `RESOLUTION_PATH.md`.

---

## 1. Base structural matrix \(M\)

$$M=\begin{pmatrix}0&1&0&0&0\\0&0&1&0&0\\3&0&0&80&0\\0&0&0&0&1\\61&0&0&-3&0\end{pmatrix}$$

Specialisation of \(T\): \((a,b,c,d,e,f)=(3,80,61,-3,0,0)\).

| item | value |
|------|-------|
| Characteristic polynomial | `x**5 + 3*x**3 - 3*x**2 - 4889` |
| Discriminant | `1781436218361244736` — **not** a square |
| Galois group | **S5TransitiveSubgroups.S5** (`odd_monodromy`) |
| Flux fingerprint | \(4889=4880+3^2\), \(80=4880/61\) |

### Block reading

| block | role |
|-------|------|
| UL \(3\times3\) | Companion-like ternary block (entry \(3\)) |
| Couplings | \(80\) and \(61\) (flux / puncture integers) |
| LR \(2\times2\) | \(\begin{pmatrix}0&1\\-3&0\end{pmatrix}\) |

**Evenness obstruction in concrete form:** full ternary + flux structure, yet disc odd \(\Rightarrow S_5\). 3-cycles present; monodromy still **odd**.

---

## 2. Structural template \(T(a,b,c,d,e,f)\)

$$T(a,b,c,d,e,f)=\begin{pmatrix}0&1&0&0&0\\0&0&1&0&0\\a&0&0&b&e\\0&0&0&0&1\\c&f&0&d&0\end{pmatrix}$$

### Characteristic polynomial (exact, verified)

$$\chi_T=x^5 - d\,x^3 - (a+ef)\,x^2 - (bf+ce)\,x + (ad-bc).$$

Symbolic identity holds: **True**

Parameters from the resonant / model lattice \(\{3,9,27,61,80,243,539,\ldots\}\) and short combinations.

### Deformation result

Restricting to the lattice and gating on square disc produces multiple explicit \(A_5\) realisations (historical: 14+). Verified example:

$$T(3,0,0,-3,1,3)\quad\Rightarrow\quad\chi=x^5+3x^3-6x^2-9,\quad\operatorname{disc}=3470769=1863^2,\quad\mathrm{Gal}=A_5.$$

Runtime check: status=`HIT_A5`, disc□=`True`.

**Lesson:** same template shape + parameters so disc is square \(\Rightarrow A_5\). **Structure alone does not force evenness; the disc gate does.**

---

## 3. BJ-embed subclass (templates \(\leftrightarrow\) pure-even theory)

Impose
$$d=0,\qquad a=-ef.$$
Then \(x^2\) terms cancel and
$$\chi=x^5-(bf+ce)\,x-bc,$$
which is Bring–Jerrard:
$$\alpha=-(bf+ce),\qquad \beta=-bc.$$

Identity verified: **True** (`-b*c - b*f*x - c*e*x + x**5`).

On this thin subclass, the **pure-even theory** applies: fix \(k=\beta/\alpha\) and run the classical envelope. Even monodromy becomes an **identity**, not a search outcome.

| property | status |
|----------|--------|
| Disc identically square on pure-even rays | **Yes** (classical BJ) |
| Forced by full unrestricted \(T\) | **No** |
| Native HQCC labelling of \((b,c,e,f)\) | Still an **ansatz** on top of the template |

---

## 4. Degree-6 enlargement \(T_6\)

- Base: `x^6 + 3*x^4 - 3*x^2 - 4889`
- Gal (approx): S4 x C2 (disc not square)
- Example A6: `x^6 - 3*x^4 + 9*x^2 ± 18*x + 9 (historical catalogue)`
- Sparse ternary specialisations produced verified A6 after disc gate; structure alone does not force even monodromy.

---

## 5. What the templates do and do not give

| claim | verdict |
|-------|---------|
| Templates encode order-3 / flux data | **Yes** |
| Templates produce many \(A_5\) (and some \(A_6\)) after disc gate | **Yes** |
| Templates force disc □ | **No** (base \(M\), base \(T_6\)) |
| BJ-embed recovers pure-even arithmetic | **Yes**, on a thin subclass |
| Templates alone yield a necessity theorem | **No** |

---

## 6. Implication for Criterion 2 (structural axioms)

Any axiom list that claims
> resonant matrix shape \(\Rightarrow\) alternating monodromy
must be **strictly stronger** than membership in \(T(a,b,c,d,e,f)\) (or \(T_6\)).

The minimal strengthening that is **known to work** is:
1. restrict to the BJ-embed (\(d=0\), \(a=-ef\)), and
2. impose a pure-even condition on \((\alpha,\beta)\),

which is exactly the **classical pure-even theory already finished** — **not** a new necessity theorem native to HQCC.

See **`NECESSITY_THEOREM.md`**.

---

## 7. Probe: larger subclasses with identically square disc?

Light scans (not exhaustive):

### `BJ-embed (d=0, a=-ef)`

- **disc_expression:** `-256*b**5*f**5 + 3125*b**4*c**4 - 1280*b**4*c*e*f**4 - 2560*b**3*c**2*e**2*f**3 - 2560*b**2*c**3*e**3*f**2 - 1280*b*c**4*e**4*f - 256*c**5*e**5`
- **identically_square:** `True`
- **note:** `Disc is the BJ form 256 A^5+3125 B^4; square iff (A,B) on pure-even locus, not for all b,c,e,f.`
- **forces_square_for_all_params:** `False`

### `e=f=0`

- **disc_simplified:** `-(a*d - b*c)*(108*a**5 - 216*a**3*d**3 + 2700*a**2*b*c*d**2 - 5625*a*b**2*c**2*d + 108*a*d**6 + 3125*b**3*c**3 - 108*b*c*d**5)`
- **sample_square_rate:** `42/756`
- **forces_square_for_all_params:** `False`

### `b=e=0 sparse`

- **sample_square_rate:** `36/108`
- **sample_hits:** `[{'a': 3, 'c': 0, 'd': -3, 'f': 0, 'disc': 3779136}, {'a': 3, 'c': 0, 'd': -3, 'f': 3, 'disc': 3779136}, {'a': 3, 'c': 0, 'd': -3, 'f': -3, 'disc': 3779136}, {'a': 3, 'c': 0, 'd': -3, 'f': 1, 'disc': 3779136}, {'a': 3, 'c': 3, 'd': -3, 'f': 0, 'disc': 3779136}]`
- **forces_square_for_all_params:** `False`

**Probe executed:** `T_SUBCLASS_IDENTICAL_SQUARE.md` / `t_subclass_identical_square.py`.

**Result:** Beyond BJ-embed, largest natural identical-square families are **1-param homogenisations** of even-disc seeds \(x^5+px^3+rx+s\) inside \(T\) (\(a=e=0\), \(d=-pt^2\), …). They do **not** exceed pure-even envelope dimension, do **not** force 3-cycles structurally, and are **not** HQCC-axiom native. **Criterion 2 necessity still open.**

```bash
python hqcc_matrix_templates.py
```

_Generated by hqcc_matrix_templates.py_