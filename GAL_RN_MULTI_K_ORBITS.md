# Does \(\mathrm{Gal}(R_n/\mathbb{Q})\) act so multi-\(k\) paths are unions of cosine orbits?

_Elapsed: 6.56s_

**Answer: NO** (for the programme’s arithmetic multi-\(k\) paths).

NO — multi-k paths (arithmetic envelope) are not unions of cosine orbits. Gal(R_n/Q) fixes Q-points (catalogue k are singletons). Cosine orbits are finite/discrete; non-constant paths have infinite rational image. Weaker true facts: (i) Q-paths are Gal-invariant as schemes; (ii) finite geometric specialisation sets over R_n can be unions of cosine orbits; (iii) Gal acts on path parameters so conjugate u hit conjugate cosine k when the path crosses the cosine locus.

---

## 0. Definitions

| object | definition |
|--------|------------|
| \(R_n\) | \(\mathbb{Q}(2\cos 2\pi/n)=\mathbb{Q}(\zeta_n)^+\) |
| \(\mathrm{Gal}(R_n/\mathbb{Q})\) | \(\cong(\mathbb{Z}/n\mathbb{Z})^\times/\{\pm1\}\), order \(\varphi(n)/2\) |
| action | \(\sigma_a:\ 2\cos\frac{2\pi}{n}\mapsto 2\cos\frac{2\pi a}{n}\) |
| multi-\(k\) path | rational curve \((m(u),k(u))\) in the pure-even envelope joining distinct ratio classes |
| cosine orbit | \(\mathrm{Gal}\cdot\bigl(2\cos\frac{2\pi p}{n}\bigr)\) — **finite** subset of \(R_n\) |

### Group orders (proxies)

| \(n\) | \([R_n:\mathbb{Q}]\) | \(\#\mathrm{Gal}\) |
|------|------------------:|---------------:|
| 5 | 2 | 2 |
| 7 | 3 | 3 |
| 11 | 5 | 5 |
| 15 | 4 | 4 |

---

## 1. Galois action on rational \(k\)

For every n, Gal(R_n/Q) acts trivially on Q. Hence every catalogue / rational pure-even ratio k ∈ Q is a Gal(R_n/Q)-fixed point: its orbit is {k}.

A multi-k path whose image lies in Q (e.g. linear k(u)∈Q(u) with rational endpoints) is a union of Gal-orbits, but each orbit is a singleton rational — not a cosine orbit of size >1.

**Catalogue \(k\)** (all in \(\mathbb{Q}\)):
`['-4', '4', '-8/5', '8/5', '4/5', '-4/5', '-12/5', '12/5', '-16/5', '16/5']`

Each has \(\mathrm{Gal}(R_n/\mathbb{Q})\)-orbit of size **1**. They are not
non-trivial cosine orbits.

---

## 2. Cosine orbits are finite and discrete

Example orbits of \(2\cos(2\pi/n)\):

- **n=5**: orbit size **2** — `['2cos(2π·1/5)', '2cos(2π·2/5)']`
- **n=7**: orbit size **3** — `['2cos(2π·1/7)', '2cos(2π·2/7)', '2cos(2π·3/7)']`
- **n=11**: orbit size **5** — `['2cos(2π·1/11)', '2cos(2π·2/11)', '2cos(2π·3/11)', '2cos(2π·4/11)', '2cos(2π·5/11)']`

A non-constant path \(k(u)\in\mathbb{Q}(u)\) takes **infinitely many** distinct
values at rational \(u\). A finite union of cosine orbits is **finite**.
Hence:

> **Obstruction.** A non-constant multi-\(k\) path cannot equal a union of
> cosine orbits as a set of field elements.

Details: A non-constant multi-k path k: P¹ ⇢ P¹ over Q has infinitely many distinct values k(u) for u ∈ Q (or u ∈ R̄). Each cosine orbit is finite. A finite union of cosine orbits is finite. Therefore a non-constant path image cannot equal a union of cosine orbits as sets of field elements.

---

## 3. Programme multi-\(k\) paths

| path | \(k(u)\) | endpoints | cosine orbit? | union of cosine orbits? |
|------|----------|-----------|:-------------:|:-----------------------:|
| flag_classical | `-8/5 + u·(4/5--8/5)` | ['-8/5', '4/5'] | **False** | **False** |
| flag_lsw | `-8/5 + u·(-4--8/5)` | ['-8/5', '-4'] | **False** | **False** |
| classical_lsw | `4/5 + u·(-4-4/5)` | ['4/5', '-4'] | **False** | **False** |

### Why not

Image of the path over Q is an infinite set of rationals (or a line segment in P¹(Q)), while every cosine orbit is finite of size ≤ φ(n)/2. An infinite set of rationals cannot be a finite union of finite cosine orbits unless the path image is finite — which a non-constant rational path is not when evaluated on infinitely many u.

Catalogue multi-k paths take values in Q. Non-rational cosine values 2cos(2π p/n) ∉ Q for n>2 (except degenerate cases). So path values are not cosine elements, and cosine orbits of size >1 consist of irrationals not hit by Q-paths.

**Catalogue vs cosine numeric collisions (n=5,7,11,15):** **none**

- n=5: hits=`[]`
- n=7: hits=`[]`
- n=11: hits=`[]`
- n=15: hits=`[]`

---

## 4. What \(\mathrm{Gal}(R_n/\mathbb{Q})\) *does* do to paths

Gal(R_n/Q) acts on R_n-points of the envelope. For a path defined over Q, σ·(m(u),k(u)) = (m(u),k(u)) when m,k ∈ Q(u). The path is fixed as a scheme over Q. Gal does not permute distinct rational k along the path into each other (each is fixed).

Weaker true statement: if a set S ⊂ P¹ of k-values is defined over Q (stable under Gal(Q̄/Q), hence under Gal(R_n/Q) after embedding), then S is a union of Gal(Q̄/Q)-orbits. For S ⊂ Q this is a union of singletons. This is Galois descent, not a cosine constraint.

### Worked example: \(R_5=\mathbb{Q}(\sqrt5)\), path flag↔classical

- Cosine orbit: `{'2cos(2π/5)': '-1/2 + sqrt(5)/2', '2cos(4π/5)': '-sqrt(5)/2 - 1/2', 'orbit_size': 2}`
- Path: `k(u) = -8/5 + (12/5)u`
- Hits cosine at u = `11/24 + 5*sqrt(5)/24`
- Hits conjugate cosine at u = `11/24 - 5*sqrt(5)/24`

The arithmetic path (over Q) can pass through a cosine value at a non-rational u ∈ R_5. The Gal-conjugate cosine is hit at the conjugate u. So: the path is NOT a union of cosine orbits; rather, the pair {u, σu} of parameters maps to a cosine orbit {k, σk}. That is Gal acting on the parameter of the path, not the path being made of cosine orbits.

**Picture:**
$$\mathrm{Gal}\curvearrowright u\quad\Longrightarrow\quadk(u)\mapsto k(\sigma u)=\sigma\bigl(k(u)\bigr)\quad\text{when }k(u)\in R_n,$$
for the linear path with \(\mathbb{Q}\)-coefficients. The path is a
\(\mathbb{Q}\)-curve; Gal permutes the **parameters** where it meets a
cosine locus, not the path into a cosine orbit.

---

## 5. When could a *weakened* cosine-orbit statement hold?

If the geometric multi-k locus is a finite Gal(R_n/Q)-stable set of special fibres (not a positive-dimensional path over Q), and those k are forced into the cyclotomic real locus by branch geometry, then that finite set can be a union of cosine orbits. That is a discrete specialisation set, not an arithmetic multi-k path in the envelope.

| scenario | multi-\(k\) path = ∪ cosine orbits? |
|----------|:-----------------------------------:|
| Arithmetic envelope paths over \(\mathbb{Q}\) | **No** |
| Finite Gal-stable set of geometric specialisations in cosine locus | **Possible** (discrete, not a path) |
| Path over \(R_n\) with image inside one cosine orbit | **No** (orbit finite, path infinite unless constant) |
| Path meets cosine locus at a Gal-orbit of parameters | **Yes as incidence**, not as path=orbit |

---

## 6. Locked answer

> **Does \(\mathrm{Gal}(R_n/\mathbb{Q})\) act so that multi-\(k\) paths are
> unions of cosine orbits?**

**No.**

1. Cosine orbits are finite; non-constant multi-\(k\) paths are not.
2. Programme paths take values in \(\mathbb{Q}\), fixed pointwise by
   \(\mathrm{Gal}(R_n/\mathbb{Q})\); each rational \(k\) is a singleton orbit,
   not a non-trivial cosine orbit.
3. Catalogue \(k\) do not coincide with \(2\cos(2\pi p/n)\) for proxy \(n\).
4. True related facts:
   - paths over \(\mathbb{Q}\) are Gal-invariant as schemes;
   - Gal acts on parameters \(u\in R_n\) along the path;
   - finite geometric \(k\)-sets cut out by cosine branch constraints
     *can* be unions of cosine orbits — that is a different object from
     arithmetic multi-\(k\) paths.

```bash
python gal_rn_multi_k_orbits.py
```

_Generated by gal_rn_multi_k_orbits.py_