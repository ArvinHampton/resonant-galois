# Largest natural subclass of \(T\) beyond BJ-embed with \(\operatorname{disc}(\chi_T)\) identically square

_Elapsed: 13.71s_

**Verdict:** Beyond BJ-embed, the largest *natural* subclass with disc identically square in free parameters is 1-parameter homogenisation of fixed even-disc seeds of shape x^5+p x^3+r x+s (realizable in T with a=e=0). It does not exceed pure-even envelope dimension; 3-cycles are not structurally forced; HQCC-axiom naming fails (homogenisation + a=e=0 are extra ansätze). Criterion 2 necessity remains open.

---

## 0. Meaning of “identically square”

After restricting parameters of \(T\) to a family with free coordinates \(u_1,\ldots,u_k\in\mathbb{Q}\), the discriminant \(\operatorname{disc}(\chi_T)\) must be a **square in the polynomial ring** \(\mathbb{Q}[u_1,\ldots,u_k]\) (or zero), not merely a square number for many integer specialisations.

BJ-embed alone (\(d=0\), \(a=-ef\)) does **not** make disc identically square in free \((b,c,e,f)\): one recovers \(256\alpha^5+3125\beta^4\), square only on the pure-even subvariety of \((\alpha,\beta)\).

---

## 1. Coordinate / sparse cuts (no identical square)

| subclass | identically square poly? | beyond BJ-embed? |
|----------|:------------------------:|:----------------:|
| unrestricted T | **False** | True |
| BJ-embed d=0,a=-ef (raw free bcef) | **False** | False |
| a=-ef only | **False** | True |
| d=0 only | **False** | True |
| e=f=0 | **False** | True |
| b=e=0 | **False** | True |
| c=f=0 | **False** | True |
| b=c=0 | **False** | True |
| a=e=0 | **False** | True |
| a=e=0,f free poly family later | **False** | True |
| LR only a=b=c=e=0 | **True** | True |
| UL only b=e=f=0 | **False** | True |

Degenerate case \(\operatorname{disc}=0\) (e.g. LR-only \(a=b=c=e=0\)) is a square but \(\chi\) is reducible — not a source of \(A_5\).

---

## 2. Homogenisation family in \(T\) **beyond** BJ-embed

### Construction

Set
$$a=0,\quad e=0,\quad d=-p\,t^2,\quad b=1,\quad f=-r\,t^4,\quad c=-s\,t^5.$$

Then (verified \(\chi_T\) match = **True**):

$$\chi = x^5 + p\, t^2 x^3 + r\, t^4 x + s\, t^5.$$

### Discriminant identity (verified = **True**)

$$\operatorname{disc}(\chi) = t^{20}\,\operatorname{disc}(x^5+p x^3+r x+s).$$

Hence disc is **identically a square in \(t\)** whenever the seed \(x^5+p x^3+r x+s\) has square (constant) discriminant.

**Beyond BJ-embed:** when \(p\neq 0\), one has \(d\neq 0\), so the family is **not** contained in \(d=0\).

| comparison | free continuous params |
|------------|------------------------:|
| Pure-even envelope (BJ-embed) | **2** \((m,k)\) |
| Pure-even fixed-\(k\) slice | **1** |
| This homogenisation (fixed even seed) | **1** \((t)\) |

**Largest beyond BJ-embed with identical-square disc:** these 1-parameter homogenisations (and the same shape with \(b=b_0\neq 0\) rescaling). No natural **2-parameter** polynomial family beyond BJ-embed was found with disc a square in \(\mathbb{Q}[u,v]\).

---

## 3. Forced 3-cycles?

On the homogenisation family, Galois behaviour is that of the seed via Hilbert specialisation: if \(\mathrm{Gal}(\mathrm{seed}/\mathbb{Q})=A_5\) and type \((3,1,1)\) appears, the same holds for many \(t\). This is **not** a structural force from \(T\)'s shape alone — it is inherited from the seed choice (disc gate + cycle gate), same as unrestricted lattice search.

Non-BJ seeds with \(p\neq 0\), disc□ found in small scan: **2** (A5 status among them: **2**).

| \(p\) | \(r\) | \(s\) | disc | status |
|----:|----:|----:|-----:|--------|
| 6 | -7 | -8 | 168792064 | HIT_A5 |
| 6 | -7 | 8 | 168792064 | HIT_A5 |

### Sample \(t\)-specialisations (disc remains □)

**Seed** \((p,r,s)=(6,-7,-8)\), HIT_A5:
- t=2: disc□=True, status=HIT_A5
- t=3: disc□=True, status=HIT_A5
- t=5: disc□=True, status=HIT_A5

**Seed** \((p,r,s)=(6,-7,8)\), HIT_A5:
- t=2: disc□=True, status=HIT_A5
- t=3: disc□=True, status=HIT_A5
- t=5: disc□=True, status=HIT_A5

---

## 4. HQCC-axiom naming?

| subclass | HQCC-native? | reason |
|----------|:------------:|--------|
| BJ-embed + pure-even | **No** | classical pure-even ansatz |
| Homogenised no-\(x^2\) (\(a=e=0\)) | **No** | homogenisation + entry zeroing not forced by ternary/flux axioms (base \(M\) has \(a=3\neq 0\)) |

a=e=0 is a structural zeroing of two entries — can be called 'decoupled companion coupling' but is not forced by ternary/flux axioms alone (base M has a=3, e=0 actually — e=0 in M!). M has a=3≠0. So a=0 is NOT the base HQCC matrix. f,c proportional to t-powers is homogenisation ansatz, not HQCC axiom.

**No subclass found that is both (i) identically disc-square in free params beyond BJ-embed pure-even/envelope dimension, and (ii) named purely by HQCC axioms without an extra evenness/homogenisation ansatz.**

---

## 5. Locked conclusion (Criterion 2)

1. **No** free multi-parameter cut of \(T\) (beyond degeneracy) makes \(\operatorname{disc}(\chi_T)\) a square polynomial without further evenness conditions.
2. **Beyond BJ-embed**, the largest natural identical-square families are **1-parameter homogenisations** of fixed even-disc seeds of shape \(x^5+px^3+rx+s\) inside \(T\) (\(a=e=0\)).
3. These do **not** beat the pure-even envelope in parameter count, do **not** force 3-cycles structurally, and are **not** HQCC-axiom native.
4. Therefore this probe **does not** produce a Criterion-2 necessity theorem.
5. Pure-even multi-\(k\) remains the finished arithmetic centre; necessity remains open.

```bash
python t_subclass_identical_square.py
```

_Generated by t_subclass_identical_square.py_