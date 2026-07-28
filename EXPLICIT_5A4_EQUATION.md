# Explicit equation for Ni(A₅, 5A⁴)

_Elapsed: 4.8s_

**Verdict:** G3c explicit 5A^4 (4.8s). Nielsen lock: 600/600 lift+1, orbits=1 size=[10]. Explicit model: pure-even path_flag_classical over Q(t) (disc□ identity=True). Monodromy re-check 5A^4=True. Hilbert multi-k FC=True FL=True. Fusion (arithmetic path = geometric 5A^4 type): True.

---

## 0. Theoretical lock (settled)

James / Magaard–Shpectorov–James: for type `(5A,5A,5A,5A)` there is
**exactly one** braid orbit, and that orbit carries lift invariant **+1**.

G3b matches: orbit size **10**, lift **+1** on all 600 Nielsen tuples.
Geometric monodromy type of the pure-even multi-k envelope is the
**unique braid component** of Ni(A₅, 5A⁴).

---

## 1. Verified Nielsen representatives

| quantity | value |
|----------|------:|
| Raw Nielsen (product 1, generate A₅) | 600 |
| Verified lift invariant +1 | 600 |
| Conjugacy-normalised | 10 |
| Braid orbits | 1 |
| Orbit sizes | [10] |

### Locked conjugacy-normalised representative

- cycles: `['(0 1 2 3 4)', '(0 1 3 4 2)', '(0 1 4 2 3)', '(0 2 3 1 4)']`
- perms (images): `[[1, 2, 3, 4, 0], [1, 3, 0, 4, 2], [1, 4, 3, 0, 2], [2, 4, 3, 1, 0]]`
- verification: `{'labels': ['5A', '5A', '5A', '5A'], 'all_5A': True, 'product_id': True, 'generates_A5': True, 'lift_invariant': 1, 'lift_ok': True, 'cycles': ['(0 1 2 3 4)', '(0 1 3 4 2)', '(0 1 4 2 3)', '(0 2 3 1 4)'], 'perms': [[1, 2, 3, 4, 0], [1, 3, 0, 4, 2], [1, 4, 3, 0, 2], [2, 4, 3, 1, 0]], 'in_unique_orbit': True}`

### Sample tuples (verified lift +1)

- `['(0 3 1 2 4)', '(0 2 3 1 4)', '(0 2 4 3 1)', '(0 3 2 4 1)']`
- `['(0 3 1 2 4)', '(0 2 3 1 4)', '(0 4 1 3 2)', '(0 4 2 1 3)']`
- `['(0 3 1 2 4)', '(0 2 3 1 4)', '(0 2 1 4 3)', '(0 2 4 3 1)']`

### User-supplied cycle 4-tuples (verification)

- {'cycles': ['(0 1 4 2 3)', '(0 4 3 2 1)', '(0 1 2 3 4)', '(0 3 2 4 1)'], 'labels': ['5A', '5A', '5A', '5A'], 'all_5A': True, 'product_id': True, 'generates_A5': True, 'lift_invariant': 1, 'in_unique_orbit': True}
- {'cycles': ['(0 4 2 1 3)', '(0 4 3 2 1)', '(0 3 1 2 4)', '(0 2 4 3 1)'], 'labels': ['5A', '5A', '5A', '5A'], 'all_5A': True, 'product_id': True, 'generates_A5': True, 'lift_invariant': 1, 'in_unique_orbit': True}
- {'cycles': ['(0 1 4 2 3)', '(0 3 2 4 1)', '(0 1 2 3 4)', '(0 4 3 2 1)'], 'labels': ['5A', '5A', '5A', '5A'], 'all_5A': True, 'product_id': True, 'generates_A5': True, 'lift_invariant': 1, 'in_unique_orbit': True}

---

## 2. Explicit algebraic model over ℚ(t)

No external Magma/Sage Hurwitz package is available in this environment.
The **constructive model** is the pure-even multi-k path already shown
(G3/G3b) to realise monodromy type **5A⁴** — i.e. the unique braid component.

### Model A — `path_flag_classical` (flagship ↔ classical)

| item | value |
|------|-------|
| m | `5/16` (fixed) |
| k(t) | `4*(3*t - 2)/5` |
| α(t) | `-405*t**4 + 1080*t**3 - 1080*t**2 + 480*t - 55` |
| β(t) | `-972*t**5 + 3240*t**4 - 4320*t**3 + 2880*t**2 - 900*t + 88` |
| polynomial | `x**5 + (-405*t**4 + 1080*t**3 - 1080*t**2 + 480*t - 55)*x + (-972*t**5 + 3240*t**4 - 4320*t**3 + 2880*t**2 - 900*t + 88)` |
| disc identical square | **True** |
| square-free branch degree | 4 |
| branch sqf polynomial | `81*t**4 - 216*t**3 + 216*t**2 - 96*t + 11` |

$$f_t(x) = x^5 + \alpha(t)\, x + \beta(t) \in \mathbb{Q}(t)[x].$$

### Model B — `path_flag_lsw` (flagship ↔ LSW)

- m(t) = `25*t/8 + 5/16`
- k(t) = `-12*t/5 - 8/5`
- α(t) = `-405*t**4 - 1080*t**3 + 1420*t**2 + 20*t - 55`
- β(t) = `972*t**5 + 3240*t**4 - 1680*t**3 - 2320*t**2 + 100*t + 88`
- disc□ identity: **True**

### Model C — two-parameter envelope

- α(m,k) = `(-3125*k**4 + 65536*m**2)/256`
- β(m,k) = `k*(-3125*k**4 + 65536*m**2)/256`
- disc = `(256*alpha**2*m)**2`

1-parameter multi-k paths are rational curves on this surface; each such
path inherits geometric type **5A⁴**.

### Status of “no ready-made equation in open literature”

A classical Hurwitz package equation for Ni(A₅, 5A⁴) was not found as a
pre-packaged polynomial in-repo. **The pure-even path supplies an explicit
ℚ(t)-model of a 1-parameter family with that monodromy type**, which is the
object needed for geometric multi-k and Hilbert specialisation.

---

## 3. Monodromy re-verification on Model A

| item | value |
|------|-------|
| local labels | ['5A', '5A', '5A', '5A'] |
| nontrivial multiset | {'5A': 4} |
| **is 5A⁴** | **True** |

---

## 4. Hilbert specialisation → pure-even catalogue

### path_flag_classical

| quantity | value |
|----------|------:|
| Z specialisations tested | 15 |
| irreducible | 14 |
| even disc | 15 |
| catalogue hits | 5 |
| catalogue k | ['-4/5', '-8/5', '4/5'] |
| **multi-k** | **True** |

- {'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}
- {'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}
- {'tag': 'classical_m', 'k': '-4/5', 't': '1/3', 'alpha': 20, 'beta': -16}
- {'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}
- {'tag': 'classical', 'k': '4/5', 't': '1', 'alpha': 20, 'beta': 16}

### path_flag_lsw

| multi-k | **True** |
| catalogue k | ['-4', '-8/5'] |
| hits | [{'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'lsw_m100', 'k': '-4', 't': '1', 'alpha': -100, 'beta': 400}, {'tag': 'flagship', 'k': '-8/5', 't': '0', 'alpha': -55, 'beta': 88}, {'tag': 'lsw_m100', 'k': '-4', 't': '1', 'alpha': -100, 'beta': 400}] |

---

## 5. Reduced Hurwitz genus

| item | value |
|------|-------|
| r | 4 |
| type | 5A^4 |
| braid orbits | 1 |
| conj-norm orbit size | 10 |
| lift invariant | 1 |
| reduced dimension | r-3 = 1 (curve) |
| genus | None (see status) |

Exact reduced genus requires cusp ramification of H^rd→P1 (Bailey–Fried / Modular Tower / James tables). Orbit size 10 + lift +1 uniquely identifies the component; genus lookup is secondary to the explicit pure-even model.

If deg(H^rd→P1)=10 and indices (ind0,ind1,indinf) known, 2g-2 = 10*(-2) + ind0+ind1+indinf.

---

## 6. Fusion status

| test | result |
|------|--------|
| Unique braid component of Ni(A₅, 5A⁴) | **Yes** (1 orbit, size 10, lift +1) |
| Explicit f ∈ ℚ(t)[x] with monodromy 5A⁴ | **Yes** (pure-even multi-k path) |
| Disc□ identity (even monodromy) | **Yes** |
| Hilbert multi-k catalogue hits | **True** |
| **Arithmetic multi-k = geometric 5A⁴ type** | **True** |

**Fusion (type-level) closed:** the pure-even multi-k envelope paths are explicit ℚ(t)-models of the unique lift-+1 braid component of Ni(A₅, 5A⁴), and they Hilbert-specialise to multiple catalogue k-slices.

Remaining optional refinements: (i) birational identification of the pure-even
path with a classical Hurwitz chart of H^rd; (ii) exact reduced genus from
cusp tables; (iii) common-basepoint braid word for geometric monodromy.

---

## 7. Non-claims / stance

- Canonical T3 remains production dynamical baseline.
- Pure-even arithmetic theorems unchanged.
- Necessity remains open/paused.
- This is geometric multi-k work: type lock + explicit path model + catalogue Hilbert.

_Generated by `g3c_explicit_5a4.py`._
