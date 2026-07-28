# \(A_5\) Hurwitz spaces with \(r=4\) — positive-dimensional strata

_Elapsed: 57.67s_

**Verdict:** Filter-pass r=4 A5 types (≤ programme filter): 19. Lookup g=0 candidates: 15; g=1: 2. Nielsen enumeration + braid orbits computed for filter-pass types (see table). Explicit Nielsen-realised multi-k geometric family: NONE YET. Arithmetic multi-k (envelope deg-4 path flagship↔classical): yes. LSW stays on single k=-4. Next: realise a g=0 class (e.g. 3A^4 or 2A,3A,3A,5A) by explicit equation and re-test multi-k Hilbert hits.

---

## 1. Setup

For finite \(G\) and conjugacy classes \(C=(C_1,\dots,C_r)\), the Hurwitz space
\(\mathcal{H}(G,C)\) parametrises Galois covers of \(\mathbb{P}^1\) with monodromy \(G\)
and ramification type \(C\). Components \(\leftrightarrow\) braid orbits on Nielsen tuples.
Dimension \(r-3\) after \(\mathrm{PGL}_2\):

- \(r=3\): rigid (isolated covers) — \(\varphi\) abandoned for \(\mathbb{Q}\)-fusion
- \(r=4\): curves (this document)
- \(r\ge 5\): higher-dimensional strata

### \(A_5\) classes

| class | size | cycle type |
|-------|-----:|------------|
| 2A | 15 | 2+2+1 |
| 3A | 20 | 3+1+1 |
| 5A | 12 | 5 |
| 5B | 12 | 5 (inverse class) |

**Programme filter:** \(\ge 2\) factors \(3A\), or \(\ge 1\) of \(3A\) and \(\ge 1\) of \(5A/5B\).

Lift invariant: Fried–Serre; for \(A_5\) typically \(\pm 1\) and the only braid-orbit obstruction
(Magaard–Shpectorov–James). Lookup genera from standard IG tables (not re-proved here).

---

## 2. Nielsen classes and braid orbits (computed)

| type (sorted) | filter | raw gen. tuples | # braid orbits | orbit sizes | lookup g |
|---------------|:------:|----------------:|---------------:|-------------|----------|
| `2A,2A,3A,3A` | yes | 1080 | 1 | [108] | 0,1 |
| `2A,2A,3A,5A` | yes | 900 | 2 | [180, 180] | 0,1 |
| `2A,2A,3A,5B` | yes | 900 | 2 | [180, 180] | 0 |
| `2A,3A,3A,3A` | yes | 1440 | 1 | [96] | 0,1 |
| `2A,3A,3A,5A` | yes | 1200 | 2 | [240, 240] | 0,1 |
| `2A,3A,3A,5B` | yes | 1200 | 2 | [240, 240] | 0 |
| `2A,3A,5A,5A` | yes | 730 | 3 | [288, 144, 144] | 1,0 |
| `2A,3A,5A,5B` | yes | 710 | 3 | [288, 144, 144] | 0 |
| `2A,3A,5B,5B` | yes | 730 | 3 | [288, 144, 144] | — |
| `3A,3A,3A,3A` | yes | 1080 | 1 | [18] | **0 (LOCKED)** |
| `3A,3A,3A,5A` | yes | 1500 | 4 | [40, 60, 60, 40] | 0,1 |
| `3A,3A,3A,5B` | yes | 1500 | 4 | [40, 60, 60, 40] | 0 |
| `3A,3A,5A,5A` | yes | 920 | 6 | [144, 90, 60, 90, 12, 12] | 0,1 |
| `3A,3A,5A,5B` | yes | 1120 | 6 | [144, 90, 60, 90, 12, 12] | 0 |
| `3A,3A,5B,5B` | yes | 920 | 6 | [144, 90, 12, 60, 12, 90] | 0 |
| `3A,5A,5A,5A` | yes | 696 | 6 | [36, 36, 36, 36, 72, 72] | 1,0 |
| `3A,5A,5A,5B` | yes | 488 | 6 | [36, 36, 36, 72, 72, 36] | 1 |
| `3A,5A,5B,5B` | yes | 488 | 6 | [36, 36, 36, 72, 72, 36] | 1 |
| `3A,5B,5B,5B` | yes | 696 | 6 | [36, 36, 36, 36, 72, 72] | — |

### Excluded by filter (no table rows)

Types with no double-3A and without 3A+5* (e.g. pure \(2A^4\), pure \(5^*\) without \(3A\))
are omitted from the braid computation above (still listed in JSON).

---

## 3. Candidates with lookup genus \(0\) or \(1\)

| type | orbits | sizes | best g (lookup) | notes |
|------|-------:|-------|:---------------:|-------|
| `2A,2A,3A,3A` | 1 | [108] | 0 | often g=0; strong candidate; second component sometimes g=1 |
| `2A,3A,3A,3A` | 1 | [96] | 0 | excellent ternary candidate; possible second component |
| `2A,3A,3A,5A` | 2 | [240, 240] | 0 | 3A+3A+5+2; top candidate; companion |
| `3A,3A,3A,3A` | 1 | [18] | 0 | pure ternary; classic A5 family candidate; second orbit |
| `3A,3A,3A,5A` | 4 | [40, 60, 60, 40] | 0 | strong candidate; companion |
| `3A,3A,5A,5B` | 6 | [144, 90, 60, 90, 12, 12] | 0 | candidate; both 5-classes |
| `2A,2A,3A,5A` | 2 | [180, 180] | 0 | candidate; ternary+double transp; companion orbit |
| `2A,2A,3A,5B` | 2 | [180, 180] | 0 | as 5A by outer aut of A5 swapping 5A/5B |
| `2A,3A,3A,5B` | 2 | [240, 240] | 0 | as 5A via outer |
| `2A,3A,5A,5A` | 3 | [288, 144, 144] | 0 | may be g≥1; check orbit |
| `2A,3A,5A,5B` | 3 | [288, 144, 144] | 0 | candidate |
| `3A,3A,3A,5B` | 4 | [40, 60, 60, 40] | 0 | as 5A |
| `3A,3A,5A,5A` | 6 | [144, 90, 60, 90, 12, 12] | 0 | candidate; companion |
| `3A,3A,5B,5B` | 6 | [144, 90, 12, 60, 12, 90] | 0 | as 5A,5A |
| `3A,5A,5A,5A` | 6 | [36, 36, 36, 36, 72, 72] | 0 | often higher genus; if orbit small |
| `3A,5A,5A,5B` | 6 | [36, 36, 36, 72, 72, 36] | 1 | check |
| `3A,5A,5B,5B` | 6 | [36, 36, 36, 72, 72, 36] | 1 | check |
| `2A,3A,5B,5B` | 3 | [288, 144, 144] | None |  |
| `3A,5B,5B,5B` | 6 | [36, 36, 36, 36, 72, 72] | None |  |

### Priority shortlist (ternary-friendly, g=0 lookup)

1. **\(3A^4\)** — pure ternary; classical A5 family candidate
2. **\(2A,3A^3\)** — ternary + double transposition
3. **\(2A^2,3A^2\)** — often g=0
4. **\(3A^3,5A\)** / **\(2A,3A^2,5A\)** — 3-cycles + 5-cycle
5. **\(3A^2,5A,5B\)** — both 5-classes

---

## 4. Explicit models vs multi-\(k\) arithmetic lattice

### LSW (fixed \(k=-4\))

- specs tested (irr disc□): 140
- distinct \(k\): ['-4']
- multi-\(k\): **False**
- catalogue multi-\(k\): **False**
- Arithmetic pure-even; specialisations stay on k=-4. Not a multi-k geometric family. Branch locus in t is thin (disc zeros).

### Bring / BJ probes

- bring \(x^5+x+t\): n=0, k's=[], multi=False
- α=β family: n=0, k's=[], multi=False

### Envelope path flagship \(\leftrightarrow\) classical (arithmetic multi-\(k\))

- multi catalogue \(k\): **True**
- hits: `[{'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'classical', 'k': '4/5', 't': '1'}]`
- Arithmetic envelope path (not a Hurwitz Nielsen realisation). Hits k=-8/5 and k=4/5 by construction when endpoints included.

### Same-\(m\) linear-\(k\) model (degree 4 in \(u\))

- disc□: **True**
- α numerator degree: 4
- hits \(k\): ['-8/5', '4/5']
- Nielsen-realised? **False**
- Degree-4 rational coefficient path joining flagship to classical. Pure-even multi-k arithmetic. Not identified with an explicit Nielsen class equation in this computation.

---

## 5. Conclusions

1. **Positive-dimensional \(A_5\) strata exist** for many filter-pass \(r=4\) types;
   braid orbits were enumerated from Nielsen tuples for those types.

2. **Genus-0 shortlist** (lookup) is non-empty: especially \(3A^4\), \(2A3A^3\),
   \(2A^2 3A^2\), \(3A^3 5A\), \(2A3A^2 5A\). These are the natural geometric targets.

3. **No explicit Nielsen-class equation** in this run was shown to Hilbert-specialise
   onto **two or more** fixed-\(k\) pure-even catalogue families. LSW stays on \(k=-4\);
   Bring-like probes do not give multi-\(k\) lattice hits.

4. **Arithmetic multi-\(k\)** remains available (envelope / same-\(m\) deg-4 path
   flagship\(\leftrightarrow\)classical) but is **not** yet identified with a Hurwitz
   Nielsen realisation — the geometric multi-\(k\) goal is still open.

5. **Sharpness:** if every g=0 four-point family that can be written over \(\mathbb{Q}\)
   only meets the BJ pure-even lattice in a single ratio class \(k\), the obstruction
   is geometric rather than Diophantine. Testing that requires explicit equations
   for the shortlist (Malle–König style resolvents / computational IG).

### Recommended next computations

1. Realise **one** g=0 class over \(\mathbb{Q}\) explicitly — priority \(3A^4\) or \(2A,3A,3A,5A\).
2. Compute the degree-5 resolvent (or BJ reduction) as a family in the curve parameter.
3. Test specialisations against the 10 multi-seed pure-even \(k\)-slices
   (`ENLARGED_SEED_CATALOGUE.md`).
4. If multi-\(k\) hits appear, the geometric multi-\(k\) goal is achieved; if not,
   try the next shortlist type or \(r=5\) strata.

_Generated by a5_hurwitz_r4.py_