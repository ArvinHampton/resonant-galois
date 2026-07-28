# Arithmetic over \(K=\mathbb{Q}(\sqrt{5})\)

_Partial run completed through relation tests; fingerprint phase interrupted (slow). Core theorems and scans are solid._

**Verdict:** Disc theorem **proved**. Over \(K=\mathbb{Q}(\sqrt{5})\), \(\operatorname{disc}(\mathrm{monic}(\varphi-t))\) is **identically a square** in \(K(t)\). Rational-\(t\) scan: **153/153** non-critical fibres even over \(K\), **152** irr over \(K\); **0** even over \(\mathbb{Q}\). Quadratic \(t=a+b\sqrt{5}\): **72/72** disc square in \(K\), all irr in sample. Disc-class vs pure-even seeds: **0** matches over \(\mathbb{Q}\), **20** over \(K\) (shared \(5\cdot\square\) class). **No descent** of evenness to \(\mathbb{Q}\). Pure-even \(m\in K\setminus\mathbb{Q}\) yields BJ coeffs in \(K\), not the HQCC \(\mathbb{Z}\)-lattice.

---

## 1. Theorem (proved)

Preferred cover: \(\varphi(y)=6y^5-15y^4+10y^3\).

Monic fibre:
\[
\mathrm{monic}(\varphi-t)=y^5-\tfrac52 y^4+\tfrac53 y^3-\tfrac t6.
\]

Discriminant identity in \(\mathbb{Q}(t)\):
\[
\operatorname{disc}_y = \frac{3125}{1296}\, t^2(t-1)^2
= 5\cdot\Bigl(\frac{25\,t(t-1)}{36}\Bigr)^2.
\]

Symbolic checks (script `k_sqrt5_even.py`):

| Identity | Status |
|----------|:------:|
| equals \(5\cdot(\mathrm{square})\) | **True** |
| equals \((\sqrt{5}\cdot\frac{25 t(t-1)}{36})^2\) | **True** |

**Over \(\mathbb{Q}\):** for \(t\in\mathbb{Q}\setminus\{0,1\}\), disc is \(5\) times a square, hence **never** a square in \(\mathbb{Q}\).

**Over \(K=\mathbb{Q}(\sqrt{5})\):** \(5=(\sqrt{5})^2\), so disc is a square in \(K\) for all non-critical \(t\in K\) (identically in \(K(t)\)).

---

## 2. Specialisations at rational \(t\)

Dense scan \(|p|\le 20\), \(q\le 6\):

| quantity | value |
|----------|------:|
| tested (non-critical) | 153 |
| even over \(K\) | **153** |
| even over \(\mathbb{Q}\) | **0** |
| irreducible over \(\mathbb{Q}\) | 152 |
| irreducible over \(K\) | 152 |
| even + irr over \(K\) | **152** |

Every non-critical rational fibre becomes even after base change to \(K\); none is even over \(\mathbb{Q}\).

---

## 3. Specialisations at \(t=a+b\sqrt{5}\)

| quantity | value |
|----------|------:|
| tested (\(b\neq 0\), \(\|a\|,\|b\|\le 4\)) | 72 |
| disc matches \(\sqrt{5}\)-square form | **72** |
| even + irr over \(K\) | **72** |

Quadratic base points behave as the theorem predicts.

---

## 4. Descent: can \(K\)-evenness reach \(\mathbb{Q}\)?

**No.**

- Z-model fibres at rational \(t\) lie in \(\mathbb{Q}[y]\).
- Their discriminants are of the form \(5\cdot(\mathrm{square})\) in \(\mathbb{Z}\) (permanent factor 5 in the square-free kernel for the geometric disc, surviving clearing of denominators up to squares).
- Base change to \(K\) makes disc a square; it does **not** produce a different \(\mathbb{Q}\)-model with disc a square.
- Same fibres remain odd over \(\mathbb{Q}\) (Gal typically in \(S_5\setminus A_5\)).

**Interpretation:** Route 2 removes the obstruction for **even monodromy of \(\varphi\)-fibres over \(K\)**, but does not create even arithmetic specialisations over \(\mathbb{Q}\).

---

## 5. Relation to pure-even BJ families

Catalogue pure-even families (LSW \(k=-4\), flagship \(k=-8/5\), classical \(k=4/5\), …) already have disc \(\square\) in \(\mathbb{Q}\) and are **not** rational fibres of \(\varphi\).

| Test | Result |
|------|--------|
| Disc-ratio square in \(\mathbb{Q}\) (fibre vs seed) | **0** fibres |
| Disc-ratio square in \(K\) | **20** (focused sample) — shared \(5\cdot\square\) vs \(\square\) becomes comparable over \(K\) |
| Number-field fingerprint (Frobenius) | interrupted; expected **no** \(\mathbb{Q}\)-isomorphism (different disc classes over \(\mathbb{Q}\)) |

### Pure-even families base-changed to \(K\)

- Parameter \(m=r+s\sqrt{5}\) (\(s\neq 0\)): disc identity still holds over \(K\).
- Typically \(\alpha,\beta\in K\setminus\mathbb{Q}\): seeds leave the HQCC \(\mathbb{Z}\)-lattice.
- Descent of such seeds to \(\mathbb{Q}\) forces constraints that collapse \(m\) back toward \(\mathbb{Q}\).

**Conclusion:** Over \(K\), \(\varphi\)-fibres and pure-even BJ families are both even, but they are not identified as the same \(\mathbb{Q}\)-objects. Relating them is at best a \(K\)-isomorphism / Tschirnhaus problem, not lattice recovery over \(\mathbb{Q}\).

---

## 6. Programme impact

| Claim | Status |
|-------|--------|
| Even monodromy of \(\varphi\)-fibres over \(K\) | **Achieved** (theorem + scan) |
| Even monodromy over \(\mathbb{Q}\) via same fibres | **Blocked** (permanent 5) |
| Recovery of HQCC \(\mathbb{Z}\)-seeds from \(\varphi\) over \(K\) | **Not obtained**; lattice interpretation breaks |
| Fusion Criterion-1 over \(\mathbb{Q}\) | **Still open** |

### Stance

1. Keep \(K=\mathbb{Q}(\sqrt{5})\) as a **side route** proving geometric evenness is possible after quadratic base change.
2. **Primary fusion fuel over \(\mathbb{Q}\)** remains multi-seed pure-even \(k\)-slices (`ENLARGED_SEED_CATALOGUE.md`): 60 A5 seeds, 10 multi-seed pure-even families.
3. Full fusion still needs either non-rigid geometry over \(\mathbb{Q}\) (Route 1) or a different construction (Route 3)—not further surgery on rigid \(\varphi/\mathbb{Q}\).

---

## 7. Code

- Script: `k_sqrt5_even.py`
- Related: `geometric_rigid_deform.py` (Q-obstruction), `enlarge_seed_catalogue.py` (pure-even slices)

_Generated from completed partial run of k_sqrt5_even.py; fingerprint section truncated._
