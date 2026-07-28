# Pure-even multi-\(k\) arithmetic — theorem-grade centre

**Resonant Galois Programme · citable core**  
**Status:** Finished arithmetic theory. HQCC enters only as lattice motivation and specialisation source.  
**Not claimed:** Forced alternating monodromy from HQCC axioms (necessity; open / paused).

Verification: `lib/lemmas.py`, `pure_even_specialisations.py`, `stage_d_density.py`, `t_subclass_identical_square.py`.

---

## Abstract

We record the finished **pure-even multi-\(k\)** theory for monic Bring–Jerrard quintics
\[
f=x^5+\alpha x+\beta\in\mathbb{Q}[x],\qquad k=\beta/\alpha\in\mathbb{Q}\setminus\{0\}.
\]
On the **pure-even locus**
\[
\alpha=256m^2-\frac{3125\,k^4}{256},\qquad\beta=k\alpha
\quad(m\in\mathbb{Q}),
\]
the discriminant is an **identical square** in the parameter ring:
\[
\operatorname{disc}(f)=\bigl(256\,\alpha^2 m\bigr)^2.
\]
Hence every irreducible specialisation has even monodromy (\(\mathrm{Gal}\le A_5\) among transitive subgroups of \(S_5\)). Combined with a Frobenius of type \((3,1,1)\), one obtains \(\mathrm{Gal}=A_5\) by the standard group-theoretic criterion.

Fixed-\(k\) rays, the two-parameter envelope \((m,k)\), homogenisation, and rational cross-\(k\) paths are theorem-grade. A lattice of HQCC / resonant integers supplies many \(Z\)-coefficient fibres (catalogue + Stage D density tables). That lattice is **motivation and specialisation data**, not a proof that every HQCC-built object has alternating monodromy.

---

## 0. Role of HQCC (what is and is not claimed)

| Claim | Status |
|-------|--------|
| BJ disc identity; pure-even \(k\)-slices; envelope; paths | **Theorem-grade** |
| Operational \(A_5\) criterion (irr + disc □ + type \((3,1,1)\)) | Classical group theory |
| Resonant / HQCC lattice as a rich source of specialisations | **Generative / empirical** |
| Alternating monodromy **forced** by HQCC axioms alone | **Not claimed** (necessity open; paused as research) |
| Four-face ternary organising principle | **Structural reading** of generative success (`TERNARY_ORGANIZING_PRINCIPLE.md`) |

**Citable result of the programme:** pure-even multi-\(k\) + organising principle + verification tables.  
**Open research:** geometric multi-\(k\) fusion (**principal next** — `GEOMETRIC_MULTI_K_FUSION.md`); necessity (Criteria 1–3, paused); arboreal dynamics.

**Contamination boundary:** Lattice integers \(\{3,9,27,61,80,243,539,\ldots\}\) may be *motivated* by the model; every identity in this document uses only their status as elements of \(\mathbb{Z}\) or \(\mathbb{Q}\). G₄/539.9 s, GW/Belle II, and 539-step dynamics are **not** hypotheses here (`CONTAMINATION_BOUNDARY.md`).

---

## 1. Classical setup

### 1.1 Discriminant of Bring–Jerrard quintics

**Theorem 1 (BJ discriminant).** For \(a,b\) in a field of characteristic not \(2\) or \(5\),
\[
\operatorname{disc}\bigl(x^5+ax+b\bigr)=256a^5+3125b^4.
\]

*Proof sketch.* Direct expansion of the Sylvester / symbolic discriminant; verified in `lib/lemmas.py` (`verify_disc_formulas`, symbolic identity).

### 1.2 Even monodromy and \(A_5\)

**Theorem 2 (operational \(A_5\)).** Let \(f\in\mathbb{Z}[x]\) be monic irreducible of degree \(5\). If \(\operatorname{disc}(f)\) is a square in \(\mathbb{Z}\) and some unramified prime has Frobenius cycle type \((3,1,1)\), then
\[
\mathrm{Gal}(f/\mathbb{Q})=A_5.
\]

*Reason.* The only transitive subgroup of \(S_5\) contained in \(A_5\) and containing a \(3\)-cycle is \(A_5\) itself.

### 1.3 Homogenisation

**Theorem 3 (homogenisation).** For any \(\alpha,\beta\) with \(\operatorname{disc}(x^5+\alpha x+\beta)\) a square, the family
\[
f_t=x^5+\alpha t^4 x+\beta t^5
\]
satisfies
\[
\operatorname{disc}(f_t)=t^{20}\operatorname{disc}(\mathrm{seed})
\]
for \(t\neq 0\). Hence disc is a square for all such \(t\), and irreducible specialisations have even monodromy.

---

## 2. Pure-even fixed-\(k\) rays

**Theorem 4 (pure-even \(k\)-slice).** Fix \(k\in\mathbb{Q}\setminus\{0\}\). Set
\[
\alpha(m)=256m^2-\frac{3125\,k^4}{256},\qquad
\beta(m)=k\cdot\alpha(m).
\]
Then, identically in \(\mathbb{Q}(m)\),
\[
\operatorname{disc}\bigl(x^5+\alpha(m)x+\beta(m)\bigr)
=\bigl(256\,\alpha(m)^2 m\bigr)^2.
\]

*Proof.* Substitute \(\beta=k\alpha\) into Theorem 1:
\[
256\alpha^5+3125k^4\alpha^4=\alpha^4\bigl(256\alpha+3125k^4\bigr).
\]
The pure-even choice \(256\alpha+3125k^4=256\cdot(256m^2)=65536\,m^2\) yields
\[
\alpha^4\cdot 65536\,m^2=(256\alpha^2 m)^2.
\]

**Corollary.** On every pure-even \(k\)-slice, evenness is an **identity**, not a search outcome. Irreducibility and the presence of type \((3,1,1)\) remain number-theoretic (Hilbert / Chebotarev).

### Flagship and classical rays (examples)

| \(k\) | name | sample \((\alpha,\beta)\) | note |
|------:|------|---------------------------|------|
| \(-8/5\) | flagship | \((-55,88)\) at suitable \(m\) | \(88=61+3^3\) lattice |
| \(4/5\) | classical | \((20,16)\), \((95,76)\), … | classical BJ ray |
| \(-4\) | LSW | \(\alpha=t^2-3125\), \(\beta=-4\alpha\) | cleared form |

---

## 3. Two-parameter envelope and cross-\(k\) paths

**Theorem 5 (envelope).** Over \(\mathbb{Q}(m,k)\) with \(k\neq 0\), the same formulae give
\[
\operatorname{disc}=\bigl(256\alpha^2 m\bigr)^2
\]
identically. This is a **two-parameter pure-even family**.

**Theorem 6 (cross-\(k\) paths).** Let \((m_0,k_0)\) and \((m_1,k_1)\) be rational points of the envelope with nonzero \(k_i\). Any rational path
\[
\bigl(m(u),k(u)\bigr)\in\mathbb{Q}(u)^2,\quad k(u)\neq 0,
\]
joining them yields a one-parameter pure-even family with disc identically square in \(\mathbb{Q}(u)\).

*Consequence.* Distinct multi-seed catalogue ratios (flagship, classical, LSW, …) lie on a single connected pure-even surface; joining them does not break evenness.

---

## 4. Lattice specialisations (HQCC as data source)

An **HQCC seed** (definitional, not a necessity claim) is a BJ polynomial \(x^5+\alpha x+\beta\in\mathbb{Z}[x]\) with:

1. \(\alpha,\beta\) from the resonant lattice \(\{3,9,27,61,80,243,539,\ldots\}\) or short combinations;
2. disc square and \(\mathrm{Gal}=A_5\) (or even + operational criteria);
3. homogenisation family as in Theorem 3.

**Examples.** \(x^5-55x+88\); \((95,\pm76)\); \((95,\pm532)\); LSW cleared forms.

**Empirical (reproducible).**

| Experiment | Result | Doc |
|------------|--------|-----|
| Pure-even \(Z\)-coeff specialisations | 1728 pts, irr=1728, even_fail=0, many \(A_5\) | `PURE_EVEN_SPECIALISATIONS.md` |
| Stage D1 irr density | irr rate \(1.0\) on large multi-seed samples | `STAGE_D_DENSITY.md` |
| Stage D2 Chebotarev proxy | cycle-type histograms near \(A_5\) class densities | same |
| Stage D3 disc height | \(\log\|\mathrm{disc}\|=10\log\|m\|+48\log 2+o(1)\) | **proved** from identity |

Model integers enter as **Hilbert specialisations of proved-even families**, not as free coefficient search on which disc happens to be square.

---

## 5. Matrix templates (context only)

The structural template
\[
T(a,b,c,d,e,f)=\begin{pmatrix}
0&1&0&0&0\\0&0&1&0&0\\a&0&0&b&e\\0&0&0&0&1\\c&f&0&d&0
\end{pmatrix}
\]
has
\[
\chi_T=x^5-d x^3-(a+ef)x^2-(bf+ce)x+(ad-bc).
\]
On the **BJ-embed** \(d=0\), \(a=-ef\), one recovers \(\chi=x^5+\alpha x+\beta\) and the pure-even theory. Base matrix \(M=T(3,80,61,-3,0,0)\) has **odd** disc (\(S_5\)). Templates alone do **not** force disc □; see `HQCC_MATRIX_TEMPLATES.md`, `TIER11_DEEPEN.md`.

---

## 6. What this paper deliberately omits

- Necessity Criteria 1–3 (open; programme stance: **paused** as citable claim).
- Geometric multi-\(k\) via Nielsen / \(3A^4\) resolvents (open research; \(P(q,w)\) has \(g>0\)).
- Physical / entanglement narratives (separate from Galois arithmetic).
- Claims that pure-even “needs” 539 or the T₃ map for the disc identity — it does **not**.

---

## 7. How to verify

```bash
python -c "from lib.lemmas import verify_disc_formulas; print(verify_disc_formulas())"
python pure_even_specialisations.py
python stage_d_density.py
```

Machine-readable tables: `build/PURE_EVEN_SPECIALISATIONS.json`, `build/STAGE_D_DATA.json`.

---

## 8. Document map

| Doc | Role |
|-----|------|
| **This file** | Standalone theorem-grade core |
| `THEOREMS.md` | Full programme ledger (thms 1–12) |
| `TERNARY_ORGANIZING_PRINCIPLE.md` | Four-face structural reading |
| `NECESSITY_THEOREM.md` | Open / paused necessity target |
| `RESEARCH_ROADMAP.md` | Priority after publish |

---

## References (internal)

1. BJ disc and pure-even identity — classical; machine check `lib/lemmas.py`.  
2. Homogenisation and multi-seed catalogue — package docs `PURE_EVEN_*`, `ENLARGED_SEED_CATALOGUE.md`.  
3. Stage D — `STAGE_D_DENSITY.md`.  
4. Negative control (rigid fibre odd) — `RIGID_FIBRE_T3.md`.

_Generated as the publishable centre of Resonant Algebra; HQCC lattice is motivation, not necessity._
