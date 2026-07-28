# RESOLUTION BUILD — closing the gap

## Programme position (locked)

> The **arithmetic** side of the resolution criteria now rests on **proved identities**
> for an infinite family whose coefficients are visibly **HQCC-lattice**.
> The **geometric cover** remains the **principal open problem**.

Ledger: **`THEOREMS.md`**. HQCC seed definition: **`HQCC_SEED.md`**.

## What a full resolution still requires (geometry)

1. **Canonical cover** (or moduli family) from HQCC/resonant *branch* data
   with **proved** alternating geometric monodromy — not only arithmetic specialisations.
2. (Optional strengthening) axioms on full structural \(T_n\) forcing disc² outside BJ embed.
3. (Optional) sign character for unrestricted monodromy representations.

Experiment (catalogues of \(A_5/A_6\)) remains necessary evidence and regression.
Arithmetic theorems: **`THEOREMS.md`**, `THEOREM_ATTACK.md`, `HQCC_STRICT.md`.

## Build status

| Step | Role |
|------|------|
| catalogue | Frozen A5/A6/D5 lists |
| criterion1–3 | Exploratory scaffolds + rates |
| **theorem_attack** | Thin-class lemmas + 1-param families + sign theorems |

Rebuild: `python build_all.py` or `python theorem_attack.py`.

## Experimental catalogues (engine room)

- Unique **A5:** 36
- Unique **A6:** 4
- D5 near-misses (sample set): 7

Details: `build/CATALOGUE.md`, `a5_brute_results/A6_T6.md`, `DEFORM_M.md`.

## Criterion 1 — Canonical HQCC object

Scaffold: Möbius/HQCC blocks, cubic resultants, BJ/near-rigid families.
See `build/CRITERION1_HQCC.md`, **`THEOREM_ATTACK.md`**, **`HQCC_NATIVE.md`**, **`HQCC_STRICT.md`**.

**Classical theorem class:** \(f_t=x^5+20t^4 x+16t^5\) — proved even; all tested
specialisations Gal \(=A_5\) (seed not HQCC-primary).

**HQCC-native theorem class (new):** strict HQCC lattice search found
**12 non-classical** BJ seeds with Gal \(A_5\), e.g.

- \(x^5 - 55x \pm 88\) (near punctures/ternary: \(88=61+27\))
- \(x^5 + 95x \pm 76\), \(x^5 + 95x \pm 532\) (\(532=539-7\), period-adjacent)
- \(x^5 - 100x \pm 400\), \(x^5 + 124x \pm 496\)

Each admits the **same homogenised theorem**:
\(\operatorname{disc}(x^5+\alpha t^4 x+\beta t^5)=t^{20}\operatorname{disc}(\mathrm{seed})\) square
for all \(t\neq 0\). Homogenised families: **10/10 A5** on tested \(t\) (see `HQCC_STRICT.md`).

**Arithmetic Crit 1:** settled at theorem grade for HQCC-lattice homogenised families.

**Principal open (geometry):** monodromy of a cover built from HQCC *branch cycle types*
(resultant cubics still prefer \(A_3\); BJ families are arithmetic, not geometric covers).

## Criterion 2 — Axioms ⇒ disc² + 3-cycles

Evenness obstruction documented (base \(M\), base \(T_6\)).
Subclass rates: `build/CRITERION2_AXIOMS.md`. Thin-class theorems:
`THEOREM_ATTACK.md`.

**Partial advance (lemma):** on the BJ class \(x^5+ax+b\),
disc \(=256a^5+3125b^4\); evenness ⇔ that integer is a square.
Homogenised A5 seed is a 1-param **theorem class**.
Self-adjoint / det±1 / ω-norm thin classes still have disc² rate \(\sim 1\%\) or all-reducible.

**Gap remaining:** no axiom list proved to force disc² for all structural \(T_n\) matrices.

## Criterion 3 — Sign character

Correlations: det(M), ternary weight vs disc² — `build/CRITERION3_SIGN.md`.

**Partial advance:** Crit 3 **solved on BJ thin class** and on the homogenised
A5 family (closed-form / proved square disc). Full T5 lattice still open
(best empirical invariants rate ≪ 1).

**Gap remaining:** no ternary invariant implying trivial sign for *all* HQCC monodromy.

## Pipeline (unchanged)

```
template → χ → irr / disc² / cycles → Gal ID on survivors
```

## Next mathematical moves (priority)

1. ~~HQCC-native analogue of homogenised A5~~ — **done** (lattice BJ seeds + lemma 5–6; see `THEOREMS.md`).
2. ~~T5 disc structure~~ — **partial**: BJ embed subclass of \(T_5\) contains all proved-even families; full \(D\) factored (~99 terms). See `T5_DISC_IDEAL.md`.
3. Geometric HQCC branch-cycle cover (true Crit 1 geometry) — still open.
4. Non-BJ thin class with \(D\equiv S^2\) identically — still open.
5. Keep catalogues as **regression tests**.

## Paths

| Path | Content |
|------|---------|
| `resonant_galois/THEOREM_ATTACK.md` | Theorem-promotion results |
| `resonant_galois/IMPLICATIONS.md` | Claims / non-claims map |
| `resonant_galois/build/` | Regenerable criterion reports |
| `a5_brute_results/` | Prior scan archives |
| `resonant_galois/build_all.py` | Full rebuild entrypoint |
| `resonant_galois/theorem_attack.py` | Criteria 1–3 theorem pass |

Full implications map: **`IMPLICATIONS.md`**.
