# Genus of P(q,w)=0 — locked (Tier 1.3)

**Status (2026-07-24).** Physical 3A⁴ eliminant is not a rational curve under this analysis. Global rational parameterisation not expected; single-valued f_s ∈ ℚ(s)[y] via param of P is blocked. Geometric multi-k remains open.

---

## Curve analysis

| Quantity | Value |
|----------|--------|
| Model | Physical 3A⁴ eliminant P(q,w)=0 |
| Degree | 6 (deg_q = deg_w = 3) |
| Arithmetic genus p_a | 10 |
| Irreducible over ℚ | Yes |
| Affine singularity | (1,1), multiplicity 3 |
| Lowest form | (Q+W)(Q²+QW+W²) — ordinary triple over ℂ |
| δ (if ordinary triple) | 3 → p_a−δ = 7 from affine alone |
| Infinity | 20q³w³; singular axis points [1:0:0], [0:1:0] |
| Genus estimate (affine + ∞ lower bound) | **g ≈ 1 (not 0)** |
| Projection from (1,1) residual degree | 3 (not birational param) |

**Genus 0?** No. Not a rational curve under this analysis. The (1,1)-line pencil does not parameterise P.

⇒ Single-valued f_s ∈ ℚ(s)[y] via rational param of P is **blocked**.

---

## Points and k sampling

| Source | Count |
|--------|-------|
| Verified ℚ-points (den ≤ 10) | 4 |
| Quadratic / special (e.g. ±1/√5) | 5 |
| Numeric real samples (subsampled) | 40 |
| Real BJ k samples | 120 |

Catalogue near-hits (tol 0.08): tags **classical, class_flip only** (k = ±4/5).

Strict geometric multi-k (≥2 multi-seed families, e.g. flagship + classical): **False**.

---

## Locked outcome

| Question | Answer |
|----------|--------|
| Genus 0? | **No** |
| Global rational param? | Not found |
| f_s ∈ ℚ(s)[y] via param of P? | **Blocked** |
| Point enum + BJ k | Done |
| Geometric multi-k | **Still open** |

---

## Geometric next options

1. Full blowup for exact geometric genus.
2. Work in the degree-3 function field over ℙ¹ (accept multi-valued coeffs in s).
3. Try other Nielsen types (2A 3A³, 2A² 3A²) for a genus-0 chart.
4. Treat arithmetic multi-k (envelope) as the completed multi-seed statement; leave Nielsen-labelled geometric multi-k as open research.

Script: `genus_p_qw.py`

See also: `EXPLICIT_3A4_EQUATION.md`, `FUSION_GAP.md`, `NECESSITY_THEOREM.md`.
