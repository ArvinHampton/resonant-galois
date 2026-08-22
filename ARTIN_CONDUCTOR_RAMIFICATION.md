# Artin conductors and ramification — executed census

**Status (2026-08-22).** Concrete disc factorisation, Frobenius cycle types, and ramification support for verified even A5 fibres. Full Artin conductor exponents (Swan) not computed; tame support and Chebotarev proxy recorded.

---

## 1. Disc factorisation census

All listed discs are perfect squares (square-free kernel = 1), consistent with image in A5.

### Flagship / PE-style

| Fibre | disc(f) | Factorisation | Support |
|-------|---------|---------------|--------|
| Flagship x^5−55x+88 | 58564000000 | 2^8 · 5^6 · 11^4 | {2,5,11} |
| LSW x^5−100x+400 | 77440000000000 | 2^16 · 5^10 · 11^2 | {2,5,11} |
| Classical x^5+20x+16 | 1024000000 | 2^16 · 5^6 | {2,5} |
| x^5+95x+76 | 2085136000000 | 2^10 · 5^6 · 19^4 | {2,5,19} |

### B-avatar P_A = x^5+75x^3+A x^2+3A

Identity verified: disc(P_A) = (18 A (A^2+84375))^2.

| A | Support (primes \| disc) |
|---|-------------------------|
| 3 | {2, 3, 293} |
| 61 | {2, 3, 61, 2753} |
| 80 | {2, 3, 5, 3631} |
| 243 | {2, 3, 83} |
| 539 | {2, 3, 7, 11, 23431} |
| −3, −61 | same as \|A\| |

Pattern: **persistent {2,3}**; moving primes from A and from A^2+84375.

### Mestre flagship P_t specialisations

| t | Support |
|---|--------|
| 0 (seed) | {2,5,11} |
| 1 | {2,5,11,17,15556511,34717435999} |
| 2 | {2,5,11,179,598966813,185712045863} |
| 3 | {2,5,11,6091,17333,9295763,595945517} |
| −1 | {2,5,11,31,367820476092638423} |

**Persistent Mestre/flagship ramification:** {2, 5, 11}.  
**Moving:** primes dividing Q(t) in disc = C · Q(t)^2.

---

## 2. Frobenius cycle types (unramified p < 100)

Degree sequence of f mod p ↔ cycle type in S5:

| Degrees | Cycle type | In A5? |
|---------|------------|--------|
| [5] | (5) | yes |
| [1,1,3] | (3,1,1) | yes |
| [1,2,2] | (2,2,1) | yes |
| [2,3] | (2,3) | no (odd) |
| [1,1,1,1,1] | id | yes |

### Flagship x^5−55x+88 (skip {2,5,11})

Summary among 22 primes: (5)×11, (3,1,1)×6, (2,2,1)×4, id×1.  
**No odd types.**

### B-avatar A=3 (skip {2,3,293})

(5)×10, (3,1,1)×9, (2,2,1)×4. No odd types.

### LSW / classical

Same even types only; frequencies in the same ballpark as A5 Chebotarev.

### A5 Chebotarev expected vs flagship sample

| Class | Expected | Flagship sample (p<100) |
|-------|----------|-------------------------|
| (5) | 0.40 | 0.50 |
| (3,1,1) | 0.33 | 0.27 |
| (2,2,1) | 0.25 | 0.18 |
| id | 0.017 | 0.045 |

Small-sample fluctuation; no contradiction.

---

## 3. Artin conductor — what is and is not computed

**Computed (support of ramification):**

- Primes dividing disc(f) = candidates for primes dividing the Artin conductor f(ρ) of the associated Artin representation (and for disc of the quintic field when monogenic).

**Not computed (full local exponents):**

- Swan conductors / wild depth at p=2,3,5.
- Exact f_p for the 4-dimensional irrep of A5 vs the permutation character.
- Proof that Z[α] = O_K (index may remove or add square factors relative to disc(K)).

**Tame heuristic:** for p > 5 unramified or tame, local contribution relates to 5 − (number of cycles of inertia). Full f(ρ) needs local resolvents or Magma/Sage NumberField conductor routines.

---

## 4. Cross-family comparison

| Family | Persistent primes | Moving primes controlled by |
|--------|-------------------|-----------------------------|
| PE / LSW / flagship seed | {2,5} or {2,5,11} | seed coefficients |
| Mestre P_t | {2,5,11} | Q(t) in disc identity |
| B-avatar | {2,3} | A and A^2+84375 |

**Monoid M0 control?**  
Lattice A ∈ {3,61,80,243,539,...} specialises B-avatar; ramification support then includes primes of A (often in M0 or small) **plus** primes of A^2+84375 (often large, outside M0).  
So M0 organises the **parameter**, not the full conductor support. Moving large primes are expected and not monoid-constrained.

---

## 5. Integrity notes

1. Square disc ⇒ even image (A5 or smaller even); Frobenius samples show types generating A5 (3-cycles + double transpositions + 5-cycles).
2. Residual D5 possible on other even fibres (known from secondary invariants) — disc□ is necessary not sufficient for A5.
3. Full Artin conductors remain a natural extension; this file locks **support and Chebotarev proxy** only.
4. Necessity unchanged: ramification fingerprints do not force image from HQCC axioms alone.

---

## 6. Optional next

1. Magma/Sage: `NumberField(f).conductor()` or Artin representation conductor for a few fibres.
2. Inertia type at p=2,3,5 via local factorisation / resolvents.
3. Tabulate radical(disc) vs height for Stage D refinement.

Script basis: inline CAS census 2026-08-22.
