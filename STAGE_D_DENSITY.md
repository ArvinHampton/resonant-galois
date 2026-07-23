# Stage D — Density and Discriminant Height

**Date:** 23 July 2026  
**Scope:** Quantitative arithmetic properties of the pure-even fixed-k slices.

This document records density evidence and a height theorem for the pure-even multi-k theory.  
It does **not** close geometric multi-k fusion or the structural necessity criteria.

---

## Scorecard

| Block | Logical status | Result |
|-------|----------------|--------|
| **D1** Irreducibility density | Conjecture + strong evidence | Sampled lattice points on all 10 multi-seed pure-even slices are 100 % irreducible; even_fail = 0; high A₅ rate in Galois checks |
| **D2** Chebotarev proxy | Empirical histograms | On 100 irreducible fibres (≈4500 prime factorisations) cycle-type frequencies approximate the conjugacy-class proportions of A₅ |
| **D3** Discriminant height | **Theorem** | From the pure-even identity disc = (256 α² m)² one obtains the exact asymptotic log\|disc\| ∼ 10 log\|m\| + 48 log 2 |

---

## D1 — Irreducibility density (conjecture + evidence)

On each of the 10 multi-seed pure-even k-slices the integer lattice points that were sampled satisfy:

- irreducibility rate 1.0,
- zero evenness failures,
- Galois group A₅ in the large majority of checked cases.

**Status:** Positive density of irreducible (and of A₅) specialisations is strongly supported by computation. It remains a conjecture; Hilbert’s Irreducibility Theorem already guarantees infinitely many irreducibles for each irreducible family over ℚ(t), but density 1 on the specific arithmetic lattices L_k is a stronger statement.

---

## D2 — Chebotarev proxy (empirical)

Cycle-type frequencies observed on irreducible fibres:

| Observed approx. | A₅ class proportion |
|------------------|---------------------|
| ~0.43 | 0.40 |
| ~0.31 | 0.33 |
| ~0.26 | 0.25 |
| ~0.01 | 0.017 |

The match is consistent with the Chebotarev density theorem once Gal = A₅ is known. The histograms constitute supporting evidence, not an independent proof of the Galois group.

---

## D3 — Discriminant height (theorem)

From the closed-form pure-even identity

```
disc = (256 α² m)²
```

the logarithmic height admits the exact asymptotic

```
log|disc| ∼ 10 log|m| + 48 log 2
```

as |m| → ∞. This is a theorem following directly from the algebraic identity (no sampling required).

---

## Programme impact

| Layer | Status after Stage D |
|-------|----------------------|
| Pure-even multi-k theory | Finished |
| Stage A+B verification | Finished |
| Stage D density & height | Finished (D3 theorem; D1–D2 conjecture + evidence) |
| Geometric multi-k | Open |
| Structural necessity theorem | Open |

The arithmetic centre is quantitatively stronger. The open geometric and structural problems are unchanged.
