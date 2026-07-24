# Non-classical resonant-field track (scaffold)

Stays in algebraic number theory but leaves f ∈ ℤ[x], Gal(f/ℚ).  
**Classical pure-even / ℤ centre is unchanged.**

## Working field

```
ℛ = ℚ(2 cos 2π/539),   [ℛ : ℚ] = 210
```

Proxies used computationally: R_n = ℚ(2 cos 2π/n), n ∈ {5,7,11,15}.

## N1 — Coefficients in R_n

BJ shapes x⁵ + α x + β with α, β ∈ ℤ[ξ], ξ = 2 cos(2π/n):

- Classical descent (e.g. −55, 88) recovers Gal/ℚ = A₅.
- Resonant deformations (3ξ−55, 88−ξ, …) need relative Gal / Weil restriction — open computation.

## N2 — Split-prime Frobenius labelling

Primes that split completely in R_n give an extra conjugacy label (cycle type, p split in R_n).  
Operational on proxies; refines Chebotarev along the resonant place.

## N3 — Pure-even over ℛ (theorem-grade, field-agnostic)

Identity (char ≠ 2, 5):

```
α = 256 m² − 3125 k⁴ / 256,   β = k α,   disc = (256 α² m)²
```

in any field F(m,k), including m, k ∈ ℛ.

## Locked distinctions

| Symbol | Meaning |
|--------|---------|
| s | Cross-ratio of four branch points (Hurwitz / M_{0,4}) |
| k | Pure-even ratio β/α of a BJ fibre |

In general **k ≠ s**.

Gal(R_n/ℚ) does **not** make multi-k paths into unions of cosine orbits (orbits finite; paths infinite; rational k are fixed points).

## Existence of f ∈ R_539[x]

- Trivial reading: YES via HQCC seeds over ℤ ⊂ R_539.
- Strong non-rational lift reducing to one fixed seed at infinitely many split primes: NO (congruence rigidity).
- Intermediate geometric readings: OPEN.

## Stance

| | |
|--|--|
| Centre | Pure-even multi-k over ℚ, ℤ-seeds |
| This track | Enrichment / research scaffold |
| Not | Replacement of classical arithmetic multi-k |
