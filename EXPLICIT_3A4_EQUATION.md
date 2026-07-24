# Explicit Ni(A₅, C₃⁴) model (3A⁴)

## Reduced Hurwitz curve

```
H^(rd) ≅ ℙ¹_s   over ℚ (genus 0, Bailey–Fried)
```

Cross-ratio chart: branch values {0, 1, ∞, s}. Infinitely many rational points: s ∈ ℚ ∖ {0,1}.

## Cover normal form

```
φ(y) = N(y)/D(y)
N = c y³ (y−1)(y−p₂)
D = y² − σ y + π
```

Triple roots of φ−1 and φ−s at parameters q, w.

## Explicit formulae

```
c = −1 / (q (6 p₂ q − 3 p₂ − 10 q² + 6 q))
σ = q (8 p₂ q − 3 p₂ − 15 q² + 8 q) / (6 p₂ q − 3 p₂ − 10 q² + 6 q)
π = q² (3 p₂ q − p₂ − 6 q² + 3 q) / (6 p₂ q − 3 p₂ − 10 q² + 6 q)
```

Eliminant physical component P(q,w) = 0 (degree-3 plane curve in each variable):

```
P(q,w) = 20 q³ w³ − 40 q³ w² + 27 q³ w − 6 q³
       − 40 q² w³ + 73 q² w² − 45 q² w + 9 q²
       + 27 q w³ − 45 q w² + 26 q w − 5 q
       − 6 w³ + 9 w² − 5 w + 1
```

## Degree-5 resolvent

```
f(y) = monic_y(N(y) − t D(y))
```

Single-valued f_s ∈ ℚ(s)[y] still **open** (H^(rd) → ℙ¹_s is multi-sheeted in this chart).

## Exact fibre s = −1 over ℚ(√5)

```
c = −√5,  p₂ = −1,  r_{1,2} = ±1/5
q, w = ±1/√5
y⁵ − y³ + (t/√5)(y² − 1/25) = 0
```

Norm to ℚ(t): degree 10 (not 5).

## Specialisation vs pure-even catalogue

| Quantity | Value |
|----------|-------|
| Catalogue seed hits (tested grid) | 0 |
| Distinct catalogue k | ∅ |
| Multi catalogue k | False |

## Status

| Item | Status |
|------|--------|
| H^(rd) ≅ ℙ¹ | Locked |
| Eliminant P(q,w)=0 | Explicit |
| c, σ, π, resolvent N−tD | Explicit |
| Fibre s=−1 | Explicit |
| f_s ∈ ℚ(s)[y] single-valued | Open |
| Geometric multi-k | **Open** |

Arithmetic multi-k (envelope / paths) remains the only completed multi-seed statement.
