# Closed-form flagship Mestre family P_t

**Seed:** P = x⁵ − 55x + 88  
**Mestre R:** R = x⁴ + 8x³ − 32x² + 33

```
P_t(z) = Res_y( P(y), z − y − t R(y) )
```

## Explicit monic form

```
P_t(z) = z⁵
  − 385 t z⁴
  − 440 t(380 t + 3) z³
  + 3520 t(18150 t² + 205 t + 3) z²
  + 55(45619200 t⁴ − 4364800 t³ − 21120 t² − 256 t − 1) z
  + 11(2269696000 t⁵ − 444928000 t⁴ + 21120000 t³ + 70400 t² + 165 t + 8)
```

## Checks

| Check | Result |
|-------|--------|
| t = 0 recovers seed | True |
| disc identically □ in ℚ(t) | True (square of a deg-10 poly) |
| Sample Gal | 10/10 A₅ for t ∈ {0, ±1, 2, 3, 5, 7, 9, 61, 80} |

Generative extension of the finished centre. Not a necessity claim.
