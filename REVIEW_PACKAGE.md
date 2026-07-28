# Review package — generative centre extensions

**Recommendation:** review **flagship Mestre lift first**, then **non-BJ B-avatar**.  
**Necessity:** remains **paused** either way. Both items enlarge the finished pure-even centre; neither is Crit-2 forcing.

---

## Review order

| Order | Object | Files | What a review checks |
|------:|--------|-------|----------------------|
| **1** | Flagship Mestre \(P_t\) | **`MESTRE_FLAGSHIP_PT.md`** | Explicit poly; disc □ identity in \(\mathbb{Q}[t]\); \(t=0\) recovery; sample Gal table |
| **2** | Non-BJ B-avatar | **`B_EMBED_LATTICE.md`** + **`EVENNESS_AVATAR.md`** | Disc identity in \(A\); \(T\)-embed \(d=-75\), \(bc=72A\); lattice \(A\to A_5\) samples |

---

## Why this order

| Criterion | Flagship \(P_t\) | Non-BJ B-avatar |
|-----------|------------------|-----------------|
| Novelty vs classical pure-even | **High** — deformation of a fixed seed into a 1-param family | **High** — beyond BJ (\(d\neq 0\)), fixed shape in \(A\) |
| Review cost | One closed form + finite Gal table | Identity + embed relations + lattice table |
| Error isolation | Single seed + single \(R\) + single resultant | Parallel track; independent of Mestre |
| Dependence | One HQCC seed + Mestre | Independent generative track |

The flagship lift is the **cleanest single theorem-shaped object**: one seed, one \(R\), one explicit \(P_t\), one disc identity, lattice \(t\) still \(A_5\). A focused review either certifies it or finds a coefficient/disc error quickly.

The B-avatar is the **better second review**: main beyond-BJ embed and second evenness avatar. Review it after the flagship so the two enlargements are checked without mixing error sources.

---

## Pass / fail criteria

### Primary — flagship \(P_t\)

| Claim | Pass if |
|-------|---------|
| Mestre condition | \(P''R-2P'R'\equiv 0\pmod{P}\) |
| Closed form | Matches \(\operatorname{Res}_y(P(y),z-y-t R(y))\), monic deg 5 |
| Seed recovery | \(P_0=P\) |
| Disc □ | \(\operatorname{disc}_z(P_t)\) square in \(\mathbb{Q}[t]\) |
| Samples | Listed \(t\) give irr + disc□ + \(A_5\) (or documented exceptions) |

### Secondary — B-avatar

| Claim | Pass if |
|-------|---------|
| Disc identity | \(\operatorname{disc}(P_A)=(18A(A^2+84375))^2\) |
| Embed | \(\chi_T=P_A\) under stated relations |
| Beyond BJ | \(d=-75\neq 0\) |
| Lattice samples | Stated \(A\in L_0\) specialisations match claimed Gal / disc□ |

---

## Out of scope for this review

- Necessity Criteria 1–3  
- Further cuts of unrestricted \(T\)  
- Geometric multi-\(k\) / genus of \(P(q,w)\)  
- Full lattice monoid / Mestre orbit graph (supporting docs only)

Supporting context (not first-pass review):

| Doc | Role |
|-----|------|
| `PURE_EVEN_MULTI_K.md` | Finished arithmetic centre |
| `L0_MESTRE_ORBIT.md` | Mestre on all lattice seeds |
| `L0_SECONDARY_INVARIANTS.md` | Secondary invariants |
| `TERNARY_LATTICE_DIRECTIONS.md` | Lattice research frame |

---

## Automated re-check

```bash
python review_flagship_b.py
python math_integrity_review.py
```

Expected: `flagship P_t: PASS` then `B-embed / avatars: PASS`; integrity suite **30/30 PASS**.  
Full write-up: **`MATH_INTEGRITY_REVIEW.md`** (independent re-derivation + verifier **VERDICT: PASS**).

---

## Stance (locked)

```
Finished centre:  pure-even multi-k + four-face organising principle
Review now:       (1) flagship Mestre Pt  (2) B-avatar
Necessity:        paused
```

## Contamination boundary

See **`CONTAMINATION_BOUNDARY.md`**.

| Risk | Locked-package status |
|------|------------------------|
| G₄ / 539.9 s as physical period in Gal proofs | Avoided — \(539\) lattice integer only |
| GW / Belle II “hits” as support for \(A_5\) | Out of scope |
| HQCC 539-step dynamics as disc-identity input | Not required |
| Motivation blurred into theorem hypothesis | Watch external write-ups |

**Rule:** Lattice integers may be motivated by the model; proofs may only use their arithmetic properties.
