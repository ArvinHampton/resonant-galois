# Tier 1.2 sharp next options — executed

**Status (2026-07-24)**

| Option | Outcome |
|--------|--------|
| A. Binary chooses \(k\), then pure-even under \(\mathcal{H}\) | Works as **composite**; not necessity |
| B. \(F\to T(\ldots)\) only, disc\(\square\)\(\to 1\) | **Failed** (rate 0.00) |
| C. New identical-square subclass of \(T\) | **None** beyond known families |

Necessity fragment: **still open**.  
Finished centre: pure-even multi-\(k\) **unchanged**.

---

## A. Binary chooses \(k\), then pure-even

\(k(n)\) from popcount / odd part / \(v_2\) / itinerary length; then pure-even monic \(\mathbb{Z}\) model.

| \(\mathcal{H}\) | Z-coeff rate | disc\(\square\) rate | \(\mathcal{H}\Rightarrow\) disc\(\square\)? |
|---------|--------------|----------------|---------------------|
| H_all | 1.00 | 1.00 | Yes (on scan) |
| H_core | 1.00 | 1.00 | Yes |
| H_small / short / bounded | 1.00 | 1.00 | Yes |

**Reading:** Once a pure-even \(\mathbb{Z}\) model exists for \(k(n)\), disc\(\square\) is classical. So \(\mathcal{H}\Rightarrow\) disc\(\square\) iff \(\mathcal{H}\Rightarrow\) “Z pure-even model for \(k(n)\)”. Composite lemma about \(F\), not HQCC necessity.

---

## B. \(F\to T(\ldots)\) only — disc\(\square\) \(\to 1\)?

Honest variants (M-deform, embed shape, full mix, binary homog seed, …):

| Best disc\(\square\) rate | Crit-2 signal (>0.5)? |
|----------------------|----------------------|
| 0.00 | **False** |

No template-only functor hits disc\(\square\) at high rate. Evenness obstruction on \(T\) stands.

---

## C. Tier 1.1 return

| Search | Result |
|--------|--------|
| New multi-param identical-square cuts | None |
| Homogenisation no-\(x^2\) (known) | 1-param, beyond BJ-embed when \(p\neq 0\) |
| Pure-even envelope | 2-param, inside BJ-embed |
| Bilinear cuts | all fail |

No new Crit-2 fragment.

---

Script: `tier12_sharp_next.py`
