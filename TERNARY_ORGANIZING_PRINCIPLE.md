# Ternary organizing principle

_The single design that sits under the generative success — structural reading, not a proved necessity theorem._

---

## Status

| Claim | Status |
|-------|--------|
| **Organizing principle** (ternary branching ↔ generators of \(A_n\)) | **Locked as structural reading** of generative success |
| **Pure-even multi-\(k\)** | **Finished + published core** (`PURE_EVEN_MULTI_K.md`) |
| **Necessity theorem** | **Open / paused** as citable claim (`NECESSITY_THEOREM.md`) |

Together with pure-even multi-\(k\), this principle is the **citable result** of the programme. Forced \(A_n\) from HQCC axioms alone remains open research.

This is **not** “classical BJ + a dense lattice.” It is that the **same ternary branching appears in four places at once**.

---

## 1. Group theory (the root fact)

For \(n\geq 3\),

\[
A_n=\langle\text{all 3-cycles}\rangle.
\]

Every even permutation is a product of 3-cycles. So any arithmetic machine that systematically produces **enough 3-cycles in Frobenius**, together with **even monodromy** (square disc), is aimed at the **generators of \(A_n\)**.

That is classical. Resonant-Gal did not invent it. It **uses** it.

Operational form for \(n=5\): irr + disc □ + Frobenius type \((3,1,1)\) \(\Rightarrow\mathrm{Gal}=A_5\).

---

## 2. Four faces of the same ternary branching

| Face | Where ternary appears | Role |
|------|----------------------|------|
| **Dynamics** | HQCC / T₃ map (residue mod 3) | Branches \(\mathbb{N}\) by \(0,1,2\bmod 3\) |
| **Lattice** | Generators \(\{3,9,27,61,80,243,539,\ldots\}\) | Coefficients built from order-3 data |
| **Matrices** | Template \(T(a,\ldots,f)\) with ternary couplings | Charpolys whose Frobenius prefer type \((3,1,1)\) |
| **Galois** | Cycle type \((3,1,1)\) in \(S_5\) | Generators of \(A_5\) once disc is square |

The pure-even multi-\(k\) theory then does one more thing: it forces the **sign** to be \(+1\) (disc □).

**Pipeline:**

```
ternary branching
    → order-3 lattice & templates
        → 3-cycles in Frobenius
            + pure-even (disc □)
                → A₅  (and, with enlargement, A₆)
```

That is why the HQCC T₃ lattice produces Resonant-Gal polynomials:

- **not** because BJ needs HQCC,  
- **not** only because the lattice is dense,  
- **but** because ternary branching is the **arithmetic avatar of the generators of \(A_n\)**.

---

## 3. Design mirror, sharpened

Classical binary Collatz-type maps are the mod-2 twin of the same idea:

| | Binary | Ternary (HQCC) |
|--|--------|----------------|
| Modulus | 2 | 3 |
| Contract | even | \(0\bmod 3\) |
| Expand | odd | \(1,2\bmod 3\) |
| Natural cycles in monodromy | 2-power structure | **3-cycles** |
| Target alternating groups | harder to force from 2-cycles alone | natural: \(A_n=\langle\text{3-cycles}\rangle\) |

The “mirror” is not mystical: **alternating groups want 3-cycles; ternary dynamics manufacture 3-cycles.** Binary dynamics do not sit in that groove as cleanly.

This is Candidate **C** in `RESONANT_ALGEBRAIC_CLOSURE.md` (design mirror), sharpened.

---

## 4. What this explains — and what it does not

### Explains

1. Why a **ternary** lattice, rather than a random integer lattice, is an efficient search space for \(A_5\).  
2. Why matrix templates with **order-3 entries** keep producing type \((3,1,1)\).  
3. Why the programme’s history (T₃ → resonant lattice → Resonant-Gal) is **one design**, not three unrelated tricks.  
4. Why Candidate C (mod-2 ↔ mod-3 design mirror) is the right language for the bigger picture.

### Does not explain (still open)

1. A **necessity theorem**: every object defined by HQCC axioms has \(\mathrm{Gal}=A_n\).  
2. Why **539** or **4880** should appear inside the disc identity (they need not; the identity is classical BJ).  
3. A **functor** from binary Collatz data to pure-even \(A_5\) polynomials.  
4. **Fusion** of the geometric Nielsen cover with the arithmetic seed lattice.

---

## 5. One-sentence bigger picture

> **HQCC T₃ produces Resonant-Gal polynomials because ternary branching is the dynamical and arithmetic form of the 3-cycles that generate \(A_n\); pure-even multi-\(k\) then locks the sign so those 3-cycles land inside \(A_n\) rather than \(S_n\).**

That is the missing organising principle. It is still a **structural reading of the generative success**, not a proved necessity theorem — but it is the right scale of picture.

---

## 6. Relation to other locks

| Doc | Relation |
|------|----------|
| `THEOREMS.md` / pure-even multi-\(k\) | **What** is proved (sign lock + operational \(A_5\)) |
| `NECESSITY_THEOREM.md` | What would make the picture **forced**, not generative |
| `RESONANT_ALGEBRAIC_CLOSURE.md` | Candidate A = closure machine; Candidate C = this mirror |
| `HQCC_MATRIX_TEMPLATES.md` | Ternary face in matrices; structure alone ⇏ disc □ |
| `T_SUBCLASS_IDENTICAL_SQUARE.md` | Crit 2 still fails beyond pure-even ansatz |

```
Group theory:  A_n = ⟨3-cycles⟩
        │
        ▼
Four faces of ternary branching  (this document)
        │
        ├── generative success  →  pure-even multi-k  (FINISHED)
        │
        └── necessity            →  Crit 1–3 still open
```

---

## 7. Locked slogans

1. **One principle under generative success:** ternary branching in dynamics, lattice, matrices, and Galois.  
2. **Pure-even multi-\(k\)** locks the sign; 3-cycles aim at the generators of \(A_n\).  
3. **Not** “BJ + dense lattice only.”  
4. **Still not** a necessity theorem.

_Generated as programme lock — ternary organizing principle._
