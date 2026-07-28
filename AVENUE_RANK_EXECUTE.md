# Avenue rank execute (1→7)

_Elapsed: 32.25s_

**Verdict:** Avenues 1–7 executed in rank order (32.25s). Geometric multi-k Nielsen hit: False. Arithmetic multi-k: True. Best geometric progress: H^rd≅P^1_s + numeric (3,1,1)^4 covers; closed form f_s∈Q(s)[x] still open. Best arithmetic: 2-param envelope + cross-k paths.

---

## Scorecard

| Rank | Avenue | Success flag | Summary |
|-----:|--------|:------------:|---------|
| 1 | Better 3A^4 resolvent | False | Symmetric-pole hits=1; p2-coordinate hits=3; p2-ansatz closed sols at s=2: 0; closed form f_s in Q(s |
| 2 | Next shortlist g=0 | True | Prior orbits: 2A3A^3 size 96, 2A^2 3A^2 size 108, g=0 lookup. Explicit Nielsen cover equations for t |
| 3 | Pure-even A5 strata | True | 2-param pure-even envelope disc_id=True; A5 hits=21/21 lattice samples; this is the positive-dim pur |
| 4 | Other rigid triples | False | Rigid triples over Q remain blocked for pure-even fibres (phi recheck even=0/7). 2A3A5* not over Q.  |
| 5 | Base change descent | False | disc=5*□ over Q (proved=True); square over Q(√5) (proved=True). No descent of evenness to Q; no HQCC |
| 6 | Higher-rank r≥5 | False | r=5 filter-pass types: 34 (dim_rd=2). Smallest class-spaces start at 414720. No explicit equations;  |
| 7 | Envelope geometric lift | False | Envelope multi-k arithmetic confirmed (hits=[{'t': '0', 'tag': 'flagship', 'k': '-8/5'}, {'t': '1',  |

---

## 1. Better rational coordinate / resolvent for \(3A^4\)

- Symmetric-pole hits: 1
- p2-coordinate hits: 3
- s(p2) polyfit: None
- p2-ansatz results: [{'ansatz': 'p2=s', 'solved_s': [-1.0], 'n': 1}, {'ansatz': 'p2=1-s', 'solved_s': [], 'n': 0}, {'ansatz': 'p2=-s', 'solved_s': [], 'n': 0}, {'ansatz': 'p2=s+1', 'solved_s': [], 'n': 0}, {'ansatz': 'p2=2-s', 'solved_s': [], 'n': 0}, {'ansatz': 'p2=-1', 'solved_s': [-3.505910655822549e-20, -9.456971697109125e-16, -0.9999999999999987, 0.9999999999999998, -1.8177838708916986e-15, -1.0000000000000002], 'n': 6}, {'ansatz': 'p2=phi2', 'solved_s': [0.9999999999999998, 0.9999999999999998, 5.033457384250216e-18, 0.9999999999999998, 1.0, 9.3390053937296e-18], 'n': 6}]
- Closed form \(f_s\in\mathbb{Q}(s)[x]\): **None**
- Envelope control multi-k: True hits=[{'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'classical', 'k': '4/5', 't': '1'}]
- **Symmetric-pole hits=1; p2-coordinate hits=3; p2-ansatz closed sols at s=2: 0; closed form f_s in Q(s)[x]: still open; envelope multi-k control=True**

---

## 2. Next shortlist genus-0 (\(2A3A^3\), \(2A^2 3A^2\))

- Prior orbits: `{'2A,3A,3A,3A': {'orbits': 1, 'sizes': [96], 'genus_lookup': 0}, '2A,2A,3A,3A': {'orbits': 1, 'sizes': [108], 'genus_lookup': 0}}`
- Family tests: `[{'id': 'path_flag_classical', 'n_even_irr': 14, 'catalogue_k': ['-8/5', '4/5'], 'multi_cat': True, 'hits': [{'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'classical', 'k': '4/5', 't': '1'}, {'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'classical', 'k': '4/5', 't': '1'}]}, {'id': 'path_flag_lsw', 'n_even_irr': 13, 'catalogue_k': ['-4', '-8/5'], 'multi_cat': True, 'hits': [{'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'lsw_m100', 'k': '-4', 't': '1'}, {'tag': 'flagship', 'k': '-8/5', 't': '0'}, {'tag': 'lsw_m100', 'k': '-4', 't': '1'}]}, {'id': 'slice_-4', 'n_even_irr': 13, 'catalogue_k': [], 'multi_cat': False, 'hits': []}, {'id': 'slice_-8/5', 'n_even_irr': 2, 'catalogue_k': [], 'multi_cat': False, 'hits': []}, {'id': 'slice_4/5', 'n_even_irr': 2, 'catalogue_k': [], 'multi_cat': False, 'hits': []}]`
- Multi-cat families: ['path_flag_classical', 'path_flag_lsw']
- **Prior orbits: 2A3A^3 size 96, 2A^2 3A^2 size 108, g=0 lookup. Explicit Nielsen cover equations for these types not closed-formed this run. Proxy pure-even paths multi_cat=['path_flag_classical', 'path_flag_lsw']. Geometric multi-k for 2A* still open.**

---

## 3. Positive-dimensional pure-even \(A_5\) strata

- Envelope 2-param: disc_id=True, A5=21/21
- Sample points: [{'m': '1', 'k': '-4', 'a': -2869, 'b': 11476}, {'m': '1', 'k': '4', 'a': -2869, 'b': -11476}, {'m': '2', 'k': '-4', 'a': -2101, 'b': 8404}, {'m': '2', 'k': '4', 'a': -2101, 'b': -8404}, {'m': '3', 'k': '-4', 'a': -821, 'b': 3284}, {'m': '3', 'k': '4', 'a': -821, 'b': -3284}]
- Even surface: `{'equation': '256*alpha**5 + 3125*beta**4 = gamma**2', 'dim_affine_cone': 2, 'k_ray_foliation': 'envelope = ruled by pure-even k-rays'}`
- **2-param pure-even envelope disc_id=True; A5 hits=21/21 lattice samples; this is the positive-dim pure-even A5 arithmetic stratum. Multi-k by construction when k varies.**

---

## 4. Other rigid triples

- **(3A,3A,5A)** (Q): even_over_Q=False — disc monic(phi-t)=5*square
- **(3A,3A,5B)** (Q): even_over_Q=False — same disc obstruction up to automorphism
- **(2A,3A,5A)** (Q(2^{1/5},3^{1/5})): even_over_Q=False — not over Q; prior numeric even scan empty
- **(2A,3A,5B)** (same): even_over_Q=False — same as 2A3A5A
- phi recheck: {'tested_irr': 7, 'even': 0}
- **Rigid triples over Q remain blocked for pure-even fibres (phi recheck even=0/7). 2A3A5* not over Q. Likelihood of multi-k via rigid triples: low.**

---

## 5. Base change + descent

- disc theorem: {'equals_5_times_square': True, 'square_in_Q_sqrt5': True}
- descent_to_Q: False, lattice_recovery: False
- **disc=5*□ over Q (proved=True); square over Q(√5) (proved=True). No descent of evenness to Q; no HQCC Z-lattice recovery. Side route only.**

---

## 6. Higher-rank rigid systems (\(r\\ge 5\))

- r=5 filter-pass types: 34
- smallest: [{'type': '3A,5A,5A,5A,5A', 'class_tuple_space': 414720, 'dim_rd': 2}, {'type': '3A,5A,5A,5A,5B', 'class_tuple_space': 414720, 'dim_rd': 2}, {'type': '3A,5A,5A,5B,5B', 'class_tuple_space': 414720, 'dim_rd': 2}, {'type': '3A,5A,5B,5B,5B', 'class_tuple_space': 414720, 'dim_rd': 2}, {'type': '3A,5B,5B,5B,5B', 'class_tuple_space': 414720, 'dim_rd': 2}]
- dim: {'r5': 2, 'r6': 3}
- **r=5 filter-pass types: 34 (dim_rd=2). Smallest class-spaces start at 414720. No explicit equations; multi-k likelihood speculative. Effort very high.**

---

## 7. Geometric lift of the envelope

- disc factors: [('81*t**4 - 216*t**3 + 216*t**2 - 96*t + 11', 4)]
- odd disc factors: 0
- catalogue hits: [{'t': '0', 'tag': 'flagship', 'k': '-8/5'}, {'t': '1', 'tag': 'classical', 'k': '4/5'}, {'t': '1/3', 'tag': 'classical_m', 'k': '-4/5'}]
- multi-k arithmetic: True
- Nielsen 3A^4?: None
- Envelope is pure-even BJ over Q(t) with multi-k Hilbert hits. Identification with Ni(A5,C_3^4) would require matching monodromy generators as 3-cycles at four branch values. Disc being a full square means all finite branch multiplicities of disc are even; cycle types need monodromy computation of the Galois closure, not done here. Speculative: envelope may be a multi-parameter pullback or a different Nielsen class with pure-even arithmetic specialisations.
- **Envelope multi-k arithmetic confirmed (hits=[{'t': '0', 'tag': 'flagship', 'k': '-8/5'}, {'t': '1', 'tag': 'classical', 'k': '4/5'}, {'t': '1/3', 'tag': 'classical_m', 'k': '-4/5'}]). Odd disc factors=0. Nielsen 3A^4 identification: open/speculative.**

---

## Global conclusions

1. **Geometric multi-\(k\) (Nielsen-labelled): still open.**
2. **Arithmetic multi-\(k\): solid** via 2-param envelope and cross-\(k\) paths.
3. **Avenue 1** remains the highest-leverage geometric attack (genus 0 is favourable)
   but closed-form \(f_s\in\mathbb{Q}(s)[x]\) needs more elimination / descent work.
4. **Avenues 4–5** are effectively closed as multi-k routes (blocked / side-only).
5. **Avenue 6** is enumerated at type-count level only; equations out of scope.
6. **Avenue 7**: envelope is the arithmetic multi-k object; geometric Nielsen ID open.

### Recommended single next move

Finish Avenue 1: resultant/Gröbner elimination of the triple-root ideal to a
plane model of \((c:p_2:r_1:r_2)\) over \(\mathbb{Q}(s)\), or a deg-5 resolvent
after a radical extension of \(\mathbb{Q}(s)\) with controlled descent.

_Generated by avenue_rank_execute.py_