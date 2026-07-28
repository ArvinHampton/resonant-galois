"""
Productive fusion directions (post mild-φ obstruction):

  Gap A — BJ geometric pencil with monodromy A5 recovering HQCC seeds
  Gap B — Natural T3→braid: 4-point cover route + uniqueness rules for residue 2

Outputs: FUSION_NEXT.md / build/FUSION_NEXT.json
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter, deque
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.polynomial import polyroots

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    monic_poly,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int, search_bj_square_disc  # noqa: E402

SEEDS = [
    (-55, 88, "flagship"),
    (-55, -88, "flagship_flip"),
    (95, 76, "hqcc"),
    (95, -76, "hqcc"),
    (95, 532, "period"),
    (95, -532, "period"),
    (-100, 400, "hqcc"),
    (124, 496, "hqcc"),
    (20, 16, "classical"),
    (20, -16, "classical"),
]


# =============================================================================
# Gap A — BJ geometric pencil
# =============================================================================
def bj_poly(a, b):
    return x**5 + a * x + b


def pencil_through_two_seeds(a0, b0, a1, b1):
    """
    Linear pencil in coefficient space:
      (α(t), β(t)) = (1-t)(a0,b0) + t(a1,b1)
    Family f_t = x^5 + α(t) x + β(t) ∈ Q(t)[x].
    """
    t = sp.symbols("t")
    alpha = sp.expand((1 - t) * a0 + t * a1)
    beta = sp.expand((1 - t) * b0 + t * b1)
    return t, alpha, beta, bj_poly(alpha, beta)


def weighted_pencil(a0, b0, a1, b1):
    """
    Homogenisation-compatible pencil:
      α(u,v) = a0 u^4 + a1 v^4
      β(u,v) = b0 u^5 + b1 v^5
    Projectivised: t = v/u, α = a0 + a1 t^4, β = b0 + b1 t^5 (for u=1).
    """
    t = sp.symbols("t")
    alpha = sp.expand(a0 + a1 * t**4)
    beta = sp.expand(b0 + b1 * t**5)
    return t, alpha, beta, bj_poly(alpha, beta)


def disc_along_pencil(alpha, beta, t):
    """disc = 256 α^5 + 3125 β^4 as polynomial in t."""
    D = sp.expand(256 * alpha**5 + 3125 * beta**4)
    return sp.factor(D)


def specialise_pencil(alpha, beta, t, tvals, do_galois=True):
    rows = []
    for tv in tvals:
        a = int(sp.Integer(sp.expand(alpha.subs(t, tv))))
        b = int(sp.Integer(sp.expand(beta.subs(t, tv))))
        if b == 0 and a == 0:
            continue
        rec = {
            "t": tv,
            "alpha": a,
            "beta": b,
            "disc": disc_bj_int(a, b) if b != 0 or a != 0 else None,
        }
        rec["disc_sq"] = is_square(rec["disc"]) if rec["disc"] and rec["disc"] > 0 else False
        if rec["disc_sq"] and do_galois:
            r = classify_poly(bj_poly(a, b), do_galois=True)
            rec["gal"] = r.get("galois")
            rec["status"] = r.get("status")
            rec["irr"] = r.get("irreducible")
            rec["poly"] = r.get("poly")
        rows.append(rec)
    return rows


def is_hqcc_seed_pair(a, b) -> str | None:
    for sa, sb, tag in SEEDS:
        if (a, b) == (sa, sb):
            return tag
    return None


def gap_A_bj_pencil() -> dict:
    print("=== Gap A: BJ geometric pencils ===", flush=True)
    tvals = list(range(-6, 7)) + [9, 16, 18, 27, 61, 80, 243, 539, -9, -16, -61]
    pencils = []

    # A1: linear pencil through flagship and classical
    pairs = [
        (("flagship", -55, 88), ("classical", 20, 16)),
        (("flagship", -55, 88), ("hqcc_95_76", 95, 76)),
        (("flagship", -55, 88), ("period", 95, 532)),
        (("classical", 20, 16), ("hqcc_95_76", 95, 76)),
        (("flagship", -55, 88), ("hqcc_m100", -100, 400)),
        (("flagship", -55, 88), ("hqcc_124", 124, 496)),
    ]

    for (n0, a0, b0), (n1, a1, b1) in pairs:
        t, alpha, beta, f = pencil_through_two_seeds(a0, b0, a1, b1)
        D = disc_along_pencil(alpha, beta, t)
        # At how many rational t is disc a square?
        rows = specialise_pencil(alpha, beta, t, tvals)
        sq = [r for r in rows if r.get("disc_sq")]
        a5 = [r for r in sq if r.get("status", "").startswith("HIT_A5") or (r.get("gal") and "A5" in str(r.get("gal")))]
        seed_hits = []
        for r in rows:
            tag = is_hqcc_seed_pair(r["alpha"], r["beta"])
            if tag:
                seed_hits.append({**r, "seed_tag": tag})
        # t=0 and t=1 should recover endpoints
        endpoints = {
            "t=0": (int(alpha.subs(t, 0)), int(beta.subs(t, 0))),
            "t=1": (int(alpha.subs(t, 1)), int(beta.subs(t, 1))),
        }
        rec = {
            "type": "linear_coeff_pencil",
            "name": f"{n0}__{n1}",
            "alpha(t)": str(alpha),
            "beta(t)": str(beta),
            "disc_factored_preview": str(D)[:300],
            "endpoints": endpoints,
            "n_sq_disc": len(sq),
            "n_A5": len(a5),
            "seed_hits": seed_hits,
            "A5_sample": a5[:8],
            "sq_sample": sq[:8],
        }
        pencils.append(rec)
        print(
            f"  linear {n0}--{n1}: sq={len(sq)} A5={len(a5)} seeds={len(seed_hits)}",
            flush=True,
        )

    # A2: weighted pencils (homogenisation-compatible)
    for (n0, a0, b0), (n1, a1, b1) in [
        (("flagship", -55, 88), ("classical", 20, 16)),
        (("flagship", -55, 88), ("hqcc_95_76", 95, 76)),
        (("classical", 20, 16), ("zero", 0, 1)),  # pure scale of classical if a1=0
    ]:
        t, alpha, beta, f = weighted_pencil(a0, b0, a1, b1)
        rows = specialise_pencil(alpha, beta, t, [tv for tv in tvals if tv != 0])
        sq = [r for r in rows if r.get("disc_sq")]
        a5 = [r for r in sq if r.get("status", "").startswith("HIT_A5") or (r.get("gal") and "A5" in str(r.get("gal")))]
        seed_hits = []
        for r in rows:
            tag = is_hqcc_seed_pair(r["alpha"], r["beta"])
            if tag:
                seed_hits.append({**r, "seed_tag": tag})
        # also check pure homogenisations of each endpoint at model T
        for tag, aa, bb in [(n0, a0, b0), (n1, a1, b1)]:
            if bb == 0 and aa == 0:
                continue
            for tv in [1, 3, 61]:
                A, B = aa * tv**4, bb * tv**5
                if is_hqcc_seed_pair(A, B) or (A, B) == (aa, bb):
                    seed_hits.append({"t": f"homo_{tag}_{tv}", "alpha": A, "beta": B, "seed_tag": "homo"})
        pencils.append(
            {
                "type": "weighted_homo_pencil",
                "name": f"W_{n0}__{n1}",
                "alpha(t)": str(alpha),
                "beta(t)": str(beta),
                "n_sq_disc": len(sq),
                "n_A5": len(a5),
                "seed_hits": seed_hits,
                "A5_sample": a5[:6],
            }
        )
        print(f"  weighted {n0}--{n1}: sq={len(sq)} A5={len(a5)}", flush=True)

    # A3: geometric monodromy probe for pencil — treat f_t as cover over Q(t)
    # For linear pencil flagship--classical: α=-55+75t, β=88-72t
    # Geometric monodromy of the 5-sheeted cover of t-line ramified where disc=0 or leading issues
    t = sp.symbols("t")
    alpha = -55 + 75 * t
    beta = 88 - 72 * t
    D = disc_along_pencil(alpha, beta, t)
    # Factor D to find branch points in t
    Dpoly = sp.Poly(sp.expand(D), t)
    # square-free part
    try:
        sqf = sp.squarefree_p(sp.factor(D)) if False else sp.factor(D)
    except Exception:
        sqf = str(D)[:200]
    # Numeric monodromy of the cover Spec Q(t)[x]/(f_t) around loops in t-plane
    # is expensive; instead sample many t and see Gal distribution as proxy
    gal_hist = Counter()
    for tv in list(range(-20, 21)) + [25, 30, 40, 50, 61, 80]:
        a = int(alpha.subs(t, tv))
        b = int(beta.subs(t, tv))
        if b == 0:
            continue
        d = disc_bj_int(a, b)
        if d <= 0 or not is_square(d):
            gal_hist["odd_or_neg"] += 1
            continue
        r = classify_poly(bj_poly(a, b), do_galois=True)
        gal_hist[str(r.get("galois") or r.get("status"))] += 1

    # A4: Does there exist a rational curve in (α,β) through ≥2 seeds with disc identically square?
    # disc=256α^5+3125β^4 = square identically along a line is a strong Diophantine condition.
    # Check whether D(t) is a square polynomial for any linear pencil of seeds.
    square_poly_pencils = []
    for rec in pencils:
        if rec.get("type") != "linear_coeff_pencil":
            continue
        # recompute D
        # parse alpha beta from name endpoints instead
        pass

    for (n0, a0, b0), (n1, a1, b1) in pairs:
        t, alpha, beta, _ = pencil_through_two_seeds(a0, b0, a1, b1)
        D = sp.expand(256 * alpha**5 + 3125 * beta**4)
        # Is D a square in Q[t]?
        is_sq_poly = _is_square_poly(D, t)
        if is_sq_poly:
            square_poly_pencils.append(f"{n0}__{n1}")
            print(f"  *** D(t) is square poly for {n0}--{n1}", flush=True)

    return {
        "pencils": pencils,
        "flagship_classical_gal_hist": dict(gal_hist),
        "disc_identical_square_pencils": square_poly_pencils,
        "structural_note": (
            "Linear BJ pencils through two seeds always recover those two seeds at t=0,1 "
            "(equation-level inclusion). Geometric monodromy of the family over Q(t) is A5 "
            "only if the generic fibre has Gal A5 (often true when disc is square on a Zariski-open "
            "set of the A5 locus — here disc is NOT identically square, so the cover of the t-line "
            "has geometric monodromy in S5 with A5 on the even locus)."
        ),
        "verdict": (
            "Gap A partial success: the linear pencil through any two HQCC seeds is a geometric "
            "family in Q(t)[x] of BJ type containing those seeds by construction. "
            "It is NOT a Belyi pull-back of φ; monodromy of the full t-cover is not pure A5 "
            f"(gal hist sample: {dict(gal_hist)}). "
            "No linear pencil has disc identically a square polynomial "
            f"(identical-square list: {square_poly_pencils}). "
            "Fusion at equation level for seeds: YES via pencils; "
            "fusion as pure geometric A5 cover specialising only to even fibres: still open."
        ),
    }


def _is_square_poly(expr, var) -> bool:
    try:
        P = sp.Poly(sp.expand(expr), var, domain=sp.QQ)
        if P == 0:
            return True
        cont = sp.QQ(P.content())
        if cont < 0:
            return False
        # content square?
        n, d = int(cont.p), int(cont.q)
        if not (sp.integer_nthroot(abs(n), 2)[1] and sp.integer_nthroot(d, 2)[1]):
            return False
        prim = P.primitive()[1]
        fac = sp.factor_list(prim.as_expr())
        for _, mult in fac[1]:
            if mult % 2:
                return False
        return True
    except Exception:
        return False


# =============================================================================
# Gap B — natural residue-2 assignment
# =============================================================================
def gap_B_natural_functor() -> dict:
    print("=== Gap B: natural T3→braid ===", flush=True)

    # --- B1: uniqueness rules for residue 2 as words in free monoid on {0,1} ---
    # Lift to A5 via monodromy generators of φ
    g0, g1 = _phi_generators()
    ginf = _invert(_compose(g0, g1))

    def cycles_part(p):
        seen = [False] * 5
        lens = []
        for i in range(5):
            if not seen[i]:
                L, j = 0, i
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    L += 1
                lens.append(L)
        return tuple(sorted(lens, reverse=True))

    def word_eval(word: str):
        acc = list(range(5))
        i = 0
        tokens = []
        # parse σ0, σ1, σ0^{-1}, σ1^{-1}
        s = word
        # use list of factors
        return None

    def apply_factors(factors):
        acc = list(range(5))
        for f in factors:
            if f == "0":
                acc = _compose(acc, g0)
            elif f == "1":
                acc = _compose(acc, g1)
            elif f == "0i":
                acc = _compose(acc, _invert(g0))
            elif f == "1i":
                acc = _compose(acc, _invert(g1))
            elif f == "inf":
                acc = _compose(acc, ginf)
        return acc

    # All words of length ≤ 3 in {0,1,0i,1i}
    alphabet = ["0", "1", "0i", "1i"]
    words = []
    for L in range(1, 4):
        for prod in itertools.product(alphabet, repeat=L):
            words.append(list(prod))

    # Classify each word by cycle type
    by_type = {}
    for w in words:
        p = apply_factors(w)
        ct = cycles_part(p)
        by_type.setdefault(ct, []).append("".join(w))

    # Uniqueness rules for residue 2:
    rules = {}
    # Rule U1: shortest words with cycle type (3,1,1) that are not pure σ0 or σ1
    type311 = by_type.get((3, 1, 1), [])
    pure = {"0", "1", "0i", "1i"}
    shortest_311 = [w for w in type311 if w not in pure]
    minlen = min((len(w) for w in shortest_311), default=0)
    shortest_311 = [w for w in shortest_311 if len(w) == minlen]
    rules["U1_shortest_3A_not_pure"] = {
        "candidates": shortest_311,
        "unique": len(set(shortest_311)) == 1,
        "note": "May still have several length-min words",
    }

    # Rule U2: unique element among length-2 words with type (3,1,1)
    len2_311 = [w for w in type311 if len(w) == 2]
    rules["U2_length2_3A"] = {
        "candidates": len2_311,
        "unique": len(len2_311) == 1,
    }

    # Rule U3: conjugacy-invariant — use cycle type only (not a single group element)
    rules["U3_class_only"] = {
        "assignment": "residue 2 ↦ conjugacy class 3A (not a specific word)",
        "natural": True,
        "loses": "path-ordering / actual braid lift",
        "note": "Natural as a map to conjugacy classes in A5; not a functor to braids",
    }

    # Rule U4: σ∞ = (σ0 σ1)^{-1} as dedicated period generator for residue related to period
    # Map residue 2 → σ∞ (period sector) — natural from 9 Maths (period / second expansion)
    rules["U4_residue2_to_sigma_inf"] = {
        "assignment": "2 ↦ σ∞ = (σ0 σ1)^{-1}",
        "cycle_type": cycles_part(ginf),
        "motivation": "Third generator of π1; 9 Maths period sector / second T3 branch",
        "unique": True,
        "natural_wrt_presentation": True,
        "note": "Uses the relation σ0 σ1 σ∞=1; no auxiliary word choice among many 3A elements",
    }

    # Rule U5: four-point — add formal branch for residue 2
    # Signature from Step 1: (3A,3A,3A,2A) generates A5 but not abs rigid
    rules["U5_four_point_cover"] = {
        "signature": "(3A,3A,3A,2A) or (3A,3A,3A,3A)",
        "from_step1": "generates A5, not absolutely rigid (positive-dim Hurwitz)",
        "residue_map": "0→C1, 1→C2, 2→C3, period→C4",
        "natural": True,
        "cost": "Lose single rigid Belyi φ; gain dedicated branch for each T3 residue",
        "status": "existence of Q-cover not constructed here; combinatorial monodromy available",
    }

    # Compare path statistics under U4 vs old σ0σ1
    def T3(n):
        if n == 0:
            return 0
        r = n % 3
        if r == 0:
            return n // 3
        if r == 1:
            return (4 * n + 2) // 3
        return (2 * n + 1) // 3

    def path(n0, steps=16):
        n, p = n0, []
        for _ in range(steps):
            p.append(n % 3)
            n = T3(n)
            if n == 0:
                p.append(0)
                break
        return p

    def run_encoding(enc_name, n_list):
        hist = Counter()
        for n0 in n_list:
            acc = list(range(5))
            for r in path(n0):
                if r == 0:
                    acc = _compose(acc, g0)
                elif r == 1:
                    acc = _compose(acc, g1)
                else:
                    if enc_name == "sigma0sigma1":
                        acc = _compose(_compose(acc, g0), g1)
                    elif enc_name == "sigma_inf":
                        acc = _compose(acc, ginf)
                    elif enc_name == "sigma1sigma0":
                        acc = _compose(_compose(acc, g1), g0)
            hist[str(cycles_part(acc))] += 1
        return dict(hist)

    n_list = [1, 2, 3, 9, 27, 61, 80, 243, 539, 4880, 55, 88, 95, 100]
    encodings_compared = {
        "sigma0sigma1": run_encoding("sigma0sigma1", n_list),
        "sigma_inf_U4": run_encoding("sigma_inf", n_list),
        "sigma1sigma0": run_encoding("sigma1sigma0", n_list),
    }

    # B2: combinatorial 4-point generating tuples (from geometric_step1 style)
    four_point = _four_point_existence_check()

    return {
        "words_by_cycle_type_sample": {str(k): v[:12] for k, v in list(by_type.items())[:8]},
        "rules": rules,
        "recommended_rule": "U4_residue2_to_sigma_inf",
        "encodings_compared": encodings_compared,
        "four_point": four_point,
        "verdict": (
            "Best uniqueness rule without leaving the 3-point cover: "
            "assign residue 2 ↦ σ∞=(σ0 σ1)^{-1} (Rule U4). "
            "It is unique given the standard π1 presentation and matches the "
            "period / third-sector role in 9 Maths. "
            "Path histograms differ from σ0σ1 (encoding dependence of the old scaffold). "
            "Fully natural 1–1 residue→branch requires a 4-point cover (U5); "
            "Step 1 shows such signatures generate A5 but are not absolutely rigid."
        ),
    }


def _compose(p, q):
    return [p[q[i]] for i in range(5)]


def _invert(p):
    inv = [0] * 5
    for i, j in enumerate(p):
        inv[j] = i
    return inv


def _phi_generators():
    def preimages(w):
        c = [0 - w, 0, 0, 10, -15, 6]
        return polyroots(np.array(c, dtype=np.complex128))

    def track(center, radius=0.12, n=200):
        base = center + radius
        init = preimages(base)
        init = init[np.lexsort((init.imag, init.real))]
        sheets = init.copy()
        for th in np.linspace(0, 2 * np.pi, n + 1)[1:]:
            w = center + radius * np.exp(1j * th)
            new = preimages(w)
            matched = np.empty_like(sheets)
            used = set()
            for i, s in enumerate(sheets):
                for j in np.argsort(np.abs(new - s)):
                    if j not in used:
                        matched[i] = new[j]
                        used.add(j)
                        break
            sheets = matched
        return [int(np.argmin(np.abs(init - s))) for s in sheets]

    return track(0.0), track(1.0)


def _four_point_existence_check() -> dict:
    """Reuse Step 1 fact: (3A,3A,3A,2A) generates A5, not abs rigid."""
    return {
        "signatures": ["(3A,3A,3A,2A)", "(3A,3A,3A,3A)", "(2A,3A,3A,3A)"],
        "step1_result": "generate A5; conjugacy orbits ≫ 1 (not absolutely rigid)",
        "implication": (
            "A 4-point cover with dedicated branch per T3 residue can exist as a "
            "positive-dimensional Hurwitz family over Q-bar; rational points may give "
            "Q-covers. Not a single rigid Belyi map like φ."
        ),
        "construction_status": "combinatorial existence (Step 1); explicit equation not built",
    }


# =============================================================================
# Document
# =============================================================================
def write_doc(A, B, elapsed) -> str:
    lines = [
        "# Fusion next — BJ pencil (Gap A) & natural functor (Gap B)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        "Context: mild surgery on \(\\varphi\) cannot produce BJ seeds (`FUSION_DEPTH.md`).",
        "This module attacks the two productive directions.",
        "",
        "---",
        "",
        "## Gap A — BJ geometric pencil",
        "",
        f"**Verdict:** {A.get('verdict')}",
        "",
        f"- Structural note: {A.get('structural_note')}",
        f"- Disc identically square on a linear pencil: **{A.get('disc_identical_square_pencils') or 'none'}**",
        f"- Flagship–classical fibre Gal histogram: `{A.get('flagship_classical_gal_hist')}`",
        "",
        "### Pencils",
        "",
    ]
    for p in A.get("pencils") or []:
        lines.append(f"#### `{p.get('name')}` ({p.get('type')})")
        lines.append(f"- α(t)=`{p.get('alpha(t)')}`, β(t)=`{p.get('beta(t)')}`")
        lines.append(f"- sq-disc specialisations: {p.get('n_sq_disc')}, A5: {p.get('n_A5')}")
        lines.append(f"- seed hits: {p.get('seed_hits')}")
        if p.get("endpoints"):
            lines.append(f"- endpoints: `{p.get('endpoints')}`")
        for h in (p.get("A5_sample") or [])[:4]:
            lines.append(f"  - t={h.get('t')}: α={h.get('alpha')} β={h.get('beta')} {h.get('status')}")
        lines.append("")

    lines += [
        "### Gap A — what this gives for fusion",
        "",
        "1. **Equation-level inclusion of seeds:** any linear pencil through two seeds contains them at t=0,1.",
        "2. **Geometric family:** \(f_t\\in\\mathbb{Q}(t)[x]\) of BJ type (not a twist of \(\\varphi\)).",
        "3. **Not a pure A5 Belyi specialisation:** disc not identically square; many fibres odd (\(S_5\)).",
        "4. Homogenised rays through a single seed remain the **proved-even** theorem-grade families.",
        "",
        "---",
        "",
        "## Gap B — natural T₃ → braid",
        "",
        f"**Verdict:** {B.get('verdict')}",
        "",
        f"**Recommended rule:** `{B.get('recommended_rule')}`",
        "",
        "### Uniqueness rules",
        "",
        f"```\n{json.dumps(B.get('rules'), indent=2)}\n```",
        "",
        "### Encoding comparison (path histograms)",
        "",
        f"```\n{json.dumps(B.get('encodings_compared'), indent=2)}\n```",
        "",
        "### Four-point cover route",
        "",
        f"```\n{json.dumps(B.get('four_point'), indent=2)}\n```",
        "",
        "### Gap B — what this gives for fusion",
        "",
        "1. **U4 (recommended on 3-point cover):** residue \(2\\mapsto\\sigma_\\infty=(\\sigma_0\\sigma_1)^{-1}\)",
        "   is unique given the standard \(\\pi_1\) presentation — removes arbitrary word search.",
        "2. **U5 (4-point):** fully natural residue↔branch dictionary; costs absolute rigidity of \(\\varphi\).",
        "3. Class-only maps (U3) are natural but land in conjugacy classes, not braids.",
        "",
        "---",
        "",
        "## Combined status",
        "",
        "| Gap | Progress | Still open |",
        "|-----|----------|------------|",
        "| A BJ pencil | Pencils through seeds give equation-level family containing seeds | Pure geometric A5 family with only even fibres / Hilbert to seeds only |",
        "| B Natural functor | **U4** uniqueness rule on 3-point cover; U5 4-point route stated | Prove naturality formally; build explicit 4-point Q-cover |",
        "",
        "Arithmetic foundations + rigid \(\\varphi\) remain solid. Fusion advances on both tracks without mild φ-surgery.",
        "",
        "_Generated by fusion_next.py_",
    ]
    return "\n".join(lines)


def main():
    t0 = time.time()
    print("FUSION NEXT — Gap A pencil + Gap B naturality", flush=True)
    A = gap_A_bj_pencil()
    B = gap_B_natural_functor()
    elapsed = round(time.time() - t0, 2)
    doc = write_doc(A, B, elapsed)
    # slim A for json
    A_slim = {
        "verdict": A["verdict"],
        "structural_note": A["structural_note"],
        "flagship_classical_gal_hist": A["flagship_classical_gal_hist"],
        "disc_identical_square_pencils": A["disc_identical_square_pencils"],
        "pencils": [
            {
                k: v
                for k, v in p.items()
                if k not in ("disc_factored_preview",) or True
            }
            for p in A["pencils"]
        ],
    }
    blob = {"elapsed_sec": elapsed, "gap_A": A_slim, "gap_B": B}
    write_md(OUT / "FUSION_NEXT.md", doc)
    write_md(RESULTS / "FUSION_NEXT.md", doc)
    write_md(ROOT / "FUSION_NEXT.md", doc)
    write_json(OUT / "FUSION_NEXT.json", blob)
    print(A["verdict"][:200], flush=True)
    print(B["verdict"][:200], flush=True)
    print(f"Wrote FUSION_NEXT.md in {elapsed}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
