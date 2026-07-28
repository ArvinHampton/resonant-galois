"""
Q1: Which pure-even k ∈ R_n arise as cross-ratios of a cover over R_n?

Precise split of objects:
  s  = cross-ratio of four branch points  (Hurwitz / M_{0,4} coordinate)
  k  = β/α  pure-even ratio class of a BJ fibre x⁵+αx+β

Answer structure:
  A. Arithmetic: every k ∈ R_n is pure-even over R_n (field-agnostic identity).
  B. Trivial geometry of s: every s ∈ R_n \\ {0,1} is a cross-ratio over R_n.
  C. Non-trivial geometry of k: which k arise from covers / fibres over R_n
     — partial answer from known models + cosine / catalogue constraints.
  D. Open: closed map s ↦ k without f_s ∈ R_n(s)[x].

Output: PURE_EVEN_K_CROSS_RATIO.md / .json
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, is_square, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

# Catalogue pure-even k over Q (multi-seed)
CATALOGUE_K = [
    Fraction(-4),
    Fraction(4),
    Fraction(-8, 5),
    Fraction(8, 5),
    Fraction(4, 5),
    Fraction(-4, 5),
    Fraction(-12, 5),
    Fraction(12, 5),
    Fraction(-16, 5),
    Fraction(16, 5),
]


def Rn_degree(n: int) -> int:
    if n < 3:
        return 1
    return int(sp.totient(n) // 2)


def xi_minpoly(n: int) -> sp.Poly | None:
    if Rn_degree(n) > 12:
        return None
    xi = 2 * sp.cos(2 * sp.pi / n)
    return sp.Poly(sp.minpoly(xi, x), x, domain=sp.ZZ)


def element_in_Rn_Q(val: Fraction, n: int) -> dict:
    """Q ⊂ R_n always — every rational k is in every R_n."""
    return {
        "k": str(val),
        "in_Rn": True,
        "reason": "Q ⊂ R_n for all n ≥ 1",
        "minimal_field_over_Q": "Q",
        "degree": 1,
    }


def cosine_values_in_Rn(n: int) -> list[dict]:
    """
    Candidates: k = 2 cos(2π p/n) ∈ R_n (generator orbit).
    Also proper subfield cosines 2 cos(2π p/d) for d|n.
    """
    out = []
    divisors = [int(d) for d in sp.divisors(n) if d >= 3]
    for d in divisors:
        # Galois orbit of 2cos(2π/d): p with gcd(p,d)=1, p ≤ d/2
        seen = set()
        for p in range(1, d):
            if sp.gcd(p, d) != 1:
                continue
            p_red = min(p, d - p)  # cos symmetry
            if p_red in seen:
                continue
            seen.add(p_red)
            # exact algebraic for small d
            expr = 2 * sp.cos(2 * sp.pi * p_red / d)
            minpoly_k = None
            deg = None
            try:
                mp = sp.minpoly(expr, x)
                minpoly_k = str(mp)
                deg = int(sp.degree(mp))
            except Exception:
                pass
            # numeric
            knum = float(sp.N(expr, 20))
            out.append(
                {
                    "form": f"2cos(2π·{p_red}/{d})",
                    "d": d,
                    "p": p_red,
                    "in_Rn": d == n or (n % d == 0),  # R_d ⊂ R_n when d|n? 
                    # Actually Q(cos 2π/d) ⊂ Q(cos 2π/n) iff d|n under standard cyclotomic inclusions
                    # more carefully: Q(ζ_d)^+ ⊂ Q(ζ_n)^+ when d|n
                    "numeric": knum,
                    "minpoly": minpoly_k,
                    "degree": deg,
                    "equals_catalogue_k": any(
                        abs(knum - float(k)) < 1e-9 for k in CATALOGUE_K
                    ),
                }
            )
    # Fix in_Rn using d|n
    for e in out:
        e["in_Rn"] = n % e["d"] == 0 or e["d"] == n
    return out


def cross_ratio(a, b, c, d):
    """λ = (a-c)/(a-d) : (b-c)/(b-d)."""
    return sp.simplify(((a - c) / (a - d)) / ((b - c) / (b - d)))


def cross_ratios_from_Rn_points(n: int) -> dict:
    """
    Four points among {0,1,∞} ∪ {ξ-orbit sample} in P¹(R_n).
    Cross-ratios that land in R_n (always if points in P¹(R_n)).
    """
    xi = sp.symbols("xi")
    # points: 0, 1, infinity handled as limit, ξ, 1-ξ, -ξ, ξ² reduced
    mp = xi_minpoly(n)
    if mp is None:
        return {"n": n, "status": "deg_too_large"}

    # Use numeric embeddings for concrete cross-ratios
    # Algebraic: s = ξ (branch points 0,1,∞,ξ) — the standard M_{0,4} chart
    # Also s = (ξ - 0)/(ξ - 1) etc.
    configs = []
    # Standard: branch points {0,1,∞,s} with s ∈ R_n ⇒ s is the cross-ratio
    for s_name, s_expr in [
        ("xi", xi),
        ("1-xi", 1 - xi),
        ("-xi", -xi),
        ("xi/(xi-1)", xi / (xi - 1)),
        ("(xi-1)/xi", (xi - 1) / xi),
        ("1/xi", 1 / xi),
    ]:
        configs.append(
            {
                "branch_chart": "{0,1,∞,s}",
                "s_form": s_name,
                "s_in_Rn": True,
                "note": "By definition every s∈R_n\\{0,1} is a cross-ratio over R_n",
            }
        )

    # Four finite points from cosine orbit → cross-ratio
    cos_pts = []
    for p in range(1, min(n, 8)):
        if sp.gcd(p, n) != 1:
            continue
        cos_pts.append(2 * sp.cos(2 * sp.pi * p / n))
        if len(cos_pts) >= 4:
            break
    cr_samples = []
    if len(cos_pts) >= 4:
        a, b, c, d = cos_pts[:4]
        try:
            lam = cross_ratio(a, b, c, d)
            cr_samples.append(
                {
                    "points": "four 2cos(2π p/n) in orbit",
                    "cross_ratio": str(sp.simplify(lam)),
                    "numeric": float(sp.N(lam, 15)),
                    "minpoly": str(sp.minpoly(lam, x)),
                }
            )
        except Exception as e:
            cr_samples.append({"error": str(e)})

    return {
        "n": n,
        "status": "ok",
        "theorem_s": (
            "Every s ∈ R_n \\ {0,1,∞} arises as the cross-ratio of the ordered "
            "4-tuple (0,1,∞,s) of R_n-rational points of P¹. Hence the Hurwitz "
            "coordinate of any 4-branch cover over R_n may be taken in R_n "
            "whenever the branch locus is R_n-rational."
        ),
        "standard_charts": configs[:4],
        "cosine_4point_samples": cr_samples,
    }


def pure_even_over_Rn_identity() -> dict:
    m, k = sp.symbols("m k")
    alpha = 256 * m**2 - 3125 * k**4 / 256
    beta = k * alpha
    D = sp.expand(256 * alpha**5 + 3125 * beta**4)
    exp = sp.expand((256 * alpha**2 * m) ** 2)
    return {
        "identity": sp.expand(D - exp) == 0,
        "statement": (
            "For every k ∈ R_n (any n) and m ∈ R_n \\ {0} with α(m,k)≠0, "
            "the BJ fibre x⁵+αx+β over R_n is pure-even: disc is a square in R_n."
        ),
    }


def catalogue_k_vs_Rn(n: int) -> list[dict]:
    return [element_in_Rn_Q(k, n) for k in CATALOGUE_K]


def known_geometric_cases() -> dict:
    """
    From programme artefacts: covers over subfields of R_n and their parameters.
    """
    return {
        "rigid_phi_over_Q": {
            "base_field": "Q",
            "branch_cross_ratio": "classical {0,1,∞} — no free s (r=3)",
            "pure_even_k": None,
            "note": "Fibres odd over Q (disc=5·□); not pure-even source",
        },
        "3A4_s_minus_1": {
            "base_field": "Q(√5) = R_5",
            "branch_cross_ratio_s": -1,
            "s_in_R5": True,
            "cover_params": {
                "c": "-√5",
                "p2": -1,
                "r1": "1/5",
                "r2": "-1/5",
            },
            "closed_form_f_s": False,
            "known_BJ_k_from_fibres": None,
            "note": (
                "s=-1 ∈ Q ⊂ R_5 is a cross-ratio over R_5. "
                "Without f_s ∈ R_5(s)[x] → BJ, no proven list of pure-even k "
                "from this cover. Numeric fibres gave even=0 over Q."
            ),
        },
        "arithmetic_envelope_over_Q": {
            "base_field": "Q",
            "pure_even_k": [str(k) for k in CATALOGUE_K],
            "geometric_cover": False,
            "note": "All catalogue k arise arithmetically; not as Nielsen cross-ratios",
        },
        "R5_cosine_k": {
            "k": "2cos(2π/5)=(-1+√5)/2",
            "in_R5": True,
            "pure_even": True,
            "is_cross_ratio_of_01inf_s": (
                "Yes if we set s=k, branch points {0,1,∞,k} — trivial chart"
            ),
            "is_BJ_ratio_from_known_cover": "unknown",
        },
    }


def necessary_conditions() -> list[str]:
    return [
        "NEC1 (arithmetic): k ∈ R_n is necessary and sufficient for the pure-even "
        "formulae to be defined over R_n with α,β ∈ R_n(m) when m ∈ R_n.",
        "NEC2 (branch locus): For a 4-point cover over R_n with R_n-rational branch "
        "locus, the cross-ratio s lies in P¹(R_n). This constrains s, not k.",
        "NEC3 (BJ link): A pure-even k arises geometrically only if some fibre of a "
        "cover over R_n is R_n-birational to a BJ quintic with β/α = k (after "
        "coordinate change over R_n).",
        "NEC4 (descent): If the cover is defined over R_n but the fibre field is a "
        "proper extension, k may lie in that extension, not in R_n.",
        "NEC5 (cosine geometry): If branch points are constrained to cyclotomic "
        "real loci (cosine relations), then s (and possibly k) lie in a thin "
        "subset of R_n — typically multi-angle / Chebyshev values — not all of R_n.",
    ]


def partial_answer() -> dict:
    return {
        "short": (
            "Every k ∈ R_n is pure-even over R_n (arithmetic). "
            "Every s ∈ R_n\\{0,1} is a branch cross-ratio over R_n (trivial). "
            "Which pure-even k arise as BJ ratios of fibres of covers over R_n "
            "is open in general; known constraints and cases are listed below."
        ),
        "classification": {
            "all_of_Rn_as_pure_even_params": {
                "status": "yes",
                "meaning": "k ∈ R_n ⇒ pure-even family over R_n exists",
            },
            "all_of_Rn_as_branch_cross_ratios_s": {
                "status": "yes",
                "meaning": "s ∈ R_n\\{0,1} ⇒ 4-branch chart over R_n",
            },
            "catalogue_Q_k_inside_Rn": {
                "status": "yes",
                "meaning": "all multi-seed k ∈ Q ⊂ R_n",
            },
            "cosine_k_as_geometric_BJ_ratio": {
                "status": "open",
                "meaning": "k=2cos(2πp/n) pure-even over R_n; geometric origin unknown",
            },
            "k_from_3A4_fibres_over_R5": {
                "status": "open",
                "meaning": "s=-1 over R_5 known; map fibre→BJ k not closed-formed",
            },
            "k_equals_s_identification": {
                "status": "not_forced",
                "meaning": (
                    "Identifying pure-even k with the Hurwitz cross-ratio s is a "
                    "coordinate choice, not a theorem. In general k = κ(s) for some "
                    "rational function / algebraic function κ of the cover moduli."
                ),
            },
        },
    }


def main():
    t0 = time.time()
    print("Q1: pure-even k ∈ R_n as cross-ratios of covers over R_n", flush=True)

    pe = pure_even_over_Rn_identity()
    print(f"  pure-even identity: {pe['identity']}", flush=True)

    known = known_geometric_cases()
    nec = necessary_conditions()
    partial = partial_answer()

    by_n = {}
    for n in (5, 7, 11, 15):
        print(f"  n={n}...", flush=True)
        by_n[str(n)] = {
            "degree": Rn_degree(n),
            "catalogue_k": catalogue_k_vs_Rn(n),
            "cosine_candidates": cosine_values_in_Rn(n),
            "cross_ratios": cross_ratios_from_Rn_points(n),
        }

    # Cosine k that accidentally match catalogue rationals (none expected)
    cosine_cat_hits = []
    for n, block in by_n.items():
        for c in block["cosine_candidates"]:
            if c.get("equals_catalogue_k"):
                cosine_cat_hits.append({"n": n, **c})

    elapsed = round(time.time() - t0, 2)

    lines = [
        r"# Which pure-even \(k\in R_n\) arise as cross-ratios of a cover over \(R_n\)?",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        r"**Question (N3 / Q1).** Among pure-even ratio classes \(k=\beta/\alpha\in R_n\),",
        r"which arise from covers defined over the real subfield",
        r"\(R_n=\mathbb{Q}(2\cos 2\pi/n)\)?",
        "",
        "---",
        "",
        r"## 0. Two different quantities (do not conflate)",
        "",
        r"| symbol | meaning | lives in |",
        r"|--------|---------|----------|",
        r"| \(s\) | **Cross-ratio** of four branch points (Hurwitz / \(M_{0,4}\)) | \(\mathbb{P}^1\) |",
        r"| \(k\) | **Pure-even ratio** \(\beta/\alpha\) of a BJ fibre \(x^5+\alpha x+\beta\) | field of coeffs |",
        "",
        r"A cover over \(R_n\) has branch cross-ratio \(s\in\mathbb{P}^1(R_n)\) when the",
        r"branch locus is \(R_n\)-rational. Separately, a BJ model of a fibre may have",
        r"ratio \(k\in R_n\). **In general \(k\neq s\)**; at best \(k=\kappa(s)\) for an",
        r"unknown algebraic map \(\kappa\) attached to the cover type.",
        "",
        "---",
        "",
        r"## 1. Arithmetic answer (complete)",
        "",
        f"**Identity:** `{pe['identity']}`",
        "",
        pe["statement"],
        "",
        r"$$\alpha(m)=256m^2-\frac{3125k^4}{256},\quad"
        r"\beta=k\alpha,\quad"
        r"\operatorname{disc}=(256\alpha^2 m)^2\quad\text{in }R_n(m).$$",
        "",
        r"**Corollary.** *Every* \(k\in R_n\setminus\{0\}\) arises as a pure-even parameter",
        r"over \(R_n\). There is no arithmetic restriction beyond \(k\in R_n\).",
        "",
        "---",
        "",
        r"## 2. Branch cross-ratio \(s\) over \(R_n\) (complete, almost trivial)",
        "",
        r"**Theorem.** For any \(s\in R_n\setminus\{0,1\}\), the ordered 4-tuple",
        r"\((0,1,\infty,s)\) consists of \(R_n\)-rational points of \(\mathbb{P}^1\), and \(s\)",
        r"is their cross-ratio. Hence **every** such \(s\) arises as the branch",
        r"cross-ratio of a 4-point configuration over \(R_n\).",
        "",
        r"So if the question is read as “which \(s\in R_n\) are cross-ratios of covers",
        r"over \(R_n\)?”, the answer is: **all \(s\in R_n\setminus\{0,1,\infty\}\)**",
        r"(subject only to the cover existing for that Nielsen type at that \(s\)).",
        "",
        r"Existence of an \(A_5\) cover of type e.g. \(3A^4\) at a given \(s\in R_n\) is a",
        r"**Hurwitz** question (genus-0 reduced space \(\cong\mathbb{P}^1_s\) over \(\mathbb{Q}\)",
        r"already known; base change to \(R_n\) is free).",
        "",
        "---",
        "",
        r"## 3. Non-trivial reading: pure-even \(k\) from geometric fibres",
        "",
        r"**Intended meaning:** which \(k\in R_n\) appear as \(\beta/\alpha\) for a BJ model",
        r"of a fibre of a cover \(X\to\mathbb{P}^1\) defined over \(R_n\)?",
        "",
        partial["short"],
        "",
        r"### Classification",
        "",
        r"| class | status | meaning |",
        r"|-------|--------|---------|",
    ]
    for key, info in partial["classification"].items():
        lines.append(
            f"| `{key}` | **{info['status']}** | {info['meaning']} |"
        )

    lines += [
        "",
        r"### Necessary conditions",
        "",
    ]
    for item in nec:
        lines.append(f"- {item}")

    lines += [
        "",
        "---",
        "",
        r"## 4. Catalogue \(k\in\mathbb{Q}\) inside every \(R_n\)",
        "",
        r"All multi-seed pure-even ratios from the HQCC catalogue lie in \(\mathbb{Q}\subset R_n\):",
        "",
        f"`{[str(k) for k in CATALOGUE_K]}`",
        "",
        r"They **do** arise as pure-even parameters over \(R_n\), but their known origin is",
        r"**arithmetic** (envelope over \(\mathbb{Q}\)), not as Nielsen-labelled cross-ratios.",
        r"Whether a cover over \(R_n\) specialises to these \(k\) remains open",
        r"(geometric multi-\(k\) problem).",
        "",
        "---",
        "",
        r"## 5. Cosine candidates \(k=2\cos(2\pi p/d)\in R_n\)",
        "",
        r"These are the \(k\) “visibly constrained by cosine relations.”",
        r"They are pure-even over \(R_n\) whenever \(d\mid n\) (so \(k\in R_d\subset R_n\)",
        r"in the cyclotomic plus tower). Matching a catalogue rational is accidental.",
        "",
    ]
    for n in ("5", "7", "11", "15"):
        cands = by_n[n]["cosine_candidates"]
        lines.append(f"### \(n={n}\) (deg \(R_n={by_n[n]['degree']}\))")
        lines.append("")
        lines.append(r"| form | deg | numeric | = catalogue \(k\)? |")
        lines.append(r"|------|----:|--------:|:------------------:|")
        for c in cands[:10]:
            lines.append(
                f"| `{c['form']}` | {c['degree']} | {c['numeric']:.6f} | "
                f"{c['equals_catalogue_k']} |"
            )
        lines.append("")

    if cosine_cat_hits:
        lines.append(f"Cosine–catalogue collisions: `{cosine_cat_hits}`")
    else:
        lines.append(
            r"**No** cosine value \(2\cos(2\pi p/d)\) in the scanned range equals a "
            r"catalogue rational \(k\) (as expected: catalogue \(k\) are rational; "
            r"non-rational cosines are irrational)."
        )

    lines += [
        "",
        "---",
        "",
        r"## 6. Known geometric cases in the programme",
        "",
        r"### Rigid \(\varphi/\mathbb{Q}\) (r=3)",
        "",
        f"- {known['rigid_phi_over_Q']['note']}",
        "",
        r"### \(3A^4\) at \(s=-1\) over \(R_5=\mathbb{Q}(\sqrt5)\)",
        "",
        f"- Branch cross-ratio \(s={known['3A4_s_minus_1']['branch_cross_ratio_s']}\) "
        f"∈ \(R_5\): **{known['3A4_s_minus_1']['s_in_R5']}**",
        f"- Cover params: `{known['3A4_s_minus_1']['cover_params']}`",
        f"- Closed form \(f_s\\to\\mathrm{{BJ}}\): **{known['3A4_s_minus_1']['closed_form_f_s']}**",
        f"- {known['3A4_s_minus_1']['note']}",
        "",
        r"### Arithmetic envelope (not geometric)",
        "",
        f"- Pure-even \(k\): `{known['arithmetic_envelope_over_Q']['pure_even_k']}`",
        f"- {known['arithmetic_envelope_over_Q']['note']}",
        "",
        r"### Cosine \(k\) as chart value",
        "",
        f"- \(k={known['R5_cosine_k']['k']}\) ∈ \(R_5\), pure-even **True**",
        f"- As branch \(s=k\) in chart \(\{{0,1,\\infty,s\}}\): trivial yes",
        f"- As BJ ratio from a known cover: **{known['R5_cosine_k']['is_BJ_ratio_from_known_cover']}**",
        "",
        "---",
        "",
        r"## 7. Partial answer (lock)",
        "",
        r"| reading of the question | answer |",
        r"|-------------------------|--------|",
        r"| Which \(k\in R_n\) admit pure-even families over \(R_n\)? | **All** \(k\in R_n\setminus\{0\}\) |",
        r"| Which \(s\in R_n\) are branch cross-ratios over \(R_n\)? | **All** \(s\in R_n\setminus\{0,1\}\) |",
        r"| Which catalogue \(k\in\mathbb{Q}\) lie in \(R_n\)? | **All** of them (\(\mathbb{Q}\subset R_n\)) |",
        r"| Which \(k\in R_n\) are BJ ratios of fibres of a cover over \(R_n\)? | **Open** — no closed \(s\mapsto k\) |",
        r"| Are cosine values distinguished? | **Yes as geometric candidates**; not forced by evenness |",
        r"| Is \(k=s\)? | **Not in general** |",
        "",
        r"### What would finish the non-trivial reading",
        "",
        r"1. Closed form \(f_s\in R_n(s)[x]\) (or over \(\mathbb{Q}(s)\)) for a Nielsen type.",
        r"2. BJ reduction: after Möbius in \(x\), read \(k(s)=\beta(s)/\alpha(s)\in R_n(s)\).",
        r"3. Image of \(k: H^{\mathrm{rd}}(R_n)\to\mathbb{P}^1\) — that image **is** the answer.",
        "",
        r"**Until then:** arithmetic multi-\(k\) over \(\mathbb{Q}\) remains the citable source of",
        r"explicit pure-even \(k\); over \(R_n\), pure-even is free in \(k\), while geometric",
        r"origin of specific \(k\) is constrained by Hurwitz data not yet converted to BJ.",
        "",
        r"---",
        "",
        r"## 8. Direct answers to the one-line question",
        "",
        r"> Which pure-even \(k\in R_n\) arise as cross-ratios of a cover over \(R_n\)?",
        "",
        r"**If “cross-ratio” means the Hurwitz parameter \(s\):**  ",
        r"every \(s\in R_n\setminus\{0,1\}\) (and the cover type must exist at that \(s\)).",
        "",
        r"**If “cross-ratio” is used loosely for pure-even \(k\) coming from such a cover:**  ",
        r"unknown list; must equal the image of the moduli map \(s\mapsto k(s)\) once a",
        r"BJ model exists. Currently known rigorously:",
        "",
        r"- all \(k\in R_n\) work **arithmetically** as pure-even parameters;",
        r"- all catalogue \(k\in\mathbb{Q}\) embed in every \(R_n\);",
        r"- \(s=-1\in R_5\) is a genuine geometric cross-ratio of a \(3A^4\) cover over \(R_5\);",
        r"- no theorem yet names a non-rational \(k\in R_n\setminus\mathbb{Q}\) as a BJ ratio",
        r"  of an \(R_n\)-cover fibre.",
        "",
        r"```bash",
        r"python pure_even_k_cross_ratio.py",
        r"```",
        "",
        r"_Generated by pure_even_k_cross_ratio.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "partial_answer": partial,
        "pure_even_identity": pe,
        "necessary_conditions": nec,
        "known_geometric_cases": known,
        "by_n": by_n,
        "cosine_catalogue_hits": cosine_cat_hits,
        "catalogue_k": [str(k) for k in CATALOGUE_K],
    }
    md = "\n".join(lines)
    write_md(ROOT / "PURE_EVEN_K_CROSS_RATIO.md", md)
    write_json(ROOT / "PURE_EVEN_K_CROSS_RATIO.json", payload)
    write_md(OUT / "PURE_EVEN_K_CROSS_RATIO.md", md)
    write_json(OUT / "PURE_EVEN_K_CROSS_RATIO.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "PURE_EVEN_K_CROSS_RATIO.md", md)
    except Exception:
        pass

    print(partial["short"], flush=True)
    print(f"Wrote PURE_EVEN_K_CROSS_RATIO.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
