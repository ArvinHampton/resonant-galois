"""
Q2: Does Gal(R_n/Q) act so that multi-k paths are unions of cosine orbits?

Objects:
  R_n = Q(ξ_n), ξ_n = 2 cos(2π/n)
  Gal(R_n/Q) ≅ (Z/nZ)*/{±1},  σ_a: ξ_n ↦ 2 cos(2π a/n)

  multi-k path: rational curve (m(u), k(u)) in pure-even envelope
                joining distinct ratio classes (programme: linear in u)

  cosine orbit: Gal·(2 cos(2π p/n)) = {2 cos(2π a p/n) : a ∈ (Z/n)*/±}
                — finite discrete subset of R_n ⊂ R̄

Answer: No in the arithmetic (catalogue) setting; only under strong geometric
hypotheses could a Galois-stable set of k-values be a union of cosine orbits.

Output: GAL_RN_MULTI_K_ORBITS.md / .json
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402

# Programme multi-k paths (from pure-even envelope)
PATHS = [
    {
        "id": "flag_classical",
        "k1": Fraction(-8, 5),
        "k2": Fraction(4, 5),
        "m1": Fraction(5, 16),
        "m2": Fraction(5, 16),
    },
    {
        "id": "flag_lsw",
        "k1": Fraction(-8, 5),
        "k2": Fraction(-4),
        "m1": Fraction(5, 16),
        "m2": Fraction(55, 16),
    },
    {
        "id": "classical_lsw",
        "k1": Fraction(4, 5),
        "k2": Fraction(-4),
        "m1": Fraction(5, 16),
        "m2": Fraction(55, 16),
    },
]

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


def gal_Rn_structure(n: int) -> dict:
    """Gal(R_n/Q) ≅ (Z/n)*/{±1}."""
    units = [a for a in range(1, n) if sp.gcd(a, n) == 1]
    # identify a ~ -a mod n
    classes = []
    seen = set()
    for a in units:
        key = frozenset({a % n, (-a) % n})
        if key in seen:
            continue
        seen.add(key)
        classes.append(sorted(key)[0])  # representative
    phi = sp.totient(n)
    deg = int(phi // 2) if n >= 3 else 1
    return {
        "n": n,
        "phi_n": int(phi),
        "degree_Rn": deg,
        "gal_order": deg,
        "isomorphism": "(Z/nZ)*/{±1}",
        "action": "σ_a: 2cos(2π/n) ↦ 2cos(2π a/n)",
        "representatives_a": classes[:20],
        "n_classes": len(classes),
    }


def cosine_orbit(n: int, p: int = 1) -> dict:
    """Orbit of 2 cos(2π p/n) under Gal(R_n/Q)."""
    units = [a for a in range(1, n) if sp.gcd(a, n) == 1]
    vals = []
    seen_num = set()
    for a in units:
        # cos is even and 2π-periodic; identify a p and -a p
        ap = (a * p) % n
        ap = min(ap, (-ap) % n)
        expr = 2 * sp.cos(2 * sp.pi * ap / n)
        num = round(float(sp.N(expr, 25)), 12)
        if num in seen_num:
            continue
        seen_num.add(num)
        try:
            mp = str(sp.minpoly(expr, sp.symbols("x")))
            deg = int(sp.degree(sp.minpoly(expr, sp.symbols("x"))))
        except Exception:
            mp, deg = None, None
        vals.append(
            {
                "form": f"2cos(2π·{ap}/{n})",
                "numeric": float(sp.N(expr, 20)),
                "minpoly": mp,
                "degree": deg,
            }
        )
    return {
        "n": n,
        "seed": f"2cos(2π·{p}/{n})",
        "orbit_size": len(vals),
        "orbit": vals,
        "finite_discrete": True,
    }


def rational_fixed_by_gal() -> dict:
    """Every k ∈ Q is fixed by Gal(R_n/Q) (Q is the base)."""
    return {
        "theorem": (
            "For every n, Gal(R_n/Q) acts trivially on Q. "
            "Hence every catalogue / rational pure-even ratio k ∈ Q is a "
            "Gal(R_n/Q)-fixed point: its orbit is {k}."
        ),
        "corollary": (
            "A multi-k path whose image lies in Q (e.g. linear k(u)∈Q(u) with "
            "rational endpoints) is a union of Gal-orbits, but each orbit is a "
            "singleton rational — not a cosine orbit of size >1."
        ),
    }


def path_analysis() -> dict:
    """Analyse programme multi-k paths under Gal and cosine."""
    u = sp.symbols("u")
    rows = []
    for path in PATHS:
        k1, k2 = path["k1"], path["k2"]
        # linear path k(u) = k1 + u(k2-k1)
        ku = k1 + u * (k2 - k1)
        # sample rational points
        samples = []
        for j in range(0, 7):
            uv = Fraction(j, 6)
            kv = k1 + uv * (k2 - k1)
            samples.append({"u": str(uv), "k": str(kv), "k_in_Q": True})

        # intermediate k that are not endpoints
        intermediate = []
        for j in range(1, 6):
            uv = Fraction(j, 6)
            kv = k1 + uv * (k2 - k1)
            intermediate.append(kv)

        rows.append(
            {
                "id": path["id"],
                "k_path": f"{k1} + u·({k2}-{k1})",
                "k_path_in_Q_u": True,
                "endpoints": [str(k1), str(k2)],
                "samples": samples,
                "all_sample_k_in_Q": True,
                "gal_orbits_of_samples": (
                    "each k∈Q has orbit {k} — path = disjoint union of singletons"
                ),
                "is_cosine_orbit": False,
                "is_union_of_cosine_orbits": False,
                "reason": (
                    "Image of the path over Q is an infinite set of rationals "
                    "(or a line segment in P¹(Q)), while every cosine orbit is "
                    "finite of size ≤ φ(n)/2. An infinite set of rationals cannot "
                    "be a finite union of finite cosine orbits unless the path "
                    "image is finite — which a non-constant rational path is not "
                    "when evaluated on infinitely many u."
                ),
            }
        )
    return {"paths": rows}


def can_path_be_union_of_cosine_orbits() -> dict:
    """
    Mathematical obstruction and the only scenarios where a weakened statement holds.
    """
    return {
        "obstruction_infinite_vs_finite": (
            "A non-constant multi-k path k: P¹ ⇢ P¹ over Q has infinitely many "
            "distinct values k(u) for u ∈ Q (or u ∈ R̄). Each cosine orbit is finite. "
            "A finite union of cosine orbits is finite. Therefore a non-constant "
            "path image cannot equal a union of cosine orbits as sets of field elements."
        ),
        "obstruction_rational_vs_cosine": (
            "Catalogue multi-k paths take values in Q. Non-rational cosine values "
            "2cos(2π p/n) ∉ Q for n>2 (except degenerate cases). "
            "So path values are not cosine elements, and cosine orbits of size >1 "
            "consist of irrationals not hit by Q-paths."
        ),
        "gal_invariance_weaker": (
            "Weaker true statement: if a set S ⊂ P¹ of k-values is defined over Q "
            "(stable under Gal(Q̄/Q), hence under Gal(R_n/Q) after embedding), "
            "then S is a union of Gal(Q̄/Q)-orbits. For S ⊂ Q this is a union of "
            "singletons. This is Galois descent, not a cosine constraint."
        ),
        "geometric_scenario_where_cosine_appears": (
            "If the geometric multi-k locus is a finite Gal(R_n/Q)-stable set of "
            "special fibres (not a positive-dimensional path over Q), and those k "
            "are forced into the cyclotomic real locus by branch geometry, then "
            "that finite set can be a union of cosine orbits. That is a discrete "
            "specialisation set, not an arithmetic multi-k path in the envelope."
        ),
        "action_on_paths": (
            "Gal(R_n/Q) acts on R_n-points of the envelope. For a path defined over Q, "
            "σ·(m(u),k(u)) = (m(u),k(u)) when m,k ∈ Q(u). The path is fixed as a "
            "scheme over Q. Gal does not permute distinct rational k along the path "
            "into each other (each is fixed)."
        ),
    }


def cosine_vs_catalogue_numeric(n: int) -> dict:
    """Check whether any catalogue k equals any cosine orbit element for this n."""
    orbit = cosine_orbit(n, 1)
    hits = []
    for k in CATALOGUE_K:
        kf = float(k)
        for o in orbit["orbit"]:
            if abs(o["numeric"] - kf) < 1e-9:
                hits.append({"k": str(k), "cosine": o["form"], "n": n})
    # also check all p
    for p in range(1, n):
        if sp.gcd(p, n) != 1:
            continue
        if p > n // 2:
            continue
        orb = cosine_orbit(n, p)
        for k in CATALOGUE_K:
            kf = float(k)
            for o in orb["orbit"]:
                if abs(o["numeric"] - kf) < 1e-9:
                    hits.append({"k": str(k), "cosine": o["form"], "n": n})
    # unique
    uniq = []
    seen = set()
    for h in hits:
        key = (h["k"], h["cosine"], h["n"])
        if key not in seen:
            seen.add(key)
            uniq.append(h)
    return {
        "n": n,
        "orbit_size_seed_p1": orbit["orbit_size"],
        "catalogue_cosine_hits": uniq,
        "any_hit": len(uniq) > 0,
    }


def worked_example_n5() -> dict:
    """R_5 = Q(√5), Gal = C2, cosine orbit of 2cos(2π/5)."""
    # 2cos(2π/5) = (-1+√5)/2,  2cos(4π/5) = (-1-√5)/2
    # Gal swaps √5 ↦ -√5, so swaps these two
    k_plus = (-1 + sp.sqrt(5)) / 2
    k_minus = (-1 - sp.sqrt(5)) / 2
    orbit = [k_plus, k_minus]
    # multi-k path flag-classical: k = -8/5 + u*(4/5 - (-8/5)) = -8/5 + u*(12/5)
    # never equals ±(√5 related) as rational function identity
    u = sp.symbols("u")
    k_path = sp.Rational(-8, 5) + u * (sp.Rational(4, 5) - sp.Rational(-8, 5))
    # Does k_path(u) = k_plus for some u in R_5?
    # -8/5 + u*12/5 = (-1+√5)/2
    # u*12/5 = (-1+√5)/2 + 8/5 = (-5 + 5√5 + 16)/10 = (11 + 5√5)/10
    # u = (5/12)*((11+5√5)/10) = (11+5√5)/24 ∈ R_5
    u_hit = sp.simplify(
        (k_plus - sp.Rational(-8, 5))
        / (sp.Rational(4, 5) - sp.Rational(-8, 5))
    )
    # conjugate
    u_hit_c = sp.simplify(
        (k_minus - sp.Rational(-8, 5))
        / (sp.Rational(4, 5) - sp.Rational(-8, 5))
    )
    return {
        "R5": "Q(√5)",
        "gal": "C2, √5 ↦ -√5",
        "cosine_orbit": {
            "2cos(2π/5)": str(k_plus),
            "2cos(4π/5)": str(k_minus),
            "orbit_size": 2,
        },
        "flag_classical_path": "k(u) = -8/5 + (12/5)u",
        "u_where_path_hits_cosine": {
            "u_for_2cos(2π/5)": str(u_hit),
            "u_for_2cos(4π/5)": str(u_hit_c),
            "are_gal_conjugates": sp.simplify(u_hit.subs(sp.sqrt(5), -sp.sqrt(5)) - u_hit_c)
            == 0
            or str(sp.simplify(u_hit.subs(sp.sqrt(5), -sp.sqrt(5)))) == str(u_hit_c),
        },
        "interpretation": (
            "The arithmetic path (over Q) can pass through a cosine value at a "
            "non-rational u ∈ R_5. The Gal-conjugate cosine is hit at the conjugate u. "
            "So: the path is NOT a union of cosine orbits; rather, the pair "
            "{u, σu} of parameters maps to a cosine orbit {k, σk}. "
            "That is Gal acting on the parameter of the path, not the path being "
            "made of cosine orbits."
        ),
    }


def main():
    t0 = time.time()
    print("Q2: Gal(R_n/Q) vs multi-k paths vs cosine orbits", flush=True)

    gal_structs = {str(n): gal_Rn_structure(n) for n in (5, 7, 11, 15)}
    fixed = rational_fixed_by_gal()
    paths = path_analysis()
    obst = can_path_be_union_of_cosine_orbits()
    n5 = worked_example_n5()

    cosine_checks = {}
    for n in (5, 7, 11, 15):
        print(f"  cosine vs catalogue n={n}", flush=True)
        cosine_checks[str(n)] = cosine_vs_catalogue_numeric(n)

    orbits_sample = {str(n): cosine_orbit(n, 1) for n in (5, 7, 11)}

    any_cat_cosine = any(c["any_hit"] for c in cosine_checks.values())

    # Verdict
    answer_no = True
    verdict = (
        "NO — multi-k paths (arithmetic envelope) are not unions of cosine orbits. "
        "Gal(R_n/Q) fixes Q-points (catalogue k are singletons). Cosine orbits are "
        "finite/discrete; non-constant paths have infinite rational image. "
        "Weaker true facts: (i) Q-paths are Gal-invariant as schemes; "
        "(ii) finite geometric specialisation sets over R_n can be unions of "
        "cosine orbits; (iii) Gal acts on path parameters so conjugate u hit "
        "conjugate cosine k when the path crosses the cosine locus."
    )

    elapsed = round(time.time() - t0, 2)
    print(verdict, flush=True)

    lines = [
        r"# Does \(\mathrm{Gal}(R_n/\mathbb{Q})\) act so multi-\(k\) paths are unions of cosine orbits?",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Answer: NO** (for the programme’s arithmetic multi-\(k\) paths).",
        "",
        verdict,
        "",
        "---",
        "",
        r"## 0. Definitions",
        "",
        r"| object | definition |",
        r"|--------|------------|",
        r"| \(R_n\) | \(\mathbb{Q}(2\cos 2\pi/n)=\mathbb{Q}(\zeta_n)^+\) |",
        r"| \(\mathrm{Gal}(R_n/\mathbb{Q})\) | \(\cong(\mathbb{Z}/n\mathbb{Z})^\times/\{\pm1\}\), order \(\varphi(n)/2\) |",
        r"| action | \(\sigma_a:\ 2\cos\frac{2\pi}{n}\mapsto 2\cos\frac{2\pi a}{n}\) |",
        r"| multi-\(k\) path | rational curve \((m(u),k(u))\) in the pure-even envelope joining distinct ratio classes |",
        r"| cosine orbit | \(\mathrm{Gal}\cdot\bigl(2\cos\frac{2\pi p}{n}\bigr)\) — **finite** subset of \(R_n\) |",
        "",
        r"### Group orders (proxies)",
        "",
        r"| \(n\) | \([R_n:\mathbb{Q}]\) | \(\#\mathrm{Gal}\) |",
        r"|------|------------------:|---------------:|",
    ]
    for n, g in gal_structs.items():
        lines.append(f"| {n} | {g['degree_Rn']} | {g['gal_order']} |")

    lines += [
        "",
        "---",
        "",
        r"## 1. Galois action on rational \(k\)",
        "",
        fixed["theorem"],
        "",
        fixed["corollary"],
        "",
        r"**Catalogue \(k\)** (all in \(\mathbb{Q}\)):",
        f"`{[str(k) for k in CATALOGUE_K]}`",
        "",
        r"Each has \(\mathrm{Gal}(R_n/\mathbb{Q})\)-orbit of size **1**. They are not",
        r"non-trivial cosine orbits.",
        "",
        "---",
        "",
        r"## 2. Cosine orbits are finite and discrete",
        "",
        r"Example orbits of \(2\cos(2\pi/n)\):",
        "",
    ]
    for n, orb in orbits_sample.items():
        lines.append(
            f"- **n={n}**: orbit size **{orb['orbit_size']}** — "
            f"`{[o['form'] for o in orb['orbit'][:6]]}`"
        )
    lines += [
        "",
        r"A non-constant path \(k(u)\in\mathbb{Q}(u)\) takes **infinitely many** distinct",
        r"values at rational \(u\). A finite union of cosine orbits is **finite**.",
        r"Hence:",
        "",
        r"> **Obstruction.** A non-constant multi-\(k\) path cannot equal a union of",
        r"> cosine orbits as a set of field elements.",
        "",
        f"Details: {obst['obstruction_infinite_vs_finite']}",
        "",
        "---",
        "",
        r"## 3. Programme multi-\(k\) paths",
        "",
        r"| path | \(k(u)\) | endpoints | cosine orbit? | union of cosine orbits? |",
        r"|------|----------|-----------|:-------------:|:-----------------------:|",
    ]
    for p in paths["paths"]:
        lines.append(
            f"| {p['id']} | `{p['k_path']}` | {p['endpoints']} | "
            f"**{p['is_cosine_orbit']}** | **{p['is_union_of_cosine_orbits']}** |"
        )

    lines += [
        "",
        r"### Why not",
        "",
    ]
    for p in paths["paths"][:1]:
        lines.append(p["reason"])
        lines.append("")

    lines += [
        obst["obstruction_rational_vs_cosine"],
        "",
        f"**Catalogue vs cosine numeric collisions (n=5,7,11,15):** "
        f"{'YES — unexpected' if any_cat_cosine else '**none**'}",
        "",
    ]
    for n, c in cosine_checks.items():
        lines.append(
            f"- n={n}: hits=`{c['catalogue_cosine_hits']}`"
        )

    lines += [
        "",
        "---",
        "",
        r"## 4. What \(\mathrm{Gal}(R_n/\mathbb{Q})\) *does* do to paths",
        "",
        obst["action_on_paths"],
        "",
        obst["gal_invariance_weaker"],
        "",
        r"### Worked example: \(R_5=\mathbb{Q}(\sqrt5)\), path flag↔classical",
        "",
        f"- Cosine orbit: `{n5['cosine_orbit']}`",
        f"- Path: `{n5['flag_classical_path']}`",
        f"- Hits cosine at u = `{n5['u_where_path_hits_cosine']['u_for_2cos(2π/5)']}`",
        f"- Hits conjugate cosine at u = `{n5['u_where_path_hits_cosine']['u_for_2cos(4π/5)']}`",
        "",
        n5["interpretation"],
        "",
        r"**Picture:**",
        r"$$\mathrm{Gal}\curvearrowright u\quad\Longrightarrow\quad"
        r"k(u)\mapsto k(\sigma u)=\sigma\bigl(k(u)\bigr)"
        r"\quad\text{when }k(u)\in R_n,$$",
        r"for the linear path with \(\mathbb{Q}\)-coefficients. The path is a",
        r"\(\mathbb{Q}\)-curve; Gal permutes the **parameters** where it meets a",
        r"cosine locus, not the path into a cosine orbit.",
        "",
        "---",
        "",
        r"## 5. When could a *weakened* cosine-orbit statement hold?",
        "",
        obst["geometric_scenario_where_cosine_appears"],
        "",
        r"| scenario | multi-\(k\) path = ∪ cosine orbits? |",
        r"|----------|:-----------------------------------:|",
        r"| Arithmetic envelope paths over \(\mathbb{Q}\) | **No** |",
        r"| Finite Gal-stable set of geometric specialisations in cosine locus | **Possible** (discrete, not a path) |",
        r"| Path over \(R_n\) with image inside one cosine orbit | **No** (orbit finite, path infinite unless constant) |",
        r"| Path meets cosine locus at a Gal-orbit of parameters | **Yes as incidence**, not as path=orbit |",
        "",
        "---",
        "",
        r"## 6. Locked answer",
        "",
        r"> **Does \(\mathrm{Gal}(R_n/\mathbb{Q})\) act so that multi-\(k\) paths are",
        r"> unions of cosine orbits?**",
        "",
        r"**No.**",
        "",
        r"1. Cosine orbits are finite; non-constant multi-\(k\) paths are not.",
        r"2. Programme paths take values in \(\mathbb{Q}\), fixed pointwise by",
        r"   \(\mathrm{Gal}(R_n/\mathbb{Q})\); each rational \(k\) is a singleton orbit,",
        r"   not a non-trivial cosine orbit.",
        r"3. Catalogue \(k\) do not coincide with \(2\cos(2\pi p/n)\) for proxy \(n\).",
        r"4. True related facts:",
        r"   - paths over \(\mathbb{Q}\) are Gal-invariant as schemes;",
        r"   - Gal acts on parameters \(u\in R_n\) along the path;",
        r"   - finite geometric \(k\)-sets cut out by cosine branch constraints",
        r"     *can* be unions of cosine orbits — that is a different object from",
        r"     arithmetic multi-\(k\) paths.",
        "",
        r"```bash",
        r"python gal_rn_multi_k_orbits.py",
        r"```",
        "",
        r"_Generated by gal_rn_multi_k_orbits.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "answer": "NO",
        "verdict": verdict,
        "gal_structures": gal_structs,
        "rational_fixed": fixed,
        "paths": paths,
        "obstructions": obst,
        "worked_example_R5": {
            "cosine_orbit": n5["cosine_orbit"],
            "path": n5["flag_classical_path"],
            "u_hits": n5["u_where_path_hits_cosine"],
            "interpretation": n5["interpretation"],
        },
        "cosine_vs_catalogue": cosine_checks,
        "orbits_sample": orbits_sample,
    }
    # JSON-serialize fractions in n5
    payload["worked_example_R5"]["u_hits"] = {
        k: str(v) for k, v in n5["u_where_path_hits_cosine"].items()
    }

    md = "\n".join(lines)
    write_md(ROOT / "GAL_RN_MULTI_K_ORBITS.md", md)
    write_json(ROOT / "GAL_RN_MULTI_K_ORBITS.json", payload)
    write_md(OUT / "GAL_RN_MULTI_K_ORBITS.md", md)
    write_json(OUT / "GAL_RN_MULTI_K_ORBITS.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "GAL_RN_MULTI_K_ORBITS.md", md)
    except Exception:
        pass

    print(f"Wrote GAL_RN_MULTI_K_ORBITS.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
