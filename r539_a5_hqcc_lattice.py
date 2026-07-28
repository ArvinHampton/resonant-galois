"""
Q4: Does there exist f ∈ R_539[x] with monodromy A5 reducing at split primes
    to the HQCC Z-lattice?

Readings:
  (T) Trivial: f ∈ Z[x] ⊂ R_539[x] an HQCC A5 seed — yes immediately.
  (N) Non-trivial: f ∈ R_539[x] not over Q, Gal(f/R_539)=A5, and at primes p
      that split completely in R_539, reductions of f recover HQCC lattice seeds
      (or their mod-p reductions / Frobenius fingerprint).

R_539 = Q(2 cos(2π/539)), [R_539:Q] = φ(539)/2 = 210. Direct Gal over R_539
is not computational; arguments are field-theoretic + proxy checks on divisors.

Output: R539_A5_HQCC_LATTICE.md / .json
"""
from __future__ import annotations

import sys
import time
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import (  # noqa: E402
    MODEL_CORE,
    OUT,
    RESULTS,
    classify_poly,
    is_square,
    write_json,
    write_md,
    x,
)
from lib.lemmas import disc_bj_int  # noqa: E402

N_PERIOD = 539  # 7^2 * 11

# HQCC Z-lattice multi-seed A5 catalogue (sample)
HQCC_SEEDS = [
    (-55, 88, "flagship"),
    (145, -232, "flagship_m"),
    (320, -512, "flagship_m2"),
    (20, 16, "classical"),
    (95, 76, "classical_m"),
    (-100, 400, "lsw"),
    (124, -496, "lsw_m"),
    (-180, 432, "s12"),
]


def r539_field_data() -> dict:
    fac = sp.factorint(N_PERIOD)
    phi = int(sp.totient(N_PERIOD))
    deg = phi // 2
    return {
        "N": N_PERIOD,
        "factorisation": {str(k): int(v) for k, v in fac.items()},
        "phi_N": phi,
        "degree_R_N": deg,
        "xi": "2*cos(2π/539)",
        "gal_order": deg,
        "computable": False,
        "proxy_divisors": [d for d in sp.divisors(N_PERIOD) if 3 <= int(d) <= 77],
    }


def xi_minpoly(n: int):
    deg = int(sp.totient(n) // 2) if n >= 3 else 1
    if deg > 10:
        return None
    return sp.Poly(sp.minpoly(2 * sp.cos(2 * sp.pi / n), x), x, domain=sp.ZZ)


def split_primes(n: int, max_p: int = 120) -> list[int]:
    mp = xi_minpoly(n)
    out = []
    if mp is None:
        for p in sp.primerange(3, max_p):
            p = int(p)
            if n % p == 0:
                continue
            if p % n in (1, n - 1):
                out.append(p)
        return out
    for p in sp.primerange(3, max_p):
        p = int(p)
        if n % p == 0:
            continue
        try:
            fac = sp.factor_list(mp.as_expr(), modulus=p)
            degs = sorted(int(sp.degree(f)) for f, m in fac[1] for _ in range(int(m)))
            if degs == [1] * mp.degree():
                out.append(p)
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------------
# (T) Trivial reading
# ---------------------------------------------------------------------------
def trivial_yes() -> dict:
    rows = []
    for a, b, name in HQCC_SEEDS:
        d = disc_bj_int(a, b)
        sq = d > 0 and is_square(d)
        rec = classify_poly(x**5 + a * x + b, do_galois=True)
        rows.append(
            {
                "name": name,
                "poly": f"x^5+({a})x+({b})",
                "a": a,
                "b": b,
                "disc_square": sq,
                "status": rec.get("status"),
                "in_R539": True,  # Z ⊂ R_539
                "monodromy_A5": (rec.get("status") or "").startswith("HIT_A5"),
            }
        )
    a5 = [r for r in rows if r["monodromy_A5"]]
    return {
        "reading": "trivial",
        "statement": (
            "Any HQCC A5 seed f ∈ Z[x] is automatically an element of R_539[x]. "
            "If the splitting field K is linearly disjoint from R_539 over Q "
            "(generic for these fields), Gal(f/R_539) ≅ Gal(f/Q) = A5. "
            "At every prime p (including those that split in R_539), f mod p is "
            "exactly the reduction of the HQCC Z-lattice polynomial."
        ),
        "answer": "YES",
        "examples": a5,
        "n_A5_examples": len(a5),
        "caveat_disjointness": (
            "If K ∩ R_539 ≠ Q, then Gal(f/R_539) ≅ Gal(K/K∩R_539) may be proper "
            "subgroup of A5. For flagship-type A5 fields this is expected to fail "
            "only for special n dividing related conductors; not checked exhaustively "
            "for n=539 (deg 210)."
        ),
    }


def compositum_heuristic_flagship() -> dict:
    """
    Heuristic: does the flagship A5 field meet R_p for small p|539?
    If disc(f) is not divisible by primes that must ramify in R_n in a way that
    forces intersection — weak check via ramified primes.
    """
    a, b = -55, 88
    d = disc_bj_int(a, b)
    # square-free kernel of disc (rough)
    sf = abs(d)
    primes_disc = []
    for p in sp.primerange(2, 500):
        p = int(p)
        if sf % p == 0:
            primes_disc.append(p)
            while sf % p == 0:
                sf //= p
        if sf == 1:
            break
    # R_539 ramifies only at primes dividing 539 = 7,11
    ram_R = [7, 11]
    meet_risk = [p for p in ram_R if p in primes_disc]
    return {
        "flagship_disc": d,
        "primes_dividing_disc_sample": primes_disc[:20],
        "R539_ramified_at": ram_R,
        "shared_ramification_7_11": meet_risk,
        "note": (
            "Shared ramification at 7 or 11 does not force K ∩ R_539 ≠ Q, but is "
            "the first place to look for non-disjointness. Full check needs the "
            "A5 field's resolvent / number field database — not done here."
        ),
    }


# ---------------------------------------------------------------------------
# (N) Non-trivial reading
# ---------------------------------------------------------------------------
def nontrivial_reading() -> dict:
    return {
        "reading": "nontrivial",
        "statement": (
            "Does there exist f ∈ R_539[x] \\ Q[x] with Gal(f/R_539)=A5 such that "
            "for a positive-density set of primes p that split completely in R_539, "
            "the reduction f mod P (P|p) is F_p-isomorphic to the reduction of some "
            "HQCC Z-lattice seed (same factorisation type / same poly up to F_p^× "
            "scaling of variable), or recovers lattice coefficients via traces?"
        ),
        "status": "OPEN — no construction; no obstruction ruling out all lifts",
        "why_hard": [
            "deg R_539 = 210 blocks direct symbolic Gal and minpoly work",
            "no closed geometric f_s over R_539 producing BJ with lattice specialisations",
            "reduction of non-rational coeffs requires integral model of O_{R_539}",
            "matching infinitely many split reductions to a fixed Z-seed is a strong "
            "rigidity condition (likely forces f over Q by approximation / Krasner)",
        ],
        "plausible_constructions": [
            {
                "name": "constant family over R_539",
                "idea": "f = HQCC seed ∈ Z[x] ⊂ R_539[x]",
                "nontrivial": False,
                "works": True,
            },
            {
                "name": "Galois conjugate twist",
                "idea": (
                    "Take α,β ∈ R_539 on a pure-even k-slice with k∈R_539, choose so "
                    "that under all embeddings R_539→R the specialised real polys "
                    "relate to lattice; reduce at split p."
                ),
                "nontrivial": True,
                "works": "unknown",
            },
            {
                "name": "norm / Weil restriction of A5 cover over R_539",
                "idea": (
                    "Cover over R_539 with monodromy A5; fibre f ∈ R_539[x]; "
                    "require specialisations at split primes match lattice Frobenius."
                ),
                "nontrivial": True,
                "works": "open (geometric multi-k)",
            },
            {
                "name": "interpolating poly with coeffs in Z[ξ_539]",
                "idea": (
                    "α = α0 + α1 ξ + … with αi ∈ Z lattice-ish; force α≡-55, β≡88 "
                    "mod all split P in a density-zero or finite set — cannot match "
                    "a fixed seed at infinitely many p unless α,β ∈ Z (approx. argument)."
                ),
                "nontrivial": True,
                "works": (
                    "Likely ONLY if f over Q: if f ≡ f0 mod P for infinitely many "
                    "split P and f0 ∈ Z[x] fixed, then f=f0 by congruence."
                ),
            },
        ],
    }


def reduction_rigidity_lemma() -> dict:
    """
    Sketch: if f ∈ O[x] (O=integers of R_539) and f ≡ f0 mod P for infinitely
    many primes P of O with f0 ∈ Z[x] fixed monic of same degree, then f = f0.
    """
    return {
        "lemma_sketch": (
            "Let O be the ring of integers of R_539 (or Z[ξ] order). Let f,f0 ∈ O[x] "
            "be monic of degree d, f0 ∈ Z[x]. Suppose f ≡ f0 (mod P) for infinitely "
            "many prime ideals P of O. Then each coefficient a_i - a_i^{(0)} lies in "
            "infinitely many P, hence is 0. Thus f = f0 ∈ Z[x]."
        ),
        "consequence": (
            "One cannot have a genuinely non-rational f ∈ R_539[x]\\Q[x] whose "
            "reduction equals a *fixed* HQCC seed at infinitely many primes of R_539. "
            "Non-trivial lifts can only match lattice seeds at finitely many split "
            "primes, or match a *varying* family of lattice seeds, or match only "
            "Frobenius cycle types (not the polynomial itself)."
        ),
        "weaker_nontrivial_goals": [
            "Match Frob cycle-type statistics of HQCC A5 seeds at split primes",
            "Match finitely many specified seeds at finitely many split P",
            "Match lattice under one fixed embedding R_539→Q_p or →R, not all reductions",
            "Have traces/norms of coeffs land in the HQCC Z-lattice",
        ],
    }


# ---------------------------------------------------------------------------
# Proxy: split reduction of flagship equals "lattice reduction"
# ---------------------------------------------------------------------------
def proxy_split_reduction_fingerprint(n: int = 7) -> dict:
    """
    For f = flagship ∈ Z[x], at p split in R_n (n|539), factorisation types —
    this IS the HQCC lattice reducing at split primes. Confirms trivial reading
    on proxies for primes dividing 539.
    """
    a, b = -55, 88
    pol = sp.Poly(x**5 + a * x + b, x, domain=sp.ZZ)
    d = disc_bj_int(a, b)
    split = split_primes(n, max_p=150)
    types = Counter()
    used = 0
    for p in split:
        if d % p == 0:
            continue
        try:
            facs = sp.factor_list(pol.as_expr(), modulus=int(p))
            degs = []
            for f, m in facs[1]:
                degs.extend([int(sp.degree(f))] * int(m))
            types[tuple(sorted(degs))] += 1
            used += 1
        except Exception:
            continue
    return {
        "proxy_n": n,
        "seed": "flagship",
        "n_split": len(split),
        "primes_used": used,
        "frob_types": {str(k): v for k, v in types.most_common()},
        "interpretation": (
            "Flagship reductions at split primes of R_n are exactly HQCC lattice "
            "reductions — trivial reading verified on proxy n|539."
        ),
    }


def lattice_coefficient_list() -> dict:
    return {
        "model_core": MODEL_CORE,
        "hqcc_seed_coeffs_sample": [(a, b, n) for a, b, n in HQCC_SEEDS],
        "period_N": N_PERIOD,
        "note": (
            "HQCC Z-lattice = integer combinations of model numbers "
            "{3,9,27,61,80,243,539,…} appearing in seeds; flagship 88=61+27."
        ),
    }


def answers_table() -> list[dict]:
    return [
        {
            "reading": "(T) Trivial — f HQCC seed over Z ⊂ R_539",
            "exists": "YES",
            "status": "proved by example",
        },
        {
            "reading": "(T') Same + Gal(f/R_539)=A5 (disjointness)",
            "exists": "YES expected / conditional",
            "status": "A5 over Q known; Gal over R_539 = A5 if K∩R_539=Q",
        },
        {
            "reading": "(N1) Non-rational f reducing to a fixed seed at ∞ many split P",
            "exists": "NO",
            "status": "ruled out by congruence rigidity",
        },
        {
            "reading": "(N2) Non-rational f matching Frob types of lattice at split p",
            "exists": "plausible / open",
            "status": "any A5 poly over R_539 with same Gal has same Chebotarev types",
        },
        {
            "reading": "(N3) Geometric f over R_539 with lattice specialisations",
            "exists": "open",
            "status": "geometric multi-k / fusion problem",
        },
    ]


def main():
    t0 = time.time()
    print("Q4: f in R_539[x], monodromy A5, split reductions → HQCC Z-lattice", flush=True)

    field = r539_field_data()
    triv = trivial_yes()
    comp = compositum_heuristic_flagship()
    nontriv = nontrivial_reading()
    rigid = reduction_rigidity_lemma()
    lattice = lattice_coefficient_list()
    answers = answers_table()

    proxies = {}
    for n in (7, 11):  # primes dividing 539
        print(f"  proxy n={n}...", flush=True)
        proxies[str(n)] = proxy_split_reduction_fingerprint(n)

    elapsed = round(time.time() - t0, 2)

    # Overall answer for the question as stated (most natural: includes trivial)
    overall = (
        "YES in the natural reading: every HQCC A5 seed f∈Z[x]⊂R_539[x] has monodromy "
        "A5 (over Q, and over R_539 if K∩R_539=Q) and reduces at split primes to the "
        "HQCC Z-lattice (it is the lattice). "
        "NO for non-rational f reducing to a *fixed* seed at infinitely many split primes "
        "(congruence rigidity). "
        "OPEN for non-rational geometric lifts matching lattice only in Frob type or "
        "at finitely many primes."
    )
    print(overall, flush=True)

    lines = [
        r"# Does there exist \(f\in R_{539}[x]\) with monodromy \(A_5\) reducing at split primes to the HQCC \(\mathbb{Z}\)-lattice?",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Overall:** {overall}",
        "",
        "---",
        "",
        r"## 0. Field \(R_{539}\)",
        "",
        f"- \(N={field['N']}={field['factorisation']}\)",
        f"- \(R_{{539}}=\\mathbb{{Q}}(2\\cos 2\\pi/539)\), degree **{field['degree_R_N']}** "
        f"\(=\\varphi(539)/2={field['phi_N']}//2\)",
        f"- Direct computation over \(R_{{539}}\): **not feasible** (deg 210)",
        f"- Proxy divisors used: `{field['proxy_divisors']}`",
        "",
        "---",
        "",
        r"## 1. Trivial reading — **YES**",
        "",
        triv["statement"],
        "",
        r"### HQCC \(A_5\) seeds (all lie in \(R_{539}[x]\))",
        "",
        r"| name | poly | disc □ | Gal status |",
        r"|------|------|:------:|------------|",
    ]
    for r in triv["examples"]:
        lines.append(
            f"| {r['name']} | `{r['poly']}` | {r['disc_square']} | {r['status']} |"
        )

    lines += [
        "",
        f"**{triv['n_A5_examples']}** verified \(A_5\) examples in the sample.",
        "",
        f"**Disjointness caveat:** {triv['caveat_disjointness']}",
        "",
        r"### Compositum heuristic (flagship)",
        "",
        f"- disc primes (sample): `{comp['primes_dividing_disc_sample']}`",
        f"- \(R_{{539}}\) ramifies at: `{comp['R539_ramified_at']}`",
        f"- shared with disc: `{comp['shared_ramification_7_11']}`",
        f"- {comp['note']}",
        "",
        r"### Proxy: flagship Frob at split primes of \(R_7\), \(R_{11}\)",
        "",
    ]
    for n, p in proxies.items():
        lines.append(
            f"- **n={n}**: split primes used={p['primes_used']}, types=`{p['frob_types']}`"
        )
    lines.append("")
    lines.append(proxies["7"]["interpretation"])

    lines += [
        "",
        "---",
        "",
        r"## 2. Non-trivial reading — rigidity vs open lifts",
        "",
        nontriv["statement"],
        "",
        f"**Status:** {nontriv['status']}",
        "",
        r"### Why hard",
        "",
    ]
    for w in nontriv["why_hard"]:
        lines.append(f"- {w}")

    lines += [
        "",
        r"### Congruence rigidity (rules out the strongest lift)",
        "",
        rigid["lemma_sketch"],
        "",
        f"**Consequence:** {rigid['consequence']}",
        "",
        r"### Weaker non-trivial goals still open",
        "",
    ]
    for g in rigid["weaker_nontrivial_goals"]:
        lines.append(f"- {g}")

    lines += [
        "",
        r"### Construction sketches",
        "",
        r"| construction | nontrivial? | works? |",
        r"|--------------|:-----------:|--------|",
    ]
    for c in nontriv["plausible_constructions"]:
        lines.append(
            f"| {c['name']} | {c['nontrivial']} | {c['works']} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## 3. HQCC \(\mathbb{Z}\)-lattice reminder",
        "",
        lattice["note"],
        "",
        f"Model core: `{lattice['model_core']}`",
        "",
        "---",
        "",
        r"## 4. Answer table",
        "",
        r"| reading | exists? | status |",
        r"|---------|---------|--------|",
    ]
    for a in answers:
        lines.append(f"| {a['reading']} | **{a['exists']}** | {a['status']} |")

    lines += [
        "",
        "---",
        "",
        r"## 5. Locked answer",
        "",
        r"> Does there exist \(f\in R_{539}[x]\) with monodromy \(A_5\) reducing at",
        r"> split primes to the HQCC \(\mathbb{Z}\)-lattice?",
        "",
        r"**Yes** — take any HQCC \(A_5\) seed \(f\in\mathbb{Z}[x]\) (e.g. flagship",
        r"\(x^5-55x+88\)). It lies in \(R_{539}[x]\), has monodromy \(A_5\) over",
        r"\(\mathbb{Q}\) (and over \(R_{539}\) if \(K\cap R_{539}=\mathbb{Q}\)), and its",
        r"reduction at every prime — including those that split in \(R_{539}\) — *is*",
        r"the HQCC lattice polynomial mod \(p\).",
        "",
        r"**No** — for a non-rational \(f\in R_{539}[x]\setminus\mathbb{Q}[x]\) that",
        r"reduces to one *fixed* lattice seed at infinitely many primes of \(R_{539}\)",
        r"(congruence rigidity forces \(f\) over \(\mathbb{Z}\)).",
        "",
        r"**Open** — non-rational geometric models over \(R_{539}\) matching the lattice",
        r"only in Frobenius type, at finitely many split primes, or under a single",
        r"embedding (fusion / geometric multi-\(k\) territory).",
        "",
        r"```bash",
        r"python r539_a5_hqcc_lattice.py",
        r"```",
        "",
        r"_Generated by r539_a5_hqcc_lattice.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "overall": overall,
        "field": field,
        "trivial": triv,
        "compositum_heuristic": comp,
        "nontrivial": nontriv,
        "rigidity": rigid,
        "answers": answers,
        "proxies": proxies,
        "lattice": lattice,
    }
    md = "\n".join(lines)
    write_md(ROOT / "R539_A5_HQCC_LATTICE.md", md)
    write_json(ROOT / "R539_A5_HQCC_LATTICE.json", payload)
    write_md(OUT / "R539_A5_HQCC_LATTICE.md", md)
    write_json(OUT / "R539_A5_HQCC_LATTICE.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "R539_A5_HQCC_LATTICE.md", md)
    except Exception:
        pass

    print(f"Wrote R539_A5_HQCC_LATTICE.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
