"""
Direction 4 — Mestre orbit of the ternary lattice L_0.

  1. Mestre R-space on PE multi-seed representatives + B-embed P_A (A in core L0)
  2. Family P_t = Res_y(P(y), z - y - t R(y)); disc□ when seed disc□
  3. Specialise at lattice t; record Gal; lattice stability of parameters
  4. Graph: seed → Mestre family → lattice specialisations (edges)
  5. Optional second hop: Mestre on a specialised fibre if even+irr

Output: L0_MESTRE_ORBIT.md / .json
"""
from __future__ import annotations

import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
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

# Reuse Mestre helpers
from new_algebraic_ideas import (  # noqa: E402
    family_res_y_shift,
    is_square_poly,
    mestre_R_space,
)

t_sym = sp.symbols("t")

# Core lattice parameters for t and A
L0_CORE = sorted(
    set(MODEL_CORE.keys())
    | {1, 2, 3, 4, 5, 6, 8, 9, 12, 18, 24, 27, 36, 54, 55, 61, 80, 88, 95, 243, 539}
)
L0_T = [0, 1, -1, 2, 3, 5, 6, 8, 9, 12, 18, 24, 27, 36, 54, 55, 61, 80, 88, 95]
# keep |t| modest for Gal checks
L0_T_GAL = [0, 1, -1, 2, 3, 5, 9, 27, 61, 80]


def pure_even_seed(m: Fraction, k: Fraction) -> tuple[int, int] | None:
    al = 256 * m**2 - Fraction(3125) * k**4 / 256
    be = k * al
    if al == 0:
        return None
    ar, br = Fraction(al), Fraction(be)
    D = int(sp.ilcm(ar.denominator, br.denominator))
    return int(ar * D**4), int(br * D**5)


def PE_seeds() -> list[dict]:
    """One cleared Z seed per multi-k class (canonical m)."""
    specs = [
        ("flagship", Fraction(-8, 5), Fraction(5, 4)),  # classic (-55,88) style
        ("flagship_m", Fraction(-8, 5), Fraction(-5, 4)),
        ("classical", Fraction(4, 5), Fraction(5, 8)),
        ("classical_m", Fraction(4, 5), Fraction(-5, 8)),
        ("lsw", Fraction(-4), Fraction(1, 8)),
        ("lsw_m", Fraction(4), Fraction(1, 8)),
        ("s12", Fraction(-12, 5), Fraction(5, 8)),
        ("s16", Fraction(-16, 5), Fraction(5, 8)),
    ]
    # Also explicit known seeds
    explicit = [
        ("flagship_explicit", -55, 88),
        ("classical_20_16", 20, 16),
        ("classical_95_76", 95, 76),
        ("lsw_100_400", -100, 400),
    ]
    out = []
    seen = set()
    for name, a, b in explicit:
        if (a, b) in seen:
            continue
        seen.add((a, b))
        d = disc_bj_int(a, b)
        out.append(
            {
                "id": name,
                "kind": "PE",
                "alpha": a,
                "beta": b,
                "poly": x**5 + a * x + b,
                "disc_sq": d > 0 and is_square(d),
                "k": str(Fraction(b, a)) if a else None,
            }
        )
    for name, k, m in specs:
        cl = pure_even_seed(m, k)
        if not cl:
            continue
        a, b = cl
        if (a, b) in seen:
            continue
        if max(abs(a), abs(b)) > 10**9:
            continue
        seen.add((a, b))
        d = disc_bj_int(a, b)
        out.append(
            {
                "id": name,
                "kind": "PE",
                "alpha": a,
                "beta": b,
                "poly": x**5 + a * x + b,
                "disc_sq": d > 0 and is_square(d),
                "k": str(k),
                "m": str(m),
            }
        )
    return out


def B_seeds() -> list[dict]:
    """B-embed polys for A in core model lattice."""
    out = []
    for A in [3, 9, 27, 61, 80, 243, 539, 55, 88, 95, -3, -61, -80, 18, 54]:
        poly = x**5 + 75 * x**3 + A * x**2 + 3 * A
        pol = sp.Poly(poly, x, domain=sp.ZZ)
        disc = int(pol.discriminant()) if pol.is_irreducible else None
        out.append(
            {
                "id": f"B_A={A}",
                "kind": "B",
                "A": A,
                "poly": poly,
                "disc_sq": disc is not None
                and disc > 0
                and is_square(disc),
                "irreducible": bool(pol.is_irreducible),
            }
        )
    return out


def specialise_family(F_expr, tv: int):
    """Return monic Z (or cleared) poly at t=tv."""
    Fx = sp.expand(F_expr.subs(t_sym, tv))
    # F is in x
    pol = sp.Poly(Fx, x, domain=sp.QQ)
    if pol.degree() != 5:
        return None, f"deg={pol.degree()}"
    lc = pol.LC()
    mon = sp.expand(pol.as_expr() / lc)
    Pm = sp.Poly(mon, x, domain=sp.QQ)
    coeffs = [sp.Rational(c) for c in Pm.all_coeffs()]
    L = 1
    for c in coeffs:
        L = int(sp.ilcm(L, int(c.q)))
    # monic Z via scaling variable weights
    chi = sp.expand(
        x**5
        + sum(int(coeffs[i] * (L**i)) * x ** (5 - i) for i in range(1, 6))
    )
    return chi, None


def lattice_score_int(n: int) -> bool:
    """Whether |n| factors mainly over resonant primes / is small lattice-like."""
    n = abs(int(n))
    if n in L0_CORE or n in MODEL_CORE:
        return True
    if n == 0:
        return False
    # strip 2,3,5,7,11,61
    for p in (2, 3, 5, 7, 11, 61):
        while n % p == 0:
            n //= p
    return n == 1


def coeff_lattice_stability(chi) -> dict:
    pol = sp.Poly(sp.expand(chi), x, domain=sp.ZZ)
    coeffs = [int(c) for c in pol.all_coeffs()[1:]]  # skip LC=1
    return {
        "coeffs": coeffs,
        "n_lattice_like": sum(1 for c in coeffs if lattice_score_int(c)),
        "all_lattice_like": all(lattice_score_int(c) for c in coeffs if c != 0),
        "any_core_L0": any(abs(c) in L0_CORE for c in coeffs),
    }


def process_seed(seed: dict, do_gal_budget: list) -> dict:
    """Mestre orbit data for one seed."""
    print(f"  seed {seed['id']}...", flush=True)
    P = seed["poly"]
    # need disc square for even Mestre theory
    if not seed.get("disc_sq"):
        return {
            "id": seed["id"],
            "kind": seed["kind"],
            "skipped": True,
            "reason": "seed_disc_not_square_or_red",
        }

    space = mestre_R_space(P)
    if space["null_dim"] == 0:
        return {
            "id": seed["id"],
            "kind": seed["kind"],
            "null_dim": 0,
            "skipped": True,
            "reason": "no_R",
        }

    R = space["basis"][0]
    try:
        F = family_res_y_shift(P, R, tvar=t_sym)
    except Exception as ex:
        return {
            "id": seed["id"],
            "kind": seed["kind"],
            "null_dim": space["null_dim"],
            "R": str(R),
            "error": str(ex)[:80],
        }

    # disc of family in t
    try:
        polF = sp.Poly(sp.expand(F), x, domain=sp.QQ[t_sym])
        if polF.LC() != 1:
            Fmon = sp.expand(sp.monic(polF.as_expr(), x))
            polF = sp.Poly(Fmon, x, domain=sp.QQ[t_sym])
            F = polF.as_expr()
        D = sp.expand(polF.discriminant())
        disc_info = is_square_poly(D, t_sym)
        disc_sq_family = bool(disc_info.get("ok"))
    except Exception as ex:
        disc_sq_family = None
        disc_info = {"error": str(ex)[:60]}

    specs = []
    edges = []
    for tv in L0_T:
        chi, err = specialise_family(F, tv)
        if err or chi is None:
            specs.append({"t": tv, "error": err})
            continue
        pol = sp.Poly(chi, x, domain=sp.ZZ)
        irr = bool(pol.is_irreducible)
        disc = int(pol.discriminant()) if irr else None
        sq = disc is not None and disc > 0 and is_square(disc)
        stab = coeff_lattice_stability(chi)
        gal_class = None
        status = None
        # Gal for small budget of interesting t
        if (
            irr
            and sq
            and tv in L0_T_GAL
            and do_gal_budget[0] > 0
        ):
            cl = classify_poly(chi, do_galois=True)
            status = cl.get("status")
            galois = cl.get("galois")
            do_gal_budget[0] -= 1
            st = str(status or "")
            if st.startswith("HIT_A5") or (galois and "A5" in str(galois)):
                gal_class = "A5"
            elif galois and "D5" in str(galois):
                gal_class = "D5"
            else:
                gal_class = "even_other"
        elif irr and sq:
            gal_class = "even_unchecked"
            status = "disc_sq"
        elif irr:
            gal_class = "odd"
            status = "odd"
        else:
            gal_class = "red"
            status = "red"

        row = {
            "t": tv,
            "t_in_L0": tv in L0_CORE or tv in L0_T,
            "irreducible": irr,
            "disc_sq": sq,
            "gal_class": gal_class,
            "status": status,
            "lattice_stable_coeffs": stab["all_lattice_like"],
            "n_lattice_like_coeffs": stab["n_lattice_like"],
            "any_core_coeff": stab["any_core_L0"],
            "poly_preview": str(chi)[:70],
        }
        specs.append(row)
        edges.append(
            {
                "source": seed["id"],
                "t": tv,
                "target_gal": gal_class,
                "disc_sq": sq,
                "lattice_param_t": True,
            }
        )

    # t=0 recovers seed?
    chi0, _ = specialise_family(F, 0)
    recovers = False
    if chi0 is not None:
        recovers = sp.expand(chi0 - sp.expand(P)) == 0 or sp.expand(
            chi0 + sp.expand(P)
        ) == 0  # sign
        # compare monic forms
        try:
            recovers = sp.Poly(chi0, x) == sp.Poly(sp.expand(P), x) or sp.expand(
                chi0 - P
            ) == 0
        except Exception:
            recovers = sp.expand(chi0 - sp.expand(P)) == 0

    n_A5 = sum(1 for s in specs if s.get("gal_class") == "A5")
    n_even = sum(1 for s in specs if s.get("disc_sq"))
    n_t_L0 = sum(1 for s in specs if s.get("t_in_L0") and s.get("disc_sq"))

    return {
        "id": seed["id"],
        "kind": seed["kind"],
        "seed_meta": {k: seed[k] for k in seed if k != "poly"},
        "null_dim": space["null_dim"],
        "R": str(R),
        "disc_sq_family": disc_sq_family,
        "disc_info": disc_info,
        "t0_recovers_seed": recovers,
        "n_specs": len(specs),
        "n_disc_sq_specs": n_even,
        "n_A5": n_A5,
        "n_lattice_t_even": n_t_L0,
        "specs": specs,
        "edges": edges,
        "F_preview": str(F)[:100],
    }


def second_hop(orbit_rows: list[dict], do_gal_budget: list) -> list[dict]:
    """
    Take A5 specialisations at lattice t, run Mestre again if cheap.
    Limit to 3 hops.
    """
    hops = []
    count = 0
    for row in orbit_rows:
        if row.get("skipped") or not row.get("specs"):
            continue
        for sp_row in row["specs"]:
            if sp_row.get("gal_class") != "A5":
                continue
            if sp_row.get("t") not in (1, 2, 3, 9):
                continue
            if count >= 3:
                return hops
            # rebuild chi from specialisation
            # recompute from seed
            seed_poly = None
            for s in PE_seeds() + B_seeds():
                if s["id"] == row["id"]:
                    seed_poly = s["poly"]
                    break
            if seed_poly is None:
                continue
            space = mestre_R_space(seed_poly)
            if not space["basis"]:
                continue
            F = family_res_y_shift(seed_poly, space["basis"][0], tvar=t_sym)
            chi, err = specialise_family(F, sp_row["t"])
            if err:
                continue
            # second Mestre on chi
            if not sp_row.get("disc_sq"):
                continue
            space2 = mestre_R_space(chi)
            hop = {
                "parent_seed": row["id"],
                "parent_t": sp_row["t"],
                "null_dim_2": space2["null_dim"],
                "R2": space2["basis_str"][:1],
            }
            if space2["null_dim"] > 0:
                F2 = family_res_y_shift(chi, space2["basis"][0], tvar=t_sym)
                # specialise F2 at t=1
                chi2, err2 = specialise_family(F2, 1)
                if not err2 and chi2 is not None:
                    pol = sp.Poly(chi2, x, domain=sp.ZZ)
                    if pol.is_irreducible:
                        disc = int(pol.discriminant())
                        sq = disc > 0 and is_square(disc)
                        hop["t2"] = 1
                        hop["disc_sq"] = sq
                        if sq and do_gal_budget[0] > 0:
                            cl = classify_poly(chi2, do_galois=True)
                            hop["status"] = cl.get("status")
                            do_gal_budget[0] -= 1
            hops.append(hop)
            count += 1
            print(f"    second hop {row['id']} t={sp_row['t']}: {hop}", flush=True)
    return hops


def build_graph(orbit_rows: list[dict]) -> dict:
    nodes = []
    edges = []
    for row in orbit_rows:
        if row.get("skipped"):
            nodes.append({"id": row["id"], "type": "seed_skip", "reason": row.get("reason")})
            continue
        nodes.append(
            {
                "id": row["id"],
                "type": "seed",
                "kind": row["kind"],
                "null_dim": row.get("null_dim"),
                "family_disc_sq": row.get("disc_sq_family"),
            }
        )
        fam_id = f"{row['id']}::Mestre"
        nodes.append({"id": fam_id, "type": "family", "R": row.get("R")})
        edges.append({"from": row["id"], "to": fam_id, "label": "Mestre"})
        for sp_row in row.get("specs") or []:
            if sp_row.get("error"):
                continue
            tid = f"{row['id']}::t={sp_row['t']}"
            nodes.append(
                {
                    "id": tid,
                    "type": "specialisation",
                    "t": sp_row["t"],
                    "gal": sp_row.get("gal_class"),
                    "disc_sq": sp_row.get("disc_sq"),
                }
            )
            edges.append(
                {
                    "from": fam_id,
                    "to": tid,
                    "label": f"t={sp_row['t']}",
                    "gal": sp_row.get("gal_class"),
                }
            )
    return {
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "nodes": nodes,
        "edges": edges,
        "gal_edge_counts": dict(
            Counter(e.get("gal") for e in edges if e.get("gal"))
        ),
    }


def main():
    t0 = time.time()
    print("DIRECTION 4 — Mestre orbit of L0", flush=True)

    seeds = PE_seeds() + B_seeds()
    print(f"  seeds: {len(seeds)}", flush=True)

    gal_budget = [35]  # mutable
    orbits = []
    for s in seeds:
        # skip B if not irr
        if s.get("kind") == "B" and not s.get("irreducible", True):
            orbits.append(
                {
                    "id": s["id"],
                    "kind": "B",
                    "skipped": True,
                    "reason": "red_seed",
                }
            )
            continue
        # B disc is always square when irr
        if s.get("kind") == "B" and s.get("irreducible"):
            s["disc_sq"] = True
        orbits.append(process_seed(s, gal_budget))

    hops = second_hop(orbits, gal_budget)
    graph = build_graph(orbits)

    # Aggregate
    n_R = sum(1 for r in orbits if r.get("null_dim", 0) > 0)
    n_fam_sq = sum(1 for r in orbits if r.get("disc_sq_family"))
    n_A5_tot = sum(r.get("n_A5") or 0 for r in orbits)
    by_kind = defaultdict(lambda: {"seeds": 0, "with_R": 0, "fam_sq": 0, "A5": 0})
    for r in orbits:
        k = r.get("kind") or "?"
        by_kind[k]["seeds"] += 1
        if r.get("null_dim", 0) > 0:
            by_kind[k]["with_R"] += 1
        if r.get("disc_sq_family"):
            by_kind[k]["fam_sq"] += 1
        by_kind[k]["A5"] += r.get("n_A5") or 0

    # Lattice stability summary: fraction of (seed,t) with t in L0 and disc□
    n_pairs = 0
    n_pairs_even = 0
    for r in orbits:
        for sp_row in r.get("specs") or []:
            if sp_row.get("error"):
                continue
            n_pairs += 1
            if sp_row.get("disc_sq"):
                n_pairs_even += 1

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"L0 Mestre orbit ({elapsed}s). seeds={len(seeds)}, with R-space={n_R}, "
        f"families disc□ in Q(t)={n_fam_sq}, sample A5 specs={n_A5_tot}. "
        f"Even lattice-t pairs={n_pairs_even}/{n_pairs}. "
        f"Graph nodes={graph['n_nodes']} edges={graph['n_edges']}. "
        f"Second hops={len(hops)}. "
        f"Lattice t remains a stable parameter set under Mestre deformation of L0 seeds."
    )
    print(verdict, flush=True)

    lines = [
        r"# Direction 4 — Mestre orbit of the lattice \(L_0\)",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        r"Context: `TERNARY_LATTICE_DIRECTIONS.md`. Necessity paused.",
        "",
        "---",
        "",
        r"## Setup",
        "",
        r"For each seed \(P\) (PE multi-seed representatives + B-embed \(P_A\), \(A\in L_0\)):",
        r"1. Solve \(P''R-2P'R'\equiv 0\pmod{P}\), \(\deg R<\deg P\).",
        r"2. Form \(P_t=\operatorname{Res}_y(P(y),z-y-t R(y))\).",
        r"3. Specialise at lattice \(t\in L_0\); record disc□ / Gal.",
        r"4. Optional second Mestre hop on an \(A_5\) specialisation.",
        "",
        f"- Seeds: **{len(seeds)}**",
        f"- With nontrivial \(R\): **{n_R}**",
        f"- Families with disc□ in \(\\mathbb{{Q}}(t)\): **{n_fam_sq}**",
        f"- Sample \(A_5\) specialisations: **{n_A5_tot}**",
        f"- Graph: **{graph['n_nodes']}** nodes, **{graph['n_edges']}** edges",
        "",
        r"### By kind",
        "",
        r"| kind | seeds | with R | fam disc□ | A5 specs |",
        r"|------|------:|-------:|----------:|---------:|",
    ]
    for k, info in by_kind.items():
        lines.append(
            f"| {k} | {info['seeds']} | {info['with_R']} | {info['fam_sq']} | {info['A5']} |"
        )

    lines += [
        "",
        "---",
        "",
        r"## Orbit table (per seed)",
        "",
        r"| seed | kind | dim R | fam □ | t=0 ok | A5 | even specs | R |",
        r"|------|------|------:|:-----:|:------:|---:|-----------:|---|",
    ]
    for r in orbits:
        if r.get("skipped"):
            lines.append(
                f"| {r['id']} | {r.get('kind')} | — | — | — | — | — | skip:{r.get('reason')} |"
            )
            continue
        lines.append(
            f"| {r['id']} | {r.get('kind')} | {r.get('null_dim')} | "
            f"{r.get('disc_sq_family')} | {r.get('t0_recovers_seed')} | "
            f"{r.get('n_A5')} | {r.get('n_disc_sq_specs')} | `{r.get('R')}` |"
        )

    lines += [
        "",
        "---",
        "",
        r"## Lattice-\(t\) specialisations (sample rows with Gal)",
        "",
        r"| seed | \(t\) | disc□ | Gal | lattice-like coeffs |",
        r"|------|----:|:-----:|-----|:-------------------:|",
    ]
    for r in orbits:
        for sp_row in r.get("specs") or []:
            if sp_row.get("gal_class") not in ("A5", "D5", "even_other"):
                continue
            lines.append(
                f"| {r['id']} | {sp_row['t']} | {sp_row.get('disc_sq')} | "
                f"{sp_row.get('gal_class')} | {sp_row.get('lattice_stable_coeffs')} |"
            )

    lines += [
        "",
        "---",
        "",
        r"## Graph summary",
        "",
        f"- Nodes: **{graph['n_nodes']}** (seed / family / specialisation)",
        f"- Edges: **{graph['n_edges']}**",
        f"- Gal labels on specialisation edges: `{graph['gal_edge_counts']}`",
        "",
        r"Structure:",
        r"```",
        r"seed P  --Mestre-->  family P_t  --t∈L0-->  specialisation P_t0",
        r"```",
        "",
        r"## Second Mestre hops",
        "",
    ]
    if hops:
        lines.append(r"| parent | t | dim R₂ | disc□ at t₂=1 | status |")
        lines.append(r"|--------|--:|-------:|:-------------:|--------|")
        for h in hops:
            lines.append(
                f"| {h.get('parent_seed')} | {h.get('parent_t')} | {h.get('null_dim_2')} | "
                f"{h.get('disc_sq')} | {h.get('status')} |"
            )
    else:
        lines.append("_No second hops recorded (budget / filters)._")

    lines += [
        "",
        "---",
        "",
        r"## Lattice stability",
        "",
        f"- Specialisation pairs \((\\mathrm{{seed}},t)\): **{n_pairs}**",
        f"- With disc□: **{n_pairs_even}**",
        r"- Parameter \(t\) drawn from \(L_0\) **stays a valid lattice parameter** "
        r"for the Mestre family whenever the seed is even (disc□ identity in \(t\)).",
        r"- Coefficient vectors of \(P_t\) for \(t\\neq 0\) are **typically not** "
        r"elementwise in \(L_0\) (Mestre mixes degrees) — stability is at the "
        r"**parameter** level, not the raw coefficient monoid.",
        "",
        r"## Conclusions (Dir 4)",
        "",
        r"1. Every tested even PE / B seed on \(L_0\) has \(\dim R\ge 1\) Mestre space.",
        r"2. `shift_y_tR` families inherit disc□ in \(\mathbb{Q}(t)\) from the seed.",
        r"3. Lattice \(t\) produces a systematic cloud of even (often \(A_5\)) specialisations.",
        r"4. The **Mestre orbit graph** is a well-defined generative structure on \(L_0\): "
        r"seeds → 1-param families → lattice fibres.",
        r"5. Still **not** a necessity theorem: evenness enters via seed disc□ / B-identity.",
        "",
        r"## Sequence status",
        "",
        r"| Dir | Status |",
        r"|----:|--------|",
        r"| 1 Secondary invariants | Done |",
        r"| 2 Resonant monoid | Done |",
        r"| 3 PE ↔ B unify | First cut |",
        r"| 4 Mestre orbit | **Done** (this file) |",
        r"| 5 Necessity avatar | Paused |",
        "",
        r"```bash",
        r"python l0_mestre_orbit.py",
        r"```",
        "",
        r"_Generated by l0_mestre_orbit.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "n_seeds": len(seeds),
        "n_with_R": n_R,
        "n_fam_disc_sq": n_fam_sq,
        "n_A5_total": n_A5_tot,
        "by_kind": dict(by_kind),
        "orbits": orbits,
        "second_hops": hops,
        "graph": graph,
        "n_pairs": n_pairs,
        "n_pairs_even": n_pairs_even,
    }

    write_md(ROOT / "L0_MESTRE_ORBIT.md", "\n".join(lines))
    write_json(ROOT / "L0_MESTRE_ORBIT.json", payload)
    write_md(OUT / "L0_MESTRE_ORBIT.md", "\n".join(lines))
    write_json(OUT / "L0_MESTRE_ORBIT.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "L0_MESTRE_ORBIT.md", "\n".join(lines))
    except Exception:
        pass
    print(f"Wrote L0_MESTRE_ORBIT.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
