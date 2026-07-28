"""
Gröbner / algebraic analysis of disc(χ_T5) on structural template parameters.

T5(a,b,c,d,e,f) companion+couplings template from the programme:

  χ = x^5 - d x^3 - (a + e f) x^2 - (b f + c e) x + (a d - b c)

Goal: express disc as a polynomial D(a,b,c,d,e,f) and find thin subclasses
where D is identically a square (or a square times a constant square).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, is_square, monic_poly, write_json, write_md, x  # noqa: E402
from lib.lemmas import disc_bj_int  # noqa: E402

a, b, c, d, e, f = sp.symbols("a b c d e f")


def chi_T5(aa=a, bb=b, cc=c, dd=d, ee=e, ff=f):
    """Characteristic polynomial of the structural T5 template."""
    # From exact matrix computation:
    # x^5 - d x^3 - (a+ef) x^2 - (bf+ce) x + (ad-bc)
    return (
        x**5
        - dd * x**3
        - (aa + ee * ff) * x**2
        - (bb * ff + cc * ee) * x
        + (aa * dd - bb * cc)
    )


def disc_of_poly_expr(expr, gens) -> sp.Expr:
    pol = sp.Poly(sp.expand(expr), x)
    return sp.expand(pol.discriminant())


def analyze_thin_ef0() -> dict:
    """
    Thin class e=f=0:
      χ = x^5 - d x^3 - a x^2 + (a d - b c)
    (no x^1 term)
    """
    print("  Thin class e=f=0...", flush=True)
    chi = chi_T5(a, b, c, d, 0, 0)
    D = disc_of_poly_expr(chi, (a, b, c, d))
    Dfac = sp.factor(D)
    # Try write as square * residual
    # Collect as polynomial in (a,b,c,d)
    polyD = sp.Poly(D, a, b, c, d)
    return {
        "class": "e=f=0",
        "chi": str(chi),
        "disc_expanded_terms": len(sp.Add.make_args(sp.expand(D))),
        "disc_factored": str(Dfac),
        "disc_factored_preview": str(Dfac)[:500],
        "total_degree": int(sp.total_degree(D)),
    }


def analyze_thin_bc0() -> dict:
    """b=c=0: χ = x^5 - d x^3 - (a+ef) x^2 + a d  (linear term 0 if b=c=0)."""
    print("  Thin class b=c=0...", flush=True)
    chi = chi_T5(a, 0, 0, d, e, f)
    D = disc_of_poly_expr(chi, (a, d, e, f))
    Dfac = sp.factor(D)
    return {
        "class": "b=c=0",
        "chi": str(chi),
        "disc_factored": str(Dfac),
        "disc_factored_preview": str(Dfac)[:500],
        "total_degree": int(sp.total_degree(D)),
    }


def analyze_thin_bj_slot() -> dict:
    """
    When does χ_T5 reduce to Bring–Jerrard x^5 + p x + q?
    Need: d=0, a+ef=0, and constant / linear match.
    Then disc = 256 p^5 + 3125 q^4 (known).
    """
    print("  Thin class → BJ form...", flush=True)
    # d=0, a = -e f, then
    # χ = x^5 - (b f + c e) x + (-e f * 0 - b c) = x^5 - (bf+ce) x - b c
    # so p = -(bf+ce), q = -bc
    p = -(b * f + c * e)
    q = -b * c
    D_form = 256 * p**5 + 3125 * q**4
    D_simp = sp.factor(sp.expand(D_form))
    # Conditions for BJ embedding in T5
    return {
        "class": "T5 → BJ (d=0, a=-ef)",
        "chi": "x**5 + p x + q with p=-(b f + c e), q=-b c",
        "conditions": ["d = 0", "a = -e*f"],
        "p": str(p),
        "q": str(q),
        "disc_as_BJ": str(D_simp),
        "theorem": (
            "On the subclass d=0, a=-e f of T5, χ is Bring–Jerrard; "
            "even monodromy ⇔ 256 p^5 + 3125 q^4 is a square, with p,q as above. "
            "This is a proved reduction of Crit 2 on a thin T5 subclass to the BJ theorem."
        ),
    }


def analyze_full_symbolic_timeout() -> dict:
    """
    Full 6-param disc — may be large. Compute and factor if possible.
    """
    print("  Full T5 disc (6 params)...", flush=True)
    t0 = time.time()
    chi = chi_T5()
    try:
        D = disc_of_poly_expr(chi, (a, b, c, d, e, f))
        elapsed = time.time() - t0
        nterms = len(sp.Add.make_args(sp.expand(D)))
        print(f"    disc has ~{nterms} terms, {elapsed:.1f}s", flush=True)
        t1 = time.time()
        Dfac = sp.factor(D)
        fac_time = time.time() - t1
        print(f"    factored in {fac_time:.1f}s", flush=True)
        # Extract square-free part / content
        # Try to see if D = S^2 * R for some polynomial S
        # Use square-free decomposition in one variable slices
        return {
            "class": "full T5",
            "chi": str(chi),
            "disc_nterms": nterms,
            "disc_total_degree": int(sp.total_degree(D)),
            "disc_factored_preview": str(Dfac)[:800],
            "disc_factored_len": len(str(Dfac)),
            "elapsed_sec": round(elapsed + fac_time, 2),
            "raw_disc_saved": False,
        }
    except Exception as ex:
        return {"class": "full T5", "error": str(ex)}


def search_square_slices() -> dict:
    """
    Fix HQCC values in most slots; free one or two params; ask when D is square.
    Also detect parametric lines where D is identically a square polynomial in t.
    """
    print("  Parametric slices for square disc...", flush=True)
    results = []

    # Line: known A5 matrix deformations already in catalogue — symbolic identity search
    # Family: a=3, b=s, c=61, d=-3, e=0, f=0 — D as poly in s
    s = sp.symbols("s")
    lines = [
        ("a=3,b=s,c=61,d=-3,e=f=0", chi_T5(3, s, 61, -3, 0, 0), s),
        ("a=s,b=3,c=61,d=-3,e=f=0", chi_T5(s, 3, 61, -3, 0, 0), s),
        ("a=3,b=80,c=s,d=-3,e=f=0", chi_T5(3, 80, s, -3, 0, 0), s),
        ("a=3,b=3,c=3,d=-3,e=s,f=0", chi_T5(3, 3, 3, -3, s, 0), s),
        ("BJ embed: d=0,a=0,b=s,c=1,e=0,f=p  wait", None, s),
        # BJ embed: d=0, a=0 (so e=0 or f=0), b=α, c=β → q=-αβ, p=0 if e=f=0 → x^5 - αβ constant only
        ("BJ: d=0,a=0,e=0,f=0,b=s,c=u fixed u=1", chi_T5(0, s, 1, 0, 0, 0), s),
        # pure BJ via e: d=0, a=0, b=0, c=1, e=s, f=0 → p=-s, q=0 → x^5 - s x
        ("BJ: d=0,a=0,b=0,c=1,e=s,f=0", chi_T5(0, 0, 1, 0, s, 0), s),
        # d=0, a=-e f with e=1, f=s, b=3, c=61: p=-(3s+61), q=-183
        ("BJ embed e=1,f=s,b=3,c=61,d=0,a=-s", chi_T5(-s, 3, 61, 0, 1, s), s),
        ("BJ embed e=3,f=s,b=61,c=80,d=0,a=-3s", chi_T5(-3 * s, 61, 80, 0, 3, s), s),
        ("BJ embed e=s,f=s,b=3,c=3,d=0,a=-s**2", chi_T5(-(s**2), 3, 3, 0, s, s), s),
    ]

    for name, chi, svar in lines:
        if chi is None:
            continue
        print(f"    line {name}", flush=True)
        D = disc_of_poly_expr(chi, (svar,))
        Dfac = sp.factor(sp.expand(D))
        # Is D a square in Q[s]?
        # Check square-free decomposition
        try:
            poly = sp.Poly(sp.expand(D), svar, domain=sp.QQ)
            # content
            cont = poly.content()
            prim = poly.primitive()[1]
            # factor over Q
            fac = sp.factor(prim.as_expr())
            # numerical: for many s, is D square?
            sq_count = 0
            irr_count = 0
            a5_count = 0
            samples = []
            for tv in list(range(-12, 13)) + [16, 27, 61, 80, 243]:
                if tv == 0 and "s**2" in name:
                    pass
                Dv = int(sp.expand(D.subs(svar, tv)))
                if Dv > 0 and is_square(Dv):
                    sq_count += 1
                    expr = sp.expand(chi.subs(svar, tv))
                    pol = monic_poly(expr)
                    if pol and pol.is_irreducible:
                        irr_count += 1
                        from lib.common import classify_poly

                        rec = classify_poly(expr, do_galois=True)
                        if (rec.get("status") or "").startswith("HIT_A5") or (
                            rec.get("galois") and "A5" in str(rec.get("galois"))
                        ):
                            a5_count += 1
                            samples.append({"t": tv, "poly": rec["poly"]})
            # identical square? check if D = const * (poly)^2
            identical_sq = False
            try:
                # sqrt as polynomial?
                # sp.sqrt only works if perfect square in the UFD
                sqf = sp.squarefree_decomposition(sp.Poly(sp.expand(D), svar))
                # if all exponents even in squarefree decomp of monic part...
                # simpler: factor and check even mults
                identical_sq = _is_polynomial_square(sp.expand(D), svar)
            except Exception:
                identical_sq = False

            results.append({
                "line": name,
                "disc_factored": str(Dfac)[:400],
                "identical_square_in_s": identical_sq,
                "n_square_at_int": sq_count,
                "n_irr_when_sq": irr_count,
                "n_A5": a5_count,
                "A5_sample": samples[:6],
            })
            print(
                f"      identical_sq={identical_sq} sq_hits={sq_count} A5={a5_count}",
                flush=True,
            )
        except Exception as ex:
            results.append({"line": name, "error": str(ex)})

    return {"lines": results}


def _is_polynomial_square(expr, var) -> bool:
    """Return True if expr is a square in Q[var] (up to units / rational content)."""
    if expr == 0:
        return True
    P = sp.Poly(sp.expand(expr), var, domain=sp.QQ)
    # content as rational square?
    cont = sp.QQ(P.content())
    # factor content numerator/denominator
    c_num, c_den = cont.p, cont.q
    if not (sp.integer_nthroot(abs(int(c_num)), 2)[1] and sp.integer_nthroot(int(c_den), 2)[1]):
        # allow rational squares: cont = ±(r/s)^2 * unit — check valuation of primes
        if not _rational_is_square(cont):
            return False
    prim = P.primitive()[1]
    fac = sp.factor_list(prim.as_expr())
    # fac = (content, [(factor, mult), ...])
    for factor, mult in fac[1]:
        if mult % 2 != 0:
            return False
    return True


def _rational_is_square(q) -> bool:
    q = sp.QQ(q)
    if q < 0:
        return False
    n, d = int(q.p), int(q.q)
    return bool(sp.integer_nthroot(n, 2)[1] and sp.integer_nthroot(d, 2)[1])


def groebner_square_conditions_bj_embed() -> dict:
    """
    On BJ-embed subclass, D = 256 p^5 + 3125 q^4.
    Seek polynomial relations that force D to be a square identically in free params.
    E.g. set p = 20 u^4, q = 16 u^5 → D square (classical).
    HQCC: p = -55 u^4, q = 88 u^5.
    """
    print("  BJ-embed parametric square families...", flush=True)
    u = sp.symbols("u")
    families = []
    seeds = [
        (20, 16, "classical"),
        (20, -16, "classical_flip"),
        (-55, 88, "hqcc_55_88"),
        (-55, -88, "hqcc_55_m88"),
        (95, 76, "hqcc_95_76"),
        (95, 532, "hqcc_95_532"),
        (-100, 400, "hqcc_100_400"),
        (124, 496, "hqcc_124_496"),
    ]
    for alpha, beta, name in seeds:
        # Realize in T5 BJ embed: p=α u^4, q=β u^5
        # Need -(b f + c e) = α u^4, -b c = β u^5, d=0, a=-e f
        # Simple realization: e=0, f=1, b = -α u^4, c = -β u^5 / b ... messy if b varies
        # Cleaner: e=0, f=0 is impossible for p≠0
        # Use: b=1, c=-β u^5, e=0, f=-α u^4 → p=-(1*(-α u^4)+0)=α u^4, q=-1*(-β u^5)=β u^5
        # a=-e f=0, d=0
        chi = chi_T5(0, 1, -beta * u**5, 0, 0, -alpha * u**4)
        # Actually p = -(b f + c e) = -(1*(-α u^4) + (-β u^5)*0) = α u^4
        # q = -b c = -1*(-β u^5) = β u^5
        D = disc_of_poly_expr(chi, (u,))
        seed_d = disc_bj_int(alpha, beta)
        expected = (u**20) * seed_d
        ok = sp.expand(D - expected) == 0
        families.append({
            "name": name,
            "alpha": alpha,
            "beta": beta,
            "T5_realization": "a=0,d=0,e=0,b=1,c=-beta u^5,f=-alpha u^4",
            "disc_identity_ok": bool(ok),
            "proved_even": bool(ok and is_square(seed_d)),
        })
        print(f"    {name}: identity_ok={ok}", flush=True)
    return {
        "families": families,
        "theorem": (
            "Every BJ homogenised family embeds in the T5 template under "
            "d=0, a=-e f (here e=0, a=0), b=1, c=-β u^5, f=-α u^4. "
            "Thus the T5 structural template *contains* all proved-even "
            "homogenised BJ families as a thin parametric subclass."
        ),
    }


def write_doc(blob: dict) -> str:
    lines = [
        "# T5 discriminant ideal / thin subclasses",
        "",
        f"_Elapsed: {blob.get('elapsed_sec')}s_",
        "",
        "## Setup",
        "",
        "Structural template \(T_5(a,b,c,d,e,f)\) has characteristic polynomial",
        "",
        "```",
        "χ = x^5 - d x^3 - (a + e f) x^2 - (b f + c e) x + (a d - b c)",
        "```",
        "",
        "Even monodromy ⇔ \(D := \\operatorname{disc}(\\chi)\) is a square in \(\\mathbb{Z}\).",
        "",
        "---",
        "",
        "## Theorem: BJ embed in T5",
        "",
    ]
    bj = blob.get("bj_embed") or {}
    lines.append(f"- conditions: `{bj.get('conditions')}`")
    lines.append(f"- p = `{bj.get('p')}`, q = `{bj.get('q')}`")
    lines.append(f"- status: {bj.get('theorem')}")
    lines.append("")

    emb = blob.get("embedding") or {}
    lines.append("### Homogenised families inside T5")
    lines.append(f"{emb.get('theorem')}")
    lines.append("")
    lines.append("| family | α | β | disc identity | proved even |")
    lines.append("|--------|--:|--:|:---:|:---:|")
    for fam in emb.get("families") or []:
        lines.append(
            f"| {fam['name']} | {fam['alpha']} | {fam['beta']} | "
            f"{fam['disc_identity_ok']} | {fam['proved_even']} |"
        )
    lines.append("")

    lines += ["---", "", "## Thin-class discriminants", ""]
    for key in ("ef0", "bc0", "full"):
        block = blob.get(key) or {}
        lines.append(f"### {block.get('class', key)}")
        if block.get("error"):
            lines.append(f"- error: {block['error']}")
        else:
            if block.get("chi"):
                lines.append(f"- χ: `{block['chi']}`")
            if block.get("total_degree") is not None:
                lines.append(f"- total degree of D: {block['total_degree']}")
            if block.get("disc_nterms") is not None:
                lines.append(f"- expanded terms: {block['disc_nterms']}")
            prev = block.get("disc_factored_preview") or block.get("disc_factored") or ""
            lines.append(f"- disc (preview): `{prev[:400]}`")
        lines.append("")

    lines += ["---", "", "## Parametric slices (when is D a square?)", ""]
    for L in (blob.get("slices") or {}).get("lines") or []:
        if L.get("error"):
            lines.append(f"- `{L.get('line')}`: error {L['error']}")
            continue
        lines.append(
            f"- **{L['line']}**: identical_square={L.get('identical_square_in_s')} "
            f"int_sq={L.get('n_square_at_int')} A5={L.get('n_A5')}"
        )
        if L.get("disc_factored"):
            lines.append(f"  - D ~ `{L['disc_factored'][:200]}`")
        for s in L.get("A5_sample") or []:
            lines.append(f"  - A5 t={s['t']}: `{s['poly']}`")
    lines.append("")

    lines += [
        "---",
        "",
        "## Status vs evenness obstruction",
        "",
        "- Base structural \(M\) (\(a=3,b=80,c=61,d=-3,e=f=0\)) remains **odd** — not in BJ embed (\(d\\neq 0\)).",
        "- All **proved-even** infinite families currently live in the **BJ embed** thin subclass of T5 "
        "(or pure BJ polynomials without matrix lift).",
        "- Full \(D(a,b,c,d,e,f)\) is a single polynomial condition; forcing it to be a square for *all* "
        "parameters is impossible (obstruction examples). Forcing it on a **subvariety** is the theorem path.",
        "",
        "## Next algebraic moves",
        "",
        "1. Square-free factorisation / Gröbner of the ideal of coefficients that make \(D\) a square "
        "in a weighted projective sense (hard Diophantine).",
        "2. Search for non-BJ thin classes (e.g. palindromic, self-adjoint) where \(D=S^2\) identically.",
        "3. Keep BJ-embed as the proved spine inside T5; extend geometric monodromy for HQCC covers separately.",
        "",
        "_Generated by t5_disc_ideal.py_",
    ]
    return "\n".join(lines)


def main():
    t0 = time.time()
    print("T5 DISC IDEAL analysis", flush=True)

    bj_embed = analyze_thin_bj_slot()
    embedding = groebner_square_conditions_bj_embed()
    ef0 = analyze_thin_ef0()
    bc0 = analyze_thin_bc0()
    full = analyze_full_symbolic_timeout()
    slices = search_square_slices()

    blob = {
        "elapsed_sec": round(time.time() - t0, 2),
        "bj_embed": bj_embed,
        "embedding": embedding,
        "ef0": ef0,
        "bc0": bc0,
        "full": full,
        "slices": slices,
    }
    doc = write_doc(blob)
    write_md(OUT / "T5_DISC_IDEAL.md", doc)
    write_md(RESULTS / "T5_DISC_IDEAL.md", doc)
    write_md(ROOT / "T5_DISC_IDEAL.md", doc)
    # JSON may be large — drop huge factored strings
    slim = {
        "elapsed_sec": blob["elapsed_sec"],
        "bj_embed": bj_embed,
        "embedding": embedding,
        "ef0": {k: v for k, v in ef0.items() if k != "disc_factored"},
        "bc0": {k: v for k, v in bc0.items() if k != "disc_factored"},
        "full": {k: v for k, v in full.items() if "factored" not in k or "preview" in k},
        "slices": slices,
    }
    write_json(OUT / "T5_DISC_IDEAL.json", slim)
    print(f"Wrote T5_DISC_IDEAL.md in {blob['elapsed_sec']}s", flush=True)
    return blob


if __name__ == "__main__":
    main()
