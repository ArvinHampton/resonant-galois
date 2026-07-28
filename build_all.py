"""
Build all: full resolution programme (Criteria 1–3) + master RESOLUTION.md.

Usage:
  python build_all.py
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402


def run_step(name: str, fn):
    print("\n" + "=" * 70, flush=True)
    print(f"BUILD STEP: {name}", flush=True)
    print("=" * 70, flush=True)
    t0 = time.time()
    try:
        fn()
        print(f"OK {name} in {time.time() - t0:.1f}s", flush=True)
        return {"step": name, "ok": True, "sec": round(time.time() - t0, 2)}
    except Exception as e:
        traceback.print_exc()
        print(f"FAIL {name}: {e}", flush=True)
        return {"step": name, "ok": False, "error": str(e), "sec": round(time.time() - t0, 2)}


def assemble_catalogues():
    """Merge A5/A6 hit lists into build/CATALOGUE.json."""
    cats = {"A5": [], "A6": [], "D5": [], "notes": []}
    # DEFORM
    p = RESULTS / "DEFORM_M.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        for h in d.get("A5") or []:
            cats["A5"].append({"poly": h["poly"], "disc": h.get("discriminant"), "src": "DEFORM_M"})
        for h in d.get("D5") or []:
            cats["D5"].append({"poly": h["poly"], "disc": h.get("discriminant"), "src": "DEFORM_M"})
    p = RESULTS / "SUMMARY_unique.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        for h in d.get("unique_A5") or []:
            cats["A5"].append({"poly": h["poly"], "disc": h.get("discriminant"), "src": "LATTICE"})
    p = RESULTS / "A6_T6.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        for h in d.get("A6") or []:
            cats["A6"].append({
                "poly": h["poly"],
                "disc": h.get("discriminant"),
                "src": "T6",
                "matrix": (h.get("meta") or {}).get("matrix"),
                "label": (h.get("meta") or {}).get("label"),
            })
    # dedupe
    for k in ("A5", "A6", "D5"):
        ded = {h["poly"]: h for h in cats[k]}
        cats[k] = list(ded.values())
    cats["counts"] = {k: len(cats[k]) for k in ("A5", "A6", "D5")}
    write_json(OUT / "CATALOGUE.json", cats)
    lines = [
        "# Consolidated catalogue",
        "",
        f"- **A5:** {cats['counts']['A5']}",
        f"- **A6:** {cats['counts']['A6']}",
        f"- **D5:** {cats['counts']['D5']}",
        "",
        "## A6 (complete T6 freeze)",
        "",
    ]
    for h in cats["A6"]:
        lines.append(f"- `{h['poly']}` disc={h.get('disc')} label={h.get('label')}")
    lines += ["", "## A5 (sample first 30)", ""]
    for h in cats["A5"][:30]:
        lines.append(f"- `{h['poly']}` src={h.get('src')}")
    if len(cats["A5"]) > 30:
        lines.append(f"- … +{len(cats['A5'])-30} more")
    write_md(OUT / "CATALOGUE.md", "\n".join(lines))
    write_md(RESULTS / "CATALOGUE.md", "\n".join(lines))
    print(f"Catalogue: {cats['counts']}", flush=True)


def write_resolution_master(step_results: list):
    # load criterion summaries if present
    def read(p: Path) -> str:
        return p.read_text(encoding="utf-8") if p.exists() else "_(missing)_"

    c1 = OUT / "CRITERION1_HQCC.md"
    c2 = OUT / "CRITERION2_AXIOMS.md"
    c3 = OUT / "CRITERION3_SIGN.md"
    cat = json.loads((OUT / "CATALOGUE.json").read_text(encoding="utf-8")) if (OUT / "CATALOGUE.json").exists() else {}

    lines = [
        "# RESOLUTION BUILD — closing the gap",
        "",
        "## What a real resolution requires",
        "",
        "1. **Canonical object** (cover / representation / moduli) from HQCC/resonant data",
        "   with **proved** alternating monodromy.",
        "2. **Axioms** on structural matrices forcing disc² + 3-cycles.",
        "3. **Sign character** linked to a ternary invariant so \(\\operatorname{sgn}\\circ\\rho=1\).",
        "",
        "Experiment (catalogues of \(A_5/A_6\)) is **necessary evidence**, not a resolution.",
        "",
        "## Build status",
        "",
        "| Step | OK | sec |",
        "|------|:--:|----:|",
    ]
    for s in step_results:
        lines.append(f"| {s['step']} | {'yes' if s.get('ok') else 'NO'} | {s.get('sec', '—')} |")

    counts = cat.get("counts") or {}
    lines += [
        "",
        "## Experimental catalogues (engine room)",
        "",
        f"- Unique **A5:** {counts.get('A5', '?')}",
        f"- Unique **A6:** {counts.get('A6', '?')}",
        f"- D5 near-misses (sample set): {counts.get('D5', '?')}",
        "",
        "Details: `build/CATALOGUE.md`, `a5_brute_results/A6_T6.md`, `DEFORM_M.md`.",
        "",
        "## Criterion 1 — Canonical HQCC object",
        "",
        "Scaffold: Möbius/HQCC blocks, cubic resultants, BJ/near-rigid families.",
        "See `build/CRITERION1_HQCC.md` and **`build/THEOREM_ATTACK.md`**.",
        "",
        "**Partial advance:** homogenised family \(f_t=x^5+20t^4 x+16t^5\) has proved even monodromy",
        "for all \(t\\neq 0\) (when irr); many specialisations are \(A_5\). Still not HQCC-native.",
        "",
        "**Gap remaining:** no proof of geometric monodromy \(A_n\) for a single *canonical HQCC* cover.",
        "",
        "## Criterion 2 — Axioms ⇒ disc² + 3-cycles",
        "",
        "Evenness obstruction documented (base \(M\), base \(T_6\)).",
        "Subclass rates: `build/CRITERION2_AXIOMS.md`. Thin-class theorems: `THEOREM_ATTACK.md`.",
        "",
        "**Partial advance (lemma):** on the BJ class \(x^5+ax+b\), disc \(=256a^5+3125b^4\);",
        "evenness ⇔ that integer is a square. Homogenised A5 seed is a 1-param theorem class.",
        "",
        "**Gap remaining:** no axiom list proved to force disc² for all structural \(T_n\) matrices.",
        "",
        "## Criterion 3 — Sign character",
        "",
        "Correlations: det(M), ternary weight vs disc² — `build/CRITERION3_SIGN.md`.",
        "",
        "**Partial advance:** Crit 3 **solved on BJ thin class** and on the homogenised A5 family",
        "(closed-form / proved square disc). Full T5 lattice still open.",
        "",
        "**Gap remaining:** no ternary invariant implying trivial sign for *all* HQCC monodromy.",
        "",
        "## Pipeline (unchanged)",
        "",
        "```",
        "template → χ → irr / disc² / cycles → Gal ID on survivors",
        "```",
        "",
        "## Next mathematical moves (priority)",
        "",
        "1. HQCC-native analogue of the homogenised A5 family (Crit 1).",
        "2. Gröbner / ideal form of disc(χ_T5) square in template parameters (Crit 2).",
        "3. Lift BJ sign theorem to a model-flux quadratic character on general T5 (Crit 3).",
        "4. Keep catalogues as **regression tests** for any proposed theorem.",
        "",
        "## Paths",
        "",
        "| Path | Content |",
        "|------|---------|",
        "| `resonant_galois/build/` | This build output |",
        "| `resonant_galois/THEOREM_ATTACK.md` | Theorem-promotion results |",
        "| `resonant_galois/IMPLICATIONS.md` | Claims / non-claims map |",
        "| `a5_brute_results/` | Prior scan archives |",
        "| `resonant_galois/build_all.py` | Rebuild entrypoint |",
        "",
        f"_Build steps: {json.dumps(step_results, indent=2)}_",
    ]
    text = "\n".join(lines)
    write_md(OUT / "RESOLUTION.md", text)
    write_md(RESULTS / "RESOLUTION.md", text)
    write_md(ROOT / "RESOLUTION.md", text)
    print(f"Wrote RESOLUTION.md", flush=True)


def write_readme():
    text = """# Resonant Galois Programme

Structural constructions linking 9 Maths / HQCC ternary–flux arithmetic
to explicit polynomials with Galois groups \(A_5\), \(A_6\), and a research path
toward a **theorem** (not only a generator).

## Quick start

```bash
cd resonant_galois
python build_all.py
```

Outputs land in `build/` and are mirrored into `../a5_brute_results/` for key docs.

## Resolution criteria

| # | Criterion | Module | Doc |
|---|-----------|--------|-----|
| 1 | Canonical HQCC monodromy object | `criterion1_hqcc.py` | `build/CRITERION1_HQCC.md` |
| 2 | Axioms ⇒ disc² + 3-cycles | `criterion2_axioms.py` | `build/CRITERION2_AXIOMS.md` |
| 3 | Sign character / ternary invariant | `criterion3_sign.py` | `build/CRITERION3_SIGN.md` |

Master synthesis: **`RESOLUTION.md`**.  
Implications: **`IMPLICATIONS.md`**.  
Theorem attack: **`THEOREM_ATTACK.md`** (`python theorem_attack.py`).

## Experimental engine (already validated)

```
matrix template (ternary) → char poly → irr / disc² / cycles → galois_group
```

Catalogues: `build/CATALOGUE.md` (A5/A6 freeze).

## Layout

```
resonant_galois/
  build_all.py
  lib/common.py
  criterion1_hqcc.py
  criterion2_axioms.py
  criterion3_sign.py
  build/           # generated
  RESOLUTION.md
  README.md
```

Prior heavy scans remain in `../a5_brute_results/` (DEFORM_M, A6_T6, fingerprints, …).
"""
    write_md(ROOT / "README.md", text)


def main():
    print("RESONANT GALOIS — BUILD ALL", flush=True)
    print(f"ROOT={ROOT}", flush=True)
    print(f"OUT={OUT}", flush=True)
    results = []

    write_readme()
    results.append(run_step("catalogue", assemble_catalogues))

    import criterion2_axioms
    results.append(run_step("criterion2_axioms", criterion2_axioms.main))

    import criterion3_sign
    results.append(run_step("criterion3_sign", criterion3_sign.main))

    import criterion1_hqcc
    results.append(run_step("criterion1_hqcc", criterion1_hqcc.main))

    import theorem_attack
    results.append(run_step("theorem_attack", theorem_attack.main))

    import hqcc_native
    results.append(run_step("hqcc_native", hqcc_native.main))

    results.append(run_step("resolution_master", lambda: write_resolution_master(results)))

    write_json(OUT / "BUILD_LOG.json", results)
    ok = all(r.get("ok") for r in results)
    print("\n" + "=" * 70, flush=True)
    print("BUILD ALL COMPLETE" if ok else "BUILD ALL FINISHED WITH FAILURES", flush=True)
    print(json.dumps(results, indent=2), flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
