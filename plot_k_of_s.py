"""
Reduce 3A^4 geometric fibres to Bring–Jerrard form and plot k(s)=β/α.

Pipeline:
  1. Newton-solve cover params (c,p2,r1,r2) at rational/real s
  2. Sample fibre parameter t; form monic N(y)-t D(y)
  3. Numeric Tschirnhaus reduction to BJ: z^5 + α z + β
  4. k = β/α (complex → real part when imag tiny; else |k| phase noted)
  5. Plot k vs s (median over t) + catalogue multi-seed horizontal lines
  6. Test whether image of k meets ≥2 catalogue ratios

Output: K_OF_S.md / .json, build/k_of_s.png, build/K_OF_S.*
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.linalg import norm

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.common import OUT, RESULTS, write_json, write_md  # noqa: E402

# ---------------------------------------------------------------------------
# Cover Newton (from build_3a4_resolvent)
# ---------------------------------------------------------------------------
def residual(v, s_val: float):
    c, p2, r1, r2, q, w = v

    def G_vals(val, pt):
        yy = pt
        N = c * yy**3 * (yy - 1) * (yy - p2)
        A, Ap, App = yy**3, 3 * yy**2, 6 * yy
        B = yy**2 - (1 + p2) * yy + p2
        Bp = 2 * yy - (1 + p2)
        Np = c * (Ap * B + A * Bp)
        Npp = c * (App * B + 2 * Ap * Bp + A * 2.0)
        D = (yy - r1) * (yy - r2)
        Dp = 2 * yy - (r1 + r2)
        return [N - val * D, Np - val * Dp, Npp - val * 2.0]

    return np.array(G_vals(1.0, q) + G_vals(s_val, w), dtype=float)


def newton(s_val, x0, niter=80):
    v = np.array(x0, dtype=float)
    for _ in range(niter):
        r = residual(v, s_val)
        if norm(r) < 1e-14:
            return v, True, float(norm(r))
        J = np.zeros((6, 6))
        eps = 1e-8
        for j in range(6):
            dv = np.zeros(6)
            dv[j] = eps
            J[:, j] = (residual(v + dv, s_val) - r) / eps
        try:
            step = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(J, -r, rcond=None)[0]
        v = v + step
    nr = float(norm(residual(v, s_val)))
    return v, nr < 1e-12, nr


def solve_at_s(s_val, warm=None, n_trials=40):
    rng = np.random.default_rng(abs(hash(round(float(s_val), 10))) % (2**32))
    cands = []
    if warm is not None:
        cands.append(np.array(warm, dtype=float))
    cands.append(
        np.array([-np.sqrt(5), -1.0, 0.2, -0.2, np.sqrt(5) / 5, -np.sqrt(5) / 5])
    )
    for _ in range(n_trials):
        cands.append(rng.normal(size=6))
    best = None
    for x0 in cands:
        v, ok, nr = newton(s_val, x0)
        if best is None or nr < best[2]:
            best = (v, ok, nr)
        if ok:
            return v, True, nr
    return best[0], False, best[2]


# ---------------------------------------------------------------------------
# Fibre poly → numeric BJ → k
# ---------------------------------------------------------------------------
def fibre_coeffs(c, p2, r1, r2, t_val):
    """
    Monic coefficients of N - t D in y, high to low.
    N = c y^3 (y-1)(y-p2) = c (y^5 - (1+p2)y^4 + p2 y^3)
    D = (y-r1)(y-r2) = y^2 - (r1+r2)y + r1 r2
    N - t D = c y^5 - c(1+p2) y^4 + c p2 y^3 - t y^2 + t(r1+r2) y - t r1 r2
    """
    # degree 5 monic: divide by c
    if abs(c) < 1e-14:
        return None
    # monic: y^5 - (1+p2)y^4 + p2 y^3 - (t/c) y^2 + (t/c)(r1+r2) y - (t/c) r1 r2
    tc = t_val / c
    return np.array(
        [
            1.0,
            -(1 + p2),
            p2,
            -tc,
            tc * (r1 + r2),
            -tc * r1 * r2,
        ],
        dtype=complex,
    )


def depress(coeffs):
    """Kill z^4: monic coeffs -> depressed monic (no z^4)."""
    # y = z - a4/5
    a4 = coeffs[1]
    shift = -a4 / 5
    # binomial expand
    n = 5
    # poly in z: sum c_k (z+shift)^{n-k} wait coeffs high to low c0=1,...,c5
    # p(y)=sum_{k=0}^5 c_k y^{5-k}, y=z+s
    s = shift
    out = np.zeros(6, dtype=complex)
    for k, ck in enumerate(coeffs):
        deg = 5 - k
        # ck (z+s)^deg
        for j in range(deg + 1):
            # binom(deg,j) z^j s^{deg-j} -> power j means index 5-j
            out[5 - j] += ck * sp_binom(deg, j) * (s ** (deg - j))
    out = out / out[0]
    out[1] = 0  # clean float noise on z^4
    return out, complex(s)


def sp_binom(n, k):
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


def _poly_from_roots(roots):
    """Monic poly coeffs high→low from roots (complex-stable product)."""
    p = np.array([1.0 + 0.0j])
    for r in roots:
        p = np.convolve(p, [1.0, -r])
    return p


def to_bring_jerrard(coeffs, n_restarts=12, rng=None):
    """
    Numeric reduction monic quintic → z^5 + α z + β over C.

    Cubic Tschirnhaus on roots: y = z^3 + p z^2 + q z + r
    (3 complex = 6 real parameters) to kill a4,a3,a2 of the monic poly in y
    after re-monicisation / depression — reaches BJ in one step for generic fibres.

    Also: scale-invariant I = β⁴/α⁵; k after α real ≤ 0 normalisation.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    coeffs = np.array(coeffs, dtype=complex)
    dep0, shift0 = depress(coeffs)
    if abs(dep0[2]) < 1e-9 and abs(dep0[3]) < 1e-9:
        return _pack_bj(
            complex(dep0[4]), complex(dep0[5]), 0.0, 0.0, 0.0, "already_BJ", shift0
        )

    # scale original variable
    # monic c5..c0; use root RMS as scale
    roots0 = np.roots(coeffs)
    rms = float(np.sqrt(np.mean(np.abs(roots0) ** 2)))
    if not np.isfinite(rms) or rms < 1e-12:
        rms = 1.0
    roots_s = roots0 / rms

    def monic_from_ys(ys):
        mon = _poly_from_roots(ys)
        if abs(mon[0]) < 1e-30 or not np.all(np.isfinite(mon)):
            return None
        return mon / mon[0]

    def apply_pqr(p, q, r):
        ys = roots_s**3 + p * roots_s**2 + q * roots_s + r
        return monic_from_ys(ys)

    def residual6(v):
        """Kill monic coeffs of y^4, y^3, y^2 (6 real conditions)."""
        p = complex(v[0], v[1])
        q = complex(v[2], v[3])
        r = complex(v[4], v[5])
        mon = apply_pqr(p, q, r)
        if mon is None:
            return np.ones(6) * 1e3
        # mon: [1, a4, a3, a2, a1, a0]
        return np.array(
            [
                mon[1].real,
                mon[1].imag,
                mon[2].real,
                mon[2].imag,
                mon[3].real,
                mon[3].imag,
            ],
            dtype=float,
        )

    def cost(v):
        return float(np.sum(residual6(v) ** 2))

    best_v, best_c = np.zeros(6), 1e300
    n_rand = 180 + 12 * n_restarts
    for i in range(n_rand):
        if i == 0:
            v = np.zeros(6)
        else:
            sc = rng.choice([0.25, 0.6, 1.2, 2.5])
            v = rng.normal(scale=sc, size=6)
        c = cost(v)
        if c < best_c:
            best_c, best_v = c, v.copy()
            if c < 1e-18:
                break

    v = best_v.copy()
    for _ in range(120):
        rvec = residual6(v)
        nr = float(norm(rvec))
        if nr < 1e-12:
            break
        J = np.zeros((6, 6))
        eps = 1e-7
        for j in range(6):
            dv = np.zeros(6)
            dv[j] = eps
            J[:, j] = (residual6(v + dv) - rvec) / eps
        try:
            step = np.linalg.solve(J + 1e-10 * np.eye(6), -rvec)
        except np.linalg.LinAlgError:
            step, *_ = np.linalg.lstsq(J, -rvec, rcond=None)
        sn = float(norm(step))
        if sn > 3:
            step *= 3 / sn
        improved = False
        base = cost(v)
        for fac in (1.0, 0.5, 0.25, 0.1, 0.03):
            v_try = v + fac * step
            if cost(v_try) < base:
                v = v_try
                improved = True
                break
        if not improved:
            break

    p = complex(v[0], v[1])
    q = complex(v[2], v[3])
    r = complex(v[4], v[5])
    mon = apply_pqr(p, q, r)
    if mon is None:
        return {"ok": False, "reason": "eval", "pq_res": best_c}

    # mon = y^5 + 0 y^4 + 0 y^3 + 0 y^2 + α_u y + β_u  (if success)
    a4, a3, a2 = float(abs(mon[1])), float(abs(mon[2])), float(abs(mon[3]))
    alpha_u, beta_u = complex(mon[4]), complex(mon[5])
    alpha, beta = alpha_u, beta_u
    pq_res = float(norm(residual6(v)))
    # encode a4 into a3 channel if needed: require all middle coeffs small
    mid = max(a4, a3, a2)
    return _pack_bj(
        alpha,
        beta,
        mid,  # a3_resid slot := max middle residual
        a2,
        pq_res,
        "cubic_tschirnhaus",
        shift0,
        p,
        q,
        rms,
    )


def _pack_bj(alpha, beta, a3, a2, pq_res, method, shift, p=0, q=0, lam=1.0):
    """Canonicalise scaling so α is real and negative; compute k and I=β⁴/α⁵."""
    if abs(alpha) < 1e-14:
        return {
            "ok": False,
            "reason": "alpha0",
            "a3_resid": a3,
            "a2_resid": a2,
            "pq_res": pq_res,
        }
    # z = μ u with α/μ^4 = -|α|  ⇒  μ^4 = α/(-|α|)
    mu4 = alpha / (-abs(alpha))
    mu = mu4**0.25  # principal 4th root
    alpha_c = alpha / mu**4  # ≈ -|α|
    beta_c = beta / mu**5
    k = beta_c / alpha_c
    I = (beta**4) / (alpha**5)
    ok = (a3 < 1e-5 and a2 < 1e-5) or (pq_res < 1e-6 and a3 < 1e-3 and a2 < 1e-3)
    return {
        "ok": bool(ok),
        "alpha": alpha_c,
        "beta": beta_c,
        "k": k,
        "I": I,
        "a3_resid": a3,
        "a2_resid": a2,
        "pq_res": pq_res,
        "method": method,
        "shift": shift,
        "lam": lam,
        "p": p,
        "q": q,
        "soft": bool(ok and not (a3 < 1e-5 and a2 < 1e-5)),
    }


def k_realish(k):
    if k is None:
        return None
    if abs(k.imag) < 1e-8 * max(1.0, abs(k)):
        return float(k.real)
    return complex(k)


# Catalogue multi-seed ratios
CATALOGUE_K = {
    "flagship": -8 / 5,
    "classical": 4 / 5,
    "lsw": -4.0,
    "s12": -12 / 5,
    "s16": -16 / 5,
    "flag_flip": 8 / 5,
    "class_flip": -4 / 5,
    "lsw_flip": 4.0,
}


def main():
    t0 = time.time()
    print("K(s) FROM 3A^4 FIBRES → BJ", flush=True)

    # s grid: dense real values avoiding 0,1
    s_vals = []
    for num, den in [
        (-5, 1), (-4, 1), (-3, 1), (-5, 2), (-2, 1), (-3, 2), (-1, 1), (-1, 2),
        (-1, 3), (1, 3), (1, 2), (2, 3), (3, 4), (5, 4), (3, 2), (5, 3),
        (2, 1), (5, 2), (3, 1), (7, 2), (4, 1), (5, 1), (-7, 2), (7, 3),
        (8, 3), (-8, 5), (9, 2), (-9, 4), (6, 5), (-6, 5),
    ]:
        sv = Fraction(num, den)
        if sv in (0, 1):
            continue
        s_vals.append(float(sv))
    # denser sample near interesting region
    s_vals += list(np.linspace(-3.5, -0.2, 18))
    s_vals += list(np.linspace(1.2, 4.5, 18))
    s_vals = sorted(set(round(s, 8) for s in s_vals if abs(s) > 1e-9 and abs(s - 1) > 1e-9))

    t_grid = [-2, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4, -3, 2.5]

    rows = []  # per (s,t)
    by_s = defaultdict(list)
    warm = None
    n_ok_s = 0
    n_bj = 0

    for si, s_val in enumerate(s_vals):
        v, ok, nr = solve_at_s(s_val, warm=warm, n_trials=35)
        if not ok:
            print(f"  s={s_val:.4f}: cover FAIL nr={nr:.2e}", flush=True)
            continue
        warm = v
        c, p2, r1, r2, q, w = v
        n_ok_s += 1
        ks_here = []
        for tv in t_grid:
            coeffs = fibre_coeffs(c, p2, r1, r2, tv)
            if coeffs is None:
                continue
            # skip near-branch t where disc degenerates
            try:
                bj = to_bring_jerrard(
                    coeffs,
                    n_restarts=8,
                    rng=np.random.default_rng(si * 17 + (hash(tv) % 1000)),
                )
            except Exception as e:
                bj = {"ok": False, "reason": str(e)[:80]}
            if not bj.get("ok") or bj.get("k") is None:
                continue
            n_bj += 1
            k = bj["k"]
            kr = k_realish(k)
            Ival = bj.get("I")
            rec = {
                "s": float(s_val),
                "t": float(tv),
                "k_real": float(kr) if isinstance(kr, float) else None,
                "k_re": float(np.real(k)),
                "k_im": float(np.imag(k)),
                "I_re": float(np.real(Ival)) if Ival is not None else None,
                "I_im": float(np.imag(Ival)) if Ival is not None else None,
                "alpha_re": float(np.real(bj["alpha"])),
                "beta_re": float(np.real(bj["beta"])),
                "a3_resid": bj.get("a3_resid"),
                "cover_c": float(c),
                "cover_p2": float(p2),
            }
            rows.append(rec)
            if rec["k_real"] is not None:
                ks_here.append(rec["k_real"])
                by_s[float(s_val)].append(rec["k_real"])
        if si % 8 == 0:
            print(
                f"  s={s_val:.4f}: cover ok, BJ fibres with k: {len(ks_here)}",
                flush=True,
            )

    # Aggregate k(s): median real k per s
    s_plot, k_med, k_q25, k_q75, k_mean = [], [], [], [], []
    for s_val in sorted(by_s.keys()):
        arr = np.array(by_s[s_val], dtype=float)
        # clip extreme outliers for plot stats
        if len(arr) == 0:
            continue
        lo, hi = np.percentile(arr, [5, 95])
        clipped = arr[(arr >= lo - 1) & (arr <= hi + 1)] if len(arr) > 4 else arr
        if len(clipped) == 0:
            clipped = arr
        s_plot.append(s_val)
        k_med.append(float(np.median(clipped)))
        k_q25.append(float(np.percentile(clipped, 25)))
        k_q75.append(float(np.percentile(clipped, 75)))
        k_mean.append(float(np.mean(clipped)))

    # Distance of median k(s) to catalogue
    cat_hits = {name: [] for name in CATALOGUE_K}
    tol = 0.08  # absolute tolerance on k for "near hit"
    near = []
    for s_val, km in zip(s_plot, k_med):
        for name, ck in CATALOGUE_K.items():
            if abs(km - ck) < tol:
                cat_hits[name].append({"s": s_val, "k_med": km, "target": ck})
                near.append({"s": s_val, "k_med": km, "catalogue": name, "target": ck})

    multi_cat = sum(1 for name, hits in cat_hits.items() if hits) >= 2
    # also check any individual fibre k near catalogue
    fibre_near = []
    for r in rows:
        if r["k_real"] is None:
            continue
        for name, ck in CATALOGUE_K.items():
            if abs(r["k_real"] - ck) < tol:
                fibre_near.append({**r, "catalogue": name, "target": ck})
                break
    fibre_cats = sorted({f["catalogue"] for f in fibre_near})

    # Plot
    plot_path = OUT / "k_of_s.png"
    plot_path_root = ROOT / "k_of_s.png"
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax = axes[0]
        # scatter all real k
        if rows:
            ss = [r["s"] for r in rows if r["k_real"] is not None]
            kk = [r["k_real"] for r in rows if r["k_real"] is not None]
            ax.scatter(ss, kk, s=8, alpha=0.25, c="C0", label="fibres (s,t)")
        if s_plot:
            ax.plot(s_plot, k_med, "k-", lw=1.8, label="median_t k(s)")
            ax.fill_between(s_plot, k_q25, k_q75, color="C0", alpha=0.15, label="IQR_t")
        for name, ck in CATALOGUE_K.items():
            ax.axhline(ck, ls="--", lw=0.9, alpha=0.7, label=f"cat {name}={ck}")
        ax.set_ylabel(r"$k=\beta/\alpha$ (BJ)")
        ax.set_title(r"3A$^4$ fibres → Bring–Jerrard: $k(s)$")
        ax.set_ylim(-6, 6)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=3, loc="upper right")

        ax2 = axes[1]
        if s_plot:
            ax2.plot(s_plot, k_med, "ko-", ms=3, label="median k(s)")
        for name, ck in list(CATALOGUE_K.items())[:5]:
            ax2.axhline(ck, ls="--", lw=0.9, alpha=0.8)
            ax2.text(s_plot[-1] if s_plot else 0, ck, f" {name}", fontsize=7, va="center")
        ax2.set_xlabel(r"branch cross-ratio $s$")
        ax2.set_ylabel(r"median $k(s)$")
        ax2.set_ylim(-6, 6)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=8)

        fig.tight_layout()
        fig.savefig(plot_path, dpi=140)
        fig.savefig(plot_path_root, dpi=140)
        plt.close(fig)
        plot_ok = True
        print(f"  wrote {plot_path}", flush=True)
    except Exception as e:
        plot_ok = False
        plot_err = str(e)
        print(f"  plot failed: {e}", flush=True)
        # ASCII fallback
        ascii_lines = ["ASCII k(s) (median):"]
        for s_val, km in zip(s_plot, k_med):
            bar = int(np.clip((km + 6) / 12 * 40, 0, 40))
            ascii_lines.append(f"  s={s_val:+6.3f} |{' ' * bar}*{km:+.3f}")
        (OUT / "k_of_s.txt").write_text("\n".join(ascii_lines), encoding="utf-8")
        (ROOT / "k_of_s.txt").write_text("\n".join(ascii_lines), encoding="utf-8")

    elapsed = round(time.time() - t0, 2)
    verdict = (
        f"k(s) from 3A^4 BJ reduction ({elapsed}s). "
        f"covers_ok={n_ok_s}/{len(s_vals)}, BJ_fibres={n_bj}, "
        f"s_with_k={len(s_plot)}. "
        f"median k near ≥2 catalogue ratios: {multi_cat}. "
        f"fibre-level catalogue tags: {fibre_cats}. "
        f"plot={'ok' if plot_ok else 'fail'}."
    )
    print(verdict, flush=True)

    lines = [
        r"# \(k(s)\) from \(3A^4\) fibres reduced to Bring–Jerrard",
        "",
        f"_Elapsed: {elapsed}s_",
        "",
        f"**Verdict:** {verdict}",
        "",
        r"## Method",
        "",
        r"1. Solve cover \(\varphi_s=N/D\) of type \((3,1,1)^4\) at real \(s\) (Newton).",
        r"2. Sample fibre \(t\); monic \(N(y)-t D(y)\).",
        r"3. Numeric cubic Tschirnhaus \(y=z^3+pz^2+qz+r\) → monic \(y^5+\alpha y+\beta\).",
        r"4. \(k=\beta/\alpha\) after canonical scaling (\(\alpha\) real \(\le 0\)); keep real \(k\) when \(\Im k\approx0\).",
        r"5. Plot fibres + median\(_t k(s)\); compare to catalogue multi-seed ratios.",
        "",
        r"## Plot",
        "",
        f"- Image: `build/k_of_s.png`" if plot_ok else f"- Plot failed: {plot_err if not plot_ok else ''}",
        f"- Covers solved: **{n_ok_s}** / {len(s_vals)}",
        f"- BJ reductions with \(k\): **{n_bj}**",
        f"- Distinct \(s\) with real \(k\): **{len(s_plot)}**",
        "",
        r"## Median \(k(s)\) table (sample)",
        "",
        r"| \(s\) | median \(k\) | IQR low | IQR high |",
        r"|------|-------------:|--------:|---------:|",
    ]
    for i in range(0, len(s_plot), max(1, len(s_plot) // 20)):
        lines.append(
            f"| {s_plot[i]:.4f} | {k_med[i]:.4f} | {k_q25[i]:.4f} | {k_q75[i]:.4f} |"
        )

    lines += [
        "",
        r"## Catalogue multi-seed test",
        "",
        f"Tolerance on \(k\): **{tol}**",
        "",
        r"| catalogue \(k\) | value | median-\(k(s)\) near-hits |",
        r"|----------------|------:|-------------------------:|",
    ]
    for name, ck in CATALOGUE_K.items():
        lines.append(f"| {name} | {ck} | {len(cat_hits[name])} |")

    lines += [
        "",
        f"**≥2 catalogue ratios met by median \(k(s)\):** **{multi_cat}**",
        f"**Fibre-level near catalogue tags:** `{fibre_cats}` (n={len(fibre_near)})",
        "",
        r"### Sample near-hits (median)",
        "",
    ]
    if near:
        for h in near[:15]:
            lines.append(
                f"- s={h['s']:.4f}: k_med={h['k_med']:.4f} ≈ {h['catalogue']} ({h['target']})"
            )
    else:
        lines.append("_None within tolerance._")

    lines += [
        "",
        r"### Sample fibre near-hits",
        "",
    ]
    if fibre_near:
        for h in fibre_near[:12]:
            lines.append(
                f"- s={h['s']:.4f}, t={h['t']}: k={h['k_real']:.4f} ≈ {h['catalogue']}"
            )
    else:
        lines.append("_None within tolerance._")

    lines += [
        "",
        r"## Interpretation",
        "",
        r"- \(k(s)\) here is the **numeric BJ ratio** of geometric fibres, not the",
        r"  arithmetic pure-even envelope parameter.",
        r"- Median \(k(s)\) near **two distinct multi-seed ratio classes** (e.g. flagship",
        r"  and classical, not only a sign-flip pair) would be geometric multi-\(k\).",
        r"- Fibre-level hits on classical \(\pm 4/5\) alone are **not** multi-slice multi-\(k\).",
        r"- Arithmetic multi-\(k\) (envelope) remains the completed fusion-level statement.",
        "",
        f"**Median multi-catalogue (≥2):** **{multi_cat}**",
        f"**Fibre tags:** `{fibre_cats}`",
        f"**Geometric multi-\(k\) (strict):** **{multi_cat}**",
        "",
        r"```bash",
        r"python plot_k_of_s.py",
        r"```",
        "",
        r"_Generated by plot_k_of_s.py_",
    ]

    payload = {
        "elapsed_s": elapsed,
        "verdict": verdict,
        "n_s_tried": len(s_vals),
        "n_covers_ok": n_ok_s,
        "n_bj_fibres": n_bj,
        "s_plot": s_plot,
        "k_median": k_med,
        "k_q25": k_q25,
        "k_q75": k_q75,
        "catalogue_k": CATALOGUE_K,
        "tol": tol,
        "median_near_hits": near,
        "fibre_near_hits_n": len(fibre_near),
        "fibre_catalogue_tags": fibre_cats,
        "multi_catalogue_median": multi_cat,
        "multi_catalogue_fibre": len(fibre_cats) >= 2,
        "plot_ok": plot_ok,
        "plot_path": str(plot_path),
        "rows_sample": rows[:40],
    }
    md = "\n".join(lines)
    write_md(ROOT / "K_OF_S.md", md)
    write_json(ROOT / "K_OF_S.json", payload)
    write_md(OUT / "K_OF_S.md", md)
    write_json(OUT / "K_OF_S.json", payload)
    try:
        if RESULTS.exists():
            write_md(RESULTS / "K_OF_S.md", md)
            if plot_ok and plot_path.exists():
                import shutil

                shutil.copy(plot_path, RESULTS / "k_of_s.png")
    except Exception:
        pass

    print(f"Wrote K_OF_S.md ({elapsed}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
