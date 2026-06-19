"""
Figure generator for the HiPPO and S4D deep-dive tutorials.

Produces publication-quality matplotlib figures (PNG) used by:
  - 01_HiPPO_Deep_Dive.md
  - 02_S4D_Deep_Dive.md
  - README.md

All figures are written to ../../images/ relative to this script
(i.e. semantic_simulation/docs/images/) with a `hippo_` / `s4d_`
prefix so they do not collide with existing images.

Run:
    python make_figures.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.polynomial import legendre as L

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "images"))
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 130,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    }
)

C = {
    "ema": "#c44e52",
    "hippo": "#4c72b0",
    "s4d": "#55a868",
    "accent": "#8172b3",
    "gray": "#7f7f7f",
}


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# 1. Legendre basis functions
# ---------------------------------------------------------------------------
def fig_legendre_basis():
    x = np.linspace(-1, 1, 400)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for n in range(6):
        coef = np.zeros(n + 1)
        coef[n] = 1.0
        # normalized so each has unit L2 norm on [-1, 1]: sqrt((2n+1)/2)
        y = L.legval(x, coef) * np.sqrt((2 * n + 1) / 2.0)
        ax.plot(x, y, label=f"P_{n}", lw=2)
    ax.set_title("Legendre polynomial basis (orthonormal on [-1, 1])")
    ax.set_xlabel("x")
    ax.set_ylabel("normalized P_n(x)")
    ax.axhline(0, color=C["gray"], lw=0.8)
    ax.legend(ncol=3, fontsize=9, loc="upper center")
    save(fig, "hippo_legendre_basis.png")


# ---------------------------------------------------------------------------
# 2. HiPPO online reconstruction with increasing order N
# ---------------------------------------------------------------------------
def _legendre_project(f_vals, x, N):
    """Project a function sampled on x in [-1,1] onto first N Legendre modes."""
    coeffs = []
    for n in range(N):
        c = np.zeros(n + 1)
        c[n] = 1.0
        Pn = L.legval(x, c) * np.sqrt((2 * n + 1) / 2.0)
        coeffs.append(np.trapz(f_vals * Pn, x))
    return np.array(coeffs)


def _legendre_reconstruct(coeffs, x):
    y = np.zeros_like(x)
    for n, cn in enumerate(coeffs):
        c = np.zeros(n + 1)
        c[n] = 1.0
        Pn = L.legval(x, c) * np.sqrt((2 * n + 1) / 2.0)
        y = y + cn * Pn
    return y


def fig_hippo_reconstruction():
    x = np.linspace(-1, 1, 500)
    # a non-trivial target "history window" signal
    f = (
        np.sin(3.0 * np.pi * (x + 1))
        * np.exp(-0.6 * (x + 1))
        + 0.3 * np.sign(np.sin(2.2 * np.pi * x))
    )
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(x, f, color="black", lw=2.4, label="signal on window", zorder=5)
    for N, col, a in [(2, "#d6e1f0", 1), (4, "#9ab4dd", 1), (8, "#4c72b0", 1),
                      (16, "#21314d", 1)]:
        coeffs = _legendre_project(f, x, N)
        rec = _legendre_reconstruct(coeffs, x)
        ax.plot(x, rec, color=col, lw=2, label=f"N = {N} coeffs", alpha=a)
    ax.set_title("HiPPO reconstruction: compress a window into N coefficients")
    ax.set_xlabel("normalized time within window  (-1 = oldest, +1 = newest)")
    ax.set_ylabel("amplitude")
    ax.legend(fontsize=9, loc="upper right")
    save(fig, "hippo_reconstruction.png")


# ---------------------------------------------------------------------------
# 3. LegT vs LegS measures
# ---------------------------------------------------------------------------
def fig_legt_legs_measures():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    t = 6.0  # current time
    s = np.linspace(0, t, 400)

    # LegT: uniform weight on a sliding window of width theta
    theta = 2.5
    legt = np.where(s >= t - theta, 1.0 / theta, 0.0)
    axes[0].fill_between(s, legt, color=C["hippo"], alpha=0.35)
    axes[0].plot(s, legt, color=C["hippo"], lw=2)
    axes[0].axvline(t - theta, color=C["gray"], ls="--", lw=1)
    axes[0].set_title("LegT: sliding window (width theta)")
    axes[0].set_xlabel("past time s  (current t = 6)")
    axes[0].set_ylabel("measure weight")
    axes[0].annotate("forgets older\nthan t - theta", xy=(t - theta, 1.0 / theta),
                     xytext=(0.4, 0.30), fontsize=9,
                     arrowprops=dict(arrowstyle="->", color=C["gray"]))

    # LegS: scaled measure, uniform on [0, t] (weight 1/t), keeps all history
    legs = np.full_like(s, 1.0 / t)
    axes[1].fill_between(s, legs, color=C["s4d"], alpha=0.35)
    axes[1].plot(s, legs, color=C["s4d"], lw=2)
    axes[1].set_title("LegS: scaled, full history (weight 1/t)")
    axes[1].set_xlabel("past time s  (current t = 6)")
    axes[1].annotate("keeps all\nhistory, rescaled", xy=(1.0, 1.0 / t),
                     xytext=(2.0, 0.28), fontsize=9,
                     arrowprops=dict(arrowstyle="->", color=C["gray"]))
    fig.suptitle("HiPPO measures: what counts as 'the window' to approximate",
                 y=1.02, fontsize=12)
    save(fig, "hippo_legt_legs_measures.png")


# ---------------------------------------------------------------------------
# 4. EMA kernels vs HiPPO-LegT kernels: redundant vs orthogonal
# ---------------------------------------------------------------------------
def fig_ema_vs_hippo_kernels():
    lag = np.linspace(0, 60, 400)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))

    # K-EMA: four exponential kernels, all monotone decay (overlapping)
    alphas = [0.0, 0.5, 0.9, 0.99]
    for a in alphas:
        if a == 0.0:
            k = np.zeros_like(lag)
            k[0] = 1.0
            axes[0].plot([0, 0], [0, 1], color=C["ema"], lw=2)
            axes[0].plot(lag, k, color=C["ema"], lw=2,
                         label=f"alpha = {a} (horizon ~1)")
        else:
            tau = -1.0 / np.log(a)
            k = (1 - a) * np.exp(-lag / tau)
            axes[0].plot(lag, k, lw=2,
                         label=f"alpha = {a} (horizon ~{tau:.0f})")
    axes[0].set_title("K-EMA bank: monotone, overlapping (redundant)")
    axes[0].set_xlabel("lag  (tokens into the past)")
    axes[0].set_ylabel("kernel weight")
    axes[0].legend(fontsize=8.5)

    # HiPPO-LegT: Legendre-modulated window kernels (orthogonal, oscillating)
    window = 60.0
    u = 1.0 - 2.0 * lag / window  # map lag in [0,window] -> u in [1,-1]
    inside = lag <= window
    for n in range(4):
        c = np.zeros(n + 1)
        c[n] = 1.0
        Pn = L.legval(u, c) * np.sqrt((2 * n + 1) / 2.0)
        Pn = np.where(inside, Pn, np.nan)
        axes[1].plot(lag, Pn, lw=2, label=f"mode n = {n}")
    axes[1].axhline(0, color=C["gray"], lw=0.8)
    axes[1].set_title("HiPPO-LegT modes: orthogonal over the window")
    axes[1].set_xlabel("lag  (tokens into the past)")
    axes[1].set_ylabel("basis weight")
    axes[1].legend(fontsize=8.5)
    fig.suptitle("Why HiPPO carries more information per vector than K-EMA",
                 y=1.02)
    save(fig, "ema_vs_hippo_kernels.png")


# ---------------------------------------------------------------------------
# HiPPO-LegT / LegS matrices (for eigenvalue plots)
# ---------------------------------------------------------------------------
def hippo_legs(N):
    """Standard HiPPO-LegS A, B (Gu et al. 2020)."""
    A = np.zeros((N, N))
    B = np.zeros(N)
    for n in range(N):
        B[n] = np.sqrt(2 * n + 1)
        for k in range(N):
            if n > k:
                A[n, k] = -np.sqrt(2 * n + 1) * np.sqrt(2 * k + 1)
            elif n == k:
                A[n, k] = -(n + 1)
            else:
                A[n, k] = 0.0
    return A, B


def hippo_legt(N, theta=1.0):
    """HiPPO-LegT A, B (translated Legendre, window theta)."""
    A = np.zeros((N, N))
    B = np.zeros(N)
    for n in range(N):
        B[n] = (2 * n + 1) * ((-1) ** n) / theta
        for k in range(N):
            sgn = 1.0 if n >= k else ((-1) ** (n - k))
            A[n, k] = -(2 * n + 1) * sgn / theta
    return A, B


# ---------------------------------------------------------------------------
# 5. S4D eigenvalues in the complex plane
# ---------------------------------------------------------------------------
def fig_s4d_eigenvalues():
    N = 16
    A_legt, _ = hippo_legt(N, theta=1.0)
    ev_legt = np.linalg.eigvals(A_legt)

    n = np.arange(N // 2)
    # S4D-Lin: A_n = -1/2 + i*pi*n
    s4d_lin = -0.5 + 1j * np.pi * n
    # S4D-Inv: A_n = -1/2 + i * (N/pi) * (N/(2n+1) - 1)
    s4d_inv = -0.5 + 1j * (N / np.pi) * (N / (2 * n + 1) - 1)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.scatter(ev_legt.real, ev_legt.imag, s=70, color=C["hippo"],
               marker="o", label="HiPPO-LegT eigenvalues (dense)", zorder=4,
               edgecolor="white")
    ax.scatter(np.r_[s4d_lin.real, s4d_lin.real],
               np.r_[s4d_lin.imag, -s4d_lin.imag], s=70, color=C["s4d"],
               marker="s", label="S4D-Lin init (diagonal)", zorder=5,
               edgecolor="white")
    ax.scatter(np.r_[s4d_inv.real, s4d_inv.real],
               np.r_[s4d_inv.imag, -s4d_inv.imag], s=55, color=C["accent"],
               marker="^", label="S4D-Inv init (diagonal)", zorder=5,
               edgecolor="white")
    ax.axvline(0, color=C["ema"], lw=1.2, ls="--", label="stability boundary (Re=0)")
    ax.set_xlabel("Re(lambda)  (decay rate)")
    ax.set_ylabel("Im(lambda)  (oscillation frequency)")
    ax.set_title("State-matrix spectra: dense HiPPO vs diagonal S4D")
    ax.set_ylim(-30, 30)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.25)
    ax.text(0.99, 0.02,
            "S4D modes cluster near Re=-1/2 (light, uniform damping);\n"
            "HiPPO-LegT eigenvalues spread to large -Re (strong, varied damping).\n"
            "A few S4D-Inv modes have |Im| > 30 and are off-axis.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color=C["gray"])
    save(fig, "s4d_eigenvalues.png")


# ---------------------------------------------------------------------------
# 6. S4D convolution kernels (damped oscillations)
# ---------------------------------------------------------------------------
def fig_s4d_kernels():
    t = np.linspace(0, 40, 500)
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    modes = [(-0.05, 0.0), (-0.08, 0.6), (-0.12, 1.4), (-0.10, 2.6)]
    for i, (re, im) in enumerate(modes):
        lam = re + 1j * im
        k = np.real(np.exp(lam * t))
        ax.plot(t, k, lw=2, label=f"lambda = {re:+.2f} {im:+.1f}i")
    ax.axhline(0, color=C["gray"], lw=0.8)
    ax.set_title("S4D continuous kernels: K_n(t) = Re(exp(lambda_n t))")
    ax.set_xlabel("t  (tokens into the past)")
    ax.set_ylabel("kernel value")
    ax.legend(fontsize=9)
    save(fig, "s4d_kernels.png")


# ---------------------------------------------------------------------------
# 7. R6 ladder tradeoff (grounding figure)
# ---------------------------------------------------------------------------
def fig_r6_tradeoff():
    cells = ["R6.h.0\nK-EMA", "R6.h.1\nK-EMA log", "R6.a\nLegT fix",
             "R6.e\nLegT learn", "R6.i\nS4D"]
    ppl = [14.78, 15.03, 19.82, 17.45, 16.85]
    keff_ratio = [0.49, 0.51, 0.87, 0.89, 0.91]
    mean_corr = [0.69, 0.67, 0.24, 0.23, 0.20]
    colors = [C["ema"], C["ema"], C["hippo"], C["hippo"], C["s4d"]]

    fig, ax1 = plt.subplots(figsize=(8.4, 4.8))
    x = np.arange(len(cells))
    ax1.bar(x, ppl, color=colors, alpha=0.85, width=0.55)
    ax1.set_ylabel("final val_ppl  (lower is better)", color="black")
    ax1.set_ylim(12, 21)
    for xi, p in zip(x, ppl):
        ax1.text(xi, p + 0.15, f"{p:.2f}", ha="center", fontsize=9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(cells, fontsize=9)

    ax2 = ax1.twinx()
    ax2.plot(x, keff_ratio, "o-", color=C["accent"], lw=2,
             label="K_eff / K (channel richness)")
    ax2.plot(x, mean_corr, "s--", color=C["gray"], lw=2,
             label="mean |corr| (redundancy)")
    ax2.set_ylabel("information diagnostics", color="black")
    ax2.set_ylim(0, 1.0)
    ax2.grid(False)
    ax2.legend(fontsize=9, loc="upper left")

    ax1.set_title("R6 ladder: more orthogonal channels, worse perplexity")
    save(fig, "r6_ladder_tradeoff.png")


if __name__ == "__main__":
    print(f"Writing figures to {OUT}")
    fig_legendre_basis()
    fig_hippo_reconstruction()
    fig_legt_legs_measures()
    fig_ema_vs_hippo_kernels()
    fig_s4d_eigenvalues()
    fig_s4d_kernels()
    fig_r6_tradeoff()
    print("done.")
