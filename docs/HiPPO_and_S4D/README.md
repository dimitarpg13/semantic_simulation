# HiPPO and S4D: tutorials for the R6 ladder

This folder contains two self-contained tutorial / deep-dive documents on the
state-space memory mechanisms that Section 15 of the Semantic Simulation
paper (the **R6 ladder**) tried as replacements for the K-EMA context pool.

| Document | Covers | R6 cells |
| -------- | ------ | -------- |
| [`01_HiPPO_Deep_Dive.md`](./01_HiPPO_Deep_Dive.md) | HiPPO online polynomial projection, the LegT / LegS measures, the HiPPO ODE and matrices, bilinear discretization, and why K-EMA is the one-pole special case | R6.a, R6.e |
| [`02_S4D_Deep_Dive.md`](./02_S4D_Deep_Dive.md) | S4D diagonal complex state-space models, the HiPPO → S4 → S4D simplification, damped-oscillation kernels, conjugate pairs, ZOH discretization, and the S4D initialisations | R6.i |

Read `01_HiPPO_Deep_Dive.md` first: S4D is initialised from the HiPPO theory
and the K-EMA-as-one-pole-SSM framing is assumed by the S4D document.

## The one-paragraph story

K-EMA, HiPPO, and S4D are all **linear state-space models** for summarising
the past; they differ only in the structure of the state matrix $A$:
diagonal-real (K-EMA), dense-real (HiPPO-LegT), and diagonal-complex (S4D).
Moving up that ladder buys strictly more information per channel and less
redundancy — and yet, in the R6 regime, every structured pool trained to
*worse* perplexity than the redundant K-EMA bank. The bottleneck is the
downstream potential $V_\theta$'s ability to fit the context, not the
information content of the context itself.

## Figures

All figures live in [`../images/`](../images/) with a `hippo_` / `s4d_` /
`ema_` / `r6_` prefix and are produced by
[`scripts/make_figures.py`](./scripts/make_figures.py):

```bash
cd scripts
python make_figures.py
```

| Figure | Used in |
| ------ | ------- |
| `hippo_legendre_basis.png` | HiPPO §4 |
| `hippo_reconstruction.png` | HiPPO §3 |
| `hippo_legt_legs_measures.png` | HiPPO §5 |
| `ema_vs_hippo_kernels.png` | HiPPO §9 |
| `s4d_eigenvalues.png` | S4D §4 |
| `s4d_kernels.png` | S4D §5 |
| `r6_ladder_tradeoff.png` | both §12 |

## Rendering

Both documents use GitHub-flavoured KaTeX and Mermaid written to comply with
[`GitHub_Markdown_LaTeX_Rendering_Cheatsheet.md`](../../semsimula-paper/companion_notes/GitHub_Markdown_LaTeX_Rendering_Cheatsheet.md).
If a symbol looks wrong on GitHub, open the file in Safari.
