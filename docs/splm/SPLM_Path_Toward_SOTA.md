# SPLM path toward SOTA: an honest roadmap

An internal memo on whether, and how, the best-performing SPLM variant
(LayerNorm-after-step, §14.16 / Q11) can be scaled toward
state-of-the-art language modelling.

Written 2026-04-24, on top of the §14.15–14.16 empirical results and
the framework argument of §14.8 (R1–R6).

---

## 1. Short answer

**No, not in the benchmark-topping sense of "SOTA". Yes, in a narrower
sense that is still a real contribution.**

The gap between those two readings is where all the interesting
research decisions live. This memo unpacks both.

---

## 2. Where we actually are

The LayerNorm-after-step result is genuinely interesting:

- **88.63** val ppl vs baseline **160.55** — a 45% relative
  improvement.
- At **zero additional parameters** (LayerNorm with
  `ln_affine=False` has no learned parameters).
- At scale $d = 128$, $L = 8$, 7.12 M params, 4000 training steps,
  321 K training tokens (Tiny Shakespeare).

That is an architectural fact about conservative-by-construction LMs,
and it holds. What it is **not** is a SOTA number on any standard
benchmark, because none of those benchmarks are measured on 7 M-param
models trained for 4 k steps on 321 K tokens. Modern SOTA on
WikiText-103, The Pile, C4, or any contemporary eval suite (MMLU,
HumanEval, GSM8K, BBH) is measured at scales three to five orders of
magnitude larger.

---

## 3. What "SOTA" currently means, in compute terms

Four rungs on the ladder, with honest cost estimates:

| Rung | Scale | Corpus | Target | Wall-clock cost |
| --- | --- | --- | --- | --- |
| 0 (now) | 7 M | Tiny Shakespeare (0.32 M tok) | beat own baseline | done |
| 1 | 125 M (GPT-2 small class) | WikiText-103 / OWT-subset (~0.1–1 B tok) | match GPT-2 small (~30–37 ppl on WT-103) | 1–4 A100-days |
| 2 | 1.3 B (Pythia / TinyLlama class) | SlimPajama / C4-subset (~100 B tok) | match Pythia-1B / TinyLlama-1.1B on HELM-lite / lm-eval-harness | 10–30 A100-weeks; ~\$20 k–\$100 k |
| 3 | 7 B–70 B (open-weights SOTA class) | RefinedWeb / Dolma (~1–10 T tok) | match Llama-3-8B / Qwen-2.5-7B on MMLU, HumanEval, GSM8K | $\gtrsim\$1$ M compute |

Rung 3 is what "SOTA" means in the current landscape. It is not a
paper follow-up; it is a year of someone's life and an industrial
compute budget.

---

## 4. Why pure SPLM likely caps before pure-transformer quality at any
given scale

The paper already contains the argument (R5 in §14.8, and the
five-negatives table). Put plainly:

1. **Conservative-by-construction is an expressivity ceiling, not a
   floor.** SPLM force fields are $f = -\nabla V$, which means the
   layer-to-layer update is a gradient field: zero curl, no
   rotational components, no antisymmetric mixing. A transformer
   attention block can produce any linear operator on tokens,
   including strongly asymmetric ones. That extra degree of freedom
   *is* expressivity.
2. **The $\xi$ causal cumulative mean is a very weak attention
   surrogate.** A single vector-valued cumulative mean cannot
   implement content-based routing the way multi-head
   key-query-value attention can. Above a certain scale the $\xi$
   bottleneck becomes limiting.
3. **LayerNorm-after-step is essentially a pre-norm gauge.** At
   $d \geq 2048$, pre-norm transformers already do this per block.
   What looks like a 45% win at $d = 128$ with a free $V_\theta$ head
   may substantially close at larger scale because a large
   transformer already enjoys the same regularization.
4. **Integration cost.** SPLM's per-step forward pass requires
   $\nabla_h V_\theta$, which means a backward pass *inside* the
   forward pass. That is roughly a 2× cost per integration step
   relative to a transformer block of the same width. At scale this
   is paid in either training time or fewer integration steps.
5. **Our own negative-result evidence.** The Gaussian-mixture variant
   (§14.16) collapsed to a context-free unigram predictor. That is
   evidence that *structurally constrained* $V_\theta$ is
   expressivity-limited. The free-$V_\theta$ case is softer, but the
   same failure mode is in the background of any scaling claim.

**Net prediction.** A *pure* SPLM almost certainly cannot match a
same-parameter-budget pure transformer at rung 2 or 3. I would bet on
a 10–30% ppl gap at 1.3 B, widening at 7 B+, if we force the
architecture to stay strictly conservative.

---

## 5. What *is* realistically achievable

Three honest targets, in increasing order of ambition.

### Target A — "SOTA for conservative-by-construction LMs at small scale"

This is what we **already have** after the LayerNorm-after-step
result. The claim is narrow and precise:

> Among conservative-by-construction language models (SPLM family,
> same-parameter budget, Tiny Shakespeare), LayerNorm-after-step with
> SARF + logfreq mass achieves 88.63 val ppl, a 45% relative
> improvement over the SARF + mass baseline (160.55) and the best
> result in the class.

This is defensible, concrete, and falsifiable. It is what §14.16 of
the paper claims. It belongs in the paper as-is.

**Status:** complete. **Cost:** zero. **Deliverable:** the current
paper.

### Target B — "Hybrid SPLM matches GPT-2 small at WikiText-103 scale"

This is the first rung that would look to an outside reader like an
LM-research result rather than an ablation study. Concrete scope:

1. Build **SPLM-ATT**, the hybrid:
   - SPLM integration with LayerNorm-after-step,
   - *plus* multi-head attention as the content-routing layer
     (attention provides the non-conservative content-based routing;
     $V_\theta$ provides the conservative force dynamics;
     LayerNorm-after-step provides the gauge).
2. Scale to $d = 768$, $L = 12$, 125 M params (GPT-2 small).
3. Train on WikiText-103 or an OpenWebText subset (~1 B tokens).
4. Target: within 15% of GPT-2 small's WT-103 perplexity.

**Cost:** 1–4 A100-days of training after a few weeks of development.
Tractable solo.

**Claim at the end:**

> A conservative + content-routing hybrid with LayerNorm-after-step
> reaches within 15% of matched-scale GPT-2 perplexity on
> WikiText-103, demonstrating that conservative-by-construction force
> dynamics can be preserved at GPT-2-small scale without catastrophic
> expressivity loss.

This is publishable as a standalone TMLR follow-up paper and would
cite the current one as the framework foundation.

### Target C — "Hybrid SPLM at 1.3 B parameters, matching Pythia-1B on zero-shot evals"

Rung 2. This is where the claim starts to look like a *model*
contribution rather than a *framework* contribution. Budget: $20 k –
$100 k of cloud compute and 2–3 months of engineering. Achievable
with a small team, a grant, or a compute-provider collaboration. Not
achievable solo in any reasonable timeframe.

**Claim at the end:**

> SPLM-ATT at Pythia-1B scale matches Pythia-1B on lm-eval-harness
> zero-shot, within $X$ % on MMLU, with interpretable per-layer
> attractor structure absent from pure transformer models.

The interpretability angle is the selling point: we get matched
*performance* plus something transformer baselines don't have (a
learned scalar potential $V_\theta$ with extractable semantic
attractors). That is the SPLM program's actual competitive edge at
scale — not raw perplexity, but perplexity-plus-interpretability at
the same parameter budget.

### Target D — "Actual SOTA"

Rung 3. Not achievable at academic-research scale. If a frontier lab
picks up the SPLM idea and bakes it into their next generation, that
is a separate conversation. As an independent research program, this
is not a realistic target.

---

## 6. Recommended sequencing

1. **Keep Target A as the current paper's empirical claim.** It is
   true, well-supported, and honestly scoped. Do not over-reach.
2. **Start Target B as the next paper.** The SPLM-ATT hybrid is a
   ~6-month project with a well-defined deliverable that plausibly
   clears TMLR on its own. The interpretability story (extractable
   attractors, content-basin analysis carried over from
   the attractor analysis pipeline) gives it a
   selling point pure transformers cannot match.
3. **Pitch Target C as a grant or collaboration proposal.** Not a
   solo project. The interpretability angle is what makes it
   fundable: a conservative-by-construction LM at the 1 B scale with
   learned potentials is a mechanistic-interpretability artefact,
   and several labs fund that direction (Anthropic, Apollo,
   MATS-adjacent groups).
4. **Do not claim Target D.** Nobody will believe it and nobody
   should.

---

## 7. The honest framing

> SPLM is not going to beat transformers on raw perplexity at any
> scale, and the paper's own R5 argument says this directly. SPLM's
> competitive claim is perplexity-close at matched budget, *plus* an
> interpretability surface transformers don't have.

Targets A–C each instantiate that claim at a different scale. SOTA in
the raw-perplexity sense is the wrong axis to chase; SOTA in the
**conservative-interpretable-LM class** is the right one, and we are
already there.

---

## 8. Open technical questions that gate Target B

These are the research questions the SPLM-ATT paper would need to
answer. They are worth enumerating here because each has a concrete
checkable answer at ~125 M scale, and the answers together determine
whether Target C is worth pursuing at all.

1. **Does the LN-after-step gain survive the introduction of
   attention?** At small scale it was a 45% ppl win *in absence* of
   attention. With attention providing its own regularization
   (implicit pre-norm, residual stream normalization), the marginal
   gain from LN-after-step may shrink or vanish.
2. **Does the gain survive at $d \geq 768$?** LayerNorm
   regularization strength scales non-trivially with width; what
   works at $d = 128$ may be near-neutral at $d = 2048$.
3. **How many integration steps $L$ are actually needed?** The
   SPLM-ATT forward pass costs roughly $L \cdot 2 \cdot C_\text{block}$
   where $C_\text{block}$ is one attention + MLP block's cost.
   Bringing $L$ down to 2–4 without loss is critical for
   economics at rung 2.
4. **What is the right mass function at scale?** Logfreq surprisal
   is a Shakespeare-scale proxy; at WT-103 scale we likely want a
   learned per-token mass head.
5. **What is the correct causal pooling for $\xi$ at long context?**
   The causal cumulative mean is $O(T)$ in context length, fine at
   256 tokens, potentially too diffuse at 2048 or 8192. Decaying or
   windowed means may be needed.
6. **Does the attractor structure carry over interpretably?** The
   attractor-extraction pipeline in
   the attractor-extraction pipeline was built for
   $d = 128$. Its K-means + silhouette machinery needs to be
   re-validated at $d = 768$ and above; otherwise the
   interpretability selling point of Target C weakens.

A SPLM-ATT paper (Target B) that answers items 1–6 honestly — even
mostly in the negative direction — is a publishable, useful
contribution. It does not need to *win* at rung 1 to be a strong
paper; it needs to *characterize* what conservative-by-construction
buys and loses at rung 1 so that rung-2 decisions are made with
evidence.
