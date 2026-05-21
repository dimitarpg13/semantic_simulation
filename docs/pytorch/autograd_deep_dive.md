# A Deep Dive into PyTorch's Autograd Computation Graph

## 1. What autograd actually is

`torch.autograd` is PyTorch's reverse-mode automatic differentiation engine. It is *not* symbolic differentiation (it does not manipulate expressions), and it is *not* finite differences (it does not perturb inputs). It is operator overloading: every differentiable operation on a tensor that "tracks gradients" records itself into a dynamic, directed acyclic graph (DAG), and a single traversal of that graph computes exact gradients of a scalar output with respect to every leaf input via the chain rule.

Three properties matter for everything that follows:

1. **The graph is built during the forward pass.** There is no `tf.Graph()`-style "define-then-run." Every Python statement that touches a tracked tensor extends the graph in place. This is what people mean by "define-by-run" or "dynamic graph."
2. **The graph is consumed (freed) by the backward pass.** By default, after `loss.backward()` returns, the saved intermediate buffers needed for differentiation are released. You opt back in with `retain_graph=True`.
3. **Only scalars (or vectors with an explicit cotangent) can be backpropagated from.** `loss.backward()` works because `loss` is a 0-dim tensor — autograd's reverse mode computes a vector-Jacobian product (VJP), and the seed vector is the implicit `1.0` corresponding to a scalar loss.

## 2. The core data structures

There are really only three objects to understand:

### 2.1 `Tensor`

A tensor that participates in autograd has three relevant attributes:

- `requires_grad: bool` — whether autograd should track operations on this tensor.
- `grad_fn: Optional[Node]` — the `Function` (graph node) that produced this tensor. `None` for leaves.
- `grad: Optional[Tensor]` — populated after `backward()` for leaves with `requires_grad=True`.

A tensor is a **leaf** if it was created directly by the user (e.g., `torch.randn(..., requires_grad=True)`) or by a non-differentiable operation. Equivalently: `t.is_leaf is True` iff `t.grad_fn is None`. Only leaves accumulate gradients into `.grad` by default; for non-leaves you need `.retain_grad()`.

### 2.2 `Function` (the node class)

Every differentiable op corresponds to a subclass of `torch.autograd.Function`. The forward method computes the output; the backward method computes the VJP given an upstream cotangent. The op also stashes whatever inputs/intermediates it needs for backward inside a `ctx` (context) object via `ctx.save_for_backward(...)`.

When you call `y = x ** 2` with `x.requires_grad=True`, PyTorch:

1. Computes `y_data = x.data ** 2`.
2. Allocates a `PowBackward0` node, stores `x` in its saved tensors and `2` as its exponent.
3. Sets `y.grad_fn = <PowBackward0>` and adds an edge from the new node back to `x.grad_fn` (or, if `x` is a leaf, to an `AccumulateGrad` node attached to `x`).

### 2.3 Edges

An edge in the graph is a `(Node, output_nr)` pair: it points to a node and selects which of that node's outputs the consumer is connected to. Multi-output ops (e.g., `torch.svd`) use `output_nr` to disambiguate.

### 2.4 `AccumulateGrad`

Leaves with `requires_grad=True` have an implicit `AccumulateGrad` node sitting at the bottom of the graph. When backward traversal reaches it, the incoming cotangent is *added* to `leaf.grad` (initializing it to zero if `None`). This is why you need `optimizer.zero_grad()` between steps — gradients accumulate by design, since the same leaf may receive contributions from multiple paths through the graph.

## 3. Watching the graph form

Inspecting `grad_fn` and walking `.next_functions` is the most direct way to see the structure.

```python
import torch

x = torch.tensor([2.0, 3.0], requires_grad=True)
w = torch.tensor([0.5, -1.0], requires_grad=True)

z = (x * w).sum()        # z is scalar, suitable for .backward()
print(z.grad_fn)          # <SumBackward0 ...>
print(z.grad_fn.next_functions)
# ((<MulBackward0 ...>, 0),)

mul = z.grad_fn.next_functions[0][0]
print(mul.next_functions)
# ((<AccumulateGrad ...>, 0), (<AccumulateGrad ...>, 0))
# These two AccumulateGrad nodes own x.grad and w.grad respectively.
```

The graph for `z = (x * w).sum()` looks like:

```
   x (leaf)         w (leaf)
       \             /
   AccumulateGrad   AccumulateGrad
       \             /
        \           /
         MulBackward0
              |
         SumBackward0
              |
              z
```

Backward traversal starts at `z`, seeds with `1.0`, computes the VJP for `Sum` (broadcasts `1.0` to shape of the mul output), then the VJP for `Mul` (which produces two gradients — `w` for the `x`-input and `x` for the `w`-input), and routes each into its `AccumulateGrad`.

## 4. The forward pass: graph construction in detail

Three things govern whether a node gets added:

1. **At least one input has `requires_grad=True`** (and the op is differentiable). Outputs inherit `requires_grad=True`.
2. **Grad mode is enabled.** The thread-local flag controlled by `torch.no_grad()`, `torch.enable_grad()`, and `torch.inference_mode()`. Inside `no_grad`, even ops on tracked tensors do not build graph nodes — `out.requires_grad` will be `False`.
3. **The op routes through a dispatcher that supports autograd.** Almost all `torch.*` and tensor methods do; raw `.data` access does not.

### 4.1 Saved tensors

The backward of many ops needs the *value* of an input or output, not just the shapes. `MulBackward0` needs both inputs (the gradient w.r.t. `a` in `a*b` is `b * grad_out`). `PowBackward0` needs the base. `SigmoidBackward0` cleverly saves only the *output* `σ(x)`, since `d/dx σ(x) = σ(x)(1−σ(x))`.

These are accessed via `ctx.saved_tensors` inside `backward`. Importantly, these references *retain memory* — saving `x` for backward keeps `x`'s storage alive at least until backward runs. For long activation chains (transformers!), this is the dominant memory cost of training, which is why activation checkpointing (`torch.utils.checkpoint`) exists: it discards saved tensors and re-runs the forward during backward.

### 4.2 Version counters

Each tensor has a `_version` counter. When a saved tensor is unpacked in backward, autograd checks the version against what it was at save time. If the storage was modified in-place after saving, you get the famous error:

```
RuntimeError: one of the variables needed for gradient computation has been
modified by an inplace operation
```

This is autograd telling you that the value it needs to differentiate is gone. The fix is almost always to remove the in-place op (`x += y` → `x = x + y`) or to clone before the modification.

## 5. The backward pass: chain rule mechanics

`loss.backward()` is sugar for `torch.autograd.backward([loss], [torch.tensor(1.0)])`. The mechanics:

1. The engine builds a topological order of nodes reachable from the seeds, traversing through `grad_fn.next_functions`.
2. It walks the graph in reverse-topological order, calling each node's `backward(*incoming_grads)` to get outgoing grads.
3. Outgoing grads are routed to consumers (parent nodes in the forward graph, which are children in the backward traversal). Multiple incoming grads at a node are summed (this is the chain-rule sum-over-paths).
4. When traversal reaches an `AccumulateGrad` node, the grad is added to the leaf's `.grad`.

For non-scalar outputs, you must provide a cotangent vector explicitly:

```python
y = model(x)           # shape (B, D), requires_grad=True
v = torch.randn_like(y)
y.backward(v)          # computes (v^T J) where J = dy/dparams
```

This computes a **vector-Jacobian product**, not the full Jacobian. The full Jacobian costs `O(output_dim)` backward passes — use `torch.func.jacrev` or `torch.autograd.functional.jacobian` when you really need it.

### 5.1 `torch.autograd.grad` vs `.backward()`

These are two different entry points:

- `.backward()` — populates `.grad` on leaves as a side effect. Returns `None`. Best for training loops.
- `torch.autograd.grad(outputs, inputs, ...)` — returns the gradients as a tuple, does not touch `.grad`. Best for higher-order derivatives, gradient-of-gradient, influence functions, and any setting where you want the grads as ordinary tensors you can differentiate again.

```python
g = torch.autograd.grad(loss, x, create_graph=True)[0]
# g is itself a tensor with a grad_fn; you can call .backward() or .grad() on it.
```

## 6. Higher-order derivatives: `create_graph` and `retain_graph`

Two flags that get confused constantly.

- `retain_graph=True` — keep the *forward* graph (saved tensors, node structure) alive after backward, so you can run `backward()` again. Needed when you want multiple backward passes from the same forward.
- `create_graph=True` — make the *backward* pass itself differentiable by building a graph over the gradient computation. Needed for double-backward, Hessian-vector products, MAML-style meta-learning, etc. Implies `retain_graph=True`.

```python
x = torch.tensor(2.0, requires_grad=True)
y = x ** 3
g = torch.autograd.grad(y, x, create_graph=True)[0]   # g = 3x² = 12, grad_fn alive
g2 = torch.autograd.grad(g, x)[0]                      # g2 = 6x = 12
```

For a Hessian-vector product without materializing the Hessian:

```python
def hvp(loss_fn, params, v):
    grads = torch.autograd.grad(loss_fn(), params, create_graph=True)
    flat = torch.cat([g.reshape(-1) for g in grads])
    hv = torch.autograd.grad(flat @ v, params, retain_graph=False)
    return hv
```

This is exact and cheap (one extra backward), and it's the backbone of conjugate-gradient natural-gradient methods and second-order optimization. Given your work computing acceleration as a derivative-of-derivative quantity along trajectories, this is the right primitive: build the gradient as a graph, then differentiate it directly rather than approximating with finite differences.

## 7. Custom autograd functions

When you need an op autograd doesn't have, or want a fused/cheaper backward, subclass `torch.autograd.Function`. The interface has two static methods.

```python
class StraightThroughSign(torch.autograd.Function):
    """Forward: sign(x). Backward: pretend we computed identity (STE)."""

    @staticmethod
    def forward(ctx, x):
        ctx.save_for_backward(x)
        return torch.sign(x)

    @staticmethod
    def backward(ctx, grad_output):
        (x,) = ctx.saved_tensors
        # Clip the surrogate gradient to a hypercube — standard STE trick.
        return grad_output * (x.abs() <= 1).to(grad_output.dtype)

# Usage
y = StraightThroughSign.apply(x)
```

A few rules:

- `forward` and `backward` are `@staticmethod`. `ctx` is passed explicitly.
- The number of tensors returned from `backward` must equal the number of tensor inputs to `forward`. Return `None` for inputs that don't need a grad (e.g., int hyperparameters).
- Save tensors via `ctx.save_for_backward(...)`, not by attaching them to `ctx` directly — the former participates in the version-counter checks.
- For non-tensor state (shapes, flags), set attributes directly: `ctx.dim = dim`.
- Always test with `torch.autograd.gradcheck(fn, inputs)` — it compares your analytical backward to a high-precision numerical Jacobian. Use `dtype=torch.double` for the inputs.

This is also where you implement **stop-gradient with structure**. A vanilla `.detach()` zeros gradient flow; sometimes you want a forward that uses one expression and a backward that pretends a different expression was computed. STE above is one case; the "reparameterization trick" with discrete samples (Gumbel-softmax STE, etc.) is another.

## 8. Stop-gradient mechanisms compared

There are four ways to prevent gradient flow, with non-obvious differences:

| Mechanism | Scope | Builds graph node? | Notes |
|---|---|---|---|
| `.detach()` | One tensor | No | Returns a view; shares storage. Use for surgical cuts (target networks, EMA teachers, JEPA target branch). |
| `with torch.no_grad():` | Block | No | Disables grad mode for everything inside. Use for inference and parameter updates. |
| `with torch.inference_mode():` | Block | No (and stronger) | Forbids ever using the tensors in autograd later. Faster than `no_grad` because it skips version-counter bookkeeping. Tensors created inside cannot be used in a grad-enabled context later — this is a hard error, not a warning. |
| `param.requires_grad_(False)` | Persistent on tensor | Existing ops still tracked through it if their other inputs require grad — but the param itself won't get a `.grad`. | Use to freeze layers. |

Worth memorizing: `inference_mode` is strictly stronger than `no_grad`. If you create a tensor inside `inference_mode` and try to use it later in an autograd-enabled computation, PyTorch raises. For permanent inference paths (deployed models, replay buffers you'll never differentiate through) it's the right tool; for code that might re-enter grad mode, stick with `no_grad`.

## 9. Hooks: instrumenting the graph

Three kinds, easy to confuse:

```python
# 1. Tensor hook — fires when gradient arrives at this tensor during backward.
h = x.register_hook(lambda grad: grad.clamp_(-1, 1))  # in-place gradient clipping per-tensor

# 2. Module forward/backward hooks — fire around nn.Module forward/backward.
def fwd_hook(module, inputs, outputs): ...
model.layer.register_forward_hook(fwd_hook)

# 3. Full backward hooks — see all grads flowing through a module.
model.layer.register_full_backward_hook(lambda mod, gin, gout: ...)
```

Tensor hooks are the right tool for **debugging vanishing/exploding gradients per-layer**, for collecting gradient statistics without changing the forward, and for surgical interventions (e.g., scaling gradients on a specific path). They have zero effect on the forward graph and are removed by calling `h.remove()` on the returned handle.

A useful pattern for diagnosing where a gradient dies:

```python
for name, p in model.named_parameters():
    p.register_hook(lambda g, n=name: print(n, g.norm().item()))
```

## 10. In-place operations and aliasing

In-place ops (those ending in `_`, plus indexed assignment `x[i] = ...`) are autograd-compatible but constrained:

- **A leaf with `requires_grad=True` cannot be modified in place after it has been used in a graph.** It can be modified before, or you can mutate `.data` (deprecated, dangerous).
- **In-place ops on non-leaf tensors are allowed**, but the version counter bumps, and any saved-for-backward reference to the old version triggers the inplace-modification error at backward time.

The safe pattern is functional: write `x = x + y` rather than `x += y` when `x` participates in autograd. The unsafe-looking-but-fine case is when the in-place op happens *before* the tensor is used in any tracked computation, or *after* backward has run and consumed everything that depended on the old value.

### 10.1 Views and aliasing

`view`, `transpose`, `permute`, `narrow`, `select`, `unsqueeze`, slicing — all return tensors that *share storage* with the source. Autograd tracks this: modifying a view in place is equivalent to modifying the base, and the version counter on the base bumps too. This is why `.detach()` followed by an in-place op silently mutates the original — the detached tensor is a view.

## 11. Memory: the real cost of the graph

For a forward pass that produces activations totaling `A` bytes, the saved-for-backward set is also `O(A)` in the worst case — every nonlinearity, every matmul, every batchnorm saves something. This is why memory profiles during training look like a triangle: activations grow linearly with depth during forward, then are freed in reverse order during backward.

Mitigations, in increasing order of disruption:

1. **`torch.utils.checkpoint.checkpoint(fn, *args)`** — runs `fn` without saving intermediates, then re-runs the forward during backward to recompute them. Trades compute for memory. The standard tool for transformer training at scale.
2. **Gradient accumulation** — split the batch, accumulate grads over micro-batches before stepping. Doesn't reduce per-step graph memory but allows large effective batch sizes.
3. **Mixed precision (`torch.amp`)** — fp16/bf16 activations are half the size; saved tensors shrink accordingly.
4. **FSDP / ZeRO** — shard the parameters and their gradients across devices. Orthogonal to the autograd graph but interacts with it through hooks.

### 11.1 Freeing the graph manually

```python
loss.backward()
# Graph is automatically freed unless retain_graph=True.
# But: anything still referenced (e.g., stored intermediate tensors with .grad_fn)
# keeps part of the graph alive. To force cleanup:
del intermediate_tensor
torch.cuda.empty_cache()  # releases cached blocks back to the allocator
```

A common leak: storing per-step losses *with their graph* in a list (`losses.append(loss)` instead of `losses.append(loss.item())`). Each entry pins its entire backward graph. Always `.item()` or `.detach()` before logging.

## 12. Worked example: implementing AdaGrad with manual gradient management

```python
import torch

torch.manual_seed(0)
x = torch.randn(100, 10)
y = torch.randn(100, 1)

W = torch.randn(10, 1, requires_grad=True)
b = torch.zeros(1, requires_grad=True)

state = {W: torch.zeros_like(W), b: torch.zeros_like(b)}
lr = 0.1

for step in range(200):
    pred = x @ W + b
    loss = ((pred - y) ** 2).mean()

    # Manual gradient computation — does not touch .grad.
    grads = torch.autograd.grad(loss, [W, b])

    # Manual update under no_grad so we don't extend the graph.
    with torch.no_grad():
        for p, g in zip([W, b], grads):
            state[p].add_(g.pow(2))
            p.sub_(lr * g / (state[p].sqrt() + 1e-8))

    if step % 50 == 0:
        print(step, loss.item())
```

Three things to notice:

- `torch.autograd.grad` returns grads without populating `.grad`, so no `zero_grad` is needed.
- The parameter update uses in-place ops (`add_`, `sub_`) inside `no_grad`. This is the standard idiom; the in-place modification of `W` is fine because by the time we update, the forward graph for this step has been fully consumed.
- If we wanted second-order info (e.g., to compute curvature), we'd pass `create_graph=True` to the `grad` call.

## 13. Gotchas worth internalizing

1. **`.grad` accumulates by default.** Forgetting `zero_grad()` adds the new step's gradient to the previous one — which is exactly what you want for gradient accumulation across micro-batches, and exactly what you don't want otherwise.
2. **`tensor.requires_grad_(True)` only works on leaves.** You cannot make a non-leaf require grad — it already does if any of its inputs did.
3. **`.grad` is `None` until backward runs at least once.** Code that does `param.grad.zero_()` before any backward will crash.
4. **Non-leaf grads are discarded.** If you need the gradient of an intermediate, call `.retain_grad()` *before* the backward.
5. **`torch.tensor(other_tensor)` copies and breaks the graph.** Use `.clone()` if you want a differentiable copy.
6. **Integer tensors cannot require grad.** Autograd is defined over floating-point and complex types only.
7. **NaN in gradients usually means a saved activation was extreme, not a bug in your backward.** Add gradient clipping (`torch.nn.utils.clip_grad_norm_`) before assuming the math is wrong.
8. **`copy.deepcopy` of a tensor with `grad_fn` will deepcopy the entire graph too.** Almost never what you want; `.detach().clone()` first.

## 14. When autograd is the wrong tool

For a few use cases there are better primitives now living in `torch.func` (the successor to `functorch`):

- **Per-sample gradients** — `torch.func.vmap(torch.func.grad(fn))` is asymptotically faster than a Python loop calling backward repeatedly.
- **Full Jacobians/Hessians** — `torch.func.jacrev`, `torch.func.jacfwd`, `torch.func.hessian`. These compose cleanly with `vmap`.
- **Functional models** — `torch.func.functional_call(model, params, args)` lets you treat parameters as inputs, which is essential for meta-learning, ensembling, and any computation where you want to differentiate through a parameter update.

These all sit on top of the same autograd engine but expose it through a JAX-style functional interface that's often cleaner for research code.

## 15. Further reading

- The autograd source itself: `torch/csrc/autograd/` in the PyTorch repo. `engine.cpp` is the topological-sort backward executor; `function.h` defines the node interface.
- Paszke et al., *Automatic differentiation in PyTorch* (NeurIPS-W 2017) — the original design doc.
- The `torch.autograd` tutorial in the official docs covers gradcheck and custom Functions in detail.
- For the theoretical foundations of reverse-mode AD: Griewank & Walther, *Evaluating Derivatives* (2008).

---

The mental model that pays off most consistently: **the graph is a record of operations, the backward pass is a single traversal of that record under the chain rule, and everything else — checkpointing, hooks, custom Functions, higher-order grads — is a way to control what gets recorded and how it gets replayed.** Once that's solid, the API stops feeling like a collection of edge cases and starts feeling like the few-dozen-line C++ engine it actually is.
