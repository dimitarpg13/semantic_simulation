# Langevin dynamics

The Langevin equation is a stochastic differential equation introduced by Paul Langevin (1908) to model a particle immersed in a heat bath. It augments deterministic Newtonian mechanics with two effects that summarize the unresolved bath degrees of freedom: a **dissipative drag** that removes energy and a **fluctuating force** that injects it. When the strengths of these two contributions are correctly balanced (the fluctuation–dissipation theorem), the system relaxes to the Gibbs equilibrium distribution at temperature $T$.

The "full" Langevin system is second-order in time and lives on phase space $(q,p)$. There is a singular limit, the **overdamped** Langevin equation, where inertia is dropped and one keeps only configuration variables $q$. The **underdamped** system is the original second-order form and is what people usually mean by "the Langevin SDE" in molecular dynamics, sampling, and (more recently) generative modelling and your own work on transformer hidden-state dynamics.

## Equations of motion

Let $q\in\mathbb{R}^d$ be configuration, $p\in\mathbb{R}^d$ momentum, $M$ a symmetric positive-definite mass matrix, $U:\mathbb{R}^d\to\mathbb{R}$ a smooth potential, $\gamma>0$ a friction coefficient (possibly a matrix $\Gamma$), and $\beta=1/(k_BT)$ inverse temperature. The **underdamped Langevin SDE** in Itô form is

$$
\begin{aligned}
dq_t &= M^{-1} p_t \ dt \\
dp_t &= -\nabla U(q_t) \ dt - \gamma M^{-1} p_t \ dt + \sqrt{2\gamma\beta^{-1}} \ dW_t,
\end{aligned}
$$

with $W_t$ a standard $d$-dimensional Wiener process. The first equation is the kinematic relation; the second is Newton's law plus a linear-in-velocity friction $-\gamma v$ and a white-noise force whose amplitude is fixed by **fluctuation–dissipation**: $\sigma\sigma^\top = 2\gamma\beta^{-1}I$. With $H(q,p)=\tfrac{1}{2}p^\top M^{-1}p + U(q)$ this reads compactly as

$$
\begin{aligned}
dq &= \partial_p H \ dt \\
dp &= -\partial_q H \ dt - \gamma \partial_p H \ dt + \sqrt{2\gamma\beta^{-1}} \ dW.
\end{aligned}
$$

The **overdamped** limit is obtained either by sending $M\to 0$ or rescaling time as $t\to \gamma t$ and letting $\gamma\to\infty$. The momentum is slaved to the force and one is left with the first-order SDE

$$
dq_t = -\nabla U(q_t) \ dt + \sqrt{2\beta^{-1}} \ dW_t.
$$

This is what's used in score-based diffusion samplers and unadjusted Langevin algorithms; the underdamped form (a.k.a. "kinetic Langevin", "second-order Langevin", or "Hamiltonian Monte Carlo with friction") is what's used when inertial mixing matters.

## Analysis

**Invariant measure.** The generator of the underdamped process is

$$
\mathcal{L} = M^{-1}p\cdot\nabla_q - \nabla U(q)\cdot\nabla_p - \gamma M^{-1}p\cdot\nabla_p + \gamma\beta^{-1}\Delta_p,
$$

which splits as $\mathcal{L} = \mathcal{L}_{\text{Ham}} + \gamma \mathcal{L}_{\text{OU}}$ — a Hamiltonian (symplectic, antisymmetric) part plus an Ornstein–Uhlenbeck (dissipative, symmetric) part acting only on momenta. The **Kramers / Fokker–Planck equation** for the density $\rho(q,p,t)$ is

$$
\partial_t\rho = -M^{-1}p\cdot\nabla_q\rho + \nabla U\cdot\nabla_p\rho + \gamma \nabla_p\cdot\big(M^{-1}p \rho + \beta^{-1}\nabla_p\rho\big).
$$

The unique invariant density on phase space is the **canonical (Gibbs) distribution**

$$
\rho_\infty(q,p) \propto e^{-\beta H(q,p)} = e^{-\beta U(q)} e^{-\beta p^\top M^{-1}p/2},
$$

which factorizes as Boltzmann in $q$ times Maxwell in $p$. This is the reason Langevin is used as a sampler: simulating the SDE long enough and time-averaging gives expectations under $e^{-\beta U}$.

**Hypoellipticity and ergodicity.** Noise enters only the momentum equation, so the diffusion matrix is degenerate. Smoothness of densities and convergence to equilibrium nevertheless hold because the commutator $[\nabla_p, M^{-1}p\cdot\nabla_q] = M^{-1}\nabla_q$ recovers the missing directions — Hörmander's condition is satisfied and the operator is **hypoelliptic**. Quantitative convergence rates (exponential in $L^2(\rho_\infty)$, or in Wasserstein-2) are obtained via **hypocoercivity** in the sense of Villani, or via coupling arguments (Eberle, Guillin–Monmarché). The optimal friction $\gamma^*$ trades off two regimes: at small $\gamma$ trajectories are nearly Hamiltonian and mix slowly through energy shells; at large $\gamma$ the dynamics becomes overdamped and configurational mixing slows as $1/\gamma$. The best mixing typically lives at $\gamma$ comparable to the slowest harmonic frequency of $U$.

**Fluctuation–dissipation.** The coefficient $\sigma=\sqrt{2\gamma\beta^{-1}}$ is not free: it is the unique choice that makes $e^{-\beta H}$ stationary. Equivalently, the **Einstein relation** $D=\beta^{-1}\gamma^{-1}$ relates diffusivity, friction, and temperature, and the Green–Kubo formula links transport coefficients to autocorrelations of fluctuating forces.

**Connections you'll recognise.** The underdamped Langevin SDE is the SDE whose drift is the Hamiltonian vector field of $H$ plus the gradient flow of the kinetic energy on momenta, with the gradient flow noised so as to preserve the Gibbs measure — this is exactly the "structure-preserving" picture that motivates BAOAB-type integrators below. In the small-noise / zero-temperature limit, large deviations are governed by the **Freidlin–Wentzell** action, and the minimizer of that action over paths connecting two wells is the analogue of a classical Lagrangian trajectory — closely related to the deceleration/restoring-force picture you've been developing for SPLM.

## Numerical methods

The integrator matters a lot. Naive schemes are stable but bias the invariant measure by $O(\Delta t)$ or $O(\Delta t^2)$ in ways that are not always benign. The state of the art uses **operator splitting** of the generator into pieces that are individually solvable.

**1. Euler–Maruyama.** Simplest, weak order 1, strong order 1/2:

$$
q_{n+1}=q_n+\Delta t \ M^{-1}p_n,\qquad p_{n+1}=p_n-\Delta t \ \nabla U(q_n)-\gamma\Delta t \ M^{-1}p_n+\sqrt{2\gamma\beta^{-1}\Delta t} \ \xi_n,
$$

with $\xi_n\sim\mathcal{N}(0,I)$. Cheap, but has substantial bias in the stationary distribution and poor energy behaviour; rarely the right choice.

**2. Brünger–Brooks–Karplus (BBK).** A symmetric trapezoidal treatment of the friction term layered on velocity Verlet. Better than Euler but still configurational bias $O(\Delta t^2)$.

**3. Stochastic velocity Verlet (van Gunsteren–Berendsen; Vanden-Eijnden–Ciccotti).** Half-kick, full drift, half-kick with the OU step solved exactly. Weak order 2.

**4. Langevin splittings (Leimkuhler–Matthews family).** Split $\mathcal{L}=A+B+O$ where

- $A$: $\dot q = M^{-1}p$ (free drift)
- $B$: $\dot p = -\nabla U(q)$ (force kick)
- $O$: $dp = -\gamma M^{-1}p \ dt + \sqrt{2\gamma\beta^{-1}} \ dW$ (Ornstein–Uhlenbeck on $p$, **solved exactly**: $p\leftarrow e^{-\gamma\Delta t/m}p + \sqrt{(1-e^{-2\gamma\Delta t/m})\beta^{-1} m} \ \xi$).

Different orderings give different schemes; **BAOAB**, i.e. $e^{(\Delta t/2)B}e^{(\Delta t/2)A}e^{\Delta t \ O}e^{(\Delta t/2)A}e^{(\Delta t/2)B}$, is the workhorse. It has the remarkable property — proved by Leimkuhler & Matthews — that the bias on **configurational averages** is $O(\Delta t^4)$ in the high-friction regime, far better than the weak order would suggest. It uses one force evaluation per step, is straightforward to implement, and is the default for sampling-quality MD. ABOBA, OBABO, etc. are alternatives with slightly different stability/accuracy tradeoffs.

**5. Metropolis-adjusted variants.** Wrapping any of the above in an accept/reject step (MALA in the overdamped case, the analogous adjustment for kinetic Langevin) makes the invariant measure exact at the cost of acceptance ratios that can collapse in high dimension. Recent work on **kinetic Langevin Monte Carlo** (Cheng, Chatterji, Bartlett, Jordan; Shen–Lee; Dalalyan–Riou-Durand; Monmarché) gives mixing-time bounds that are quadratically better in dimension than overdamped Langevin Monte Carlo when $U$ is strongly convex and smooth — one of the practical reasons the underdamped form is preferred in ML sampling pipelines.

**6. Exponential / geometric integrators.** When $\nabla U$ has stiff linear part $-Kq$, integrate the harmonic + OU block exactly (Ornstein–Uhlenbeck transition kernels on the affine system are Gaussian) and treat only the nonlinear remainder by a kick. These are highly accurate near minima of $U$ and useful for noise-dominated regimes.

**Practical pointers.** Stability requires $\Delta t \lesssim \omega_{\max}^{-1}$ where $\omega_{\max}$ is the highest harmonic frequency of $U$. The OU step should always be done with the exact Gaussian transition, not an Euler discretization — this is essentially free and removes a large source of bias. For sampling, choose $\gamma$ near the slowest mode; for nonequilibrium / kinetic studies, use the physical friction. For configurational averages BAOAB is usually the right default; for joint $(q,p)$ averages OBABO is often preferred.
