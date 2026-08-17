---
date: 2026-03-10
title: "GEM vs full EM: LL decrease under linked CTMC specification"
category: "data-quality"
language: "R"
tags: [em-algorithm, gem, generalised-em, linked-specification, ctmc, log-likelihood, monotonicity, mixture-model]
root-cause: "Closed-form M-step for Markov parameters only maximises the Markov block of Q; CTMC link propagates theta -> lambda into EMG emission terms not covered by that maximisation, breaking the full EM guarantee."
severity: "P1"
---

# GEM vs full EM: LL decrease under linked CTMC specification

## Problem

An EM algorithm for a 3-wave binary panel model with tenure durations produced
**decreasing log-likelihood** on early iterations:

```
iter 1 | ll = -477.5
iter 2 | ll = -479.2   ← decreased by 1.7
iter 3 | ll = -478.8   ← recovered slightly
iter 4 | ll = -478.1
...
iter N | converged
```

The algorithm still converged to a reasonable solution, but the LL violation
triggered warnings and undermined confidence in the implementation.

## Root Cause

### Model structure

The model has 6 free parameters: (α, θ₁, θ₀, π, σ²_g, σ²_d). The exponential
duration rates are **linked** deterministically to the transition probabilities
via the continuous-time Markov chain (CTMC) formula:

```
λ_g = -log(1 - θ₁) / 3
λ_d = -log(1 - θ₀) / 3
```

### Why standard EM guarantees fail

The standard EM monotonicity guarantee requires the M-step to maximise (or at
least not decrease) the complete-data expected log-likelihood:

```
Q(Θ | Θ^old) = E_{H|Y,Θ^old}[log p(Y, H | Θ)]
```

Q decomposes additively:

```
Q = Q_prior(α, θ₁, θ₀) + Q_misclass(π) + Q_emission(σ²_g, σ²_d, λ_g, λ_d)
```

The closed-form M-step for (θ₁, θ₀) correctly maximises **Q_prior**:

```r
theta1 <- T11 / D1
theta0 <- T01 / D0
```

However, under the linked specification, updating θ also changes λ via
`ctmc_lambda_from_theta(theta)`. The **Q_emission** block (EMG density terms)
depends non-trivially on λ:

```
log f_EMG(x; λ, σ²) = log(λ) + λ²σ²/2 - λx + log(erfc((λσ² - x)/(√2 σ))) - log(2)
```

The M-step for θ does **not** account for this λ-dependence in Q_emission. When
θ changes substantially (as on early iterations), the induced change in λ can
decrease Q_emission by more than Q_prior increases, so:

```
Q(Θ^new | Θ^old) < Q(Θ^old | Θ^old)   # overall Q decreases
→ L(Θ^new) < L(Θ^old)                  # observed-data LL decreases
```

This makes the algorithm a **Generalised EM (GEM)**, not full EM. GEM converges
to a stationary point but does not guarantee monotone ascent at every iteration.

### TeX documentation error

The original TeX document claimed: *"Monotone ascent is guaranteed because every
M-step update maximizes Q in closed form."* This statement is **incorrect** for
the linked specification. It would only hold if λ were a free parameter with its
own independent M-step.

## Solution

### Immediate: tighten monitoring threshold

Change the LL-decrease warning threshold to flag real bugs while tolerating GEM
behaviour:

```r
# Was: -0.1 * abs(ll_prev)  — too loose, hides real errors on large datasets
# Now: -1e-4 * abs(ll_prev) — tight enough to catch bugs, loose enough for GEM
if (ll_change < -1e-4 * abs(ll_vec[iter - 1])) {
  warning(sprintf("EM iter %d: log-likelihood decreased by %.6e", iter, ll_change))
}
```

The empirical drops from GEM behaviour are typically `< 1e-3 × |LL|`, so
`-1e-4` catches genuine errors while silently accepting the GEM drops.

### Proper fix A: Safeguarded M-step (line search)

After the closed-form θ update, check if Q increased. If not, line-search
between old and new θ until Q does not decrease:

```r
# Pseudocode
theta_new <- closed_form_theta_update(suff)
Q_new     <- compute_Q(theta_new, gamma, df)
Q_old     <- compute_Q(theta_old, gamma, df)
if (Q_new < Q_old) {
  theta_new <- line_search(theta_old, theta_new, Q_old, gamma, df)
}
```

This restores full EM monotonicity at cost of ~1 extra Q evaluation per iteration.

### Proper fix B: Profile optimisation for θ

Replace closed-form update with numerical maximisation of
`Q_prior(θ) + Q_emission(σ², λ(θ))` jointly over (θ₁, θ₀). This is a 2D
bounded optimisation (fast: 2 parameters, smooth objective):

```r
opt <- optim(
  par = c(theta1_old, theta0_old),
  fn  = function(th) -q_full(th, suff, sigma2_g, sigma2_d),
  method = "L-BFGS-B",
  lower  = c(1e-6, 1e-6),
  upper  = c(0.999, 0.999)
)
```

### Proper fix C: Unlinked specification

Treat (λ_g, λ_d) as free parameters with their own M-step. The EMG M-step for
λ has a known closed form:

```
λ̂ = N_emg / S_emg_adj
```

where `S_emg_adj` is the weighted sum of (EMG mean – observed duration) terms.
This gives true EM with guaranteed monotonicity, at the cost of 2 extra free
parameters and loss of the theoretical CTMC link.

## Prevention

**Rule**: When using a **linked specification** where one parameter set
determines another deterministically (θ → λ, or any similar link), the M-step
must maximise Q over the full linked parameter space, not just the "source"
parameter block. Failing to do this produces a GEM, not full EM.

**Pattern to watch for**: Model has parameter Θ = (A, B) where B = f(A) is
deterministic. The Q function includes terms from both A and B. If the M-step
only maximises Q_A (terms depending directly on A), it ignores ∂Q_B/∂A which
is non-zero because B depends on A. → Always verify monotonicity empirically
with synthetic data before deploying.

**Diagnostic test**:
```r
# After any EM implementation, run on synthetic data and plot LL path
df  <- simulate_panel(n = 200, seed = 1)
fit <- em_fit_tenure(df, max_iter = 50, verbose = 0)
plot(fit$history$loglik, type = "l", main = "EM log-likelihood path")
# Should be non-decreasing after warm-up; sustained decreases = M-step bug
```

## Related

- `EM-tenure/CONCERNS.md §13` — detailed formal analysis of the GEM issue
- `EM-tenure/R/em_driver.R` — current monitoring threshold implementation
- `2026-03-10-vectorise-r-estep-over-obs.md` — performance fix applied in same session
- Wu (1983) "On the Convergence Properties of the EM Algorithm" — formal GEM convergence theory
- Dempster, Laird & Rubin (1977) — original EM paper, Section 4 on generalised EM
