---
date: 2026-05-05
title: "Baum-Welch free-alpha restores EM monotonicity for short AR(2) panels"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, hmm, baum-welch, ar2, stationary-distribution, monotonicity, initial-distribution, free-alpha, short-panel]
root-cause: "Forcing the initial pair distribution alpha(h1,h2) to equal the stationary distribution of the current theta estimates at every M-step is an *approximation*. The stationarity constraint is not maintained by the closed-form M-step, so Q can decrease."
severity: "P1"
---

# Baum-Welch Free-Alpha Restores EM Monotonicity for Short AR(2) Panels

## Problem

The AR(2) EM driver (`EM-AR2/R/em_driver.R`) produced a monotonically *decreasing*
log-likelihood across iterations. The test:

```r
test_that("em_fit_ar2 log-likelihood is non-decreasing across iterations", {
  diffs <- diff(fit$history$loglik)
  expect_true(all(diffs > -1e-4))
})
```

...failed with LL decreases of several nats per iteration from the very first step.

## Root Cause

The standard "initialise α from the stationary distribution" approach for HMMs on
long sequences implicitly assumes that the stationarity constraint is maintained
throughout EM. On short (4-wave) panels this approximation breaks:

1. The M-step updates θ (transition parameters) via `p̂_{jk→1} = T_{jk}/D_{jk}`.
2. The new θ has a *different* stationary distribution than α was computed from.
3. The prior `prior_over_histories_ar2(alpha=stationary(θ_new))` evaluated at the
   *new* θ can assign very different weights to histories, decreasing the observed-data
   LL even though the M-step maximised Q.

In short, the constraint `alpha = stationary(theta)` is not maintained by the
unconstrained closed-form M-step. This violates the EM ascent guarantee.

The TeX documentation noted this explicitly:
> "we adopt the standard approximation for short panels: update the transition
> parameters from the transition counts alone, then set the initial distribution by
> stationarity."

This is an **approximation** — it does NOT guarantee monotone LL.

## Solution

Use the **Baum-Welch free-alpha** approach: treat `alpha(h1,h2)` as a free
length-4 parameter updated from data via sufficient statistics at each M-step.

### E-step change

Add a 2×2 matrix `S_{jk}` to the sufficient statistics:

```r
# S_{jk} = sum_i sum_h gamma_{ih} * w_i * 1(h1=j, h2=k)
# = weighted responsibility of all histories starting with (h1=j, h2=k)
S <- matrix(0, 2, 2, dimnames=list(c("0","1"), c("0","1")))
for (j in 0:1) for (k in 0:1) {
  ind <- (hmat[,1] == j) & (hmat[,2] == k)
  S[j+1, k+1] <- sum(col_wgamma * ind)  # col_wgamma = colSums(wgamma)
}
```

### M-step change

Update alpha directly from `S`:

```r
alpha_hat <- S / sum(S)   # 2x2 matrix → normalise → use as named 4-vector
```

No call to `stationary_ar2()` during EM iterations.

### Driver change

Include `alpha` in `params` from the start, initialised from stationary:

```r
# init_params_ar2: one-time init from stationary
alpha <- tryCatch(
  stationary_ar2(theta0, theta01, theta1, theta10),
  error = function(e) c(`00`=0.25, `10`=0.25, `01`=0.25, `11`=0.25)
)
```

Pass `params$alpha` to `prior_over_histories_ar2()` during every EM iteration.
After convergence, recompute the stationary distribution from the converged `theta`
for **inference** only (not during iterations).

### Why this works

The free-alpha M-step update `alpha_hat = S/sum(S)` is the exact maximiser of the
Q-function term for alpha (it maximises the expected log-prior over initial pairs).
Since both the alpha update and the theta update each maximise their respective Q
terms, the full M-step is the standard Baum-Welch M-step and guarantees
`Q(θ_new, α_new | θ_old, α_old) ≥ Q(θ_old, α_old | θ_old, α_old)`, which implies
non-decreasing LL.

## Inference After Convergence

After EM converges, the estimated `alpha` is the sample-based initial pair
distribution (not necessarily stationary). For publishable implied transitions and
GoF, call `implied_transitions_ar2(params)` with `alpha=NULL` to use the stationary
distribution derived from the converged theta — this gives model-implied long-run
quantities independent of the particular sample composition.

```r
# During EM:
prior <- prior_over_histories_ar2(hmat, ..., alpha = params$alpha)  # free alpha

# After convergence, for inference:
trans <- implied_transitions_ar2(params_without_alpha)  # uses stationary
```

## Prevention

- **Never** force `alpha = stationary(theta)` at every M-step for short panels with
  a closed-form transition M-step. Use free-alpha Baum-Welch instead.
- If you want to report results under stationarity, compute the stationary
  distribution from the **converged** theta **once after EM**, not during iterations.
- Always test LL monotonicity with `diffs > -1e-4` tolerance on ≥ 50 iterations.
- The stationary-approximation approach from the TeX is a reasonable starting point
  for initialisation but not for the iterative updates.

## Related

- [2026-03-20-em-monotonicity-guard.md](../testing-patterns/2026-03-20-em-monotonicity-guard.md) — monotonicity guard for constrained (stationary=TRUE) EM
- [2026-04-30-em-guard-unconditional-for-approximate-mstep.md](../testing-patterns/2026-04-30-em-guard-unconditional-for-approximate-mstep.md) — guard for approximate M-step
- `EM-AR2/R/estep.R` — S matrix computation
- `EM-AR2/R/mstep.R` — alpha_hat = S/sum(S)
- `EM-AR2/R/em_driver.R` — params$alpha threaded through iterations
