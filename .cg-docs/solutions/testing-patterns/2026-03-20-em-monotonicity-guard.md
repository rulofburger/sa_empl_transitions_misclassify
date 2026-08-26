---
date: 2026-03-20
title: "EM ascent property violation under constrained M-step"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, monotonicity, stationarity, constrained-mstep, convergence, log-likelihood]
root-cause: "Stationarity constraint forces alpha = theta0/(theta0+1-theta1), which can overshoot the unconstrained maximum and decrease LL."
severity: "P1"
---

# EM Ascent Property Violation Under Constrained M-step

## Problem

The EM algorithm guarantees that the log-likelihood is non-decreasing at every
iteration (`Q(θ_new | θ_old) ≥ Q(θ_old | θ_old)`). When `stationary = TRUE` the
M-step imposes the extra constraint:

```
alpha = theta0 / (theta0 + 1 - theta1)
```

This constraint can cause the M-step to move *away* from the unconstrained maximum,
violating the ascent property. Observed symptoms:

- Test suite emitted **9 warnings** of the form `"EM iter N: LL decreased by X"`.
- Worst case: stationarity test, iter 4, decrease of **−0.108** (not just numerical
  noise).
- After the decrease the algorithm still reported `converged = TRUE`, meaning results
  were silently at a sub-optimal constrained point.

## Root Cause

The standard EM M-step is derived under the assumption that each parameter is
maximised independently. When a stationarity constraint ties `alpha` to `theta0` and
`theta1`, the joint update no longer guarantees an increase in the observed-data
log-likelihood. The constrained alpha can be higher or lower than the unconstrained
optimum, and if it is sufficiently different it causes a net LL decrease.

## Solution

Add a **monotonicity guard** immediately after the M-step in `em_driver.R`:

```r
# --- Monotonicity guard ---
# When the M-step uses constrained updates (e.g., stationary alpha), the
# parameter update can violate the EM ascent property. We check by
# evaluating LL at the new params. If LL decreased, revert to previous
# params and declare convergence (the algorithm is at a constrained
# stationary point).
estep_check <- e_step(df, params, discrete_timegap = discrete_timegap)
if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
  if (verbose >= 2) {
    message(sprintf(
      "EM iter %d: M-step decreased LL (%.6e -> %.6e); reverting to previous params.",
      iter, ll_new, estep_check$loglik
    ))
  }
  params <- params_prev
  converged <- TRUE
  history[[iter]] <- unlist(params)
  break
}
```

**Effect**: Warnings dropped from **9 → 0**. The algorithm now treats a LL decrease as
a signal that it has reached the constrained stationary point and stops cleanly.

### Cost

One extra E-step evaluation per iteration. In practice the guard triggers only near
convergence (typically 0–2 extra evaluations per run), so the runtime impact is
negligible.

### Tolerance choice

`1e-8 * abs(ll_new)` is a relative tolerance. It ignores decreases that are purely
floating-point noise (typically < 1e-12 in absolute terms) while catching genuine
constraint-induced overshoots (≥ 1e-4 in absolute terms in the cases observed here).

## Prevention

- Any time a constrained M-step is added (stationarity, simplex constraints, etc.),
  immediately add a post-M-step monotonicity guard.
- Write an integration test that explicitly asserts `all(diff(ll) >= -tol)` for the
  constrained mode:

```r
test_that("em_fit_tenure stationary: monotone LL", {
  df <- simulate_panel(n = 200, seed = 88, discrete_timegap = TRUE)
  fit <- em_fit_tenure(df, stationary = TRUE, max_iter = 200,
                       verbose = 0, discrete_timegap = TRUE)
  ll <- fit$history$loglik[!is.na(fit$history$loglik)]
  expect_true(all(diff(ll) >= -1e-6 * max(abs(ll))))
})
```

- Treat **any warning about LL decrease** as a P1 bug, not acceptable noise.

## Related

- `.cg-docs/solutions/testing-patterns/2026-03-20-discrete-mode-test-coverage-gap.md`
  — the coverage work that produced the stationarity integration test which catches
  this pattern.
- [2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md)
  — **performance follow-up**: guard was unnecessarily running on every unconstrained
  iteration, doubling E-step cost. Fixed by gating on `stationary || linked`.
- [EM algorithm convergence theory (Dempster, Laird & Rubin 1977)](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x)
