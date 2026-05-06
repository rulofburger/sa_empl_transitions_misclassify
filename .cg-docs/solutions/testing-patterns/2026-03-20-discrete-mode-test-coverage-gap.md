---
date: 2026-03-20
title: "Closing discrete-mode test coverage gaps in EM estimation"
category: "testing-patterns"
language: "R"
tags: [testthat, em-algorithm, discrete-mode, integration-tests, unit-tests, coverage]
root-cause: "All e_step/m_step and em_fit_tenure tests used discrete_timegap=FALSE; the new discrete path had zero test coverage."
severity: "P2"
---

# Closing Discrete-Mode Test Coverage Gaps in EM Estimation

## Problem

After implementing a new `discrete_timegap = TRUE` code path across `estep.R`,
`mstep.R`, `emissions.R`, and `em_driver.R`, the full test suite (247 tests at the
time) only exercised `discrete_timegap = FALSE`. The new path was untested at both
the unit-test and integration-test level:

- `test-estep-mstep.R` — all 8 tests hard-coded `discrete_timegap = FALSE`
- `test-integration.R` — all 6 `em_fit_tenure` tests used the legacy continuous path

Symptoms: high confidence in CI "green" but silent regressions possible in the default
code path (`discrete_timegap = TRUE` is the project default).

## Root Cause

Test suite was written alongside the original (continuous-timegap) implementation and
never updated when the discrete path was added. The two modes share the same function
signatures, so calls with the old default `discrete_timegap = FALSE` continued to pass
without ever exercising the new code.

## Solution

### Unit tests added to `test-estep-mstep.R`

Four `e_step` tests for discrete mode:
1. **Structure** — `gamma` is `n × 8`, rows sum to 1, `loglik` is finite.
2. **Suff stats content** — `cat_d_marginal_c`, `cat_d_marginal_w` present; `Sd`/`Nd`
   absent (those are continuous-mode sufficient stats).
3. **Perfect-classification data** — with `pi = 0` the diagonal histories dominate.
4. **NA rejection** — a single `NA` in `tenure1` triggers an informative `stop()`.

Four `m_step` tests for discrete mode:
1. **Required params / no sigma2_d** — all of `alpha, theta1, theta0, pi, sigma2_g,
   lambda_g` present; `sigma2_d` is `NULL`.
2. **No-misclassification** — `pi = 0` when `misclassification = FALSE`.
3. **Stationarity constraint** — `alpha == theta0 / (theta0 + 1 - theta1)`.
4. **Lambda-theta consistency** — `lambda_g` equals `ctmc_lambda_from_theta(theta1)`.

### Integration tests added to `test-integration.R`

Six `em_fit_tenure` tests for discrete mode:
1. Convergence on synthetic data (parameter recovery within tolerance).
2. Monotone log-likelihood (no decrease > 1e-6 × |LL|).
3. No-misclassification model (`misclassification = FALSE`).
4. Stationarity mode convergence + alpha constraint satisfied + monotone LL.
5. Gamma dimensions (n × 8, row-sums = 1).
6. Legacy stationarity with monotonicity guard (regression test for P1 fix).

### Pattern

```r
# Always add a parallel discrete-mode block when testing a function
# that has a discrete_timegap argument.

# --- Legacy (continuous) ---
test_that("e_step foo: legacy", {
  df <- simulate_panel(n = 80, seed = 1, discrete_timegap = FALSE)
  params <- init_params(df, discrete_timegap = FALSE)
  out <- e_step(df, params, discrete_timegap = FALSE)
  # ... assertions ...
})

# --- Discrete (default) ---
test_that("e_step foo: discrete", {
  df <- simulate_panel(n = 80, seed = 1, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  out <- e_step(df, params, discrete_timegap = TRUE)
  # ... same assertions ...
})
```

## Prevention

- When adding a new `mode` / `method` / `engine` argument to an existing function,
  **immediately** add parallel tests for every mode, not just the new one.
- Prefer a helper `make_test_df(discrete = TRUE/FALSE)` to reduce boilerplate.
- In code review, flag any test file where all calls share the same value of a
  boolean-dispatch argument — that is almost always a coverage gap.

## Related

- `.cg-docs/solutions/testing-patterns/2026-03-20-em-monotonicity-guard.md` — the P1
  fix that the new stationarity integration test now guards against.
