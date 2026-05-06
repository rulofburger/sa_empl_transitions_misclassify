---
date: 2026-05-05
title: "validate=FALSE pattern: skip per-call validation in iterative EM loops"
category: "performance-issues"
language: "R"
tags: [em-algorithm, validation, hot-path, performance, estep, design-pattern]
root-cause: "Running full data validation on every EM iteration is O(N × iterations) work; data does not change between iterations, so validation need only run once."
severity: "P2"
---

# `validate=FALSE` pattern: skip per-call validation in iterative EM loops

## Problem

Iterative algorithms like EM call the same inner function hundreds of times
with identical data. If that inner function validates its data inputs on every
call, the validation overhead accumulates across all iterations and becomes
a significant fraction of total runtime — especially for large N.

Example (before fix):

```r
# em_fit_baseline() — old pattern
for (iter in seq_len(max_iter)) {
  estep_out <- e_step(df, params, model_type)   # validates df every iteration
  params    <- m_step(estep_out$suff, ...)
}
```

With `validate = TRUE` (the default), `e_step` runs binary checks, weight
checks, and column existence checks on `df` on every iteration. `df` is
immutable across iterations.

## Root Cause

The original `e_step()` had no `validate` parameter — validation was
unconditional. Adding the parameter creates a toggle, but the pattern only
becomes correct when the caller (the EM driver) understands the invariant:
*data quality is guaranteed once before the loop starts*.

## Solution

Add a `validate` parameter to the inner function (default `TRUE`) and
have the driver call it with `validate = FALSE` after an explicit up-front check:

```r
# e_step() signature (estep.R)
e_step <- function(df, params, model_type = "none", validate = TRUE) {
  if (validate) {
    # Full data quality checks (binary, weights, columns, params)
    .check_df_binary(df)
    .check_weights(df)
    for (nm in c("alpha", "theta0", "theta1")) {
      v <- params[[nm]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 || v >= 1)
        stop(sprintf("params$%s must be a single finite numeric in (0, 1); got: %s",
                     nm, deparse(v)))
    }
  }
  # ... rest of e_step ...
}

# em_fit_baseline() — correct pattern (em_driver.R)
e_step(df, params0, model_type, validate = TRUE)   # warm-up: validate once

for (iter in seq_len(max_iter)) {
  estep_out <- e_step(df, params, model_type, validate = FALSE)  # skip validation
  params    <- m_step(estep_out$suff, ...)
}
# Final E-step after convergence also uses validate = FALSE
final_estep <- e_step(df, params, model_type, validate = FALSE)
```

**Key rule**: the driver calls `e_step(..., validate = TRUE)` exactly once — for the
warm-up / initial E-step before the loop — and `validate = FALSE` on every
subsequent call. This is safe because `df` is never modified inside the loop.

## Prevention

Apply this pattern whenever:
1. A function validates its inputs
2. That function is called inside a loop
3. The validated inputs are immutable across loop iterations

```r
# Pattern template
outer_function <- function(data, ...) {
  inner_function(data, ..., validate = TRUE)  # once, before loop
  for (i in seq_len(n_iter)) {
    result <- inner_function(data, ..., validate = FALSE)  # every iteration
  }
}
```

**Roxygen documentation convention** — document both sides:

In the inner function:
```r
#' @param validate Logical (default TRUE). If FALSE, skips per-call data
#'   checks. Only set to FALSE when data has been validated once upstream.
```

In the driver `@details`:
```r
#' @details **Validation strategy**: data quality is validated *once* before
#'   entering the EM loop. \code{e_step} is then called with
#'   \code{validate = FALSE} on every iteration and on the final E-step to
#'   avoid redundant checks on immutable data.
```

### Tests to write for this pattern

```r
# 1. validate=FALSE produces same result as validate=TRUE on valid data
test_that("e_step validate=FALSE produces same result as validate=TRUE", {
  r1 <- e_step(df, params, model_type = "symmetric", validate = TRUE)
  r2 <- e_step(df, params, model_type = "symmetric", validate = FALSE)
  expect_equal(r1$loglik, r2$loglik, tolerance = 1e-12)
})

# 2. validate=FALSE skips the data check (invalid data reaches the next guard)
test_that("e_step validate=FALSE skips data check", {
  # validate=TRUE: binary check fires first
  expect_error(e_step(bad_df, params, validate = TRUE), regexp = "binary")
  # validate=FALSE: binary check skipped; next guard fires instead
  expect_error(e_step(bad_df, params, validate = FALSE), regexp = "zero probability")
})

# 3. validate=TRUE catches bad params
test_that("e_step validate=TRUE catches invalid params", {
  params_bad <- list(alpha = 2.0, ...)
  expect_error(e_step(df, params_bad, validate = TRUE), regexp = "params\\$alpha")
})
```

## Related

- [2026-03-15-vectorise-em-estep-over-observations.md](../performance-issues/2026-03-15-vectorise-em-estep-over-observations.md) — vectorised E-step (related hot-path optimization)
- [2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md) — gating expensive guards in the EM loop
- [2026-05-05-nan-guard-test-needs-truly-impossible-data.md](../testing-patterns/2026-05-05-nan-guard-test-needs-truly-impossible-data.md) — how to test the NaN guard that fires when validate=FALSE
