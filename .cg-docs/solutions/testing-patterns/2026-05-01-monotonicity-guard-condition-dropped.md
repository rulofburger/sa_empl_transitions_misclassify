---
date: 2026-05-01
title: "Monotonicity guard condition accidentally dropped during refactor — silent always-revert bug"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, monotonicity-guard, refactor, parse-error, always-revert, convergence, code-review]
root-cause: "During a refactor the `if (estep_check$loglik < threshold)` guard condition was dropped, leaving the revert block executing unconditionally after every M-step — EM always stops at iteration 1 with converged=TRUE"
severity: "P1"
---

# Monotonicity Guard Condition Accidentally Dropped in Refactor

## Problem

After refactoring `em_driver_eps.R`, the EM algorithm appeared to "converge"
immediately at iteration 1 on every run, returning `converged = TRUE` with
a single-iteration result. The test suite still passed because the monotonicity
tolerance test (`all(diffs >= -1e-10)`) is trivially satisfied with only one
data point.

The file also produced a parse error:

```
Error: unexpected '}' at line 198
```

## Root Cause

During the refactor, the `if (estep_check$loglik < ...)` condition was
accidentally removed from the monotonicity guard block, leaving only the
inner verbose block and the revert statements:

```r
# BROKEN: guard condition missing — always reverts at iteration 1
estep_check <- e_step_eps(df, params, check_df = FALSE)
  if (verbose >= 2) {                    # ← orphaned inner if, one level too deep
    message(sprintf("reverting..."))
  }
  params    <- params_prev
  converged <- TRUE                      # ← runs unconditionally every iteration
  history[[iter]] <- unlist(params)
  break                                  # ← always breaks at iter 1
}                                        # ← unmatched: outer if {} never opened
```

The unmatched `}` at the end caused an R parse error. Without the parse error
(e.g., if the orphaned braces happened to balance), this bug would be completely
silent — `converged = TRUE` is the expected outcome, just not at iteration 1.

## Solution

Restore the guard condition wrapping the entire revert block:

```r
# CORRECT: guard condition present — only reverts when LL actually decreases
estep_check <- e_step_eps(df, params, check_df = FALSE)
if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
  if (verbose >= 2) {
    message(sprintf(
      "EM-eps iter %d: M-step decreased LL (%.6e -> %.6e); reverting.",
      iter, ll_new, estep_check$loglik
    ))
  }
  params    <- params_prev
  converged <- TRUE
  history[[iter]] <- unlist(params)
  break
}
```

## Prevention

**Test for multi-iteration convergence, not just final-state convergence.**

A test that only checks `fit$converged == TRUE` or `fit$loglik < Inf` will not
catch this bug. The following test patterns catch it:

```r
# 1. Verify EM ran more than 1 iteration on non-trivial data
test_that("em_fit_tenure_eps: runs multiple iterations before converging", {
  fit <- em_fit_tenure_eps(.make_eps_data(n = 200), max_iter = 50)
  expect_true(fit$iterations > 1,
              info = "EM converged at iteration 1 — guard may always be firing")
})

# 2. Verify LL changes across iterations (monotone + non-trivial)
test_that("em_fit_tenure_eps: LL is strictly non-decreasing across iterations", {
  fit   <- em_fit_tenure_eps(.make_eps_data(n = 200), max_iter = 50)
  diffs <- diff(fit$history$loglik)
  expect_true(all(diffs >= -1e-10))
  expect_true(length(diffs) > 0,
              info = "Need at least 2 iterations to test monotonicity")
})

# 3. Verify parameter estimates change from initial values
test_that("em_fit_tenure_eps: params change from initialization", {
  p0  <- init_params_eps(.make_eps_data())
  fit <- em_fit_tenure_eps(.make_eps_data(n = 200), params0 = p0, max_iter = 50)
  expect_false(identical(fit$params$eps, p0$eps),
               info = "eps unchanged — EM may not have run")
})
```

**Code review checklist for monotonicity guards:**
- [ ] The `if (check$loglik < threshold)` condition wraps the *entire* revert block
- [ ] The revert block is not reachable by any path that bypasses the condition
- [ ] `verbose` inner blocks are *inside* the outer `if`, not adjacent to it
- [ ] Parse check: brace depth returns to for-loop level after the guard block

## Related

- [2026-04-30-em-guard-unconditional-for-approximate-mstep.md](../testing-patterns/2026-04-30-em-guard-unconditional-for-approximate-mstep.md) — guard must be unconditional for approximate M-steps (correctness)
- [2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md) — guard caching strategy to avoid double E-step (performance)
- [2026-03-20-em-monotonicity-guard.md](2026-03-20-em-monotonicity-guard.md) — original guard diagnosis
- [2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md](2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md) — third guard failure mode: guard present and syntactically correct, but mathematically stale after emission contract change
