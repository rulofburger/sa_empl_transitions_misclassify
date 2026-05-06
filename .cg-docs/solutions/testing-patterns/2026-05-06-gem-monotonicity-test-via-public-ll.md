---
date: 2026-05-06
title: "Testing GEM monotonicity via public API instead of private Q function"
category: "testing-patterns"
language: "R"
tags: [GEM, EM-algorithm, private-function-testing, monotonicity, public-API, m-step]
root-cause: "A monotonicity test directly called the private .eval_q_miscl() function, coupling the test to an implementation detail and bypassing the public observed-data log-likelihood as the canonical correctness check"
severity: "P2"
---

# Testing GEM monotonicity via public API instead of private Q function

## Problem

The NR step in `m_step_inconsistency()` is guaranteed not to decrease the
Q-function for the misclassification block (`Q_miscl`). A test was written to
verify this by directly calling the private helper `.eval_q_miscl()`:

```r
# FRAGILE: direct call to private implementation detail
q_old <- .eval_q_miscl(params$delta,  p_mat, Y_age, Y_edu, w)
q_new <- .eval_q_miscl(out_m$delta,   p_mat, Y_age, Y_edu, w)
expect_true(q_new >= q_old - 1e-10)
```

**Problems with this approach**:
1. `.eval_q_miscl` is a private (`.`-prefixed) function with no stability
   guarantee. If it is renamed, inlined, or removed during a refactor, the test
   breaks for the wrong reason.
2. The test does not verify that the *observed-data log-likelihood* improves —
   only that an internal Q proxy does. A bug that made Q increase while LL
   decreased would pass undetected.
3. Private R functions are in `.GlobalEnv` due to sourcing, so they are
   technically accessible from tests — but this is an implementation accident,
   not a contract.

## Root Cause

The GEM guarantee (Q non-decreasing ⇒ LL non-decreasing) was tested at the
wrong level of abstraction. The private Q function was used because it was
convenient and because the public LL function is more expensive to evaluate.

## Solution

Test the GEM guarantee using two consecutive E-step log-likelihoods evaluated
under the public API:

```r
test_that("m_step_inconsistency: NR step does not decrease LL (public API)", {
  df      <- .make_incons_panel(n = 200L)
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()

  # E-step at initial params
  out_e1  <- e_step_inconsistency(df, inc_mat, params)
  # M-step produces updated params
  params1 <- m_step_inconsistency(out_e1$suff, inc_mat, params)
  # E-step at updated params — this evaluates the observed-data LL
  out_e2  <- e_step_inconsistency(df, inc_mat, params1)

  expect_true(out_e2$loglik >= out_e1$loglik - 1e-8,
              label = "GEM iteration must not decrease observed-data LL")
})
```

**Why this works**:
- The GEM guarantee `Q(θ_new | θ_old) ≥ Q(θ_old | θ_old)` implies
  `LL(θ_new) ≥ LL(θ_old)` under regularity conditions (Dempster et al. 1977).
- Testing the LL directly verifies the guarantee as observed by the caller,
  not by an internal proxy.
- The test uses only exported functions and is robust to any internal refactoring
  of `.eval_q_miscl`.

## Prevention

- **Rule**: Tests must not call `.`-prefixed (private) functions. If you need
  to verify an internal property, find its observable consequence in the
  public API.
- **GEM/EM monotonicity pattern**: always test via two consecutive E-step
  log-likelihoods rather than Q values. This gives a stronger guarantee
  (full LL, not a partial Q block) and uses only the public interface.
- **Sign that a test is too deep**: if a test imports variables by name from
  `suff` and passes them directly to a private helper, it is testing
  implementation rather than contract.

## Related

- [2026-03-20-em-monotonicity-guard.md](2026-03-20-em-monotonicity-guard.md)
- [2026-04-30-em-guard-unconditional-for-approximate-mstep.md](2026-04-30-em-guard-unconditional-for-approximate-mstep.md)
- [2026-05-01-monotonicity-guard-condition-dropped.md](2026-05-01-monotonicity-guard-condition-dropped.md)
