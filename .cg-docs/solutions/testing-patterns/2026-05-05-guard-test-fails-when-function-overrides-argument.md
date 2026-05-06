---
date: 2026-05-05
title: "Testing guards on functions that pass a fixed internal argument tests nothing"
category: "testing-patterns"
language: "R"
tags: [guard-testing, hardcoded-argument, null-argument, inference, em-algorithm, degenerate-prior, test-design]
root-cause: "A function that unconditionally passes a fixed value (e.g. `alpha=NULL`) to a lower-level helper makes any test of the *guard triggered by non-NULL alpha* untestable — the guard path is structurally unreachable."
severity: "P2"
---

# Testing Guards on Functions That Override an Argument Internally

## Problem

`implied_transitions_ar2()` contained a comment:
```r
alpha = NULL  # force stationary distribution for inference
```

Two tests were written to verify that the function raises errors for degenerate
`alpha` values:

```r
test_that("implied_transitions_ar2 errors on degenerate prior (no mass on h1=0)", {
  degenerate_params <- list(
    theta0=0.10, theta01=0.15, theta1=0.08, theta10=0.12, pi=0,
    alpha = c(`00`=0, `10`=0, `01`=0.5, `11`=0.5)
  )
  expect_error(implied_transitions_ar2(degenerate_params), regexp="no prior mass")
})
```

Both tests failed: no error was thrown. The `params$alpha` was never passed to
`prior_over_histories_ar2()` — the function unconditionally passed `NULL`.

## Root Cause

The function hardcoded `alpha = NULL` regardless of `params$alpha`:

```r
prior <- prior_over_histories_ar2(hmat, ..., alpha = NULL)
```

The degenerate guard in `prior_over_histories_ar2` requires a non-NULL alpha.
Since `alpha=NULL` always flows through, the guard is **structurally unreachable**
via `implied_transitions_ar2`. The tests correctly described the expected behaviour
but the implementation never exercised the code path.

There was an additional error: the `alpha` names did not match the convention. The
AR(2) parameterisation names pairs as `"XY"` where X=h1, Y=h2. So h1=0 states are
`"00"` and `"01"` (not `"00"` and `"10"`).

## Solution

### Fix the function

Thread `params$alpha` through to the helper, falling back to stationary when
`params$alpha` is NULL:

```r
implied_transitions_ar2 <- function(params) {
  hmat  <- latent_histories_ar2()
  # Use params$alpha if provided (e.g. during EM diagnostics).
  # For standard post-estimation inference, params$alpha = NULL → stationary.
  prior <- prior_over_histories_ar2(hmat,
    theta0  = params$theta0,
    theta01 = params$theta01,
    theta1  = params$theta1,
    theta10 = params$theta10,
    alpha   = params$alpha  # NULL → stationary; non-NULL → use as-is
  )
  ...
}
```

### Fix the tests — get the alpha names right

The alpha naming convention for AR(2) is `"XY"` where X=h1, Y=h2:

| Name | Meaning |
|------|---------|
| `"00"` | h1=0, h2=0 |
| `"01"` | h1=0, h2=1 |
| `"10"` | h1=1, h2=0 |
| `"11"` | h1=1, h2=1 |

"h1=0 states" are `"00"` and `"01"`, **not** `"00"` and `"10"`. Setting `"00"` and
`"10"` to zero actually removes mass from one h1=0 state and one h1=1 state, so the
guard may or may not fire depending on the remaining weights:

```r
# ❌ Wrong: removes "00" (h1=0) and "10" (h1=1)
alpha = c(`00`=0, `10`=0, `01`=0.5, `11`=0.5)

# ✅ Correct: removes both h1=0 states
alpha = c(`00`=0, `01`=0, `10`=0.5, `11`=0.5)
```

## Prevention

1. **Before writing a guard test, trace how the argument flows.** If the public
   function unconditionally overrides the argument before passing it to the internal
   function, the guard is unreachable and the test will always pass vacuously.

2. **Fix the function first, then write the test.** If the guard tests a path that
   only makes sense when the argument is user-controlled, make the function accept it:

   > Public function: thread `params$alpha` through.  
   > Tests: supply degenerate alpha in `params`.

3. **Verify alpha name conventions match the helper's expectation.** Print or
   inspect the expected names with `names(stationary_ar2(...))` before constructing
   test fixtures.

4. **Run tests immediately after adding them** — don't batch guard tests with other
   changes. A guard test that passes without triggering the error is a false positive.

## Related

- `EM-AR2/R/inference.R` — `implied_transitions_ar2` (fixed)
- `EM-AR2/tests/testthat/test-inference.R` — the tests (fixed)
- [2026-05-05-nan-guard-test-needs-truly-impossible-data.md](2026-05-05-nan-guard-test-needs-truly-impossible-data.md) — related false-positive guard test pattern
