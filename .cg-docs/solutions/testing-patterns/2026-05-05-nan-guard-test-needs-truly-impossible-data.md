---
date: 2026-05-05
title: "NaN guard tests must use truly impossible data, not low-probability data"
category: "testing-patterns"
language: "R"
tags: [nan-guard, estep, em-algorithm, numerical, testing, log-sum-exp]
root-cause: "Low-probability observations still produce finite log_joint values via floating-point arithmetic; only genuinely impossible observations (no valid history) produce all-Inf rows"
severity: "P1"
---

# NaN guard tests must use truly impossible data, not low-probability data

## Problem

A test designed to trigger the EM E-step NaN guard failed silently. The test used
an extreme parameter combination (`alpha=0.001, theta1=0.999`) expecting the observation
`(1,1,0)` to be numerically impossible. Instead, `e_step()` computed a valid (if
extremely skewed) responsibility matrix and returned successfully.

```r
# ❌ This did NOT trigger the NaN guard — (1,1,0) is compatible with several
# histories and has small but nonzero probability
df     <- data.frame(y1 = 1L, y2 = 1L, y3 = 0L, weight = 1)
params <- list(alpha = 0.001, theta0 = 0.001, theta1 = 0.999)
expect_error(
  e_step(df, params, model_type = "none", validate = FALSE),
  regexp = "zero probability"
)
# Result: test passed silently (no error thrown), test failed
```

## Root Cause

Floating-point log-probability arithmetic keeps "nearly impossible" observations
alive. As long as at least one of the 8 latent histories produces a finite
`log_joint` value, `logsumexp` normalises successfully and the NaN guard never
fires. The guard fires only when **every** history produces `-Inf` log_joint —
i.e., when the observation is logically incompatible with the model's support.

With `model_type = "none"` (exact-match misclassification), history `(1,1,0)` has
a small but nonzero prior under `alpha=0.001`:

```
log(prior for h=(1,1,0)) = log(0.001) + log(0.001) + log(0.999) ≈ -13.8  (finite!)
```

So `log_joint[h=(1,1,0)]` is finite and the guard does not fire.

## Solution

Use a value outside the model's support — a non-binary integer (`y1 = 2L`) that
cannot match any of the 8 binary histories `{0,1}^3`. This guarantees every row of
`log_joint` is `-Inf`, making `logsumexp` return `NaN` and triggering the guard.

```r
# ✅ y1=2 is outside the binary support of model_type='none'
# All 8 histories have h1 ∈ {0,1}; none matches y1=2 → log_joint all -Inf
df     <- data.frame(y1 = 2L, y2 = 0L, y3 = 0L, weight = 1)
params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
expect_error(
  e_step(df, params, model_type = "none", validate = FALSE),
  regexp = "zero probability"
)
```

Note: `validate = FALSE` is required to bypass the upstream binary check, which
would fire first and mask the NaN guard. The goal of this test is specifically to
verify the NaN guard, not the input validation layer.

## Prevention

When testing numeric guards (NaN, Inf, overflow):

- **Use structurally impossible data**, not probabilistically extreme data.
- **Disable upstream guards** (`validate = FALSE`) to isolate the specific guard
  being tested.
- **Document why** `validate = FALSE` is intentional in the test comment.
- **Separate guard tests by layer**: one test for the validation guard, one for the
  NaN guard — they should test independent code paths.

### Checklist for NaN guard tests

```r
# ✅ Correct: obs is outside model support (all log_joint = -Inf)
df <- data.frame(y1 = 2L, ...)  # non-binary value

# ❌ Wrong: obs is just unlikely (some log_joint values are tiny but finite)
df <- data.frame(y1 = 1L, ...)
params <- list(theta1 = 0.999)   # just extreme params
```

## Related

- [2026-03-20-em-monotonicity-guard.md](../testing-patterns/2026-03-20-em-monotonicity-guard.md) — related guard testing in the EM driver
- [2026-05-01-monotonicity-guard-condition-dropped.md](../testing-patterns/2026-05-01-monotonicity-guard-condition-dropped.md) — guard condition logic pitfall
