---
date: 2026-05-05
title: "pmax(0, inv_logit(x)) is dead code — inv_logit always returns (0,1)"
category: "data-quality"
language: "R"
tags: [inv_logit, logistic, dead-code, numerical, clamp, lower-bound, underflow]
root-cause: "inv_logit (logistic function) is strictly positive by definition; pmax(0, ...) can never activate. The intended lower bound (preventing IEEE underflow) requires max(eps, ...) instead."
severity: "P1"
---

# `pmax(0, inv_logit(x))` is dead code — `inv_logit` always returns (0,1)

## Problem

A parameter clamping helper used `pmax(0, ...)` as a lower bound on an
`inv_logit`-transformed value, intending to prevent the clamped parameter
from reaching zero:

```r
# ❌ Dead code — pmax(0, inv_logit(x)) is always equal to inv_logit(x)
.clamp_pi <- function(p, cap = PI_CAP, margin = 0.01) {
  min(cap - margin, pmax(0, .perturb_param(p)))
}
# .perturb_param() returns inv_logit(x + noise), which is always in (0,1)
```

The bug is silent: the function runs without error, returns a plausible value,
but provides **no protection against IEEE underflow** — a sufficiently small
perturbed `p` could produce a value like `1e-320` that causes confusing `NaN`
downstream rather than a clear error at the clamping site.

## Root Cause

The logistic (inverse logit) function $\sigma(x) = \frac{1}{1 + e^{-x}}$ is
**strictly positive for all finite inputs**:

$$\sigma(x) > 0 \quad \forall x \in \mathbb{R}$$

Therefore `pmax(0, sigma(x)) == sigma(x)` always — the comparison never
activates. The developer intended a lower bound to guard against underflow, but
used a bound that the function can never violate.

## Solution

Replace `pmax(0, ...)` with `max(eps, ...)` where `eps` is a small positive
constant (e.g., `1e-6`) that sits above IEEE underflow territory:

```r
# ✅ eps provides a genuine lower bound above IEEE underflow
.clamp_pi <- function(p, cap = PI_CAP, margin = 0.01, eps = 1e-6) {
  min(cap - margin, max(eps, .perturb_param(p)))
}
```

Why `1e-6`? IEEE double precision underflows below ~`5e-324`. The `eps = 1e-6`
bound sits well above this, catching pathological perturbations while still
allowing the parameter to be very small. Adjust if the model's numerical context
requires a different floor.

Also note: `pmax` is a vectorised function; `max` is correct here since `p` is
scalar. Using `pmax` on a scalar is not wrong, but it signals intent incorrectly.

## Prevention

**General rule**: When you want a lower bound on a value derived from `inv_logit`,
`exp`, or any other function that is already bounded below by zero, you must use
`max(eps, ...)` (not `max(0, ...)`), because:

- `max(0, positive_value) == positive_value` — the 0 floor is never reached
- The real danger is underflow to subnormal territory, which needs `eps > 0`

```r
# Anti-patterns to avoid
pmax(0, inv_logit(x))      # ❌ dead code
max(0,  inv_logit(x))      # ❌ dead code (0 is never triggered)
max(0,  exp(x))            # ❌ dead code (exp is always positive)

# Correct patterns
max(1e-6, inv_logit(x))    # ✅ genuine lower bound above underflow
max(1e-6, exp(x))          # ✅ genuine lower bound above underflow
pmin(1 - 1e-6, inv_logit(x))  # ✅ genuine upper bound below 1
```

**Code review heuristic**: any `pmax(0, f(x))` or `max(0, f(x))` where `f` maps
to positive reals is a dead-code smell. Flag it in review.

## Related

- [2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md) — related numerical guard patterns in the EM pipeline
