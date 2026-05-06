---
date: 2026-05-05
title: "any(x <= 0) silently returns NA when x contains NA, bypassing the guard"
category: "data-quality"
language: "R"
tags: [na-handling, input-validation, weights, guard, any, logical-na, edge-case]
root-cause: "`any(NA <= 0)` returns `NA` (not TRUE), so a validation guard `if (any(w <= 0)) stop(...)` silently passes when weights include NA values."
severity: "P1"
---

# `any(x <= 0)` Returns NA When x Contains NA, Bypassing the Guard

## Problem

The EM driver had an input validation guard for non-positive weights:

```r
if (any(w <= 0)) stop("all weights must be positive")
```

When `w` contained `NA`, the guard silently passed (no error thrown), allowing the
EM loop to proceed. The algorithm then produced `NaN` log-likelihood values without
a diagnostic message.

## Root Cause

In R, `NA <= 0` evaluates to `NA`, not `TRUE`. And `any(c(1, NA, 2) <= 0)` returns
`NA` (not `TRUE`), because R cannot rule out that the `NA` element might be ≤ 0.
The `if (NA)` branch is not taken.

```r
any(c(1, NA, 2) <= 0)   # NA
if (any(c(1, NA, 2) <= 0)) cat("triggered")  # silent — no output
```

This means the positivity guard provides **zero protection** when the weight vector
contains `NA`.

## Solution

Add a separate, explicit `any(is.na(w))` check **before** the positivity check:

```r
if (any(is.na(w)))   stop("em_fit_ar2: weights must not contain NA values")
if (any(w <= 0))     stop("em_fit_ar2: all weights must be strictly positive")
```

The NA check comes first so the positivity guard can safely assume no NA values.
Do **not** use `any(w <= 0, na.rm = TRUE)` — that silently treats NA as non-positive
and converts what should be a hard error into a silent ignore.

## Prevention

- For every `any(x < threshold)` or `any(x == value)` guard in an input validation
  block, ask: "can x contain NA?" If yes, add `any(is.na(x))` first.
- Never use `na.rm = TRUE` inside a validation guard — that turns errors into silence.
- Pattern to follow for input validation of numeric vectors:

```r
# ✅ Correct guard sequence
if (any(is.na(x)))   stop("x must not be NA")
if (any(x <= 0))     stop("x must be positive")
```

## Related

- `EM-AR2/R/em_driver.R` — the file where this fix was applied
- [2026-04-30-na-structural-missingness-vs-error.md](../data-quality/2026-04-30-na-structural-missingness-vs-error.md) — related NA-handling patterns
