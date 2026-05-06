---
date: 2026-03-15
title: "EM converged flag wrong when loglik_trace is trimmed before comparison"
category: "data-quality"
language: "R"
tags: [em-algorithm, convergence, flag, bug, length-comparison, trimming]
root-cause: "converged = length(loglik_trace) < max_iter evaluated AFTER the trace was trimmed to the actual iteration count, so the condition was always FALSE (trimmed length is always < max_iter whether or not convergence occurred)"
severity: "P2"
---

# EM Convergence Flag Wrong After Trace Trimming

## Problem

`em_fit()` returned `converged = FALSE` for runs that had genuinely converged.
The `converged` field in the output list was always `FALSE` regardless of
whether the relative-change tolerance was satisfied.

## Root Cause

The driver pre-allocated `loglik_trace <- numeric(max_iter)` and then, upon
convergence, trimmed it:

```r
if (rel_change < tol) {
  loglik_trace <- loglik_trace[seq_len(it)]  # trim here
  break
}

# ... later:
converged = length(loglik_trace) < max_iter  # WRONG
```

After trimming, `length(loglik_trace)` is the convergence iteration count
(e.g. 47), which is **always** less than `max_iter` (500).  So the expression
evaluates to `TRUE` when converged *and* when not converged (the non-convergent
path also trims: `loglik_trace <- loglik_trace[seq_len(it)]` at `it ==
max_iter`).

Even if the non-convergent path did not trim, `length < max_iter` tests the
wrong thing: it conflates "we trimmed" with "we converged".

## Solution

Introduce an explicit boolean flag and set it inside the convergence block:

```r
converged <- FALSE   # initialise before loop

for (it in seq_len(max_iter)) {
  ...
  if (rel_change < tol) {
    loglik_trace <- loglik_trace[seq_len(it)]
    converged <- TRUE           # set explicitly
    if (verbose) message(...)
    break
  }
}

# output:
list(
  ...
  converged = converged,   # not length(loglik_trace) < max_iter
  ...
)
```

## Prevention

- **Never infer convergence from trace length.**  Always use an explicit
  flag initialised to `FALSE` and set to `TRUE` only inside the `break` block.
- The same pattern applies to any iterative algorithm: optimisers, MCMC chains,
  fixed-point iterations.
- Add a unit test that checks `fit$converged == TRUE` for a run known to
  converge quickly, and `fit$converged == FALSE` for a run capped at 1
  iteration.

## Related

- See `em-algorithm/R/em_driver.R` for the corrected implementation.
