---
date: 2026-05-06
title: "<<- in inner functions leaks intent and breaks testability"
category: "data-quality"
language: "R"
tags: [super-assignment, closures, <<-, side-effects, functional-style, prepare-covariates]
root-cause: "Inner helper functions used <<- to accumulate results into parent-frame vectors instead of returning their outputs, making the code harder to test and reason about"
severity: "P2"
---

# `<<-` in inner functions leaks intent and breaks testability

## Problem

`prepare_covariate_matrix()` contained two inner helper functions, `.std()` and
`.onehot()`, that silently mutated `center_vals` and `scale_vals` in the
enclosing function frame via `<<-`:

```r
# PROBLEMATIC PATTERN
.std <- function(x, name) {
  mu  <- mean(x)
  sig <- sd(x)
  center_vals[name] <<- mu   # mutates enclosing frame
  scale_vals[name]  <<- sig  # mutates enclosing frame
  (x - mu) / sig
}

age_std <- .std(df$age1, "age")   # side-effect hidden inside call
```

The same pattern appeared in `estimate_extensions_pipeline.R` with `.add_row()`:
```r
.add_row <- function(label, ..., fit) {
  run_rows[[label]] <<- data.frame(...)   # mutates enclosing list
}
.add_row("cov_s1_sym_stat", ...)   # reader can't see what changes
```

**Consequences**:
- A reader of the call site cannot tell that `center_vals`, `scale_vals`, or
  `run_rows` are being modified — the effect is hidden inside the function name.
- Unit-testing `.std()` or `.onehot()` in isolation requires setting up the
  parent-frame state, which is fragile.
- If the helper is ever extracted to a different scope (e.g., moved to a shared
  helpers file), `<<-` will walk up the wrong call stack and pollute `.GlobalEnv`.

## Root Cause

`<<-` (super-assignment) assigns to the first matching name found by walking
up the call stack. When used inside a closure defined within a function, it
modifies the *enclosing function's* local variables. This works correctly as
long as the inner function is never called from outside that exact enclosing
scope — but it creates fragile, implicit coupling between the inner function
and its context.

## Solution

Return the outputs from the inner function as a named list; accumulate
explicitly in the caller.

**`.std()` refactored**:
```r
# CORRECT: return list, accumulate in caller
.std <- function(x, name) {
  mu  <- mean(x)
  sig <- sd(x)
  if (sig < 1e-10) { warning(...); sig <- 1 }
  list(mu = mu, sig = sig, val = (x - mu) / sig)
}

age_res            <- .std(df$age1, "age")
center_vals["age"] <- age_res$mu    # explicit accumulation
scale_vals["age"]  <- age_res$sig   # visible at call site
col_parts[["age"]] <- matrix(age_res$val, ncol = 1L, dimnames = list(NULL, "age"))
```

**`.add_row()` → `.make_row()` refactored**:
```r
# CORRECT: return the row, let the caller assign it
.make_row <- function(label, family, model_type, stationary, fit) {
  data.frame(label = label, loglik = fit$loglik, ...)
}

run_rows[["cov_s1_sym_stat"]] <- .make_row("cov_s1_sym_stat", ...)
```

## Prevention

- **Rule**: Inner helper functions must not use `<<-`. They must return their
  output as a value (scalar, vector, or named list). The caller is responsible
  for accumulation.
- **Exception**: Module-level memoisation via a private environment is
  acceptable (e.g., `.hmat_env$.cache <- ...`), because the cache object
  explicitly models a persistent store, not a hidden side-effect on a caller's
  local variable.
- **Code review signal**: Any `<<-` in an inner function is a flag for review.
  Ask: "Could this function be pure (return-value only)?" Almost always yes.

## Related

- [2026-05-06-em-hot-path-constant-matrix-caching.md](../performance-issues/2026-05-06-em-hot-path-constant-matrix-caching.md) — the
  memoisation exception where a private environment is appropriate
