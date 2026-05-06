---
date: 2026-05-06
title: "EM hot-path: memoised constant matrix and lapply logsumexp"
category: "performance-issues"
language: "R"
tags: [EM-algorithm, memoization, logsumexp, hot-path, latent-histories, matrix-allocation]
root-cause: "Two unnecessary allocations per E-step iteration: (1) latent_histories() reconstructed its 8x3 constant matrix from scratch on every call; (2) as.data.frame(log_joint) allocated a new data frame just to enable Reduce(pmax,...) column iteration"
severity: "P2"
---

# EM hot-path: memoised constant matrix and lapply logsumexp

## Problem

In a 3-wave binary EM model with N observations, the E-step runs 200–1000
times per estimation. Two allocation patterns inside the E-step had unnecessary
overhead:

1. **`latent_histories()`** was called at the top of every E-step. It always
   returns the same constant 8×3 integer matrix (all 8 binary histories of 3
   waves). Rebuilding it 1000× per estimation was pure waste.

2. **`Reduce(pmax, as.data.frame(log_joint))`** is the standard trick for a
   row-wise maximum over an N×8 matrix without materialising an N-length apply.
   But `as.data.frame()` converts each column to a list element with attribute
   copying — needlessly slow for a numeric matrix already in memory.

## Root Cause

- `latent_histories()` computes `expand.grid(0:1, 0:1, 0:1)` and coerces to
  integer on every call — cheap per call but adds up over thousands of
  iterations across 18 models × 5 starts each.
- `as.data.frame()` on a matrix copies column data and sets column names,
  class, and row.names attributes. `lapply()` on column indices avoids all of
  this.

## Solution

### 1. Memoised history matrix via module-level cache

```r
# helpers_ext.R (or any shared module file)
.hmat_env <- new.env(parent = emptyenv())
.hmat_env$.cache <- NULL

.hmat_cache <- function() {
  if (is.null(.hmat_env$.cache))
    .hmat_env$.cache <- latent_histories()
  .hmat_env$.cache
}
```

Then in every E-step replace:
```r
hmat <- latent_histories()      # OLD: recomputes every call
hmat <- .hmat_cache()           # NEW: returns cached matrix
```

**Why a private environment?** Storing the cache in `.hmat_env` avoids
polluting `.GlobalEnv` with a mutable variable while keeping the cache
persistent across calls within a session.

### 2. `lapply` instead of `as.data.frame` for logsumexp column iteration

```r
# OLD: allocates data.frame with attribute copying
log_max <- Reduce(pmax, as.data.frame(log_joint))

# NEW: zero-copy column iteration via index function
log_max <- Reduce(pmax, lapply(seq_len(ncol(log_joint)), function(j) log_joint[, j]))
```

For a fixed 8-column matrix, `seq_len(8L)` can be hardcoded:
```r
log_max <- Reduce(pmax, lapply(seq_len(8L), function(j) log_joint[, j]))
```

This extracts each column as a plain numeric vector (no attribute copy) and
reduces them with `pmax` — identical result, lower allocation pressure.

## Prevention

- **Rule**: Any constant that is computed inside a function called in a tight
  loop (>10× per estimation) should be precomputed outside the loop or memoised.
  The 8×3 history matrix and any static indicator vectors derived from it
  (e.g., `h1 == 1L`, `h1 == 1L & h2 == 1L`) are candidates.
- **Pattern for N×8 mismatch indicator**: also precompute outside the driver loop:
  ```r
  # Driver, before EM loop
  hmat_s  <- .hmat_cache()
  mm_static <- outer(as.integer(df$y1), hmat_s[,1L], "!=") +
               outer(as.integer(df$y2), hmat_s[,2L], "!=") +
               outer(as.integer(df$y3), hmat_s[,3L], "!=")
  # Pass mm_static to e_step as mm_precomp argument
  ```
- **Avoid `as.data.frame()` for column iteration on numeric matrices**. Use
  `lapply(seq_len(ncol(m)), function(j) m[, j])` instead.

## Related

- [2026-03-10-vectorise-r-estep-over-obs.md](2026-03-10-vectorise-r-estep-over-obs.md)
- [2026-05-05-validate-false-pattern-em-hot-path.md](2026-05-05-validate-false-pattern-em-hot-path.md)
- [2026-05-01-tapply-suff-stats-aggregation-before-brent.md](2026-05-01-tapply-suff-stats-aggregation-before-brent.md)
