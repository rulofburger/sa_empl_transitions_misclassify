---
date: 2026-03-15
title: "Vectorise EM E-step over N observations using N×H matrix operations"
category: "performance-issues"
language: "R"
tags: [em-algorithm, estep, vectorisation, data.table, matrix, HMM, performance]
root-cause: "E-step implemented as a for-loop over N observations; each iteration extracted a single row from a data.frame and looped over H=8 latent histories, making the algorithm O(N) in interpreted R"
severity: "P1"
---

# Vectorise EM E-step over N Observations

## Problem

The E-step of the EM algorithm for the 3-wave employment HMM was implemented
as `for (i in seq_len(N))` loop.  Inside each iteration:

- A single data row was extracted via `df[i, c("y1","y2","y3")]`
- Log-emission contributions were accumulated into a length-8 vector
- Sufficient statistics were accumulated into scalars

With N = 10 000+ observations this is catastrophically slow: R's interpreter
overhead dominates and the code cannot exploit BLAS/LAPACK or data.table's
multi-thread engine.

## Root Cause

The original scalar-centric design was correct but not vectorised.  Every
intermediate result (`lq`, `ld`, `wg`) was length-8 (one per latent history)
rather than N×8 (one per observation × history).  Sufficient statistics were
accumulated inside the loop with `suff_x <- suff_x + ...`, preventing any
vectorised reduction.

## Solution

Lift all computation to **N×H matrices** (`N = nrow(df)`, `H = 8`).

### 1. Extract column vectors once before any computation

```r
s1 <- as.numeric(df[["y1"]]); s2 <- as.numeric(df[["y2"]]); s3 <- as.numeric(df[["y3"]])
g1 <- as.numeric(df[["tenure1"]]); ...
w  <- as.numeric(df[["weight"]])
```

`df[[col]]` is faster than `df[, col]` for both `data.frame` and
`data.table`; it avoids S3 dispatch and returns the vector directly.

### 2. Build N×H emission matrices via `outer()`

```r
# Misclassification: N x H logical outer product
lq <- matrix(0, nrow = N, ncol = H)
lq <- lq +
  ifelse(outer(s1, h1, "=="), log1p(-pi_b), log(pi_b)) +
  ifelse(outer(s2, h2, "=="), log1p(-pi_b), log(pi_b)) +
  ifelse(outer(s3, h3, "=="), log1p(-pi_b), log(pi_b))
```

`outer(s1, h1, "==")` produces an N×8 logical matrix in one C call.

### 3. Fill emission matrix by subsetting rows and columns

For each (obs-type, history-type) cell:

```r
hcols    <- which(h1 == 1)          # length-H subset
obs_rows <- which(s1 == 1)          # length-N subset
val      <- log_emission_start_emg_g(g1[obs_rows], lambda_g, sigma2_g)
ld[obs_rows, hcols] <- ld[obs_rows, hcols] + val
```

This avoids allocating the full N×H matrix for each case: most cases apply to
a small subset of rows and columns.

### 4. Vectorised posterior and log-likelihood

```r
log_post_kernel <- sweep(lq + ld, 2L, log_prior, "+")  # broadcast length-H prior
log_denom       <- apply(log_post_kernel, 1L, .logsumexp)   # length-N
gamma_mat       <- exp(log_post_kernel - log_denom)          # N x H
ll              <- sum(w * log_denom)                        # scalar
```

### 5. Vectorised sufficient statistics

```r
wg <- gamma_mat * w        # N x H, weight broadcast over rows

C1  <- sum(wg * outer(rep(1, N), as.numeric(h1 == 1)))
T11 <- sum(wg * outer(rep(1, N), as.numeric(h1==1 & h2==1))) +
       sum(wg * outer(rep(1, N), as.numeric(h2==1 & h3==1)))

# Variance sums (continuations only)
hcols_g2 <- which(h2 == 1 & h1 == 1)
Sg <- sum(wg[, hcols_g2] * (g2 - g1 - increment)^2) + ...
```

## Prevention

When writing any EM or HMM E-step for tabular panel data:

- **Design with matrices from the start**: `gamma` is N×H, `wg` is N×H.
  Sufficient statistics are sums over this matrix, not loop accumulators.
- **Never extract single rows** inside iteration loops for vectorisable
  operations.  Use `df[[col]]` to extract entire columns before the loop.
- **Use `outer()` for cross-products** of observation-level and
  history-level indicator vectors.
- **Test with N=1 and N=1000** to catch both correctness and performance
  regressions.

## Related

- See `em-algorithm/R/estep.R` for the full vectorised implementation.
- See `testing-patterns/2026-03-15-em-estep-mstep-tests.md` for the
  accompanying test suite that validates correctness of the refactor.
