---
date: 2026-03-10
title: "Vectorising a per-observation EM E-step in R"
category: "performance-issues"
language: "R"
tags: [em-algorithm, vectorisation, e-step, matrix-ops, mixture-model, latent-variable]
root-cause: "Per-observation for loop in R is O(N) interpreted overhead; E-step is vectorisable because all observations share the same parameter set and latent structure."
severity: "P1"
---

# Vectorising a per-observation EM E-step in R

## Problem

The E-step of an EM algorithm looped over all N observations:

```r
for (i in seq_len(N)) {
  s  <- as.numeric(df[i, c("y1", "y2", "y3")])
  g  <- as.numeric(df[i, c("tenure1", "tenure2", "tenure3")])
  # ... compute responsibilities for observation i ...
  gamma_mat[i, ] <- gamma_i
}
```

At N = 100 with 30 EM iterations this timed out. N = 10 with 3 iterations took ~1 second. For N = 5,000 (typical survey dataset) the algorithm was completely unusable.

## Root Cause

Per-observation `for` loops in interpreted R incur a ~constant interpreter overhead per iteration regardless of body complexity. With H = 8 latent histories and T = 3 waves, the inner body is simple enough that the loop overhead dominated the actual computation. The key insight: **all observations share the same `hmat` (8 × 3), the same parameters, and the same emission function structure** — only the data vectors (s, g, d) differ per observation. This is exactly the structure R's vectorised operations are designed for.

## Solution

### 1. Change the data layout

Instead of accessing `df[i, col]` inside a loop, extract N-vectors once before the loop:

```r
s1 <- df$y1;       s2 <- df$y2;       s3 <- df$y3
g1 <- df$tenure1;  g2 <- df$tenure2;  g3 <- df$tenure3
d1 <- df$timegap1; d2 <- df$timegap2; d3 <- df$timegap3
wi <- df$weight
```

### 2. Build N × H matrices instead of length-H vectors

Replace the per-observation `ld` vector with an `N × H` matrix:

```r
ld <- matrix(0, nrow = N, ncol = H)
lq <- matrix(0, nrow = N, ncol = H)
```

### 3. Loop over H (fixed = 8), not over N

The outer loop becomes over histories (H = 8, fixed), not over observations:

```r
for (j in seq_len(H)) {
  # mask is a length-N logical vector
  mask <- (s1 == 1) & (h1[j] == 1)
  if (any(mask)) {
    ld[mask, j] <- ld[mask, j] + log_emg(g1[mask], lambda_g, sigma2_g)
  }
}
```

All emission functions (log_emg, log_emission_increment_g, etc.) already accept vectors, so they work directly on `g1[mask]`.

### 4. Responsibilities and log-likelihood in one matrix op

```r
log_kernel <- sweep(lq + ld, 2, log_prior, "+")   # N x H
row_max    <- apply(log_kernel, 1, max)
log_denom  <- row_max + log(rowSums(exp(log_kernel - row_max)))
gamma_mat  <- exp(log_kernel - log_denom)
ll         <- sum(wi * log_denom)
```

### 5. Sufficient statistics without the obs loop

Markov statistics reduce to column sums of the weighted responsibility matrix:

```r
wg    <- gamma_mat * wi          # N x H (broadcast wi as column)
wg_cs <- colSums(wg)             # length-H
C1    <- sum(wg_cs * h1)         # scalar
D1    <- sum(wg_cs * (h1 == 1)) + sum(wg_cs * (h2 == 1))
```

Variance statistics still loop over H=8 but apply vectorised `sum()` over masked N-subsets:

```r
for (j in seq_len(H)) {
  wj   <- wg[, j]
  mask <- (s2 == 1) & (s1 == 1) & (h1[j] == 1) & (h2[j] == 1)
  if (any(mask)) {
    dg <- g2[mask] - g1[mask] - .QUARTER_YEARS
    Sg <- Sg + sum(wj[mask] * dg^2)
    Ng <- Ng + sum(wj[mask])
  }
}
```

## Performance

| N | Before (loop) | After (vectorised) | Speedup |
|---|---|---|---|
| 10 | ~1 s / iteration | < 0.01 s / iteration | ~100× |
| 100 | timeout (>30 s) | ~0.1 s / iteration | >300× |
| 500 | timeout | ~0.5 s / iteration | usable |

## Prevention

**Pattern**: In R, any computation of the form "apply the same function to each row of a data frame with fixed parameters" should be vectorised. The E-step of mixture/latent-variable models always has this structure:

```
gamma[i, h] = f(data[i, ], params)   for fixed params
```

This is always vectorisable by transposing the loop structure: iterate over the fixed-size latent dimension (H), applying R vectorised ops over the N-sized data dimension.

**Anti-pattern to avoid**:
```r
for (i in seq_len(N)) {
  row_data <- df[i, ]   # R data frame row access is slow
  result[i] <- compute(row_data, params)
}
```

**Preferred pattern**:
```r
col_data <- df$col      # extract N-vectors once
for (j in seq_len(H)) {  # fixed small dimension
  mask <- (col_data == condition[j])
  result[mask, j] <- compute(col_data[mask], params)
}
```

## Related

- See `EM-tenure/R/estep.R` for the complete vectorised implementation.
- If N > 100k and H grows larger, consider porting the H-loop to C++ via Rcpp (`Rcpp::NumericMatrix` + `Rcpp::LogicalVector`).
- `2026-03-10-em-gem-linked-ctmc.md` (data-quality) — GEM issue discovered in the same EM; performance and numerical correctness are related.
- `2026-03-10-em-integration-test-skip-guards.md` (testing-patterns) — test size tiering that became feasible after vectorisation.
