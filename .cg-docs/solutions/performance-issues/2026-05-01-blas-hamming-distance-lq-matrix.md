---
date: 2026-05-01
title: "Replace per-history ifelse() loop for misclassification matrix with BLAS Hamming distance"
category: "performance-issues"
language: "R"
tags: [em-algorithm, estep, misclassification, lq-matrix, blas, matrix-multiply, hamming-distance, ifelse, vectorisation, performance]
root-cause: "Misclassification log-probability matrix (N×H) was computed with a per-history loop calling ifelse() 3× per column, causing 24 N-length allocations and 24 full-vector traversals per E-step"
severity: "P1"
---

# BLAS Hamming Distance Replaces Per-History ifelse() Loop for lq Matrix

## Problem

The misclassification log-probability matrix `lq` (N×H, H=8) was built with a
`for (j in seq_len(H))` loop using `ifelse()`:

```r
# SLOW: 24 N-length allocations per E-step (3 ifelse per history × 8 histories)
for (j in seq_len(H)) {
  lq[, j] <- ifelse(s1 == h1[j], lp_match, lp_mismatch) +
              ifelse(s2 == h2[j], lp_match, lp_mismatch) +
              ifelse(s3 == h3[j], lp_match, lp_mismatch)
}
```

At N=402K and 200+ EM iterations this is ~24 × 402K × 8 bytes × 200 ≈ **15 GB**
of allocation traffic in the hot path.

## Root Cause

The pattern treats each of the 8 histories independently, recomputing wave-by-wave
matches via `ifelse()`. But `lq[i, j] = sum_{t} log_p(s_{it}, h_{jt})` is a linear
combination of Hamming distance `d(s_i, h_j) = # waves where s_i ≠ h_j`, which
can be computed for all (i, j) pairs in a single BLAS matrix multiply.

## Solution

Pre-compute the N×H Hamming distance matrix via BLAS `%*%`, then derive `lq`
in one vectorised operation:

```r
# 1. Compute H_mat (3×H binary history matrix) — column j is h_j
H_mat <- rbind(h1, h2, h3)                        # 3 × H, fixed across iterations

# 2. Compute s_full (N×3 logical → numeric 0/1) — built once above for n_mis_mat
#    s_full <- cbind(s1 == 1L, s2 == 1L, s3 == 1L)

# 3. Hamming distance via BLAS: d(s_i, h_j) = sum(s_i) + sum(h_j) - 2*dot(s_i, h_j)
#    Using (s1+s2+s3) instead of rowSums(s_full) avoids logical→integer coercion.
n_mis_mat <- (s1 + s2 + s3) +
             matrix(colSums(H_mat), N, H, byrow = TRUE) -
             2L * (s_full %*% H_mat)               # BLAS dgemm: N×H in one pass

# 4. Derive lq: scalar operations on an already-computed N×H matrix
if (pi_par < .Machine$double.eps) {
  lq <- matrix(-Inf, nrow = N, ncol = H)
  lq[n_mis_mat == 0L] <- 0                         # exact match: log(1) = 0
} else {
  pi_b <- .bound01(pi_par)
  lq   <- n_mis_mat * log(pi_b) + (3L - n_mis_mat) * log1p(-pi_b)
}
```

Key points:
- `n_mis_mat` is also needed later for the M-step `M_count` sufficient statistic —
  compute once and reuse both times.
- `(s1 + s2 + s3)` is marginally faster than `rowSums(s_full)` because it avoids
  R's internal logical→integer coercion in `rowSums`.
- For `pi_par ≈ 0`, using `n_mis_mat == 0L` directly (no log) avoids `log(0)`.

## Prevention

When computing a matrix of pairwise functions `f(row_i, col_j)` over fixed small
sets of columns (e.g., 8 latent histories), check whether `f` decomposes into
linear algebra. For binary inputs, the identity:

```
Hamming(a, b) = sum(a) + sum(b) - 2 * dot(a, b)
```

enables a single BLAS `%*%` call instead of H per-column `ifelse()` traversals.
This is especially impactful when the function is called inside an EM loop
(200+ iterations × N=400K).

Also applies to: `n_mis_mat` in the sufficient statistics block (replaces triple
`outer(s1, h1, "!=") + outer(s2, h2, "!=") + outer(s3, h3, "!=")`).

## Related

- [2026-03-10-vectorise-r-estep-over-obs.md](2026-03-10-vectorise-r-estep-over-obs.md) — general E-step vectorisation strategy
- [2026-04-29-discrete-emission-lookup-tables.md](2026-04-29-discrete-emission-lookup-tables.md) — lookup table strategy for discrete emissions
