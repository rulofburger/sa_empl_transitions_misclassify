---
date: 2026-05-01
title: "Aggregate N-length sufficient statistics to K-bin histograms before passing to Brent FOC evaluator"
category: "performance-issues"
language: "R"
tags: [em-algorithm, mstep, brent, lambda-d, sufficient-statistics, tapply, histogram, aggregation, performance, discrete-timegap]
root-cause: "Lambda_d Brent FOC evaluator received O(8N) length vectors of raw weights because suff-stats were accumulated per-history then unlisted, causing ~100 full-vector dot-products per M-step despite discrete categories having at most 7 or 49 distinct inputs"
severity: "P1"
---

# Aggregate Sufficient Statistics to K-Bin Histograms Before Brent Solver

## Problem

The lambda_d M-step uses Brent's method to find the root of a first-order
condition (FOC). The FOC is evaluated ~100 times per M-step call. Each
evaluation computes:

```r
# Called ~100 times; cat_marg has O(8N) elements
sum(w_marg * interval_grad_lambda_d(cat_marg, lam))
```

The vectors `cat_marg` and `w_marg` were assembled by unlisting per-history
accumulation lists:

```r
# After E-step loop:
cat_d_marginal_c <- unlist(cat_d_marginal_c_list)  # length ≈ 8N
cat_d_marginal_w <- unlist(cat_d_marginal_w_list)  # length ≈ 8N
```

At N=402K this means each of the ~100 Brent FOC calls processes ~3.2M elements
when only 7 distinct category values exist (marginal) or 49 (transition pairs).

Additionally, the per-history accumulation itself used repeated `c()` appends:

```r
# Inside H loop: O(copy) per append
ec <- c(ec, c2[m]); ew <- c(ew, wj[m])
```

With up to 6 appends per history and N=402K, this is O(6×8×N) copy traffic
just to build the input vectors.

## Root Cause

The E-step accumulated raw observation-level weights for every history,
preserving the full N-length detail even though the Brent FOC only depends on
the weighted count per discrete category level. Discrete timegap categories
take integer values 1–7 (7 marginal bins) or pairs in 1:7 × 1:7 (49 transition
bins). All observations in the same category bin have identical FOC contributions.

## Solution

After accumulating per-history suff-stats, aggregate to bin-level histograms
with `tapply()` before storing:

```r
# --- At the end of per-j marginal accumulation ---
# (ec, ew are the raw N-length category / weight vectors for history j)

# Aggregate: O(N) tapply → ≤7 elements
wt <- tapply(ew, ec, sum)
cat_d_marginal_c_list[[j]] <- as.integer(names(wt))    # ≤7 distinct levels
cat_d_marginal_w_list[[j]] <- as.double(wt)

# --- For transition pairs ---
# (tc, tp: current/previous category; tw: weight vector)
key  <- tc + 7L * (tp - 1L)           # unique key in 1:49
wt2  <- tapply(tw, key, sum)
keys <- as.integer(names(wt2))
cat_d_trans_curr_list[[j]] <- (keys - 1L) %% 7L + 1L
cat_d_trans_prev_list[[j]] <- (keys - 1L) %/% 7L + 1L
cat_d_trans_w_list[[j]]    <- as.double(wt2)
```

After `unlist()`, the combined suff-stats are:
- **Marginal**: ≤ 7 × H = 56 elements (vs O(8N) before)
- **Transition**: ≤ 49 × H = 392 elements (vs O(8N) before)

Each Brent FOC evaluation now does ~56 multiplications instead of ~3.2M.
Total M-step cost drops from O(100 × 8N) to O(100 × 56) — a ~57,000× reduction
for N=402K.

## Prevention

**Before passing vectors to an iterative numerical solver (Brent, Newton, etc.),
check whether the input has bounded discrete structure.** If the function
evaluated inside the solver depends only on category membership (not on
individual observation identity), pre-aggregate to a weighted histogram.

Pattern to follow:

```r
# Before Brent: aggregate O(N) → O(K_levels)
wt <- tapply(weights, categories, sum)
cat_bins    <- as.integer(names(wt))
weight_bins <- as.double(wt)

# Brent FOC now operates on K_levels elements, not N
foc <- function(lam) sum(weight_bins * f(cat_bins, lam))
```

This applies whenever:
1. Categories are bounded integers (e.g., discrete timegap: 7 levels)
2. The solver evaluates the same function with varying scalar parameter
3. N >> K_levels (which is always true at panel survey scales)

## Related

- [2026-04-29-discrete-emission-lookup-tables.md](2026-04-29-discrete-emission-lookup-tables.md) — lookup table strategy for the emission functions themselves (complementary: this doc covers the suff-stat aggregation *before* the solver; that doc covers vectorising *inside* the emission function)
- [2026-05-01-blas-hamming-distance-lq-matrix.md](2026-05-01-blas-hamming-distance-lq-matrix.md) — BLAS strategy for the E-step lq matrix
