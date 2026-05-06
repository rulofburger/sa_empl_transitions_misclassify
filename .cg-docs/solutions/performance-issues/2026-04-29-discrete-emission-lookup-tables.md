---
date: 2026-04-29
title: "Replace per-observation emission loops with lookup tables for discrete categories"
category: "performance-issues"
language: "R"
tags: [em-algorithm, emissions, discrete-timegap, lookup-table, vectorisation, mstep, brent, performance]
root-cause: "Discrete category emissions looped over N observations despite having at most 7 or 49 distinct inputs; called O(800N) times per EM iteration due to nesting inside 8-history E-step and 100-call Brent M-step"
severity: "P2"
---

# Lookup Tables for Discrete Emission Functions (7 or 49 Distinct Inputs)

## Problem

Four emission functions for the `discrete_timegap = TRUE` model looped over N
observations one at a time:

```r
# SLOW: O(N) interpreted iterations per call
log_emission_interval_d <- function(cat, lambda_d) {
  result <- rep(-Inf, length(cat))
  for (i in which(!is.na(cat) & cat >= 1L)) {
    iv <- .timegap_interval(cat[i])
    result[i] <- .log_interval_prob(iv[1], iv[2], lambda_d)
  }
  result
}
```

These functions are called:
- **Inside the 8-history E-step loop**: effectively `O(8N)` per wave per E-step call
- **Inside the Brent M-step FOC evaluator**: `~100` function evaluations per M-step,
  each calling `log_emission_interval_d` on an N-vector

Total: approximately **O(800N)** interpreted loop iterations per EM iteration —
dominated by R's interpreter overhead, not arithmetic.

## Root Cause

The category inputs are integers in `1:7` (marginal emissions) or pairs from
`1:7 × 1:7` (transition emissions). There are exactly **7** and **49** distinct
inputs respectively, regardless of N. Yet the original implementations re-computed
the same 7 (or 49) values for every observation in every call.

## Solution

### Pattern: build lookup table, then vectorised index

**For 7-input functions** (`log_emission_interval_d`, `interval_grad_lambda_d`):

```r
log_emission_interval_d <- function(cat, lambda_d) {
  lut <- vapply(seq_len(.N_TIMEGAP_CATS), function(k) {
    iv <- .timegap_interval(k)
    .log_interval_prob(iv[1], iv[2], lambda_d)
  }, numeric(1))
  result        <- rep(-Inf, length(cat))
  valid         <- !is.na(cat) & cat >= 1L & cat <= .N_TIMEGAP_CATS
  result[valid] <- lut[cat[valid]]   # single vectorised index
  result
}
```

**For 49-input functions** (`log_emission_transition_d`, `transition_grad_lambda_d`):

```r
log_emission_transition_d <- function(cat_curr, cat_prev, lambda_d) {
  K <- .N_TIMEGAP_CATS
  tmat <- matrix(-Inf, K, K)          # 7x7 table
  for (j in seq_len(K)) {
    log_denom <- .log_cat_mass(j, lambda_d)
    if (is.infinite(log_denom) && log_denom < 0) next
    for (k in seq_len(K)) {
      # ... compute intersection interval [L, U) ...
      if (L < U) tmat[j, k] <- .log_interval_prob(L, U, lambda_d) - log_denom
    }
  }
  result <- rep(-Inf, length(cat_curr))
  valid  <- !is.na(cat_prev) & !is.na(cat_curr) &
            cat_prev >= 1L & cat_prev <= K & cat_curr >= 1L & cat_curr <= K
  result[valid] <- tmat[cbind(cat_prev[valid], cat_curr[valid])]  # matrix index
  result
}
```

Key: `tmat[cbind(row_idx, col_idx)]` is a vectorised matrix subscript that
extracts N values in a single C operation — no R loop over N.

### Performance impact

| Function | Old cost per call | New cost per call |
|----------|------------------|------------------|
| `log_emission_interval_d` | O(N) R iterations | O(7) + O(N) C index |
| `log_emission_transition_d` | O(N) R iterations | O(49) + O(N) C index |
| `interval_grad_lambda_d` | O(N) R iterations | O(7) + O(N) C index |
| `transition_grad_lambda_d` | O(N) R iterations | O(49) + O(N) C index |

Net EM iteration speedup is especially large in the M-step (Brent calls each
function ~100×).

## Prevention

**Checklist** for any new emission function that takes a category code:
1. Count distinct inputs. If ≤ 7 inputs (or ≤ K² for pairs), use a LUT.
2. Build the LUT with `vapply()` over the small domain.
3. Use vectorised indexing: `lut[cat[valid]]` or `mat[cbind(r, c)]`.
4. Never `for (i in which(valid))` over N observations.

## Related

- [2026-04-29-log-emg-erfc-underflow.md](2026-04-29-log-emg-erfc-underflow.md) — related numerical stability fix in emissions
- [2026-03-15-vectorise-em-estep-over-observations.md](2026-03-15-vectorise-em-estep-over-observations.md) — vectorising the E-step outer loop
- [2026-04-29-em-monotonicity-guard-gating.md](2026-04-29-em-monotonicity-guard-gating.md) — other EM performance fix from same session
- [2026-05-01-tapply-suff-stats-aggregation-before-brent.md](2026-05-01-tapply-suff-stats-aggregation-before-brent.md) — complementary: aggregating the sufficient statistics *before* passing to the Brent solver (this doc covers vectorising *inside* the emission function)
