---
date: 2026-03-15
title: "Testing a vectorised EM E-step and M-step without external data files"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, estep, mstep, testthat, HMM, inline-data, monotonicity, data.table]
root-cause: "No test suite existed for the EM E-step/M-step; adding one required constructing minimal self-contained fixtures that exercise all emission branches"
severity: "P2"
---

# Testing a Vectorised EM E-step and M-step

## Problem

The `em-algorithm` package had zero tests.  Without them, the vectorisation
refactor of `e_step()` could silently break correctness (e.g. wrong row/column
indexing, mismatched obs-type and history-type subsets).

## Pattern

### Fixture design

Keep test data **inline** — no external files.  Use `data.table` for
performance consistency with the production path.

```r
make_df <- function() {
  data.table(
    y1 = c(1L, 0L, 1L), y2 = c(1L, 1L, 0L), y3 = c(1L, 1L, 0L),
    tenure1  = c(0.25, 0.10, 0.30), tenure2 = c(0.50, 0.35, 0.20),
    tenure3  = c(0.75, 0.60, 0.10),
    timegap1 = c(0.10, 0.25, 0.10), timegap2 = c(0.10, 0.50, 0.35),
    timegap3 = c(0.10, 0.75, 0.60),
    weight   = c(1.0,  1.0,  1.0)
  )
}
make_params <- function() {
  list(alpha=0.5, theta0=0.2, theta1=0.8, pi=0.05,
       sigma2_g=0.04, sigma2_h=0.04, lambda_g=2.0, lambda_h=2.0)
}
```

Choose values that exercise **all 4×3 emission branches per wave**:
- Mix of employed/nonemployed observations
- Some continuations, some transitions
- Non-trivial `pi` so misclassification branch fires

### Core invariants to test

| Test | What it catches |
|------|----------------|
| `rowSums(gamma) == 1` | Normalisation bug in log-sum-exp |
| `gamma >= 0` | Sign error in log-posterior |
| `C1 + C0 == N` | Sufficient stat accumulation error |
| `Sg >= 0`, `Sh >= 0` | Variance sum sign flip |
| `ll` is finite scalar | NaN/-Inf in any emission path |

### pi=0 strict-matching test

```r
# With pi=0, only histories matching obs on ALL three waves get weight
for (i in seq_len(nrow(df))) {
  s <- c(df$y1[i], df$y2[i], df$y3[i])
  for (h in seq_len(8L)) {
    if (any(hmat[h, ] != s))
      expect_equal(result$gamma[i, h], 0, tolerance = 1e-14)
  }
}
```

This is the most targeted regression test for the `outer()` misclassification
matrix: any wrong column/row indexing shows up as non-zero weight on a
mismatched history.

### One-cycle monotonicity test

EM guarantees non-decreasing log-likelihood.  One E→M cycle is cheap to test:

```r
r1 <- e_step(df, params)
params2 <- m_step(r1$suff, sum_w, params)
r2 <- e_step(df, params2)
expect_gte(r2$loglik, r1$loglik - 1e-8)
```

The `1e-8` tolerance absorbs floating-point rounding; a genuine M-step bug
produces differences of order 0.01+.

### Linear weight scaling test

```r
df_w2 <- copy(df_w1); df_w2[, weight := 2.0]
expect_equal(e_step(df_w2, p)$loglik, 2 * e_step(df_w1, p)$loglik, tol=1e-10)
```

Catches any hard-coded `sum(weight)` vs weight-per-obs confusion.

### Input validation tests

```r
df_bad <- make_df(); df_bad$tenure1[1] <- -0.1
expect_error(e_step(df_bad, make_params()), "non-negative")
```

### `cbind()` names in matrix row extraction

`clocks_from_histories()` uses `cbind(g1, g2, g3)` which retains column
names.  Single-row extraction returns a named vector:

```r
clks$Gstar[row, ]
# g1   g2   g3
# 0.25  0  0.25
```

Comparing this to `c(0.25, 0, 0.25)` fails `expect_equal()` on the `names`
attribute.  Use `unname()`:

```r
expect_equal(unname(clks$Gstar[row, ]), c(0.25, 0, 0.25))
```

## Prevention

- Always test `rowSums(gamma)` and `C1 + C0 == N` after any refactor of the
  E-step.
- The pi=0 strict-matching test is a cheap and powerful smoke test: run it
  first when debugging responsibility anomalies.
- Use `unname()` when comparing named-vector outputs to bare `c()` literals.

## Related

- `performance-issues/2026-03-15-vectorise-em-estep-over-observations.md` —
  the refactor that motivated this test suite.
- `data-quality/2026-03-15-em-converged-flag-trace-trimming-bug.md` —
  companion bug found in the same review pass.
