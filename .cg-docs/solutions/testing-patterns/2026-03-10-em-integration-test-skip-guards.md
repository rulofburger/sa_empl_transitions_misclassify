---
date: 2026-03-10
title: "Testing slow EM/statistical algorithms: skip guards and size tiers"
category: "testing-patterns"
language: "R"
tags: [testthat, integration-tests, em-algorithm, slow-tests, skip-guards, ci, synthetic-data]
root-cause: "Integration tests with large synthetic datasets (n=500+) timeout in both CI and interactive test_dir() runs when the algorithm under test has not yet been vectorised."
severity: "P2"
---

# Testing slow EM/statistical algorithms: skip guards and size tiers

## Problem

Integration tests for an EM algorithm used `n = 500` and `n = 800` to check
parameter recovery. These tests timed out:

- `n = 100`, 30 iterations: **timeout** in interactive `test_dir()`
- `n = 500`, single run: **timeout** in CI

This prevented the test suite from running at all, masking other failures.

## Root Cause

The algorithm had a per-observation R loop making it O(N) in interpreted code.
Even after vectorisation, large-N tests are inherently slower and should not run
unconditionally in every test_dir() invocation.

Two separate issues:

1. **Algorithm not yet vectorised** — immediate fix is vectorise (see
   `2026-03-10-vectorise-r-estep-over-obs.md`).
2. **Test size not tiered** — structural issue that persists even after
   vectorisation: parameter-recovery tests genuinely need large N, but unit
   tests should run fast.

## Solution

### 1. Reduce N for the "always-run" tier

For correctness checks (convergence, gamma dimensions, column validation),
small N is sufficient. Reduce to the minimum that exercises the code path:

```r
# Before
df <- simulate_panel(n = 500, seed = 42)
# After: n=100 is enough to test convergence logic
df <- simulate_panel(n = 100, seed = 42)
```

Rule of thumb for EM integration tests:
- Structural tests (output shape, convergence flag): n = 50–100
- Parameter recovery tests (θ estimates close to truth): n = 300–800

### 2. Gate large-N tests behind an environment variable

```r
test_that("larger sample recovers parameters more precisely", {
  skip_if_not(
    identical(Sys.getenv("EM_FULL_TESTS"), "true"),
    "Set EM_FULL_TESTS=true to run large-N tests"
  )
  df  <- simulate_panel(n = 800, seed = 1234, ...)
  fit <- em_fit_tenure(df, max_iter = 200, verbose = 0)
  # ... assertions ...
})
```

Run locally when needed:
```bash
EM_FULL_TESTS=true Rscript -e "testthat::test_dir('EM-tenure/tests/testthat')"
```

Add to CI only on scheduled runs, not every push:
```yaml
# .github/workflows/full-tests.yml
schedule:
  - cron: '0 2 * * 0'   # Sundays at 2 AM
env:
  EM_FULL_TESTS: "true"
```

### 3. Add a canonical test runner

Create `tests/testthat.R` so tests can be found by standard tooling:

```r
# tests/testthat.R
library(testthat)
test_dir(file.path(dirname(sys.frame(1)$ofile), "testthat"))
```

## Test Size Guidelines

| Test type | Recommended n | Notes |
|---|---|---|
| Output structure (gamma dims, suff stats) | 20–50 | Just needs valid data |
| Convergence flag | 50–100 | Needs enough data to converge |
| Parameter sign/direction | 100–200 | Rough recovery only |
| Quantitative parameter recovery | 500–1000 | Gate with `EM_FULL_TESTS` |
| Bootstrap / SE validation | 1000+ | Always gate |

## Prevention

**Checklist for any statistical algorithm test suite:**

- [ ] All unconditional tests complete in < 10 seconds total
- [ ] Large-N / long-running tests are behind `skip_if_not()` with env var
- [ ] Synthetic data generator accepts `seed` for reproducibility
- [ ] `tests/testthat.R` runner exists alongside `tests/testthat/`
- [ ] Helper files in `tests/testthat/` are named `helper-*.R` (auto-loaded by testthat)

**Pattern**: Use `simulate_*()` functions that accept `seed` so tests are
deterministic:

```r
# ✓ reproducible
df <- simulate_panel(n = 100, seed = 42)

# ✗ non-deterministic, flaky tests
df <- simulate_panel(n = 100)
```

## Related

- `2026-03-10-vectorise-r-estep-over-obs.md` — the performance fix that made fast tests feasible
- `testthat::skip_if_not()` — R docs for conditional test skipping
- `testthat::skip_on_ci()` — alternative for CI-only skipping (less flexible than env var approach)
