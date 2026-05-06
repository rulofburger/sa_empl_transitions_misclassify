---
date: 2026-05-01
title: "Stale E-step accumulation guard after emission refactor — K=1/offset>0 silently excluded from eps sufficient stats"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, e-step, m-step, eps-model, accumulation-guard, emission-contract, monotonicity, sufficient-statistics, k1-singleton, refactor]
root-cause: "The E-step used a hardcoded `K >= 2L` guard to decide which spells contribute to eps sufficient statistics. After the emission function was refactored to make K=1/offset>0 spells eps-informative (with tau_sum > 0), the guard was never updated. The E-step and M-step then optimised different Q-functions, breaking EM monotonicity for eps."
severity: "P1"
---

# Stale E-step Accumulation Guard After Emission Refactor

## Problem

After refactoring `emissions_eps.R` to correctly handle K=1 spells where
the observed wave is not the spell start (`t_1 > a`, offset > 0), the EM
algorithm's eps estimate was systematically biased.

Specifically:

- The emission function now correctly returns `tau_sum > 0` for K=1/offset>0
  spells (they carry eps information because the clean and contaminated
  branches differ by a shift of `(t_1 - a) * Delta`).
- But the E-step in `estep_eps.R` still used the old guard:
  ```r
  mask_ge2 <- out$K >= 2L
  ```
  and silently discarded K=1/offset>0 tau_sum values from `Eps_num`/`Eps_den`.

The emission Q-function now included eps terms for K=1/offset>0 observations,
while the M-step denominator `Eps_den` did not count them. The E-step and
M-step were no longer maximising the same objective — **EM monotonicity was
broken for eps**.

No error was raised. The algorithm ran to convergence but produced an eps
estimate from K>=2 spells only, discarding information from K=1/offset>0
singletons.

## Root Cause

The `K >= 2L` guard was derived from the two-branch ("all-or-nothing")
approximation where K=1 spells *genuinely* have tau_sum = 0 (both patterns
give the same density). The guard was correct for that model.

When the 2^K exact enumeration was introduced, K=1/offset>0 became
eps-informative, but the guard was inherited verbatim without re-derivation.
The mathematical condition "this spell informs eps" changed, but the code
condition did not.

**Pattern: guards gated on model properties are brittle across emission
changes.** Any guard of the form `if (some_structural_property)` must be
re-derived from first principles whenever the emission formula changes what
it returns.

## Solution

### Step 1: Add `eps_informative` to the emission return list

In `emissions_eps.R`, introduce a logical vector `eps_informative` set to
`TRUE` for every spell where the emission varies with eps:

```r
eps_informative <- logical(N)          # initialise FALSE

# K=1, offset>0: eps is identifiable (clean/contaminated branches differ)
eps_informative[idx_sh] <- TRUE

# K>=2: always eps-informative
eps_informative[mask2]  <- TRUE
eps_informative[mask3]  <- TRUE

# Return it
list(loglik = loglik, tau_sum = tau_sum, K = K,
     eps_informative = eps_informative,
     lambda_count = lambda_count, lambda_xsum = lambda_xsum)
```

The decision of "is this spell informative for eps?" now lives inside the
emission function, where the logic is, not in the E-step.

### Step 2: Drive the E-step guard from `eps_informative`

In `estep_eps.R`, replace the hardcoded structural guard:

```r
# BEFORE (stale — fails after emission refactor)
mask_ge2 <- out$K >= 2L
if (any(mask_ge2)) {
  eps_num_mat[mask_ge2, j] <- eps_num_mat[mask_ge2, j] + out$tau_sum[mask_ge2]
  eps_den_mat[mask_ge2, j] <- eps_den_mat[mask_ge2, j] + out$K[mask_ge2]
}

# AFTER (driven by emission contract)
mask_eps <- out$eps_informative
if (any(mask_eps)) {
  eps_num_mat[mask_eps, j] <- eps_num_mat[mask_eps, j] + out$tau_sum[mask_eps]
  eps_den_mat[mask_eps, j] <- eps_den_mat[mask_eps, j] + out$K[mask_eps]
}
```

### Step 3: Add integration test

```r
test_that("e_step_eps: K=1 offset>0 spells contribute to Eps_den", {
  df    <- .make_eps_data(n = 300L, seed = 123L)
  df$y1 <- 0L; df$y3 <- 0L    # force s=(0,1,0) patterns
  df$tenure1 <- 0.5; df$tenure3 <- 0.5

  out <- e_step_eps(df, .make_eps_params())

  # Under h=(1,1,0) or h=(0,1,1), wave 2 is K=1/offset>0 → must contribute
  expect_true(out$suff$Eps_den > 0)
  expect_true(out$suff$Eps_num >= 0)
  expect_true(out$suff$Eps_num <= out$suff$Eps_den + 1e-10)
})
```

## Prevention

**Rule: emission return contracts must drive downstream guards.**

When an emission function changes what it returns (new fields, changed
semantics of existing fields), audit every consumer of that function for
guards that may now be stale.

Concrete checklist after any emission refactor:
1. Does the emission now return non-zero values for rows it previously
   returned zero for? → Check all accumulation guards in the E-step.
2. Does the emission's "informativeness" condition change? → Add an
   `_informative` boolean to the return list and drive all guards from it.
3. Are there comments in the E-step that refer to specific K values or
   other structural properties? → Treat them as stale-guard red flags.

**Anti-pattern**: Hardcoding structural properties (`K >= 2L`, `stationary`,
`linked`) in the E-step accumulation guard. These properties reflect the
model structure at a point in time; when the emission changes, the guard
must be re-derived, not inherited.

**Canonical guard placement**: The "is this spell informative for parameter
θ?" decision must be made by the function that computes the θ-relevant
density (the emission), not by its consumer (the E-step).

## Related

- [2026-05-01-monotonicity-guard-condition-dropped.md](2026-05-01-monotonicity-guard-condition-dropped.md)
  — Different failure mode: guard condition accidentally *dropped* entirely
  (unconditional revert). This document: guard condition is present but
  *stale* after emission refactor.
- [2026-04-30-em-guard-unconditional-for-approximate-mstep.md](2026-04-30-em-guard-unconditional-for-approximate-mstep.md)
  — Guards gated on model properties (`stationary || linked`) break for
  approximate M-steps; related class of "structural property used as proxy
  for guard condition" failures.
- [2026-03-20-em-monotonicity-guard.md](2026-03-20-em-monotonicity-guard.md)
  — Original EM ascent property guard pattern.
- [2026-03-15-em-estep-mstep-tests.md](2026-03-15-em-estep-mstep-tests.md)
  — Testing E-step/M-step consistency patterns; E-step sufficient stat
  accumulation tests.
