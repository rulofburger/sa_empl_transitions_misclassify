---
date: 2026-04-30
title: "EM monotonicity guard must be unconditional when M-step is approximate"
category: "testing-patterns"
language: "R"
tags: [em-algorithm, monotonicity, convergence, brent-solver, mstep, approximate-mstep, eps-model, ascent-property]
root-cause: "Gating the monotonicity guard on `stationary || linked` assumed unconstrained EM always guarantees ascent, but the eps model's custom mixture M-step is approximate and can violate monotonicity even in the free variant."
severity: "P1"
---

# EM Monotonicity Guard Must Be Unconditional for Approximate M-steps

## Problem

After gating the guard on `stationary || linked` (per
[2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md)),
the eps (Spec I) driver `em_fit_tenure_eps()` was written with the same gating
pattern:

```r
if (stationary || linked) {
  estep_check <- e_step_eps(df, params)
  if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) { ... }
}
```

Test 17 (`em_fit_tenure_eps: LL non-decreasing across iterations`) failed across
all three seeds with genuine LL decreases of **−0.9 to −6.0 nats** — well outside
numerical noise, occurring in the **free** (`stationary=FALSE, linked=FALSE`) variant.

## Root Cause

The standard EM ascent guarantee holds when the M-step is the **exact maximizer**
of Q(θ | θ_old). It requires:

```
Q(θ_new | θ_old) ≥ Q(θ_old | θ_old)    [condition for L(θ+) ≥ L(θ)]
```

This is guaranteed for standard exponential family models (closed-form M-step) and
for **exact** Brent-solver optima. However, the eps model emission is:

```
ell_spell(g) = (1-ε)^K · f_clean(g) + (1-(1-ε)^K) · f_contam(g)
```

This is **not a standard exponential family** in ε. The M-step update:

```r
eps_hat <- Eps_num / Eps_den   # closed-form derived from expected sufficient stats
```

is correct *in expectation under the assumed sufficient statistics*, but it is not
derived as the exact argmax of Q over ε — it is an approximation that follows from
treating the eps contribution as separable. When the eps M-step result is approximate,
the standard EM guarantee is lost and the guard must run unconditionally.

The earlier diagnosis ([2026-04-29 doc](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md))
was correct for the **base/rho models** where every M-step update is an exact
closed-form solution of the Q-function. It cannot be applied without verification to
models with custom, non-exponential-family M-step components.

## Solution

Remove the `stationary || linked` gating; run the guard unconditionally. Also: cache
the guard result as the next iteration's E-step to eliminate the double-evaluation
cost (this benefit now applies to all variants, not just constrained ones):

```r
# --- Monotonicity guard (all variants) ---
# Runs unconditionally because the eps M-step is an approximation and cannot
# rely on the algebraic EM ascent guarantee. The guard result is cached and
# reused as the next iteration's E-step to avoid double-evaluation cost.
estep_check <- e_step_eps(df, params)
if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
  params    <- params_prev
  converged <- FALSE   # revert: best-known params retained, not converged
  history[[iter]] <- unlist(params)
  break
}
cached_estep <- estep_check   # reuse next iteration
```

**Cost**: One E-step per iteration (same as the gated case for constrained models).
For the free variant this adds one E-step per iteration vs. the gated pattern —
but this cost was already implied by the monotonicity test failures; we cannot
safely skip the guard.

## Test Pattern

Strengthen LL monotonicity tests to use multiple seeds and a larger sample to give
real statistical power:

```r
test_that("em_fit_tenure_eps: LL non-decreasing across iterations (free)", {
  for (seed in c(1L, 42L, 999L)) {
    df  <- .make_eps_data(n = 300L, seed = seed)   # pass seed to helper
    fit <- em_fit_tenure_eps(df, max_iter = 50L, verbose = 0L)

    ll_hist <- fit$history$loglik
    diffs   <- diff(ll_hist)
    expect_true(all(diffs >= -1e-6),
                info = sprintf("Seed %d: min diff = %.3e", seed, min(diffs)))
  }
})
```

**Note**: pass seed to the data-generation helper rather than calling `set.seed()`
before it — the helper's internal `set.seed()` overrides any external call.

## Checklist: When to Gate the EM Guard

| M-step type | Gate on `stationary || linked`? |
|---|---|
| Exact closed-form (exponential family) | ✅ Gate — algebraic guarantee holds |
| Exact numerical (true argmax via Brent) | ✅ Gate — guarantee holds |
| Approximate (custom sufficient-stat formula for non-expfam mixture) | ❌ Unconditional — guarantee may fail |
| Any M-step with unknown guarantee | ❌ Unconditional until verified |

## Prevention

Before gating the monotonicity guard on `stationary || linked` in a new EM driver,
verify that **every parameter update** in the unconstrained M-step is either:
(a) a closed-form exact maximizer of Q, or
(b) a numerical exact maximizer (Brent with tolerance ≪ convergence tolerance).

If any update is an approximation derived from expected sufficient statistics (not
from direct Q maximization), the guard must run unconditionally.

## Related

- [2026-04-29-em-monotonicity-guard-gating.md](../performance-issues/2026-04-29-em-monotonicity-guard-gating.md) — original gating fix (valid for base/rho models with exact M-steps)
- [2026-03-20-em-monotonicity-guard.md](2026-03-20-em-monotonicity-guard.md) — original guard implementation
- [2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md](2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md) — E-step accumulation guard stale after emission contract change (third failure mode)
- EM-tenure/R/em_driver_eps.R — fix applied here
