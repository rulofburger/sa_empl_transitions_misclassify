---
date: 2026-04-29
title: "EM monotonicity guard doubles iteration cost when M-step is unconstrained"
category: "performance-issues"
language: "R"
tags: [em-algorithm, monotonicity-guard, mstep, stationary, linked, performance, convergence]
root-cause: "Verification E-step ran unconditionally after every M-step, doubling E-step cost even though unconstrained EM algebraically guarantees ascent"
severity: "P2"
---

# EM Monotonicity Guard Doubles Iteration Cost (Unconstrained Case)

## Problem

The monotonicity guard in `em_fit_tenure()` and `em_fit_tenure_rho()` ran a
full E-step verification after every M-step iteration:

```r
# SLOW: runs every iteration regardless of whether M-step is constrained
estep_check <- e_step(df, params, discrete_timegap = discrete_timegap)
if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
  params <- params_prev; converged <- TRUE; break
}
```

For the **default specification** (`stationary = FALSE`, `linked = FALSE`),
this guard always passes and costs **one full E-step per iteration**, effectively
doubling total EM runtime.

## Root Cause

The standard EM M-step with all-free parameters algebraically guarantees
log-likelihood ascent (`EM guarantee: Q(θ+) ≥ Q(θ)`). The guard was only added
to handle the `stationary = TRUE` case where `alpha` is tied to `theta0` and
`theta1` via a constraint that can violate ascent. It was mistakenly applied
unconditionally.

See also [2026-03-20-em-monotonicity-guard.md](../testing-patterns/2026-03-20-em-monotonicity-guard.md)
for the original diagnosis.

## Solution

Gate the guard on `stationary || linked`:

```r
# --- Monotonicity guard ---
# Only needed when the M-step is constrained (stationary alpha or CTMC link);
# unconstrained EM guarantees ascent algebraically, so skip the check.
if (stationary || linked) {
  estep_check <- e_step(df, params, discrete_timegap = discrete_timegap)
  if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
    params <- params_prev; converged <- TRUE; break
  }
}
```

Apply identically in `em_fit_tenure_rho()` (replace `e_step` with `e_step_rho`).

## Impact

| Specification | Before | After |
|---|---|---|
| `stationary=FALSE, linked=FALSE` (default) | 2× E-step cost/iter | 1× (guard skipped) |
| `stationary=TRUE` or `linked=TRUE` | 2× E-step cost/iter | 2× (guard needed) |

**~2× speedup** for the most common use case with no change in results.

## Prevention

When adding a post-M-step check that requires a full forward pass, always ask:
> "Is this check algebraically necessary, or only for certain constrained modes?"

Gate expensive checks on the condition that makes them necessary. Document
explicitly which M-step variants require the check and why.

## Related

- [2026-03-20-em-monotonicity-guard.md](../testing-patterns/2026-03-20-em-monotonicity-guard.md) — original guard diagnosis and implementation
- [2026-04-29-discrete-emission-lookup-tables.md](2026-04-29-discrete-emission-lookup-tables.md) — other EM performance fix from same session
- [2026-04-29-log-emg-erfc-underflow.md](2026-04-29-log-emg-erfc-underflow.md) — P1 numerical stability fix from same session
- ⚠️ [2026-04-30-em-guard-unconditional-for-approximate-mstep.md](../testing-patterns/2026-04-30-em-guard-unconditional-for-approximate-mstep.md) — **EXCEPTION**: the eps (Spec I) model has a custom non-exponential-family M-step for ε that breaks the algebraic ascent guarantee even in the free variant. For models with approximate M-steps, the guard must run unconditionally (do NOT apply this gating pattern).
