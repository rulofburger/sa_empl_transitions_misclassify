---
date: 2026-05-01
title: "K=1 singleton emission incorrectly collapsed to plain Exp when spell start is unobserved"
category: "data-quality"
language: "R"
tags: [em-algorithm, emission, k1-singleton, eps-model, shifted-exp, contamination, two-pattern-mixture, offset, spell-start]
root-cause: "The K=1 emission handler in emissions_eps.R treated all K=1 spells as plain Exp(lambda_g) evaluations, regardless of whether the observed wave was the spell start (offset=0) or a later wave (offset>0). When offset>0, the clean branch gives a *shifted* Exp (T_g = g - d*Delta) while the contaminated branch gives unshifted Exp(g). These differ, so eps is identifiable and a 2-pattern mixture is required."
severity: "P1"
---

# K=1 Singleton Emission Incorrectly Collapsed to Plain Exp

## Problem

The `log_emission_spell_g()` function in `EM-tenure/R/emissions_eps.R`
handled all K=1 spells (exactly one observed tenure in the maximal E-spell)
with a single plain Exp evaluation:

```r
# WRONG: applied unconditionally for all K=1 rows
loglik[mask1] <- log(lambda_g) - lambda_g * g_obs
tau_sum[mask1] <- 0  # eps drops out (assumed)
```

This is correct only when the observed wave is the spell start (`t_1 = a`,
offset = 0). When `t_1 > a` (the spell started in an earlier unobserved
wave), the two contamination patterns give **different densities**:

- Clean branch (`m=0`): `T_g = g - d*Delta` is pinned → shifted Exp density
  `lambda_g * exp(-lambda_g * (g - d*Delta))` for `g >= d*Delta`
- Contaminated branch (`m=1`): `T_g` is free → unshifted Exp density
  `lambda_g * exp(-lambda_g * g)`

By treating both as identical, the emission was systematically wrong for
K=1/offset>0 spells, and eps was treated as unidentifiable from singletons
even when it is.

## Root Cause

The TeX spec (Section 3.4, Singletons paragraph) correctly states:

> When `t_1 > a` (the spell start is unobserved), the clean and contaminated
> branches differ by a shift of `(t_1 - a)*Delta`, and eps remains active.

But the code comment said only:
```r
# K=1: plain Exp; eps drops out
```
— conflating the `t_1 = a` case (where eps truly drops out) with the general
K=1 case. The offset condition was never checked.

This class of spell arises in a 3-wave panel whenever an E-spell spans
waves 1–2 and only wave 2 is observed (`s=(0,1,*)` with `h=(1,1,*)`),
or similarly for the 2–3 bridge.

## Solution

Split the K=1 handler into two branches based on offset:

```r
mask1 <- K == 1L
if (any(mask1)) {
  # ... extract g_obs, d_obs per row ...
  at_start <- d_obs == 0L   # t_1 = a: offset = 0
  idx_at   <- which(mask1)[at_start]
  idx_sh   <- which(mask1)[!at_start]

  # Branch 1: t_1 = a — eps drops out, plain Exp
  if (length(idx_at) > 0L) {
    loglik[idx_at]       <- log_lam - lambda_g * g_obs[at_start]
    lambda_count[idx_at] <- 1L
    lambda_xsum[idx_at]  <- g_obs[at_start]
    # eps_informative stays FALSE; tau_sum stays 0
  }

  # Branch 2: t_1 > a — 2-pattern mixture (shifted vs unshifted Exp)
  if (length(idx_sh) > 0L) {
    g1s <- g_obs[!at_start]; d1s <- d_obs[!at_start]
    T_g <- g1s - d1s * Delta           # T_g from clean branch
    v   <- T_g > 0                     # clean branch valid iff T_g > 0

    lp_clean  <- rep(-Inf, length(g1s))
    if (any(v)) lp_clean[v] <- log_1me + log_lam - lambda_g * T_g[v]
    lp_contam <- log_e + log_lam - lambda_g * g1s

    lp_mx  <- pmax(lp_clean, lp_contam)
    ll_k1  <- lp_mx + log(exp(lp_clean - lp_mx) + exp(lp_contam - lp_mx))
    nu     <- exp(lp_contam - ll_k1)   # P(contaminated | g)
    T_safe <- pmax(T_g, 0)

    loglik[idx_sh]          <- ll_k1
    tau_sum[idx_sh]         <- nu
    eps_informative[idx_sh] <- TRUE    # signals E-step to accumulate
    lambda_count[idx_sh]    <- 1L
    lambda_xsum[idx_sh]     <- nu * g1s + (1 - nu) * T_safe
  }
}
```

Key points:
- Use **precomputed indices** (`idx_at`, `idx_sh`) rather than double-subscript
  `loglik[mask1][at_start]` which creates intermediate copies.
- Set `eps_informative = TRUE` for offset>0 rows so the E-step guard can
  be driven by the emission contract (see companion solution
  `2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md`).
- `lambda_xsum` for the offset>0 branch is the responsibility-weighted mean
  of the two Exp arguments: `nu * g + (1-nu) * T_g` — this feeds the
  lambda_g M-step correctly.

## Tests

Verify all three K=1 sub-cases analytically:

```r
# offset = 0: eps drops out; tau_sum = 0; loglik = log(lambda_g) - lambda_g*g
test_that("K=1 offset=0: plain Exp, tau_sum=0, eps_informative=FALSE", { ... })

# offset > 0, T_g > 0: 2-pattern mixture; tau_sum = nu; eps_informative=TRUE
test_that("K=1 offset=1: 2-pattern mixture matches analytical", { ... })

# offset > 0, T_g <= 0 (impossible clean): fully contaminated; tau_sum=1
test_that("K=1 offset=1, g < Delta: tau_sum = 1, eps_informative=TRUE", { ... })

# Vectorised (N > 1): all rows independently correct
test_that("K=1 offset>0 vectorized over N rows", { ... })

# Mixed batch (offset=0 and offset>0 in same call): no cross-contamination
test_that("K=1 mixed at_start and offset>0 in same call", { ... })
```

## Prevention

**Rule: always condition on offset when handling K=1 spells.**

Any emission function for a spell model where the spell start may be
unobserved must distinguish:
- `t_observed = spell_start` (offset = 0): T_g is pinned by the observation;
  clean and contaminated patterns may collapse.
- `t_observed > spell_start` (offset > 0): T_g is shifted by the offset;
  the two patterns are structurally different.

The shift `d * Delta` is small (0.25 years for quarterly data) but not
negligible relative to typical tenure values (mean ~0.5 years), so the
approximation error is not numerically insignificant.

**Checklist when writing per-wave emission for a spell model:**
1. For each K, enumerate all 2^K contamination patterns symbolically.
2. Identify which patterns are structurally identical (density is the same)
   and which are distinct.
3. Only collapse identical patterns. Never assume K=1 collapses without
   checking whether the first observed wave is the spell start.

## Related

- [2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md](../testing-patterns/2026-05-01-stale-estep-accumulation-guard-after-emission-refactor.md)
  — Companion fix: after this emission correction, the E-step guard also
  needed updating.
- `documents/EM tenure epsilon.tex`, Section 3.4 (Singletons paragraph) —
  Mathematical justification for the K=1 case split.
- [2026-04-30-na-structural-missingness-vs-error.md](2026-04-30-na-structural-missingness-vs-error.md)
  — Related: distinguishing structural zeros (offset=0, eps drops out) from
  genuine missing information.
