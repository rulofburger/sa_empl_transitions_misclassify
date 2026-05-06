---
date: 2026-04-29
title: "log(erfc(z)) silently underflows for large z — replace with pnorm(log.p=TRUE)"
category: "performance-issues"
language: "R"
tags: [em-algorithm, emg, emissions, numerical-stability, erfc, pnorm, underflow, log-density]
root-cause: "erfc(z) reaches exact 0.0 in double precision for z > ~27, so log(erfc(z)) returns -Inf silently, collapsing responsibilities for short-tenure observations"
severity: "P1"
---

# `log(erfc(z))` Silently Underflows for Large z in EMG Density

## Problem

The ExGaussian (EMG) log-density was implemented as:

```r
log_density <- log(lambda) + (lambda^2 * sigma2 / 2) - lambda * x +
  log(erfc(z)) - log(2)
```

For large positive `z` — which occurs for short spells (small `x`) or high
`lambda` — `erfc(z)` underflows to exactly `0.0` in double precision *before*
`log()` is applied. The result is `-Inf`, not an approximation error but a
silent exact collapse. This zeroes out the responsibility `gamma_{ih}` for
the affected history, corrupting the E-step.

**Symptom**: Responsibilities can collapse to 0 for certain observation-history
combinations at realistic parameter values (e.g., high retention rates and
short observed tenures), causing the EM algorithm to produce incorrect parameter
estimates without error.

## Root Cause

`erfc(z)` uses the complementary error function which, for `z > ~27`, is below
double-precision machine epsilon (~2.2e-16) and rounds to 0. The R `erfc()`
function (from `pracma`) does not have a log-scale variant.

The identity `log(erfc(z)/2) = pnorm(-z*sqrt(2), log.p=TRUE)` provides the
numerically stable alternative: `pnorm` with `log.p=TRUE` uses the
Wichura-Hill algorithm which computes the log of the tail probability
directly without first computing the probability itself.

## Solution

Replace the erfc-based computation with `pnorm(log.p=TRUE)`:

```r
log_emg <- function(x, lambda, sigma2) {
  sigma <- sqrt(.pos(sigma2))
  z     <- (lambda * sigma2 - x) / (sqrt(2) * sigma)

  log_density <- log(lambda) +
    (lambda^2 * sigma2 / 2) -
    lambda * x +
    pnorm(-z * sqrt(2), log.p = TRUE)   # replaces log(erfc(z)) - log(2)

  ifelse(x > 0, log_density, -Inf)
}
```

**Identity used**: `log(erfc(z)/2) = pnorm(-z*sqrt(2), log.p=TRUE)`

The gradient function `log_emg_grad_lambda` already used `pnorm(log.p=TRUE)`
for the inverse Mills ratio — this fix makes the density consistent with it.

## Prevention

- **Always use `log.p=TRUE` variants of distribution functions** when computing
  log-densities or log-probabilities. Never compute `log(p_fun(...))` — the
  `p_fun` may underflow to 0 first.
- Equivalences to remember:
  - `log(pnorm(u)) = pnorm(u, log.p=TRUE)`
  - `log(erfc(z)/2) = pnorm(-z*sqrt(2), log.p=TRUE)`
  - `log(1 - pnorm(u)) = pnorm(u, lower.tail=FALSE, log.p=TRUE)`
- The same pattern applies in `log_emg` calls for `lambda_d` (nonemployment).

## Related

- [2026-03-15-vectorise-em-estep-over-observations.md](2026-03-15-vectorise-em-estep-over-observations.md) — vectorisation of E-step
- [2026-03-20-em-monotonicity-guard.md](../testing-patterns/2026-03-20-em-monotonicity-guard.md) — EM numerical stability checks
