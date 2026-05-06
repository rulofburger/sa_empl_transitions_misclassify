# EM-baseline

EM algorithm for three-wave binary employment transitions with misclassification.

**Created**: 2026-05-05  
**Companion paper**: `documents/EM baseline.tex`  
**Related module**: `EM-tenure/` (extends this model with duration emissions)

---

## Overview

This module implements the Expectation–Maximisation (EM) algorithm described
in `EM baseline.tex` for a three-wave panel model of binary employment status
with measurement error.

**True employment** follows a first-order Markov chain with parameters
(θ₀, θ₁, α). **Observed employment** is measured with misclassification
probability π (symmetric) or (π₀, π₁) (asymmetric). The EM algorithm exploits
the finite, enumerable set of 8 latent histories h ∈ {0,1}³ to obtain
**closed-form M-step updates** — no numerical optimisation is needed.

---

## Model Configurations

| Label | `model_type` | `stationary` | Free parameters |
|---|---|---|---|
| `none_stat` | `"none"` | `TRUE` | θ₀, θ₁ |
| `none_free` | `"none"` | `FALSE` | α, θ₀, θ₁ |
| `sym_stat` | `"symmetric"` | `TRUE` | θ₀, θ₁, π |
| `sym_free` | `"symmetric"` | `FALSE` | α, θ₀, θ₁, π |
| `asym_stat` | `"asymmetric"` | `TRUE` | θ₀, θ₁, π₀, π₁ |
| `asym_free` | `"asymmetric"` | `FALSE` | α, θ₀, θ₁, π₀, π₁ |

**Stationary restriction** (TeX Eq 4): α = θ₀ / (θ₀ + 1 − θ₁). When
`stationary = TRUE`, α is derived from the M-step transition estimates.
When `stationary = FALSE`, α is estimated freely from the initial-state
sufficient statistic C₁ / (C₁ + C₀).

---

## File Structure

```
EM-baseline/
├── R/
│   ├── utils.R              # .bound01, .logsumexp, %||%
│   ├── transforms.R         # logit, inv_logit
│   ├── latent_histories.R   # latent_histories(), prior_over_histories(),
│   │                        #   stationary_alpha()           [TeX Eq 6, 4]
│   ├── estep.R              # e_step()                       [TeX Eqs 12-16, 23]
│   ├── mstep.R              # m_step()                       [TeX Eqs 18-19, 23]
│   ├── em_driver.R          # em_fit_baseline(), init_params(),
│   │                        #   compute_observed_loglik()    [TeX Sec 2.4]
│   └── source_all.R         # Sources files in dependency order
├── output/
│   ├── results/             # fit_{label}.rds, run_summary.csv
│   ├── tables/              # table_baseline_results.tex
│   └── figures/             # (reserved for diagnostic plots)
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── helper-source.R
│       ├── test-utils.R
│       ├── test-transforms.R
│       ├── test-latent-histories.R
│       ├── test-estep.R
│       ├── test-mstep.R
│       └── test-em-driver.R
├── estimate_baseline_pipeline.R  # End-to-end orchestrator
└── README.md                     # This file
```

---

## Usage

### Run the full pipeline

```r
# From project root
source("EM-baseline/estimate_baseline_pipeline.R")
```

### Source the module interactively

```r
library(here)
source(here::here("EM-baseline", "R", "source_all.R"))
```

### Fit a single model

```r
source("EM-baseline/R/source_all.R")

# Load data however you like — needs y1, y2, y3 (integer 0/1) and weight
source("scripts/ingest_data_3waves_SA.R")
df <- df_qlfs[, c("y1", "y2", "y3", "weight")]

# Symmetric misclassification, stationary alpha
fit <- em_fit_baseline(
  df         = df,
  model_type = "symmetric",
  stationary = TRUE,
  max_iter   = 1000L,
  tol        = 1e-8,
  verbose    = 1L
)

# Results
fit$params     # final estimates: alpha, theta0, theta1, pi
fit$loglik     # observed-data log-likelihood
fit$converged  # TRUE if converged
fit$iterations # number of EM iterations
fit$gamma      # N x 8 responsibility matrix
```

### Run the tests

```r
# From project root
library(testthat)
test_dir("EM-baseline/tests/testthat")
```

---

## Data Requirements

| Column | Type | Description |
|--------|------|-------------|
| `y1`   | integer (0/1) | Observed employment status at wave 1 |
| `y2`   | integer (0/1) | Observed employment status at wave 2 |
| `y3`   | integer (0/1) | Observed employment status at wave 3 |
| `weight` | numeric > 0 | Survey weight |

---

## Key Function Signatures

```r
em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                params0 = NULL, max_iter = 1000L, tol = 1e-8,
                theta_cap = 0.999, pi_cap = 0.49, verbose = 1L)

e_step(df, params, model_type = "symmetric")
# Returns: list(gamma [N x 8], loglik, suff)

m_step(suff, model_type = "symmetric", stationary = TRUE,
       theta_cap = 0.999, pi_cap = 0.49)
# Returns: list(alpha, theta0, theta1, [pi | pi0 + pi1])

init_params(model_type = "symmetric", stationary = TRUE)
compute_observed_loglik(df, params, model_type = "symmetric")
stationary_alpha(theta0, theta1)
```

---

## Relationship to EM-tenure

The `EM-tenure/` module extends this baseline by adding **duration emission
terms** (tenure and nonemployment duration) to the E-step. The Markov prior
(`prior_over_histories`) and misclassification kernel (`e_step` normalisation)
are structurally identical. Key differences:

| Feature | EM-baseline | EM-tenure |
|---------|------------|---------|
| Duration data | None | Tenure + timegap |
| M-step | Fully closed-form | Closed-form + Brent 1D |
| Emission model | Bernoulli (misclass only) | EMG / Exponential / Gaussian |
| CTMC link | Not applicable | Optional λ ↔ θ link |
| Variants | 6 (2 stat × 3 miscl) | 12 (2 stat × 2 linked × 3 models) |

---

## Numerical Notes

- The **stationarity approximation** (TeX Sec 2.3) updates θ₀, θ₁ from
  transition counts alone, ignoring the initial-state contribution to Q.
  This can cause tiny per-step log-likelihood decreases (< 0.01 per iteration
  in typical data). The free-alpha model (`stationary = FALSE`) satisfies
  strict EM monotonicity.
- **Multiple starts**: the pipeline runs 5 random starts per model and
  retains the solution with the highest log-likelihood, guarding against
  local optima.
- **Boundary constraints**: θ₀, θ₁ ∈ (1e-10, theta_cap), π ∈ [0, pi_cap),
  α ∈ (1e-10, 1 − 1e-10).
