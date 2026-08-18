# EM-baseline

EM algorithm for three-wave binary employment transitions with misclassification.

**Created**: 2026-05-05  
**Companion paper**: `documents/EM baseline.tex`  
**Related module**: `EM-tenure/` (extends this model with duration emissions)

---

## Overview

This module implements the three-wave binary employment model with measurement
error. The EM algorithm remains available for latent-history responsibilities
and starting values. Reported estimates use an exact direct maximisation of the
observed likelihood collapsed to the eight possible observed histories.

**True employment** follows a first-order Markov chain with parameters
(θ₀, θ₁, α). **Observed employment** is measured with misclassification
probability π (symmetric) or (π₀, π₁) (asymmetric). The EM algorithm exploits
the finite, enumerable set of 8 latent histories h ∈ {0,1}³ to obtain
For free initial conditions the EM M-step is closed form. Under stationarity,
the initial-state term couples the transition parameters, so the reporting
pipeline uses the exact stationarity-constrained observed-data MLE.

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
│   ├── implied_quantities.R # implied_baseline() — post-estimation transforms
│   ├── mle_baseline.R       # exact eight-cell observed-data MLE + diagnostics
│   ├── analytical_se_baseline.R # weighted sandwich + delta-method inference
│   ├── bootstrap_utils.R    # bootstrap_resample(), bootstrap_one_baseline(),
│   │                        #   summarise_bootstrap()
│   └── source_all.R         # Sources files in dependency order
├── output/
│   ├── results/             # fit_{label}.rds, run_summary.csv
│   ├── tables/              # generated Tables 2 and 3
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
fit <- fit_baseline_mle(
  df         = df,
  model_type = "symmetric",
  stationary = TRUE,
  verbose    = 1L
)

# Results
fit$params     # final estimates: alpha, theta0, theta1, pi
fit$loglik     # observed-data log-likelihood
fit$converged  # TRUE if converged
fit$iterations # number of EM iterations
fit$diagnostics # score, information, and optimizer checks
```

### Compute implied probabilities

```r
# After fitting, compute implied quantities
source("EM-baseline/R/source_all.R")
fit <- readRDS("EM-baseline/output/results/fit_sym_stat.rds")
q   <- implied_baseline(fit$params, model_type = "symmetric")
# q$entry_rate, q$exit_rate, q$employment_rate, q$pi
```

---

## Analytical Standard Errors

`estimate_baseline_pipeline.R` calculates individual-level survey-weighted
sandwich covariance matrices for all six models. Delta-method standard errors
are produced for entry, exit, steady-state employment, and misclassification
probabilities. Results are saved as
`output/results/analytical_se_{label}.rds`. Strata and PSU clustering are not
included because those design identifiers are unavailable in the extract.

## Optional Bootstrap Standard Errors

Bootstrap SEs for all 6 baseline models are computed by `bootstrap_pipeline.R`
(project root). Each model's B=200 results are saved to:

```
EM-baseline/output/results/bootstrap/boot_{label}_B200.rds
```

Each file contains:
- `$boot_results` — B-element list of replicate parameter estimates and implied quantities
- `$summary` — tidy data frame of SEs and 95% CIs for all scalar quantities
- `$n_ok` / `$B` — number of successful reps / total attempted

### Run the bootstrap pipeline

```bash
# From project root (requires point estimates to exist first)
Rscript bootstrap_pipeline.R
```

**To increase B** (e.g., for journal submission): edit `B <- 200L` at the top of
`bootstrap_pipeline.R`, delete the existing per-model `.rds` files, and re-run.

**To re-run a single model**: delete
`EM-baseline/output/results/bootstrap/boot_{label}_B200.rds` and re-run the
pipeline — it skips models whose output file already exists.

### Build LaTeX tables

```bash
Rscript build_tables.R
```

Outputs:
- `EM-baseline/output/tables/table_baseline_params.tex` — parameter estimates with SEs
- `EM-baseline/output/tables/table_baseline_implied.tex` — implied probabilities with SEs

---

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
| Reported estimator | Exact eight-cell observed-data MLE | EM with closed-form + Brent 1D updates |
| Emission model | Bernoulli (misclass only) | EMG / Exponential / Gaussian |
| CTMC link | Not applicable | Optional λ ↔ θ link |
| Variants | 6 (2 stat × 3 miscl) | 12 (2 stat × 2 linked × 3 models) |

---

## Numerical Notes

- The observed likelihood is evaluated using eight survey-weighted history
  cells, which is exactly equivalent to the row-level likelihood.
- Probabilities are optimized on unconstrained logit scales; misclassification
  probabilities use a half-logit transform to remain below 0.5.
- Multiple starts and deterministic nested-model starts guard against local
  optima. The pipeline requires small scores, positive observed information,
  and all model-nesting likelihood inequalities before saving results.
