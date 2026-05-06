# EM-AR2 — Second-Order Markov EM for Misclassified Employment

Created: 2026-05-05

## Overview

This module implements an Expectation-Maximisation (EM) algorithm for the
AR(2) misclassified employment model described in `documents/EM AR2.tex`.
True employment follows a second-order Markov chain; observed employment is
measured with misclassification probability π. Uses 4-wave panel data from
the South Africa QLFS.

## Model Variants

| Variant | Params | Description |
|---------|--------|-------------|
| Symmetric ME | θ₀, θ₀₁, θ₁₀, θ₁, π | Symmetric misclassification |
| No ME | θ₀, θ₀₁, θ₁₀, θ₁ | π fixed at 0 |
| Asymmetric ME | θ₀, θ₀₁, θ₁₀, θ₁, π₀, π₁ | False positive ≠ false negative |

### Parameters

- **θ₀**: Baseline job-finding probability P(y*_t=1 | y*_{t-2}=0, y*_{t-1}=0)
- **θ₀₁**: Duration-dependence increment for job-finding P(y*_t=1 | y*_{t-2}=1, y*_{t-1}=0) − θ₀
- **θ₁**: Baseline separation probability P(y*_t=0 | y*_{t-2}=1, y*_{t-1}=1)
- **θ₁₀**: Duration-dependence increment for separation P(y*_t=0 | y*_{t-2}=0, y*_{t-1}=1) − θ₁
- **π**: P(observed ≠ true)
- **π₀**: P(observed=0 | true=1) (false non-employment)
- **π₁**: P(observed=1 | true=0) (false employment)

## Folder Layout

```
EM-AR2/
├── R/
│   ├── source_all.R        # Sources all modules in order
│   ├── utils.R             # .bound01(), .logsumexp()
│   ├── latent_histories.R  # 16×4 history matrix, stationary dist, prior
│   ├── estep.R             # E-step: responsibilities + sufficient stats
│   ├── mstep.R             # M-step: closed-form parameter updates
│   ├── em_driver.R         # init_params_ar2(), em_fit_ar2() driver
│   └── inference.R         # Cell probs, implied transitions, GoF
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── helper-source.R
│       ├── test-utils.R
│       ├── test-latent-histories.R
│       ├── test-estep.R
│       ├── test-mstep.R
│       ├── test-em-driver.R
│       └── test-inference.R
├── output/
│   ├── results/            # .rds results + run_summary.csv
│   ├── tables/             # LaTeX tables via stargazer
│   └── figures/
├── estimate_pipeline.R     # End-to-end runner
└── README.md
```

## Usage

From the project root:

```r
# Run all three EM variants
source("EM-AR2/estimate_pipeline.R")
```

## Dependencies

- `here` — path management
- `tidyverse` — pipeline filtering and output formatting
- `stargazer` — LaTeX tables (pipeline only)

No external optimisation packages required — M-step is closed-form.

## Methodology

See `documents/EM AR2.tex` for full derivation. Key references:
- Latent histories: 2⁴ = 16 binary sequences over 4 waves
- Prior: AR(2) stationary distribution α(h₁,h₂) with Φ normalisation
- E-step: standard log-sum-exp responsibilities; sufficient stats D_{jk}, T_{jk→l}, M
- M-step: closed-form p̂_{jk→l} = T_{jk→l}/D_{jk}, θ recovery via linear transform
