# EM-tenure

EM algorithm for estimating employment transition probabilities with tenure and nonemployment durations, allowing for symmetric misclassification.

## Model

This implements the **linked specification** from *"A Unified EM Framework for Misclassified Employment Dynamics with Tenure and Nonemployment Durations"* (`documents/EM tenure.tex`).

Two emission modes are supported via the `discrete_timegap` flag (default `TRUE`):

- **Discrete mode** (default): Nonemployment duration D is observed as an integer category code 1–7 corresponding to intervals [0, 3), [3, 6), [6, 9), [9, 12), [12, 36), [36, 60), [60, ∞) months. The model treats the true D as Exp(λ_d) and the observed code as interval-censored. `sigma2_d` is not estimated.
- **Legacy continuous mode** (`discrete_timegap = FALSE`): D is observed as a continuous duration in years. The emission density is EMG(λ_d, σ²_d).

### Parameters

| Parameter | Description | Discrete mode | Legacy mode |
|-----------|-------------|:---:|:---:|
| α | Initial employment probability | ✓ | ✓ |
| θ₁ | P(employed \| employed last period) | ✓ | ✓ |
| θ₀ | P(employed \| nonemployed last period) | ✓ | ✓ |
| π | Symmetric misclassification probability | ✓ | ✓ |
| σ²_g | Employment tenure measurement variance | ✓ | ✓ |
| σ²_d | Nonemployment duration measurement variance | — | ✓ |

The exponential rates λ_g and λ_d are derived from θ₁ and θ₀ via the CTMC link:

```
λ_g = -log(1 - θ₁) / 3
λ_d = -log(1 - θ₀) / 3
```

### Key features

- **Symmetric misclassification**: Estimates π, the probability of observing the wrong employment state
- **Stationarity**: Set `stationary = TRUE` to impose α = θ₀/(θ₀ + 1 - θ₁)
- **Discrete timegap**: Set `discrete_timegap = TRUE` (default) for interval-censored Exp(λ_d); set `FALSE` for legacy continuous EMG
- **Pooled variance estimation**: Combines continuation increments and within-panel starts (tenure only in discrete mode)
- **Linked CTMC**: θ₁, θ₀ updated via Brent's method (joint FOC; full EM, not GEM)

### Timegap category codes

| Code | Interval (months) | Interval (years) |
|:----:|:-----------------:|:----------------:|
| 1 | [0, 3) | [0, 0.25) |
| 2 | [3, 6) | [0.25, 0.50) |
| 3 | [6, 9) | [0.50, 0.75) |
| 4 | [9, 12) | [0.75, 1.00) |
| 5 | [12, 36) | [1.00, 3.00) |
| 6 | [36, 60) | [3.00, 5.00) |
| 7 | [60, ∞) | [5.00, ∞) |

## File structure

```
EM-tenure/
├── R/
│   ├── utils.R              # Constants, .bound01, .logsumexp, timegap helpers
│   ├── transforms.R         # Logit, CTMC link functions
│   ├── latent_histories.R   # 8 latent histories, clocks, prior
│   ├── emissions.R          # EMG, Normal increment, Normal start, discrete interval
│   ├── estep.R              # E-step: responsibilities + sufficient statistics
│   ├── mstep.R              # M-step: Brent theta update + closed-form others
│   ├── em_driver.R          # em_fit_tenure() main driver + init_params()
│   ├── simulate.R           # simulate_panel() synthetic data generator
│   └── source_all.R         # Source all files in dependency order
├── estimate_pipeline.R      # End-to-end estimation script
├── tests/
│   └── testthat/
│       ├── helper-source.R
│       ├── test-utils.R
│       ├── test-transforms.R
│       ├── test-latent_histories.R
│       ├── test-emissions.R
│       ├── test-emissions-discrete.R
│       ├── test-estep-mstep.R
│       └── test-integration.R
├── README.md
└── CONCERNS.md
```

## Usage

### Quick start with synthetic data

```r
source("EM-tenure/R/source_all.R")

# Simulate data (discrete mode, default)
df <- simulate_panel(n = 1000, alpha = 0.6, theta1 = 0.9, theta0 = 0.1,
                     pi = 0.05, sigma2_g = 0.01, seed = 42,
                     discrete_timegap = TRUE)

# Fit model (discrete mode)
fit <- em_fit_tenure(df, discrete_timegap = TRUE, verbose = 2)
fit$params
```

### Legacy continuous mode

```r
df <- simulate_panel(n = 1000, alpha = 0.6, theta1 = 0.9, theta0 = 0.1,
                     pi = 0.05, sigma2_g = 0.01, sigma2_d = 0.01, seed = 42,
                     discrete_timegap = FALSE)
fit <- em_fit_tenure(df, discrete_timegap = FALSE, verbose = 2)
```

### With real data

```r
source("EM-tenure/estimate_pipeline.R")
```

### Running tests

From the project root:

```r
testthat::test_dir("EM-tenure/tests/testthat")
```

## Data requirements

The input data frame must have columns:

| Column | Description | Required |
|--------|-------------|:--------:|
| y1, y2, y3 | Observed employment (0/1) at waves 1-3 | Always |
| tenure1, tenure2, tenure3 | Observed tenure in **years** | Always |
| weight | Survey weight | Always |
| timegap_cat1, timegap_cat2, timegap_cat3 | Nonemployment category codes (1–7) | Discrete mode |
| timegap1, timegap2, timegap3 | Nonemployment duration in **years** | Legacy mode |

### Imputation of wrong-state durations (DIAGNOSIS.md Issue 1)

The EM-tenure model evaluates duration emissions for **both** states at each
wave (the observed state and the hypothetical misclassified state). Raw data
often contain zeros for wrong-state durations (e.g., timegap = 0 for employed
waves). These structural zeros are imputed in the ingest script
(`scripts/ingest_data_3waves_SA.R`) before estimation:

**Timegap category imputation** (for employed waves, where the timegap code is
set to 0 and maps to `NA` after code validation):
- Look for the nearest wave with a valid category code (1–7) for the same individual.
- "Nearest" means |t' − t| = 1 is preferred over |t' − t| = 2.
- **Q6 assumption**: Only valid category codes 1–7 are eligible donors. Never-worked
  waves have `timegap_cat = NA` and are excluded from the donor pool. If an
  individual has no valid donor across all 3 waves (e.g., always employed), assign
  category 1 (`[0, 3)` months) — the floor.
- The never-worked filter in `estimate_pipeline.R` removes fully never-worked
  individuals from estimation, so this floor assumption does not affect EM estimates
  in practice.

**Tenure imputation** (for nonemployed waves, where tenure is set to 0):
- Look for the nearest wave with `tenure > 0` for the same individual.
- If no valid donor: assign `.DURATION_FLOOR = 0.25/12 ≈ 0.021 years ≈ 1 week`.

## TeX reference

All formulas reference the TeX document at `documents/EM tenure.tex`:
- Prior: Eq (1)
- CTMC link: Eq (3)
- Emission densities: Section 2.7 (E-step itemized list)
- Sufficient statistics: Eqs (14)-(18)
- M-step updates: Eqs (19)-(22)
