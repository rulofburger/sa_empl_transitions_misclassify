# EM-tenure

EM algorithm for estimating employment transition probabilities with tenure and nonemployment durations, allowing for symmetric misclassification.

## Model

This implements the **linked specification** from *"A Unified EM Framework for Misclassified Employment Dynamics with Tenure and Nonemployment Durations"* (`documents/EM tenure.tex`).

### Parameters (6 free)

| Parameter | Description |
|-----------|-------------|
| α | Initial employment probability |
| θ₁ | P(employed \| employed last period) |
| θ₀ | P(employed \| nonemployed last period) |
| π | Symmetric misclassification probability |
| σ²_g | Employment measurement variance |
| σ²_d | Nonemployment measurement variance |

The exponential rates λ_g and λ_d are derived from θ₁ and θ₀ via the CTMC link:

```
λ_g = -log(1 - θ₁) / 3
λ_d = -log(1 - θ₀) / 3
```

### Key features

- **With/without misclassification**: Set `misclassification = TRUE/FALSE`
- **Stationarity**: Set `stationary = TRUE` to impose α = θ₀/(θ₀ + 1 - θ₁)
- **Pooled variance estimation**: Combines continuation increments and within-panel starts
- **Linked CTMC**: θ₁, θ₀ updated via Newton-Raphson (joint FOC; full EM, not GEM)

## File structure

```
EM-tenure/
├── R/
│   ├── utils.R              # Numerical utilities (.bound01, .logsumexp, erfc)
│   ├── transforms.R         # Logit, CTMC link functions
│   ├── latent_histories.R   # 8 latent histories, clocks, prior
│   ├── emissions.R          # EMG, Normal increment, Normal start densities
│   ├── estep.R              # E-step: responsibilities + sufficient statistics
│   ├── mstep.R              # M-step: Newton-Raphson theta update + closed-form others
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
│       ├── test-estep-mstep.R
│       └── test-integration.R
├── README.md
└── CONCERNS.md
```

## Usage

### Quick start with synthetic data

```r
source("EM-tenure/R/source_all.R")

# Simulate data
df <- simulate_panel(n = 1000, alpha = 0.6, theta1 = 0.9, theta0 = 0.1,
                     pi = 0.05, sigma2_g = 0.01, sigma2_d = 0.01, seed = 42)

# Fit with misclassification
fit <- em_fit_tenure(df, misclassification = TRUE, verbose = 2)
fit$params
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

| Column | Description |
|--------|-------------|
| y1, y2, y3 | Observed employment (0/1) at waves 1-3 |
| tenure1, tenure2, tenure3 | Observed tenure in **years** |
| timegap1, timegap2, timegap3 | Observed nonemployment duration in **years** |
| weight | Survey weight |

Duration columns must be in years (divide months by 12 in the ingest script).

## TeX reference

All formulas reference the TeX document at `documents/EM tenure.tex`:
- Prior: Eq (1)
- CTMC link: Eq (3)
- Emission densities: Section 2.7 (E-step itemized list)
- Sufficient statistics: Eqs (14)-(18)
- M-step updates: Eqs (19)-(22)
