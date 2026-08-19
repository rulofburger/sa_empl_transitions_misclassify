# EM-tenure

EM algorithm for estimating employment transition probabilities with tenure and nonemployment durations, allowing for symmetric misclassification. This folder implements three model specifications; the **contamination model (epsilon)** is the primary specification of interest.

## Overview for reviewers

If you are reviewing this code alongside the paper, start here:

1. **Theory**: `documents/EM tenure epsilon.tex` — the self-contained companion paper describing the contamination model.
2. **Estimation**: `estimate_pipeline.R` — end-to-end script that ingests data and runs all model variants (base EM-tenure, epsilon). The epsilon section starts at the comment block `# EPSILON-AUGMENTED ESTIMATION`.
3. **Epsilon model code**: The five `*_eps.R` files in `R/` (see file structure below).
4. **Tests**: `tests/testthat/test-eps-model.R`.

## Why the contamination (epsilon) model?

The base EM-tenure model (`documents/EM tenure.tex`) uses Gaussian measurement error (σ²_g) for tenure and an EMG spell-level distribution. Three empirical facts in the QLFS data rule this out:

1. **Tenure is integer-month, not continuous Gaussian.** Reported tenures are discrete month values.
2. **63.2% of E→E continuation increments are exactly 0.25 years.** The Gaussian MLE for σ²_g collapses to zero. The remaining 36.8% include wildly off values (negatives, year-shifted) — a contamination distribution, not Gaussian noise.
3. **36.7% of within-panel new job entrants report tenure > 12 months.** A genuine job entry cannot have 12+ months of tenure — these must be misclassifications or tenure reporting errors.

The base EM-tenure model cannot accommodate these facts: σ²_g collapses to zero, and all inconsistency is routed through π (inflating it to over 25%). The contamination model replaces the Gaussian + EMG measurement model with a **point-mass + Exponential contamination mixture**, eliminating σ²_g entirely, and introduces a **spell-pair joint emission** that uses tenure flanks to identify misclassification.

## Contamination model (epsilon) — specification

### Parameters

| Parameter | Description |
|-----------|-------------|
| α | Initial employment probability |
| θ₁ | P(employed \| employed last period) |
| θ₀ | P(employed \| nonemployed last period) |
| π | Symmetric misclassification probability |
| ε | Per-wave tenure contamination probability |
| λ_g | Exponential rate for employment spells |
| λ_d | Exponential rate for nonemployment spells |

Note: σ²_g is **eliminated** — there is no Gaussian component. The parameter ε replaces it.

### Tenure measurement model

Each observed tenure report is independently either:
- **Clock-consistent** (probability 1 − ε): the observation equals the deterministic clock value exactly (point mass).
- **Contaminated** (probability ε): an independent draw from Exp(λ_g), unrelated to the clock.

For a maximal latent employment spell with K observed tenure waves, the joint tenure emission sums over all 2^K contamination patterns. This is the **spell-pair joint emission** — the key innovation over the base EM-tenure model. When a spell spans an intermediate misclassified wave (e.g., h = (1,1,1) but s = (1,0,1)), the tenure at waves 1 and 3 are "flanks" of the same spell. If both are clock-consistent, they must satisfy g₃ = g₁ + 2Δ. This structural constraint lets tenure directly inform misclassification identification.

### Nonemployment duration

Nonemployment duration D is observed as an integer category code 1–7 (interval-censored). The model treats the true D as Exp(λ_d). This component is shared with the base EM-tenure model.

| Code | Interval (months) | Interval (years) |
|:----:|:-----------------:|:----------------:|
| 1 | [0, 3) | [0, 0.25) |
| 2 | [3, 6) | [0.25, 0.50) |
| 3 | [6, 9) | [0.50, 0.75) |
| 4 | [9, 12) | [0.75, 1.00) |
| 5 | [12, 36) | [1.00, 3.00) |
| 6 | [36, 60) | [3.00, 5.00) |
| 7 | [60, ∞) | [5.00, ∞) |

### Model variants

Each specification is estimated in four variants via two orthogonal switches:

| Variant | Stationarity | CTMC link | Free parameters |
|---------|:---:|:---:|:---:|
| Free, non-stationary | — | — | 7 (α, θ₀, θ₁, π, ε, λ_g, λ_d) |
| Free, stationary | α = θ₀/(θ₀ + 1 − θ₁) | — | 6 |
| Linked, non-stationary | — | λ_g = −log(θ₁)/Δ, λ_d = −log(1−θ₀)/Δ | 5 |
| Linked, stationary | both | both | 4 (θ₀, θ₁, π, ε) |

The **linked + stationary** variant is the most parsimonious. `estimate_pipeline.R` runs all four and performs nested LR tests.

## File structure

```
EM-tenure/
├── R/
│   ├── utils.R              # Constants, .bound01, .logsumexp, timegap helpers
│   ├── transforms.R         # Logit, CTMC link functions
│   ├── latent_histories.R   # 8 latent histories, clocks, prior
│   │
│   │  ── Base EM-tenure model ──
│   ├── emissions.R          # EMG, Normal increment, Normal start, discrete interval
│   ├── estep.R              # E-step (base): responsibilities + sufficient statistics
│   ├── mstep.R              # M-step (base): Brent theta update + closed-form others
│   ├── em_driver.R          # em_fit_tenure() driver + init_params()
│   │
│   │  ── Contamination (epsilon) model ──
│   ├── emissions_eps.R      # Spell-pair point-mass + Exp emission (2^K enumeration)
│   ├── estep_eps.R          # E-step (eps): responsibilities + tau / spell statistics
│   ├── mstep_eps.R          # M-step (eps): closed-form eps, lambda_g; Brent for linked
│   ├── init_params_eps.R    # init_params_eps(): starting values (no sigma2_g)
│   ├── em_driver_eps.R      # em_fit_tenure_eps() driver
│   │
│   │  ── Shared utilities ──
│   ├── simulate.R           # simulate_panel() synthetic data generator
│   ├── diagnostics.R        # Post-estimation: multistart, profile LL, Hessian
│   ├── compare_distributions.R  # Real vs simulated comparison helpers
│   └── source_all.R         # Source all files in dependency order
├── estimate_pipeline.R      # End-to-end estimation (base + eps, all 4 variants each)
├── simulation_diagnostic.R  # Distributional comparison: real vs simulated data
├── tests/
│   └── testthat/
│       ├── helper-source.R
│       ├── test-utils.R
│       ├── test-transforms.R
│       ├── test-latent_histories.R
│       ├── test-emissions.R
│       ├── test-emissions-discrete.R
│       ├── test-estep-mstep.R
│       ├── test-eps-model.R       # ← Tests for the contamination model
│       ├── test-free-spec.R
│       └── test-integration.R
├── CONCERNS.md
├── SIMULATION_DIAGNOSTIC_REPORT.md
└── README.md
```

### Where to look for the epsilon model

| What | File(s) |
|------|---------|
| Theory / derivation | `documents/EM tenure epsilon.tex` |
| Spell-pair emission densities | `R/emissions_eps.R` |
| E-step (responsibilities, tau, sufficient stats) | `R/estep_eps.R` |
| M-step (ε, λ_g closed-form; Brent for linked θ) | `R/mstep_eps.R` |
| Starting values | `R/init_params_eps.R` |
| EM loop / convergence / monotonicity guard | `R/em_driver_eps.R` |
| Estimation pipeline (data → fit → save) | `estimate_pipeline.R` (search for `EPSILON`) |
| Tests | `tests/testthat/test-eps-model.R` |

## Usage

### Epsilon model quick start (synthetic data)

```r
source("EM-tenure/R/source_all.R")

# Initialise parameters
p0 <- init_params_eps(df, linked = FALSE, eps_init = 0.20)

# Fit the contamination model (free, non-stationary)
fit <- em_fit_tenure_eps(df, params0 = p0, stationary = FALSE,
                         linked = FALSE, verbose = 2L)
fit$params
# Returns: alpha, theta0, theta1, pi, eps, lambda_g, lambda_d
```

### With real data (full pipeline)

```r
source("EM-tenure/estimate_pipeline.R")
```

For epsilon point estimates only (three starts per specification, exact-cell
collapsing, weights normalized to the original sample size, and no bootstrap),
run `EM-tenure/estimate_point_tenure_contamination.R`.

This runs all model variants (base EM-tenure × 4 variants, epsilon × 4 variants), saves each fit as a timestamped `.rds` in `output/results/`, appends a summary row to `output/results/run_summary.csv`, and performs LR tests.

### Running tests

From the project root:

```r
testthat::test_dir("EM-tenure/tests/testthat")
```

## Key epsilon model functions

| Function | File | Description |
|----------|------|-------------|
| `em_fit_tenure_eps()` | `em_driver_eps.R` | Top-level EM driver (mirrors `em_fit_tenure()`) |
| `init_params_eps()` | `init_params_eps.R` | Build starting parameter vector (no σ²_g; adds ε) |
| `e_step_eps()` | `estep_eps.R` | E-step: responsibilities γ, spell-pair tau, sufficient stats |
| `m_step_eps()` | `mstep_eps.R` | M-step: closed-form ε, λ_g; Brent solver for linked θ |
| `log_spell_emission_eps()` | `emissions_eps.R` | Spell-pair joint tenure emission (2^K enumeration) |

The point-only runner `estimate_duration_hazard_tenure_contamination.R`
compares the constant-hazard linked epsilon model with the nested hazards
`h(x) = lambda * (1 + x)^beta`. It uses three starts, integrates transition
risk over interval-censored timegap categories, prints duration-specific and
posterior-risk-weighted rates, and does not run the bootstrap.

The follow-on runner `estimate_piecewise_hazard_tenure_contamination.R`
estimates positive piecewise-constant entry and exit hazards over 0--3 months,
3--12 months, 1--3 years, 3--5 years, and 5+ years. The final positive hazard
is maintained over the open tail, guaranteeing a proper distribution and a
finite mean duration. It reuses the corrected constant and power-law fits,
runs multiple starts plus a tighter refinement, and does not bootstrap.

### Key arguments to `em_fit_tenure_eps()`

| Argument | Default | Description |
|----------|---------|-------------|
| `df` | — | Data frame (see data requirements below) |
| `params0` | `NULL` | Starting parameters; if NULL, calls `init_params_eps()` |
| `stationary` | `FALSE` | Impose α = θ₀/(θ₀ + 1 − θ₁) |
| `linked` | `FALSE` | Impose CTMC link: λ_g = −log(θ₁)/Δ, λ_d = −log(1−θ₀)/Δ |
| `max_iter` | `500` | Maximum EM iterations |
| `tol` | `1e-8` | Convergence tolerance (relative ΔLL) |
| `verbose` | `0` | 0 = silent, 1 = final summary, 2 = per-iteration |

## Data requirements

The input data frame must have columns:

| Column | Description | Required |
|--------|-------------|:--------:|
| y1, y2, y3 | Observed employment (0/1) at waves 1-3 | Always |
| tenure1, tenure2, tenure3 | Observed tenure in **years** | Always |
| timegap_cat1, timegap_cat2, timegap_cat3 | Nonemployment category codes (1–7) | Always |
| weight | Survey weight | Always |

### Missing wrong-state durations

The shared ingest script historically fills durations belonging to the
unreported state for compatibility with older estimators. The epsilon model no
longer treats these constructed values as data. `prepare_eps_estimation_data()`
sets tenure to missing at reported-nonemployment waves and timegap to missing
at reported-employment waves. Emissions integrate these unavailable clocks out.
In the duration-dependent transition prior, risk is averaged over the
model-implied duration distribution whenever the latent state differs from the
reported state. Actually observed tenure and timegap reports retain the exact
clock-consistency restrictions.

## TeX references

| Specification | Document |
|---------------|----------|
| Base EM-tenure model | `documents/EM tenure.tex` |
| **Contamination (epsilon) model** | `documents/EM tenure epsilon.tex` |

**Epsilon model TeX cross-references**:
- Contamination model definition: Section 3.2 (Eq. for ε branch)
- Spell-pair joint emission: Section 3.4 (Eq. spell_emission)
- Per-pattern density (clean/contaminated decomposition): Section 3.4 (Eq. per_pattern_density)
- E-step sufficient statistics: Section 4
- M-step updates (ε, λ_g): Section 4.3
