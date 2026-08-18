# EM-baseline-ext: Baseline Model Extensions

<!-- last-reviewed: 2026-05-06 -->

## Overview

This subfolder implements three extensions of the baseline symmetric EM model
(`EM-baseline/`) for the 3-wave binary employment panel. Each extension relaxes
a different assumption of the baseline model.

| Extension | TeX Section | Model | Key addition |
|-----------|-------------|-------|-------------|
| I  | Section 5 | Covariate heterogeneity | Probit transition probs via observable covariates |
| III | Section 7 | 2-type FMM | Latent type mixture (stable vs unstable employment) |
| IV | Section 8 | Inconsistency-augmented | Wave-specific π tied to age/education inconsistencies |

The duration covariates here are observed origin-wave controls. A structural
model of endogenous duration dependence remains separate future work.

## Model Descriptions

### Extension I: Observable heterogeneity (covariate model)

Transition probabilities θ₀ᵢ and θ₁ᵢ are individual-specific via probit link
functions: θ₁ᵢₜ = Φ(X₁ᵢₜ β₁), θ₀ᵢₜ = Φ(X₀ᵢₜ β₀). Under stationarity, the initial
employment term couples β₀ and β₁. The M-step therefore applies limited joint
BFGS updates to the full Q-function, with Q-function and observed-likelihood
backtracking safeguards.

Four nested covariate sets:

| Set | Variables | p |
|-----|-----------|---|
| 1 | intercept, age (std), age² (std), educ (std) | 4 |
| 2 | Set 1 + race dummies + female | 4 + r |
| 3 | Set 2 + log tenure in persistence; log time since work and never worked in entry | 12 on the current sample |
| 4 | Set 3 + origin-wave contract type + informal-sector dummy in persistence only | 14 on the current sample |

Sets 1 and 2 use symmetric and no-error models under stationary and free
initial conditions. Sets 3 and 4 use free α because their origin-wave
covariates change between transitions → **12 models** in total.

### Extension III: 2-type finite mixture model (FMM)

Individuals belong to one of two latent employment types (A = stable/high
persistence, B = unstable/low persistence) with unknown mixing weight φ. All
free-initial-condition M-step updates are in closed form. Label switching is
resolved post-convergence by enforcing θ₁ᴬ ≥ θ₁ᴮ. Under stationarity, however,
the initial-state likelihood couples each type's transition parameters, so the
legacy closed-form transition update is not an exact stationary M-step.

Symmetric and no-error model × stationary and free α → **4 models**.

**Identification warning.** The direct eight-cell likelihood audit in
`estimate_fmm_table5.R` finds local Jacobian ranks 4/5 (no error, stationary),
4/6 (symmetric, stationary), 7/7 (no error, free), and 7/8 (symmetric, free).
Thus three binary waves do not identify three of the four Table 5 columns; in
particular, the symmetric-error FMM cannot separately identify π from type
heterogeneity. Numerical values in those columns are points on likelihood
ridges and do not have unique structural interpretations or standard errors.

### Extension IV: Inconsistency-augmented misclassification

Wave-specific misclassification probability:
πᵢₜ = ½ σ(δ₀ + δ₁ Yᵢₜᵃᵍᵉ + δ₂ Yᵢₜᵉᵈᵘ)

where σ(·) is the logistic function, ensuring πᵢₜ ∈ (0, ½). Inconsistency
indicators are computed from age and education reports across waves (see
`compute_inconsistencies.R`). The production Table 6 estimates the exact
observed-data likelihood on collapsed outcome/indicator cells, uses a free
initial employment probability, and reports analytical survey-weighted
sandwich/delta-method standard errors. The EM implementation remains available;
its stationary transition block now jointly maximises the complete-data
objective because stationarity couples the initial probability to both
transition parameters.

Robustness specifications allow the true transition process to differ between
records with and without an inconsistency and allow misclassification to differ
between mild and severe inconsistencies. The indicators are interpreted as
linked-record reliability measures that encompass response and matching errors.

## File Structure

```
EM-baseline-ext/
├── estimate_extensions_pipeline.R   # Main estimation script
├── estimate_analytical_se_table4.R  # Sandwich/delta SEs for Table 4
├── estimate_fmm_table5.R             # Direct-MLE replication + rank audit
├── replicate_table6.R                # Direct MLE, analytical SEs, robustness
├── README.md                        # This file
├── R/
│   ├── source_all.R                 # Sources all extension modules
│   ├── helpers_ext.R                # Shared extension helpers
│   ├── compute_inconsistencies.R    # Wave-specific age/edu inconsistency indicators
│   ├── prepare_covariates.R         # Build covariate design matrices (3 sets)
│   ├── estep_covariates.R           # E-step: covariate extension
│   ├── mstep_covariates.R           # M-step: joint monotone GEM for β
│   ├── em_driver_covariates.R       # Driver + init + loglik: covariate
│   ├── analytical_se_covariates.R   # Weighted sandwich + delta method
│   ├── estep_fmm.R                  # E-step: 2-type FMM
│   ├── mstep_fmm.R                  # M-step: closed form for FMM + label switch
│   ├── em_driver_fmm.R              # Driver + init + loglik: FMM
│   ├── mle_fmm.R                    # Exact eight-cell FMM MLE + rank diagnostics
│   ├── estep_inconsistency.R        # E-step: inconsistency model
│   ├── mstep_inconsistency.R        # M-step: NR GEM for δ
│   ├── em_driver_inconsistency.R    # Driver + init + loglik: inconsistency
│   ├── implied_quantities_ext.R     # implied_covariates(), implied_fmm(),
│   │                                #   implied_inconsistency()
│   └── bootstrap_utils_ext.R        # bootstrap_one_covariates/fmm/inconsistency(),
│                                    #   summarise_bootstrap_ame()
├── tests/
│   └── testthat/
│       ├── helper-source.R          # Sources modules before tests
│       ├── test-compute-inconsistencies.R
│       ├── test-prepare-covariates.R
│       ├── test-covariates.R
│       ├── test-fmm.R
│       └── test-inconsistency.R
└── output/
    ├── results/                     # fit_<label>.rds + run_summary.csv
    └── tables/                      # LaTeX tables (generated by pipeline)
```

## Usage

> ⚠️ **All scripts must be run from the project root directory**, not from within `EM-baseline-ext/`. The pipeline uses `here::here()` for all paths.

Run from the **project root**:

```r
# Interactive
source("EM-baseline-ext/estimate_extensions_pipeline.R")

# Command line
Rscript EM-baseline-ext/estimate_extensions_pipeline.R
```

Each of the 18 model configurations is estimated with `N_STARTS = 5`
random starts. The best-fitting start (highest log-likelihood) is retained.
Results are saved to `output/results/fit_<label>.rds`.

## Running Tests

```r
# From project root (individual test files — see Pester safety rules)
Rscript -e "testthat::test_file('EM-baseline-ext/tests/testthat/test-fmm.R')"
```

## Model Labels

| Label | Extension | Set | Model type | Stationary |
|-------|-----------|-----|-----------|-----------|
| `cov_s1_sym_stat` | I | 1 | symmetric | ✓ |
| `cov_s1_sym_free` | I | 1 | symmetric | |
| `cov_s1_non_stat` | I | 1 | none | ✓ |
| `cov_s1_non_free` | I | 1 | none | |
| `cov_s2_sym_stat` | I | 2 | symmetric | ✓ |
| `cov_s2_sym_free` | I | 2 | symmetric | |
| `cov_s2_non_stat` | I | 2 | none | ✓ |
| `cov_s2_non_free` | I | 2 | none | |
| `cov_s3_sym_stat` | I | 3 | symmetric | ✓ |
| `cov_s3_sym_free` | I | 3 | symmetric | |
| `cov_s3_non_stat` | I | 3 | none | ✓ |
| `cov_s3_non_free` | I | 3 | none | |
| `fmm_sym_stat` | III | — | symmetric | ✓ |
| `fmm_sym_free` | III | — | symmetric | |
| `fmm_non_stat` | III | — | none | ✓ |
| `fmm_non_free` | III | — | none | |
| `incons_sym_stat` | IV | — | symmetric | ✓ |
| `incons_sym_free` | IV | — | symmetric | |

## Bootstrap Standard Errors

Bootstrap SEs for all 18 extension models are computed by `bootstrap_pipeline.R`
(project root). Each model's B=200 results are saved to:

```
EM-baseline-ext/output/results/bootstrap/boot_{label}_B200.rds
```

Each file contains:
- `$boot_results` — B-element list of replicate parameter + implied quantity estimates
- `$summary` — tidy SE table for all scalar quantities
- `$ame_summary` — AME SE table (covariate models only)
- `$n_ok` / `$B` — successful reps / total

### Run the bootstrap pipeline

```bash
# From project root (requires point estimates to exist first)
Rscript bootstrap_pipeline.R
```

**To increase B**: edit `B <- 200L` in `bootstrap_pipeline.R`, delete per-model
`.rds` files, and re-run.

**To re-run a single model**: delete its bootstrap file and re-run — the pipeline
skips models whose output file already exists.

### Build LaTeX tables

```bash
Rscript build_tables.R
```

Outputs (in `EM-baseline-ext/output/tables/`):
- `table_cov_risk_weighted.tex` — main risk-set and survey-weighted transitions
- `table_cov_coefficients_appendix.tex` — raw coefficients and hazard distributions
- `table_cov_ame.tex` — Average Marginal Effects (all 3 sets side-by-side)
- `table_fmm.tex` — FMM implied probabilities
- `table_inconsistency.tex` — inconsistency model results

---

## Dependencies

- R ≥ 4.1
- `here` (path resolution)
- `dplyr` (pipeline data manipulation)
- Shared utilities from `EM-baseline/R/`: `utils.R`, `transforms.R`,
  `latent_histories.R` (sourced automatically)
- `EM-baseline/R/em_driver.R` (for `em_fit_baseline` used in tests)

## TeX References

All model derivations are in `documents/EM baseline.tex`:

- Extension I: Sections 5–6, Eqs (10)–(14)  
- Extension III: Section 7, Eqs (19)–(27)
- Extension IV: Section 8, Eqs (26)–(35)

## Data Availability

The pipeline requires `data/raw/df_qlfs_A.rds` (South Africa QLFS 3-wave balanced
panel) and `data/raw/QLFSmerged_mapped.rds` (the upstream long-format file used
to restore wave-1 sector). These files are **not committed** to the repository
as they contain confidential microdata. To obtain them:

1. Contact the DECDG team at the World Bank for access to the project data share.
2. Place both files in `data/raw/` relative to the project root.
3. The ingest script `scripts/ingest_data_3waves_SA.R` will be called automatically
   by the pipeline.

Unit tests do **not** require these files — they use synthetic panel data generated
by helper functions within each test file.
