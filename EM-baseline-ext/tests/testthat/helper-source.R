# ==============================================================================
# EM-baseline-ext: Tests — helper to source all module files
# Created: 2026-05-06
#
# Sources EM-baseline shared utilities (utils, transforms, latent_histories)
# and all EM-baseline-ext extension modules. Called automatically by testthat
# before running test files in this directory.
# ==============================================================================
library(testthat)

if (requireNamespace("here", quietly = TRUE)) {
  .baseline_r <- here::here("EM-baseline", "R")
  .ext_r      <- here::here("EM-baseline-ext", "R")
} else {
  # Fallback: resolve relative to testthat directory (tests/testthat/ -> project root)
  .root       <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  .baseline_r <- file.path(.root, "EM-baseline", "R")
  .ext_r      <- file.path(.root, "EM-baseline-ext", "R")
}

if (!dir.exists(.baseline_r))
  stop("helper-source.R: Cannot find EM-baseline/R/. Run tests from project root.")
if (!dir.exists(.ext_r))
  stop("helper-source.R: Cannot find EM-baseline-ext/R/. Run tests from project root.")

# Shared baseline utilities + baseline EM driver (needed for em_fit_baseline in tests)
# Note: unlike source_all.R (which skips estep/mstep to avoid symbol collision),
# we deliberately source EM-baseline's estep.R, mstep.R, and em_driver.R here
# because test fixtures compare extension outputs against em_fit_baseline results.
source(file.path(.baseline_r, "utils.R"),            local = FALSE)
source(file.path(.baseline_r, "transforms.R"),       local = FALSE)
source(file.path(.baseline_r, "latent_histories.R"), local = FALSE)
source(file.path(.baseline_r, "estep.R"),            local = FALSE)
source(file.path(.baseline_r, "mstep.R"),            local = FALSE)
source(file.path(.baseline_r, "em_driver.R"),        local = FALSE)

# Extension preprocessing
source(file.path(.ext_r, "compute_inconsistencies.R"), local = FALSE)
source(file.path(.ext_r, "prepare_covariates.R"),      local = FALSE)

# Shared extension helpers (must precede all extension E-steps)
source(file.path(.ext_r, "helpers_ext.R"),             local = FALSE)

# Extension I: Covariates
source(file.path(.ext_r, "estep_covariates.R"),     local = FALSE)
source(file.path(.ext_r, "mstep_covariates.R"),     local = FALSE)
source(file.path(.ext_r, "em_driver_covariates.R"), local = FALSE)

# Extension III: FMM
source(file.path(.ext_r, "estep_fmm.R"),     local = FALSE)
source(file.path(.ext_r, "mstep_fmm.R"),     local = FALSE)
source(file.path(.ext_r, "em_driver_fmm.R"), local = FALSE)
source(file.path(.baseline_r, "mle_baseline.R"), local = FALSE)
source(file.path(.ext_r, "mle_fmm.R"),       local = FALSE)

# Extension IV: Inconsistency
source(file.path(.ext_r, "estep_inconsistency.R"),     local = FALSE)
source(file.path(.ext_r, "mstep_inconsistency.R"),     local = FALSE)
source(file.path(.ext_r, "em_driver_inconsistency.R"), local = FALSE)

# Publication transforms and bootstrap wrappers are tested separately but rely
# on all model families above.
source(file.path(.ext_r, "implied_quantities_ext.R"), local = FALSE)
source(file.path(.baseline_r, "bootstrap_utils.R"), local = FALSE)
source(file.path(.ext_r, "bootstrap_utils_ext.R"), local = FALSE)
source(file.path(.ext_r, "analytical_se_covariates.R"), local = FALSE)
source(file.path(.ext_r, "covariate_reliability.R"), local = FALSE)

rm(.baseline_r, .ext_r)
