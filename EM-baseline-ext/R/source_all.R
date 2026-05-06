# ==============================================================================
# EM-baseline-ext: Source all module files in correct dependency order
# Created: 2026-05-06
#
# Sources shared utilities from EM-baseline/R/ (utils, transforms,
# latent_histories) then loads extension-specific modules.
#
# IMPORTANT: This file does NOT source EM-baseline's estep.R, mstep.R, or
# em_driver.R. Extensions define their own functions with distinct names
# (e_step_covariates, e_step_fmm, e_step_inconsistency, etc.) to avoid
# symbol collisions in the same R session.
#
# Usage (from project root):
#   source("EM-baseline-ext/R/source_all.R")
# ==============================================================================

# Locate the EM-baseline utilities directory
if (requireNamespace("here", quietly = TRUE)) {
  .baseline_r <- here::here("EM-baseline", "R")
  .ext_r      <- here::here("EM-baseline-ext", "R")
} else {
  .baseline_r <- normalizePath(file.path(getwd(), "EM-baseline", "R"), mustWork = FALSE)
  .ext_r      <- normalizePath(file.path(getwd(), "EM-baseline-ext", "R"), mustWork = FALSE)
}

if (!dir.exists(.baseline_r))
  stop(paste0(
    "EM-baseline-ext/source_all.R: Cannot find EM-baseline/R/ at: ", .baseline_r, "\n",
    "Run from the project root directory."
  ))

if (!dir.exists(.ext_r))
  stop(paste0(
    "EM-baseline-ext/source_all.R: Cannot find EM-baseline-ext/R/ at: ", .ext_r, "\n",
    "Run from the project root directory."
  ))

# ---- 1. Shared utilities from EM-baseline (utils, transforms, latent histories)
# NOTE: We deliberately do NOT source estep.R, mstep.R, or em_driver.R from
# EM-baseline. Extension functions have their own distinct names.
# Guard: skip re-sourcing if already loaded (idempotent when pipeline sources
# EM-baseline/R/source_all.R first).
# Prefer the explicit flag set by EM-baseline's source_all.R; fall back to
# checking for the latent_histories symbol for older sessions that did not use
# the flagged version.
if (!exists(".em_baseline_utils_loaded", envir = .GlobalEnv, inherits = FALSE) &&
    !exists("latent_histories", envir = .GlobalEnv, inherits = FALSE)) {
  source(file.path(.baseline_r, "utils.R"),            local = FALSE)
  source(file.path(.baseline_r, "transforms.R"),       local = FALSE)
  source(file.path(.baseline_r, "latent_histories.R"), local = FALSE)
}

# ---- 2. Shared extension helpers (must precede all extension E-steps)
source(file.path(.ext_r, "helpers_ext.R"),             local = FALSE)

# ---- 3. Extension preprocessing helpers
source(file.path(.ext_r, "compute_inconsistencies.R"), local = FALSE)
source(file.path(.ext_r, "prepare_covariates.R"),      local = FALSE)

# ---- 3. Extension I: Observable heterogeneity (covariates + probit GEM)
source(file.path(.ext_r, "estep_covariates.R"),    local = FALSE)
source(file.path(.ext_r, "mstep_covariates.R"),    local = FALSE)
source(file.path(.ext_r, "em_driver_covariates.R"), local = FALSE)

# ---- 4. Extension III: Unobserved heterogeneity (2-type FMM)
source(file.path(.ext_r, "estep_fmm.R"),    local = FALSE)
source(file.path(.ext_r, "mstep_fmm.R"),    local = FALSE)
source(file.path(.ext_r, "em_driver_fmm.R"), local = FALSE)

# ---- 5. Extension IV: Inconsistency-augmented misclassification (GEM)
source(file.path(.ext_r, "estep_inconsistency.R"),    local = FALSE)
source(file.path(.ext_r, "mstep_inconsistency.R"),    local = FALSE)
source(file.path(.ext_r, "em_driver_inconsistency.R"), local = FALSE)

rm(.baseline_r, .ext_r)
