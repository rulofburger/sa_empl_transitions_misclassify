# ==============================================================================
# EM-baseline: Source all module files in correct dependency order
# Created: 2026-05-05
#
# NOTE: If EM-tenure is also loaded in the same session, sourcing this file
# will overwrite shared symbols (e_step, m_step, init_params, logit, etc.)
# with the EM-baseline versions. Load only one module per session, or call
# functions via a dedicated environment (see architecture review notes).
# ==============================================================================
# Usage (from project root):
#   source("EM-baseline/R/source_all.R")
# ==============================================================================

# Runtime guard: abort if EM-baseline symbols already exist in .GlobalEnv.
# This prevents silent overwrite of EM-tenure (or other module) functions.
#
# PROBLEM: sourcing two modules in the same R session silently overwrites
# shared symbols (e_step, m_step, init_params, logit, ...) with whichever
# module was loaded last, causing wrong results with no error signal.
#
# BYPASS (force reload): remove the conflicting symbol first, then source again:
#   rm(list = c("e_step", "m_step", "init_params"), envir = .GlobalEnv)
#   source("EM-baseline/R/source_all.R")
if (exists("e_step", envir = .GlobalEnv, inherits = FALSE)) {
  stop(paste0(
    "EM-baseline/source_all.R: 'e_step' already exists in .GlobalEnv. ",
    "Another module may be loaded. Restart R or remove the conflicting ",
    "symbols (see BYPASS comment above) to avoid silent overwrite."
  ))
}

.em_baseline_root <- here::here("EM-baseline")
.em_baseline_r    <- file.path(.em_baseline_root, "R")

# Source in dependency order
source(file.path(.em_baseline_r, "utils.R"))
source(file.path(.em_baseline_r, "transforms.R"))
source(file.path(.em_baseline_r, "latent_histories.R"))
source(file.path(.em_baseline_r, "estep.R"))
source(file.path(.em_baseline_r, "mstep.R"))
source(file.path(.em_baseline_r, "em_driver.R"))
source(file.path(.em_baseline_r, "implied_quantities.R"))
source(file.path(.em_baseline_r, "bootstrap_utils.R"))

rm(.em_baseline_root, .em_baseline_r)

# Signal to extension source_all.R that shared utilities are loaded.
.em_baseline_utils_loaded <- TRUE
