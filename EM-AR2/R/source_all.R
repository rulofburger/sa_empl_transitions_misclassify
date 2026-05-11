# ==============================================================================
# EM-AR2: Source all module files in correct dependency order
# Created: 2026-05-05
# ==============================================================================
# Usage from the project root:
#   source("EM-AR2/R/source_all.R")
# ==============================================================================

.em_ar2_root <- here::here("EM-AR2")
.em_ar2_r    <- file.path(.em_ar2_root, "R")

# Source in dependency order:
#   utils      — no dependencies
#   latent_histories — depends on utils (.bound01)
#   estep      — depends on utils, latent_histories
#   mstep      — depends on utils
#   em_driver  — depends on utils, latent_histories, estep, mstep
#   inference  — depends on latent_histories, estep (for cell prob computation)
source(file.path(.em_ar2_r, "utils.R"),             local = FALSE)
source(file.path(.em_ar2_r, "latent_histories.R"),   local = FALSE)
source(file.path(.em_ar2_r, "estep.R"),              local = FALSE)
source(file.path(.em_ar2_r, "mstep.R"),              local = FALSE)
source(file.path(.em_ar2_r, "em_driver.R"),              local = FALSE)
source(file.path(.em_ar2_r, "inference.R"),              local = FALSE)
# implied_quantities_AR2 depends on: latent_histories (stationary_ar2), utils (%||%)
source(file.path(.em_ar2_r, "implied_quantities_AR2.R"), local = FALSE)
# bootstrap_utils_AR2 depends on: em_driver (em_fit_ar2), implied_quantities_AR2
# (implied_ar2), and EM-baseline::bootstrap_resample [external — must be
# sourced before this file: source("EM-baseline/R/source_all.R") first].
if (!exists("bootstrap_resample", mode = "function"))
  stop("EM-AR2/R/source_all.R: bootstrap_resample() not found. ",
       "Source EM-baseline/R/source_all.R before EM-AR2/R/source_all.R.")
source(file.path(.em_ar2_r, "bootstrap_utils_AR2.R"),    local = FALSE)

rm(.em_ar2_root, .em_ar2_r)
