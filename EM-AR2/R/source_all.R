# ==============================================================================
# EM-AR2: Source all module files in correct dependency order
# Created: 2026-05-05
# ==============================================================================
# Usage from the project root:
#   source("EM-AR2/R/source_all.R")
# ==============================================================================

.em_ar2_cursor <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
.em_ar2_project <- NULL
repeat {
  if (file.exists(file.path(.em_ar2_cursor, "EM-AR2", "R", "utils.R"))) {
    .em_ar2_project <- .em_ar2_cursor
    break
  }
  .em_ar2_parent <- dirname(.em_ar2_cursor)
  if (identical(.em_ar2_parent, .em_ar2_cursor)) break
  .em_ar2_cursor <- .em_ar2_parent
}
if (is.null(.em_ar2_project))
  stop("Cannot locate the project root containing EM-AR2/R/utils.R.")
.em_ar2_root <- file.path(.em_ar2_project, "EM-AR2")
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
source(file.path(.em_ar2_r, "analytical_se_AR2.R"),      local = FALSE)
# implied_quantities_AR2 depends on: latent_histories (stationary_ar2), utils (%||%)
source(file.path(.em_ar2_r, "implied_quantities_AR2.R"), local = FALSE)
# Bootstrap helpers are optional. The core estimator and analytical inference
# must remain independently sourceable; bootstrap_pipeline_AR2.R sources the
# baseline resampling helper before requesting these functions.
if (exists("bootstrap_resample", mode = "function")) {
  source(file.path(.em_ar2_r, "bootstrap_utils_AR2.R"), local = FALSE)
}

rm(list = intersect(c(".em_ar2_cursor", ".em_ar2_parent", ".em_ar2_project",
                      ".em_ar2_root", ".em_ar2_r"), ls()))
