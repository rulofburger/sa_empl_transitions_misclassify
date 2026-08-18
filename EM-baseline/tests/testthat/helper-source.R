# ==============================================================================
# EM-baseline: Tests — helper to source all module files
# Created: 2026-05-05
# ==============================================================================
library(testthat)

if (requireNamespace("here", quietly = TRUE)) {
  .r_dir <- here::here("EM-baseline", "R")
} else {
  .r_dir <- normalizePath(file.path(getwd(), "..", "..", "R"), mustWork = FALSE)
  if (!dir.exists(.r_dir)) {
    .r_dir <- normalizePath("EM-baseline/R", mustWork = FALSE)
  }
}
if (!dir.exists(.r_dir)) {
  stop("Cannot find EM-baseline/R directory. Run tests from project root or EM-baseline/tests/testthat/.")
}

source(file.path(.r_dir, "utils.R"))
source(file.path(.r_dir, "transforms.R"))
source(file.path(.r_dir, "latent_histories.R"))
source(file.path(.r_dir, "estep.R"))
source(file.path(.r_dir, "mstep.R"))
source(file.path(.r_dir, "em_driver.R"))
source(file.path(.r_dir, "implied_quantities.R"))
source(file.path(.r_dir, "mle_baseline.R"))
source(file.path(.r_dir, "analytical_se_baseline.R"))
source(file.path(.r_dir, "bootstrap_utils.R"))

rm(.r_dir)
