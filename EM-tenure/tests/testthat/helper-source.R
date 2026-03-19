# ==============================================================================
# EM-tenure: Tests — helper to source all module files
# ==============================================================================
library(testthat)

# Navigate up from tests/testthat/ to EM-tenure/R/
.r_dir <- normalizePath(file.path(getwd(), "..", "..", "R"), mustWork = FALSE)
if (!dir.exists(.r_dir)) {

  # Fallback: try from project root
  .r_dir <- normalizePath("EM-tenure/R", mustWork = FALSE)
}
if (!dir.exists(.r_dir)) {
  stop("Cannot find EM-tenure/R directory. Run tests from project root or EM-tenure/tests/testthat/.")
}

source(file.path(.r_dir, "utils.R"))
source(file.path(.r_dir, "transforms.R"))
source(file.path(.r_dir, "latent_histories.R"))
source(file.path(.r_dir, "emissions.R"))
source(file.path(.r_dir, "estep.R"))
source(file.path(.r_dir, "mstep.R"))
source(file.path(.r_dir, "em_driver.R"))
source(file.path(.r_dir, "simulate.R"))

rm(.r_dir)
