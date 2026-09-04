# ==============================================================================
# EM-tenure: Tests — helper to source all module files
# ==============================================================================
library(testthat)

# Use here::here() for working-directory-agnostic path resolution (same
# approach as source_all.R). Falls back to normalizePath-based discovery.
if (requireNamespace("here", quietly = TRUE)) {
  .r_dir <- here::here("EM-tenure", "R")
} else {
  .r_dir <- normalizePath(file.path(getwd(), "..", "..", "R"), mustWork = FALSE)
  if (!dir.exists(.r_dir)) {
    .r_dir <- normalizePath("EM-tenure/R", mustWork = FALSE)
  }
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
# eps (Spec I) model files
source(file.path(.r_dir, "emissions_eps.R"))
source(file.path(.r_dir, "tenure_local_gross_eps.R"))
source(file.path(.r_dir, "job_change_eps.R"))
source(file.path(.r_dir, "tenure_monthly_eps.R"))
source(file.path(.r_dir, "estep_eps.R"))
source(file.path(.r_dir, "four_wave_eps.R"))
source(file.path(.r_dir, "timegap_clock_joint.R"))
source(file.path(.r_dir, "piecewise_risk_exact.R"))
source(file.path(.r_dir, "four_wave_fast.R"))
source(file.path(.r_dir, "four_wave_ar2.R"))
source(file.path(.r_dir, "simulate_four_wave_ar2.R"))
source(file.path(.r_dir, "fit_four_wave_ar2.R"))
source(file.path(.r_dir, "mstep_eps.R"))
source(file.path(.r_dir, "init_params_eps.R"))
source(file.path(.r_dir, "em_driver_eps.R"))
source(file.path(.r_dir, "data_prep_eps.R"))
source(file.path(.r_dir, "diagnostics_eps.R"))
source(file.path(.r_dir, "duration_hazard_eps.R"))
source(file.path(.r_dir, "piecewise_hazard_eps.R"))
source(file.path(.r_dir, "simulate.R"))

rm(.r_dir)
