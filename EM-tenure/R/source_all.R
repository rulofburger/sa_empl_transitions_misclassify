# ==============================================================================
# EM-tenure: Source all module files in correct dependency order
# ==============================================================================
# Usage:
#   source("EM-tenure/R/source_all.R")
# ==============================================================================

# NOTE: sys.frame(1)$ofile only works when invoked via source(). If using
# Rscript, ensure the working directory is the project root or EM-tenure/.
.em_tenure_root <- here::here("EM-tenure")

.em_tenure_r <- file.path(.em_tenure_root, "R")

# Source in dependency order
source(file.path(.em_tenure_r, "utils.R"),            local = FALSE)
source(file.path(.em_tenure_r, "transforms.R"),        local = FALSE)
source(file.path(.em_tenure_r, "latent_histories.R"),  local = FALSE)
source(file.path(.em_tenure_r, "emissions.R"),         local = FALSE)
source(file.path(.em_tenure_r, "estep.R"),             local = FALSE)
source(file.path(.em_tenure_r, "mstep.R"),             local = FALSE)
source(file.path(.em_tenure_r, "em_driver.R"),         local = FALSE)
# eps (Spec I) model: depends on emissions.R (interval-censored d emissions),
# mstep.R (.m_step_lambda_d_brent, .m_step_theta0_brent_discrete), and
# latent_histories.R (latent_histories, prior_over_histories).
source(file.path(.em_tenure_r, "emissions_eps.R"),     local = FALSE)
source(file.path(.em_tenure_r, "estep_eps.R"),         local = FALSE)
source(file.path(.em_tenure_r, "mstep_eps.R"),         local = FALSE)
source(file.path(.em_tenure_r, "init_params_eps.R"),   local = FALSE)
source(file.path(.em_tenure_r, "em_driver_eps.R"),     local = FALSE)
source(file.path(.em_tenure_r, "simulate.R"),               local = FALSE)
source(file.path(.em_tenure_r, "diagnostics.R"),            local = FALSE)
source(file.path(.em_tenure_r, "compare_distributions.R"),  local = FALSE)
# Contamination model: implied quantities and bootstrap utilities
source(file.path(.em_tenure_r, "implied_quantities_tenure_contamination.R"),  local = FALSE)
source(file.path(.em_tenure_r, "bootstrap_utils_tenure_contamination.R"),     local = FALSE)

rm(.em_tenure_root, .em_tenure_r)
