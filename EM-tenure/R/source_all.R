# ==============================================================================
# EM-tenure: Source all module files in correct dependency order
# ==============================================================================
# Usage:
#   source("EM-tenure/R/source_all.R")
# ==============================================================================

# NOTE: sys.frame(1)$ofile only works when invoked via source(). If using
# Rscript, ensure the working directory is the project root or EM-tenure/.
.em_tenure_root <- normalizePath(file.path(dirname(
  if (sys.nframe() > 0) sys.frame(1)$ofile else "."
), ".."), mustWork = FALSE)

.em_tenure_r <- file.path(.em_tenure_root, "R")

# Source in dependency order
source(file.path(.em_tenure_r, "utils.R"),            local = FALSE)
source(file.path(.em_tenure_r, "transforms.R"),        local = FALSE)
source(file.path(.em_tenure_r, "latent_histories.R"),  local = FALSE)
source(file.path(.em_tenure_r, "emissions.R"),         local = FALSE)
source(file.path(.em_tenure_r, "estep.R"),             local = FALSE)
source(file.path(.em_tenure_r, "mstep.R"),             local = FALSE)
source(file.path(.em_tenure_r, "em_driver.R"),         local = FALSE)
source(file.path(.em_tenure_r, "simulate.R"),          local = FALSE)

rm(.em_tenure_root, .em_tenure_r)
