# METADATA =====================================================================
# DESCRIPTION: Estimate 3-wave AR(1) simple models for SA using the same
#              never-worked filter applied in EM-tenure/estimate_pipeline.R.
#              This removes individuals ever recorded as "never worked" at any
#              wave, raising the employment rate from ~50% to ~70%. The results
#              are directly comparable with the EM-tenure estimates.
# CREATED: 2026-04-14
# BASED ON: scripts/estimate_models_3waves_SA.R (Section 2 only)
#           EM-tenure/estimate_pipeline.R (never-worked filter)

# NOTE: No period filter and no weight normalization are applied, matching
#       the EM-tenure pipeline's data preparation.

#-------------------------------------------------------------------------------
# 1) INITIALISE ----------------------------------------------------------------
#-------------------------------------------------------------------------------

library(tidyverse)
library(haven)
library(data.table)
library(fastverse)

#> Load estimation functions ----
source("scripts/define_estimation_functions_3waves_mle_ar1.R")

#> Define output table helper ----
# Latex output table for simple models without covariates
create_stargazer_table <- function(model_object, df, formula) {
  if (!requireNamespace("maxLik", quietly = TRUE) ||
      !requireNamespace("stargazer", quietly = TRUE)) {
    stop("Required packages 'maxLik' and 'stargazer' are not installed.")
  }

  model_se <- sqrt(((exp(model_object$estimate) / ((1 + exp(model_object$estimate))^2))^2) *
                     diag(vcov(model_object)))

  formula <- as.formula(formula)
  lm_model <- lm(data = df, formula)

  coefficients_transformed <- unlist(logit_inverse(model_object$estimate)) * 100
  std_errors_transformed <- unlist(model_se * 100)

  lm_model$coefficients <- coefficients_transformed
  obs <- nrow(df)
  lm_model$residuals <- rnorm(obs)
  lm_model$se_list <- std_errors_transformed
  lm_model$ll <- model_object$maximum
  return(lm_model)
}

#-------------------------------------------------------------------------------
# 2) INGEST DATA & APPLY NEVER-WORKED FILTER -----------------------------------
#-------------------------------------------------------------------------------

# Run script that loads 3-wave SA data as df_qlfs
source("scripts/ingest_data_3waves_SA.R")

# Convert haven-labelled columns to numeric (identical to EM-tenure pipeline)
df_qlfs <- df_qlfs |>
  mutate(across(where(haven::is.labelled), \(x) as.numeric(x)))

# Remove individuals ever recorded as "never worked" at any wave
# (identical filter to EM-tenure/estimate_pipeline.R lines 29-35)
df_qlfs <- df_qlfs |>
  filter(!(!is.na(neverworked1) & as.numeric(neverworked1) == 1)) |>
  filter(!(!is.na(neverworked2) & as.numeric(neverworked2) == 1)) |>
  filter(!(!is.na(neverworked3) & as.numeric(neverworked3) == 1))

message(sprintf("Observations after never-worked filter: %d", nrow(df_qlfs)))
message(sprintf("Employment rate wave 1: %.1f%%", mean(df_qlfs$y1) * 100))
message(sprintf("Employment rate wave 2: %.1f%%", mean(df_qlfs$y2) * 100))
message(sprintf("Employment rate wave 3: %.1f%%", mean(df_qlfs$y3) * 100))

#-------------------------------------------------------------------------------
# 3) THREE WAVES AR(1) — NO COVARIATES ----------------------------------------
#-------------------------------------------------------------------------------

df_template <- data.table::CJ(y1      = c(0, 1),
                               y1_star = c(0, 1),
                               y2      = c(0, 1),
                               y2_star = c(0, 1),
                               y3      = c(0, 1),
                               y3_star = c(0, 1))
df_estimate <- df_qlfs

#>> No ME ----

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085)
param_init_transformed <- logit_transform(param_init)
model_mle_3w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_pi0,
  grad = calc_mle_derivatives_3waves_ar1_pi0,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_pi0,
  grad = calc_mle_derivatives_3waves_ar1_pi0,
  start = param_init_transformed,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
print(logit_inverse(model_mle_3w_ar1_pi0$estimate))

#>> ME ----

param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, pi = 0.03)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1 <- maxLik::maxLik(
  calc_mle_3waves_ar1,
  grad = calc_mle_derivatives_3waves_ar1,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1 <- maxLik::maxLik(
  calc_mle_3waves_ar1,
  grad = calc_mle_derivatives_3waves_ar1,
  start = model_mle_3w_ar1$estimate,
  method = "NM",
  reltol = 0,
  gradtol = 0
)

print(logit_inverse(model_mle_3w_ar1$estimate))

#>> Asymmetric ME ----

param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, pi_0 = 0.03, pi_1 = 0.03)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1_asymmetric <- maxLik::maxLik(
  calc_mle_3waves_ar1_asymmetric,
  grad = calc_mle_derivatives_3waves_ar1_asymmetric,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_asymmetric <- maxLik::maxLik(
  calc_mle_3waves_ar1_asymmetric,
  grad = calc_mle_derivatives_3waves_ar1_asymmetric,
  start = model_mle_3w_ar1_asymmetric$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
print(logit_inverse(model_mle_3w_ar1_asymmetric$estimate))

#-------------------------------------------------------------------------------
# 4) OUTPUT --------------------------------------------------------------------
#-------------------------------------------------------------------------------

#>> Create stargazer table objects ----

# AR(1) model, 3 waves, no ME
lm_model_3w_ar1_pi0 <- create_stargazer_table(model_object = model_mle_3w_ar1_pi0, df = df_estimate, formula = "y3 ~ y1 + y2")

# AR(1) model, 3 waves, ME
lm_model_3w_ar1 <- create_stargazer_table(model_object = model_mle_3w_ar1, df = df_estimate, formula = "y3 ~ y1 + y2 + y3")

# AR(1) model, 3 waves, Asymmetric ME
lm_model_3w_ar1_asymmetric <- create_stargazer_table(model_object = model_mle_3w_ar1_asymmetric, df = df_estimate, formula = "y3 ~ y1 + y2 + y3 + age1")


# Generate the stargazer table
table_simple_implied <- stargazer::stargazer(
  lm_model_3w_ar1_pi0, lm_model_3w_ar1, lm_model_3w_ar1_asymmetric,
  type = "latex",
  no.space = TRUE,
  label = "table_simple_implied_neverworked_filter",
  digits = 2,
  order = c("theta_1", "theta_10", "theta_0", "theta_01", "pi", "pi_0", "pi_1"),
  se = list(lm_model_3w_ar1_pi0$se_list, lm_model_3w_ar1$se_list, lm_model_3w_ar1_asymmetric$se_list),
  keep.stat = c("n"),
  covariate.labels = c("Entry rate", "Exit rate", "Misclassification rate", "Misclass. rate: non-employed", "Misclass. rate: employed"),
  add.lines = list(
    c("Misclassification", "No", "Symmetric", "Asymmetric"),
    c("LL", round(lm_model_3w_ar1_pi0$ll, 1), round(lm_model_3w_ar1$ll, 1), round(lm_model_3w_ar1_asymmetric$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)

output_file <- "./output/tables/SA/table_simple_implied_neverworked_filter.tex"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
cat(table_simple_implied, file = output_file, sep = "\n")
message(sprintf("LaTeX table saved: %s", output_file))
