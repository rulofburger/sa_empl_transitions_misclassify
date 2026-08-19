# Four-wave, one-type AR(2) model with tenure and non-employment duration.
if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source EM-AR2/estimate_duration_covariates_4w.R from the project root")
source("EM-AR2/R/source_all.R")

data <- prepare_ar2_duration_4w()
message(sprintf("AR(2) duration sample: N = %s; exact cells = %s",
  format(data$n_original, big.mark=","), format(data$n_cells, big.mark=",")))
fit <- fit_ar2_duration_4w(data)
inference <- analytical_se_ar2_duration_4w(data, fit)
fit$analytical_inference <- inference

results_dir <- "EM-AR2/output/results"
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
saveRDS(fit, file.path(results_dir, "ar2_duration_covariates_4w_latest.rds"))
write.csv(inference$coefficient_summary,
          file.path(results_dir, "ar2_duration_coefficients_latest.csv"), row.names=FALSE)
write.csv(inference$probability_summary,
          file.path(results_dir, "ar2_duration_probabilities_latest.csv"), row.names=FALSE)

cat("\nFour-wave one-type AR(2), with observed duration controls\n")
cat(sprintf("N = %s; log likelihood = %.3f; converged = %s; max |score| = %.3g\n\n",
  format(data$n_original,big.mark=","), fit$loglik, fit$converged, fit$max_abs_score))
print(inference$probability_summary, row.names=FALSE, digits=6)
cat("\nTransition coefficients (probit index)\n")
coef_keep <- grepl("^(entry|persistence)_", inference$coefficient_summary$term)
print(inference$coefficient_summary[coef_keep,], row.names=FALSE, digits=6)
cat("\nInformation diagnostics\n")
print(as.data.frame(inference$diagnostics), row.names=FALSE, digits=6)
