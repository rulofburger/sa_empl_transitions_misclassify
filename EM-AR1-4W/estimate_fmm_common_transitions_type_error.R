if (!file.exists("EM-AR1-4W/R/source_all.R")) stop("Run from project root")
source("EM-AR1-4W/R/source_all.R")
# Use the covariate-complete Table 6 sample for exact comparability, but do not
# put any observed covariate in this model's likelihood.
table6_data <- prepare_fmm_covariates_inconsistency_4w(collapse = FALSE)
df <- as.data.frame(table6_data$y)
names(df) <- paste0("y", 1:4)
df$weight <- table6_data$weight
random_starts <- as.integer(Sys.getenv("FMM_RANDOM_STARTS", "120"))

cat("Estimating two-type four-wave FMM with common transitions and type-specific error...\n")
fit <- fit_fmm_ctte_4w(df, random_starts = random_starts)
inference <- analytical_se_fmm_ctte_4w(df, fit)
fit$analytical_inference <- inference

out <- "EM-AR1-4W/output/results/fmm_common_transitions_type_error_4w_latest.rds"
saveRDS(fit, out)
write.csv(inference$summary,
  "EM-AR1-4W/output/results/fmm_common_transitions_type_error_4w_latest.csv",
  row.names = FALSE)

cat("\nEstimated probabilities (percent)\n")
print(transform(inference$summary, estimate = 100 * estimate, se = 100 * se),
      row.names = FALSE, digits = 6)
cat("\nDiagnostics\n")
print(c(loglik = fit$loglik, fit$diagnostics), digits = 7)
