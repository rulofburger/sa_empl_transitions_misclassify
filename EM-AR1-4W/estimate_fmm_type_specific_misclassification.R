if(!file.exists("EM-AR1-4W/R/source_all.R"))stop("Run from the project root")
source("EM-AR1-4W/R/source_all.R")

cat("Preparing four-wave panel for the type-specific misclassification model...\n")
data <- prepare_fmm_covariates_inconsistency_4w()
common_file <- "EM-AR1-4W/output/results/fmm_covariates_inconsistency_4w_latest.rds"
if(!file.exists(common_file))stop("Estimate the common-intercept model first")
common_fit <- readRDS(common_file)
start <- common_fit$params
start$delta <- setNames(c(start$delta[1],start$delta[1],start$delta[2:3]),
  c("intercept_A","intercept_B","age_inconsistent","education_inconsistent"))

cat(sprintf("Analysis sample: %s people; %s exact likelihood cells\n",
  format(data$n_original,big.mark=","),format(data$n_cells,big.mark=",")))
fit <- fit_fmm_covariates_inconsistency_4w(data,starts=list(start),random_starts=4L)
summary <- summarize_fmm_covariates_inconsistency_4w(data,fit)
inference <- analytical_se_fmm_covinc_4w(data,fit)
derived <- derived_se_fmm_covinc_4w(data,fit,inference)
fit$summary <- summary;fit$analytical_inference <- inference;fit$derived_inference <- derived
fit$model <- "two_type_4w_ar1_common_covariate_slopes_type_specific_pi_intercepts"
fit$comparison <- list(common_loglik=common_fit$loglik,
  loglik_gain=fit$loglik-common_fit$loglik,additional_parameters=1L)
out <- "EM-AR1-4W/output/results/fmm_type_specific_misclassification_4w_latest.rds"
saveRDS(fit,out)

cat("\nRisk-set weighted transition rates (percent)\n")
print(transform(summary$type_summary,share=100*share,
  initial_employment=100*initial_employment,risk_weighted_entry=100*risk_weighted_entry,
  risk_weighted_exit=100*risk_weighted_exit),row.names=FALSE,digits=5)
cat("\nMisclassification coefficients\n");print(round(summary$misclassification_coefficients,6))
cat("\nImplied misclassification probabilities (percent)\n")
print(transform(summary$misclassification_probabilities,probability=100*probability),
  row.names=FALSE,digits=5)
cat("\nDerived estimates and analytical delta-method standard errors (percent)\n")
print(transform(derived$summary,estimate=100*estimate,se=100*se),
  row.names=FALSE,digits=5)
cat("\nFit and identification diagnostics\n")
print(data.frame(loglik=fit$loglik,common_loglik=common_fit$loglik,
  gain=fit$comparison$loglik_gain,converged=fit$converged,
  max_abs_score=fit$fixed_point_error),row.names=FALSE,digits=8)
print(as.data.frame(inference$diagnostics),row.names=FALSE,digits=6)
cat("\nTop refined starts\n");print(fit$candidates,row.names=FALSE,digits=10)
