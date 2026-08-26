if(!file.exists("EM-AR1-4W/R/source_all.R"))stop("Run from the project root")
source("EM-AR1-4W/R/source_all.R")
cat("Preparing four-wave panel, covariates, and inconsistency indicators...\n")
data <- prepare_fmm_covariates_inconsistency_4w()
cat(sprintf("Analysis sample: %s people; %s exact likelihood cells\n",
  format(data$n_original,big.mark=","),format(data$n_cells,big.mark=",")))
base_files <- list.files("EM-AR1-4W/output/results",
  pattern="^fmm_ar1_4w_sym_free_.*\\.rds$",full.names=TRUE)
base_fit <- if(length(base_files))readRDS(base_files[which.max(file.info(base_files)$mtime)])else NULL
old_file <- "EM-AR1-4W/output/results/fmm_covariates_inconsistency_4w_latest.rds"
old_fit <- if(file.exists(old_file))readRDS(old_file)else NULL
start <- if(!is.null(old_fit)) .expand_fmm_covinc_start_4w(data,old_fit$params) else
  .initial_fmm_covinc_4w(data,base_fit)
random_starts <- as.integer(Sys.getenv("FMM_RANDOM_STARTS", "6"))
max_iter <- as.integer(Sys.getenv("FMM_MAX_ITER", "600"))
fit <- fit_fmm_covariates_inconsistency_4w(
  data,starts=list(start),random_starts=random_starts,max_iter=max_iter)
summary <- summarize_fmm_covariates_inconsistency_4w(data,fit)
inference <- analytical_se_fmm_covinc_4w(data,fit)
derived_inference <- derived_se_fmm_covinc_4w(data,fit,inference)
outdir <- "EM-AR1-4W/output/results";dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
fit$summary <- summary
fit$analytical_inference <- inference
fit$derived_inference <- derived_inference
fit$model <- "two_type_4w_ar1_common_slopes_duration_job_covariates_inconsistency_pi"
saveRDS(fit,file.path(outdir,"fmm_covariates_inconsistency_4w_latest.rds"))
cat("\nRisk-set weighted transition rates (percent)\n")
print(transform(summary$type_summary,share=100*share,
  initial_employment=100*initial_employment,risk_weighted_entry=100*risk_weighted_entry,
  risk_weighted_exit=100*risk_weighted_exit),row.names=FALSE,digits=5)
cat("\nEntry probit coefficients (type-specific intercepts, common slopes)\n")
print(round(summary$entry_coefficients,5))
cat("\nPersistence probit coefficients (type-specific intercepts, common slopes)\n")
print(round(summary$persistence_coefficients,5))
cat("\nMisclassification equation coefficients and implied probabilities\n")
print(round(summary$misclassification_coefficients,5))
print(transform(summary$misclassification_probabilities,probability=100*probability),
  row.names=FALSE,digits=5)
cat("\nFit diagnostics\n")
print(data.frame(loglik=fit$loglik,converged=fit$converged,iterations=fit$iterations,
  fixed_point_error=fit$fixed_point_error,n=fit$sample$n,cells=fit$sample$cells),row.names=FALSE)
cat("\nTop refined starts\n");print(fit$candidates,row.names=FALSE,digits=10)
cat("\nIdentification and analytical-score diagnostics\n")
print(as.data.frame(inference$diagnostics),row.names=FALSE,digits=6)
cat("\nDerived estimates and analytical delta-method standard errors (percent)\n")
print(transform(derived_inference$summary,estimate=100*estimate,se=100*se),
  row.names=FALSE,digits=5)
