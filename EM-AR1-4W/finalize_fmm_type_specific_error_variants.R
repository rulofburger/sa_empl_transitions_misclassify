if (!file.exists("EM-AR1-4W/R/source_all.R")) stop("Run from project root")
source("EM-AR1-4W/R/source_all.R")

configs <- list(
  type_error_intercept_only="intercept_only",
  type_error_table3_column1="table3_column1")
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) configs <- configs[names(configs) %in% args]
if (!length(configs)) stop("No matching variant requested")
outdir <- "EM-AR1-4W/output/results"

for (stem in names(configs)) {
  cat("\nFinalizing ",stem,"...\n",sep="")
  data <- prepare_fmm_covariates_inconsistency_4w(error_design=configs[[stem]])
  path <- file.path(outdir,paste0("fmm_",stem,"_4w_latest.rds"))
  fit <- readRDS(path)
  fit$e <- e_step_fmm_covinc_4w(data,fit$params)
  fit$loglik <- fit$e$loglik
  fit$summary <- summarize_fmm_covariates_inconsistency_4w(data,fit)
  fit$analytical_inference <- analytical_se_fmm_covinc_4w(data,fit)
  fit$derived_inference <- derived_se_fmm_covinc_4w(
    data,fit,fit$analytical_inference)
  diag <- fit$analytical_inference$diagnostics
  fit$diagnostically_accepted <- fit$fixed_point_error < 1e-5 &&
    diag$rank == diag$parameters && diag$minimum_eigenvalue > 0
  saveRDS(fit,path)
  write.csv(fit$derived_inference$summary,
    file.path(outdir,paste0("fmm_",stem,"_derived_latest.csv")),row.names=FALSE)
  cat(sprintf("ll=%.3f score=%.3e rank=%d/%d minEig=%.3e condition=%.3g accepted=%s\n",
    fit$loglik,fit$fixed_point_error,diag$rank,diag$parameters,
    diag$minimum_eigenvalue,diag$condition,fit$diagnostically_accepted))
  wanted <- c("type_A_share","type_B_share","A_risk_weighted_entry",
    "B_risk_weighted_entry","A_risk_weighted_exit","B_risk_weighted_exit",
    "A_posterior_mean_misclassification","B_posterior_mean_misclassification",
    "A_initial_employment","B_initial_employment")
  tab <- fit$derived_inference$summary[match(wanted,
    fit$derived_inference$summary$quantity),]
  print(transform(tab,estimate=100*estimate,se=100*se),row.names=FALSE,digits=6)
}
