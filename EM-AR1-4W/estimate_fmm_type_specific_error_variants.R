if (!file.exists("EM-AR1-4W/R/source_all.R")) stop("Run from the project root")
source("EM-AR1-4W/R/source_all.R")

outdir <- "EM-AR1-4W/output/results"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
old_file <- file.path(outdir, "fmm_type_specific_misclassification_4w_latest.rds")
common_file <- file.path(outdir, "fmm_covariates_inconsistency_4w_latest.rds")
old_fit <- if (file.exists(old_file)) readRDS(old_file) else NULL
common_fit <- if (file.exists(common_file)) readRDS(common_file) else NULL

adapt_start <- function(data, slope_names) {
  base_fit <- old_fit %||% common_fit
  start <- .initial_fmm_covinc_4w(data)
  if (!is.null(base_fit)) {
    p <- base_fit$params
    start$alpha <- p$alpha; start$phi <- p$phi
    for (block in c("beta0", "beta1")) {
      common <- intersect(names(start[[block]]), names(p[[block]]))
      start[[block]][common] <- p[[block]][common]
    }
    old_delta <- p$delta
    intercepts <- if (.fmm_covinc_type_specific_delta(old_delta))
      old_delta[c("intercept_A", "intercept_B")] else rep(old_delta[1], 2L)
  } else {
    intercepts <- rep(start$delta[1], 2L)
    old_delta <- start$delta
  }
  slopes <- setNames(rep(0, length(slope_names)), slope_names)
  common_slopes <- intersect(names(old_delta), slope_names)
  slopes[common_slopes] <- old_delta[common_slopes]
  start$delta <- c(setNames(unname(intercepts), c("intercept_A", "intercept_B")),
                   slopes)
  # With no error covariates, initialize the two intercepts at the posterior
  # mean class-specific rates from the encompassing controlled model. These
  # values only seed the optimizer; both intercepts remain freely estimated.
  if (!length(slope_names))
    start$delta[c("intercept_A", "intercept_B")] <- qlogis(2 * c(.0068, .0460))
  start
}

estimate_variant <- function(name, error_design, random_starts, max_iter) {
  cat(sprintf("\nPreparing %s model...\n", name))
  data <- prepare_fmm_covariates_inconsistency_4w(error_design=error_design)
  slopes <- colnames(data$Z[[1]])[-1L]
  start <- adapt_start(data, slopes)
  cat(sprintf("Sample: %s people; %s cells; %d error slopes\n",
    format(data$n_original,big.mark=","),format(data$n_cells,big.mark=","),length(slopes)))
  fit <- fit_fmm_covariates_inconsistency_4w(data, starts=list(start),
    random_starts=random_starts, max_iter=max_iter, em_warm_start=FALSE)
  summary <- summarize_fmm_covariates_inconsistency_4w(data,fit)
  inference <- analytical_se_fmm_covinc_4w(data,fit)
  derived <- derived_se_fmm_covinc_4w(data,fit,inference)
  fit$summary <- summary; fit$analytical_inference <- inference
  fit$derived_inference <- derived; fit$error_design <- error_design
  fit$model <- paste0("two_type_4w_controlled_type_specific_intercepts_",error_design)
  saveRDS(fit,file.path(outdir,paste0("fmm_",name,"_4w_latest.rds")))
  write.csv(derived$summary,file.path(outdir,paste0("fmm_",name,"_derived_latest.csv")),
            row.names=FALSE)
  cat("\nType-specific transition estimates\n")
  print(transform(summary$type_summary,share=100*share,
    initial_employment=100*initial_employment,
    risk_weighted_entry=100*risk_weighted_entry,
    risk_weighted_exit=100*risk_weighted_exit),row.names=FALSE,digits=6)
  wanted <- c("A_posterior_mean_misclassification",
              "B_posterior_mean_misclassification")
  cat("\nPosterior mean misclassification rates (percent)\n")
  print(transform(derived$summary[match(wanted,derived$summary$quantity),],
                  estimate=100*estimate,se=100*se),row.names=FALSE,digits=6)
  cat("\nMisclassification coefficients\n"); print(round(summary$misclassification_coefficients,6))
  cat("\nDiagnostics\n")
  print(data.frame(loglik=fit$loglik,n=fit$sample$n,cells=fit$sample$cells,
    converged=fit$converged,max_score=fit$fixed_point_error,
    rank=inference$diagnostics$rank,parameters=inference$diagnostics$parameters,
    min_eigenvalue=inference$diagnostics$minimum_eigenvalue,
    condition=inference$diagnostics$condition),row.names=FALSE,digits=7)
  fit
}

random_starts <- as.integer(Sys.getenv("FMM_RANDOM_STARTS", "2"))
max_iter <- as.integer(Sys.getenv("FMM_MAX_ITER", "350"))
intercept_only <- estimate_variant("type_error_intercept_only", "intercept_only",
                                   random_starts, max_iter)
table3_column1 <- estimate_variant("type_error_table3_column1", "table3_column1",
                                   random_starts, max_iter)

cat("\nModel comparison\n")
print(data.frame(model=c("Type-specific intercepts only","Table 3 column (1) slopes"),
  loglik=c(intercept_only$loglik,table3_column1$loglik),
  parameters=c(intercept_only$analytical_inference$diagnostics$parameters,
    table3_column1$analytical_inference$diagnostics$parameters)),row.names=FALSE)
