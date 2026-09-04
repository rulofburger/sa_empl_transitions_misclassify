# Continue the Panel A separate-dispersion model until the four-parameter
# contamination/reliability block converges. Completed stages are cached.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month",1:3))

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability")
fit_file <- file.path(outdir,"fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate separate duration reliability first")
saved <- readRDS(fit_file)
incumbent <- saved$extension
workers <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_REFINE_WORKERS","8"))
maxit <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_REFINE_MAXIT","15"))
chunks <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_REFINE_CHUNKS","3"))
resume <- identical(tolower(Sys.getenv("SEPARATE_RELIABILITY_REFINE_RESUME",
  "true")),"true")
free_names <- c("eps","eps_d","tenure_reliability_shift",
  "timegap_reliability_shift")

fit_stage <- function(label,start,scheme="central",iterations=maxit) {
  path <- file.path(outdir,paste0("refine_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Separate-reliability longer refinement: ",label)
  p <- start$params
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    tenure_reliability_shift_start=p$tenure_reliability_shift,
    timegap_reliability_shift_start=p$timegap_reliability_shift,
    q_start=p$job_change_prob,maxit=iterations,reltol=1e-11,pgtol=1e-8,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=2e-4,
    free_names=free_names,gradient_scheme=scheme)
  fit$objective_function <- NULL
  fit$stage <- label
  evaluation <- e_step_eps(df_fit,fit$params,check_df=FALSE,suff_stats=FALSE)
  fit$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  saveRDS(fit,path)
  fit
}

central <- fit_stage("four_parameter_central",incumbent)
central_fits <- list(central_1=central)
if (chunks>1L) for (chunk in 2:chunks) {
  previous <- central_fits[[length(central_fits)]]
  if (identical(previous$convergence,0L)) break
  central_fits[[paste0("central_",chunk)]] <- fit_stage(
    paste0("four_parameter_central_",chunk),previous)
}
candidates <- c(list(incumbent=incumbent),central_fits)
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

summary <- do.call(rbind,lapply(names(candidates),function(label) {
  f <- candidates[[label]]
  data.frame(stage=label,loglik=f$loglik,
    gain_over_incumbent=f$loglik-incumbent$loglik,
    convergence=f$convergence,iterations=f$iterations,
    tenure_midpoint=f$params$eps,timegap_midpoint=f$params$eps_d,
    tenure_shift=f$params$tenure_reliability_shift,
    timegap_shift=f$params$timegap_reliability_shift)
}))
if (best$loglik>=saved$extension$loglik-1e-7) {
  saved$fits$long_four_parameter_refinement <- best
  saved$extension <- best
}
saved$long_refinement <- summary
saveRDS(saved,fit_file)
write.csv(summary,file.path(outdir,
  "long_four_parameter_refinement_latest.csv"),row.names=FALSE)
cat("\nSeparate-reliability longer four-parameter refinement\n")
print(summary,row.names=FALSE,digits=12)
cat("\nRetained stage\n")
print(summary[which.max(summary$loglik),],row.names=FALSE,digits=12)
