# Conditional common-dispersion polish followed by a full restricted score
# check. This resolves whether the restricted/unrestricted LR comparison is
# driven by the equality restriction or by an incompletely optimized common
# reliability shift.

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
saved <- readRDS(fit_file)
workers <- as.integer(Sys.getenv("COMMON_POLISH_WORKERS","8"))
resume <- identical(tolower(Sys.getenv("COMMON_POLISH_RESUME","true")),"true")

fit_stage <- function(path,start,free_names,scheme,maxit,label) {
  if (resume && file.exists(path)) return(readRDS(path))
  p <- start$params
  message("Common-reliability polish: ",label)
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    reliability_shift_start=p$duration_reliability_shift,
    q_start=p$job_change_prob,maxit=maxit,reltol=1e-11,pgtol=1e-8,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=2e-4,
    free_names=free_names,gradient_scheme=scheme)
  fit$objective_function <- NULL; fit$stage <- label
  saveRDS(fit,path); fit
}

delta_path <- file.path(outdir,"refine_common_delta_polished_latest.rds")
delta_fit <- fit_stage(delta_path,saved$common,
  "duration_reliability_shift","central",40L,"common_delta_polished")
full_path <- file.path(outdir,"refine_common_full_checked_latest.rds")
full_fit <- fit_stage(full_path,delta_fit,NULL,"forward",1L,
  "common_full_checked")
candidates <- list(saved=saved$common,delta=delta_fit,full=full_fit)
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]
saved$common <- best
saved$fits$common_delta_polished <- delta_fit
saved$fits$common_full_checked <- full_fit

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  p <- fit$params; shifts <- .duration_reliability_shifts(p)
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination_midpoint=p$eps,
    timegap_contamination_midpoint=p$eps_d,job_change_prob=p$job_change_prob,
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    local_anchor_revision_prob=p$tenure_local_revision_prob,
    tenure_reliability_shift=unname(shifts["tenure"]),
    timegap_reliability_shift=unname(shifts["timegap"]),
    convergence=fit$convergence,iterations=fit$iterations)
}
stage_summary <- do.call(rbind,lapply(names(candidates),function(label) {
  fit <- candidates[[label]]
  data.frame(stage=label,loglik=fit$loglik,
    gain_over_saved=fit$loglik-candidates$saved$loglik,
    common_shift=fit$params$duration_reliability_shift,
    convergence=fit$convergence,iterations=fit$iterations)
}))
comparison <- rbind(
  summarise_fit("Common reliability dispersion",best,32L),
  summarise_fit("Separate reliability dispersions",saved$extension,33L))
LR <- max(0,2*(saved$extension$loglik-best$loglik))
lr <- data.frame(comparison=
    "Equal tenure and timegap reliability dispersions",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  note="Equality restriction after conditional common-shift polish and full restricted score check")
saved$comparison <- comparison; saved$lr <- lr
saved$common_polish <- stage_summary
saveRDS(saved,fit_file)
write.csv(stage_summary,file.path(outdir,"common_polish_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
cat("\nCommon-dispersion conditional and full polish\n")
print(stage_summary,row.names=FALSE,digits=12)
cat("\nFinal restricted/unrestricted comparison\n")
print(comparison,row.names=FALSE,digits=12)
cat("\nEquality-restriction test\n")
print(lr,row.names=FALSE,digits=12)
