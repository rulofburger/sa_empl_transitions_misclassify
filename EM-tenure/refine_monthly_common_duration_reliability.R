# Fair restricted-model refinement for the equal tenure/timegap reliability
# dispersion test. The converged unrestricted nuisance estimates are projected
# onto a common dispersion, then the full 32-parameter restricted likelihood is
# refined in resumable one-step chunks.

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
workers <- as.integer(Sys.getenv("COMMON_JOINT_REFINE_WORKERS","8"))
chunks <- as.integer(Sys.getenv("COMMON_JOINT_REFINE_CHUNKS","4"))
gain_tolerance <- as.numeric(Sys.getenv(
  "COMMON_JOINT_REFINE_GAIN_TOLERANCE","1"))
resume <- identical(tolower(Sys.getenv("COMMON_JOINT_REFINE_RESUME","true")),
  "true")

projected_params <- saved$extension$params
projected_params$duration_reliability_shift <- mean(c(
  projected_params$tenure_reliability_shift,
  projected_params$timegap_reliability_shift))
projected_params$tenure_reliability_shift <- NULL
projected_params$timegap_reliability_shift <- NULL
projected_eval <- e_step_eps(df_fit,projected_params,check_df=FALSE,
  suff_stats=FALSE)
projected <- list(params=projected_params,loglik=projected_eval$loglik,
  gamma=projected_eval$gamma,
  job_change_posterior=projected_eval$job_change_posterior,
  convergence=NA_integer_,iterations=0L,stage="projected_unrestricted")
current <- if (projected$loglik>saved$common$loglik) projected else saved$common
base <- current

for (chunk in seq_len(chunks)) {
  path <- file.path(outdir,sprintf("refine_common_joint_%02d_latest.rds",chunk))
  cached <- resume && file.exists(path)
  if (cached) candidate <- readRDS(path) else {
    p <- current$params
    message("Common-reliability joint refinement chunk ",chunk)
    candidate <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
      heaping_start=max(p$tenure_heaping_prob,1e-8),
      revision_start=p$tenure_year_revision_prob,
      clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
      start_month_probs_start=p$tenure_start_month_probs,
      exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
      local_revision_start=p$tenure_local_revision_prob,
      reliability_shift_start=p$duration_reliability_shift,
      q_start=p$job_change_prob,maxit=1L,reltol=1e-10,pgtol=1e-7,
      workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
      gradient_scheme="forward")
    candidate$objective_function <- NULL
    candidate$stage <- paste0("common_joint_refinement_",chunk)
    saveRDS(candidate,path)
  }
  gain <- candidate$loglik-current$loglik
  if (candidate$loglik>=current$loglik-1e-7) current <- candidate
  saved$common <- current
  saved$fits[[paste0("common_joint_refinement_",chunk)]] <- candidate
  saveRDS(saved,fit_file)
  if (!cached && is.finite(gain) && gain<gain_tolerance) break
}

cached_paths <- file.path(outdir,sprintf(
  "refine_common_joint_%02d_latest.rds",seq_len(chunks)))
cached_paths <- cached_paths[file.exists(cached_paths)]
cached_fits <- lapply(cached_paths,readRDS)
chain_loglik <- c(base$loglik,vapply(cached_fits,`[[`,numeric(1),"loglik"))
history <- data.frame(chunk=0:(length(chain_loglik)-1L),
  start=if (base$loglik==projected$loglik) "projected unrestricted" else
    "saved restricted",loglik=chain_loglik,
  gain=c(NA_real_,diff(chain_loglik)),
  convergence=c(base$convergence,vapply(cached_fits,`[[`,integer(1),
    "convergence")),
  iterations=c(base$iterations,vapply(cached_fits,function(x)
    as.integer(x$iterations),integer(1))))

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  p <- fit$params; shifts <- .duration_reliability_shifts(p)
  jp <- fit$job_change_posterior
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination_midpoint=p$eps,
    timegap_contamination_midpoint=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=sum(df_fit$weight*jp$expected_changes)/
      sum(df_fit$weight*jp$opportunities),
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    local_anchor_revision_prob=p$tenure_local_revision_prob,
    tenure_reliability_shift=unname(shifts["tenure"]),
    timegap_reliability_shift=unname(shifts["timegap"]),
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Common reliability dispersion",current,32L),
  summarise_fit("Separate reliability dispersions",saved$extension,33L))
LR <- max(0,2*(saved$extension$loglik-current$loglik))
lr <- data.frame(comparison=
    "Equal tenure and timegap reliability dispersions",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  note="Equality restriction after symmetric restricted/unrestricted nuisance refinement")
saved$common <- current
saved$comparison <- comparison
saved$lr <- lr
saved$common_joint_refinement <- history
saveRDS(saved,fit_file)
write.csv(history,file.path(outdir,"common_joint_refinement_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
cat("\nCommon-dispersion joint refinement\n")
print(history,row.names=FALSE,digits=12)
cat("\nFair restricted/unrestricted comparison\n")
print(comparison,row.names=FALSE,digits=12)
cat("\nEquality-restriction test\n")
print(lr,row.names=FALSE,digits=12)
