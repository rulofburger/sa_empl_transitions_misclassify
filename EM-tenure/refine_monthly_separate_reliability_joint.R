# Resumable one-step joint continuation for the 33-parameter model with
# separate tenure and timegap reliability dispersions. Short cached chunks
# avoid losing a long numerical gradient if a session is interrupted.

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
current <- saved$extension
base_loglik <- if (!is.null(saved$block_refinement) &&
    any(saved$block_refinement$stage=="full"))
  saved$block_refinement$loglik[saved$block_refinement$stage=="full"][1L] else
  current$loglik
workers <- as.integer(Sys.getenv("SEPARATE_JOINT_REFINE_WORKERS","8"))
chunks <- as.integer(Sys.getenv("SEPARATE_JOINT_REFINE_CHUNKS","3"))
gain_tolerance <- as.numeric(Sys.getenv(
  "SEPARATE_JOINT_REFINE_GAIN_TOLERANCE","1"))
resume <- identical(tolower(Sys.getenv("SEPARATE_JOINT_REFINE_RESUME","true")),
  "true")
history <- data.frame(chunk=0L,loglik=base_loglik,gain=NA_real_,
  convergence=current$convergence,iterations=current$iterations)

for (chunk in seq_len(chunks)) {
  path <- file.path(outdir,sprintf("refine_joint_%02d_latest.rds",chunk))
  prior_loglik <- current$loglik
  cached <- resume && file.exists(path)
  if (cached) {
    candidate <- readRDS(path)
  } else {
    p <- current$params
    message("Separate-reliability joint continuation chunk ",chunk)
    candidate <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
      heaping_start=max(p$tenure_heaping_prob,1e-8),
      revision_start=p$tenure_year_revision_prob,
      clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
      start_month_probs_start=p$tenure_start_month_probs,
      exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
      local_revision_start=p$tenure_local_revision_prob,
      tenure_reliability_shift_start=p$tenure_reliability_shift,
      timegap_reliability_shift_start=p$timegap_reliability_shift,
      q_start=p$job_change_prob,maxit=1L,reltol=1e-10,pgtol=1e-7,
      workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
      gradient_scheme="forward")
    candidate$objective_function <- NULL
    candidate$stage <- paste0("joint_continuation_",chunk)
    candidate$duration_reliability_posterior <- e_step_eps(df_fit,
      candidate$params,check_df=FALSE,suff_stats=FALSE)$
        duration_reliability_posterior
    saveRDS(candidate,path)
  }
  gain <- candidate$loglik-prior_loglik
  history <- rbind(history,data.frame(chunk=chunk,loglik=candidate$loglik,
    gain=gain,convergence=candidate$convergence,
    iterations=candidate$iterations))
  if (candidate$loglik>=current$loglik-1e-7) current <- candidate
  saved$extension <- current
  saved$fits[[paste0("joint_continuation_",chunk)]] <- candidate
  saved$joint_continuation <- history
  saveRDS(saved,fit_file)
  # A cached chunk may equal the incumbent because the preceding run already
  # saved it. Continue to the first missing chunk rather than treating that
  # zero bookkeeping gain as numerical convergence.
  if (!cached && is.finite(gain) && gain<gain_tolerance) break
}

# Reconstruct the reported chain from the cached files. This keeps gains
# meaningful when a completed script is rerun from its final incumbent.
cached_paths <- file.path(outdir,sprintf("refine_joint_%02d_latest.rds",
  seq_len(chunks)))
cached_paths <- cached_paths[file.exists(cached_paths)]
if (length(cached_paths)) {
  cached_fits <- lapply(cached_paths,readRDS)
  chain_loglik <- c(base_loglik,
    vapply(cached_fits,`[[`,numeric(1),"loglik"))
  history <- data.frame(chunk=0:(length(chain_loglik)-1L),
    loglik=chain_loglik,gain=c(NA_real_,diff(chain_loglik)),
    convergence=c(NA_integer_,vapply(cached_fits,`[[`,integer(1),
      "convergence")),
    iterations=c(NA_integer_,vapply(cached_fits,function(x)
      as.integer(x$iterations),integer(1))))
}

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
  summarise_fit("Common reliability dispersion",saved$common,32L),
  summarise_fit("Separate reliability dispersions",current,33L))
component_probabilities <- do.call(rbind,lapply(c("reliable","unreliable"),
  function(class) {
    p <- .duration_reliability_component_params(current$params,class)
    data.frame(class=class,prior_probability=.5,
      tenure_contamination=p$eps,timegap_contamination=p$eps_d)
  }))
posterior <- current$duration_reliability_posterior
if (is.null(posterior)) posterior <- e_step_eps(df_fit,current$params,
  check_df=FALSE,suff_stats=FALSE)$duration_reliability_posterior
.weighted_quantile <- function(value,weight,probability) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  value <- value[keep]; weight <- weight[keep]
  ordering <- order(value); value <- value[ordering]; weight <- weight[ordering]
  value[findInterval(probability,cumsum(weight)/sum(weight))+1L]
}
posterior_summary <- data.frame(
  statistic=c("weighted mean","minimum","p10","p25","median","p75",
    "p90","maximum"),
  posterior_unreliable=c(weighted.mean(posterior,df_fit$weight),min(posterior),
    vapply(c(.10,.25,.50,.75,.90),function(probability)
      .weighted_quantile(posterior,df_fit$weight,probability),numeric(1)),
    max(posterior)))
LR <- max(0,2*(current$loglik-saved$common$loglik))
lr <- data.frame(comparison=
    "Equal tenure and timegap reliability dispersions",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  note="Equality restriction; chi-square reference remains conditional on a stable selected likelihood mode")
saved$comparison <- comparison
saved$component_probabilities <- component_probabilities
saved$posterior_summary <- posterior_summary
saved$lr <- lr
saved$joint_continuation <- history
saveRDS(saved,fit_file)
write.csv(history,file.path(outdir,"joint_continuation_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(component_probabilities,file.path(outdir,
  "component_probabilities_latest.csv"),row.names=FALSE)
write.csv(posterior_summary,file.path(outdir,
  "posterior_unreliable_summary_latest.csv"),row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
cat("\nSeparate-reliability joint continuation\n")
print(history,row.names=FALSE,digits=12)
cat("\nRetained model comparison\n")
print(comparison,row.names=FALSE,digits=12)
cat("\nRetained component probabilities\n")
print(component_probabilities,row.names=FALSE,digits=12)
