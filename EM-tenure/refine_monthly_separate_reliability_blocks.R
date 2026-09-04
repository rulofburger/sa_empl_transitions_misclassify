# Test broader parameter blocks around the converged four-parameter
# separate-duration-reliability fit. Stages are deliberately short, chained,
# cached, and resumable. A final full forward-gradient step tests whether the
# blockwise result leaves a material joint direction unexplored.

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
if (!file.exists(fit_file)) stop("Run the four-parameter refinement first")
saved <- readRDS(fit_file)
incumbent <- saved$extension
workers <- as.integer(Sys.getenv("SEPARATE_BLOCK_REFINE_WORKERS","8"))
reporting_maxit <- as.integer(Sys.getenv("SEPARATE_BLOCK_REPORTING_MAXIT","3"))
transition_maxit <- as.integer(Sys.getenv("SEPARATE_BLOCK_TRANSITION_MAXIT","5"))
full_maxit <- as.integer(Sys.getenv("SEPARATE_BLOCK_FULL_MAXIT","1"))
resume <- identical(tolower(Sys.getenv("SEPARATE_BLOCK_REFINE_RESUME","true")),
  "true")

reporting_names <- c("pi","eps","eps_d","job_change","tenure_heaping",
  "tenure_year_revision","tenure_clean_anchor_revision",
  paste0("start_month_logit",1:11),"tenure_exact_anchor_retention",
  "tenure_local_revision","tenure_reliability_shift",
  "timegap_reliability_shift")
transition_names <- c("alpha",paste0("log_hg",1:5),paste0("log_hd",1:5))

fit_stage <- function(label,start,free_names,scheme,iterations) {
  path <- file.path(outdir,paste0("refine_block_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Separate-reliability block refinement: ",label)
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
    q_start=p$job_change_prob,maxit=iterations,reltol=1e-10,pgtol=1e-7,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
    free_names=free_names,gradient_scheme=scheme)
  fit$objective_function <- NULL
  fit$stage <- label
  evaluation <- e_step_eps(df_fit,fit$params,check_df=FALSE,suff_stats=FALSE)
  fit$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  saveRDS(fit,path)
  fit
}

reporting <- fit_stage("reporting",incumbent,reporting_names,"forward",
  reporting_maxit)
transition <- fit_stage("transition",reporting,transition_names,"central",
  transition_maxit)
full <- fit_stage("full",transition,NULL,"forward",full_maxit)
candidates <- list(incumbent=incumbent,reporting=reporting,
  transition=transition,full=full)
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

summarise_stage <- function(label,fit) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  data.frame(stage=label,loglik=fit$loglik,
    gain_over_incumbent=fit$loglik-incumbent$loglik,
    incremental_gain=NA_real_,convergence=fit$convergence,
    iterations=fit$iterations,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,status_misclassification=fit$params$pi,
    tenure_midpoint=fit$params$eps,timegap_midpoint=fit$params$eps_d,
    tenure_shift=fit$params$tenure_reliability_shift,
    timegap_shift=fit$params$timegap_reliability_shift,
    job_change_prob=fit$params$job_change_prob,
    exact_anchor_retention=fit$params$tenure_exact_anchor_retention_prob,
    local_revision=fit$params$tenure_local_revision_prob)
}
block_summary <- do.call(rbind,Map(summarise_stage,names(candidates),candidates))
block_summary$incremental_gain <- c(NA,diff(block_summary$loglik))

if (best$loglik>=saved$extension$loglik-1e-7) {
  saved$fits$broader_block_refinement <- best
  saved$extension <- best
}
saved$block_refinement <- block_summary

component_probabilities <- do.call(rbind,lapply(
  c("reliable","unreliable"),function(class) {
    p <- .duration_reliability_component_params(saved$extension$params,class)
    data.frame(class=class,prior_probability=.5,
      tenure_contamination=p$eps,timegap_contamination=p$eps_d)
  }))
posterior <- saved$extension$duration_reliability_posterior
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
  summarise_fit("Separate reliability dispersions",saved$extension,33L))
LR <- max(0,2*(saved$extension$loglik-saved$common$loglik))
lr <- data.frame(comparison="Equal tenure and timegap reliability dispersions",
  LR=LR,df=1L,p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  note="Equality restriction; final reference distribution remains conditional on a stable selected likelihood mode")
saved$comparison <- comparison
saved$component_probabilities <- component_probabilities
saved$posterior_summary <- posterior_summary
saved$lr <- lr
saveRDS(saved,fit_file)
write.csv(block_summary,file.path(outdir,"block_refinement_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(component_probabilities,file.path(outdir,
  "component_probabilities_latest.csv"),row.names=FALSE)
write.csv(posterior_summary,file.path(outdir,
  "posterior_unreliable_summary_latest.csv"),row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
cat("\nSeparate-reliability broader block refinement\n")
print(block_summary,row.names=FALSE,digits=12)
cat("\nRetained model comparison\n")
print(comparison,row.names=FALSE,digits=12)
cat("\nRetained component probabilities\n")
print(component_probabilities,row.names=FALSE,digits=12)
