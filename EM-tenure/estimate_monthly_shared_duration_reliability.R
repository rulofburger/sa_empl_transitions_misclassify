# Add a shared person-level duration-reporting reliability mixture to the
# local-anchor model. The two equally weighted latent classes shift both the
# tenure- and timegap-contamination logits by -delta and +delta. The preceding
# model is reproduced exactly at delta=0. The fixed class share and nonnegative
# shift avoid an unnecessary class-weight parameter and label switching.

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

base_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_local_anchor_revision/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate local anchor revisions first")
base <- readRDS(base_file)$extension
base$params$duration_reliability_shift <- 0
base_eval <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_eval$loglik
base$gamma <- base_eval$gamma
base$job_change_posterior <- base_eval$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "shared_duration_reliability")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("SHARED_RELIABILITY_WORKERS","8"))
maxit <- as.integer(Sys.getenv("SHARED_RELIABILITY_MAXIT","2"))
screen_maxit <- as.integer(Sys.getenv("SHARED_RELIABILITY_SCREEN_MAXIT","1"))
resume <- identical(tolower(Sys.getenv("SHARED_RELIABILITY_RESUME","true")),
  "true")

evaluate_shift <- function(shift,reference=base) {
  p <- reference$params
  p$duration_reliability_shift <- shift
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
profile_grid <- c(0,.10,.25,.50,.75,1,1.5,2,3,4,5)
profile_file <- file.path(outdir,"conditional_profile_latest.csv")
if (resume && file.exists(profile_file)) {
  profile <- read.csv(profile_file)
} else {
  profile <- data.frame(reliability_shift=profile_grid,loglik=NA_real_)
  for (j in seq_along(profile_grid)) {
    message("Conditional shared-reliability profile: delta=",profile_grid[j])
    profile$loglik[j] <- evaluate_shift(profile_grid[j])
    write.csv(profile,profile_file,row.names=FALSE)
  }
}
conditional_file <- file.path(outdir,"conditional_fit_latest.rds")
if (resume && file.exists(conditional_file)) {
  conditional_fit <- readRDS(conditional_file)
} else {
  conditional <- optimize(function(shift) -evaluate_shift(shift),
    interval=c(1e-6,5.5),tol=1e-4)
  conditional_params <- base$params
  conditional_params$duration_reliability_shift <- conditional$minimum
  conditional_eval <- e_step_eps(df_fit,conditional_params,check_df=FALSE,
    suff_stats=FALSE)
  conditional_fit <- structure(list(params=conditional_params,
    loglik=-conditional$objective,gamma=conditional_eval$gamma,
    job_change_posterior=conditional_eval$job_change_posterior,
    duration_reliability_posterior=
      conditional_eval$duration_reliability_posterior,
    convergence=0L,iterations=NA_integer_),
    class="conditional_shared_duration_reliability_fit")
  saveRDS(conditional_fit,conditional_file)
}

make_start <- function(shift,label) {
  value <- conditional_fit
  value$params$duration_reliability_shift <- shift
  evaluation <- e_step_eps(df_fit,value$params,check_df=FALSE,
    suff_stats=FALSE)
  value$loglik <- evaluation$loglik
  value$gamma <- evaluation$gamma
  value$job_change_posterior <- evaluation$job_change_posterior
  value$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  value$stage <- label
  value
}
free_names <- c("eps","eps_d","duration_reliability_shift")
fit_stage <- function(label,start,iterations) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Shared-reliability refinement: ",label)
  p <- start$params
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    reliability_shift_start=p$duration_reliability_shift,
    q_start=p$job_change_prob,maxit=iterations,reltol=1e-9,pgtol=1e-7,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
    free_names=free_names,gradient_scheme="forward")
  fit$objective_function <- NULL
  fit$stage <- label
  evaluation <- e_step_eps(df_fit,fit$params,check_df=FALSE,
    suff_stats=FALSE)
  fit$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  saveRDS(fit,path)
  fit
}

main <- fit_stage("conditional",conditional_fit,maxit)
low <- fit_stage("low",make_start(.10,"low"),screen_maxit)
high <- fit_stage("high",make_start(3,"high"),screen_maxit)
candidates <- list(main=main,low=low,high=high)
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

# Re-polish delta exactly at the best estimated baseline contamination rates.
polish_objective <- function(shift) {
  p <- best$params; p$duration_reliability_shift <- shift
  -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
polished_file <- file.path(outdir,"fit_delta_polished_latest.rds")
if (resume && file.exists(polished_file)) {
  polished <- readRDS(polished_file)
} else {
  polish <- optimize(polish_objective,c(1e-6,5.5),tol=1e-5)
  polished_params <- best$params
  polished_params$duration_reliability_shift <- polish$minimum
  polished_eval <- e_step_eps(df_fit,polished_params,check_df=FALSE,
    suff_stats=FALSE)
  polished <- best
  polished$params <- polished_params
  polished$par_unconstrained["duration_reliability_shift"] <- polish$minimum
  polished$loglik <- polished_eval$loglik
  polished$gamma <- polished_eval$gamma
  polished$job_change_posterior <- polished_eval$job_change_posterior
  polished$duration_reliability_posterior <-
    polished_eval$duration_reliability_posterior
  polished$stage <- "delta_polished"
  saveRDS(polished,polished_file)
}
if (polished$loglik>=best$loglik) best <- polished

component_probabilities <- do.call(rbind,lapply(
  c("reliable","unreliable"),function(class) {
    p <- .duration_reliability_component_params(best$params,class)
    data.frame(class=class,prior_probability=.5,
      tenure_contamination=p$eps,timegap_contamination=p$eps_d)
  }))
posterior <- best$duration_reliability_posterior
.weighted_quantile <- function(value,weight,probability) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  value <- value[keep]; weight <- weight[keep]
  ordering <- order(value)
  value <- value[ordering]; weight <- weight[ordering]
  value[findInterval(probability,cumsum(weight)/sum(weight))+1L]
}
posterior_summary <- data.frame(
  statistic=c("weighted mean","minimum","p10","p25","median","p75",
    "p90","maximum"),
  posterior_unreliable=c(weighted.mean(posterior,df_fit$weight),
    min(posterior),vapply(c(.10,.25,.50,.75,.90),function(probability)
      .weighted_quantile(posterior,df_fit$weight,probability),numeric(1)),
    max(posterior)))

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior; p <- fit$params
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=sum(df_fit$weight*jp$expected_changes)/
      sum(df_fit$weight*jp$opportunities),
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    local_anchor_revision_prob=p$tenure_local_revision_prob,
    duration_reliability_shift=p$duration_reliability_shift,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Independent duration reliability",base,31L),
  summarise_fit("Conditional shared reliability only",conditional_fit,32L),
  summarise_fit("Shared duration reliability",best,32L))
refinement <- do.call(rbind,lapply(names(candidates),function(label) {
  f <- candidates[[label]]
  data.frame(stage=label,loglik=f$loglik,gain_over_independent=
      f$loglik-base$loglik,convergence=f$convergence,
    iterations=f$iterations,shift=f$params$duration_reliability_shift,
    tenure_contamination=f$params$eps,
    timegap_contamination=f$params$eps_d)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="No shared duration reliability",LR=LR,df=1L,
  p_value_chisq1_heuristic=pchisq(LR,1,lower.tail=FALSE),
  note="Nonstandard homogeneity boundary; use parametric bootstrap for final inference")
zero_loglik <- evaluate_shift(0)
nesting <- data.frame(independent_loglik=base$loglik,
  zero_shift_loglik=zero_loglik,difference=zero_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("delta=0 does not reproduce the local-anchor model")

cat("\nShared-reliability nesting check\n")
print(nesting,row.names=FALSE,digits=12)
cat("\nConditional profile\n")
print(profile,row.names=FALSE,digits=10)
cat("\nMultistart refinement\n")
print(refinement,row.names=FALSE,digits=10)
cat("\nModel comparison\n")
print(comparison,row.names=FALSE,digits=10)
cat("\nReliability-class contamination probabilities\n")
print(component_probabilities,row.names=FALSE,digits=10)
cat("\nPosterior probability of unreliable class\n")
print(posterior_summary,row.names=FALSE,digits=10)
cat("\nLikelihood-ratio diagnostic\n")
print(lr,row.names=FALSE,digits=10)

write.csv(profile,profile_file,row.names=FALSE)
write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(refinement,file.path(outdir,"refinement_summary_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(component_probabilities,file.path(outdir,
  "component_probabilities_latest.csv"),row.names=FALSE)
write.csv(posterior_summary,file.path(outdir,
  "posterior_unreliable_summary_latest.csv"),row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
saveRDS(list(base=base,conditional=conditional_fit,extension=best,
  fits=candidates,profile=profile,nesting=nesting,refinement=refinement,
  comparison=comparison,component_probabilities=component_probabilities,
  posterior_summary=posterior_summary,lr=lr),
  file.path(outdir,"fits_latest.rds"))
