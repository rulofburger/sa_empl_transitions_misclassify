# Relax the shared duration-reliability model by allowing separate logit
# dispersions for tenure and timegap contamination. The latent class remains
# person-level and equally weighted. Equality of the two dispersions exactly
# reproduces the preceding common-shift model.

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

common_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "shared_duration_reliability/fits_latest.rds")
if (!file.exists(common_file)) stop("Estimate shared duration reliability first")
common <- readRDS(common_file)$extension
common_shift <- common$params$duration_reliability_shift
common_params <- common$params
common_params$duration_reliability_shift <- 0
common_params$tenure_reliability_shift <- common_shift
common_params$timegap_reliability_shift <- common_shift
common_eval <- e_step_eps(df_fit,common_params,check_df=FALSE,suff_stats=FALSE)
common$params <- common_params
common$loglik <- common_eval$loglik
common$gamma <- common_eval$gamma
common$job_change_posterior <- common_eval$job_change_posterior
common$duration_reliability_posterior <-
  common_eval$duration_reliability_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_WORKERS","8"))
maxit <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_MAXIT","3"))
screen_maxit <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_SCREEN_MAXIT","1"))
resume <- identical(tolower(Sys.getenv("SEPARATE_RELIABILITY_RESUME","true")),
  "true")

evaluate_shifts <- function(shifts,reference=common) {
  p <- reference$params
  p$duration_reliability_shift <- 0
  p$tenure_reliability_shift <- shifts[1]
  p$timegap_reliability_shift <- shifts[2]
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}

# A conditional two-dimensional surface is cheap relative to joint fitting and
# reveals whether the equality restriction is locally plausible.
tenure_grid <- c(.5,1,1.5,2,2.5,3)
timegap_grid <- c(.25,.75,1.25,1.75,2.25,2.75)
profile_file <- file.path(outdir,"conditional_surface_latest.csv")
if (resume && file.exists(profile_file)) {
  surface <- read.csv(profile_file)
} else {
  surface <- expand.grid(tenure_shift=tenure_grid,
    timegap_shift=timegap_grid)
  surface$loglik <- NA_real_
}
for (j in which(!is.finite(surface$loglik))) {
  message("Conditional separate-reliability surface: tenure=",
    surface$tenure_shift[j],", timegap=",surface$timegap_shift[j])
  surface$loglik[j] <- evaluate_shifts(unname(unlist(surface[j,1:2])))
  write.csv(surface,profile_file,row.names=FALSE)
}

conditional_file <- file.path(outdir,"conditional_fit_latest.rds")
if (resume && file.exists(conditional_file)) {
  conditional_fit <- readRDS(conditional_file)
} else {
  conditional_opt <- optim(c(common_shift,common_shift),
    function(shifts) -evaluate_shifts(shifts),method="L-BFGS-B",
    lower=c(1e-6,1e-6),upper=c(5.5,5.5),
    control=list(factr=1e8,pgtol=1e-6,maxit=40))
  conditional_params <- common$params
  conditional_params$tenure_reliability_shift <- conditional_opt$par[1]
  conditional_params$timegap_reliability_shift <- conditional_opt$par[2]
  conditional_eval <- e_step_eps(df_fit,conditional_params,check_df=FALSE,
    suff_stats=FALSE)
  conditional_fit <- structure(list(params=conditional_params,
    loglik=conditional_eval$loglik,gamma=conditional_eval$gamma,
    job_change_posterior=conditional_eval$job_change_posterior,
    duration_reliability_posterior=
      conditional_eval$duration_reliability_posterior,
    convergence=conditional_opt$convergence,
    iterations=unname(conditional_opt$counts["function"]),
    par_unconstrained=setNames(conditional_opt$par,
      c("tenure_reliability_shift","timegap_reliability_shift"))),
    class="conditional_separate_duration_reliability_fit")
  saveRDS(conditional_fit,conditional_file)
}

make_start <- function(tenure_shift,timegap_shift,label) {
  value <- conditional_fit
  value$params$tenure_reliability_shift <- tenure_shift
  value$params$timegap_reliability_shift <- timegap_shift
  evaluation <- e_step_eps(df_fit,value$params,check_df=FALSE,suff_stats=FALSE)
  value$loglik <- evaluation$loglik
  value$gamma <- evaluation$gamma
  value$job_change_posterior <- evaluation$job_change_posterior
  value$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  value$stage <- label
  value
}
free_names <- c("eps","eps_d","tenure_reliability_shift",
  "timegap_reliability_shift")
fit_stage <- function(label,start,iterations) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Separate-reliability refinement: ",label)
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
    q_start=p$job_change_prob,maxit=iterations,reltol=1e-9,pgtol=1e-7,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
    free_names=free_names,gradient_scheme="forward")
  fit$objective_function <- NULL
  fit$stage <- label
  evaluation <- e_step_eps(df_fit,fit$params,check_df=FALSE,suff_stats=FALSE)
  fit$duration_reliability_posterior <-
    evaluation$duration_reliability_posterior
  saveRDS(fit,path)
  fit
}

main <- fit_stage("conditional",conditional_fit,maxit)
common_start <- fit_stage("common",common,screen_maxit)
tenure_high <- fit_stage("tenure_high",
  make_start(2.75,.50,"tenure_high"),screen_maxit)
timegap_high <- fit_stage("timegap_high",
  make_start(.50,2.75,"timegap_high"),screen_maxit)
candidates <- list(main=main,common=common_start,tenure_high=tenure_high,
  timegap_high=timegap_high)
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

# Finish with a two-dimensional conditional polish at the best midpoint
# contamination rates. This makes the equality comparison numerically sharp.
polished_file <- file.path(outdir,"fit_shifts_polished_latest.rds")
if (resume && file.exists(polished_file)) {
  polished <- readRDS(polished_file)
} else {
  polish <- optim(c(best$params$tenure_reliability_shift,
      best$params$timegap_reliability_shift),
    function(shifts) {
      p <- best$params
      p$tenure_reliability_shift <- shifts[1]
      p$timegap_reliability_shift <- shifts[2]
      -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
    },method="L-BFGS-B",lower=c(1e-6,1e-6),upper=c(5.5,5.5),
    control=list(factr=1e7,pgtol=1e-7,maxit=40))
  polished_params <- best$params
  polished_params$tenure_reliability_shift <- polish$par[1]
  polished_params$timegap_reliability_shift <- polish$par[2]
  polished_eval <- e_step_eps(df_fit,polished_params,check_df=FALSE,
    suff_stats=FALSE)
  polished <- best
  polished$params <- polished_params
  polished$par_unconstrained["tenure_reliability_shift"] <- polish$par[1]
  polished$par_unconstrained["timegap_reliability_shift"] <- polish$par[2]
  polished$loglik <- polished_eval$loglik
  polished$gamma <- polished_eval$gamma
  polished$job_change_posterior <- polished_eval$job_change_posterior
  polished$duration_reliability_posterior <-
    polished_eval$duration_reliability_posterior
  polished$stage <- "shifts_polished"
  polished$shift_polish_convergence <- polish$convergence
  polished$shift_polish_iterations <- unname(polish$counts["function"])
  saveRDS(polished,polished_file)
}
# A cached polished fit may predate the explicit distinction between the
# converged conditional polish and the capped joint reporting block.
polished$shift_polish_convergence <- if (is.null(
    polished$shift_polish_convergence)) polished$convergence else
  polished$shift_polish_convergence
polished$shift_polish_iterations <- if (is.null(
    polished$shift_polish_iterations)) polished$iterations else
  polished$shift_polish_iterations
polished$convergence <- best$convergence
polished$iterations <- best$iterations
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
  shifts <- .duration_reliability_shifts(p)
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination_midpoint=p$eps,
    timegap_contamination_midpoint=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=sum(df_fit$weight*jp$expected_changes)/
      sum(df_fit$weight*jp$opportunities),
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    local_anchor_revision_prob=p$tenure_local_revision_prob,
    tenure_reliability_shift=shifts["tenure"],
    timegap_reliability_shift=shifts["timegap"],
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Common reliability dispersion",common,32L),
  summarise_fit("Separate reliability dispersions",best,33L))
refinement <- do.call(rbind,lapply(names(candidates),function(label) {
  f <- candidates[[label]]; shifts <- .duration_reliability_shifts(f$params)
  data.frame(stage=label,loglik=f$loglik,
    gain_over_common=f$loglik-common$loglik,convergence=f$convergence,
    iterations=f$iterations,tenure_shift=shifts["tenure"],
    timegap_shift=shifts["timegap"],
    tenure_contamination_midpoint=f$params$eps,
    timegap_contamination_midpoint=f$params$eps_d)
}))
LR <- max(0,2*(best$loglik-common$loglik))
lr <- data.frame(comparison="Equal tenure and timegap reliability dispersions",
  LR=LR,df=1L,p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  note="Regular equality restriction conditional on the fitted interior mixture and selected likelihood mode")
nesting <- data.frame(common_loglik=readRDS(common_file)$extension$loglik,
  equal_separate_shifts_loglik=common$loglik,
  difference=common$loglik-readRDS(common_file)$extension$loglik)
if (abs(nesting$difference)>1e-6)
  stop("equal separate shifts do not reproduce the common-shift model")

cat("\nSeparate-reliability nesting check\n")
print(nesting,row.names=FALSE,digits=12)
cat("\nConditional surface\n")
print(surface,row.names=FALSE,digits=10)
cat("\nMultistart refinement\n")
print(refinement,row.names=FALSE,digits=10)
cat("\nModel comparison\n")
print(comparison,row.names=FALSE,digits=10)
cat("\nReliability-class contamination probabilities\n")
print(component_probabilities,row.names=FALSE,digits=10)
cat("\nPosterior probability of unreliable class\n")
print(posterior_summary,row.names=FALSE,digits=10)
cat("\nEquality-restriction LR test\n")
print(lr,row.names=FALSE,digits=10)

write.csv(surface,profile_file,row.names=FALSE)
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
saveRDS(list(common=common,conditional=conditional_fit,extension=best,
  fits=candidates,surface=surface,nesting=nesting,refinement=refinement,
  comparison=comparison,component_probabilities=component_probabilities,
  posterior_summary=posterior_summary,lr=lr),
  file.path(outdir,"fits_latest.rds"))
