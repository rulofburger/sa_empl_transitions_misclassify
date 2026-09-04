# Add an exact reported-start-date retention channel to the flexible baseline
# calendar model. Conditional on a gross current tenure report, rho is the
# probability of retaining the preceding reported start date exactly. The
# existing whole-year revision and redraw probabilities apply conditional on
# not retaining it, so the preceding model is nested at rho=0.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
interview_cols <- paste0("interview_month",1:3)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=interview_cols)

base_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_start_month_baseline/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the flexible baseline model first")
base <- readRDS(base_file)$seasonal
base$params$tenure_exact_anchor_retention_prob <- 0
base_eval <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_eval$loglik
base$gamma <- base_eval$gamma
base$job_change_posterior <- base_eval$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("EXACT_ANCHOR_WORKERS","8"))
maxit <- as.integer(Sys.getenv("EXACT_ANCHOR_MAXIT","8"))
resume <- identical(tolower(Sys.getenv("EXACT_ANCHOR_RESUME","true")),
  "true")
total_weight <- sum(df_fit$weight)

# First profile only the new probability. This both diagnoses identification
# and supplies a stable start for joint optimization.
profile_grid <- unique(c(0,seq(.01,.20,by=.01),seq(.25,.80,by=.05)))
evaluate_rho <- function(rho) {
  p <- base$params
  p$tenure_exact_anchor_retention_prob <- rho
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
profile_file <- file.path(outdir,"conditional_profile_latest.csv")
if (resume && file.exists(profile_file)) {
  profile <- read.csv(profile_file)
} else {
  profile <- data.frame(rho=profile_grid,loglik=NA_real_)
  for (j in seq_along(profile_grid)) {
    message("Conditional exact-anchor profile: rho=",profile_grid[j])
    profile$loglik[j] <- evaluate_rho(profile_grid[j])
    write.csv(profile,profile_file,row.names=FALSE)
  }
}
conditional <- optimize(function(rho) -evaluate_rho(rho),
  interval=c(1e-7,.95),tol=1e-6)
conditional_rho <- conditional$minimum
conditional_loglik <- -conditional$objective
conditional_params <- base$params
conditional_params$tenure_exact_anchor_retention_prob <- conditional_rho
conditional_eval <- e_step_eps(df_fit,conditional_params,check_df=FALSE,
  suff_stats=FALSE)

fit_joint <- function(label,start_params,scheme,iterations) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Joint exact-anchor refinement: ",label)
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start_params,
    heaping_start=start_params$tenure_heaping_prob,
    revision_start=start_params$tenure_year_revision_prob,
    clean_anchor_revision_start=
      start_params$tenure_clean_anchor_revision_prob,
    start_month_probs_start=start_params$tenure_start_month_probs,
    exact_anchor_retention_start=
      start_params$tenure_exact_anchor_retention_prob,
    q_start=start_params$job_change_prob,maxit=iterations,reltol=1e-9,
    pgtol=1e-7,workers=max(1L,min(workers,8L)),verbose=1L,
    gradient_step=1e-3,gradient_scheme=scheme)
  fit$objective_function <- NULL
  saveRDS(fit,path)
  fit
}
fit_forward <- fit_joint("full_forward",conditional_params,"forward",maxit)
fit_central <- fit_joint("full_central",fit_forward$params,"central",
  min(maxit,3L))
joint_fits <- list(forward=fit_forward,central=fit_central)
best <- joint_fits[[which.max(vapply(joint_fits,`[[`,numeric(1),"loglik"))]]

# The short joint run is deliberately bounded because each 30-dimensional
# numerical gradient is expensive. Finish with an exact one-dimensional
# conditional maximization of rho at the best joint nuisance parameters.
polish_objective <- function(rho) {
  p <- best$params; p$tenure_exact_anchor_retention_prob <- rho
  -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
polish <- optimize(polish_objective,c(1e-7,.95),tol=1e-6)
polished_params <- best$params
polished_params$tenure_exact_anchor_retention_prob <- polish$minimum
polished_eval <- e_step_eps(df_fit,polished_params,check_df=FALSE,
  suff_stats=FALSE)
polished <- best
polished$params <- polished_params
polished$loglik <- polished_eval$loglik
polished$gamma <- polished_eval$gamma
polished$job_change_posterior <- polished_eval$job_change_posterior
polished$rho_polished <- TRUE
joint_fits$rho_polished <- polished
if (polished$loglik>=best$loglik) best <- polished
saveRDS(polished,file.path(outdir,"fit_rho_polished_latest.rds"))

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior
  posterior_q <- sum(df_fit$weight*jp$expected_changes)/
    sum(df_fit$weight*jp$opportunities)
  p <- fit$params
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
conditional_fit <- structure(list(params=conditional_params,
  loglik=conditional_loglik,gamma=conditional_eval$gamma,
  job_change_posterior=conditional_eval$job_change_posterior,
  convergence=0L,iterations=NA_integer_),
  class="conditional_exact_anchor_fit")
comparison <- rbind(
  summarise_fit("Flexible baseline; no exact retention",base,29L),
  summarise_fit("Conditional exact retention only",conditional_fit,30L),
  summarise_fit("Joint exact-anchor retention",best,30L))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="No exact-anchor retention",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  p_value_boundary_mixture=.5*pchisq(LR,1,lower.tail=FALSE))
nesting <- data.frame(base_loglik=base$loglik,
  zero_retention_loglik=evaluate_rho(0),difference=
    evaluate_rho(0)-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("rho=0 does not reproduce the flexible baseline model")

cat("\nExact nesting check\n")
print(nesting,row.names=FALSE,digits=12)
cat("\nConditional profile\n")
print(profile,row.names=FALSE,digits=9)
cat("\nExact-anchor retention comparison\n")
print(comparison,row.names=FALSE,digits=9)
cat("\nLikelihood-ratio test\n")
print(lr,row.names=FALSE,digits=9)

write.csv(profile,profile_file,row.names=FALSE)
write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
saveRDS(list(base=base,conditional=conditional_fit,extension=best,
  fits=joint_fits,profile=profile,nesting=nesting,comparison=comparison,
  lr=lr),file.path(outdir,"fits_latest.rds"))
