# Add a short-range reported-start-date revision channel to the exact-anchor
# model. Given a gross current tenure report and a usable preceding anchor, the
# reporting sequence is: retain the anchor exactly with probability rho;
# otherwise make a non-zero local revision with probability kappa; otherwise
# use the existing whole-year-revision/redraw mixture. The fixed local kernel
# has support +/- 1,...,+/- 6 months and discrete-Laplace scale three months.
# The exact-anchor model is nested at kappa=0.

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
  "calendar_exact_anchor_retention/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the exact-anchor model first")
base <- readRDS(base_file)$extension
base$params$tenure_local_revision_prob <- 0
base_eval <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_eval$loglik
base$gamma <- base_eval$gamma
base$job_change_posterior <- base_eval$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_local_anchor_revision")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("LOCAL_ANCHOR_WORKERS","8"))
maxit <- as.integer(Sys.getenv("LOCAL_ANCHOR_MAXIT","6"))
screen_maxit <- as.integer(Sys.getenv("LOCAL_ANCHOR_SCREEN_MAXIT","2"))
resume <- identical(tolower(Sys.getenv("LOCAL_ANCHOR_RESUME","true")),
  "true")

evaluate_kappa <- function(kappa,reference=base) {
  p <- reference$params
  p$tenure_local_revision_prob <- kappa
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}

# A sparse conditional profile is sufficient to expose a boundary solution or
# a broad/weakly identified interior maximum before expensive joint fitting.
profile_grid <- c(0,.01,.025,.05,.10,.20,.35,.50,.70,.90)
profile_file <- file.path(outdir,"conditional_profile_latest.csv")
if (resume && file.exists(profile_file)) {
  profile <- read.csv(profile_file)
} else {
  profile <- data.frame(kappa=profile_grid,loglik=NA_real_)
  for (j in seq_along(profile_grid)) {
    message("Conditional local-revision profile: kappa=",profile_grid[j])
    profile$loglik[j] <- evaluate_kappa(profile_grid[j])
    write.csv(profile,profile_file,row.names=FALSE)
  }
}
conditional <- optimize(function(kappa) -evaluate_kappa(kappa),
  interval=c(1e-7,.98),tol=1e-6)
conditional_params <- base$params
conditional_params$tenure_local_revision_prob <- conditional$minimum
conditional_eval <- e_step_eps(df_fit,conditional_params,check_df=FALSE,
  suff_stats=FALSE)
conditional_fit <- structure(list(params=conditional_params,
  loglik=-conditional$objective,gamma=conditional_eval$gamma,
  job_change_posterior=conditional_eval$job_change_posterior,
  convergence=0L,iterations=NA_integer_),
  class="conditional_local_anchor_fit")

make_start <- function(kappa,label) {
  value <- conditional_fit
  value$params$tenure_local_revision_prob <- kappa
  evaluation <- e_step_eps(df_fit,value$params,check_df=FALSE,
    suff_stats=FALSE)
  value$loglik <- evaluation$loglik
  value$gamma <- evaluation$gamma
  value$job_change_posterior <- evaluation$job_change_posterior
  value$stage <- label
  value
}
reporting_names <- c("pi","eps","eps_d","job_change","tenure_heaping",
  "tenure_year_revision","tenure_clean_anchor_revision",
  "tenure_exact_anchor_retention","tenure_local_revision")
transition_names <- c("alpha",paste0("log_hg",1:5),paste0("log_hd",1:5))
fit_stage <- function(label,start,scheme,iterations,
    free_names=reporting_names) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Joint local-revision refinement: ",label)
  p <- start$params
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,p,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    q_start=p$job_change_prob,maxit=iterations,reltol=1e-9,pgtol=1e-7,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=5e-4,
    free_names=free_names,gradient_scheme=scheme)
  fit$objective_function <- NULL
  fit$stage <- label
  saveRDS(fit,path)
  fit
}

# Screen well-separated values of kappa. The conditional optimum gets the
# longer joint run; the tails get short runs that can reveal another basin.
main <- fit_stage("conditional_forward",conditional_fit,"forward",maxit)
low <- fit_stage("low_forward",make_start(.02,"low"),"forward",screen_maxit)
high <- fit_stage("high_forward",make_start(.65,"high"),"forward",screen_maxit)
screened <- list(main=main,low=low,high=high)
screen_best_name <- names(screened)[which.max(vapply(screened,`[[`,numeric(1),
  "loglik"))]
screen_best <- screened[[screen_best_name]]
central <- fit_stage(paste0(screen_best_name,"_central"),screen_best,
  "central",min(3L,maxit))
candidates <- c(screened,list(central=central))
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

# One transition-hazard block verifies that the new reporting channel is not
# borrowing its fit merely because the inherited transition rates were held
# fixed. It is intentionally followed by the exact scalar kappa polish below.
transition <- fit_stage("transition_check",best,"forward",1L,
  free_names=transition_names)
candidates$transition <- transition
best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]

# Polish the new mixing probability conditionally and synchronize both stored
# parameter representations so resumed inference cannot use a stale value.
polish_objective <- function(kappa) {
  p <- best$params; p$tenure_local_revision_prob <- kappa
  -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
polish <- optimize(polish_objective,c(1e-7,.98),tol=1e-7)
polished_params <- best$params
polished_params$tenure_local_revision_prob <- polish$minimum
polished_eval <- e_step_eps(df_fit,polished_params,check_df=FALSE,
  suff_stats=FALSE)
polished <- best
polished$params <- polished_params
polished$par_unconstrained["tenure_local_revision"] <-
  qlogis(polish$minimum)
polished$loglik <- polished_eval$loglik
polished$gamma <- polished_eval$gamma
polished$job_change_posterior <- polished_eval$job_change_posterior
polished$stage <- "kappa_polished"
saveRDS(polished,file.path(outdir,"fit_kappa_polished_latest.rds"))
if (polished$loglik>=best$loglik) best <- polished

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
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    local_anchor_revision_prob=p$tenure_local_revision_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Exact-anchor model; no local revision",base,30L),
  summarise_fit("Conditional local revision only",conditional_fit,31L),
  summarise_fit("Joint local anchor revision",best,31L))
refinement <- do.call(rbind,lapply(names(candidates),function(label) {
  f <- candidates[[label]]
  data.frame(stage=label,loglik=f$loglik,
    gain_over_exact=f$loglik-base$loglik,convergence=f$convergence,
    iterations=f$iterations,
    local_anchor_revision=f$params$tenure_local_revision_prob,
    exact_anchor_retention=f$params$tenure_exact_anchor_retention_prob,
    tenure_contamination=f$params$eps,
    gross_anchor_revision=f$params$tenure_year_revision_prob,
    clean_anchor_revision=f$params$tenure_clean_anchor_revision_prob)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="No local anchor revision",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  p_value_boundary_mixture=.5*pchisq(LR,1,lower.tail=FALSE))
zero_loglik <- evaluate_kappa(0)
nesting <- data.frame(exact_model_loglik=base$loglik,
  zero_local_revision_loglik=zero_loglik,
  difference=zero_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("kappa=0 does not reproduce the exact-anchor model")

cat("\nLocal-revision nesting check\n")
print(nesting,row.names=FALSE,digits=12)
cat("\nConditional profile\n")
print(profile,row.names=FALSE,digits=10)
cat("\nMultistart refinement\n")
print(refinement,row.names=FALSE,digits=10)
cat("\nLocal-anchor revision comparison\n")
print(comparison,row.names=FALSE,digits=10)
cat("\nLikelihood-ratio test\n")
print(lr,row.names=FALSE,digits=10)

write.csv(profile,profile_file,row.names=FALSE)
write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(refinement,file.path(outdir,"refinement_summary_latest.csv"),
  row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
saveRDS(list(base=base,conditional=conditional_fit,extension=best,
  fits=candidates,profile=profile,nesting=nesting,refinement=refinement,
  comparison=comparison,lr=lr),file.path(outdir,"fits_latest.rds"))
