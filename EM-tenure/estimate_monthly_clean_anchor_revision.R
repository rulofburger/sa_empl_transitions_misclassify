# Extend the preferred monthly calendar model so a current gross tenure report
# may revise a preceding clean reported start date. The prior gross-to-gross
# revision model is nested at eta=0.

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
  "calendar_revision/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the calendar-revision model first")
base <- readRDS(base_file)$revision
base$params$tenure_clean_anchor_revision_prob <- 0
base_eval <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_eval$loglik
base$gamma <- base_eval$gamma
base$job_change_posterior <- base_eval$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_clean_anchor_revision")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("CLEAN_ANCHOR_REVISION_WORKERS","1"))
maxit <- as.integer(Sys.getenv("CLEAN_ANCHOR_REVISION_MAXIT","80"))
resume <- identical(tolower(Sys.getenv("CLEAN_ANCHOR_REVISION_RESUME","true")),
  "true")

nested_loglik <- e_step_eps(df_fit,base$params,check_df=FALSE,
  suff_stats=FALSE)$loglik
nesting <- data.frame(base_loglik=base$loglik,
  clean_anchor_zero_loglik=nested_loglik,
  difference=nested_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("eta=0 does not reproduce the gross-anchor revision likelihood")

fit_one <- function(label,eta_start,start_params=base$params) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Estimating clean-anchor start %s (eta = %.2f)",
    label,eta_start))
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start_params,
    heaping_start=start_params$tenure_heaping_prob,
    revision_start=start_params$tenure_year_revision_prob,
    clean_anchor_revision_start=eta_start,
    q_start=start_params$job_change_prob,maxit=maxit,reltol=1e-9,
    pgtol=1e-7,workers=workers,verbose=1L)
  fit$objective_function <- NULL
  if (!is.finite(fit$loglik)) stop("Non-finite fit for start ",label)
  saveRDS(fit,path)
  fit
}

starts <- c(low=.03,main=.25,high=.65)
fit_main <- fit_one("main",starts[["main"]])
fit_high <- fit_one("high",starts[["high"]],fit_main$params)
fit_low <- fit_one("low",starts[["low"]],fit_main$params)
fits <- list(low=fit_low,main=fit_main,high=fit_high)
best_positive <- fits[[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]]
best <- if (base$loglik>=best_positive$loglik) base else best_positive

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
    clean_anchor_revision_prob=if (is.null(
      p$tenure_clean_anchor_revision_prob)) 0 else
      p$tenure_clean_anchor_revision_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Gross-anchor revisions only",base,17L),
  summarise_fit("Gross and clean-anchor revisions",best,18L))
start_comparison <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,clean_anchor_revision_prob=
      p$tenure_clean_anchor_revision_prob,gross_anchor_revision_prob=
      p$tenure_year_revision_prob,gross_january_heaping_prob=
      p$tenure_heaping_prob,job_change_prob=p$job_change_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="Positive clean-anchor revision probability",
  LR=LR,df=1L,p_chibar2=if (LR==0) 1 else
    .5*pchisq(LR,1,lower.tail=FALSE))

cat("\nExact nesting check\n"); print(nesting,row.names=FALSE,digits=12)
cat("\nClean-anchor revision comparison\n")
print(comparison,row.names=FALSE,digits=8)
cat("\nStarting-value comparison\n")
print(start_comparison,row.names=FALSE,digits=8)
cat("\nBoundary likelihood-ratio test\n"); print(lr,row.names=FALSE,digits=8)

write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison,file.path(outdir,"starting_values_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"boundary_lr_latest.csv"),row.names=FALSE)
saveRDS(list(base=base,extension=best,best_positive=best_positive,fits=fits,
  nesting=nesting,comparison=comparison,start_comparison=start_comparison,
  lr=lr),file.path(outdir,"fits_latest.rds"))
