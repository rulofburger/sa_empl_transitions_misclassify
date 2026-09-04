# Estimate the whole-year job-start-date revision extension on Panel A.
# January heaping remains in the gross-report distribution; omega additionally
# links consecutive gross reports by preserving start month but changing year.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
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
  "calendar_heaping/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the January-heaping model first")
base <- readRDS(base_file)$heaping
base$params$tenure_measurement_model <- "monthly"
base$params$tenure_report_persistence <- 0
base$params$tenure_year_revision_prob <- 0
base_estep <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_estep$loglik
base$gamma <- base_estep$gamma
base$job_change_posterior <- base_estep$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_revision")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("REVISION_TENURE_WORKERS","1"))
maxit <- as.integer(Sys.getenv("REVISION_TENURE_MAXIT","100"))
resume <- identical(tolower(Sys.getenv("REVISION_TENURE_RESUME","true")),
  "true")

nested_params <- base$params
nested_params$tenure_year_revision_prob <- 0
nested_loglik <- e_step_eps(df_fit,nested_params,check_df=FALSE,
  suff_stats=FALSE)$loglik
nesting <- data.frame(base_loglik=base$loglik,
  revision_zero_loglik=nested_loglik,
  difference=nested_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("omega=0 does not exactly reproduce the January-heaping likelihood")

fit_one <- function(label,omega_start,start_params=base$params) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Estimating revision start %s (omega = %.2f)",
    label,omega_start))
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start_params,
    heaping_start=base$params$tenure_heaping_prob,
    revision_start=omega_start,q_start=start_params$job_change_prob,
    maxit=maxit,reltol=1e-9,pgtol=1e-7,workers=workers,verbose=1L)
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
    whole_year_revision_prob=p$tenure_year_revision_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(summarise_fit("January-heaped gross reports",base,16L),
  summarise_fit("January heaping + whole-year revisions",best,17L))
start_comparison <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,whole_year_revision_prob=
      p$tenure_year_revision_prob,gross_january_heaping_prob=
      p$tenure_heaping_prob,job_change_prob=p$job_change_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="Positive whole-year revision probability",
  LR=LR,df=1L,p_chibar2=if (LR==0) 1 else
    .5*pchisq(LR,1,lower.tail=FALSE))

observed_pair <- function(first,second) {
  keep <- df_full[[paste0("y",first)]]==1L &
    df_full[[paste0("y",second)]]==1L
  prior <- round(12*df_full[[paste0("tenure",first)]][keep])
  current <- round(12*df_full[[paste0("tenure",second)]][keep])
  gap <- second-first
  revised <- (current-prior-gap*.TENURE_MONTHS_PER_WAVE)%%12L==0L &
    current!=prior+gap*.TENURE_MONTHS_PER_WAVE
  weighted.mean(revised,df_full$weight[keep])
}
observed <- data.frame(pair=c("Waves 1-2","Waves 2-3","Waves 1-3"),
  observed_whole_year_revision_share=c(observed_pair(1,2),
    observed_pair(2,3),observed_pair(1,3)))

cat("\nExact nesting check\n"); print(nesting,row.names=FALSE,digits=12)
cat("\nCalendar-revision comparison\n")
print(comparison,row.names=FALSE,digits=8)
cat("\nStarting-value comparison\n")
print(start_comparison,row.names=FALSE,digits=8)
cat("\nBoundary likelihood-ratio test\n"); print(lr,row.names=FALSE,digits=8)
cat("\nObserved same-month/different-year revision shares\n")
print(observed,row.names=FALSE,digits=8)

write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison,file.path(outdir,"starting_values_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"boundary_lr_latest.csv"),row.names=FALSE)
write.csv(observed,file.path(outdir,"observed_revision_share_latest.csv"),
  row.names=FALSE)
saveRDS(list(base=base,revision=best,best_positive=best_positive,fits=fits,
  nesting=nesting,comparison=comparison,start_comparison=start_comparison,
  lr=lr,observed=observed),file.path(outdir,"fits_latest.rds"))
