# Estimate a one-parameter calendar-heaping extension of the corrected
# discrete-month tenure model on Panel A.  Conditional on a gross tenure
# report, h is the probability that the reported current-job start month is
# drawn from the marginal duration law restricted to January.  The existing
# independent-gross model is nested exactly at h=0.  No bootstrap is run.

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
df_full <- add_nominal_interview_months(df_qlfs)
df_full <- prepare_eps_estimation_data(df_full,allow_zero_tenure=TRUE)
interview_cols <- paste0("interview_month",1:3)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=interview_cols)

base_file <- "EM-tenure/output/results/job_change_monthly/fits_latest.rds"
if (!file.exists(base_file)) stop("Estimate the corrected Panel A model first")
base <- readRDS(base_file)$reset
base$params$tenure_measurement_model <- "monthly"
base$params$tenure_report_persistence <- 0
base$params$tenure_heaping_prob <- 0

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_heaping")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("HEAPING_TENURE_WORKERS","1"))
maxit <- as.integer(Sys.getenv("HEAPING_TENURE_MAXIT","100"))
resume <- identical(tolower(Sys.getenv("HEAPING_TENURE_RESUME","true")),
  "true")

nested_estep <- e_step_eps(df_fit,base$params,check_df=FALSE,
  suff_stats=FALSE)
nested_loglik <- nested_estep$loglik
# The saved base posterior is indexed by the earlier collapse that omitted
# interview month.  Recompute it on the calendar-aware cells before any
# weighted posterior summaries; the likelihood and parameters are unchanged.
base$gamma <- nested_estep$gamma
base$job_change_posterior <- nested_estep$job_change_posterior
nesting <- data.frame(base_loglik=base$loglik,
  heaping_zero_loglik=nested_loglik,difference=nested_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("h=0 does not exactly reproduce the corrected monthly likelihood")

fit_one <- function(label,h_start,start_params=base$params) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Estimating heaping start %s (h = %.2f)",label,h_start))
  fit <- fit_eps_piecewise_heaping_monthly(df_fit,start_params,
    heaping_start=h_start,q_start=base$params$job_change_prob,
    maxit=maxit,reltol=1e-9,pgtol=1e-7,workers=workers,verbose=1L)
  fit$objective_function <- NULL
  if (!is.finite(fit$loglik)) stop("Non-finite fit for start ",label)
  saveRDS(fit,path)
  fit
}

starts <- c(low=.05,main=.25,high=.60)
fit_low <- fit_one("low",starts[["low"]])
fit_main <- fit_one("main",starts[["main"]],fit_low$params)
fit_high <- fit_one("high",starts[["high"]],fit_main$params)
fits <- list(low=fit_low,main=fit_main,high=fit_high)
best_positive <- fits[[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]]
best <- if (base$loglik>=best_positive$loglik) base else best_positive
best$params$tenure_measurement_model <- "monthly"
best$params$tenure_report_persistence <- 0
best$params$tenure_heaping_prob <- if (base$loglik>=best_positive$loglik) 0 else
  best_positive$params$tenure_heaping_prob

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior
  posterior_q <- sum(df_fit$weight*jp$expected_changes)/
    sum(df_fit$weight*jp$opportunities)
  p <- fit$params
  heap <- if (is.null(p$tenure_heaping_prob)) 0 else p$tenure_heaping_prob
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    gross_january_heaping_prob=heap,
    prior_heaped_tenure_report_prob=p$eps*heap,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(summarise_fit("Independent gross reports",base,15L),
  summarise_fit("January-heaped gross reports",best,16L))
start_comparison <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,gross_january_heaping_prob=
      p$tenure_heaping_prob,job_change_prob=p$job_change_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="Positive January-heaping probability",
  LR=LR,df=1L,p_chibar2=if (LR==0) 1 else
    .5*pchisq(LR,1,lower.tail=FALSE))

observed_january <- numeric(3)
for (wave in 1:3) {
  keep <- df_full[[paste0("y",wave)]]==1L
  tenure_month <- round(12*df_full[[paste0("tenure",wave)]][keep])
  month <- df_full[[paste0("interview_month",wave)]][keep]
  observed_january[wave] <- weighted.mean(
    tenure_month%%12L==month-1L,df_full$weight[keep])
}
observed <- data.frame(quantity=c("Wave-specific minimum","Wave-specific mean",
  "Wave-specific maximum"),observed_january_share=c(min(observed_january),
  mean(observed_january),max(observed_january)))

cat("\nExact nesting check\n"); print(nesting,row.names=FALSE,digits=12)
cat("\nCalendar-heaping comparison\n")
print(comparison,row.names=FALSE,digits=8)
cat("\nStarting-value comparison\n")
print(start_comparison,row.names=FALSE,digits=8)
cat("\nBoundary likelihood-ratio test\n"); print(lr,row.names=FALSE,digits=8)
cat("\nObserved January start-date shares\n")
print(observed,row.names=FALSE,digits=8)

write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison,file.path(outdir,"starting_values_latest.csv"),
  row.names=FALSE)
write.csv(lr,file.path(outdir,"boundary_lr_latest.csv"),row.names=FALSE)
write.csv(observed,file.path(outdir,"observed_january_share_latest.csv"),
  row.names=FALSE)
saveRDS(list(base=base,heaping=best,best_positive=best_positive,fits=fits,
  nesting=nesting,comparison=comparison,start_comparison=start_comparison,
  lr=lr,observed=observed),file.path(outdir,"fits_latest.rds"))
