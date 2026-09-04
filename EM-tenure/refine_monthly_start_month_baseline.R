# Continue the full-sample central-difference refinement of the flexible
# baseline start-month model without repeating the subsample mode search.

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
  "calendar_start_month_baseline")
fit_file <- file.path(outdir,"fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate the start-month model first")
saved <- readRDS(fit_file)
start <- saved$seasonal
workers <- as.integer(Sys.getenv("START_MONTH_REFINE_WORKERS","8"))
maxit <- as.integer(Sys.getenv("START_MONTH_REFINE_MAXIT","8"))

fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start$params,
  heaping_start=start$params$tenure_heaping_prob,
  revision_start=start$params$tenure_year_revision_prob,
  clean_anchor_revision_start=start$params$tenure_clean_anchor_revision_prob,
  start_month_probs_start=start$params$tenure_start_month_probs,
  q_start=start$params$job_change_prob,maxit=maxit,reltol=1e-9,pgtol=1e-7,
  workers=workers,verbose=1L,gradient_step=1e-3,
  gradient_scheme="central")
fit$objective_function <- NULL
saveRDS(fit,file.path(outdir,"fit_full_central_refined_latest.rds"))
if (fit$loglik>=saved$seasonal$loglik) {
  saved$fits$central_refined <- fit
  saved$seasonal <- fit
}

summarise_fit <- function(model,value,k) {
  rates <- duration_weighted_transition_rates(df_fit,value)[1,]
  jp <- value$job_change_posterior
  posterior_q <- sum(df_fit$weight*jp$expected_changes)/
    sum(df_fit$weight*jp$opportunities)
  p <- value$params
  data.frame(model=model,loglik=value$loglik,parameters=k,
    AIC=-2*value$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    convergence=value$convergence,iterations=value$iterations)
}
saved$comparison <- rbind(
  summarise_fit("Uniform baseline start months",saved$base,18L),
  summarise_fit("Flexible baseline start months",saved$seasonal,29L))
saved$month_distribution$estimated_baseline <-
  saved$seasonal$params$tenure_start_month_probs
LR <- max(0,2*(saved$seasonal$loglik-saved$base$loglik))
saved$lr <- data.frame(
  comparison="Equal baseline start-month probabilities",LR=LR,df=11L,
  p_value=pchisq(LR,11,lower.tail=FALSE))
saveRDS(saved,fit_file)
write.csv(saved$comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(saved$month_distribution,file.path(outdir,
  "start_month_distribution_latest.csv"),row.names=FALSE)
write.csv(saved$lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)

cat("\nExtended full-sample central refinement\n")
print(saved$comparison,row.names=FALSE,digits=9)
cat("\nBaseline start-month distribution\n")
print(saved$month_distribution,row.names=FALSE,digits=9)
cat("\nLikelihood-ratio test\n")
print(saved$lr,row.names=FALSE,digits=9)
