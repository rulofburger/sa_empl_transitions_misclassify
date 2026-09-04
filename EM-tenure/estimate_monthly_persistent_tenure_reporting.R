# Estimate the first persistent reported-start-date extension of the corrected
# discrete-month model on Panel A. Conditional on consecutive gross tenure
# reports within a latent job, rho is the probability that the preceding
# reported start-date anchor is retained rather than redrawn. The model nests
# the independent-gross specification exactly at rho=0. No bootstrap is run.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(df_qlfs,allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE)

base_file <- "EM-tenure/output/results/job_change_monthly/fits_latest.rds"
if (!file.exists(base_file)) stop("Estimate the corrected Panel A model first")
base <- readRDS(base_file)$reset
base$params$tenure_measurement_model <- "monthly"
base$params$tenure_report_persistence <- 0

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "persistent_reporting")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("PERSISTENT_TENURE_WORKERS","1"))
maxit <- as.integer(Sys.getenv("PERSISTENT_TENURE_MAXIT","80"))
resume <- identical(tolower(Sys.getenv("PERSISTENT_TENURE_RESUME","true")),
  "true")

nested_loglik <- e_step_eps(df_fit,base$params,check_df=FALSE,
  suff_stats=FALSE)$loglik
nesting <- data.frame(base_loglik=base$loglik,
  rho_zero_loglik=nested_loglik,difference=nested_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("rho=0 does not exactly reproduce the corrected monthly likelihood")

fit_one <- function(label,rho_start,start_params=base$params) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Estimating persistence start %s (rho = %.2f)",
    label,rho_start))
  fit <- fit_eps_piecewise_persistent_monthly(df_fit,start_params,
    persistence_start=rho_start,q_start=base$params$job_change_prob,
    maxit=maxit,reltol=1e-9,pgtol=1e-7,workers=workers,verbose=1L)
  fit$objective_function <- NULL
  if (!is.finite(fit$loglik)) stop("Non-finite fit for start ",label)
  saveRDS(fit,path)
  fit
}

starts <- c(low=.10,main=.50,high=.90)
fit_low <- fit_one("low",starts[["low"]])
fit_main <- fit_one("main",starts[["main"]])
# Reset rho to the deliberately high value, but retain the converged interior
# nuisance estimates so this checks the rho basin without repeating a long
# nuisance-parameter traverse from the nested model.
fit_high <- fit_one("high",starts[["high"]],fit_main$params)
fits <- list(low=fit_low,main=fit_main,high=fit_high)
best_positive <- fits[[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]]
# rho is boundary-nested. If all numerical positive-rho fits remain below the
# exact rho=0 likelihood, the constrained MLE is the nested fit itself.
best <- if (base$loglik >= best_positive$loglik) base else best_positive
best$params$tenure_measurement_model <- "monthly"
best$params$tenure_report_persistence <- if (base$loglik >=
  best_positive$loglik) 0 else best_positive$params$tenure_report_persistence

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior
  posterior_q <- sum(df_fit$weight*jp$expected_changes)/
    sum(df_fit$weight*jp$opportunities)
  p <- fit$params
  rho <- if (is.null(p$tenure_report_persistence)) 0 else
    p$tenure_report_persistence
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,tenure_report_persistence=rho,
    prior_retained_gross_pair_prob=p$eps^2*rho,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(summarise_fit("Independent gross reports",base,15L),
  summarise_fit("Persistent gross anchor",best,16L))
start_comparison <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,tenure_report_persistence=
      p$tenure_report_persistence,job_change_prob=p$job_change_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi)
}))
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="Positive tenure-report persistence",
  LR=LR,df=1L,p_chibar2=if (LR==0) 1 else .5*pchisq(LR,1,
    lower.tail=FALSE))

cat("\nExact nesting check\n"); print(nesting,row.names=FALSE,digits=12)
cat("\nPersistent tenure-reporting comparison\n")
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
saveRDS(list(base=base,persistent=best,best_positive=best_positive,fits=fits,
  nesting=nesting,
  comparison=comparison,start_comparison=start_comparison,lr=lr),
  file.path(outdir,"fits_latest.rds"))
