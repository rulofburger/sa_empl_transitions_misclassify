# Estimate the first start-date improvement: a one-parameter job-to-job reset
# process inside latent employment spells. The established Table 7 estimator
# and outputs are not overwritten. No bootstrap is run.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
baseline_file <-
  "EM-tenure/output/results/timegap_contamination_robustness/fits_latest.rds"
if (!file.exists(baseline_file)) stop("Run the Table 7 robustness estimator first")
saved <- readRDS(baseline_file)
baseline <- saved$joint$fit
if (nrow(baseline$gamma) != nrow(df_fit)) stop("Saved fit does not match data")

outdir <- "EM-tenure/output/results/job_change"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Exact nesting is a non-negotiable acceptance check.
p0 <- baseline$params
p0$job_change_prob <- 0
ll0 <- e_step_eps(df_fit, p0, check_df=FALSE, suff_stats=FALSE)$loglik
if (!identical(ll0, baseline$loglik)) stop(
  sprintf("q=0 failed exact nesting: %.16g versus %.16g", ll0, baseline$loglik))

evaluate_q <- function(q, posterior=FALSE) {
  p <- baseline$params
  p$job_change_prob <- q
  e <- e_step_eps(df_fit, p, check_df=FALSE, suff_stats=FALSE)
  rate <- NA_real_
  if (posterior) {
    jp <- e$job_change_posterior
    rate <- sum(df_fit$weight * jp$expected_changes) /
      sum(df_fit$weight * jp$opportunities)
  }
  c(loglik=e$loglik, posterior_rate=rate)
}

profile_file <- file.path(outdir, "conditional_profile_latest.csv")
resume_profile <- identical(tolower(
  Sys.getenv("JOB_CHANGE_PROFILE_RESUME", "true")), "true")
if (resume_profile && file.exists(profile_file)) {
  message("Loading saved conditional q profile")
  profile <- read.csv(profile_file)
  q_conditional <- profile$job_change_prob[which.max(profile$loglik)]
} else {
  profile_grid <- c(0, .001, .0025, .005, .01, .015, .02, .03, .05, .08, .12)
  message("Evaluating conditional q profile")
  profile <- do.call(rbind, lapply(profile_grid, function(q) {
    value <- evaluate_q(q, posterior=TRUE)
    data.frame(job_change_prob=q, loglik=value["loglik"],
      loglik_difference=value["loglik"]-baseline$loglik,
      posterior_rate=value["posterior_rate"])
  }))
  conditional <- optimize(function(q) -evaluate_q(q)["loglik"],
    interval=c(1e-7, .20), tol=2e-5)
  q_conditional <- conditional$minimum
  ll_conditional <- -conditional$objective
  profile <- rbind(profile, data.frame(job_change_prob=q_conditional,
    loglik=ll_conditional, loglik_difference=ll_conditional-baseline$loglik,
    posterior_rate=evaluate_q(q_conditional, posterior=TRUE)["posterior_rate"]))
  profile <- profile[order(profile$job_change_prob), ]
  write.csv(profile, profile_file, row.names=FALSE)
}

fit_file <- file.path(outdir, "fit_latest.rds")
resume <- identical(tolower(Sys.getenv("JOB_CHANGE_RESUME", "true")), "true")
if (resume && file.exists(fit_file)) {
  message("Loading saved job-change fit")
  fit <- readRDS(fit_file)
} else {
  maxit <- as.integer(Sys.getenv("JOB_CHANGE_MAXIT", "100"))
  workers <- as.integer(Sys.getenv("JOB_CHANGE_WORKERS", "1"))
  message(sprintf("Jointly optimizing from conditional q=%.6f", q_conditional))
  fit <- fit_eps_piecewise_job_change(df_fit, baseline$params,
    q_start=q_conditional, maxit=maxit, reltol=1e-9, pgtol=1e-7,
    method="L-BFGS-B", verbose=1L, workers=workers)
  fit$objective_function <- NULL
  saveRDS(fit, fit_file)
}

# Verify that the accepted fit is finite, nested, and internally coherent.
if (!is.finite(fit$loglik) || fit$loglik < baseline$loglik - 1e-6)
  stop("Job-change fit is invalid or worse than its nested model")
jp <- fit$job_change_posterior
posterior_rate <- sum(df_fit$weight * jp$expected_changes) /
  sum(df_fit$weight * jp$opportunities)
rates0 <- duration_weighted_transition_rates(df_fit, baseline)[1, ]
rates1 <- duration_weighted_transition_rates(df_fit, fit)[1, ]

comparison <- data.frame(
  model=c("Preferred Table 7", "Plus E-E job change"),
  parameters=c(14L, 15L), loglik=c(baseline$loglik, fit$loglik),
  AIC=c(-2*baseline$loglik+28, -2*fit$loglik+30),
  job_change_prob=c(0, fit$params$job_change_prob),
  posterior_job_change_rate=c(0, posterior_rate),
  alpha=c(baseline$params$alpha, fit$params$alpha),
  status_misclassification=c(baseline$params$pi, fit$params$pi),
  tenure_contamination=c(baseline$params$eps, fit$params$eps),
  timegap_contamination=c(baseline$params$eps_d, fit$params$eps_d),
  weighted_exit=c(rates0$exit_rate, rates1$exit_rate),
  weighted_entry=c(rates0$entry_rate, rates1$entry_rate),
  convergence=c(0L, fit$convergence), iterations=c(NA_integer_, fit$iterations))
lr <- data.frame(comparison="Job-change probability greater than zero",
  LR=max(0, 2*(fit$loglik-baseline$loglik)), df=1L,
  p_chibar2=if (fit$loglik <= baseline$loglik) 1 else
    .5*pchisq(2*(fit$loglik-baseline$loglik), 1, lower.tail=FALSE))

cat("\nConditional job-change profile\n")
print(profile, row.names=FALSE, digits=8)
cat("\nJoint model comparison\n")
print(comparison, row.names=FALSE, digits=8)
cat("\nBoundary likelihood-ratio test\n")
print(lr, row.names=FALSE, digits=8)

write.csv(comparison, file.path(outdir, "model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(lr, file.path(outdir, "boundary_lr_latest.csv"), row.names=FALSE)
saveRDS(list(fit=fit, baseline=baseline, comparison=comparison,
  profile=profile, lr=lr), file.path(outdir, "fits_latest.rds"))
