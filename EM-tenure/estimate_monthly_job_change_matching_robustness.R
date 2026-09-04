# Re-estimate the corrected discrete-month job-change/reset model on the
# stricter B and C matched panels. Panel A remains the main specification.
# Zero-month tenure is retained and tenure contributions are monthly masses.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

main_fit_file <- "EM-tenure/output/results/job_change_monthly/fits_latest.rds"
if (!file.exists(main_fit_file)) stop(
  "Run estimate_monthly_job_change_tenure_contamination.R first")
fit_a_all <- readRDS(main_fit_file)
fit_a <- fit_a_all$reset

outdir <- "EM-tenure/output/results/job_change_monthly/matching_robustness"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("MONTHLY_MATCHING_WORKERS", "1"))
maxit <- as.integer(Sys.getenv("MONTHLY_MATCHING_MAXIT", "80"))
resume <- identical(tolower(Sys.getenv("MONTHLY_MATCHING_RESUME", "true")),
  "true")

load_panel <- function(letter) {
  old <- options(
    sa_empl_transitions.qlfs_3wave_panel=paste0("df_qlfs_", letter, ".rds"),
    sa_empl_transitions.preserve_zero_tenure=TRUE)
  on.exit(options(old), add=TRUE)
  ingest <- new.env(parent=globalenv())
  sys.source("scripts/ingest_data_3waves_SA.R", envir=ingest)
  prepared <- prepare_eps_estimation_data(ingest$df_qlfs,
    allow_zero_tenure=TRUE)
  grid <- unlist(prepared[paste0("tenure", 1:3)])
  if (any(abs(12*grid - round(12*grid)) > 1e-6, na.rm=TRUE))
    stop("Panel ", letter, " contains tenure observations off the monthly grid")
  list(df=collapse_eps_cells(prepared, allow_zero_tenure=TRUE),
    observations=nrow(prepared))
}

fit_one_start <- function(letter, dat, label, q_start) {
  path <- file.path(outdir,
    paste0("panel_", letter, "_fit_", label, "_latest.rds"))
  if (resume && file.exists(path)) {
    message("Loading Panel ", letter, " start ", label)
    return(readRDS(path))
  }
  message(sprintf("Estimating Panel %s, start %s (q = %.4f)",
    letter, label, q_start))
  fit <- fit_eps_piecewise_job_change(dat$df, fit_a$params,
    q_start=q_start, maxit=maxit, reltol=1e-9, pgtol=1e-7,
    method="L-BFGS-B", verbose=1L, workers=workers,
    tenure_measurement_model="monthly")
  fit$objective_function <- NULL
  if (!is.finite(fit$loglik))
    stop("Panel ", letter, ", start ", label, " has non-finite likelihood")
  saveRDS(fit, path)
  fit
}

fit_panel <- function(letter) {
  message("Preparing Panel ", letter)
  dat <- load_panel(letter)
  message(sprintf("Panel %s: %s observations, %s collapsed cells",
    letter, format(dat$observations, big.mark=","),
    format(nrow(dat$df), big.mark=",")))
  q_a <- fit_a$params$job_change_prob
  starts <- c(low=max(.002, q_a/2), main=q_a, high=min(.05, 3*q_a))
  fits <- lapply(names(starts), function(label)
    fit_one_start(letter, dat, label, starts[[label]]))
  names(fits) <- names(starts)
  best <- fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
  list(panel=letter, fit=best, fits=fits, df=dat$df,
    observations=dat$observations, collapsed_cells=nrow(dat$df))
}

results <- lapply(c("B", "C"), fit_panel)
names(results) <- c("B", "C")
dat_a <- load_panel("A")
results <- c(list(A=list(panel="A", fit=fit_a,
  fits=fit_a_all$reset_fits, df=dat_a$df,
  observations=dat_a$observations, collapsed_cells=nrow(dat_a$df))), results)

summarise_panel <- function(x) {
  fit <- x$fit
  p <- fit$params
  jp <- fit$job_change_posterior
  posterior_q <- sum(x$df$weight * jp$expected_changes) /
    sum(x$df$weight * jp$opportunities)
  rates <- duration_weighted_transition_rates(x$df, fit)[1, ]
  data.frame(panel=x$panel, observations=x$observations,
    collapsed_cells=x$collapsed_cells, loglik=fit$loglik,
    parameters=length(fit$par_unconstrained),
    AIC=-2*fit$loglik + 2*length(fit$par_unconstrained),
    entry_rate=rates$entry_rate, exit_rate=rates$exit_rate,
    initial_employment=p$alpha, status_misclassification=p$pi,
    tenure_contamination=p$eps, timegap_contamination=p$eps_d,
    job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    convergence=fit$convergence, iterations=fit$iterations)
}
comparison <- do.call(rbind, lapply(results, summarise_panel))

start_comparison <- do.call(rbind, lapply(results, function(x) {
  do.call(rbind, lapply(names(x$fits), function(label) {
    fit <- x$fits[[label]]
    data.frame(panel=x$panel, start=label, loglik=fit$loglik,
      convergence=fit$convergence, iterations=fit$iterations,
      job_change_prob=fit$params$job_change_prob,
      status_misclassification=fit$params$pi,
      tenure_contamination=fit$params$eps,
      timegap_contamination=fit$params$eps_d)
  }))
}))

cat("\nCorrected discrete-month model: matching-panel robustness\n")
print(comparison, row.names=FALSE, digits=8)
cat("\nStarting-value comparison\n")
print(start_comparison, row.names=FALSE, digits=8)
write.csv(comparison, file.path(outdir, "comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison, file.path(outdir, "starting_values_latest.csv"),
  row.names=FALSE)
saveRDS(list(results=lapply(results, function(x) x$fit),
  all_starts=lapply(results, function(x) x$fits),
  comparison=comparison, start_comparison=start_comparison),
  file.path(outdir, "fits_latest.rds"))
