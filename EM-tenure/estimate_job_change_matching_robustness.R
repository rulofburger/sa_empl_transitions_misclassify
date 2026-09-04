# Re-estimate the preferred job-change/reset extension on the stricter B and C
# matched panels. Panel A remains the main specification. These fits change
# only the input panel and are written to a separate ignored output directory.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

main_fit_file <- "EM-tenure/output/results/job_change/fit_latest.rds"
if (!file.exists(main_fit_file)) stop(
  "Run estimate_job_change_tenure_contamination.R on Panel A first")
fit_a <- readRDS(main_fit_file)

outdir <- "EM-tenure/output/results/job_change/matching_robustness"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("JOB_CHANGE_WORKERS", "1"))
maxit <- as.integer(Sys.getenv("JOB_CHANGE_MAXIT", "60"))
resume <- identical(tolower(Sys.getenv("JOB_CHANGE_RESUME", "true")), "true")

load_panel <- function(letter) {
  old <- options(sa_empl_transitions.qlfs_3wave_panel =
    paste0("df_qlfs_", letter, ".rds"))
  on.exit(options(old), add=TRUE)
  ingest <- new.env(parent=globalenv())
  sys.source("scripts/ingest_data_3waves_SA.R", envir=ingest)
  prepared <- prepare_eps_estimation_data(ingest$df_qlfs)
  list(df=collapse_eps_cells(prepared), observations=nrow(prepared))
}

fit_panel <- function(letter) {
  fit_file <- file.path(outdir, paste0("panel_", letter, "_fit.rds"))
  if (resume && file.exists(fit_file)) {
    message("Loading saved Panel ", letter, " fit")
    return(readRDS(fit_file))
  }
  message("Preparing Panel ", letter)
  dat <- load_panel(letter)
  message(sprintf("Estimating Panel %s: %s observations, %s collapsed cells",
    letter, format(dat$observations, big.mark=","),
    format(nrow(dat$df), big.mark=",")))
  fit <- fit_eps_piecewise_job_change(dat$df, fit_a$params,
    q_start=fit_a$params$job_change_prob, maxit=maxit,
    reltol=1e-9, pgtol=1e-7, method="L-BFGS-B", verbose=1L,
    workers=workers)
  fit$objective_function <- NULL
  if (!is.finite(fit$loglik)) stop("Panel ", letter, " fit is not finite")
  ans <- list(panel=letter, fit=fit, df=dat$df,
    observations=dat$observations, collapsed_cells=nrow(dat$df))
  saveRDS(ans, fit_file)
  ans
}

results <- lapply(c("B", "C"), fit_panel)
names(results) <- c("B", "C")

# Recreate Panel A's prepared data only for comparable posterior and flow rates.
dat_a <- load_panel("A")
results <- c(list(A=list(panel="A", fit=fit_a, df=dat_a$df,
  observations=dat_a$observations, collapsed_cells=nrow(dat_a$df))), results)

comparison <- do.call(rbind, lapply(results, function(x) {
  fit <- x$fit
  p <- fit$params
  jp <- fit$job_change_posterior
  posterior_q <- sum(x$df$weight * jp$expected_changes) /
    sum(x$df$weight * jp$opportunities)
  rates <- duration_weighted_transition_rates(x$df, fit)[1, ]
  data.frame(panel=x$panel, observations=x$observations,
    collapsed_cells=x$collapsed_cells, loglik=fit$loglik,
    job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    status_misclassification=p$pi, tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,
    weighted_entry=rates$entry_rate, weighted_exit=rates$exit_rate,
    convergence=fit$convergence, iterations=fit$iterations)
}))

cat("\nJob-change model: matching-panel robustness\n")
print(comparison, row.names=FALSE, digits=8)
write.csv(comparison, file.path(outdir, "comparison_latest.csv"),
  row.names=FALSE)
saveRDS(list(results=lapply(results, function(x) x$fit),
  comparison=comparison), file.path(outdir, "fits_latest.rds"))
