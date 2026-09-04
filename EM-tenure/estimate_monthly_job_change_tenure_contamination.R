# Estimate the job-start-month version of the three-wave tenure-contamination
# model on Panel A. Zero-month reports are retained and tenure contributions
# are probability masses on the monthly grid. Existing Table 7 results are not
# overwritten. No bootstrap is run.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(df_qlfs, allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full, allow_zero_tenure=TRUE)

g <- unlist(df_full[paste0("tenure", 1:3)])
if (any(abs(12*g - round(12*g)) > 1e-6, na.rm=TRUE))
  stop("Panel A tenure contains observations off the monthly grid")

old <- readRDS("EM-tenure/output/results/job_change/fit_latest.rds")
start <- old$params
start$tenure_measurement_model <- "monthly"

outdir <- "EM-tenure/output/results/job_change_monthly"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("MONTHLY_ESTIMATION_WORKERS", "1"))
maxit <- as.integer(Sys.getenv("MONTHLY_ESTIMATION_MAXIT", "80"))
resume <- identical(tolower(Sys.getenv("MONTHLY_ESTIMATION_RESUME", "true")),
  "true")

fit_cached <- function(path, expression) {
  if (resume && file.exists(path)) return(readRDS(path))
  value <- force(expression)
  value$objective_function <- NULL
  saveRDS(value, path)
  value
}

message("Estimating discrete-month model without within-employment resets")
no_reset <- fit_cached(file.path(outdir, "fit_no_reset_latest.rds"),
  fit_eps_piecewise_job_change(df_fit, start, q_fixed=0, maxit=maxit,
    reltol=1e-9, pgtol=1e-7, workers=workers,
    tenure_measurement_model="monthly"))

message("Estimating discrete-month model with within-employment resets")
starts <- c(low=.008, main=max(.02, start$job_change_prob), high=.06)
reset_fits <- lapply(names(starts), function(label) fit_cached(
  file.path(outdir, paste0("fit_reset_", label, "_latest.rds")),
  fit_eps_piecewise_job_change(df_fit, no_reset$params,
    q_start=starts[[label]], maxit=maxit, reltol=1e-9, pgtol=1e-7,
    workers=workers, tenure_measurement_model="monthly")))
names(reset_fits) <- names(starts)
reset <- reset_fits[[which.max(vapply(reset_fits, `[[`, numeric(1),
  "loglik"))]]

summarise_fit <- function(label, fit) {
  rates <- duration_weighted_transition_rates(df_fit, fit)[1, ]
  jp <- fit$job_change_posterior
  posterior_q <- if (is.null(jp)) 0 else
    sum(df_fit$weight * jp$expected_changes) /
      sum(df_fit$weight * jp$opportunities)
  data.frame(model=label, loglik=fit$loglik,
    parameters=length(fit$par_unconstrained), AIC=-2*fit$loglik +
      2*length(fit$par_unconstrained), convergence=fit$convergence,
    iterations=fit$iterations, entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate, initial_employment=fit$params$alpha,
    status_misclassification=fit$params$pi,
    tenure_contamination=fit$params$eps,
    timegap_contamination=fit$params$eps_d,
    job_change_prob=fit$params$job_change_prob,
    posterior_job_change_rate=posterior_q)
}
comparison <- rbind(summarise_fit("No job reset", no_reset),
  summarise_fit("Job reset", reset))
start_comparison <- do.call(rbind, lapply(names(reset_fits), function(label) {
  fit <- reset_fits[[label]]
  data.frame(start=label, loglik=fit$loglik,
    convergence=fit$convergence, iterations=fit$iterations,
    job_change_prob=fit$params$job_change_prob,
    status_misclassification=fit$params$pi,
    tenure_contamination=fit$params$eps)
}))
lr <- data.frame(comparison="Positive within-employment reset probability",
  LR=max(0, 2*(reset$loglik-no_reset$loglik)), df=1L,
  p_chibar2=if (reset$loglik <= no_reset$loglik) 1 else
    .5*pchisq(2*(reset$loglik-no_reset$loglik), 1, lower.tail=FALSE))

cat("\nPanel A discrete-month model comparison\n")
print(comparison, row.names=FALSE, digits=8)
cat("\nReset-model starting-value comparison\n")
print(start_comparison, row.names=FALSE, digits=8)
cat("\nBoundary likelihood-ratio test\n")
print(lr, row.names=FALSE, digits=8)

write.csv(comparison, file.path(outdir, "model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison,
  file.path(outdir, "starting_values_latest.csv"), row.names=FALSE)
write.csv(lr, file.path(outdir, "boundary_lr_latest.csv"), row.names=FALSE)
saveRDS(list(no_reset=no_reset, reset=reset, reset_fits=reset_fits,
  comparison=comparison, start_comparison=start_comparison, lr=lr),
  file.path(outdir, "fits_latest.rds"))
