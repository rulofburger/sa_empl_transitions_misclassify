# Reoptimize one fixed reset-probability profile point for the corrected
# monthly model. The target is supplied through MONTHLY_PROFILE_Q; the nearest
# converged saved profile point is used as the start.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs,
  allow_zero_tenure=TRUE), allow_zero_tenure=TRUE)

q <- as.numeric(Sys.getenv("MONTHLY_PROFILE_Q", ""))
if (length(q)!=1L || !is.finite(q) || q<=0 || q>=.5)
  stop("Set MONTHLY_PROFILE_Q to an interior probability")
workers <- as.integer(Sys.getenv("MONTHLY_PROFILE_WORKERS", "1"))
maxit <- as.integer(Sys.getenv("MONTHLY_PROFILE_MAXIT", "35"))
outdir <- "EM-tenure/output/results/job_change_monthly/inference"

main <- readRDS("EM-tenure/output/results/job_change_monthly/fits_latest.rds")$reset
paths <- Sys.glob(file.path(outdir, "profile_q_*.rds"))
candidates <- list(main=main)
if (length(paths)) for (path in paths) {
  value <- readRDS(path)
  if (identical(value$convergence, 0L)) candidates[[basename(path)]] <- value
}
distance <- vapply(candidates, function(x)
  abs(x$params$job_change_prob-q), numeric(1))
start_name <- names(candidates)[which.min(distance)]
message(sprintf("Profiling q=%.8f from %s", q, start_name))
fit <- fit_eps_piecewise_job_change(df_fit, candidates[[start_name]]$params,
  q_fixed=q, maxit=maxit, reltol=1e-10, pgtol=1e-7, workers=workers,
  tenure_measurement_model="monthly")
fit$objective_function <- NULL
path <- file.path(outdir,
  paste0("profile_q_", formatC(q,digits=8,format="f"), ".rds"))
saveRDS(fit, path)
cat("\nRefined monthly profile point\n")
print(data.frame(q=q, loglik=fit$loglik, convergence=fit$convergence,
  iterations=fit$iterations, start=start_name), row.names=FALSE, digits=9)
if (fit$convergence!=0L) stop("Refined profile point did not converge")
