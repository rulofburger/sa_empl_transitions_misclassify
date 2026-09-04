# Repeated recovery validation for the discrete-month job-change model.
# Runs a resumable sample-size grid, checks a dispersed start at the first
# replication of each size, and reports bias and RMSE. No empirical fit is
# overwritten and no bootstrap is run.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")

old <- readRDS("EM-tenure/output/results/job_change/fit_latest.rds")
truth <- old$params
truth$tenure_measurement_model <- "monthly"

sizes <- as.integer(strsplit(Sys.getenv("MONTHLY_RECOVERY_SIZES",
  "10000,40000"), ",", fixed=TRUE)[[1]])
reps <- as.integer(Sys.getenv("MONTHLY_RECOVERY_REPS", "3"))
workers <- as.integer(Sys.getenv("MONTHLY_RECOVERY_WORKERS", "1"))
maxit <- as.integer(Sys.getenv("MONTHLY_RECOVERY_MAXIT", "45"))
resume <- identical(tolower(Sys.getenv("MONTHLY_RECOVERY_RESUME", "true")),
  "true")
if (any(!is.finite(sizes) | sizes < 1000) || reps < 1L)
  stop("Invalid recovery grid")

outdir <- "EM-tenure/output/results/job_change_monthly/recovery_grid"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

fit_cached <- function(path, expression) {
  if (resume && file.exists(path)) return(readRDS(path))
  value <- force(expression)
  value$objective_function <- NULL
  saveRDS(value, path)
  value
}

evaluate_score <- function(df, name, step) {
  lo <- hi <- truth
  lo[[name]] <- truth[[name]] - step
  hi[[name]] <- truth[[name]] + step
  (e_step_eps(df, hi, suff_stats=FALSE)$loglik -
    e_step_eps(df, lo, suff_stats=FALSE)$loglik) / (2*step*nrow(df))
}

rows <- list()
row_id <- 0L
for (n in sizes) for (r in seq_len(reps)) {
  seed <- 27000L + n + r
  message(sprintf("Recovery n=%d, replication=%d/%d", n, r, reps))
  df <- simulate_eps_piecewise_job_change(n, truth, seed=seed)

  start <- truth
  start$job_change_prob <- .015
  start$pi <- truth$pi * 1.15
  start$eps <- truth$eps * .90
  main_path <- file.path(outdir,
    sprintf("fit_n%d_rep%d_main.rds", n, r))
  main <- fit_cached(main_path, fit_eps_piecewise_job_change(df, start,
    q_start=start$job_change_prob, maxit=maxit, reltol=1e-9, pgtol=1e-7,
    workers=workers, tenure_measurement_model="monthly"))

  candidates <- list(main=main)
  if (r == 1L) {
    dispersed <- truth
    dispersed$job_change_prob <- .06
    dispersed$pi <- min(.10, truth$pi*1.6)
    dispersed$eps <- min(.60, truth$eps*1.35)
    high_path <- file.path(outdir,
      sprintf("fit_n%d_rep%d_high.rds", n, r))
    candidates$high <- fit_cached(high_path,
      fit_eps_piecewise_job_change(df, dispersed,
        q_start=dispersed$job_change_prob, maxit=maxit, reltol=1e-9,
        pgtol=1e-7, workers=workers,
        tenure_measurement_model="monthly"))
  }
  best_name <- names(candidates)[which.max(vapply(candidates, `[[`,
    numeric(1), "loglik"))]
  best <- candidates[[best_name]]
  actual_q <- (sum(df$reset12) + sum(df$reset23)) /
    sum((df$h1 == 1L & df$h2 == 1L) +
      (df$h2 == 1L & df$h3 == 1L))

  row_id <- row_id + 1L
  rows[[row_id]] <- data.frame(n=n, replication=r, seed=seed,
    best_start=best_name, convergence=best$convergence,
    iterations=best$iterations, loglik=best$loglik, actual_reset=actual_q,
    job_change_prob=best$params$job_change_prob, pi=best$params$pi,
    eps=best$params$eps, eps_d=best$params$eps_d,
    alpha=best$params$alpha,
    lambda_g_rmse=sqrt(mean((best$params$lambda_g-truth$lambda_g)^2)),
    lambda_d_rmse=sqrt(mean((best$params$lambda_d-truth$lambda_d)^2)),
    score_q=evaluate_score(df, "job_change_prob", .00025),
    score_pi=evaluate_score(df, "pi", .0005),
    score_eps=evaluate_score(df, "eps", .002),
    score_eps_d=evaluate_score(df, "eps_d", .002))
  write.csv(do.call(rbind, rows), file.path(outdir, "replications_latest.csv"),
    row.names=FALSE)
}

results <- do.call(rbind, rows)
parameters <- c(job_change_prob=truth$job_change_prob, pi=truth$pi,
  eps=truth$eps, eps_d=truth$eps_d, alpha=truth$alpha)
summary_rows <- lapply(sort(unique(results$n)), function(n) {
  z <- results[results$n == n, ]
  do.call(rbind, lapply(names(parameters), function(p) data.frame(
    n=n, replications=nrow(z), parameter=p, truth=parameters[[p]],
    mean=mean(z[[p]]), bias=mean(z[[p]])-parameters[[p]],
    rmse=sqrt(mean((z[[p]]-parameters[[p]])^2)),
    minimum=min(z[[p]]), maximum=max(z[[p]]))))
})
recovery_summary <- do.call(rbind, summary_rows)
score_summary <- aggregate(results[c("score_q", "score_pi", "score_eps",
  "score_eps_d")], list(n=results$n), function(x)
    c(mean=mean(x), sd=sd(x), rms=sqrt(mean(x^2))))

cat("\nRepeated monthly-model recovery\n")
print(results, row.names=FALSE, digits=8)
cat("\nBias and RMSE by sample size\n")
print(recovery_summary, row.names=FALSE, digits=8)
cat("\nScore per observation at the generating parameters\n")
print(score_summary, row.names=FALSE, digits=8)

write.csv(recovery_summary, file.path(outdir, "summary_latest.csv"),
  row.names=FALSE)
write.csv(score_summary, file.path(outdir, "score_summary_latest.csv"),
  row.names=FALSE)
saveRDS(list(truth=truth, results=results, summary=recovery_summary,
  score_summary=score_summary), file.path(outdir, "recovery_grid_latest.rds"))
