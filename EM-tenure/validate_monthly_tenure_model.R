# Validate the discrete-month tenure likelihood against its own data-generating
# process. This is deliberately separate from the empirical estimator so a
# failed recovery check cannot be mistaken for a substantive result.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")

old <- readRDS("EM-tenure/output/results/job_change/fit_latest.rds")
truth <- old$params
truth$tenure_measurement_model <- "monthly"

outdir <- "EM-tenure/output/results/job_change_monthly/validation"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
n <- as.integer(Sys.getenv("MONTHLY_VALIDATION_N", "50000"))
workers <- as.integer(Sys.getenv("MONTHLY_VALIDATION_WORKERS", "1"))

sim <- simulate_eps_piecewise_job_change(n, truth, seed=25041L)
validate_df_eps(sim, allow_zero_tenure=TRUE)

evaluate <- function(p) e_step_eps(sim, p, suff_stats=FALSE)$loglik
central_score <- function(name, step) {
  lower <- upper <- truth
  lower[[name]] <- truth[[name]] - step
  upper[[name]] <- truth[[name]] + step
  (evaluate(upper) - evaluate(lower)) / (2 * step)
}
score <- data.frame(
  parameter=c("job_change_prob", "eps", "pi"),
  truth=unlist(truth[c("job_change_prob", "eps", "pi")]),
  step=c(.00025, .002, .0005))
score$score <- mapply(central_score, score$parameter, score$step)
score$score_per_observation <- score$score / n

# Estimate from starts on either side of the truth. The better likelihood is
# retained; agreement between starts is reported rather than assumed.
fit_one <- function(q_start) {
  start <- truth
  start$job_change_prob <- q_start
  start$eps <- if (q_start < truth$job_change_prob) .18 else .30
  start$pi <- .03
  fit_eps_piecewise_job_change(sim, start, q_start=q_start, maxit=60L,
    reltol=1e-9, pgtol=1e-7, workers=workers,
    tenure_measurement_model="monthly")
}
fits <- list(low=fit_one(.008), high=fit_one(.055))
fit_summary <- do.call(rbind, lapply(names(fits), function(label) {
  f <- fits[[label]]
  data.frame(start=label, loglik=f$loglik, convergence=f$convergence,
    iterations=f$iterations, job_change_prob=f$params$job_change_prob,
    pi=f$params$pi, eps=f$params$eps, eps_d=f$params$eps_d)
}))
best <- fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
recovery <- data.frame(parameter=c("job_change_prob", "pi", "eps", "eps_d"),
  truth=unlist(truth[c("job_change_prob", "pi", "eps", "eps_d")]),
  estimate=unlist(best$params[c("job_change_prob", "pi", "eps", "eps_d")]))
recovery$bias <- recovery$estimate - recovery$truth

cat("\nLarge-sample score at generating parameters\n")
print(score, row.names=FALSE, digits=8)
cat("\nFull-likelihood recovery starts\n")
print(fit_summary, row.names=FALSE, digits=8)
cat("\nBest-fit recovery\n")
print(recovery, row.names=FALSE, digits=8)

write.csv(score, file.path(outdir, "score_latest.csv"), row.names=FALSE)
write.csv(fit_summary, file.path(outdir, "fits_latest.csv"), row.names=FALSE)
write.csv(recovery, file.path(outdir, "recovery_latest.csv"), row.names=FALSE)
saveRDS(list(truth=truth, score=score, fit_summary=fit_summary,
  recovery=recovery), file.path(outdir, "validation_latest.rds"))
