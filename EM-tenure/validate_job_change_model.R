# Numerical stability, nuisance-adjusted profile, and full-likelihood recovery
# checks for the Panel A job-change/reset model. The established fit is read but
# never overwritten. No bootstrap is run.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
main_file <- "EM-tenure/output/results/job_change/fit_latest.rds"
if (!file.exists(main_file)) stop("Estimate the Panel A job-change model first")
main <- readRDS(main_file)
outdir <- "EM-tenure/output/results/job_change/validation"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("JOB_CHANGE_WORKERS", "1"))
resume <- identical(tolower(Sys.getenv("JOB_CHANGE_VALIDATION_RESUME", "true")),
  "true")

cached_fit <- function(path, expression) {
  if (resume && file.exists(path)) return(readRDS(path))
  value <- force(expression)
  value$objective_function <- NULL
  saveRDS(value, path)
  value
}

# Two deliberately dispersed starts alter q, error probabilities, and hazards.
starts <- list(low=main$params, high=main$params)
starts$low$job_change_prob <- .005
starts$low$pi <- .012
starts$low$eps <- .15
starts$low$eps_d <- .10
starts$low$lambda_g <- starts$low$lambda_g * c(.80, 1.20, .85, 1.15, .90)
starts$low$lambda_d <- starts$low$lambda_d * c(1.20, .80, 1.15, .85, 1.10)
starts$high$job_change_prob <- .08
starts$high$pi <- .05
starts$high$eps <- .38
starts$high$eps_d <- .28
starts$high$lambda_g <- starts$high$lambda_g * c(1.25, .75, 1.20, .80, 1.10)
starts$high$lambda_d <- starts$high$lambda_d * c(.75, 1.25, .80, 1.20, .90)

message("Running dispersed Panel A starts")
fits <- lapply(names(starts), function(label) cached_fit(
  file.path(outdir, paste0("start_", label, ".rds")),
  fit_eps_piecewise_job_change(df_fit, starts[[label]],
    q_start=starts[[label]]$job_change_prob, maxit=50L,
    reltol=1e-9, pgtol=1e-7, gradient_step=1e-4,
    workers=workers, verbose=1L)))
names(fits) <- names(starts)

start_summary <- do.call(rbind, c(list(data.frame(start="established",
  loglik=main$loglik, convergence=main$convergence,
  iterations=main$iterations, job_change_prob=main$params$job_change_prob,
  pi=main$params$pi, eps=main$params$eps, eps_d=main$params$eps_d)),
  lapply(names(fits), function(label) {
    f <- fits[[label]]
    data.frame(start=label, loglik=f$loglik, convergence=f$convergence,
      iterations=f$iterations, job_change_prob=f$params$job_change_prob,
      pi=f$params$pi, eps=f$params$eps, eps_d=f$params$eps_d)
  })))
best_pool <- c(list(main), fits)
best <- best_pool[[which.max(vapply(best_pool, `[[`, numeric(1), "loglik"))]]

message("Tight refinement of the best solution")
refined <- cached_fit(file.path(outdir, "tight_refinement.rds"),
  fit_eps_piecewise_job_change(df_fit, best$params,
    q_start=best$params$job_change_prob, maxit=50L,
    reltol=1e-11, pgtol=1e-8, gradient_step=5e-5,
    workers=workers, verbose=1L))
if (!is.finite(refined$loglik) || refined$loglik < best$loglik - 1e-5)
  stop("Tight refinement is invalid or inferior")

# Nuisance-adjusted local profile. Each noncentral point reoptimizes all 14
# nuisance parameters at fixed q and starts from the tight optimum.
message("Profiling the reset probability with nuisance reoptimization")
q_hat <- refined$params$job_change_prob
offsets <- c(-.0015, -.001, -.0005, 0, .0005, .001, .0015)
q_grid <- q_hat + offsets
profile_fits <- lapply(seq_along(q_grid), function(k) {
  if (offsets[k] == 0) return(refined)
  tag <- sprintf("%+.4f", offsets[k])
  tag <- gsub("\\+", "p", gsub("-", "m", tag))
  cached_fit(file.path(outdir, paste0("profile_", tag, ".rds")),
    fit_eps_piecewise_job_change(df_fit, refined$params,
      q_fixed=q_grid[k], maxit=35L, reltol=1e-10, pgtol=1e-7,
      gradient_step=7.5e-5, workers=workers, verbose=0L))
})
profile <- data.frame(q=q_grid,
  loglik=vapply(profile_fits, `[[`, numeric(1), "loglik"),
  convergence=vapply(profile_fits, `[[`, integer(1), "convergence"),
  iterations=vapply(profile_fits, function(z) as.integer(z$iterations), integer(1)))
profile$lr <- 2 * (max(profile$loglik) - profile$loglik)
x <- profile$q - q_hat
quad <- lm(loglik ~ x + I(x^2), data=profile)
co <- coef(quad)
target <- max(profile$loglik) - qchisq(.95, 1) / 2
roots <- polyroot(c(co[1] - target, co[2], co[3]))
roots <- sort(Re(roots[abs(Im(roots)) < 1e-8])) + q_hat
if (length(roots) != 2L || roots[1] < min(q_grid) || roots[2] > max(q_grid))
  stop("Profile grid does not bracket the quadratic 95% interval")
profile_interval <- data.frame(estimate=q_hat, lower=roots[1], upper=roots[2],
  profile_points=nrow(profile), quadratic_R2=summary(quad)$r.squared)

# Full observed-likelihood recovery. A modest number of medium-sized replicates
# tests all parameters while keeping the validation reproducible and resumable.
message("Running full-likelihood recovery simulations")
sim_n <- as.integer(Sys.getenv("JOB_CHANGE_SIM_N", "20000"))
sim_reps <- as.integer(Sys.getenv("JOB_CHANGE_SIM_REPS", "3"))
sim_results <- lapply(seq_len(sim_reps), function(r) {
  d <- simulate_eps_piecewise_job_change(sim_n, refined$params,
    seed=17030L + r)
  actual_reset_rate <-
    (sum(d$reset12) + sum(d$reset23)) /
    sum((d$h1 == 1L & d$h2 == 1L) + (d$h2 == 1L & d$h3 == 1L))
  fit_one <- function(label, q_start) {
    cache <- file.path(outdir,
      paste0("simulation_", r, "_", label, "_fit.rds"))
    if (resume && file.exists(cache)) return(readRDS(cache))
    # Reuse the completed legacy run when it used this starting side.
    legacy <- file.path(outdir, paste0("simulation_", r, ".rds"))
    legacy_label <- if (r %% 2L) "low" else "high"
    if (resume && label == legacy_label && file.exists(legacy)) {
      value <- readRDS(legacy)$fit
      saveRDS(value, cache)
      return(value)
    }
    start <- refined$params
    start$job_change_prob <- q_start
    start$pi <- min(.10, refined$params$pi * 1.35)
    start$eps <- min(.60, refined$params$eps * .80)
    start$eps_d <- min(.50, refined$params$eps_d * 1.20)
    value <- fit_eps_piecewise_job_change(d, start,
      q_start=q_start, maxit=50L, reltol=1e-9, pgtol=1e-7,
      gradient_step=1e-4, workers=workers, verbose=0L)
    value$objective_function <- NULL
    saveRDS(value, cache)
    value
  }
  candidates <- list(low=fit_one("low", .012),
    high=fit_one("high", .045))
  best_name <- names(candidates)[which.max(vapply(candidates, `[[`,
    numeric(1), "loglik"))]
  ans <- list(fit=candidates[[best_name]], best_start=best_name,
    candidates=candidates, actual_reset_rate=actual_reset_rate)
  saveRDS(ans, file.path(outdir,
    paste0("simulation_", r, "_multistart.rds")))
  ans
})
simulation <- do.call(rbind, lapply(seq_along(sim_results), function(r) {
  z <- sim_results[[r]]; p <- z$fit$params
  data.frame(replicate=r, n=sim_n, best_start=z$best_start,
    convergence=z$fit$convergence,
    loglik=z$fit$loglik, actual_reset_rate=z$actual_reset_rate,
    job_change_prob=p$job_change_prob, pi=p$pi, eps=p$eps, eps_d=p$eps_d,
    lambda_g_rmse=sqrt(mean((p$lambda_g-refined$params$lambda_g)^2)),
    lambda_d_rmse=sqrt(mean((p$lambda_d-refined$params$lambda_d)^2)))
}))
truth <- data.frame(job_change_prob=refined$params$job_change_prob,
  pi=refined$params$pi, eps=refined$params$eps, eps_d=refined$params$eps_d)
recovery_summary <- data.frame(parameter=names(truth),
  truth=as.numeric(unlist(truth[1, ], use.names=FALSE)),
  mean_estimate=vapply(names(truth), function(nm)
    mean(simulation[[nm]]), numeric(1)),
  bias=vapply(names(truth), function(nm)
    mean(simulation[[nm]]) - truth[[nm]], numeric(1)),
  min_estimate=vapply(names(truth), function(nm)
    min(simulation[[nm]]), numeric(1)),
  max_estimate=vapply(names(truth), function(nm)
    max(simulation[[nm]]), numeric(1)))

# A large-sample score check distinguishes finite-sample recovery noise from a
# systematic mismatch between the intended DGP and the implemented likelihood.
score_file <- file.path(outdir, "large_sample_score_latest.csv")
if (resume && file.exists(score_file)) {
  score_diagnostic <- read.csv(score_file)
} else {
  score_n <- as.integer(Sys.getenv("JOB_CHANGE_SCORE_N", "100000"))
  score_data <- simulate_eps_piecewise_job_change(score_n, refined$params,
    seed=18001L)
  evaluate_score <- function(q, eps) {
    p <- refined$params
    p$job_change_prob <- q
    p$eps <- eps
    e_step_eps(score_data, p, check_df=FALSE, suff_stats=FALSE)$loglik
  }
  hq <- .00025; he <- .002
  qm <- evaluate_score(q_hat-hq, refined$params$eps)
  qp <- evaluate_score(q_hat+hq, refined$params$eps)
  em <- evaluate_score(q_hat, refined$params$eps-he)
  ep <- evaluate_score(q_hat, refined$params$eps+he)
  score_diagnostic <- data.frame(n=score_n,
    q_score=(qp-qm)/(2*hq), eps_score=(ep-em)/(2*he),
    q_score_per_observation=(qp-qm)/(2*hq*score_n),
    eps_score_per_observation=(ep-em)/(2*he*score_n))
  write.csv(score_diagnostic, score_file, row.names=FALSE)
}

cat("\nDispersed starts\n")
print(start_summary, row.names=FALSE, digits=9)
cat("\nTight optimum\n")
print(data.frame(loglik=refined$loglik,
  job_change_prob=refined$params$job_change_prob,
  convergence=refined$convergence, iterations=refined$iterations),
  row.names=FALSE, digits=9)
cat("\nNuisance-adjusted profile\n")
print(profile, row.names=FALSE, digits=9)
print(profile_interval, row.names=FALSE, digits=9)
cat("\nFull-likelihood simulation estimates\n")
print(simulation, row.names=FALSE, digits=8)
cat("\nRecovery summary\n")
print(recovery_summary, row.names=FALSE, digits=8)
cat("\nLarge-sample score diagnostic at the generating parameters\n")
print(score_diagnostic, row.names=FALSE, digits=8)

write.csv(start_summary, file.path(outdir, "dispersed_starts_latest.csv"),
  row.names=FALSE)
write.csv(profile, file.path(outdir, "profile_latest.csv"), row.names=FALSE)
write.csv(profile_interval, file.path(outdir, "profile_interval_latest.csv"),
  row.names=FALSE)
write.csv(simulation, file.path(outdir, "simulation_latest.csv"), row.names=FALSE)
write.csv(recovery_summary, file.path(outdir, "recovery_summary_latest.csv"),
  row.names=FALSE)
saveRDS(list(main=main, starts=fits, refined=refined, profile=profile,
  profile_interval=profile_interval, simulation=simulation,
  recovery_summary=recovery_summary, score_diagnostic=score_diagnostic),
  file.path(outdir, "validation_latest.rds"))
