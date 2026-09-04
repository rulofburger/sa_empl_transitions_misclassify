# Joint likelihood continuation only. No three-wave comparator is estimated.
source("EM-tenure/R/source_all.R")
source("EM-tenure/R/four_wave_fast.R")
outdir <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
validation <- readRDS(file.path(outdir, "compiled_equivalence_latest.rds"))
fingerprint <- four_wave_fast_source_fingerprint()
if (!identical(validation$source_md5, fingerprint))
  stop("Run validate_four_wave_fast.R for the current compiled source first")
recovery <- readRDS(file.path(outdir,
  "recovery_continuous_clock/recovery_status_latest.rds"))
if (!identical(validation$timegap_clock, "continuous_joint") ||
    !identical(recovery$source_md5, fingerprint) ||
    !isTRUE(recovery$ready_for_recovery_fits) ||
    !isTRUE(recovery$recovery_fits_all_converged))
  stop("Current-source repaired-clock validation and recovery must pass first")
df <- readRDS(file.path(outdir, "prepared_cells_latest.rds"))$df4
data <- prepare_four_wave_kernel_data(df)
workers <- as.integer(Sys.getenv("FOUR_WAVE_FAST_WORKERS", "8"))
workers <- max(1L, min(8L, workers))
maxit <- as.integer(Sys.getenv("FOUR_WAVE_FAST_MAXIT", "80"))
label <- Sys.getenv("FOUR_WAVE_FAST_LABEL", "clock_joint1")
gradient_scheme <- Sys.getenv("FOUR_WAVE_FAST_GRADIENT", "central")
if (!gradient_scheme %in% c("central", "forward")) stop("Invalid gradient scheme")
checkpoint <- file.path(outdir, paste0("checkpoint_fast_", label, ".rds"))
output <- file.path(outdir, paste0("fit_fast_", label, "_latest.rds"))
start_path <- Sys.getenv("FOUR_WAVE_FAST_START",
  file.path(outdir, "fit_fast_start_plus_latest.rds"))
start <- readRDS(start_path)$params
z0 <- .pack_four_wave_preferred(start)
is_resuming <- file.exists(checkpoint) &&
  identical(Sys.getenv("FOUR_WAVE_FAST_RESUME", "true"), "true")
if (is_resuming) {
  resumed <- readRDS(checkpoint)
  if (!identical(resumed$source_md5,fingerprint))
    stop("Checkpoint sources differ: use a new label and an explicit starting fit")
  saveRDS(resumed, file.path(outdir, paste0("resume_fast_", label, "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")))
  z0 <- resumed$par
}
perturb <- as.numeric(Sys.getenv("FOUR_WAVE_FAST_PERTURB", "0"))
if (is_resuming) perturb <- 0 # Do not perturb an already-refined checkpoint again.
bounds <- .four_wave_parameter_bounds(z0)
z0 <- pmin(bounds$upper, pmax(bounds$lower, z0 + perturb * sin(seq_along(z0))))
weight <- sum(df$weight)
information_path <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/inference_opg/opg_inference_latest.rds")
information <- readRDS(information_path)$information
# A diagonal numerical preconditioner only: it is not four-wave inference.
# Cap nearly flat reporting directions instead of inverting their tiny curvature.
parscale <- pmin(10, pmax(.5, sqrt(weight / diag(information)[names(z0)])))
if (any(!is.finite(parscale))) stop("Invalid saved-information parameter scaling")
score_tolerance <- as.numeric(Sys.getenv("FOUR_WAVE_FAST_SCORE_TOL", "1e-5"))
if (!is.finite(score_tolerance) || score_tolerance <= 0) stop("Invalid score tolerance")
# optim applies pgtol to g * parscale. This sufficient scaled threshold
# guarantees the predeclared unscaled criterion for every coordinate.
optimizer_pgtol <- score_tolerance * min(parscale)
worker_objective <- function(z) {
  p <- .piecewise_calendar_revision_monthly_unpack(z, "joint_marginal")
  v <- evaluate_four_wave_monthly_fast(data, p, posterior = FALSE, threads = 1L)$loglik
  if (!is.finite(v)) stop("Nonfinite likelihood in joint refinement")
  -v / weight
}
cluster <- parallel::makePSOCKcluster(workers)
# Shutdown is handled by tryCatch(finally=...) below.
initialize_worker <- function(path) {
  setwd(path)
  source("EM-tenure/R/source_all.R")
  source("EM-tenure/R/four_wave_fast.R")
  load_four_wave_monthly_kernel()
  NULL
}
# Rcpp's shared build-cache metadata must not be rewritten concurrently.
tryCatch(for (node in seq_along(cluster))
  parallel::clusterCall(cluster[node], initialize_worker, getwd()),
  error = function(e) {parallel::stopCluster(cluster); stop(e)})
parallel::clusterExport(cluster, c("data", "weight", "worker_objective"),
  envir = environment())
best <- Inf
evaluations <- 0L
gradient_calls <- 0L
last_z <- last_value <- NULL
history <- data.frame()
objective <- function(z) {
  if (identical(z, last_z)) return(last_value)
  p <- .piecewise_calendar_revision_monthly_unpack(z, "joint_marginal")
  value <- -evaluate_four_wave_monthly_fast(data, p, posterior = FALSE,
    threads = 8L)$loglik / weight
  if (!is.finite(value)) stop("Nonfinite likelihood in joint refinement")
  evaluations <<- evaluations + 1L
  last_z <<- z
  last_value <<- value
  if (value < best) {
    best <<- value
    history <<- rbind(history, data.frame(evaluation = evaluations,
      gradient_calls = gradient_calls, loglik = -value * weight,
      timestamp = as.character(Sys.time())))
    saveRDS(list(params = p, par = z, loglik = -value * weight,
      convergence = 1L, converged = FALSE, history = history,
      source_md5 = fingerprint), checkpoint)
    message(sprintf("%s %s: evaluation %d, LL %.6f", Sys.time(), label,
      evaluations, -value * weight))
  }
  value
}
gradient <- function(z, central = identical(gradient_scheme, "central")) {
  step <- 1e-5 * pmax(1, abs(z))
  plus <- pmin(step, bounds$upper - z)
  minus <- pmin(step, z - bounds$lower)
  if (central) {
    points <- vector("list", 2L * length(z))
    for (j in seq_along(z)) {
      points[[2*j-1L]] <- points[[2*j]] <- z
      points[[2*j-1L]][j] <- z[j] + plus[j]
      points[[2*j]][j] <- z[j] - minus[j]
    }
  } else {
    step <- ifelse(plus >= step, step, -minus)
    points <- lapply(seq_along(z), function(j) {
      value <- z; value[j] <- z[j] + step[j]; value
    })
  }
  values <- unlist(parallel::parLapplyLB(cluster, points, worker_objective,
    chunk.size = 1L),
    use.names = FALSE)
  ans <- if (central) (values[seq(1, length(values), 2)] -
    values[seq(2, length(values), 2)]) / (plus + minus) else
    (values - objective(z)) / step
  if (any(!is.finite(ans))) stop("Nonfinite numerical score")
  gradient_calls <<- gradient_calls + 1L
  message(sprintf("%s %s: gradient %d, max |score/N| %.8g", Sys.time(),
    label, gradient_calls, max(abs(ans))))
  ans
}
result <- tryCatch({
  initial <- objective(z0)
  opt <- optim(z0, objective, gr = gradient, method = "L-BFGS-B",
    lower = bounds$lower, upper = bounds$upper,
    control = list(maxit = maxit, factr = 1e5, pgtol = optimizer_pgtol, lmm = 33,
      parscale = parscale))
  score <- gradient(opt$par, central = TRUE)
  projected <- score
  projected[opt$par <= bounds$lower + 1e-8 & score > 0] <- 0
  projected[opt$par >= bounds$upper - 1e-8 & score < 0] <- 0
  p <- .piecewise_calendar_revision_monthly_unpack(opt$par, "joint_marginal")
  p$timegap_clock_model <- "continuous_joint"
  final <- evaluate_four_wave_monthly_fast(data, p, threads = 8L)
  fit <- c(list(params = p, par_unconstrained = opt$par,
    integration_method = "exact_piecewise",
    convergence = opt$convergence, message = opt$message,
    independent_score_pass = max(abs(projected)) <= score_tolerance,
    converged = opt$convergence == 0L && max(abs(projected)) <= score_tolerance,
    central_score = score, projected_score = projected,
    gradient_max = max(abs(projected)), counts = opt$counts,
    loglik_gain = final$loglik + initial * weight,
    history = history, source_md5 = fingerprint,
    controls = list(maxit = maxit, workers = workers, step = 1e-5,
      gradient_scheme = gradient_scheme,
      score_tolerance = score_tolerance, pgtol = optimizer_pgtol,
      parscale = parscale, scaling_source_md5 = tools::md5sum(information_path)),
    session_info = sessionInfo()), final)
  saveRDS(fit, output)
  write.csv(history, file.path(outdir, paste0("history_fast_", label, ".csv")),
    row.names = FALSE)
  print(c(loglik = fit$loglik, gain = fit$loglik_gain,
    convergence = fit$convergence, projected_score = fit$gradient_max))
  print(duration_weighted_transition_rates_4w(df, fit))
  fit
}, finally = parallel::stopCluster(cluster))
