# Table 5 robustness: Set 2 AR(2) model with piecewise-constant duration hazards.
if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source EM-AR2/estimate_set2_piecewise_duration_4w.R from the project root")
source("EM-AR1-4W/R/source_all.R")
source("EM-AR2/R/source_all.R")

set.seed(20260828L)
results_dir <- "EM-AR2/output/results/set2_piecewise_duration"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
data <- prepare_ar2_set2_piecewise_4w()
cat(sprintf("Piecewise-duration Table 5 sample: N=%s; exact cells=%s; K=%d\n",
  format(data$n_original, big.mark = ","), format(data$n_cells, big.mark = ","),
  length(ar2_reliability_names(data))))

linear_fit <- readRDS(paste0(
  "EM-AR2/output/results/set4_reliability/",
  "fit_ar2_set4_constant_latest.rds"))

piecewise_start <- function() {
  out <- setNames(rep(0, length(ar2_reliability_names(data))),
                  ar2_reliability_names(data))
  common <- intersect(names(linear_fit$eta), names(out))
  out[common] <- linear_fit$eta[common]
  representative_months <- c(4.5, 9, 18, 42, 90)
  add_duration_warm <- function(block, old_term, prefix, scaling) {
    slope <- unname(linear_fit$eta[old_term])
    reference <- (log1p(1.5 / 12) - scaling$center) / scaling$scale
    midpoint <- (log1p(representative_months / 12) - scaling$center) /
      scaling$scale
    out[paste0(block, "_intercept")] <<-
      out[paste0(block, "_intercept")] + slope * reference
    bin_terms <- grep(paste0("^", block, "_", prefix, "_"), names(out), value = TRUE)
    out[bin_terms] <<- slope * (midpoint - reference)
  }
  add_duration_warm("entry", "entry_log_time_since_work", "timegap",
                    data$scaling$log_time_since_work)
  add_duration_warm("persistence", "persistence_log_tenure", "tenure",
                    data$scaling$log_tenure)
  out
}

base_start <- piecewise_start()
n_starts <- as.integer(Sys.getenv("AR2_PIECEWISE_STARTS", "2"))
if (!is.finite(n_starts) || n_starts < 1L || n_starts > 5L)
  stop("AR2_PIECEWISE_STARTS must be between 1 and 5")
starts <- list(linear_warm = base_start)
if (n_starts >= 2L) for (j in 2:n_starts) {
  scale <- ifelse(grepl("timegap_|tenure_", names(base_start)), .12,
                  ifelse(grepl("lag2|intercept", names(base_start)), .08, .04))
  starts[[paste0("perturbed_", j - 1L)]] <-
    base_start + rnorm(length(base_start), 0, scale)
}

resume <- identical(Sys.getenv("AR2_PIECEWISE_RESUME"), "1")
prelim_maxit <- as.integer(Sys.getenv("AR2_PIECEWISE_PRELIM_MAXIT", "250"))
fits <- vector("list", length(starts)); names(fits) <- names(starts)
elapsed <- numeric(length(starts))
for (j in seq_along(starts)) {
  checkpoint <- file.path(results_dir,
    sprintf("fit_piecewise_start_%s.rds", names(starts)[j]))
  if (resume && file.exists(checkpoint)) {
    fits[[j]] <- readRDS(checkpoint)
    cat(sprintf("resumed %-14s ll=%.3f score=%.3e\n", names(starts)[j],
                fits[[j]]$loglik, fits[[j]]$max_abs_score))
  } else {
    tm <- system.time(fits[[j]] <- fit_ar2_set4_reliability(
      data, starts[[j]], maxit = prelim_maxit, reltol = 1e-7))
    elapsed[j] <- unname(tm["elapsed"])
    saveRDS(fits[[j]], checkpoint)
    cat(sprintf("%-22s ll=%.3f code=%d eval=%d score=%.3e time=%.1fs\n",
      names(starts)[j], fits[[j]]$loglik, fits[[j]]$optimizer_code,
      fits[[j]]$iterations, fits[[j]]$max_abs_score, elapsed[j]))
  }
}
multistart <- data.frame(start = names(fits),
  loglik = vapply(fits, `[[`, numeric(1), "loglik"),
  code = vapply(fits, `[[`, integer(1), "optimizer_code"),
  evaluations = vapply(fits, `[[`, numeric(1), "iterations"),
  max_abs_score = vapply(fits, `[[`, numeric(1), "max_abs_score"),
  elapsed_seconds = elapsed)
write.csv(multistart, file.path(results_dir, "multistart_piecewise.csv"),
          row.names = FALSE)
fit <- fits[[which.max(multistart$loglik)]]

cat("Computing analytical inference for the piecewise-duration model...\n")
inference <- analytical_se_ar2_set4_reliability(data, fit)
correction <- inference$diagnostics$newton_correction
if (max(abs(correction)) > 1e-3) {
  fractions <- c(1, .5, .25, .125, 0)
  trials <- lapply(fractions, function(frac) fit$eta + frac * correction)
  likelihood <- vapply(trials, function(z)
    sum(data$weight * .ar2r_components(z, data)$loglik_i), numeric(1))
  best <- which.max(likelihood)
  if (likelihood[best] > fit$loglik + 1e-5) {
    cat(sprintf("Newton polish fraction %.3f; gain %.3f\n",
                fractions[best], likelihood[best] - fit$loglik))
    saveRDS(list(eta = trials[[best]], fraction = fractions[best],
      loglik = likelihood[best], source_fit = fit$eta, inference = inference),
      file.path(results_dir, "newton_start_piecewise_latest.rds"))
    fit <- fit_ar2_set4_reliability(data, trials[[best]], 500L, 1e-10)
    inference <- analytical_se_ar2_set4_reliability(data, fit)
  }
}
if (fit$max_abs_score > 2e-5)
  stop("Final piecewise-duration normalized score remains too large")
if (inference$diagnostics$information_rank <
    inference$diagnostics$parameter_count)
  stop("Piecewise-duration model is not locally identified")

fit$analytical_inference <- inference
fit$multistart <- multistart
saveRDS(fit, file.path(results_dir, "fit_ar2_set2_piecewise_latest.rds"))
write.csv(inference$coefficient_summary,
  file.path(results_dir, "coefficients_ar2_set2_piecewise.csv"), row.names = FALSE)
write.csv(inference$probability_summary,
  file.path(results_dir, "probabilities_ar2_set2_piecewise.csv"), row.names = FALSE)

constant <- linear_fit$analytical_inference$probability_summary
comparison <- rbind(
  transform(constant, model = "Set 2, log-linear duration"),
  transform(inference$probability_summary,
            model = "Set 2, piecewise-constant duration"))
comparison <- comparison[, c("model", "quantity", "estimate", "se")]
write.csv(comparison, file.path(results_dir, "table5_piecewise_comparison.csv"),
          row.names = FALSE)

cat(sprintf(paste0("accepted piecewise duration: ll=%.3f score=%.3e ",
  "rank=%d/%d cond=%.3g\n"), fit$loglik, fit$max_abs_score,
  inference$diagnostics$information_rank,
  inference$diagnostics$parameter_count,
  inference$diagnostics$information_condition))
print(inference$probability_summary, row.names = FALSE, digits = 6)
