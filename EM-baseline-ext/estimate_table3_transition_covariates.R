# Estimate the Table 3 transition-covariate specifications. The preferred model
# adds origin-wave inconsistencies to the Set 2 entry and persistence equations;
# it is estimated with both free and conditionally stationary initial states.

library(here)
library(dplyr)
library(ggplot2)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

results_dir <- here::here("EM-baseline-ext", "output", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260827L)

# Use exactly the current Table 3 sample.
keep <- complete.cases(df_qlfs[, c(
  "y1", "y2", "y3", "weight", "age1", "age2", "age3",
  "educ1", "educ2", "educ3", "race1", "female1")]) &
  df_qlfs$weight > 0
df <- as.data.frame(df_qlfs[keep, , drop = FALSE])
df$y1 <- as.integer(df$y1); df$y2 <- as.integer(df$y2)
df$y3 <- as.integer(df$y3); df$weight <- as.numeric(df$weight)
stopifnot(nrow(df) == 394823L)

inc <- add_inconsistency_count_dummies(
  compute_demographic_inconsistencies(compute_inconsistencies(df)))
characteristics <- c("age", "education", "race", "gender")
inc_prefix <- c(age = "age", education = "edu", race = "race", gender = "gender")
Z <- lapply(1:3, function(tt) {
  out <- cbind(error_intercept = 1,
    age_inconsistency = inc[[paste0("Y_age_", tt)]],
    education_inconsistency = inc[[paste0("Y_edu_", tt)]],
    race_inconsistency = inc[[paste0("Y_race_", tt)]],
    gender_inconsistency = inc[[paste0("Y_gender_", tt)]],
    two_inconsistencies = inc[[paste0("Y_exactly_2_", tt)]],
    three_inconsistencies = inc[[paste0("Y_exactly_3_", tt)]],
    four_inconsistencies = inc[[paste0("Y_exactly_4_", tt)]])
  storage.mode(out) <- "double"
  out
})

X2 <- prepare_covariate_matrix(df, covariate_set = 2L)$X_transition

append_transition_inconsistencies <- function(X, inc) {
  old_entry <- attr(X$X12, "entry_active")
  old_persistence <- attr(X$X12, "persistence_active")
  add <- lapply(1:2, function(tt) {
    out <- sapply(names(inc_prefix), function(nm)
      inc[[paste0("Y_", inc_prefix[[nm]], "_", tt)]])
    colnames(out) <- paste0("origin_", names(inc_prefix), "_inconsistency")
    storage.mode(out) <- "double"
    out
  })
  out <- list(X12 = cbind(X$X12, add[[1L]]),
              X23 = cbind(X$X23, add[[2L]]))
  for (tt in 1:2) {
    attr(out[[tt]], "entry_active") <- c(old_entry, rep(TRUE, 4L))
    attr(out[[tt]], "persistence_active") <- c(old_persistence, rep(TRUE, 4L))
  }
  out
}
X2_inc <- append_transition_inconsistencies(X2, inc)

cells4 <- collapse_covariate_reliability(df, X2, Z)
cells5 <- collapse_covariate_reliability(df, X2_inc, Z)
cat(sprintf("Table 3 Set 2 specification: N=%s; %s cells; K=25\n",
  format(cells4$n_original, big.mark = ","),
  format(nrow(cells4$df), big.mark = ",")))
cat(sprintf("Table 3 transition-inconsistency specification: N=%s; %s cells; K=33\n",
  format(cells5$n_original, big.mark = ","),
  format(nrow(cells5$df), big.mark = ",")))

table3 <- readRDS(file.path(results_dir, "fit_table6_inconsistency_audit.rds"))
set2 <- readRDS(file.path(results_dir, "fit_cov_s2_sym_free.rds"))
direct_delta_names <- c(
  error_intercept = "delta0", age_inconsistency = "delta_age",
  education_inconsistency = "delta_education",
  race_inconsistency = "delta_race", gender_inconsistency = "delta_gender",
  two_inconsistencies = "delta_two_inconsistencies",
  three_inconsistencies = "delta_three_inconsistencies",
  four_inconsistencies = "delta_four_inconsistencies")
delta <- setNames(unname(table3$free$params$delta[direct_delta_names]),
                  names(direct_delta_names))
start4 <- list(beta0 = setNames(set2$params$beta0, colnames(X2$X12)),
               beta1 = setNames(set2$params$beta1, colnames(X2$X12)),
               alpha = set2$params$alpha, delta = delta)

fit_multistart <- function(data, center, stem, nested_loglik,
                           requested = 4L, transition_sd = .08,
                           delta_sd = .15, stationary = FALSE) {
  eta <- pack_covariate_reliability(
    center, data$X, stationary = stationary)
  q <- length(center$delta)
  transition_k <- length(eta) - q
  starts <- c(list(warm_start = eta), setNames(lapply(seq_len(requested - 1L),
    function(i) eta + rnorm(length(eta), 0,
      c(rep(transition_sd, transition_k), rep(delta_sd, q)))),
    paste0("perturbed_", seq_len(requested - 1L))))
  fits <- vector("list", length(starts)); names(fits) <- names(starts)
  elapsed <- numeric(length(starts))
  resume <- identical(Sys.getenv("TABLE3_COV_RESUME"), "1")
  for (i in seq_along(starts)) {
    checkpoint <- file.path(results_dir,
      paste0(stem, "_start_", names(starts)[i], ".rds"))
    cat(sprintf("\n%s: start %d/%d (%s)\n", stem, i, length(starts), names(starts)[i]))
    if (resume && file.exists(checkpoint)) {
      fits[[i]] <- readRDS(checkpoint)
      cat(sprintf("  resumed ll=%.3f code=%d max|score|=%.3e\n",
        fits[[i]]$loglik, fits[[i]]$optimizer_code, fits[[i]]$max_abs_score))
      next
    }
    timing <- system.time(fits[[i]] <- fit_covariate_reliability(
      data, starts[[i]], maxit = 2000L, reltol = 1e-9,
      stationary = stationary))
    elapsed[i] <- unname(timing["elapsed"])
    saveRDS(fits[[i]], checkpoint)
    cat(sprintf("  ll=%.3f code=%d evaluations=%d max|score|=%.3e elapsed=%.1fs\n",
      fits[[i]]$loglik, fits[[i]]$optimizer_code, fits[[i]]$iterations,
      fits[[i]]$max_abs_score, elapsed[i]))
  }
  best <- which.max(vapply(fits, `[[`, numeric(1L), "loglik"))
  fit <- fits[[best]]
  candidates <- data.frame(start = names(starts),
    loglik = vapply(fits, `[[`, numeric(1L), "loglik"),
    optimizer_code = vapply(fits, `[[`, integer(1L), "optimizer_code"),
    evaluations = vapply(fits, `[[`, numeric(1L), "iterations"),
    max_abs_score = vapply(fits, `[[`, numeric(1L), "max_abs_score"),
    elapsed_seconds = elapsed)
  # Conditional stationarity can be very flat near the optimum. Permit a
  # slightly looser preliminary score here; the Hessian-polishing and final
  # Newton-correction checks below remain strict.
  preliminary_score_limit <- if (stationary) 1e-4 else 2e-5
  if (fit$optimizer_code != 0L || fit$max_abs_score > preliminary_score_limit)
    stop(stem, ": best fit failed convergence/score checks")
  if (!is.null(nested_loglik) && fit$loglik < nested_loglik - 1e-4)
    stop(stem, ": likelihood does not dominate the nested specification")
  list(fit = fit, candidates = candidates)
}

requested_starts <- as.integer(Sys.getenv("TABLE3_COV_STARTS", "4"))
if (!requested_starts %in% 1:8)
  stop("TABLE3_COV_STARTS must be between 1 and 8")

reuse_audited_fits <- identical(
  Sys.getenv("TABLE3_COV_REUSE_AUDITED_FITS"), "1")
audited_path <- file.path(results_dir, "fit_table3_transition_covariates.rds")
if (reuse_audited_fits) {
  if (!file.exists(audited_path)) stop("Missing audited transition fit: ", audited_path)
  audited <- readRDS(audited_path)
  fit4 <- audited$column4
  fit4_run <- list(fit = fit4, candidates = audited$column4_candidates)
} else {
  fit4_run <- fit_multistart(cells4, start4, "table3_col4_covariates",
    nested_loglik = set2$loglik, requested = requested_starts)
  fit4 <- fit4_run$fit
}

start5 <- list(
  beta0 = c(fit4$params$beta0,
    setNames(rep(0, 4L), setdiff(colnames(X2_inc$X12), colnames(X2$X12)))),
  beta1 = c(fit4$params$beta1,
    setNames(rep(0, 4L), setdiff(colnames(X2_inc$X12), colnames(X2$X12)))),
  alpha = fit4$params$alpha, delta = fit4$params$delta)
if (reuse_audited_fits) {
  fit5 <- audited$column5
  fit5_run <- list(fit = fit5, candidates = audited$column5_candidates)
} else {
  fit5_run <- fit_multistart(cells5, start5,
    "table3_col5_covariates_inconsistency",
    nested_loglik = fit4$loglik, requested = requested_starts)
  fit5 <- fit5_run$fit
}

# Conditional stationarity sets each person's wave-1 latent employment
# probability to theta0_i / (theta0_i + 1 - theta1_i), using the wave-1 Set 2
# covariates and origin-wave inconsistency indicators. Later hazards may differ
# because those origin-wave indicators can change between transitions.
start5_stationary <- fit5$params
start5_stationary$alpha <- NULL
if (reuse_audited_fits) {
  fit5_stationary <- audited$column5_stationary
  fit5_stationary_run <- list(
    fit = fit5_stationary,
    candidates = audited$column5_stationary_candidates)
} else {
  fit5_stationary_run <- fit_multistart(
    cells5, start5_stationary, "table3_col5_stationary",
    nested_loglik = NULL, requested = requested_starts, stationary = TRUE)
  fit5_stationary <- fit5_stationary_run$fit
}

cat("\nCalculating analytical sandwich/delta inference for column (4)...\n")
inf4 <- analytical_se_covariate_reliability(cells4, fit4)
cat("Calculating analytical sandwich/delta inference for column (5)...\n")
inf5 <- analytical_se_covariate_reliability(cells5, fit5)
cat("Calculating analytical sandwich/delta inference for conditionally stationary column (5)...\n")
inf5_stationary <- analytical_se_covariate_reliability(
  cells5, fit5_stationary, stationary = TRUE)

polish_fit <- function(data, fit, inference, label, max_rounds = 3L,
                       stationary = fit$stationary %||% FALSE) {
  for (round in seq_len(max_rounds)) {
    correction <- inference$diagnostics$newton_correction
    size <- max(abs(correction))
    cat(sprintf("%s: polishing round %d, max Newton correction %.3e\n",
                label, round, size))
    if (size <= .002) break
    fractions <- c(1, .5, .25, .125, 0)
    trials <- lapply(fractions, function(frac) fit$eta + frac * correction)
    trial_ll <- vapply(trials, function(z)
      sum(data$df$weight * .covrel_components(
        z, data, stationary = stationary)$loglik_i), numeric(1L))
    best <- which.max(trial_ll)
    cat(sprintf("  line search fraction %.3f; likelihood gain %.3f\n",
                fractions[best], trial_ll[best] - fit$loglik))
    fit_new <- fit_covariate_reliability(
      data, trials[[best]], maxit = 400L, reltol = 1e-9,
      stationary = stationary)
    if (fit_new$loglik < fit$loglik - 1e-4)
      stop(label, ": polishing reduced the likelihood")
    fit <- fit_new
    inference <- analytical_se_covariate_reliability(
      data, fit, stationary = stationary)
  }
  list(fit = fit, inference = inference)
}
polished4 <- polish_fit(cells4, fit4, inf4, "Table 3 column (4)")
fit4 <- polished4$fit; inf4 <- polished4$inference
polished5 <- polish_fit(cells5, fit5, inf5, "Table 3 column (5)")
fit5 <- polished5$fit; inf5 <- polished5$inference
polished5_stationary <- polish_fit(
  cells5, fit5_stationary, inf5_stationary,
  "Table 3 column (5), conditional stationarity", stationary = TRUE)
fit5_stationary <- polished5_stationary$fit
inf5_stationary <- polished5_stationary$inference

for (obj in list(inf4, inf5, inf5_stationary)) {
  d <- obj$diagnostics
  if (d$information_rank != d$K)
    stop("Table 3 transition-covariate model is rank deficient")
  if (d$max_abs_newton_correction > .01)
    stop("Table 3 transition-covariate model needs material Newton correction")
}

append_polished <- function(candidates, fit) {
  candidates <- candidates[candidates$start != "hessian_polished", , drop = FALSE]
  rbind(candidates, data.frame(
    start = "hessian_polished", loglik = fit$loglik,
    optimizer_code = fit$optimizer_code, evaluations = fit$iterations,
    max_abs_score = fit$max_abs_score, elapsed_seconds = NA_real_))
}
fit4_run$candidates <- append_polished(fit4_run$candidates, fit4)
fit5_run$candidates <- append_polished(fit5_run$candidates, fit5)
fit5_stationary_run$candidates <- append_polished(
  fit5_stationary_run$candidates, fit5_stationary)

out <- list(
  column4 = fit4, column4_inference = inf4,
  column4_candidates = fit4_run$candidates,
  column5 = fit5, column5_inference = inf5,
  column5_candidates = fit5_run$candidates,
  column5_stationary = fit5_stationary,
  column5_stationary_inference = inf5_stationary,
  column5_stationary_candidates = fit5_stationary_run$candidates,
  design = list(column4 = colnames(X2$X12),
                column5 = colnames(X2_inc$X12),
                error = colnames(Z[[1L]])),
  N = nrow(df), cells = c(column4 = nrow(cells4$df), column5 = nrow(cells5$df)))
saveRDS(out, file.path(results_dir, "fit_table3_transition_covariates.rds"))
write.csv(fit4_run$candidates,
  file.path(results_dir, "table3_col4_multistart.csv"), row.names = FALSE)
write.csv(fit5_run$candidates,
  file.path(results_dir, "table3_col5_multistart.csv"), row.names = FALSE)
write.csv(fit5_stationary_run$candidates,
  file.path(results_dir, "table3_col5_stationary_multistart.csv"),
  row.names = FALSE)
write.csv(inf4$summary,
  file.path(results_dir, "table3_col4_summary.csv"), row.names = FALSE)
write.csv(inf5$summary,
  file.path(results_dir, "table3_col5_summary.csv"), row.names = FALSE)
write.csv(inf5_stationary$summary,
  file.path(results_dir, "table3_col5_stationary_summary.csv"),
  row.names = FALSE)

report <- function(inference, fit, label) {
  keep <- c("mean_entry_rate", "mean_exit_rate", "mean_misclassification_rate",
            "initial_employment", "pi_base", "pi_age_mild", "pi_education_mild",
            "pi_race", "pi_gender", paste0("pi_count_", 0:4))
  z <- inference$summary[inference$summary$quantity %in% keep, ]
  z$estimate_percent <- 100 * z$estimate
  z$std_error_percent <- 100 * z$std_error
  cat("\n", label, " (percent):\n", sep = "")
  print(z[, c("quantity", "estimate_percent", "std_error_percent")],
        row.names = FALSE, digits = 4)
  cat(sprintf("Log-likelihood %.3f; rank %d/%d; condition %.3g; max score %.3e\n",
    fit$loglik, inference$diagnostics$information_rank,
    inference$diagnostics$K, inference$diagnostics$information_condition,
    inference$diagnostics$max_abs_score))
}
report(inf4, fit4, "Table 3 column (4)")
report(inf5, fit5, "Table 3 column (5)")
report(inf5_stationary, fit5_stationary,
       "Table 3 column (5), conditional stationarity")
