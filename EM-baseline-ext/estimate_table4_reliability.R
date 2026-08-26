# Estimate the Table 4 Set 4 model with reliability-dependent symmetric error.

library(here)
library(dplyr)
library(ggplot2)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

results_dir <- here::here("EM-baseline-ext", "output", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260826L)

# Reproduce the Table 4 analysis sample and time-varying Set 4 design.
sector_source <- readRDS(here::here("data", "raw", "QLFSmerged_mapped.rds"))
df_qlfs <- attach_transition_informal_sector(df_qlfs, sector_source)
rm(sector_source)
keep <- complete.cases(df_qlfs[, c("y1", "y2", "y3", "weight",
  "age1", "age2", "age3", "educ1", "educ2", "educ3")]) &
  df_qlfs$weight > 0
df <- as.data.frame(df_qlfs[keep, , drop = FALSE])
df$y1 <- as.integer(df$y1); df$y2 <- as.integer(df$y2)
df$y3 <- as.integer(df$y3); df$weight <- as.numeric(df$weight)
for (nm in c("contracttype1", "contracttype2"))
  df[[nm]] <- ifelse(is.na(df[[nm]]), 0L, as.integer(df[[nm]]))
X4 <- prepare_covariate_matrix(df, covariate_set = 4L)$X_transition
rm(df_qlfs); gc(verbose = FALSE)

# Reproduce the wave-specific reliability variables in Table 3, columns 5--6.
inc <- compute_inconsistency_extent(df)
table2_person_wave_keys <- function(panel) {
  raw <- readRDS(here::here("data", "raw", paste0("df_qlfs_", panel, ".rds")))
  num <- function(x) as.numeric(unclass(x))
  weight <- (num(raw$weight1) * num(raw$weight2) * num(raw$weight3))^(1 / 3)
  use <- num(raw$age1) > 17 & num(raw$age1) < 56 &
    complete.cases(raw[paste0("employed", 1:3)]) &
    is.finite(weight) & weight > 0
  raw <- raw[use, c("hhnr", "pnr", paste0("period", 1:3)), drop = FALSE]
  unique(unlist(lapply(1:3, function(tt)
    paste(raw$hhnr, raw$pnr, raw[[paste0("period", tt)]], sep = "|")),
    use.names = FALSE))
}
panel_B_keys <- table2_person_wave_keys("B")
panel_C_keys <- table2_person_wave_keys("C")
panel_A_keys <- sapply(1:3, function(tt)
  paste(df$hhnr, df$pnr, df[[paste0("period", tt)]], sep = "|"))
in_B <- matrix(panel_A_keys %in% panel_B_keys, nrow(df), 3L)
in_C <- matrix(panel_A_keys %in% panel_C_keys, nrow(df), 3L)
B_not_C <- 1L * (in_B & !in_C)
A_not_B <- 1L * !in_B
stopifnot(sum(B_not_C * A_not_B) == 0L)
Z <- lapply(1:3, function(tt) cbind(
  error_intercept = 1,
  age_inconsistency = inc[[paste0("Y_age_", tt)]],
  education_inconsistency = inc[[paste0("Y_edu_", tt)]],
  large_age_inconsistency = as.integer(inc[[paste0("extent_age_", tt)]] >= 2),
  large_education_inconsistency =
    as.integer(inc[[paste0("extent_edu_", tt)]] >= 2),
  panel_B_not_C = B_not_C[, tt], panel_A_not_B = A_not_B[, tt]))
rm(inc, panel_B_keys, panel_C_keys, panel_A_keys, in_B, in_C,
   B_not_C, A_not_B)
gc(verbose = FALSE)

cells <- collapse_covariate_reliability(df, X4, Z)
cat(sprintf("Combined Table 4 sample: N=%s; %s exact likelihood cells; K=30\n",
  format(cells$n_original, big.mark = ","),
  format(nrow(cells$df), big.mark = ",")))

# Nested transition and error-equation estimates provide deterministic starts.
set4 <- readRDS(file.path(results_dir, "fit_cov_s4_sym_free.rds"))
table3 <- readRDS(file.path(results_dir, "fit_table6_inconsistency_audit.rds"))
delta_names <- colnames(Z[[1L]])
start_combined <- list(beta0 = set4$params$beta0, beta1 = set4$params$beta1,
  alpha = set4$params$alpha,
  delta = setNames(table3$matching_free$params$delta, delta_names))
start_constant <- start_combined
start_constant$delta[] <- 0
start_constant$delta["error_intercept"] <- qlogis(2 * set4$params$pi)
eta_combined <- pack_covariate_reliability(start_combined, X4)
eta_constant <- pack_covariate_reliability(start_constant, X4)
starts <- c(list(nested_constant_error = eta_constant,
                 combined_warm_start = eta_combined),
  setNames(lapply(1:3, function(i) eta_combined +
    rnorm(length(eta_combined), 0, c(rep(.08, 23L), rep(.15, 7L)))),
    paste0("perturbed_", 1:3)))
requested_starts <- as.integer(Sys.getenv("TABLE4_RELIABILITY_STARTS", "5"))
if (!is.finite(requested_starts) || requested_starts < 1L || requested_starts > 5L)
  stop("TABLE4_RELIABILITY_STARTS must be between 1 and 5")
starts <- starts[seq_len(requested_starts)]

fits <- vector("list", length(starts)); names(fits) <- names(starts)
elapsed <- numeric(length(starts))
resume_checkpoint <- identical(Sys.getenv("TABLE4_RELIABILITY_RESUME"), "1")
for (i in seq_along(starts)) {
  cat(sprintf("\nStart %d/%d: %s\n", i, length(starts), names(starts)[i]))
  checkpoint <- file.path(results_dir,
    paste0("fit_cov_s4_reliability_start_", names(starts)[i], ".rds"))
  if (resume_checkpoint && file.exists(checkpoint)) {
    fits[[i]] <- readRDS(checkpoint)
    cat(sprintf("  resumed ll=%.3f code=%d evaluations=%d max|score|=%.3e\n",
      fits[[i]]$loglik, fits[[i]]$optimizer_code, fits[[i]]$iterations,
      fits[[i]]$max_abs_score))
    if (identical(Sys.getenv("TABLE4_RELIABILITY_POLISH"), "1")) {
      preliminary_path <- file.path(results_dir,
        "analytical_se_cov_s4_reliability_preliminary.rds")
      if (!file.exists(preliminary_path)) stop("Missing preliminary inference for polishing")
      preliminary <- readRDS(preliminary_path)
      correction <- preliminary$diagnostics$newton_correction
      fractions <- c(1, .5, .25, .125, 0)
      trial_eta <- lapply(fractions, function(frac) fits[[i]]$eta + frac * correction)
      trial_ll <- vapply(trial_eta, function(z) sum(cells$df$weight *
        .covrel_components(z, cells)$loglik_i), numeric(1L))
      best_trial <- which.max(trial_ll)
      cat(sprintf("  Newton line search: fraction=%.3f ll=%.3f (gain %.3f)\n",
        fractions[best_trial], trial_ll[best_trial],
        trial_ll[best_trial] - fits[[i]]$loglik))
      timing <- system.time(fits[[i]] <- fit_covariate_reliability(
        cells, trial_eta[[best_trial]], maxit = 1500L, reltol = 1e-11))
      elapsed[i] <- unname(timing["elapsed"])
      saveRDS(fits[[i]], checkpoint)
      cat(sprintf("  polished ll=%.3f code=%d evaluations=%d max|score|=%.3e elapsed=%.1fs\n",
        fits[[i]]$loglik, fits[[i]]$optimizer_code, fits[[i]]$iterations,
        fits[[i]]$max_abs_score, elapsed[i]))
    }
    next
  }
  timing <- system.time(fits[[i]] <- fit_covariate_reliability(
    cells, starts[[i]], maxit = 1500L, reltol = 1e-9))
  elapsed[i] <- unname(timing["elapsed"])
  saveRDS(fits[[i]], checkpoint)
  cat(sprintf("  ll=%.3f code=%d evaluations=%d max|score|=%.3e elapsed=%.1fs\n",
    fits[[i]]$loglik, fits[[i]]$optimizer_code, fits[[i]]$iterations,
    fits[[i]]$max_abs_score, elapsed[i]))
}
best_index <- which.max(vapply(fits, `[[`, numeric(1L), "loglik"))
fit <- fits[[best_index]]
candidate_table <- data.frame(start = names(starts),
  loglik = vapply(fits, `[[`, numeric(1L), "loglik"),
  optimizer_code = vapply(fits, `[[`, integer(1L), "optimizer_code"),
  evaluations = vapply(fits, `[[`, numeric(1L), "iterations"),
  max_abs_score = vapply(fits, `[[`, numeric(1L), "max_abs_score"),
  elapsed_seconds = elapsed)

if (fit$loglik < set4$loglik - 1e-4)
  stop("Combined model failed to dominate the nested constant-error Set 4 model")
if (fit$optimizer_code != 0L)
  stop("Best combined model did not pass the optimizer convergence check")

cat("\nCalculating analytical sandwich/delta inference...\n")
inference_time <- system.time(
  inference <- analytical_se_covariate_reliability(cells, fit))
diag <- inference$diagnostics
if (diag$information_rank != diag$K)
  stop("Combined model information matrix is rank deficient")
saveRDS(inference,
  file.path(results_dir, "analytical_se_cov_s4_reliability_preliminary.rds"))
cat(sprintf("Preliminary information condition: %.3g; max Newton correction: %.3e\n",
  diag$information_condition, diag$max_abs_newton_correction))
if (diag$max_abs_newton_correction > 0.01)
  stop("Combined model requires a material post-optimization Newton correction")

fit$analytical_inference <- inference
fit$candidate_table <- candidate_table
fit$estimation_elapsed_seconds <- sum(elapsed)
fit$inference_elapsed_seconds <- unname(inference_time["elapsed"])
fit$n_obs <- cells$n_original
fit$n_cells <- nrow(cells$df)
saveRDS(fit, file.path(results_dir, "fit_cov_s4_reliability_free.rds"))
saveRDS(inference,
  file.path(results_dir, "analytical_se_cov_s4_reliability_free.rds"))
write.csv(candidate_table,
  file.path(results_dir, "cov_s4_reliability_multistart.csv"), row.names = FALSE)
write.csv(inference$summary,
  file.path(results_dir, "cov_s4_reliability_summary.csv"), row.names = FALSE)

reported <- inference$summary[inference$summary$quantity %in% c(
  "mean_entry_rate", "mean_exit_rate", "mean_misclassification_rate",
  "mean_employment_rate", "initial_employment", "pi_base", "pi_age_mild",
  "pi_age_severe", "pi_education_mild", "pi_education_severe",
  "pi_B_not_C", "pi_A_not_B"), ]
reported$estimate_percent <- 100 * reported$estimate
reported$std_error_percent <- 100 * reported$std_error
cat("\nCombined Set 4 / reliability-dependent-error results (percent):\n")
print(reported[, c("quantity", "estimate_percent", "std_error_percent")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nLog-likelihood: %.3f (nested Set 4: %.3f; gain: %.3f)\n",
  fit$loglik, set4$loglik, fit$loglik - set4$loglik))
cat(sprintf("Information rank: %d/%d; condition number: %.3g; inference: %.1fs\n",
  diag$information_rank, diag$K, diag$information_condition,
  fit$inference_elapsed_seconds))
cat(sprintf("Max normalized score: %.3e; max Newton correction: %.3e\n",
  diag$max_abs_score, diag$max_abs_newton_correction))
