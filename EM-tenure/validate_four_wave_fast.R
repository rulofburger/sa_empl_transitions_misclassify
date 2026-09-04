source("EM-tenure/R/source_all.R")
source("EM-tenure/R/four_wave_fast.R")
outdir <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
df <- readRDS(file.path(outdir, "prepared_cells_latest.rds"))$df4
fit_path <- Sys.getenv("FOUR_WAVE_FAST_VALIDATE_FIT",
  file.path(outdir, "fit_reporting2_latest.rds"))
fit <- readRDS(fit_path)
clock_model <- match.arg(Sys.getenv("FOUR_WAVE_TIMEGAP_CLOCK","continuous_joint"),
  c("continuous_joint","legacy_pairwise"))
fit$params$timegap_clock_model <- clock_model
set.seed(3903)
code <- as.vector(as.matrix(df[paste0("y", 1:4)]) %*% 2^(0:3))
rows <- unique(c(which(!duplicated(code)), sample(nrow(df), 250)))
small <- df[rows, ]
data <- prepare_four_wave_kernel_data(small)
z <- .pack_four_wave_preferred(fit$params)
perturb <- .piecewise_calendar_revision_monthly_unpack(
  z + .12 * sin(seq_along(z)), "joint_marginal")
perturb$timegap_clock_model <- clock_model
simple <- fit$params
simple$tenure_start_month_probs <- rep(1/12, 12)
simple$tenure_heaping_prob <- simple$tenure_year_revision_prob <-
  simple$tenure_clean_anchor_revision_prob <-
  simple$tenure_exact_anchor_retention_prob <-
  simple$tenure_local_revision_prob <- simple$job_change_prob <- 0
simple$tenure_reliability_shift <- simple$timegap_reliability_shift <- 0
exact_status <- fit$params
exact_status$pi <- 0
parameters <- list(current = fit$params, perturbed = perturb,
  simple = simple, exact_status = exact_status)
results <- lapply(names(parameters), function(label) {
  p <- parameters[[label]]
  a <- e_step_eps_4w(small, p, check_df = FALSE, exact_risk = TRUE)
  b <- evaluate_four_wave_monthly_fast(data, p, threads = 1L)
  parallel <- evaluate_four_wave_monthly_fast(data, p, threads = 8L)
  result <- data.frame(case = label,
    max_row_difference = max(abs(a$row_loglik - b$row_loglik)),
    max_gamma_difference = max(abs(a$gamma - b$gamma)),
    max_job_difference = max(abs(a$job_change_posterior$expected_changes -
      b$job_change_posterior$expected_changes)),
    max_opportunities_difference = max(abs(a$job_change_posterior$opportunities -
      b$job_change_posterior$opportunities)),
    max_parallel_difference = max(abs(b$gamma - parallel$gamma)))
  stopifnot(all(is.finite(unlist(result[-1]))), all(unlist(result[-1]) < 1e-10))
  result
})
checks <- do.call(rbind, results)
data <- prepare_four_wave_kernel_data(df)
timing <- system.time(full <- evaluate_four_wave_monthly_fast(data, fit$params,
  threads = 8L))
# Validate the computational kernel against the saved legacy R result too.
legacy <- evaluate_four_wave_monthly_fast(data, fit$params,
  threads = 8L, exact_risk = FALSE)
exact_reference <- clock_model=="continuous_joint" ||
  identical(Sys.getenv("FOUR_WAVE_FAST_REFERENCE_FULL", "false"),
    "true") || identical(fit$integration_method, "exact_piecewise")
reference <- if (exact_reference)
  e_step_eps_4w(df, fit$params, check_df = FALSE, exact_risk = TRUE) else fit
target <- if (exact_reference) full else legacy
full_check <- c(loglik_difference = target$loglik - reference$loglik,
  gamma_max_difference = max(abs(target$gamma - reference$gamma)),
  job_max_difference = max(abs(target$job_change_posterior$expected_changes -
    reference$job_change_posterior$expected_changes)))
integration_difference <- c(loglik = full$loglik - legacy$loglik,
  max_posterior = max(abs(full$gamma - legacy$gamma)))
stopifnot(abs(full_check[1]) < 1e-6, all(abs(full_check[-1]) < 1e-10))
print(checks)
print(full_check)
print(integration_difference)
print(timing)
saveRDS(list(checks = checks, full_check = full_check, timing = timing,
  timegap_clock=clock_model,
  reference_fit = fit_path, exact_reference = exact_reference,
  integration_difference = integration_difference,
  source_md5 = four_wave_fast_source_fingerprint(), session_info = sessionInfo()),
  file.path(outdir, "compiled_equivalence_latest.rds"))
write.csv(checks, file.path(outdir, "compiled_equivalence_latest.csv"), row.names=FALSE)
