# Optional compiled evaluator for the preferred monthly, joint-timegap model.
# The R evaluator remains the reference and the default for existing callers.
load_four_wave_monthly_kernel <- function() {
  if (exists("four_wave_monthly_kernel", mode = "function")) return(invisible())
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required")
  cache <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1/rcpp-cache"
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp("EM-tenure/src/four_wave_monthly.cpp", cacheDir = cache,
    showOutput = FALSE, verbose = FALSE)
  invisible()
}

prepare_four_wave_kernel_data <- function(df) {
  validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  tenure <- as.matrix(df[paste0("tenure", 1:4)])
  observed <- as.matrix(df[paste0("y", 1:4)]) == 1L
  if (any(observed & !matrix(.tenure_month_index(tenure)$valid,
      nrow(tenure), ncol(tenure))))
    stop("Observed tenure is not on the non-negative monthly grid")
  month <- round(12 * tenure)
  month[!is.finite(month)] <- -1
  y <- as.matrix(df[paste0("y", 1:4)])
  category <- as.matrix(df[paste0("timegap_cat", 1:4)])
  category[is.na(category)] <- 0
  interview <- as.matrix(df[paste0("interview_month", 1:4)])
  if (any(!interview %in% 1:12)) stop("Invalid interview month")
  storage.mode(month) <- storage.mode(y) <- storage.mode(category) <-
    storage.mode(interview) <- "integer"
  list(y = y, month = month, category = category, interview = interview,
    weight = df$weight, max_month = max(month) + 9L)
}

evaluate_four_wave_monthly_fast <- function(data, params, posterior = TRUE,
    threads = 4L, exact_risk = TRUE,
    timegap_clock = params$timegap_clock_model %||% "continuous_joint") {
  timegap_clock <- match.arg(timegap_clock,c("continuous_joint","legacy_pairwise"))
  load_four_wave_monthly_kernel()
  if (!identical(params$tenure_measurement_model, "monthly") ||
      !identical(params$timegap_contamination_model, "joint_marginal") ||
      length(params$lambda_g) != 5L || length(params$lambda_d) != 5L ||
      (params$tenure_report_persistence %||% 0) != 0 ||
      (params$beta_g %||% 0) != 0 || (params$beta_d %||% 0) != 0)
    stop("Compiled evaluator supports the preferred monthly piecewise model only")
  components <- lapply(c("reliable", "unreliable"), function(cl)
    .duration_reliability_component_params(params, cl))
  if (any(vapply(components, function(p)
      p$eps <= 0 || p$eps >= 1 || p$eps_d < 0 || p$eps_d > 1, logical(1))))
    stop("Compiled duration contamination probabilities must be interior")
  m <- 0:data$max_month
  probs <- .validate_start_month_probs(params$tenure_start_month_probs)
  norm <- .log_january_duration_normalizers(params$lambda_g)
  pair <- as.matrix(expand.grid(previous = 1:7, current = 1:7))
  cat3 <- as.matrix(expand.grid(rep(list(1:7), 3)))
  cat4 <- as.matrix(expand.grid(rep(list(1:7), 4)))
  repaired <- identical(timegap_clock,"continuous_joint")
  clock_tables <- if (repaired) .four_wave_timegap_clock_tables(
    params$lambda_d,components) else rep(list(matrix(0,0,2)),4)
  tables <- list(
    timegap_clock_joint=repaired,
    timegap_clock1=clock_tables[[1]],timegap_clock2=clock_tables[[2]],
    timegap_clock3=clock_tables[[3]],timegap_clock4=clock_tables[[4]],
    raw = .log_duration_month_mass(m, params$lambda_g), normalizers = norm,
    calendar = vapply(1:12, function(im)
      .log_calendar_duration_month_mass(m, rep(im, length(m)),
        params$lambda_g, start_month_probs = probs, log_residue_masses = norm),
      numeric(length(m))),
    within = vapply(1:12, function(im)
      .log_within_interval_start_month_mass(0:2, rep(im, 3), probs), numeric(3)),
    exit = .duration_transition_probability(m / 12, params$lambda_g),
    entry = if (exact_risk) .piecewise_category_transition_risk_exact(1:7,
      params$lambda_d) else .duration_category_transition_probability(1:7, params$lambda_d),
    exit_missing = if (exact_risk) .piecewise_mean_transition_risk_exact(params$lambda_g)
      else .duration_marginal_transition_probability(params$lambda_g),
    entry_missing = if (exact_risk) .piecewise_mean_transition_risk_exact(params$lambda_d)
      else .duration_marginal_transition_probability(params$lambda_d),
    timegap_marginal = log_emission_interval_d(1:7, params$lambda_d),
    timegap_pair = if (repaired) matrix(0,0,2) else vapply(components, function(p)
      log_emission_transition_d_contaminated(pair[, 2], pair[, 1],
        params$lambda_d, eps_d = 1 - (1 - p$eps_d)^2), numeric(49)),
    timegap_joint3 = if (repaired) matrix(0,0,2) else vapply(components, function(p)
      log_emission_timegap_spell_joint(cat3, params$lambda_d, eps_d = p$eps_d),
      numeric(343)),
    timegap_joint4 = if (repaired) matrix(0,0,2) else vapply(components, function(p)
      log_emission_timegap_spell_joint(cat4, params$lambda_d, eps_d = p$eps_d),
      numeric(2401)))
  p <- c(alpha = params$alpha, pi = params$pi,
    job_change = params$job_change_prob, heaping = params$tenure_heaping_prob,
    year = params$tenure_year_revision_prob,
    clean_revision = params$tenure_clean_anchor_revision_prob,
    exact = params$tenure_exact_anchor_retention_prob,
    local = params$tenure_local_revision_prob,
    eps_reliable = components[[1]]$eps, eps_unreliable = components[[2]]$eps)
  threads <- max(1L, min(8L, as.integer(threads)))
  four_wave_monthly_kernel(data, tables, p, posterior, threads)
}

four_wave_fast_source_fingerprint <- function() tools::md5sum(c(
  "EM-tenure/src/four_wave_monthly.cpp", "EM-tenure/R/four_wave_fast.R",
  "EM-tenure/R/piecewise_risk_exact.R", "EM-tenure/R/four_wave_eps.R",
  "EM-tenure/R/emissions.R", "EM-tenure/R/tenure_monthly_eps.R",
  "EM-tenure/R/job_change_eps.R", "EM-tenure/R/timegap_clock_joint.R"))
