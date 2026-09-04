# Four-wave observed likelihood for the corrected discrete-month duration
# model.  This module deliberately leaves the established three-wave E-step
# unchanged.  It supplies the AR(1) baseline that must be validated before a
# second-order transition law is introduced.

latent_histories_eps_4w <- function() {
  out <- as.matrix(expand.grid(rep(list(0:1), 4L)))
  storage.mode(out) <- "integer"
  colnames(out) <- paste0("h", 1:4)
  out
}

.maximal_state_spells <- function(history, state) {
  out <- list()
  start <- which(history == state &
    c(TRUE, history[-length(history)] != state))
  if (!length(start)) return(out)
  end <- which(history == state &
    c(history[-1L] != state, TRUE))
  Map(seq.int, start, end)
}

validate_df_eps_4w <- function(df, allow_zero_tenure = TRUE) {
  y_cols <- paste0("y", 1:4)
  tenure_cols <- paste0("tenure", 1:4)
  timegap_cols <- paste0("timegap_cat", 1:4)
  required <- c(y_cols, tenure_cols, timegap_cols, "weight")
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("four-wave eps data require: ",
    paste(missing, collapse = ", "))
  y <- as.matrix(df[y_cols])
  tenure <- as.matrix(df[tenure_cols])
  timegap <- as.matrix(df[timegap_cols])
  if (anyNA(y) || any(!y %in% 0:1))
    stop("y1-y4 must be observed and binary")
  invalid_tenure <- y == 1L & (!is.finite(tenure) |
    if (allow_zero_tenure) tenure < 0 else tenure <= 0)
  if (any(invalid_tenure))
    stop(sum(invalid_tenure), " employed-wave tenure reports are invalid")
  invalid_timegap <- y == 0L &
    (is.na(timegap) | !timegap %in% seq_len(.N_TIMEGAP_CATS))
  if (any(invalid_timegap))
    stop(sum(invalid_timegap), " nonemployment-wave timegap reports are invalid")
  if (any(!is.finite(df$weight)) || any(df$weight <= 0))
    stop("weights must be finite and positive")
  invisible(NULL)
}

prepare_eps_estimation_data_4w <- function(df, drop_incomplete = TRUE) {
  if (requireNamespace("haven", quietly = TRUE)) {
    labelled <- vapply(df, haven::is.labelled, logical(1))
    df[labelled] <- lapply(df[labelled], as.numeric)
  }
  y_cols <- paste0("y", 1:4)
  tenure_cols <- paste0("tenure", 1:4)
  raw_timegap_cols <- paste0("timegap", 1:4)
  timegap_cols <- paste0("timegap_cat", 1:4)
  required <- c(y_cols, tenure_cols, raw_timegap_cols,
    paste0("age", 1:4), paste0("neverworked", 1:4))
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("four-wave preparation requires: ",
    paste(missing, collapse = ", "))
  neverworked_category <- function(age_years) as.integer(cut(
    pmax((age_years - 16) * 12, 0),
    breaks = c(0, 3, 6, 9, 12, 36, 60, Inf), right = FALSE,
    include.lowest = TRUE))
  for (wave in 1:4) {
    y <- as.integer(df[[y_cols[wave]]])
    tenure <- as.numeric(df[[tenure_cols[wave]]])
    tenure[tenure < 0] <- NA_real_
    df[[tenure_cols[wave]]] <- tenure / 12
    raw_timegap <- as.numeric(df[[raw_timegap_cols[wave]]])
    category <- ifelse(raw_timegap %in% 1:7,
      as.integer(raw_timegap), NA_integer_)
    neverworked <- as.numeric(df[[paste0("neverworked", wave)]]) == 1
    replace <- !is.na(y) & y == 0L & !is.na(neverworked) & neverworked
    category[replace] <- neverworked_category(
      as.numeric(df[[paste0("age", wave)]][replace]))
    df[[timegap_cols[wave]]] <- as.integer(category)
    df[[tenure_cols[wave]]][y == 0L] <- NA_real_
    df[[timegap_cols[wave]]][y == 1L] <- NA_integer_
  }
  if (!"weight" %in% names(df)) {
    weight <- as.matrix(df[paste0("weight", 1:4)])
    df$weight <- apply(weight, 1L, prod)^(1/4)
  }
  complete <- rep(TRUE, nrow(df))
  for (wave in 1:4) {
    y <- df[[y_cols[wave]]]
    complete <- complete & is.finite(y) & y %in% 0:1 &
      (y == 0L | is.finite(df[[tenure_cols[wave]]])) &
      (y == 1L | !is.na(df[[timegap_cols[wave]]]))
  }
  complete <- complete & is.finite(df$weight) & df$weight > 0
  attr(df, "complete_duration_reports") <- complete
  if (drop_incomplete) df <- df[complete, , drop = FALSE]
  validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  df
}

add_nominal_interview_months_4w <- function(df) {
  period_cols <- paste0("period", 1:4)
  missing <- setdiff(period_cols, names(df))
  if (length(missing)) stop("nominal interview months require period1-period4")
  for (wave in 1:4) {
    period <- as.integer(df[[period_cols[wave]]])
    if (any(!is.finite(period) | period < 1L))
      stop("period columns must contain positive integers")
    df[[paste0("interview_month", wave)]] <-
      3L * ((period - 1L) %% 4L + 1L)
  }
  df
}

collapse_eps_cells_4w <- function(df, normalize_weights = TRUE,
    extra_cols = character()) {
  validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  missing <- setdiff(extra_cols, names(df))
  if (length(missing)) stop("unknown extra collapse columns: ",
    paste(missing, collapse = ", "))
  cols <- c(paste0("y", 1:4), paste0("tenure", 1:4),
    paste0("timegap_cat", 1:4), extra_cols)
  key <- do.call(paste, c(df[cols], sep = "\r"))
  first <- !duplicated(key)
  group <- match(key, key[first])
  out <- df[first, c(cols, "weight"), drop = FALSE]
  out$weight <- as.vector(rowsum(df$weight, group, reorder = FALSE))
  if (normalize_weights) out$weight <- nrow(df) * out$weight / sum(out$weight)
  attr(out, "n_original") <- nrow(df)
  attr(out, "n_cells") <- nrow(out)
  out
}

.log_duration_history_prior_eps_4w <- function(hmat, alpha, s, g, d,
    lambda_g, beta_g, lambda_d, beta_d, exact_risk = FALSE) {
  n <- length(s[[1L]])
  out <- matrix(0, n, nrow(hmat))
  if (exact_risk && (length(lambda_g) != 5L || length(lambda_d) != 5L ||
      beta_g != 0 || beta_d != 0)) stop("Exact risks require piecewise hazards")
  p_exit_missing <- if (exact_risk) .piecewise_mean_transition_risk_exact(lambda_g)
    else .duration_marginal_transition_probability(lambda_g, beta_g)
  p_entry_missing <- if (exact_risk) .piecewise_mean_transition_risk_exact(lambda_d)
    else .duration_marginal_transition_probability(lambda_d, beta_d)
  p_exit <- p_entry <- vector("list", 3L)
  for (wave in 1:3) {
    p_exit[[wave]] <- .duration_transition_probability(
      g[[wave]], lambda_g, beta_g)
    missing_exit <- s[[wave]] == 0L | !is.finite(p_exit[[wave]])
    p_exit[[wave]][missing_exit] <- p_exit_missing
    p_entry[[wave]] <- if (exact_risk)
      .piecewise_category_transition_risk_exact(d[[wave]], lambda_d) else
      .duration_category_transition_probability(d[[wave]], lambda_d, beta_d)
    missing_entry <- s[[wave]] == 1L | !is.finite(p_entry[[wave]])
    p_entry[[wave]][missing_entry] <- p_entry_missing
    p_exit[[wave]] <- pmin(pmax(p_exit[[wave]], 1e-12), 1 - 1e-12)
    p_entry[[wave]] <- pmin(pmax(p_entry[[wave]], 1e-12), 1 - 1e-12)
  }
  for (history in seq_len(nrow(hmat))) {
    h <- hmat[history, ]
    out[, history] <- if (h[1L] == 1L) log(alpha) else log1p(-alpha)
    for (wave in 1:3) {
      p_change <- if (h[wave] == 1L) p_exit[[wave]] else p_entry[[wave]]
      out[, history] <- out[, history] +
        if (h[wave + 1L] == h[wave]) log1p(-p_change) else log(p_change)
    }
  }
  out
}

# Per-report contamination patterns for a fully observed latent
# nonemployment spell of arbitrary length.  At length three this is exactly
# the established triplet likelihood.
log_emission_timegap_spell_joint <- function(categories, lambda_d,
    beta_d = 0, eps_d) {
  if (!is.matrix(categories)) categories <- as.matrix(categories)
  if (!ncol(categories)) return(numeric(nrow(categories)))
  if (!is.finite(eps_d) || eps_d <= 0 || eps_d >= 1)
    stop("eps_d must be in (0,1) for a joint timegap spell")
  patterns <- as.matrix(expand.grid(rep(list(0:1), ncol(categories))))
  marginal <- lapply(seq_len(ncol(categories)), function(wave)
    log_emission_interval_d(categories[, wave], lambda_d, beta_d))
  terms <- matrix(-Inf, nrow(categories), nrow(patterns))
  for (pattern in seq_len(nrow(patterns))) {
    contaminated <- which(patterns[pattern, ] == 1L)
    clean <- which(patterns[pattern, ] == 0L)
    value <- rep(length(contaminated) * log(eps_d) +
      length(clean) * log1p(-eps_d), nrow(categories))
    if (length(contaminated)) for (wave in contaminated)
      value <- value + marginal[[wave]]
    if (length(clean)) {
      value <- value + marginal[[clean[1L]]]
      if (length(clean) > 1L) for (position in 2:length(clean)) {
        current <- clean[position]
        previous <- clean[position - 1L]
        value <- value + log_emission_transition_d(
          categories[, current], categories[, previous], lambda_d,
          beta_d, delta = (current - previous) * .QUARTER_YEARS)
      }
    }
    terms[, pattern] <- value
  }
  .row_logsumexp(terms)
}

e_step_eps_4w <- function(df, params, check_df = TRUE,
    suff_stats = FALSE, exact_risk = FALSE,
    timegap_clock = params$timegap_clock_model %||% "continuous_joint") {
  timegap_clock <- match.arg(timegap_clock,c("continuous_joint","legacy_pairwise"))
  if (suff_stats)
    stop("the four-wave direct-likelihood estimator does not use EM statistics")
  shifts <- .duration_reliability_shifts(params)
  if (any(shifts > 0)) {
    reliable_params <- .duration_reliability_component_params(params,
      "reliable")
    unreliable_params <- .duration_reliability_component_params(params,
      "unreliable")
    reliable <- e_step_eps_4w(df, reliable_params, check_df = check_df,
      exact_risk = exact_risk, timegap_clock=timegap_clock)
    unreliable <- e_step_eps_4w(df, unreliable_params, check_df = FALSE,
      exact_risk = exact_risk, timegap_clock=timegap_clock)
    row_max <- pmax(reliable$row_loglik, unreliable$row_loglik)
    row_loglik <- row_max + log(.5 * exp(reliable$row_loglik - row_max) +
      .5 * exp(unreliable$row_loglik - row_max))
    posterior_unreliable <- .5 * exp(unreliable$row_loglik - row_loglik)
    posterior_reliable <- 1 - posterior_unreliable
    gamma <- reliable$gamma * posterior_reliable +
      unreliable$gamma * posterior_unreliable
    job_change_posterior <- list(
      expected_changes = reliable$job_change_posterior$expected_changes *
        posterior_reliable +
        unreliable$job_change_posterior$expected_changes *
        posterior_unreliable,
      opportunities = reliable$job_change_posterior$opportunities *
        posterior_reliable +
        unreliable$job_change_posterior$opportunities *
        posterior_unreliable)
    return(list(gamma = gamma, loglik = sum(df$weight * row_loglik),
      row_loglik = row_loglik, suff = NULL,
      job_change_posterior = job_change_posterior,
      duration_reliability_posterior = posterior_unreliable,
      duration_reliability_component_probabilities = data.frame(
        class = c("reliable", "unreliable"),
        tenure_contamination = c(reliable_params$eps,
          unreliable_params$eps),
        timegap_contamination = c(reliable_params$eps_d,
          unreliable_params$eps_d))))
  }
  if (check_df) validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  if (!identical(params$tenure_measurement_model, "monthly"))
    stop("the four-wave corrected estimator requires monthly tenure")
  alpha <- params$alpha
  pi_par <- params$pi
  eps <- params$eps
  eps_d <- params$eps_d %||% 0
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d
  beta_g <- params$beta_g %||% 0
  beta_d <- params$beta_d %||% 0
  job_change <- params$job_change_prob %||% 0
  heaping <- params$tenure_heaping_prob %||% 0
  year_revision <- params$tenure_year_revision_prob %||% 0
  clean_revision <- params$tenure_clean_anchor_revision_prob %||% 0
  exact_retention <- params$tenure_exact_anchor_retention_prob %||% 0
  local_revision <- params$tenure_local_revision_prob %||% 0
  start_month_probs <- .validate_start_month_probs(
    params$tenure_start_month_probs)
  timegap_model <- params$timegap_contamination_model %||% "joint_marginal"
  local_decay <- params$timegap_local_decay %||% 1
  probabilities <- c(alpha, pi_par, eps, eps_d, job_change, heaping,
    year_revision, clean_revision, exact_retention, local_revision)
  if (any(!is.finite(probabilities)) || alpha <= 0 || alpha >= 1 ||
      pi_par < 0 || pi_par >= 1 || eps <= 0 || eps >= 1 || eps_d < 0 ||
      eps_d >= 1 || any(probabilities[-1L] < 0))
    stop("invalid four-wave probability parameter")
  if (any(!is.finite(c(lambda_g, lambda_d))) ||
      any(c(lambda_g, lambda_d) <= 0)) stop("duration hazards must be positive")
  if (!timegap_model %in% c("marginal", "local", "joint_marginal"))
    stop("unknown timegap contamination model")
  interview_cols <- paste0("interview_month", 1:4)
  if (length(setdiff(interview_cols, names(df))) ||
      any(!as.matrix(df[interview_cols]) %in% 1:12))
    stop("four-wave calendar likelihood requires interview_month1-4")

  hmat <- latent_histories_eps_4w()
  n <- nrow(df)
  histories <- nrow(hmat)
  s <- lapply(1:4, function(wave) df[[paste0("y", wave)]])
  g <- lapply(1:4, function(wave) df[[paste0("tenure", wave)]])
  d <- lapply(1:4, function(wave) df[[paste0("timegap_cat", wave)]])
  s_mat <- as.matrix(df[paste0("y", 1:4)]) == 1L
  g_mat <- as.matrix(df[paste0("tenure", 1:4)])
  interview <- as.matrix(df[interview_cols])
  h_by_wave <- t(hmat)
  mismatch <- matrix(rowSums(s_mat), n, histories) +
    matrix(colSums(h_by_wave), n, histories, byrow = TRUE) -
    2L * (s_mat %*% h_by_wave)
  log_prior <- .log_duration_history_prior_eps_4w(hmat, alpha, s, g, d,
    lambda_g, beta_g, lambda_d, beta_d, exact_risk = exact_risk)
  if (pi_par < .Machine$double.eps) {
    status_loglik <- matrix(-Inf, n, histories)
    status_loglik[mismatch == 0L] <- 0
  } else {
    pi_b <- .bound01(pi_par)
    status_loglik <- mismatch * log(pi_b) +
      (4L - mismatch) * log1p(-pi_b)
  }
  emission <- matrix(0, n, histories)
  job_num <- matrix(0, n, histories)
  job_den <- matrix(0, n, histories)
  seasonal <- max(abs(start_month_probs - 1/12)) > 1e-12
  pair_eps_d <- if (identical(timegap_model, "joint_marginal"))
    1 - (1 - eps_d)^2 else eps_d
  conditional_model <- if (identical(timegap_model, "local"))
    "local" else "marginal"
  spell_cache <- list()
  segment_cache <- new.env(parent = emptyenv())
  joint_cache <- list()
  singleton_tenure <- matrix(0, n, 4L)
  timegap_marginal <- lapply(1:4, function(wave)
    log_emission_interval_d(d[[wave]], lambda_d, beta_d))
  timegap_conditional <- vector("list", 4L)
  for (wave in 2:4)
    timegap_conditional[[wave]] <- log_emission_transition_d_contaminated(
      d[[wave]], d[[wave - 1L]], lambda_d, beta_d, pair_eps_d,
      conditional_model, local_decay)
  for (wave in 1:4) {
    mask <- s[[wave]] == 1L
    if (!any(mask)) next
    marginal <- if (seasonal)
      .log_calendar_duration_month_mass(
        round(12 * g[[wave]][mask]), interview[mask, wave],
        lambda_g, beta_g, start_month_probs) else
      .log_duration_month_mass(round(12 * g[[wave]][mask]),
        lambda_g, beta_g)
    if (heaping > 0) {
      heaped <- .log_january_duration_month_mass(
        round(12 * g[[wave]][mask]), interview[mask, wave],
        lambda_g, beta_g)
      marginal <- .log_probability_mixture(marginal, heaped, heaping)
    }
    singleton_tenure[mask, wave] <- marginal
  }

  for (history in seq_len(histories)) {
    h <- hmat[history, ]
    for (spell in .maximal_state_spells(h, 1L)) {
      key <- paste(spell, collapse = "_")
      out <- spell_cache[[key]]
      if (is.null(out)) {
        out <- log_emission_spell_g_monthly(
        g_mat[, spell, drop = FALSE],
        s_mat[, spell, drop = FALSE],
        as.integer(spell - spell[1L]), lambda_g, eps, job_change,
        beta_g, initial_model = if (spell[1L] > 1L)
          "within_interval" else "marginal",
        tenure_heaping_prob = heaping,
        tenure_year_revision_prob = year_revision,
        tenure_clean_anchor_revision_prob = clean_revision,
        tenure_exact_anchor_retention_prob = exact_retention,
        tenure_local_revision_prob = local_revision,
        tenure_start_month_probs = start_month_probs,
        interview_month_mat = interview[, spell, drop = FALSE],
        segment_cache = segment_cache, segment_ids = spell)
        spell_cache[[key]] <- out
      }
      emission[, history] <- emission[, history] + out$loglik
      job_num[, history] <- job_num[, history] + out$job_changes
      job_den[, history] <- job_den[, history] +
        out$job_change_opportunities
    }
    for (wave in which(h == 0L))
      emission[, history] <- emission[, history] + singleton_tenure[, wave]

    joint_covered <- matrix(FALSE, n, 4L)
    if (identical(timegap_clock,"continuous_joint")) {
      if (!identical(timegap_model,"joint_marginal"))
        stop("Continuous joint clock requires joint_marginal contamination")
      for (spell in .maximal_state_spells(h,0L)) {
        key <- paste(spell,collapse="_")
        joint <- joint_cache[[key]]
        if (is.null(joint)) {
          joint <- log_emission_timegap_clock_joint(
            as.matrix(df[paste0("timegap_cat",spell)]),lambda_d,beta_d,
            eps_d,observed=!s_mat[,spell,drop=FALSE])
          joint_cache[[key]] <- joint
        }
        emission[,history] <- emission[,history]+joint
        joint_covered[,spell] <- TRUE
      }
    } else if (identical(timegap_model, "joint_marginal") && eps_d > 0) {
      for (spell in .maximal_state_spells(h, 0L)) if (length(spell) >= 3L) {
        key <- paste(spell, collapse = "_")
        joint <- joint_cache[[key]]
        if (is.null(joint)) {
          use <- rowSums(!s_mat[, spell, drop = FALSE]) == length(spell)
          joint <- list(use = use, loglik = numeric(sum(use)))
          if (any(use)) joint$loglik <- log_emission_timegap_spell_joint(
            as.matrix(df[use, paste0("timegap_cat", spell), drop = FALSE]),
            lambda_d, beta_d, eps_d)
          joint_cache[[key]] <- joint
        }
        use <- joint$use
        if (any(use)) {
          emission[use, history] <- emission[use, history] +
            joint$loglik
          joint_covered[use, spell] <- TRUE
        }
      }
    }
    for (wave in 1:4) {
      observed_nonemployment <- s[[wave]] == 0L &
        !joint_covered[, wave]
      if (!any(observed_nonemployment)) next
      continuation <- wave > 1L && h[wave - 1L] == 0L && h[wave] == 0L
      if (continuation) {
        conditional <- observed_nonemployment & s[[wave - 1L]] == 0L &
          !joint_covered[, wave - 1L]
        if (any(conditional))
          emission[conditional, history] <- emission[conditional, history] +
            timegap_conditional[[wave]][conditional]
        marginal <- observed_nonemployment & !conditional
      } else marginal <- observed_nonemployment
      if (any(marginal))
        emission[marginal, history] <- emission[marginal, history] +
          timegap_marginal[[wave]][marginal]
    }
  }

  kernel <- status_loglik + emission + log_prior
  if (pi_par < .Machine$double.eps) {
    match_index <- 1L + s[[1L]] + 2L * s[[2L]] + 4L * s[[3L]] +
      8L * s[[4L]]
    gamma <- matrix(0, n, histories)
    gamma[cbind(seq_len(n), match_index)] <- 1
    row_loglik <- kernel[cbind(seq_len(n), match_index)]
  } else {
    row_loglik <- .row_logsumexp(kernel)
    gamma <- exp(kernel - row_loglik)
    gamma[!is.finite(gamma)] <- 0
  }
  list(gamma = gamma, loglik = sum(df$weight * row_loglik),
    row_loglik = row_loglik, suff = NULL,
    job_change_posterior = list(
      expected_changes = rowSums(gamma * job_num),
      opportunities = rowSums(gamma * job_den)))
}

.pack_four_wave_preferred <- function(params) {
  shifts <- .duration_reliability_shifts(params)
  separate <- !is.null(params$tenure_reliability_shift) ||
    !is.null(params$timegap_reliability_shift)
  common <- params$duration_reliability_shift %||% 0
  .piecewise_calendar_revision_monthly_pack(params,
    heaping_start = params$tenure_heaping_prob,
    revision_start = params$tenure_year_revision_prob,
    q_start = params$job_change_prob,
    clean_anchor_revision_start = params$tenure_clean_anchor_revision_prob,
    start_month_probs_start = params$tenure_start_month_probs,
    exact_anchor_retention_start =
      params$tenure_exact_anchor_retention_prob,
    local_revision_start = params$tenure_local_revision_prob,
    reliability_shift_start = if (!separate) common else NULL,
    tenure_reliability_shift_start = if (separate)
      unname(shifts["tenure"]) else NULL,
    timegap_reliability_shift_start = if (separate)
      unname(shifts["timegap"]) else NULL)
}

.four_wave_parameter_bounds <- function(z_full0) {
  lower_full <- rep(-Inf, length(z_full0))
  upper_full <- rep(Inf, length(z_full0))
  names(lower_full) <- names(upper_full) <- names(z_full0)
  hazards <- c(paste0("log_hg", 1:5), paste0("log_hd", 1:5))
  lower_full[hazards] <- log(1e-4)
  upper_full[hazards] <- log(20)
  probability_logits <- intersect(c("job_change", "tenure_heaping",
    "tenure_year_revision", "tenure_clean_anchor_revision",
    "tenure_exact_anchor_retention", "tenure_local_revision"),
    names(z_full0))
  lower_full[probability_logits] <- -12
  upper_full[probability_logits] <- 12
  shift_names <- intersect(c("duration_reliability_shift",
    "tenure_reliability_shift", "timegap_reliability_shift"),
    names(z_full0))
  lower_full[shift_names] <- 0
  upper_full[shift_names] <- 6
  month_names <- intersect(paste0("start_month_logit", 1:11),
    names(z_full0))
  lower_full[month_names] <- -8
  upper_full[month_names] <- 8
  list(lower = lower_full, upper = upper_full)
}

fit_eps_piecewise_calendar_revision_monthly_4w <- function(df, start,
    maxit = 10L, reltol = 1e-9, pgtol = 1e-7, workers = 1L,
    gradient_step = 1e-4, free_names = NULL,
    gradient_scheme = c("forward", "central"), verbose = 0L) {
  gradient_scheme <- match.arg(gradient_scheme)
  validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  z_full0 <- .pack_four_wave_preferred(start)
  bounds <- .four_wave_parameter_bounds(z_full0)
  lower_full <- bounds$lower
  upper_full <- bounds$upper
  if (is.null(free_names)) free_names <- names(z_full0)
  if (!length(free_names) || anyDuplicated(free_names) ||
      length(setdiff(free_names, names(z_full0))))
    stop("free_names must uniquely identify packed parameters")
  z0 <- z_full0[free_names]
  lower <- lower_full[free_names]
  upper <- upper_full[free_names]
  total_weight <- sum(df$weight)
  expand <- function(z) {
    out <- z_full0
    out[names(z)] <- z
    out
  }
  objective <- function(z) {
    params <- .piecewise_calendar_revision_monthly_unpack(expand(z),
      "joint_marginal")
    value <- tryCatch(e_step_eps_4w(df, params, check_df = FALSE)$loglik,
      error = function(error) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value / total_weight
  }
  cluster <- NULL
  gradient <- NULL
  workers <- max(1L, min(8L, as.integer(workers)))
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster, function(path) {
      setwd(path)
      source("EM-tenure/R/source_all.R")
      NULL
    }, worker_wd)
    df_worker <- df
    z_full_worker <- z_full0
    free_worker <- free_names
    weight_worker <- total_weight
    parallel::clusterExport(cluster, c("df_worker", "z_full_worker",
      "free_worker", "weight_worker"), envir = environment())
    worker_objective <- function(z) {
      full <- z_full_worker
      full[free_worker] <- z
      params <- .piecewise_calendar_revision_monthly_unpack(full,
        "joint_marginal")
      value <- tryCatch(e_step_eps_4w(df_worker, params,
        check_df = FALSE)$loglik, error = function(error) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value / weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step * pmax(1, abs(z))
      plus_step <- pmin(step, upper - z)
      minus_step <- pmin(step, z - lower)
      if (gradient_scheme == "forward") {
        step <- ifelse(plus_step >= step, step, -minus_step)
        points <- lapply(seq_along(z), function(index) {
          point <- z
          point[index] <- point[index] + step[index]
          point
        })
        plus <- unlist(parallel::parLapplyLB(cluster, points,
          worker_objective), use.names = FALSE)
        centre <- objective(z)
        out <- (plus - centre) / step
      } else {
        points <- vector("list", 2L * length(z))
        for (index in seq_along(z)) {
          points[[2L * index - 1L]] <- points[[2L * index]] <- z
          points[[2L * index - 1L]][index] <- z[index] + plus_step[index]
          points[[2L * index]][index] <- z[index] - minus_step[index]
        }
        values <- unlist(parallel::parLapplyLB(cluster, points,
          worker_objective), use.names = FALSE)
        out <- (values[seq(1L, length(values), 2L)] -
          values[seq(2L, length(values), 2L)]) / (plus_step + minus_step)
      }
      if (any(!is.finite(out)) || any(abs(out) > 1e90))
        stop("Invalid finite-difference likelihood in four-wave optimizer")
      pmin(pmax(out, -1e6), 1e6)
    }
  }
  control <- list(maxit = maxit, factr = reltol / .Machine$double.eps,
    pgtol = pgtol, trace = as.integer(verbose))
  opt <- optim(z0, objective, gr = gradient, method = "L-BFGS-B",
    lower = lower, upper = upper, control = control)
  full <- expand(opt$par)
  params <- .piecewise_calendar_revision_monthly_unpack(full,
    "joint_marginal")
  evaluation <- e_step_eps_4w(df, params, check_df = FALSE)
  list(params = params, loglik = evaluation$loglik,
    gamma = evaluation$gamma,
    job_change_posterior = evaluation$job_change_posterior,
    duration_reliability_posterior =
      evaluation$duration_reliability_posterior,
    convergence = opt$convergence, message = opt$message,
    iterations = unname(opt$counts["function"]), objective = opt$value,
    par_unconstrained = full, objective_function = objective,
    free_names = free_names)
}

duration_weighted_transition_rates_4w <- function(df, fit) {
  hmat <- latent_histories_eps_4w()
  params <- fit$params
  exact_risk <- identical(fit$integration_method, "exact_piecewise")
  rows <- lapply(1:3, function(wave) {
    employed <- as.vector(fit$gamma %*% hmat[, wave])
    nonemployed <- 1 - employed
    exit <- .duration_transition_probability(
      df[[paste0("tenure", wave)]], params$lambda_g, params$beta_g)
    entry <- if (exact_risk) .piecewise_category_transition_risk_exact(
      df[[paste0("timegap_cat", wave)]], params$lambda_d) else
      .duration_category_transition_probability(
        df[[paste0("timegap_cat", wave)]], params$lambda_d, params$beta_d)
    exit[!is.finite(exit)] <- if (exact_risk)
      .piecewise_mean_transition_risk_exact(params$lambda_g) else
      .duration_marginal_transition_probability(params$lambda_g, params$beta_g)
    entry[!is.finite(entry)] <- if (exact_risk)
      .piecewise_mean_transition_risk_exact(params$lambda_d) else
      .duration_marginal_transition_probability(params$lambda_d, params$beta_d)
    data.frame(wave = wave,
      exit_rate = weighted.mean(exit, df$weight * employed),
      entry_rate = weighted.mean(entry, df$weight * nonemployed),
      exit_numerator = sum(df$weight * employed * exit),
      exit_denominator = sum(df$weight * employed),
      entry_numerator = sum(df$weight * nonemployed * entry),
      entry_denominator = sum(df$weight * nonemployed))
  })
  out <- do.call(rbind, rows)
  pooled <- data.frame(wave = 0L,
    exit_rate = sum(out$exit_numerator) / sum(out$exit_denominator),
    entry_rate = sum(out$entry_numerator) / sum(out$entry_denominator),
    exit_numerator = sum(out$exit_numerator),
    exit_denominator = sum(out$exit_denominator),
    entry_numerator = sum(out$entry_numerator),
    entry_denominator = sum(out$entry_denominator))
  rbind(pooled, out)
}

# A bounded alternative to an optimizer-controlled line search. One numerical
# gradient is evaluated, followed by at most `length(backtrack)` likelihood
# evaluations along a normalized descent direction. This makes long four-wave
# continuation stages predictable and independently cacheable.
four_wave_gradient_step <- function(df, start, free_names, workers = 2L,
    gradient_step = 1e-4, max_parameter_move = .25,
    backtrack = c(1, .5, .25, .125), scheme = c("forward", "central")) {
  scheme <- match.arg(scheme)
  validate_df_eps_4w(df, allow_zero_tenure = TRUE)
  full0 <- .pack_four_wave_preferred(start)
  if (!length(free_names) || length(setdiff(free_names, names(full0))))
    stop("unknown or empty four-wave gradient block")
  z0 <- full0[free_names]
  bounds <- .four_wave_parameter_bounds(full0)
  lower <- bounds$lower[free_names]
  upper <- bounds$upper[free_names]
  if (any(z0 < lower | z0 > upper)) stop("Starting parameters violate bounds")
  total_weight <- sum(df$weight)
  objective <- function(z) {
    full <- full0
    full[free_names] <- z
    params <- .piecewise_calendar_revision_monthly_unpack(full,
      "joint_marginal")
    value <- tryCatch(e_step_eps_4w(df, params,
      check_df = FALSE)$loglik, error = function(error) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value / total_weight
  }
  centre <- objective(z0)
  if (!is.finite(centre) || centre >= 1e99)
    stop("Starting likelihood is not finite")
  step <- gradient_step * pmax(1, abs(z0))
  plus_step <- pmin(step, upper - z0)
  minus_step <- pmin(step, z0 - lower)
  points <- vector("list", if (scheme == "central") 2L * length(z0) else
    length(z0))
  if (scheme == "central") {
    for (index in seq_along(z0)) {
      points[[2L * index - 1L]] <- points[[2L * index]] <- z0
      points[[2L * index - 1L]][index] <- z0[index] + plus_step[index]
      points[[2L * index]][index] <- z0[index] - minus_step[index]
    }
  } else for (index in seq_along(z0)) {
    points[[index]] <- z0
    step[index] <- if (plus_step[index] >= step[index]) step[index] else
      -minus_step[index]
    points[[index]][index] <- z0[index] + step[index]
  }
  workers <- max(1L, min(8L, as.integer(workers)))
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster, function(path) {
      setwd(path)
      source("EM-tenure/R/source_all.R")
      NULL
    }, worker_wd)
    df_worker <- df
    full_worker <- full0
    free_worker <- free_names
    weight_worker <- total_weight
    parallel::clusterExport(cluster, c("df_worker", "full_worker",
      "free_worker", "weight_worker"), envir = environment())
    worker_objective <- function(z) {
      full <- full_worker
      full[free_worker] <- z
      params <- .piecewise_calendar_revision_monthly_unpack(full,
        "joint_marginal")
      value <- tryCatch(e_step_eps_4w(df_worker, params,
        check_df = FALSE)$loglik, error = function(error) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value / weight_worker
    }
    values <- unlist(parallel::parLapplyLB(cluster, points,
      worker_objective), use.names = FALSE)
  } else values <- vapply(points, objective, numeric(1))
  gradient <- if (scheme == "central")
    (values[seq(1L, length(values), 2L)] -
      values[seq(2L, length(values), 2L)]) / (plus_step + minus_step) else
    (values - centre) / step
  if (any(!is.finite(gradient)) || any(values >= 1e99))
    stop("Invalid finite-difference likelihood; gradient step was not accepted")
  names(gradient) <- names(z0)
  if (!any(abs(gradient) > 0)) stop("four-wave numerical gradient is zero")
  direction <- -gradient / max(abs(gradient)) * max_parameter_move
  proposals <- lapply(backtrack, function(scale)
    pmin(upper, pmax(lower, z0 + scale * direction)))
  proposal_values <- vapply(proposals, objective, numeric(1))
  candidates <- c(centre, proposal_values)
  selected <- which.min(candidates)
  z_best <- if (selected == 1L) z0 else proposals[[selected - 1L]]
  full_best <- full0
  full_best[free_names] <- z_best
  params <- .piecewise_calendar_revision_monthly_unpack(full_best,
    "joint_marginal")
  evaluation <- e_step_eps_4w(df, params, check_df = FALSE)
  list(params = params, loglik = evaluation$loglik,
    gamma = evaluation$gamma,
    job_change_posterior = evaluation$job_change_posterior,
    duration_reliability_posterior =
      evaluation$duration_reliability_posterior,
    convergence = 1L, converged = FALSE, step_accepted = selected != 1L,
    message = if (selected == 1L) "No improving bounded step" else
      "Accepted bounded gradient step",
    iterations = 1L, objective = candidates[selected],
    par_unconstrained = full_best, free_names = free_names,
    gradient = gradient, gradient_max = max(abs(gradient)),
    gradient_at = full0, gradient_scheme = scheme,
    backtrack = data.frame(scale = c(0, backtrack),
      objective = candidates, selected = seq_along(candidates) == selected))
}
