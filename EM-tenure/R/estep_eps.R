# ==============================================================================
# EM-tenure: E-step for the eps (Spec I) model
# ==============================================================================
# Created: 2026-04-30
#
# Computes responsibilities gamma_{ih} = P(H = h | y_i, Theta) and sufficient
# statistics under Spec I (point-mass + Exp contamination, no sigma_g).
#
# Key differences from base e_step():
#   1. Tenure emissions decompose by maximal latent E-spell (per history h).
#      For a spell with K observed tenures (K = sum_t s_t in spell):
#        - K = 0: no contribution.
#        - K = 1: log Exp(g_obs).
#        - K >= 2: log-sum-exp of clean (clock-consistent point mass) and
#          contaminated (iid Exp) branches.
#   2. Miscl-employed (h_t = 0, s_t = 1) waves contribute a singleton Exp
#      tenure emission.
#   3. Timegap (d) emissions reuse the base interval-censored Exp(lambda_d)
#      machinery (marginal + transition), unchanged from the base / rho models.
#   4. No sigma_g sufficient statistics (sigma is dropped).
#   5. Per-spell posterior contamination probability tau and spell length K
#      are accumulated into Eps_num and Eps_den for the eps M-step.
#   6. Exp-emission counts (Exp_count) and tenure x-sums (Exp_xsum) are
#      accumulated for the lambda_g M-step (handles clean vs contaminated
#      decomposition; see emissions_eps.R for derivation).
#
# Companion: documents/EM tenure epsilon.tex (Section 4).
# Only supports discrete_timegap = TRUE.
# ==============================================================================

#' E-step for the eps (Spec I) model
#'
#' Computes responsibilities and sufficient statistics for the spell-pair
#' point-mass + Exp contamination model.
#'
#' @param df Data frame with columns: y1-y3, tenure1-tenure3,
#'   timegap_cat1-timegap_cat3, weight.
#' @param params Named list: alpha, theta0, theta1, pi, eps, lambda_g,
#'   lambda_d. (No sigma2_g.)
#' @return Named list with:
#'   \describe{
#'     \item{gamma}{N x 8 matrix of responsibilities.}
#'     \item{loglik}{Weighted observed-data log-likelihood.}
#'     \item{suff}{Sufficient statistics for \code{m_step_eps()}.}
#'   }
#' @references TeX: \emph{EM tenure epsilon.tex}, Section 4.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y1 = 1L, y2 = 1L, y3 = 0L,
#'   tenure1 = 2.0, tenure2 = 2.25, tenure3 = NA_real_,
#'   timegap_cat1 = 4L, timegap_cat2 = 4L, timegap_cat3 = 2L,
#'   weight = 1.0
#' )
#' params <- list(alpha = 0.6, theta0 = 0.1, theta1 = 0.9, pi = 0.05,
#'                eps = 0.20, lambda_g = 0.15, lambda_d = 0.5)
#' e_step_eps(df, params)
#' }
#' Validate data-frame inputs for the eps E-step
#'
#' Checks column presence, binary y, employed-wave tenure NA/non-positive,
#' timegap category range, and weight validity. Called once per EM run from
#' \code{em_fit_tenure_eps()} so that \code{e_step_eps()} can skip these
#' O(N) traversals on every iteration (via \code{check_df = FALSE}).
#'
#' @param df Data frame; same as for \code{e_step_eps()}.
#' @return Invisibly NULL; stops with an informative error on failure.
#' @keywords internal
validate_df_eps <- function(df, allow_zero_tenure = FALSE) {
  cat_cols <- c("timegap_cat1", "timegap_cat2", "timegap_cat3")
  missing_cats <- setdiff(cat_cols, names(df))
  if (length(missing_cats) > 0)
    stop("e_step_eps requires columns: ", paste(missing_cats, collapse = ", "))
  # Check y before tenure: NA in y would propagate through y == 1L comparisons.
  na_y <- is.na(df$y1) | is.na(df$y2) | is.na(df$y3)
  if (any(na_y))
    stop(sprintf("e_step_eps: %d obs have NA in y1/y2/y3.", sum(na_y)))
  bad_y <- !all(df$y1 %in% 0:1) || !all(df$y2 %in% 0:1) || !all(df$y3 %in% 0:1)
  if (bad_y) stop("e_step_eps: y1/y2/y3 must be binary (0 or 1).")
  na_timegap_nonemp <- (df$y1 == 0L & is.na(df$timegap_cat1)) |
                       (df$y2 == 0L & is.na(df$timegap_cat2)) |
                       (df$y3 == 0L & is.na(df$timegap_cat3))
  if (any(na_timegap_nonemp))
    stop(sprintf("e_step_eps: %d obs have NA timegap at a nonemployed wave.",
                 sum(na_timegap_nonemp)))
  bad_cats <- any(!is.na(df$timegap_cat1) & !df$timegap_cat1 %in% 1:7) ||
              any(!is.na(df$timegap_cat2) & !df$timegap_cat2 %in% 1:7) ||
              any(!is.na(df$timegap_cat3) & !df$timegap_cat3 %in% 1:7)
  if (bad_cats) stop("e_step_eps: observed timegap categories must be integers 1-7.")
  na_tenure_emp <- (df$y1 == 1L & is.na(df$tenure1)) |
                   (df$y2 == 1L & is.na(df$tenure2)) |
                   (df$y3 == 1L & is.na(df$tenure3))
  if (any(na_tenure_emp))
    stop(sprintf("e_step_eps: %d obs have NA tenure at an employed wave.",
                 sum(na_tenure_emp)))
  bad_tenure <- if (allow_zero_tenure) {
    (df$y1 == 1 & df$tenure1 < 0) |
      (df$y2 == 1 & df$tenure2 < 0) |
      (df$y3 == 1 & df$tenure3 < 0)
  } else {
    (df$y1 == 1 & df$tenure1 <= 0) |
      (df$y2 == 1 & df$tenure2 <= 0) |
      (df$y3 == 1 & df$tenure3 <= 0)
  }
  if (any(bad_tenure))
    stop(sprintf("e_step_eps: %d obs have invalid tenure for employed state.",
      sum(bad_tenure)))
  if (is.null(df$weight)) stop("e_step_eps: df must contain a 'weight' column.")
  if (any(is.na(df$weight)))
    stop(sprintf("e_step_eps: %d obs have NA weight.", sum(is.na(df$weight))))
  if (any(df$weight <= 0))
    stop(sprintf("e_step_eps: %d obs have non-positive weight.", sum(df$weight <= 0)))
  invisible(NULL)
}

.log_duration_history_prior_eps <- function(hmat, alpha, s_list, g_list, c_list,
                                             lambda_g, beta_g,
                                             lambda_d, beta_d) {
  N <- length(g_list[[1]])
  H <- nrow(hmat)
  out <- matrix(0, N, H)
  p_entry <- lapply(c_list[1:2], function(z)
    .duration_category_transition_probability(z, lambda_d, beta_d))
  p_exit_missing <- .duration_marginal_transition_probability(
    lambda_g, beta_g)
  p_entry_missing <- .duration_marginal_transition_probability(
    lambda_d, beta_d)
  for (j in seq_len(H)) {
    h <- hmat[j, ]
    out[, j] <- if (h[1] == 1L) log(alpha) else log1p(-alpha)
    for (t in 1:2) {
      if (h[t] == 1L) {
        # Tenure is observed only when employment is reported.  Under a
        # latent-employed/reported-nonemployed mismatch, integrate the
        # unavailable spell age out rather than using an imputed clock.
        p_change <- .duration_transition_probability(
          g_list[[t]], lambda_g, beta_g)
        p_change[s_list[[t]] == 0L | !is.finite(p_change)] <- p_exit_missing
      } else {
        p_change <- p_entry[[t]]
        p_change[s_list[[t]] == 1L | !is.finite(p_change)] <- p_entry_missing
      }
      p_change <- pmin(pmax(p_change, 1e-12), 1 - 1e-12)
      out[, j] <- out[, j] + if (h[t + 1L] == h[t]) {
        log1p(-p_change)
      } else {
        log(p_change)
      }
    }
  }
  out
}

#' @export
e_step_eps <- function(df, params, check_df = TRUE, suff_stats = TRUE) {
  reliability_shifts <- .duration_reliability_shifts(params)

  # A fixed equal-weight two-point mixture links tenure and timegap reporting
  # reliability across all waves of a record. The reliable class shifts the
  # contamination logits down and the unreliable class shifts them up. The
  # shifts may be common or clock-specific. At zero shifts the two components
  # coincide exactly with the preceding model.
  if (any(reliability_shifts>0)) {
    if (suff_stats)
      stop("duration-reliability mixtures require suff_stats=FALSE")
    reliable_params <- .duration_reliability_component_params(params,
      "reliable")
    unreliable_params <- .duration_reliability_component_params(params,
      "unreliable")
    reliable <- e_step_eps(df,reliable_params,check_df=check_df,
      suff_stats=FALSE)
    unreliable <- e_step_eps(df,unreliable_params,check_df=FALSE,
      suff_stats=FALSE)
    row_max <- pmax(reliable$row_loglik,unreliable$row_loglik)
    row_loglik <- row_max+log(.5*exp(reliable$row_loglik-row_max)+
      .5*exp(unreliable$row_loglik-row_max))
    posterior_unreliable <- .5*exp(unreliable$row_loglik-row_loglik)
    posterior_reliable <- 1-posterior_unreliable
    gamma <- reliable$gamma*posterior_reliable+
      unreliable$gamma*posterior_unreliable
    job_change_posterior <- list(
      expected_changes=reliable$job_change_posterior$expected_changes*
        posterior_reliable+
        unreliable$job_change_posterior$expected_changes*
          posterior_unreliable,
      opportunities=reliable$job_change_posterior$opportunities*
        posterior_reliable+
        unreliable$job_change_posterior$opportunities*
          posterior_unreliable)
    return(list(gamma=gamma,loglik=sum(df$weight*row_loglik),
      row_loglik=row_loglik,suff=NULL,
      job_change_posterior=job_change_posterior,
      duration_reliability_posterior=posterior_unreliable,
      duration_reliability_component_probabilities=data.frame(
        class=c("reliable","unreliable"),
        tenure_contamination=c(reliable_params$eps,
          unreliable_params$eps),
        timegap_contamination=c(reliable_params$eps_d,
          unreliable_params$eps_d))))
  }
  # --- Validate df inputs (skipped from the EM loop for efficiency) ---
  monthly_tenure <- identical(params$tenure_measurement_model, "monthly")
  if (check_df) validate_df_eps(df, allow_zero_tenure = monthly_tenure)
  # --- Unpack parameters ---
  alpha    <- params$alpha
  theta0   <- params$theta0
  theta1   <- params$theta1
  pi_par   <- params$pi
  eps      <- params$eps
  tenure_model <- if (is.null(params$tenure_contamination_model))
    "marginal" else params$tenure_contamination_model
  eps_local <- if (is.null(params$eps_local)) 0 else params$eps_local
  eps_gross <- if (is.null(params$eps_gross)) eps else params$eps_gross
  eps_d    <- if (is.null(params$eps_d)) 0 else params$eps_d
  job_change_prob <- if (is.null(params$job_change_prob)) 0 else
    params$job_change_prob
  tenure_report_persistence <- if (is.null(params$tenure_report_persistence))
    0 else params$tenure_report_persistence
  tenure_heaping_prob <- if (is.null(params$tenure_heaping_prob)) 0 else
    params$tenure_heaping_prob
  tenure_year_revision_prob <-
    if (is.null(params$tenure_year_revision_prob)) 0 else
      params$tenure_year_revision_prob
  tenure_clean_anchor_revision_prob <-
    if (is.null(params$tenure_clean_anchor_revision_prob)) 0 else
      params$tenure_clean_anchor_revision_prob
  tenure_exact_anchor_retention_prob <-
    if (is.null(params$tenure_exact_anchor_retention_prob)) 0 else
      params$tenure_exact_anchor_retention_prob
  tenure_local_revision_prob <-
    if (is.null(params$tenure_local_revision_prob)) 0 else
      params$tenure_local_revision_prob
  tenure_start_month_probs <- if (is.null(params$tenure_start_month_probs))
    rep(1/12,12L) else .validate_start_month_probs(
      params$tenure_start_month_probs)
  seasonal_start_month <- max(abs(tenure_start_month_probs-1/12))>1e-12
  timegap_model <- if (is.null(params$timegap_contamination_model))
    "marginal" else params$timegap_contamination_model
  local_decay <- if (is.null(params$timegap_local_decay)) 1 else
    params$timegap_local_decay
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d
  duration_dependent <- !is.null(params$beta_g) || !is.null(params$beta_d)
  beta_g <- if (is.null(params$beta_g)) 0 else params$beta_g
  beta_d <- if (is.null(params$beta_d)) 0 else params$beta_d

  if (!is.finite(eps) || eps <= 0 || eps >= 1) {
    stop(sprintf("e_step_eps: params$eps must be in (0, 1); got %.4g", eps))
  }
  if (!tenure_model %in% c("marginal", "local_gross"))
    stop("e_step_eps: unknown tenure contamination model")
  if (identical(tenure_model, "local_gross") &&
      (!is.finite(eps_local) || !is.finite(eps_gross) || eps_local < 0 ||
       eps_gross < 0 || abs(eps - eps_local - eps_gross) > 1e-8))
    stop("e_step_eps: invalid local/gross tenure-error probabilities")
  if (!is.finite(eps_d) || eps_d < 0 || eps_d >= 1) {
    stop(sprintf("e_step_eps: params$eps_d must be in [0, 1); got %.4g", eps_d))
  }
  if (!is.finite(job_change_prob) || job_change_prob < 0 ||
      job_change_prob >= 1)
    stop("e_step_eps: params$job_change_prob must be in [0,1)")
  if (!is.finite(tenure_report_persistence) ||
      tenure_report_persistence < 0 || tenure_report_persistence >= 1)
    stop("e_step_eps: params$tenure_report_persistence must be in [0,1)")
  if (tenure_report_persistence > 0 && !monthly_tenure)
    stop("e_step_eps: tenure-report persistence requires monthly tenure")
  if (!is.finite(tenure_heaping_prob) || tenure_heaping_prob < 0 ||
      tenure_heaping_prob >= 1)
    stop("e_step_eps: params$tenure_heaping_prob must be in [0,1)")
  if (tenure_heaping_prob > 0 && !monthly_tenure)
    stop("e_step_eps: tenure heaping requires monthly tenure")
  if (tenure_heaping_prob > 0 && tenure_report_persistence > 0)
    stop("e_step_eps: persistence and calendar heaping are not jointly implemented")
  if (!is.finite(tenure_year_revision_prob) ||
      tenure_year_revision_prob < 0 || tenure_year_revision_prob >= 1)
    stop("e_step_eps: params$tenure_year_revision_prob must be in [0,1)")
  if (tenure_year_revision_prob > 0 && !monthly_tenure)
    stop("e_step_eps: whole-year revisions require monthly tenure")
  if (tenure_year_revision_prob > 0 && tenure_report_persistence > 0)
    stop("e_step_eps: persistence and whole-year revisions are not jointly implemented")
  if (!is.finite(tenure_clean_anchor_revision_prob) ||
      tenure_clean_anchor_revision_prob < 0 ||
      tenure_clean_anchor_revision_prob >= 1)
    stop(paste0("e_step_eps: params$tenure_clean_anchor_revision_prob ",
      "must be in [0,1)"))
  if (tenure_clean_anchor_revision_prob > 0 && !monthly_tenure)
    stop("e_step_eps: clean-anchor revisions require monthly tenure")
  if (tenure_clean_anchor_revision_prob > 0 &&
      tenure_report_persistence > 0)
    stop("e_step_eps: persistence and clean-anchor revisions are not jointly implemented")
  if (!is.finite(tenure_exact_anchor_retention_prob) ||
      tenure_exact_anchor_retention_prob < 0 ||
      tenure_exact_anchor_retention_prob >= 1)
    stop(paste0("e_step_eps: params$tenure_exact_anchor_retention_prob ",
      "must be in [0,1)"))
  if (tenure_exact_anchor_retention_prob > 0 && !monthly_tenure)
    stop("e_step_eps: exact-anchor retention requires monthly tenure")
  if (tenure_exact_anchor_retention_prob > 0 &&
      tenure_report_persistence > 0)
    stop("e_step_eps: persistence and exact-anchor retention are not jointly implemented")
  if (!is.finite(tenure_local_revision_prob) ||
      tenure_local_revision_prob < 0 || tenure_local_revision_prob >= 1)
    stop(paste0("e_step_eps: params$tenure_local_revision_prob ",
      "must be in [0,1)"))
  if (tenure_local_revision_prob > 0 && !monthly_tenure)
    stop("e_step_eps: local revisions require monthly tenure")
  if (tenure_local_revision_prob > 0 && tenure_report_persistence > 0)
    stop("e_step_eps: persistence and local revisions are not jointly implemented")
  interview_cols <- paste0("interview_month", 1:3)
  if ((tenure_heaping_prob > 0 || tenure_year_revision_prob > 0 ||
      tenure_clean_anchor_revision_prob > 0 ||
      tenure_exact_anchor_retention_prob > 0 ||
      tenure_local_revision_prob > 0 || seasonal_start_month) &&
      (length(setdiff(interview_cols, names(df))) ||
       any(!as.matrix(df[interview_cols]) %in% 1:12)))
    stop("e_step_eps: calendar revisions require interview_month1-3 in 1:12")
  if (job_change_prob > 0 && identical(tenure_model, "local_gross"))
    stop("e_step_eps: job changes are currently implemented for marginal tenure contamination only")
  if (!timegap_model %in% c("marginal","local","joint_marginal"))
    stop("e_step_eps: unknown timegap contamination model")
  if (!is.finite(local_decay) || local_decay<=0)
    stop("e_step_eps: timegap_local_decay must be positive")
  if (any(!is.finite(lambda_g)) || any(lambda_g <= 0)) {
    stop("e_step_eps: all params$lambda_g hazards must be finite and positive")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop(sprintf("e_step_eps: params$alpha must be in (0, 1); got %.4g", alpha))
  }
  if (!duration_dependent && (!is.finite(theta0) || theta0 <= 0 || theta0 >= 1)) {
    stop(sprintf("e_step_eps: params$theta0 must be in (0, 1); got %.4g", theta0))
  }
  if (!duration_dependent && (!is.finite(theta1) || theta1 <= 0 || theta1 >= 1)) {
    stop(sprintf("e_step_eps: params$theta1 must be in (0, 1); got %.4g", theta1))
  }
  if (!is.finite(pi_par) || pi_par < 0 || pi_par >= 1) {
    stop(sprintf("e_step_eps: params$pi must be in [0, 1); got %.4g", pi_par))
  }
  if (any(!is.finite(lambda_d)) || any(lambda_d <= 0)) {
    stop("e_step_eps: all params$lambda_d hazards must be finite and positive")
  }
  if (!all(is.finite(c(beta_g, beta_d))) || beta_g <= -1 || beta_d <= -1) {
    stop("e_step_eps: beta_g and beta_d must be finite and greater than -1")
  }

  # --- Latent structure (H = 8) ---
  hmat      <- latent_histories()

  N  <- nrow(df)
  H  <- nrow(hmat)
  h1 <- hmat[, 1]; h2 <- hmat[, 2]; h3 <- hmat[, 3]

  # --- Extract data as N-vectors ---
  s1 <- df$y1; s2 <- df$y2; s3 <- df$y3

  # --- Pre-compute mismatch matrix (N×H) for lq and M-step suff-stats ---
  # H_mat: 3×H binary history matrix; s_full: N×3 logical observation matrix.
  # n_mis_mat[i,j] = # waves where latent history h_j differs from obs s_i.
  # (s1+s2+s3) avoids logical→integer coercion of rowSums(s_full). [P3-24]
  H_mat     <- rbind(h1, h2, h3)
  s_full    <- cbind(s1 == 1L, s2 == 1L, s3 == 1L)
  n_mis_mat <- (s1 + s2 + s3) +
               matrix(colSums(H_mat), N, H, byrow = TRUE) -
               2L * (s_full %*% H_mat)   # BLAS dgemm: N×H in one pass [P1-1]
  g1 <- df$tenure1; g2 <- df$tenure2; g3 <- df$tenure3
  c1 <- df$timegap_cat1; c2 <- df$timegap_cat2; c3 <- df$timegap_cat3
  wi <- df$weight  # already validated by validate_df_eps when check_df = TRUE

  s_list <- list(s1, s2, s3)
  g_list <- list(g1, g2, g3)
  c_list <- list(c1, c2, c3)

  if (duration_dependent) {
    log_prior <- .log_duration_history_prior_eps(
      hmat, alpha, s_list, g_list, c_list,
      lambda_g, beta_g, lambda_d, beta_d)
  } else {
    prior_h <- prior_over_histories(hmat, theta1, theta0, alpha)
    log_prior <- log(.bound01(prior_h))
  }

  # --- Misclassification log-probability: N x H ---
  # Vectorised via n_mis_mat: no per-history loop or ifelse() allocation. [P1-1]
  if (pi_par < .Machine$double.eps) {
    lq <- matrix(-Inf, nrow = N, ncol = H)
    lq[n_mis_mat == 0L] <- 0
  } else {
    pi_b <- .bound01(pi_par)
    lq   <- n_mis_mat * log(pi_b) + (3L - n_mis_mat) * log1p(-pi_b)
  }

  # --- Per-history emission accumulators ---
  ld               <- matrix(0, nrow = N, ncol = H)
  lambda_count_mat <- matrix(0, nrow = N, ncol = H)
  lambda_xsum_mat  <- matrix(0, nrow = N, ncol = H)
  eps_num_mat      <- matrix(0, nrow = N, ncol = H)  # tau_sum (E[# contam waves], eps-informative spells)
  eps_den_mat      <- matrix(0, nrow = N, ncol = H)  # K (# obs waves, eps-informative spells)
  job_num_mat      <- matrix(0, nrow = N, ncol = H)
  job_den_mat      <- matrix(0, nrow = N, ncol = H)

  # Pre-compute full N x 3 tenure matrix (s_full is already computed above).
  g_full <- cbind(g1, g2, g3)
  interview_month_full <- if (all(interview_cols %in% names(df)))
    as.matrix(df[interview_cols]) else NULL

  # --- Loop over histories ---
  for (j in seq_len(H)) {
    h_j    <- as.integer(hmat[j, ])
    spells <- .maximal_e_spells(h_j)

    # ---- (a) Tenure spell emissions (per maximal E-spell in h_j) ----
    for (spell in spells) {
      L         <- length(spell)
      g_mat     <- g_full[, spell, drop = FALSE]
      s_mat     <- s_full[, spell, drop = FALSE]
      t_offsets <- as.integer(spell - spell[1L])

      out <- if (monthly_tenure)
        log_emission_spell_g_monthly(g_mat, s_mat, t_offsets, lambda_g,
          eps, job_change_prob, beta_g = beta_g,
          initial_model = if (spell[1L] > 1L) "within_interval" else
            "marginal",
          tenure_report_persistence=tenure_report_persistence,
          tenure_heaping_prob=tenure_heaping_prob,
          tenure_year_revision_prob=tenure_year_revision_prob,
          tenure_clean_anchor_revision_prob=
            tenure_clean_anchor_revision_prob,
          tenure_exact_anchor_retention_prob=
            tenure_exact_anchor_retention_prob,
          tenure_local_revision_prob=tenure_local_revision_prob,
          tenure_start_month_probs=tenure_start_month_probs,
          interview_month_mat=if (is.null(interview_month_full)) NULL else
            interview_month_full[, spell, drop=FALSE]) else if
        (job_change_prob > 0)
        log_emission_spell_g_job_change(g_mat, s_mat, t_offsets, lambda_g,
          eps, job_change_prob, beta_g = beta_g) else if
        (identical(tenure_model, "local_gross"))
        log_emission_spell_g_local_gross(g_mat, s_mat, t_offsets, lambda_g,
          eps_local, eps_gross, beta_g = beta_g) else
        log_emission_spell_g(g_mat, s_mat, t_offsets, lambda_g, eps,
          beta_g = beta_g)

      ld[, j]               <- ld[, j]               + out$loglik
      lambda_count_mat[, j] <- lambda_count_mat[, j] + out$lambda_count
      lambda_xsum_mat[, j]  <- lambda_xsum_mat[, j]  + out$lambda_xsum
      if (!is.null(out$job_changes)) {
        job_num_mat[, j] <- job_num_mat[, j] + out$job_changes
        job_den_mat[, j] <- job_den_mat[, j] + out$job_change_opportunities
      }

      # eps stats: accumulate for all eps-informative spells (K >= 2, or
      # K = 1 with offset > 0 where clean/contaminated branches differ).
      mask_eps <- out$eps_informative
      if (any(mask_eps)) {
        eps_num_mat[mask_eps, j] <- eps_num_mat[mask_eps, j] + out$tau_sum[mask_eps]
        eps_den_mat[mask_eps, j] <- eps_den_mat[mask_eps, j] + out$K[mask_eps]
      }
    }

    # ---- (b) Miscl-employed singleton tenure emissions (h_t=0, s_t=1) ----
    for (t in 1:3) {
      if (h_j[t] == 0L) {
        s_t  <- s_list[[t]]
        g_t  <- g_list[[t]]
        mask <- (s_t == 1L)
        if (any(mask)) {
          marginal_g <- if (monthly_tenure && seasonal_start_month)
            .log_calendar_duration_month_mass(round(12*g_t[mask]),
              df[[interview_cols[t]]][mask],lambda_g,beta_g,
              tenure_start_month_probs) else if (monthly_tenure)
            .log_duration_month_mass(round(12*g_t[mask]),lambda_g,
              beta_g) else
            .log_duration_density(g_t[mask], lambda_g, beta_g)
          if (monthly_tenure && tenure_heaping_prob > 0) {
            heaped_g <- .log_january_duration_month_mass(
              round(12 * g_t[mask]), df[[interview_cols[t]]][mask],
              lambda_g, beta_g)
            marginal_g <- .log_probability_mixture(marginal_g, heaped_g,
              tenure_heaping_prob)
          }
          ld[mask, j] <- ld[mask, j] + marginal_g
          lambda_count_mat[mask, j] <- lambda_count_mat[mask, j] + 1
          lambda_xsum_mat[mask, j]  <- lambda_xsum_mat[mask, j]  + g_t[mask]
        }
      }
    }

    # ---- (c) Timegap emissions: marginal + transition (mirror base model) ----
    # The joint model differs only for a fully observed three-wave latent
    # nonemployment spell. Two-wave blocks reduce to a conditional mixture with
    # contamination probability 1-(1-eps_d)^2.
    joint3 <- identical(timegap_model,"joint_marginal") && all(h_j==0L) &
      eps_d>0
    joint3_mask <- if (joint3) s1==0L & s2==0L & s3==0L else rep(FALSE,N)
    if (any(joint3_mask)) ld[joint3_mask,j] <- ld[joint3_mask,j] +
      log_emission_timegap_triplet_joint(c1[joint3_mask],c2[joint3_mask],
        c3[joint3_mask],lambda_d,beta_d,eps_d)
    pair_eps_d <- if (identical(timegap_model,"joint_marginal"))
      1-(1-eps_d)^2 else eps_d
    conditional_model <- if (identical(timegap_model,"local")) "local" else
      "marginal"
    # Wave 1: marginal interval whenever s_1 = 0.
    mask_w1 <- (s1 == 0L) & !joint3_mask
    if (any(mask_w1)) {
      ld[mask_w1, j] <- ld[mask_w1, j] +
        log_emission_interval_d(c1[mask_w1], lambda_d, beta_d)
    }
    # Waves 2 and 3: inlined timegap transitions (avoids 2× rep(0,N) alloc) [P2-9]
    # --- Transition (2 ← 1) ---
    {
      hp12 <- h_j[1L]; hc12 <- h_j[2L]
      if (hp12 == 0L && hc12 == 0L) {
        m12 <- (s2 == 0L) & (s1 == 0L) & !joint3_mask
        if (any(m12))  ld[m12, j] <- ld[m12, j] + log_emission_transition_d_contaminated(c2[m12], c1[m12], lambda_d, beta_d, pair_eps_d,conditional_model,local_decay)
        m12m <- (s2 == 0L) & (s1 == 1L)
        if (any(m12m)) ld[m12m, j] <- ld[m12m, j] + log_emission_interval_d(c2[m12m], lambda_d, beta_d)
      } else if (hp12 == 1L && hc12 == 0L) {
        m12 <- (s2 == 0L)
        if (any(m12)) ld[m12, j] <- ld[m12, j] + log_emission_interval_d(c2[m12], lambda_d, beta_d)
      } else if (hc12 == 1L) {
        m12 <- (s2 == 0L)
        if (any(m12)) ld[m12, j] <- ld[m12, j] + log_emission_interval_d(c2[m12], lambda_d, beta_d)
      }
    }
    # --- Transition (3 ← 2) ---
    {
      hp23 <- h_j[2L]; hc23 <- h_j[3L]
      if (hp23 == 0L && hc23 == 0L) {
        m23 <- (s3 == 0L) & (s2 == 0L) & !joint3_mask
        if (any(m23))  ld[m23, j] <- ld[m23, j] + log_emission_transition_d_contaminated(c3[m23], c2[m23], lambda_d, beta_d, pair_eps_d,conditional_model,local_decay)
        m23m <- (s3 == 0L) & (s2 == 1L)
        if (any(m23m)) ld[m23m, j] <- ld[m23m, j] + log_emission_interval_d(c3[m23m], lambda_d, beta_d)
      } else if (hp23 == 1L && hc23 == 0L) {
        m23 <- (s3 == 0L)
        if (any(m23)) ld[m23, j] <- ld[m23, j] + log_emission_interval_d(c3[m23], lambda_d, beta_d)
      } else if (hc23 == 1L) {
        m23 <- (s3 == 0L)
        if (any(m23)) ld[m23, j] <- ld[m23, j] + log_emission_interval_d(c3[m23], lambda_d, beta_d)
      }
    }
  }

  # --- Posterior responsibilities: N x H ---
  if (pi_par < .Machine$double.eps) {
    # Encode (s1, s2, s3) as a 1-based integer index into the 8 history rows:
    # index = 1 + s1*2^0 + s2*2^1 + s3*2^2  (binary-to-row mapping)
    match_idx <- 1L + s1 + 2L * s2 + 4L * s3
    gamma_mat <- matrix(0, nrow = N, ncol = H)
    gamma_mat[cbind(seq_len(N), match_idx)] <- 1
    lp_match <- if (is.matrix(log_prior)) {
      log_prior[cbind(seq_len(N), match_idx)]
    } else log_prior[match_idx]
    log_lik_i <- lp_match + ld[cbind(seq_len(N), match_idx)]
    ll <- sum(wi * log_lik_i)
  } else {
    log_kernel <- if (is.matrix(log_prior)) lq + ld + log_prior else
      sweep(lq + ld, 2, log_prior, "+")
    # log-sum-exp normalization (H=8 is fixed so explicit-column pmax avoids the
    # as.data.frame() copy overhead of do.call(pmax, ...))
    row_max   <- pmax(log_kernel[,1], log_kernel[,2], log_kernel[,3], log_kernel[,4],
                      log_kernel[,5], log_kernel[,6], log_kernel[,7], log_kernel[,8])
    log_denom <- row_max + log(rowSums(exp(log_kernel - row_max)))
    gamma_mat <- exp(log_kernel - log_denom)
    log_lik_i <- log_denom
    ll        <- sum(wi * log_denom)
  }

  if (!suff_stats) {
    return(list(gamma = gamma_mat, loglik = ll, row_loglik = log_lik_i,
      suff = NULL,
      job_change_posterior = list(
        expected_changes = rowSums(gamma_mat * job_num_mat),
        opportunities = rowSums(gamma_mat * job_den_mat))))
  }

  # --- Sufficient statistics ---
  wg    <- gamma_mat * wi
  wg_cs <- colSums(wg)

  # Markov stats (identical to base)
  C1  <- sum(wg_cs * h1)
  C0  <- sum(wg_cs * (1 - h1))
  D1  <- sum(wg_cs * (h1 == 1)) + sum(wg_cs * (h2 == 1))
  D0  <- sum(wg_cs * (h1 == 0)) + sum(wg_cs * (h2 == 0))
  T11 <- sum(wg_cs * (h1 == 1 & h2 == 1)) +
         sum(wg_cs * (h2 == 1 & h3 == 1))
  T01 <- sum(wg_cs * (h1 == 0 & h2 == 1)) +
         sum(wg_cs * (h2 == 0 & h3 == 1))

  # H_mat and n_mis_mat already computed in the lq block above (reused here).
  M_count <- sum(wg * n_mis_mat)

  # eps stats
  Eps_num <- sum(wg * eps_num_mat)
  Eps_den <- sum(wg * eps_den_mat)

  # lambda_g stats
  Lg_count <- sum(wg * lambda_count_mat)
  Lg_xsum  <- sum(wg * lambda_xsum_mat)

  # --- Discrete lambda_d sufficient stats (mirror base / rho) ---
  cat_d_marginal_c_list <- vector("list", H)
  cat_d_marginal_w_list <- vector("list", H)
  cat_d_trans_curr_list <- vector("list", H)
  cat_d_trans_prev_list <- vector("list", H)
  cat_d_trans_w_list    <- vector("list", H)

  for (j in seq_len(H)) {
    wj <- wg[, j]

    # Marginal stats: every observation contributing a marginal interval-d emission.
    m1 <- (s1 == 0)
    ec <- c1[m1]; ew <- wj[m1]

    # Wave 2 marginals
    if (h1[j] == 0 && h2[j] == 0) {
      m <- (s2 == 0) & (s1 == 1)
      if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
    }
    if (h1[j] == 1 && h2[j] == 0) {
      m <- (s2 == 0)
      if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
    }
    if (h2[j] == 1) {
      m <- (s2 == 0)
      if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
    }

    # Wave 3 marginals
    if (h2[j] == 0 && h3[j] == 0) {
      m <- (s3 == 0) & (s2 == 1)
      if (any(m)) { ec <- c(ec, c3[m]); ew <- c(ew, wj[m]) }
    }
    if (h2[j] == 1 && h3[j] == 0) {
      m <- (s3 == 0)
      if (any(m)) { ec <- c(ec, c3[m]); ew <- c(ew, wj[m]) }
    }
    if (h3[j] == 1) {
      m <- (s3 == 0)
      if (any(m)) { ec <- c(ec, c3[m]); ew <- c(ew, wj[m]) }
    }

    # Aggregate marginal stats to ≤7 category bins (eliminates O(N) Brent FOC) [P1-2]
    if (length(ec) > 0) {
      wt <- tapply(ew, ec, sum)
      cat_d_marginal_c_list[[j]] <- as.integer(names(wt))
      cat_d_marginal_w_list[[j]] <- as.double(wt)
    } else {
      cat_d_marginal_c_list[[j]] <- integer(0)
      cat_d_marginal_w_list[[j]] <- double(0)
    }

    # Transition stats
    if (h1[j] == 0 && h2[j] == 0) {
      m <- (s2 == 0) & (s1 == 0)
      if (any(m)) {
        cat_d_trans_curr_list[[j]] <- c(cat_d_trans_curr_list[[j]], c2[m])
        cat_d_trans_prev_list[[j]] <- c(cat_d_trans_prev_list[[j]], c1[m])
        cat_d_trans_w_list[[j]]    <- c(cat_d_trans_w_list[[j]],    wj[m])
      }
    }
    if (h2[j] == 0 && h3[j] == 0) {
      m <- (s3 == 0) & (s2 == 0)
      if (any(m)) {
        cat_d_trans_curr_list[[j]] <- c(cat_d_trans_curr_list[[j]], c3[m])
        cat_d_trans_prev_list[[j]] <- c(cat_d_trans_prev_list[[j]], c2[m])
        cat_d_trans_w_list[[j]]    <- c(cat_d_trans_w_list[[j]],    wj[m])
      }
    }
    # Aggregate transition stats to ≤49 unique (curr,prev) pairs [P1-2]
    tw_j <- cat_d_trans_w_list[[j]]
    if (length(tw_j) > 0) {
      tc_j <- cat_d_trans_curr_list[[j]]; tp_j <- cat_d_trans_prev_list[[j]]
      key  <- tc_j + 7L * (tp_j - 1L)
      wt2  <- tapply(tw_j, key, sum)
      ks   <- as.integer(names(wt2))
      cat_d_trans_curr_list[[j]] <- (ks - 1L) %% 7L + 1L
      cat_d_trans_prev_list[[j]] <- (ks - 1L) %/% 7L + 1L
      cat_d_trans_w_list[[j]]    <- as.double(wt2)
    }
  }

  suff <- list(
    C1 = C1, C0 = C0,
    D1 = D1, D0 = D0,
    T11 = T11, T01 = T01,
    M = M_count,
    Eps_num = Eps_num, Eps_den = Eps_den,
    Lg_count = Lg_count, Lg_xsum = Lg_xsum,
    cat_d_marginal_c = unlist(cat_d_marginal_c_list, use.names = FALSE),
    cat_d_marginal_w = unlist(cat_d_marginal_w_list, use.names = FALSE),
    cat_d_trans_curr = unlist(cat_d_trans_curr_list, use.names = FALSE),
    cat_d_trans_prev = unlist(cat_d_trans_prev_list, use.names = FALSE),
    cat_d_trans_w    = unlist(cat_d_trans_w_list,    use.names = FALSE)
  )

  list(gamma = gamma_mat, loglik = ll, row_loglik = log_lik_i, suff = suff)
}
