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
validate_df_eps <- function(df) {
  cat_cols <- c("timegap_cat1", "timegap_cat2", "timegap_cat3")
  missing_cats <- setdiff(cat_cols, names(df))
  if (length(missing_cats) > 0)
    stop("e_step_eps requires columns: ", paste(missing_cats, collapse = ", "))
  na_timegap <- is.na(df$timegap_cat1) | is.na(df$timegap_cat2) | is.na(df$timegap_cat3)
  if (any(na_timegap))
    stop(sprintf("e_step_eps: %d obs have NA in timegap_cat columns.", sum(na_timegap)))
  bad_cats <- !all(df$timegap_cat1 %in% 1:7) || !all(df$timegap_cat2 %in% 1:7) ||
              !all(df$timegap_cat3 %in% 1:7)
  if (bad_cats) stop("e_step_eps: timegap_cat1/2/3 must contain only integers 1-7.")
  # Check y before tenure: NA in y would propagate through y == 1L comparisons.
  na_y <- is.na(df$y1) | is.na(df$y2) | is.na(df$y3)
  if (any(na_y))
    stop(sprintf("e_step_eps: %d obs have NA in y1/y2/y3.", sum(na_y)))
  bad_y <- !all(df$y1 %in% 0:1) || !all(df$y2 %in% 0:1) || !all(df$y3 %in% 0:1)
  if (bad_y) stop("e_step_eps: y1/y2/y3 must be binary (0 or 1).")
  na_tenure_emp <- (df$y1 == 1L & is.na(df$tenure1)) |
                   (df$y2 == 1L & is.na(df$tenure2)) |
                   (df$y3 == 1L & is.na(df$tenure3))
  if (any(na_tenure_emp))
    stop(sprintf("e_step_eps: %d obs have NA tenure at an employed wave.",
                 sum(na_tenure_emp)))
  bad_tenure <- (df$y1 == 1 & df$tenure1 <= 0) |
                (df$y2 == 1 & df$tenure2 <= 0) |
                (df$y3 == 1 & df$tenure3 <= 0)
  if (any(bad_tenure))
    stop(sprintf("e_step_eps: %d obs have non-positive tenure for employed state.",
                 sum(bad_tenure)))
  if (is.null(df$weight)) stop("e_step_eps: df must contain a 'weight' column.")
  if (any(is.na(df$weight)))
    stop(sprintf("e_step_eps: %d obs have NA weight.", sum(is.na(df$weight))))
  if (any(df$weight <= 0))
    stop(sprintf("e_step_eps: %d obs have non-positive weight.", sum(df$weight <= 0)))
  invisible(NULL)
}

.log_duration_history_prior_eps <- function(hmat, alpha, g_list, c_list,
                                             lambda_g, beta_g,
                                             lambda_d, beta_d) {
  N <- length(g_list[[1]])
  H <- nrow(hmat)
  out <- matrix(0, N, H)
  p_entry <- lapply(c_list[1:2], function(z)
    .duration_category_transition_probability(z, lambda_d, beta_d))
  for (j in seq_len(H)) {
    h <- hmat[j, ]
    out[, j] <- if (h[1] == 1L) log(alpha) else log1p(-alpha)
    for (t in 1:2) {
      if (h[t] == 1L) {
        p_change <- .duration_transition_probability(
          g_list[[t]], lambda_g, beta_g)
      } else {
        p_change <- p_entry[[t]]
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
  # --- Validate df inputs (skipped from the EM loop for efficiency) ---
  if (check_df) validate_df_eps(df)
  # --- Unpack parameters ---
  alpha    <- params$alpha
  theta0   <- params$theta0
  theta1   <- params$theta1
  pi_par   <- params$pi
  eps      <- params$eps
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d
  duration_dependent <- !is.null(params$beta_g) || !is.null(params$beta_d)
  beta_g <- if (is.null(params$beta_g)) 0 else params$beta_g
  beta_d <- if (is.null(params$beta_d)) 0 else params$beta_d

  if (!is.finite(eps) || eps <= 0 || eps >= 1) {
    stop(sprintf("e_step_eps: params$eps must be in (0, 1); got %.4g", eps))
  }
  if (!is.finite(lambda_g) || lambda_g <= 0) {
    stop(sprintf("e_step_eps: params$lambda_g must be > 0; got %.4g", lambda_g))
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
  if (!is.finite(lambda_d) || lambda_d <= 0) {
    stop(sprintf("e_step_eps: params$lambda_d must be > 0; got %.4g", lambda_d))
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
      hmat, alpha, g_list, c_list, lambda_g, beta_g, lambda_d, beta_d)
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

  # Pre-compute full N x 3 tenure matrix (s_full is already computed above).
  g_full <- cbind(g1, g2, g3)

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

      out <- log_emission_spell_g(g_mat, s_mat, t_offsets, lambda_g, eps,
                                  beta_g = beta_g)

      ld[, j]               <- ld[, j]               + out$loglik
      lambda_count_mat[, j] <- lambda_count_mat[, j] + out$lambda_count
      lambda_xsum_mat[, j]  <- lambda_xsum_mat[, j]  + out$lambda_xsum

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
          ld[mask, j] <- ld[mask, j] +
            .log_duration_density(g_t[mask], lambda_g, beta_g)
          lambda_count_mat[mask, j] <- lambda_count_mat[mask, j] + 1
          lambda_xsum_mat[mask, j]  <- lambda_xsum_mat[mask, j]  + g_t[mask]
        }
      }
    }

    # ---- (c) Timegap emissions: marginal + transition (mirror base model) ----
    # Wave 1: marginal interval whenever s_1 = 0.
    mask_w1 <- (s1 == 0L)
    if (any(mask_w1)) {
      ld[mask_w1, j] <- ld[mask_w1, j] +
        log_emission_interval_d(c1[mask_w1], lambda_d, beta_d)
    }
    # Waves 2 and 3: inlined timegap transitions (avoids 2× rep(0,N) alloc) [P2-9]
    # --- Transition (2 ← 1) ---
    {
      hp12 <- h_j[1L]; hc12 <- h_j[2L]
      if (hp12 == 0L && hc12 == 0L) {
        m12 <- (s2 == 0L) & (s1 == 0L)
        if (any(m12))  ld[m12, j] <- ld[m12, j] + log_emission_transition_d(c2[m12], c1[m12], lambda_d, beta_d)
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
        m23 <- (s3 == 0L) & (s2 == 0L)
        if (any(m23))  ld[m23, j] <- ld[m23, j] + log_emission_transition_d(c3[m23], c2[m23], lambda_d, beta_d)
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
    ll        <- sum(wi * log_denom)
  }

  if (!suff_stats) {
    return(list(gamma = gamma_mat, loglik = ll, suff = NULL))
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

  list(gamma = gamma_mat, loglik = ll, suff = suff)
}
