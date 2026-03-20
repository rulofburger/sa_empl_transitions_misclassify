# ==============================================================================
# EM-tenure: E-step — responsibilities and sufficient statistics (vectorised)
# ==============================================================================
# Computes gamma_{ih} = P(H=h | y_i, Theta^old) for each observation i and
# latent history h, plus all weighted sufficient statistics needed by the M-step.
#
# KEY DESIGN DECISIONS (all from TeX):
# 1. Wave 1 matched durations use EMG (left-censored spells).
# 2. Continuations use increment likelihood ONLY when the previous wave is
#    also correctly observed (s_{t-1} matches h_{t-1}), so the lag is available.
# 3. When h_{t-1}=h_t but s_{t-1}!=h_{t-1} (continuation with previous
#    misclassified), the increment is undefined — use EMG instead.
# 4. Within-panel starts (h_{t-1}!=h_t, s_t matches h_t) use Normal(0.25, sigma^2).
# 5. Misclassified durations (s_t != h_t) use EMG with the "inappropriate" clock.
# 6. Sufficient stats for sigma^2 come from (a) increments and (b) starts,
#    pooled in the M-step. Both require s_{t-1} conditioning.
#
# VECTORISATION: All operations are fully vectorised over observations. Each
# observation i has its own (s, g, d) but shares the same hmat and parameters.
# We work with N x H matrices where N = #observations, H = 8 latent histories.
# The only explicit loop is over the 8 histories (fixed size).
#
# TeX reference: Section 2.7 (E-step), Eqs (8)-(18)
# ==============================================================================

# --- Internal helper: vectorised wave t >= 2 emission log-densities ----------
#
# Returns an N x H matrix of log-density contributions for wave t.
# All 8 emission cases from the TeX are handled.
#
# @param s_t,s_prev Length-N vectors: observed state at wave t and t-1.
# @param g_t,g_prev Length-N vectors: observed tenure at wave t and t-1.
# @param d_t,d_prev Length-N vectors: observed timegap (continuous, years) at
#   wave t and t-1. Used only when discrete_timegap = FALSE.
# @param cat_t,cat_prev Length-N integer vectors: timegap category codes (1-7)
#   at wave t and t-1. Used only when discrete_timegap = TRUE.
# @param h_prev,h_curr Length-H vectors: latent states at t-1 and t.
# @param lambda_g,lambda_d Scalar exponential rates.
# @param sigma2_g,sigma2_d Scalar measurement variances (sigma2_d used only
#   when discrete_timegap = FALSE).
# @param discrete_timegap Logical; if TRUE use interval-censored discrete
#   emissions for d; if FALSE use legacy continuous EMG emissions.
# @return N x H matrix of log-density contributions.
# @keywords internal
.wave_emission_vec <- function(s_t, s_prev, g_t, g_prev, d_t, d_prev,
                               cat_t, cat_prev,
                               h_prev, h_curr,
                               lambda_g, lambda_d, sigma2_g, sigma2_d,
                               discrete_timegap = TRUE) {
  N <- length(s_t)
  H <- length(h_prev)
  ld <- matrix(0, nrow = N, ncol = H)

  # Pre-compute per-observation scalars (N-vectors)
  delta_g <- g_t - g_prev - .QUARTER_YEARS
  if (!discrete_timegap) {
    delta_d <- d_t - d_prev - .QUARTER_YEARS
  }

  # Loop over histories only (H = 8, fixed)
  for (j in seq_len(H)) {
    hp <- h_prev[j]
    hc <- h_curr[j]

    # --- CORRECTLY CLASSIFIED at wave t ---
    if (hp == 1 && hc == 1) {
      # Employment continuation, previous observed: s_t=1, s_prev=1
      mask_obs <- (s_t == 1) & (s_prev == 1)
      if (any(mask_obs)) {
        ld[mask_obs, j] <- ld[mask_obs, j] +
          log_emission_increment_g(delta_g[mask_obs], sigma2_g)
      }
      # Employment continuation, previous misclassified: s_t=1, s_prev=0
      mask_mis <- (s_t == 1) & (s_prev == 0)
      if (any(mask_mis)) {
        ld[mask_mis, j] <- ld[mask_mis, j] +
          log_emg(g_t[mask_mis], lambda_g, sigma2_g)
      }
    } else if (hp == 0 && hc == 1) {
      # Within-panel employment start: s_t=1
      mask <- (s_t == 1)
      if (any(mask)) {
        ld[mask, j] <- ld[mask, j] +
          log_emission_start_g(g_t[mask], sigma2_g)
      }
    } else if (hp == 0 && hc == 0) {
      if (discrete_timegap) {
        # Case 3 (discrete): interval-censored transition emission
        # Nonemployment continuation, previous observed: s_t=0, s_prev=0
        mask_obs <- (s_t == 0) & (s_prev == 0)
        if (any(mask_obs)) {
          ld[mask_obs, j] <- ld[mask_obs, j] +
            log_emission_transition_d(cat_t[mask_obs], cat_prev[mask_obs], lambda_d)
        }
        # Nonemployment continuation, previous misclassified: s_t=0, s_prev=1
        # (Case 2: interval-censored marginal, since lag is unavailable)
        mask_mis <- (s_t == 0) & (s_prev == 1)
        if (any(mask_mis)) {
          ld[mask_mis, j] <- ld[mask_mis, j] +
            log_emission_interval_d(cat_t[mask_mis], lambda_d)
        }
      } else {
        # Legacy continuous EMG
        mask_obs <- (s_t == 0) & (s_prev == 0)
        if (any(mask_obs)) {
          ld[mask_obs, j] <- ld[mask_obs, j] +
            log_emission_increment_d(delta_d[mask_obs], sigma2_d)
        }
        mask_mis <- (s_t == 0) & (s_prev == 1)
        if (any(mask_mis)) {
          ld[mask_mis, j] <- ld[mask_mis, j] +
            log_emg(d_t[mask_mis], lambda_d, sigma2_d)
        }
      }
    } else if (hp == 1 && hc == 0) {
      # Within-panel nonemployment start: s_t=0
      mask <- (s_t == 0)
      if (any(mask)) {
        if (discrete_timegap) {
          # Case 4: interval-censored marginal (fresh start, no lag)
          ld[mask, j] <- ld[mask, j] +
            log_emission_interval_d(cat_t[mask], lambda_d)
        } else {
          ld[mask, j] <- ld[mask, j] +
            log_emission_start_d(d_t[mask], sigma2_d)
        }
      }
    }

    # --- MISCLASSIFIED at wave t ---
    if (hc == 0) {
      # Latent: nonemployed; observed: employed (s_t=1). Use tenure clock.
      mask <- (s_t == 1)
      if (any(mask)) {
        ld[mask, j] <- ld[mask, j] +
          log_emg(g_t[mask], lambda_g, sigma2_g)
      }
    } else {
      # Latent: employed; observed: nonemployed (s_t=0). Use timegap clock.
      mask <- (s_t == 0)
      if (any(mask)) {
        if (discrete_timegap) {
          # Case 6: interval-censored marginal (misclassified nonemployed)
          ld[mask, j] <- ld[mask, j] +
            log_emission_interval_d(cat_t[mask], lambda_d)
        } else {
          ld[mask, j] <- ld[mask, j] +
            log_emg(d_t[mask], lambda_d, sigma2_d)
        }
      }
    }
  }

  return(ld)
}


#' E-step: compute responsibilities and sufficient statistics (vectorised)
#'
#' Fully vectorised over observations: all N individuals are processed
#' simultaneously via matrix operations. The only explicit loop is over
#' the 8 latent histories (fixed size).
#'
#' @param df Data frame with columns: y1, y2, y3 (observed employment 0/1),
#'   tenure1, tenure2, tenure3 (observed tenure in years),
#'   timegap1, timegap2, timegap3 (observed nonemployment duration in years,
#'     used only when discrete_timegap = FALSE),
#'   timegap_cat1, timegap_cat2, timegap_cat3 (integer category codes 1-7,
#'     used only when discrete_timegap = TRUE),
#'   weight (survey weights).
#' @param params Named list: alpha, theta0, theta1, pi, sigma2_g, lambda_g,
#'   lambda_d. When discrete_timegap = FALSE, must also include sigma2_d.
#'   When discrete_timegap = TRUE, sigma2_d is not used and may be omitted.
#' @param discrete_timegap Logical (default TRUE). If TRUE, use interval-
#'   censored discrete Exp(lambda_d) emissions for nonemployment durations
#'   (columns timegap_cat1-3 must be present and in 1:7). If FALSE, use
#'   legacy continuous EMG(lambda_d, sigma2_d) emissions (columns timegap1-3
#'   must be present).
#' @return List with:
#'   - gamma: N x 8 matrix of responsibilities
#'   - loglik: weighted observed-data log-likelihood
#'   - suff: list of sufficient statistics for M-step, including:
#'       emg_g_x, emg_g_w: duration values and weights for all EMG-lambda_g obs
#'       When discrete_timegap = FALSE:
#'         emg_d_x, emg_d_w: duration values and weights for EMG-lambda_d obs
#'       When discrete_timegap = TRUE:
#'         cat_d_marginal_c: integer vector of category codes for marginal d obs
#'         cat_d_marginal_w: corresponding weights
#'         cat_d_trans_curr: integer vector of current-wave categories (Case 3)
#'         cat_d_trans_prev: integer vector of previous-wave categories (Case 3)
#'         cat_d_trans_w: corresponding weights
#' @export
e_step <- function(df, params, discrete_timegap = TRUE) {
  # --- Unpack parameters ---
  alpha    <- params$alpha
  theta0   <- params$theta0
  theta1   <- params$theta1
  pi_par   <- params$pi
  sigma2_g <- params$sigma2_g
  sigma2_d <- if (discrete_timegap) NA_real_ else params$sigma2_d
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d

  # --- Validate required columns ---
  if (discrete_timegap) {
    cat_cols <- c("timegap_cat1", "timegap_cat2", "timegap_cat3")
    missing_cats <- setdiff(cat_cols, names(df))
    if (length(missing_cats) > 0) {
      stop("e_step discrete_timegap=TRUE requires columns: ",
           paste(missing_cats, collapse = ", "))
    }
    bad_cats <- !all(df$timegap_cat1 %in% 1:7, na.rm = TRUE) ||
                !all(df$timegap_cat2 %in% 1:7, na.rm = TRUE) ||
                !all(df$timegap_cat3 %in% 1:7, na.rm = TRUE)
    if (bad_cats) {
      stop("e_step: timegap_cat1/2/3 must contain only integers 1-7 (no NA, no 0/8/99).")
    }
  }

  # --- Latent structure (H = 8) ---
  hmat    <- latent_histories()
  prior_h <- prior_over_histories(hmat, theta1, theta0, alpha)
  log_prior <- log(.bound01(prior_h))    # length-8

  N <- nrow(df)
  H <- nrow(hmat)
  h1 <- hmat[, 1]; h2 <- hmat[, 2]; h3 <- hmat[, 3]

  # --- Validate no NA durations ---
  # NA in tenure or timegap would propagate silently through log-density
  # computations, producing NaN in the log-likelihood.
  na_tenure <- is.na(df$tenure1) | is.na(df$tenure2) | is.na(df$tenure3)
  if (discrete_timegap) {
    na_timegap <- is.na(df$timegap_cat1) | is.na(df$timegap_cat2) | is.na(df$timegap_cat3)
  } else {
    na_timegap <- is.na(df$timegap1) | is.na(df$timegap2) | is.na(df$timegap3)
  }
  n_na <- sum(na_tenure | na_timegap)
  if (n_na > 0) {
    stop(sprintf(
      "E-step: %d observation(s) have NA in tenure or timegap columns. Remove these rows first.",
      n_na
    ))
  }

  # --- Validate duration positivity ---
  # When observed as employed (y=1), tenure must be > 0; when nonemployed,
  # timegap must be > 0 (continuous) or timegap_cat must be in 1:7 (discrete).
  if (discrete_timegap) {
    bad <- (df$y1 == 1 & df$tenure1 <= 0) |
           (df$y2 == 1 & df$tenure2 <= 0) |
           (df$y3 == 1 & df$tenure3 <= 0)
    # NB: timegap_cat validity already checked above
  } else {
    bad <- (df$y1 == 1 & df$tenure1 <= 0) | (df$y1 == 0 & df$timegap1 <= 0) |
           (df$y2 == 1 & df$tenure2 <= 0) | (df$y2 == 0 & df$timegap2 <= 0) |
           (df$y3 == 1 & df$tenure3 <= 0) | (df$y3 == 0 & df$timegap3 <= 0)
  }
  if (any(bad)) {
    stop(sprintf(
      paste0("E-step: %d observation(s) have non-positive duration for their ",
             "observed state. Filter these rows before calling e_step()."),
      sum(bad)
    ))
  }

  # --- Extract data as N-vectors ---
  s1 <- df$y1;       s2 <- df$y2;       s3 <- df$y3
  g1 <- df$tenure1;  g2 <- df$tenure2;  g3 <- df$tenure3
  wi <- df$weight

  if (discrete_timegap) {
    c1 <- df$timegap_cat1; c2 <- df$timegap_cat2; c3 <- df$timegap_cat3
    d1 <- d2 <- d3 <- NULL  # not used in discrete mode
  } else {
    d1 <- df$timegap1; d2 <- df$timegap2; d3 <- df$timegap3
    c1 <- c2 <- c3 <- NULL  # not used in continuous mode
  }

  # --- Misclassification log-probability: N x H matrix ---
  # lq[i, j] = sum_t log P(s_it | h_jt, pi)
  if (pi_par == 0) {
    lq <- matrix(-Inf, nrow = N, ncol = H)
    for (j in seq_len(H)) {
      match_all <- (s1 == h1[j]) & (s2 == h2[j]) & (s3 == h3[j])
      lq[match_all, j] <- 0
    }
  } else {
    pi_b <- .bound01(pi_par)
    lp_match    <- log1p(-pi_b)
    lp_mismatch <- log(pi_b)
    lq <- matrix(0, nrow = N, ncol = H)
    for (j in seq_len(H)) {
      lq[, j] <- ifelse(s1 == h1[j], lp_match, lp_mismatch) +
                  ifelse(s2 == h2[j], lp_match, lp_mismatch) +
                  ifelse(s3 == h3[j], lp_match, lp_mismatch)
    }
  }

  # --- Duration emission log-densities: N x H matrix ---
  ld <- matrix(0, nrow = N, ncol = H)

  # ============ WAVE 1 (left-censored: all durations use EMG / cat marginal) ============
  for (j in seq_len(H)) {
    if (h1[j] == 1) {
      # h1=1: any correctly-classified or misclassified scenario at wave 1
      # Observed emp (s1=1): EMG on g1 (tenure clock)
      mask_emp <- (s1 == 1)
      if (any(mask_emp)) {
        ld[mask_emp, j] <- ld[mask_emp, j] +
          log_emg(g1[mask_emp], lambda_g, sigma2_g)
      }
      # Observed nonemp (s1=0): timegap clock (Case 1: interval marginal or EMG)
      mask_non <- (s1 == 0)
      if (any(mask_non)) {
        if (discrete_timegap) {
          # Case 1: interval-censored marginal (wave 1 = left-censored)
          ld[mask_non, j] <- ld[mask_non, j] +
            log_emission_interval_d(c1[mask_non], lambda_d)
        } else {
          ld[mask_non, j] <- ld[mask_non, j] +
            log_emg(d1[mask_non], lambda_d, sigma2_d)
        }
      }
    } else {
      # h1=0: observed nonemp (s1=0): EMG / interval on d1
      mask_non <- (s1 == 0)
      if (any(mask_non)) {
        if (discrete_timegap) {
          ld[mask_non, j] <- ld[mask_non, j] +
            log_emission_interval_d(c1[mask_non], lambda_d)
        } else {
          ld[mask_non, j] <- ld[mask_non, j] +
            log_emg(d1[mask_non], lambda_d, sigma2_d)
        }
      }
      # Observed emp (s1=1): EMG on g1 (tenure clock)
      mask_emp <- (s1 == 1)
      if (any(mask_emp)) {
        ld[mask_emp, j] <- ld[mask_emp, j] +
          log_emg(g1[mask_emp], lambda_g, sigma2_g)
      }
    }
  }

  # ============ WAVES 2-3 ============
  ld <- ld + .wave_emission_vec(
    s_t = s2, s_prev = s1, g_t = g2, g_prev = g1, d_t = d2, d_prev = d1,
    cat_t = c2, cat_prev = c1,
    h_prev = h1, h_curr = h2,
    lambda_g = lambda_g, lambda_d = lambda_d,
    sigma2_g = sigma2_g, sigma2_d = sigma2_d,
    discrete_timegap = discrete_timegap
  )
  ld <- ld + .wave_emission_vec(
    s_t = s3, s_prev = s2, g_t = g3, g_prev = g2, d_t = d3, d_prev = d2,
    cat_t = c3, cat_prev = c2,
    h_prev = h2, h_curr = h3,
    lambda_g = lambda_g, lambda_d = lambda_d,
    sigma2_g = sigma2_g, sigma2_d = sigma2_d,
    discrete_timegap = discrete_timegap
  )

  # --- Posterior responsibilities: N x H ---
  log_kernel <- sweep(lq + ld, 2, log_prior, "+")

  # Row-wise log-sum-exp for denominator (numerically stable)
  row_max <- apply(log_kernel, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_kernel - row_max)))

  gamma_mat <- exp(log_kernel - log_denom)

  # Weighted log-likelihood
  ll <- sum(wi * log_denom)

  # --- Sufficient statistics (vectorised over observations) ---
  wg <- gamma_mat * wi  # N x H

  # Column sums of weighted responsibilities: length-H
  wg_cs <- colSums(wg)

  # Initial state (Eq 14)
  C1 <- sum(wg_cs * h1)
  C0 <- sum(wg_cs * (1 - h1))

  # Transition counts (Eq 14)
  D1  <- sum(wg_cs * (h1 == 1)) + sum(wg_cs * (h2 == 1))
  D0  <- sum(wg_cs * (h1 == 0)) + sum(wg_cs * (h2 == 0))
  T11 <- sum(wg_cs * (h1 == 1 & h2 == 1)) +
         sum(wg_cs * (h2 == 1 & h3 == 1))
  T01 <- sum(wg_cs * (h1 == 0 & h2 == 1)) +
         sum(wg_cs * (h2 == 0 & h3 == 1))

  # Misclassification count (Eq 14)
  M_count <- 0
  for (j in seq_len(H)) {
    n_mis <- (h1[j] != s1) + (h2[j] != s2) + (h3[j] != s3)  # N-vector
    M_count <- M_count + sum(wg[, j] * n_mis)
  }

  # --- Variance sufficient stats (Eqs 15-18, tenure only) ---
  # sigma2_d stats are only needed in the legacy continuous mode.
  Sg <- 0; Ng <- 0
  Sg_start <- 0; Ng_start <- 0

  if (!discrete_timegap) {
    Sd <- 0; Nd <- 0
    Sd_start <- 0; Nd_start <- 0
  }

  for (j in seq_len(H)) {
    wj <- wg[, j]

    # ---- t = 2 ----
    if (h1[j] == 1 && h2[j] == 1) {
      mask <- (s2 == 1) & (s1 == 1)
      if (any(mask)) {
        dg <- g2[mask] - g1[mask] - .QUARTER_YEARS
        Sg <- Sg + sum(wj[mask] * dg^2)
        Ng <- Ng + sum(wj[mask])
      }
    }
    if (!discrete_timegap) {
      if (h1[j] == 0 && h2[j] == 0) {
        mask <- (s2 == 0) & (s1 == 0)
        if (any(mask)) {
          dd <- d2[mask] - d1[mask] - .QUARTER_YEARS
          Sd <- Sd + sum(wj[mask] * dd^2)
          Nd <- Nd + sum(wj[mask])
        }
      }
    }
    if (h1[j] == 0 && h2[j] == 1) {
      mask <- (s2 == 1)
      if (any(mask)) {
        Sg_start <- Sg_start + sum(wj[mask] * (g2[mask] - .QUARTER_YEARS)^2)
        Ng_start <- Ng_start + sum(wj[mask])
      }
    }
    if (!discrete_timegap) {
      if (h1[j] == 1 && h2[j] == 0) {
        mask <- (s2 == 0)
        if (any(mask)) {
          Sd_start <- Sd_start + sum(wj[mask] * (d2[mask] - .QUARTER_YEARS)^2)
          Nd_start <- Nd_start + sum(wj[mask])
        }
      }
    }

    # ---- t = 3 ----
    if (h2[j] == 1 && h3[j] == 1) {
      mask <- (s3 == 1) & (s2 == 1)
      if (any(mask)) {
        dg <- g3[mask] - g2[mask] - .QUARTER_YEARS
        Sg <- Sg + sum(wj[mask] * dg^2)
        Ng <- Ng + sum(wj[mask])
      }
    }
    if (!discrete_timegap) {
      if (h2[j] == 0 && h3[j] == 0) {
        mask <- (s3 == 0) & (s2 == 0)
        if (any(mask)) {
          dd <- d3[mask] - d2[mask] - .QUARTER_YEARS
          Sd <- Sd + sum(wj[mask] * dd^2)
          Nd <- Nd + sum(wj[mask])
        }
      }
    }
    if (h2[j] == 0 && h3[j] == 1) {
      mask <- (s3 == 1)
      if (any(mask)) {
        Sg_start <- Sg_start + sum(wj[mask] * (g3[mask] - .QUARTER_YEARS)^2)
        Ng_start <- Ng_start + sum(wj[mask])
      }
    }
    if (!discrete_timegap) {
      if (h2[j] == 1 && h3[j] == 0) {
        mask <- (s3 == 0)
        if (any(mask)) {
          Sd_start <- Sd_start + sum(wj[mask] * (d3[mask] - .QUARTER_YEARS)^2)
          Nd_start <- Nd_start + sum(wj[mask])
        }
      }
    }
  }

  # --- EMG gradient data for M-step ---
  # lambda_g: unchanged in both modes (always EMG)
  # lambda_d (discrete): gradient FOC uses (cat, w) pairs for all d emissions
  # lambda_d (continuous): gradient FOC uses (x, w) pairs for EMG emissions

  # --- EMG-lambda_g sufficient stats (same in both modes) ---
  # Observation types (TeX Section 2.7):
  #   (A+B) Wave 1: s1=1 for ANY h1 (matched or misclassified).
  #   (C) t>=2: hp=1, hc=1, s_t=1, s_{t-1}=0  (continuation, prev misclassified)
  #   (D) t>=2: hc=0, s_t=1                     (misclassified as employed)
  emg_g_x_list <- vector("list", H)
  emg_g_w_list <- vector("list", H)

  if (discrete_timegap) {
    # Discrete mode: collect (cat, w) for all d-emissions (Cases 1–6)
    # Case 1 (wave 1), Case 2 (prev miscl.), Case 4 (within-panel start),
    # Case 6 (misclassified nonemployed) — all use log_emission_interval_d
    # Case 3 (continuation, prev observed) — uses log_emission_transition_d
    cat_d_marginal_c_list <- vector("list", H)
    cat_d_marginal_w_list <- vector("list", H)
    cat_d_trans_curr_list <- vector("list", H)
    cat_d_trans_prev_list <- vector("list", H)
    cat_d_trans_w_list    <- vector("list", H)
  } else {
    emg_d_x_list <- vector("list", H)
    emg_d_w_list <- vector("list", H)
  }

  for (j in seq_len(H)) {
    wj <- wg[, j]

    # ---- EMG-lambda_g ----
    mask <- (s1 == 1)
    ex_g <- g1[mask];  ew_g <- wj[mask]
    if (h1[j] == 1 && h2[j] == 1) {
      m <- (s2 == 1) & (s1 == 0)
      if (any(m)) { ex_g <- c(ex_g, g2[m]); ew_g <- c(ew_g, wj[m]) }
    }
    if (h2[j] == 0) {
      m <- (s2 == 1)
      if (any(m)) { ex_g <- c(ex_g, g2[m]); ew_g <- c(ew_g, wj[m]) }
    }
    if (h2[j] == 1 && h3[j] == 1) {
      m <- (s3 == 1) & (s2 == 0)
      if (any(m)) { ex_g <- c(ex_g, g3[m]); ew_g <- c(ew_g, wj[m]) }
    }
    if (h3[j] == 0) {
      m <- (s3 == 1)
      if (any(m)) { ex_g <- c(ex_g, g3[m]); ew_g <- c(ew_g, wj[m]) }
    }
    emg_g_x_list[[j]] <- ex_g
    emg_g_w_list[[j]] <- ew_g

    if (discrete_timegap) {
      # ---- Discrete lambda_d sufficient stats ----
      # Marginal emissions (Cases 1, 2, 4, 6):
      #   Wave 1: s1=0 for any h1 (Case 1)
      #   t>=2: hp=0, hc=0, s_t=0, s_{t-1}=1 (Case 2: prev miscl.)
      #   t>=2: hp=1, hc=0, s_t=0 (Case 4: within-panel start)
      #   t>=2: hc=1, s_t=0 (Case 6: misclassified nonemployed)
      m1 <- (s1 == 0)
      ec <- c1[m1]; ew <- wj[m1]
      # Case 2: hp=0, hc=0, s2=0, s1=1
      if (h1[j] == 0 && h2[j] == 0) {
        m <- (s2 == 0) & (s1 == 1)
        if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
      }
      # Case 4: hp=1, hc=0, s2=0
      if (h1[j] == 1 && h2[j] == 0) {
        m <- (s2 == 0)
        if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
      }
      # Case 6: hc=1, s2=0
      if (h2[j] == 1) {
        m <- (s2 == 0)
        if (any(m)) { ec <- c(ec, c2[m]); ew <- c(ew, wj[m]) }
      }
      # Same for wave 3
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
      cat_d_marginal_c_list[[j]] <- ec
      cat_d_marginal_w_list[[j]] <- ew

      # Transition emissions (Case 3): hp=0, hc=0, s_t=0, s_{t-1}=0
      if (h1[j] == 0 && h2[j] == 0) {
        m <- (s2 == 0) & (s1 == 0)
        cat_d_trans_curr_list[[j]] <- c(cat_d_trans_curr_list[[j]], c2[m])
        cat_d_trans_prev_list[[j]] <- c(cat_d_trans_prev_list[[j]], c1[m])
        cat_d_trans_w_list[[j]]    <- c(cat_d_trans_w_list[[j]],    wj[m])
      }
      if (h2[j] == 0 && h3[j] == 0) {
        m <- (s3 == 0) & (s2 == 0)
        cat_d_trans_curr_list[[j]] <- c(cat_d_trans_curr_list[[j]], c3[m])
        cat_d_trans_prev_list[[j]] <- c(cat_d_trans_prev_list[[j]], c2[m])
        cat_d_trans_w_list[[j]]    <- c(cat_d_trans_w_list[[j]],    wj[m])
      }

    } else {
      # ---- Legacy continuous EMG-lambda_d sufficient stats ----
      # (A+B) Wave 1: s1=0 for any h1
      mask <- (s1 == 0)
      ex_d <- d1[mask];  ew_d <- wj[mask]
      if (h1[j] == 0 && h2[j] == 0) {
        m <- (s2 == 0) & (s1 == 1)
        if (any(m)) { ex_d <- c(ex_d, d2[m]); ew_d <- c(ew_d, wj[m]) }
      }
      if (h2[j] == 1) {
        m <- (s2 == 0)
        if (any(m)) { ex_d <- c(ex_d, d2[m]); ew_d <- c(ew_d, wj[m]) }
      }
      if (h2[j] == 0 && h3[j] == 0) {
        m <- (s3 == 0) & (s2 == 1)
        if (any(m)) { ex_d <- c(ex_d, d3[m]); ew_d <- c(ew_d, wj[m]) }
      }
      if (h3[j] == 1) {
        m <- (s3 == 0)
        if (any(m)) { ex_d <- c(ex_d, d3[m]); ew_d <- c(ew_d, wj[m]) }
      }
      emg_d_x_list[[j]] <- ex_d
      emg_d_w_list[[j]] <- ew_d
    }
  }

  emg_g_x <- unlist(emg_g_x_list, use.names = FALSE)
  emg_g_w <- unlist(emg_g_w_list, use.names = FALSE)

  # --- Assemble sufficient stats list ---
  suff <- list(
    C1 = C1, C0 = C0,
    D1 = D1, D0 = D0,
    T11 = T11, T01 = T01,
    M = M_count,
    Sg = Sg, Ng = Ng,
    Sg_start = Sg_start, Ng_start = Ng_start,
    emg_g_x = emg_g_x, emg_g_w = emg_g_w
  )

  if (discrete_timegap) {
    suff$cat_d_marginal_c <- unlist(cat_d_marginal_c_list, use.names = FALSE)
    suff$cat_d_marginal_w <- unlist(cat_d_marginal_w_list, use.names = FALSE)
    suff$cat_d_trans_curr <- unlist(cat_d_trans_curr_list, use.names = FALSE)
    suff$cat_d_trans_prev <- unlist(cat_d_trans_prev_list, use.names = FALSE)
    suff$cat_d_trans_w    <- unlist(cat_d_trans_w_list,    use.names = FALSE)
  } else {
    suff$Sd       <- Sd;       suff$Nd       <- Nd
    suff$Sd_start <- Sd_start; suff$Nd_start <- Nd_start
    suff$emg_d_x  <- unlist(emg_d_x_list, use.names = FALSE)
    suff$emg_d_w  <- unlist(emg_d_w_list, use.names = FALSE)
  }

  return(list(gamma = gamma_mat, loglik = ll, suff = suff))
}
