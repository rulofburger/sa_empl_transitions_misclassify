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
# @param d_t,d_prev Length-N vectors: observed timegap at wave t and t-1.
# @param h_prev,h_curr Length-H vectors: latent states at t-1 and t.
# @param lambda_g,lambda_d Scalar exponential rates.
# @param sigma2_g,sigma2_d Scalar measurement variances.
# @return N x H matrix of log-density contributions.
# @keywords internal
.wave_emission_vec <- function(s_t, s_prev, g_t, g_prev, d_t, d_prev,
                               h_prev, h_curr,
                               lambda_g, lambda_d, sigma2_g, sigma2_d) {
  N <- length(s_t)
  H <- length(h_prev)
  ld <- matrix(0, nrow = N, ncol = H)

  # Pre-compute per-observation scalars (N-vectors)
  delta_g <- g_t - g_prev - .QUARTER_YEARS
  delta_d <- d_t - d_prev - .QUARTER_YEARS

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
      # Nonemployment continuation, previous observed: s_t=0, s_prev=0
      mask_obs <- (s_t == 0) & (s_prev == 0)
      if (any(mask_obs)) {
        ld[mask_obs, j] <- ld[mask_obs, j] +
          log_emission_increment_d(delta_d[mask_obs], sigma2_d)
      }
      # Nonemployment continuation, previous misclassified: s_t=0, s_prev=1
      mask_mis <- (s_t == 0) & (s_prev == 1)
      if (any(mask_mis)) {
        ld[mask_mis, j] <- ld[mask_mis, j] +
          log_emg(d_t[mask_mis], lambda_d, sigma2_d)
      }
    } else if (hp == 1 && hc == 0) {
      # Within-panel nonemployment start: s_t=0
      mask <- (s_t == 0)
      if (any(mask)) {
        ld[mask, j] <- ld[mask, j] +
          log_emission_start_d(d_t[mask], sigma2_d)
      }
    }

    # --- MISCLASSIFIED at wave t ---
    if (hc == 0) {
      mask <- (s_t == 1)
      if (any(mask)) {
        ld[mask, j] <- ld[mask, j] +
          log_emg(g_t[mask], lambda_g, sigma2_g)
      }
    } else {
      mask <- (s_t == 0)
      if (any(mask)) {
        ld[mask, j] <- ld[mask, j] +
          log_emg(d_t[mask], lambda_d, sigma2_d)
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
#'   timegap1, timegap2, timegap3 (observed nonemployment duration in years),
#'   weight (survey weights).
#' @param params Named list: alpha, theta0, theta1, pi, sigma2_g, sigma2_d,
#'   lambda_g, lambda_d.
#' @return List with:
#'   - gamma: N x 8 matrix of responsibilities
#'   - loglik: weighted observed-data log-likelihood
#'   - suff: list of sufficient statistics for M-step, including:
#'       emg_g_x, emg_g_w: duration values and weights for all EMG-lambda_g obs
#'       emg_d_x, emg_d_w: duration values and weights for all EMG-lambda_d obs
#' @export
e_step <- function(df, params) {
  # --- Unpack parameters ---
  alpha    <- params$alpha
  theta0   <- params$theta0
  theta1   <- params$theta1
  pi_par   <- params$pi
  sigma2_g <- params$sigma2_g
  sigma2_d <- params$sigma2_d
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d

  # --- Latent structure (H = 8) ---
  hmat    <- latent_histories()
  prior_h <- prior_over_histories(hmat, theta1, theta0, alpha)
  log_prior <- log(.bound01(prior_h))    # length-8

  N <- nrow(df)
  H <- nrow(hmat)
  h1 <- hmat[, 1]; h2 <- hmat[, 2]; h3 <- hmat[, 3]

  # --- Validate duration positivity ---
  # When observed as employed (y=1), tenure must be > 0; when nonemployed,
  # timegap must be > 0. Zero or negative durations cause log_emg -> -Inf for
  # ALL 8 histories, producing NaN responsibilities.
  bad <- (df$y1 == 1 & df$tenure1 <= 0) | (df$y1 == 0 & df$timegap1 <= 0) |
         (df$y2 == 1 & df$tenure2 <= 0) | (df$y2 == 0 & df$timegap2 <= 0) |
         (df$y3 == 1 & df$tenure3 <= 0) | (df$y3 == 0 & df$timegap3 <= 0)
  if (any(bad, na.rm = TRUE)) {
    stop(sprintf(
      paste0("E-step: %d observation(s) have non-positive duration for their ",
             "observed state. Filter these rows before calling e_step()."),
      sum(bad, na.rm = TRUE)
    ))
  }

  # --- Extract data as N-vectors ---
  s1 <- df$y1;       s2 <- df$y2;       s3 <- df$y3
  g1 <- df$tenure1;  g2 <- df$tenure2;  g3 <- df$tenure3
  d1 <- df$timegap1; d2 <- df$timegap2; d3 <- df$timegap3
  wi <- df$weight

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

  # ============ WAVE 1 (left-censored: all durations use EMG) ============
  for (j in seq_len(H)) {
    if (h1[j] == 1) {
      # h1=1: observed emp (s1=1) -> EMG on g; observed nonemp (s1=0) -> EMG on d
      mask_emp <- (s1 == 1)
      if (any(mask_emp)) {
        ld[mask_emp, j] <- ld[mask_emp, j] +
          log_emg(g1[mask_emp], lambda_g, sigma2_g)
      }
      mask_non <- (s1 == 0)
      if (any(mask_non)) {
        ld[mask_non, j] <- ld[mask_non, j] +
          log_emg(d1[mask_non], lambda_d, sigma2_d)
      }
    } else {
      # h1=0: observed nonemp (s1=0) -> EMG on d; observed emp (s1=1) -> EMG on g
      mask_non <- (s1 == 0)
      if (any(mask_non)) {
        ld[mask_non, j] <- ld[mask_non, j] +
          log_emg(d1[mask_non], lambda_d, sigma2_d)
      }
      mask_emp <- (s1 == 1)
      if (any(mask_emp)) {
        ld[mask_emp, j] <- ld[mask_emp, j] +
          log_emg(g1[mask_emp], lambda_g, sigma2_g)
      }
    }
  }

  # ============ WAVES 2-3 ============
  ld <- ld + .wave_emission_vec(s2, s1, g2, g1, d2, d1, h1, h2,
                                lambda_g, lambda_d, sigma2_g, sigma2_d)
  ld <- ld + .wave_emission_vec(s3, s2, g3, g2, d3, d2, h2, h3,
                                lambda_g, lambda_d, sigma2_g, sigma2_d)

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

  # --- Variance sufficient stats (Eqs 15-18) ---
  Sg <- 0; Ng <- 0; Sd <- 0; Nd <- 0
  Sg_start <- 0; Ng_start <- 0; Sd_start <- 0; Nd_start <- 0

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
    if (h1[j] == 0 && h2[j] == 0) {
      mask <- (s2 == 0) & (s1 == 0)
      if (any(mask)) {
        dd <- d2[mask] - d1[mask] - .QUARTER_YEARS
        Sd <- Sd + sum(wj[mask] * dd^2)
        Nd <- Nd + sum(wj[mask])
      }
    }
    if (h1[j] == 0 && h2[j] == 1) {
      mask <- (s2 == 1)
      if (any(mask)) {
        Sg_start <- Sg_start + sum(wj[mask] * (g2[mask] - .QUARTER_YEARS)^2)
        Ng_start <- Ng_start + sum(wj[mask])
      }
    }
    if (h1[j] == 1 && h2[j] == 0) {
      mask <- (s2 == 0)
      if (any(mask)) {
        Sd_start <- Sd_start + sum(wj[mask] * (d2[mask] - .QUARTER_YEARS)^2)
        Nd_start <- Nd_start + sum(wj[mask])
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
    if (h2[j] == 0 && h3[j] == 0) {
      mask <- (s3 == 0) & (s2 == 0)
      if (any(mask)) {
        dd <- d3[mask] - d2[mask] - .QUARTER_YEARS
        Sd <- Sd + sum(wj[mask] * dd^2)
        Nd <- Nd + sum(wj[mask])
      }
    }
    if (h2[j] == 0 && h3[j] == 1) {
      mask <- (s3 == 1)
      if (any(mask)) {
        Sg_start <- Sg_start + sum(wj[mask] * (g3[mask] - .QUARTER_YEARS)^2)
        Ng_start <- Ng_start + sum(wj[mask])
      }
    }
    if (h2[j] == 1 && h3[j] == 0) {
      mask <- (s3 == 0)
      if (any(mask)) {
        Sd_start <- Sd_start + sum(wj[mask] * (d3[mask] - .QUARTER_YEARS)^2)
        Nd_start <- Nd_start + sum(wj[mask])
      }
    }
  }

  # --- EMG gradient data for joint theta M-step (TeX Eq. suff_emg_grad) ------
  # Collect (x, w) pairs for every (i, t, h) triple that contributes
  # log_emg(x; lambda_g, sigma2_g) or log_emg(x; lambda_d, sigma2_d).
  # The Newton solver in the M-step re-evaluates log_emg_grad_lambda over
  # these vectors at each candidate theta value (TeX Section 2.7).
  #
  # EMG-lambda_g observation types (TeX Section 2.7):
  #   (A+B) Wave 1: s1=1 for ANY h1 (matched or misclassified).
  #   (C) t>=2: hp=1, hc=1, s_t=1, s_{t-1}=0  (continuation, prev misclassified)
  #   (D) t>=2: hc=0, s_t=1                     (misclassified as employed)
  # EMG-lambda_d observation types (symmetric):
  #   (A+B) Wave 1: s1=0 for ANY h1.
  #   (C) t>=2: hp=0, hc=0, s_t=0, s_{t-1}=1  (continuation, prev misclassified)
  #   (D) t>=2: hc=1, s_t=0                     (misclassified as nonemployed)

  emg_g_x_list <- vector("list", H)
  emg_g_w_list <- vector("list", H)
  emg_d_x_list <- vector("list", H)
  emg_d_w_list <- vector("list", H)

  for (j in seq_len(H)) {
    wj <- wg[, j]

    # ---- EMG-lambda_g ----
    # (A+B) Wave 1: s1=1 regardless of h1
    mask <- (s1 == 1)
    ex_g <- g1[mask]
    ew_g <- wj[mask]

    # (C) t=2: hp=1, hc=1, s2=1, s1=0
    if (h1[j] == 1 && h2[j] == 1) {
      m <- (s2 == 1) & (s1 == 0)
      if (any(m)) { ex_g <- c(ex_g, g2[m]); ew_g <- c(ew_g, wj[m]) }
    }
    # (D) t=2: hc=0, s2=1
    if (h2[j] == 0) {
      m <- (s2 == 1)
      if (any(m)) { ex_g <- c(ex_g, g2[m]); ew_g <- c(ew_g, wj[m]) }
    }
    # (C) t=3: hp=1, hc=1, s3=1, s2=0
    if (h2[j] == 1 && h3[j] == 1) {
      m <- (s3 == 1) & (s2 == 0)
      if (any(m)) { ex_g <- c(ex_g, g3[m]); ew_g <- c(ew_g, wj[m]) }
    }
    # (D) t=3: hc=0, s3=1
    if (h3[j] == 0) {
      m <- (s3 == 1)
      if (any(m)) { ex_g <- c(ex_g, g3[m]); ew_g <- c(ew_g, wj[m]) }
    }

    emg_g_x_list[[j]] <- ex_g
    emg_g_w_list[[j]] <- ew_g

    # ---- EMG-lambda_d ----
    # (A+B) Wave 1: s1=0 regardless of h1
    mask <- (s1 == 0)
    ex_d <- d1[mask]
    ew_d <- wj[mask]

    # (C) t=2: hp=0, hc=0, s2=0, s1=1
    if (h1[j] == 0 && h2[j] == 0) {
      m <- (s2 == 0) & (s1 == 1)
      if (any(m)) { ex_d <- c(ex_d, d2[m]); ew_d <- c(ew_d, wj[m]) }
    }
    # (D) t=2: hc=1, s2=0
    if (h2[j] == 1) {
      m <- (s2 == 0)
      if (any(m)) { ex_d <- c(ex_d, d2[m]); ew_d <- c(ew_d, wj[m]) }
    }
    # (C) t=3: hp=0, hc=0, s3=0, s2=1
    if (h2[j] == 0 && h3[j] == 0) {
      m <- (s3 == 0) & (s2 == 1)
      if (any(m)) { ex_d <- c(ex_d, d3[m]); ew_d <- c(ew_d, wj[m]) }
    }
    # (D) t=3: hc=1, s3=0
    if (h3[j] == 1) {
      m <- (s3 == 0)
      if (any(m)) { ex_d <- c(ex_d, d3[m]); ew_d <- c(ew_d, wj[m]) }
    }

    emg_d_x_list[[j]] <- ex_d
    emg_d_w_list[[j]] <- ew_d
  }

  emg_g_x <- unlist(emg_g_x_list, use.names = FALSE)
  emg_g_w <- unlist(emg_g_w_list, use.names = FALSE)
  emg_d_x <- unlist(emg_d_x_list, use.names = FALSE)
  emg_d_w <- unlist(emg_d_w_list, use.names = FALSE)

  suff <- list(
    C1 = C1, C0 = C0,
    D1 = D1, D0 = D0,
    T11 = T11, T01 = T01,
    M = M_count,
    Sg = Sg, Ng = Ng,
    Sd = Sd, Nd = Nd,
    Sg_start = Sg_start, Ng_start = Ng_start,
    Sd_start = Sd_start, Nd_start = Nd_start,
    emg_g_x = emg_g_x, emg_g_w = emg_g_w,
    emg_d_x = emg_d_x, emg_d_w = emg_d_w
  )

  return(list(gamma = gamma_mat, loglik = ll, suff = suff))
}
