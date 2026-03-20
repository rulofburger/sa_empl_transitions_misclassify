# ==============================================================================
# EM-tenure: Synthetic data generator (vectorised)
# ==============================================================================
# Simulates panel data from the linked-specification model for testing and
# validation. The DGP follows the TeX exactly:
#   1. Draw true history from Markov chain (alpha, theta1, theta0)
#   2. Misclassify with probability pi
#   3. Generate durations:
#      - Wave 1: EMG(lambda, sigma^2) for all
#      - Continuation (h_{t-1}=h_t): increment ~ N(0.25, 2*sigma^2)
#      - Start (h_{t-1}!=h_t): g_t ~ N(0.25, sigma^2)
#      - Misclassified: EMG with inappropriate lambda/sigma
#
# VECTORISATION: All draws are vectorised across individuals.
# ==============================================================================

#' Simulate panel data from the EM-tenure model
#'
#' @param n Number of individuals to simulate.
#' @param alpha Initial employment probability.
#' @param theta1 P(employed | employed last period).
#' @param theta0 P(employed | nonemployed last period).
#' @param pi Misclassification probability.
#' @param sigma2_g Employment measurement variance.
#' @param sigma2_d Nonemployment measurement variance. Ignored (and set to NA
#'   in output) when \code{discrete_timegap = TRUE}.
#' @param discrete_timegap Logical (default \code{TRUE}). When \code{TRUE},
#'   nonemployment durations are represented as integer category codes (1–7)
#'   matching the QLFS survey instrument, and the timegap columns are
#'   populated via the discrete interval-censored model. When \code{FALSE},
#'   the legacy continuous EMG/Gaussian model is used.
#' @param seed Optional random seed.
#' @return Data frame with columns:
#'   \itemize{
#'     \item y1–y3: observed employment (0/1)
#'     \item tenure1–tenure3: observed tenure in years
#'     \item timegap1–timegap3: observed nonemployment duration (midpoint in
#'       years when \code{discrete_timegap = FALSE}; category midpoint in years
#'       when \code{discrete_timegap = TRUE})
#'     \item timegap_cat1–timegap_cat3: integer category codes 1–7 (only when
#'       \code{discrete_timegap = TRUE})
#'     \item weight: survey weight (all 1)
#'     \item h1–h3: true latent employment states
#'   }
#' @export
simulate_panel <- function(n,
                           alpha = 0.6,
                           theta1 = 0.9,
                           theta0 = 0.1,
                           pi = 0.05,
                           sigma2_g = 0.01,
                           sigma2_d = 0.01,
                           discrete_timegap = TRUE,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  lambda_g <- ctmc_lambda_from_theta(theta1)
  lambda_d <- ctmc_lambda_from_theta(theta0)
  sigma_g  <- sqrt(sigma2_g)
  sigma_d  <- if (discrete_timegap) NA_real_ else sqrt(sigma2_d)

  # Pre-allocate
  h <- matrix(NA_integer_, nrow = n, ncol = 3)
  s <- matrix(NA_integer_, nrow = n, ncol = 3)
  g <- matrix(0, nrow = n, ncol = 3)
  d <- matrix(0, nrow = n, ncol = 3)
  d_cat <- matrix(NA_integer_, nrow = n, ncol = 3)  # category codes

  # --- Wave 1: all vectorised ---
  h[, 1] <- as.integer(runif(n) < alpha)
  s[, 1] <- ifelse(runif(n) < pi, 1L - h[, 1], h[, 1])

  if (discrete_timegap) {
    # Wave 1: draw true durations from Exp(lambda_d), map to category codes
    true_d1 <- rexp(n, rate = lambda_d)
    d_cat[, 1] <- .continuous_to_cat(true_d1)
    d[, 1] <- .TIMEGAP_MIDPOINTS_MONTHS[d_cat[, 1]] / 12  # midpoint in years

    # Wave 1 tenure: EMG (unchanged)
    emg_g1 <- abs(rexp(n, rate = lambda_g) + rnorm(n, 0, sigma_g))
    emp1 <- (s[, 1] == 1)
    g[emp1,  1] <- emg_g1[emp1]
    # Wrong-state tenure for nonemployed at wave 1: use DURATION_FLOOR
    # (will be imputed by nearest-non-zero in the real data pipeline;
    # here we assign a floor for simulator internal consistency)
    g[!emp1, 1] <- .DURATION_FLOOR

    # Wrong-state timegap category for employed at wave 1: category 1 (floor)
    d_cat[emp1, 1] <- 1L
    d[emp1, 1] <- .TIMEGAP_MIDPOINTS_MONTHS[1L] / 12

  } else {
    # Legacy continuous path (unchanged from original)
    emg_g1 <- abs(rexp(n, rate = lambda_g) + rnorm(n, 0, sigma_g))
    emg_d1 <- abs(rexp(n, rate = lambda_d) + rnorm(n, 0, sigma_d))
    emp1 <- (s[, 1] == 1)
    g[emp1,  1] <- emg_g1[emp1]
    d[!emp1, 1] <- emg_d1[!emp1]
  }

  # --- Waves 2-3: vectorised ---
  for (tt in 2:3) {
    was_emp <- (h[, tt - 1] == 1)
    p_emp <- ifelse(was_emp, theta1, theta0)
    h[, tt] <- as.integer(runif(n) < p_emp)

    s[, tt] <- ifelse(runif(n) < pi, 1L - h[, tt], h[, tt])

    is_emp      <- (h[, tt] == 1)
    was_emp_t   <- (h[, tt - 1] == 1)
    obs_emp     <- (s[, tt] == 1)
    prev_obs_emp <- (s[, tt - 1] == 1)
    obs_non     <- !obs_emp

    # --- Observed as employed: tenure generation (same for both modes) ---
    # Case 1: Truly employed, continuation, previous observed
    mask <- obs_emp & is_emp & was_emp_t & prev_obs_emp
    if (any(mask)) {
      nm <- sum(mask)
      g[mask, tt] <- abs(g[mask, tt - 1] + .QUARTER_YEARS +
                         rnorm(nm, 0, sqrt(2) * sigma_g))
    }
    # Case 2: Truly employed, continuation, previous misclassified
    mask <- obs_emp & is_emp & was_emp_t & !prev_obs_emp
    if (any(mask)) {
      nm <- sum(mask)
      g[mask, tt] <- abs(rexp(nm, lambda_g) + rnorm(nm, 0, sigma_g))
    }
    # Case 3: Truly employed, start
    mask <- obs_emp & is_emp & !was_emp_t
    if (any(mask)) {
      nm <- sum(mask)
      g[mask, tt] <- abs(.QUARTER_YEARS + rnorm(nm, 0, sigma_g))
    }
    # Case 4: Misclassified as employed (truly nonemployed)
    mask <- obs_emp & !is_emp
    if (any(mask)) {
      nm <- sum(mask)
      g[mask, tt] <- abs(rexp(nm, lambda_g) + rnorm(nm, 0, sigma_g))
    }

    if (discrete_timegap) {
      # --- Observed as nonemployed: discrete category generation ---

      # Case 5: Truly nonemp, continuation, prev observed (s_{t-1}=0 → cat known)
      mask <- obs_non & !is_emp & !was_emp_t & !prev_obs_emp
      if (any(mask)) {
        # True duration = D_{t-1} + 0.25; D_{t-1} in [a_{c_{t-1}}, b_{c_{t-1}})
        # We use midpoint of prev category as point estimate for D_{t-1}
        prev_d_approx <- .TIMEGAP_MIDPOINTS_MONTHS[d_cat[mask, tt - 1]] / 12
        true_d_new <- prev_d_approx + .QUARTER_YEARS
        d_cat[mask, tt] <- .continuous_to_cat(true_d_new)
        d[mask, tt] <- .TIMEGAP_MIDPOINTS_MONTHS[d_cat[mask, tt]] / 12
      }
      # Case 6: Truly nonemp, continuation, prev misclassified (prev cat unknown)
      mask <- obs_non & !is_emp & !was_emp_t & prev_obs_emp
      if (any(mask)) {
        nm <- sum(mask)
        true_d <- rexp(nm, lambda_d)
        d_cat[mask, tt] <- .continuous_to_cat(true_d)
        d[mask, tt] <- .TIMEGAP_MIDPOINTS_MONTHS[d_cat[mask, tt]] / 12
      }
      # Case 7: Truly nonemp, start (new spell: duration < 0.25 years → cat 1)
      mask <- obs_non & !is_emp & was_emp_t
      if (any(mask)) {
        d_cat[mask, tt] <- 1L
        d[mask, tt] <- .TIMEGAP_MIDPOINTS_MONTHS[1L] / 12
      }
      # Case 8: Misclassified as nonemp (truly employed)
      mask <- obs_non & is_emp
      if (any(mask)) {
        nm <- sum(mask)
        true_d <- rexp(nm, lambda_d)
        d_cat[mask, tt] <- .continuous_to_cat(true_d)
        d[mask, tt] <- .TIMEGAP_MIDPOINTS_MONTHS[d_cat[mask, tt]] / 12
      }

      # Wrong-state tenure for nonemployed: DURATION_FLOOR
      d_floor <- (s[, tt] == 0)
      g[d_floor, tt] <- ifelse(g[d_floor, tt] == 0, .DURATION_FLOOR, g[d_floor, tt])

      # Wrong-state timegap_cat for employed: category 1 (floor)
      emp_mask <- (s[, tt] == 1)
      d_cat[emp_mask, tt] <- 1L
      d[emp_mask, tt] <- .TIMEGAP_MIDPOINTS_MONTHS[1L] / 12

    } else {
      # --- Legacy continuous path ---
      # Case 5: Truly nonemp, continuation, prev observed
      mask <- obs_non & !is_emp & !was_emp_t & !prev_obs_emp
      if (any(mask)) {
        nm <- sum(mask)
        d[mask, tt] <- abs(d[mask, tt - 1] + .QUARTER_YEARS +
                           rnorm(nm, 0, sqrt(2) * sigma_d))
      }
      # Case 6: Truly nonemp, continuation, prev misclassified
      mask <- obs_non & !is_emp & !was_emp_t & prev_obs_emp
      if (any(mask)) {
        nm <- sum(mask)
        d[mask, tt] <- abs(rexp(nm, lambda_d) + rnorm(nm, 0, sigma_d))
      }
      # Case 7: Truly nonemp, start
      mask <- obs_non & !is_emp & was_emp_t
      if (any(mask)) {
        nm <- sum(mask)
        d[mask, tt] <- abs(.QUARTER_YEARS + rnorm(nm, 0, sigma_d))
      }
      # Case 8: Misclassified as nonemp
      mask <- obs_non & is_emp
      if (any(mask)) {
        nm <- sum(mask)
        d[mask, tt] <- abs(rexp(nm, lambda_d) + rnorm(nm, 0, sigma_d))
      }
    }
  }

  out <- data.frame(
    y1 = s[, 1], y2 = s[, 2], y3 = s[, 3],
    tenure1 = g[, 1], tenure2 = g[, 2], tenure3 = g[, 3],
    timegap1 = d[, 1], timegap2 = d[, 2], timegap3 = d[, 3],
    weight = rep(1, n),
    h1 = h[, 1], h2 = h[, 2], h3 = h[, 3]
  )

  if (discrete_timegap) {
    out$timegap_cat1 <- d_cat[, 1]
    out$timegap_cat2 <- d_cat[, 2]
    out$timegap_cat3 <- d_cat[, 3]
  }

  return(out)
}
