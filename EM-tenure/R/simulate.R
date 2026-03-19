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
#' @param sigma2_d Nonemployment measurement variance.
#' @param seed Optional random seed.
#' @return Data frame with columns: y1-y3, tenure1-tenure3, timegap1-timegap3,
#'   weight, h1-h3 (true latent states).
#' @export
simulate_panel <- function(n,
                           alpha = 0.6,
                           theta1 = 0.9,
                           theta0 = 0.1,
                           pi = 0.05,
                           sigma2_g = 0.01,
                           sigma2_d = 0.01,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  lambda_g <- ctmc_lambda_from_theta(theta1)
  lambda_d <- ctmc_lambda_from_theta(theta0)
  sigma_g  <- sqrt(sigma2_g)
  sigma_d  <- sqrt(sigma2_d)

  # Pre-allocate

  h <- matrix(NA_integer_, nrow = n, ncol = 3)
  s <- matrix(NA_integer_, nrow = n, ncol = 3)
  g <- matrix(0, nrow = n, ncol = 3)
  d <- matrix(0, nrow = n, ncol = 3)

  # --- Wave 1: all vectorised ---
  h[, 1] <- as.integer(runif(n) < alpha)
  s[, 1] <- ifelse(runif(n) < pi, 1L - h[, 1], h[, 1])

  # Wave 1 durations: all EMG
  emg_g1 <- abs(rexp(n, rate = lambda_g) + rnorm(n, 0, sigma_g))
  emg_d1 <- abs(rexp(n, rate = lambda_d) + rnorm(n, 0, sigma_d))

  # Assign based on observed state
  emp1 <- (s[, 1] == 1)
  g[emp1, 1]  <- emg_g1[emp1]
  d[!emp1, 1] <- emg_d1[!emp1]

  # --- Waves 2-3: vectorised ---
  for (tt in 2:3) {
    # State transitions (vectorised)
    was_emp <- (h[, tt - 1] == 1)
    p_emp <- ifelse(was_emp, theta1, theta0)
    h[, tt] <- as.integer(runif(n) < p_emp)

    # Misclassification (vectorised)
    s[, tt] <- ifelse(runif(n) < pi, 1L - h[, tt], h[, tt])

    # Duration generation (vectorised by case)
    is_emp    <- (h[, tt] == 1)
    was_emp_t <- (h[, tt - 1] == 1)
    obs_emp   <- (s[, tt] == 1)
    prev_obs_emp <- (s[, tt - 1] == 1)

    # --- Observed as employed ---
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
    # Case 3: Truly employed, start (h_{t-1}=0)
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

    # --- Observed as nonemployed ---
    obs_non <- !obs_emp

    # Case 5: Truly nonemployed, continuation, previous observed
    mask <- obs_non & !is_emp & !was_emp_t & !prev_obs_emp
    if (any(mask)) {
      nm <- sum(mask)
      d[mask, tt] <- abs(d[mask, tt - 1] + .QUARTER_YEARS +
                         rnorm(nm, 0, sqrt(2) * sigma_d))
    }
    # Case 6: Truly nonemployed, continuation, previous misclassified
    mask <- obs_non & !is_emp & !was_emp_t & prev_obs_emp
    if (any(mask)) {
      nm <- sum(mask)
      d[mask, tt] <- abs(rexp(nm, lambda_d) + rnorm(nm, 0, sigma_d))
    }
    # Case 7: Truly nonemployed, start (h_{t-1}=1)
    mask <- obs_non & !is_emp & was_emp_t
    if (any(mask)) {
      nm <- sum(mask)
      d[mask, tt] <- abs(.QUARTER_YEARS + rnorm(nm, 0, sigma_d))
    }
    # Case 8: Misclassified as nonemployed (truly employed)
    mask <- obs_non & is_emp
    if (any(mask)) {
      nm <- sum(mask)
      d[mask, tt] <- abs(rexp(nm, lambda_d) + rnorm(nm, 0, sigma_d))
    }
  }

  return(data.frame(
    y1 = s[, 1], y2 = s[, 2], y3 = s[, 3],
    tenure1 = g[, 1], tenure2 = g[, 2], tenure3 = g[, 3],
    timegap1 = d[, 1], timegap2 = d[, 2], timegap3 = d[, 3],
    weight = rep(1, n),
    h1 = h[, 1], h2 = h[, 2], h3 = h[, 3]
  ))
}
