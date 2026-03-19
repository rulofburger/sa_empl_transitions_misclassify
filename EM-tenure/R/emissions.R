# ==============================================================================
# EM-tenure: Emission log-densities
# ==============================================================================
# Each function computes the log-density contribution of observed durations
# given a latent history. The emission type depends on:
#   (1) Whether the observation is correctly classified or misclassified
#   (2) Whether it is wave 1 (left-censored), a continuation, or a start
#   (3) Whether the previous wave was correctly observed (for increments)
#
# TeX reference: Section 2.7 (E-step emission densities)
# ==============================================================================

# --- Misclassification probabilities -----------------------------------------

#' Log-probability of observed state sequence given latent history
#'
#' Computes sum_t log P(s_t | h_t, pi) for the symmetric misclassification model:
#'   P(s_t = h_t) = 1 - pi
#'   P(s_t != h_t) = pi
#'
#' When pi = 0: returns 0 for exact matches, -Inf otherwise.
#'
#' TeX ref: Eq (5)
#'
#' @param s_vec Length-3 observed employment sequence (0/1).
#' @param hmat 8 x 3 matrix of latent histories.
#' @param pi Misclassification probability in [0, 1].
#' @return Length-8 vector of log-probabilities.
#' @export
log_misclass_prob <- function(s_vec, hmat, pi) {
  if (pi == 0) {
    match_all <- (hmat[, 1] == s_vec[1]) &
                 (hmat[, 2] == s_vec[2]) &
                 (hmat[, 3] == s_vec[3])
    return(ifelse(match_all, 0, -Inf))
  }

  pi <- .bound01(pi)
  lp <- 0
  for (t in 1:3) {
    lp <- lp + ifelse(hmat[, t] == s_vec[t], log1p(-pi), log(pi))
  }
  lp
}

# --- ExGaussian (EMG) density ------------------------------------------------

#' Log-density of the Exponentially Modified Gaussian (EMG)
#'
#' The EMG models the convolution of Exp(lambda) and N(0, sigma^2):
#'   X = E + epsilon, where E ~ Exp(lambda), epsilon ~ N(0, sigma^2)
#'
#' Used for:
#'   - Wave 1 matched durations (left-censored spells)
#'   - Continuation with previous misclassified (lagged duration unobserved)
#'   - Misclassified durations (inappropriate clock)
#'
#' Implementation uses the numerically stable form:
#'   log f_EMG(x) = log(lambda) + lambda^2*sigma^2/2 - lambda*x
#'                  + log(erfc((lambda*sigma^2 - x) / (sqrt(2)*sigma)))
#'                  - log(2)
#'
#' Returns -Inf for x <= 0 (strict support).
#'
#' TeX ref: Used in multiple emission cases in Section 2.7
#'
#' @param x Observed duration (scalar or vector, in years).
#' @param lambda Exponential rate (per month, from CTMC link).
#' @param sigma2 Gaussian variance (in years^2).
#' @return Log-density value(s). -Inf for x <= 0.
#' @export
log_emg <- function(x, lambda, sigma2) {
  sigma <- sqrt(.pos(sigma2))
  # Compute the argument of erfc
  z <- (lambda * sigma2 - x) / (sqrt(2) * sigma)

  log_density <- log(lambda) +
    (lambda^2 * sigma2 / 2) -
    lambda * x +
    log(erfc(z)) -
    log(2)

  ifelse(x > 0, log_density, -Inf)
}

# --- EMG gradient w.r.t. lambda ----------------------------------------------

#' Gradient of the EMG log-density with respect to the exponential rate lambda
#'
#' Computes \eqn{\partial \log f_{\mathrm{EMG}}(x; \lambda, \sigma^2) / \partial \lambda}
#' for use in the joint M-step first-order condition for \eqn{\theta}.
#'
#' The EMG log-density is (TeX Eq log_emg):
#'   log lambda + lambda^2*sigma^2/2 - lambda*x + log Phi(u)
#' where u = (x - lambda*sigma^2) / sigma.
#'
#' Differentiating with respect to lambda (TeX Eq demg_dlambda):
#'   d/dlambda log f_EMG = 1/lambda + lambda*sigma^2 - x - sigma * phi(u)/Phi(u)
#'
#' The inverse Mills ratio phi(u)/Phi(u) is computed in log-space as
#'   exp(dnorm(u, log=TRUE) - pnorm(u, log.p=TRUE))
#' to avoid underflow when u is very negative (where Phi(u) -> 0).
#'
#' Returns NA for x <= 0 (outside the strict support of the EMG).
#'
#' TeX ref: Eqs (log_emg), (demg_dlambda), (suff_emg_grad)
#'
#' @param x Observed duration (scalar or vector, in years). Must be > 0.
#' @param lambda Exponential rate (per month, from CTMC link). Must be > 0.
#' @param sigma2 Gaussian variance (in years^2). Must be > 0.
#' @return Gradient value(s); NA for x <= 0.
#' @export
log_emg_grad_lambda <- function(x, lambda, sigma2) {
  sigma <- sqrt(.pos(sigma2))
  u     <- (x - lambda * sigma2) / sigma

  # Inverse Mills ratio: phi(u) / Phi(u) computed in log-space for stability
  log_mills <- dnorm(u, log = TRUE) - pnorm(u, log.p = TRUE)
  mills      <- exp(log_mills)

  grad <- 1 / lambda + lambda * sigma2 - x - sigma * mills

  ifelse(x > 0, grad, NA_real_)
}

# --- Normal increment density (continuations) --------------------------------

#' Log-density for employment continuation increment
#'
#' When h_{t-1} = h_t = 1, s_t = 1, s_{t-1} = 1 (both waves correctly observed),
#' the tenure increment has the distribution:
#'   Delta^g_t = (g_t - g_{t-1}) - 0.25 ~ N(0, 2*sigma_g^2)
#'
#' Variance is 2*sigma_g^2 because Var(eps_t - eps_{t-1}) = 2*sigma_g^2.
#'
#' TeX ref: Eq (6)
#'
#' @param delta_g Observed increment minus 0.25 (scalar or vector).
#' @param sigma2_g Employment measurement variance.
#' @return Log-density value(s).
#' @export
log_emission_increment_g <- function(delta_g, sigma2_g) {
  v <- 2 * .pos(sigma2_g)
  -0.5 * log(2 * base::pi * v) - 0.5 * delta_g^2 / v
}

#' Log-density for nonemployment continuation increment
#'
#' When h_{t-1} = h_t = 0, s_t = 0, s_{t-1} = 0:
#'   Delta^d_t = (d_t - d_{t-1}) - 0.25 ~ N(0, 2*sigma_d^2)
#'
#' TeX ref: Eq (7)
#'
#' @param delta_d Observed increment minus 0.25 (scalar or vector).
#' @param sigma2_d Nonemployment measurement variance.
#' @return Log-density value(s).
#' @export
log_emission_increment_d <- function(delta_d, sigma2_d) {
  v <- 2 * .pos(sigma2_d)
  -0.5 * log(2 * base::pi * v) - 0.5 * delta_d^2 / v
}

# --- Normal level density (within-panel starts) ------------------------------

#' Log-density for within-panel employment start
#'
#' When t >= 2, h_{t-1} = 0, h_t = 1, s_t = 1:
#'   g_t - 0.25 ~ N(0, sigma_g^2)
#'
#' The person started employment this quarter; true tenure is exactly 0.25 years.
#' Observed tenure deviates by measurement error only.
#'
#' TeX ref: Section 2.7, "Within-panel employment start"
#'
#' @param g Observed tenure at the start wave.
#' @param sigma2_g Employment measurement variance.
#' @return Log-density value(s).
#' @export
log_emission_start_g <- function(g, sigma2_g) {
  v <- .pos(sigma2_g)
  -0.5 * log(2 * base::pi * v) - 0.5 * (g - .QUARTER_YEARS)^2 / v
}

#' Log-density for within-panel nonemployment start
#'
#' When t >= 2, h_{t-1} = 1, h_t = 0, s_t = 0:
#'   d_t - 0.25 ~ N(0, sigma_d^2)
#'
#' TeX ref: Section 2.7, "Within-panel nonemployment start"
#'
#' @param d Observed nonemployment duration at the start wave.
#' @param sigma2_d Nonemployment measurement variance.
#' @return Log-density value(s).
#' @export
log_emission_start_d <- function(d, sigma2_d) {
  v <- .pos(sigma2_d)
  -0.5 * log(2 * base::pi * v) - 0.5 * (d - .QUARTER_YEARS)^2 / v
}
