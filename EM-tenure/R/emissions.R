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

# ==============================================================================
# Discrete interval-censored timegap emissions (Issue 3 remedy)
# ==============================================================================
# These functions replace the Gaussian/EMG emissions on the nonemployment
# (timegap) side when discrete_timegap = TRUE. The underlying model is
# D ~ Exp(lambda_d) where lambda_d = -log(1 - theta0) / 3 (CTMC link).
#
# Instead of evaluating a continuous density at a midpoint, we compute the
# PROBABILITY that the true duration falls in the observed category interval.
# This is the correct likelihood for interval-censored data.
#
# TeX reference: DIAGNOSIS.md, Issue 3 Remedy; brainstorm
# .cg-docs/brainstorms/2026-03-19-em-tenure-diagnosis-expansion.md
# ==============================================================================

# --- Helper: log interval probability of Exp(lambda) -------------------------

#' Log interval probability: log P(D in [a, b) | D ~ Exp(lambda))
#'
#' Numerically stable implementation:
#'   - Finite b: log(exp(-lambda*a) - exp(-lambda*b))
#'             = -lambda*a + log1p(-exp(-lambda*(b-a)))
#'   - Infinite b: -lambda*a  (since P(D >= a) = exp(-lambda*a))
#'
#' @param a Left endpoint (>= 0).
#' @param b Right endpoint (> a; may be Inf).
#' @param lambda Exponential rate (> 0).
#' @return Scalar log probability.
#' @keywords internal
.log_interval_prob <- function(a, b, lambda) {
  if (is.infinite(b)) {
    return(-lambda * a)
  }
  # Stable: -lambda*a + log(1 - exp(-lambda*(b-a)))
  span <- lambda * (b - a)
  if (span > 700) {
    # exp(-span) underflows to 0, log1p(-0) = 0
    return(-lambda * a)
  }
  -lambda * a + log1p(-exp(-span))
}

# --- Helper: log denominator for conditional transition ----------------------

#' Log of the Exp(lambda) mass in category j: log P(D in [a_j, b_j))
#' @keywords internal
.log_cat_mass <- function(j, lambda) {
  iv <- .timegap_interval(j)
  .log_interval_prob(iv[1], iv[2], lambda)
}


# --- Discrete timegap emission functions -------------------------------------

#' Log interval probability for a nonemployment category (discrete model)
#'
#' Computes log P(D in [a_k, b_k) | D ~ Exp(lambda_d)) for category k,
#' the natural likelihood for interval-censored exponential data.
#'
#' This replaces log_emg() in six emission cases when discrete_timegap = TRUE:
#'   - Wave 1 matched nonemployment (h1=0, s1=0)
#'   - Wave 1 misclassified as nonemployed (h1=1, s1=0)
#'   - Nonemployment continuation, previous misclassified (hp=0, hc=0, s_prev=1)
#'   - Misclassified as nonemployed, t >= 2 (hc=1, s_t=0)
#'
#' @param cat Integer vector of category codes in 1:7.
#' @param lambda_d Exponential rate for nonemployment durations (> 0).
#' @return Log-probability vector. -Inf for cat outside 1:7.
#' @export
log_emission_interval_d <- function(cat, lambda_d) {
  result <- rep(-Inf, length(cat))
  valid <- !is.na(cat) & cat >= 1L & cat <= .N_TIMEGAP_CATS
  for (i in which(valid)) {
    iv <- .timegap_interval(cat[i])
    result[i] <- .log_interval_prob(iv[1], iv[2], lambda_d)
  }
  result
}

#' Log conditional transition probability for nonemployment category (discrete)
#'
#' Computes log P(c_t = cat_curr | c_{t-1} = cat_prev, lambda_d) for an
#' individual who was nonemployed at both waves t-1 and t.
#'
#' The transition probability uses the intersection interval formula:
#'   P(c_t | c_{t-1}) = P(D_{t-1} + 0.25 in [a_k, b_k) | D_{t-1} in [a_j, b_j))
#'
#' This is computed as:
#'   [exp(-lambda*L) - exp(-lambda*U)] / [exp(-lambda*a_j) - exp(-lambda*b_j)]
#' where L = max(a_j, a_k - 0.25), U = min(b_j, b_k - 0.25).
#'
#' The 7x7 transition matrix is sparse (at most 2 non-zero entries per row).
#' Categories 1-4 are DETERMINISTIC (the narrow 3-month intervals shift exactly
#' one category forward with a 0.25-year increment). Categories 5-6 are
#' PROBABILISTIC (two possible destinations).
#'
#' This replaces log_emission_increment_d() for Case 3:
#'   Nonemployment continuation, previous correctly observed (hp=0, hc=0,
#'   s_t=0, s_{t-1}=0).
#'
#' @param cat_curr Integer vector: current category codes (1:7).
#' @param cat_prev Integer vector: previous category codes (1:7), same length.
#' @param lambda_d Exponential rate for nonemployment durations (> 0).
#' @return Log-probability vector. -Inf for invalid or unreachable transitions.
#' @export
log_emission_transition_d <- function(cat_curr, cat_prev, lambda_d) {
  n <- length(cat_curr)
  stopifnot(length(cat_prev) == n)
  result <- rep(-Inf, n)

  for (i in seq_len(n)) {
    j <- cat_prev[i]
    k <- cat_curr[i]
    if (is.na(j) || is.na(k) || j < 1L || j > .N_TIMEGAP_CATS ||
        k < 1L || k > .N_TIMEGAP_CATS) {
      next
    }

    iv_j <- .timegap_interval(j)
    iv_k <- .timegap_interval(k)
    a_j <- iv_j[1]; b_j <- iv_j[2]
    a_k <- iv_k[1]; b_k <- iv_k[2]

    # Intersection interval: the part of [a_j, b_j) that, after adding 0.25,
    # maps into [a_k, b_k). Equivalently, D_{t-1} in [a_k-0.25, b_k-0.25)
    # intersected with [a_j, b_j).
    L <- max(a_j, a_k - .QUARTER_YEARS)
    U <- min(b_j, if (is.infinite(b_k)) Inf else b_k - .QUARTER_YEARS)

    if (L >= U) next  # unreachable transition

    log_numerator   <- .log_interval_prob(L, U, lambda_d)
    log_denominator <- .log_cat_mass(j, lambda_d)

    if (is.infinite(log_denominator) && log_denominator < 0) next

    result[i] <- log_numerator - log_denominator
  }

  result
}

#' Log emission for within-panel nonemployment start (discrete model)
#'
#' When t >= 2, h_{t-1} = 1, h_t = 0, s_t = 0: a new nonemployment spell
#' began this quarter. Its duration is at most 0.25 years (one quarter), so
#' it MUST fall in category 1 ([0, 0.25) years = [0, 3) months).
#' Any other category is structurally impossible.
#'
#' This replaces log_emission_start_d() when discrete_timegap = TRUE.
#'
#' @param cat Integer vector of category codes.
#' @return Log-probability vector: 0 if cat == 1, -Inf otherwise.
#' @export
log_emission_start_d_cat <- function(cat) {
  ifelse(!is.na(cat) & cat == 1L, 0, -Inf)
}

#' Gradient of log interval probability w.r.t. lambda_d (discrete model)
#'
#' Computes d/d(lambda_d) log P(D in [a_k, b_k) | D ~ Exp(lambda_d)).
#'
#' For finite b_k:
#'   d/d(lambda) log(exp(-lambda*a) - exp(-lambda*b))
#'   = (-a * exp(-lambda*a) + b * exp(-lambda*b)) / (exp(-lambda*a) - exp(-lambda*b))
#'
#' For b_k = Inf:
#'   d/d(lambda) log(exp(-lambda*a)) = -a
#'
#' Used in the M-step FOC for theta_0 when discrete_timegap = TRUE.
#'
#' @param cat Integer vector of category codes in 1:7.
#' @param lambda_d Exponential rate (> 0).
#' @return Gradient vector. NA for cat outside 1:7.
#' @export
interval_grad_lambda_d <- function(cat, lambda_d) {
  result <- rep(NA_real_, length(cat))
  valid <- !is.na(cat) & cat >= 1L & cat <= .N_TIMEGAP_CATS
  for (i in which(valid)) {
    iv <- .timegap_interval(cat[i])
    a  <- iv[1]; b <- iv[2]
    if (is.infinite(b)) {
      result[i] <- -a
    } else {
      ea <- exp(-lambda_d * a)
      eb <- exp(-lambda_d * b)
      denom <- ea - eb
      if (abs(denom) < 1e-15) {
        # Near-zero denominator: use limit (L'Hopital) -> -(a+b)/2
        result[i] <- -(a + b) / 2
      } else {
        result[i] <- (-a * ea + b * eb) / denom
      }
    }
  }
  result
}

#' Gradient of log conditional transition probability w.r.t. lambda_d
#'
#' Computes d/d(lambda_d) log P(c_t | c_{t-1}, lambda_d), i.e., the gradient
#' of the log of the conditional transition probability used in Case 3.
#'
#' This equals: interval_grad on [L, U) - interval_grad on [a_j, b_j),
#' where [L, U) is the intersection interval from log_emission_transition_d().
#'
#' Used in the M-step FOC for theta_0 when discrete_timegap = TRUE.
#'
#' @param cat_curr Integer vector: current category codes (1:7).
#' @param cat_prev Integer vector: previous category codes (1:7), same length.
#' @param lambda_d Exponential rate (> 0).
#' @return Gradient vector. NA for invalid or unreachable transitions.
#' @export
transition_grad_lambda_d <- function(cat_curr, cat_prev, lambda_d) {
  n <- length(cat_curr)
  stopifnot(length(cat_prev) == n)
  result <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    j <- cat_prev[i]
    k <- cat_curr[i]
    if (is.na(j) || is.na(k) || j < 1L || j > .N_TIMEGAP_CATS ||
        k < 1L || k > .N_TIMEGAP_CATS) {
      next
    }

    iv_j <- .timegap_interval(j)
    iv_k <- .timegap_interval(k)
    a_j <- iv_j[1]; b_j <- iv_j[2]
    a_k <- iv_k[1]; b_k <- iv_k[2]

    L <- max(a_j, a_k - .QUARTER_YEARS)
    U <- min(b_j, if (is.infinite(b_k)) Inf else b_k - .QUARTER_YEARS)

    if (L >= U) next  # unreachable transition

    # grad = d/dlambda log P([L,U)) - d/dlambda log P([a_j, b_j))
    # Each term is computed using the same formula as interval_grad_lambda_d
    grad_num <- if (is.infinite(U)) {
      -L
    } else {
      eL <- exp(-lambda_d * L); eU <- exp(-lambda_d * U)
      dLU <- eL - eU
      if (abs(dLU) < 1e-15) -(L + U) / 2 else (-L * eL + U * eU) / dLU
    }

    ea <- exp(-lambda_d * a_j)
    eb <- if (is.infinite(b_j)) 0 else exp(-lambda_d * b_j)
    denom <- ea - eb
    grad_den <- if (abs(denom) < 1e-15) {
      -(a_j + (if (is.infinite(b_j)) a_j else b_j)) / 2
    } else {
      (-a_j * ea + (if (is.infinite(b_j)) 0 else b_j * eb)) / denom
    }

    result[i] <- grad_num - grad_den
  }

  result
}
