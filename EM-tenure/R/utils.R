# ==============================================================================
# EM-tenure: Internal utility functions
# ==============================================================================
# These are low-level helpers used across all modules. They are not exported
# for direct use but are essential for numerical stability and convenience.
# ==============================================================================

#' Quarterly increment in years. All modules derive the inter-wave
#' time step from this single constant (TeX: Delta = 3 months = 0.25 years).
#' @keywords internal
.QUARTER_YEARS <- 0.25

#' Null coalescing operator
#' Returns `a` if non-NULL, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

#' Clip probabilities to (eps, 1-eps)
#'
#' Prevents log(0) and log(1) in likelihood computations.
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive bound (default 1e-10).
#' @return Clipped probabilities in (eps, 1-eps).
#' @keywords internal
.bound01 <- function(p, eps = 1e-10) {
  pmin(1 - eps, pmax(eps, p))
}

#' Enforce strict positivity with a floor
#'
#' Used for variances, rates, and denominators to avoid division by zero.
#'
#' @param x Numeric vector.
#' @param eps Small positive floor (default 1e-12).
#' @return Values >= eps.
#' @keywords internal
.pos <- function(x, eps = 1e-12) {
  pmax(eps, x)
}

#' Safe logit transform
#'
#' Returns NA at boundaries 0 and 1 (where logit is undefined).
#'
#' @param p Probability in (0, 1).
#' @return log(p/(1-p)), or NA if p is on boundary.
#' @keywords internal
.safe_logit <- function(p) {
  ifelse(p <= 0 | p >= 1, NA_real_, log(p / (1 - p)))
}

#' Numerically stable log-sum-exp
#'
#' Computes log(sum(exp(x))) without overflow/underflow.
#' Shifts by max(x) before exponentiating.
#'
#' @param x Numeric vector (may contain -Inf).
#' @return Scalar: log(sum(exp(x))).
#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)  # all -Inf
  m + log(sum(exp(x - m)))
}

#' Complementary error function
#'
#' Computes erfc(x) = 1 - erf(x) = 2 * Phi(-sqrt(2) * x),
#' where Phi is the standard Normal CDF. Numerically stable for all finite x.
#'
#' @param x Numeric vector.
#' @return erfc(x).
#' @keywords internal
erfc <- function(x) {
  2 * pnorm(-sqrt(2) * x)
}
