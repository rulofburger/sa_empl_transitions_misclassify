# ==============================================================================
# EM-baseline: Implied probability transforms
# Created: 2026-05-06
#
# Functions that map raw EM parameter lists to economically interpretable
# implied quantities: job entry rate, job exit rate, employment rate, and
# misclassification rate(s).  These are used both in the point-estimate
# pipeline and in each bootstrap replicate so that bootstrap SEs cover the
# implied quantities directly.
# ==============================================================================

#' Compute implied probabilities from a baseline EM parameter list
#'
#' Converts raw EM parameter estimates into economically interpretable
#' quantities.  For the symmetric and no-error variants the transformation is
#' trivial; the asymmetric variant returns separate false-positive and
#' false-negative rates.  The employment rate is always derived from the
#' stationarity formula \eqn{\alpha = \theta_0 / (\theta_0 + 1 - \theta_1)},
#' regardless of whether the model was estimated with a free \eqn{\alpha},
#' so that it is comparable across stationary and free-\eqn{\alpha} models.
#'
#' @param params  Named list of EM parameters.  Must contain \code{theta0} and
#'   \code{theta1}.  May contain \code{pi} (symmetric), \code{pi0} and
#'   \code{pi1} (asymmetric), or neither (no-error variant).
#' @param model_type  Character scalar: \code{"symmetric"}, \code{"asymmetric"},
#'   or \code{"none"}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{entry_rate}{Job-finding probability per period (\eqn{\theta_0}).}
#'     \item{exit_rate}{Job-separation probability per period
#'       (\eqn{1 - \theta_1}).}
#'     \item{employment_rate}{Steady-state employment rate derived from the
#'       stationarity formula.}
#'     \item{pi}{Symmetric misclassification probability.  \code{NA} for
#'       asymmetric and no-error models.}
#'     \item{pi0}{False-positive rate (\eqn{\pi_0}).  \code{NA} except for
#'       asymmetric models.}
#'     \item{pi1}{False-negative rate (\eqn{\pi_1}).  \code{NA} except for
#'       asymmetric models.}
#'   }
#'
#' @examples
#' p <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
#' implied_baseline(p, "symmetric")
#'
#' p_asym <- list(theta0 = 0.10, theta1 = 0.90, pi0 = 0.03, pi1 = 0.07)
#' implied_baseline(p_asym, "asymmetric")
implied_baseline <- function(params, model_type = "symmetric") {
  stopifnot(
    is.list(params),
    is.numeric(params$theta0), length(params$theta0) == 1L, !is.na(params$theta0),
    params$theta0 >= 0, params$theta0 <= 1,
    is.numeric(params$theta1), length(params$theta1) == 1L, !is.na(params$theta1),
    params$theta1 >= 0, params$theta1 <= 1,
    model_type %in% c("symmetric", "asymmetric", "none")
  )

  theta0 <- params$theta0
  theta1 <- params$theta1

  # Steady-state employment rate from stationarity formula, regardless of
  # whether alpha was a free parameter or constrained.
  denom <- theta0 + (1 - theta1)
  employment_rate <- if (denom > 0) theta0 / denom else NA_real_

  pi  <- NA_real_
  pi0 <- NA_real_
  pi1 <- NA_real_

  if (model_type == "symmetric") {
    pi <- params$pi %||% NA_real_
    if (!is.na(pi) && (pi < 0 || pi > 1))
      stop(sprintf("implied_baseline: pi out of [0,1]: %g", pi))
    if (is.na(pi))
      warning("implied_baseline: pi is NA for symmetric model — params may be incomplete.")
  } else if (model_type == "asymmetric") {
    pi0 <- params$pi0 %||% NA_real_
    pi1 <- params$pi1 %||% NA_real_
    if (!is.na(pi0) && (pi0 < 0 || pi0 > 1))
      stop(sprintf("implied_baseline: pi0 out of [0,1]: %g", pi0))
    if (!is.na(pi1) && (pi1 < 0 || pi1 > 1))
      stop(sprintf("implied_baseline: pi1 out of [0,1]: %g", pi1))
  }
  # model_type == "none": all NA, which is correct.

  list(
    entry_rate       = theta0,
    exit_rate        = 1 - theta1,
    employment_rate  = employment_rate,
    pi               = pi,
    pi0              = pi0,
    pi1              = pi1
  )
}
