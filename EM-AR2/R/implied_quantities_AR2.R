# ==============================================================================
# EM-AR2: Implied probability transforms
# Created: 2026-05-07
#
# Maps raw EM-AR2 parameter lists to economically interpretable implied
# quantities:
#   - The four AR(2) transition probabilities p_{jk} (eq. 4 of the paper)
#   - Baseline entry and exit rates (theta0, 1-theta1)
#   - Steady-state employment rate (from stationary_ar2())
#   - Misclassification rate(s)
#
# The p_{jk} are exact algebraic functions of theta — no simulation required.
# Employment rate uses the AR(2) stationary distribution formula (eqs. 9–12).
# ==============================================================================

# Null-coalescing operator — defined locally so this file can be sourced
# independently without relying on utils.R being loaded first.
if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

#' Compute implied quantities from an AR(2) EM parameter list
#'
#' Converts raw EM-AR2 parameter estimates into economically interpretable
#' quantities. The four transition probabilities \eqn{p_{jk}} are defined in
#' eq. (4) of the paper:
#'
#' \deqn{p_{00} = \theta_0}
#' \deqn{p_{10} = \theta_0 + \theta_{01}}
#' \deqn{p_{01} = 1 - \theta_1 - \theta_{10}}
#' \deqn{p_{11} = 1 - \theta_1}
#'
#' The steady-state employment rate is \eqn{\alpha(0,1) + \alpha(1,1)} from
#' \code{stationary_ar2()}.
#'
#' @param params Named list of EM parameters. Must contain \code{theta0},
#'   \code{theta01}, \code{theta1}, \code{theta10}. May contain \code{pi}
#'   (symmetric), \code{pi0} and \code{pi1} (asymmetric), or neither
#'   (no-error variant).
#' @param model_type Character scalar: \code{"symmetric"}, \code{"asymmetric"},
#'   or \code{"none"}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{entry_rate}{Baseline job-finding probability \eqn{\theta_0}.}
#'     \item{exit_rate}{Baseline separation probability \eqn{1 - \theta_1}.}
#'     \item{p_00}{Transition probability \eqn{p_{00} = \theta_0}: entry from
#'       two periods of non-employment.}
#'     \item{p_10}{Transition probability \eqn{p_{10} = \theta_0 + \theta_{01}}:
#'       entry from employed then non-employed.}
#'     \item{p_01}{Transition probability \eqn{p_{01} = 1 - \theta_1 - \theta_{10}}:
#'       retention from non-employed then employed.}
#'     \item{p_11}{Transition probability \eqn{p_{11} = 1 - \theta_1}:
#'       retention from two periods of employment.}
#'     \item{employment_rate}{Steady-state employment rate
#'       \eqn{\alpha(0,1) + \alpha(1,1)} from the stationary distribution.
#'       \code{NA} if stationary distribution is undefined (degenerate params).}
#'     \item{pi}{Symmetric misclassification probability. \code{NA} for
#'       asymmetric and no-error models.}
#'     \item{pi0}{False-positive rate \eqn{\pi_0}. \code{NA} except for
#'       asymmetric models.}
#'     \item{pi1}{False-negative rate \eqn{\pi_1}. \code{NA} except for
#'       asymmetric models.}
#'   }
#'
#' @examples
#' \dontrun{
#' p <- list(theta0=0.10, theta01=0.05, theta1=0.08, theta10=0.03, pi=0.05)
#' implied_ar2(p, "symmetric")
#' }
#' @export
implied_ar2 <- function(params, model_type = "symmetric") {
  stopifnot(
    is.list(params),
    is.numeric(params$theta0),  length(params$theta0)  == 1L, !is.na(params$theta0),
    is.numeric(params$theta01), length(params$theta01) == 1L, !is.na(params$theta01),
    is.numeric(params$theta1),  length(params$theta1)  == 1L, !is.na(params$theta1),
    is.numeric(params$theta10), length(params$theta10) == 1L, !is.na(params$theta10),
    params$theta0  >= 0, params$theta0  <= 1,
    params$theta1  >= 0, params$theta1  <= 1,
    params$theta01 >= 0,  # p_10 = theta0 + theta01 must be a valid probability
    params$theta10 >= 0,  # p_01 = 1 - theta1 - theta10 must be non-negative
    params$theta0 + params$theta01 <= 1,   # p_10 ∈ [0, 1]
    params$theta1 + params$theta10 <= 1,   # p_01 = 1 - theta1 - theta10 ≥ 0
    model_type %in% c("symmetric", "asymmetric", "none")
  )

  theta0  <- params$theta0
  theta01 <- params$theta01
  theta1  <- params$theta1
  theta10 <- params$theta10

  # Four AR(2) transition probabilities (paper eq. 4)
  p_00 <- theta0
  p_10 <- theta0 + theta01
  p_01 <- 1 - theta1 - theta10
  p_11 <- 1 - theta1

  # Baseline (AR(1)-comparable) entry and exit rates
  entry_rate <- theta0
  exit_rate  <- 1 - theta1

  # Steady-state employment rate from stationary distribution
  employment_rate <- tryCatch({
    alpha_stat <- stationary_ar2(theta0, theta01, theta1, theta10)
    # Employment rate = P(h=1) = alpha(0,1) + alpha(1,1)
    unname(alpha_stat["01"] + alpha_stat["11"])
  }, error = function(e) NA_real_)

  # Misclassification quantities
  pi  <- NA_real_
  pi0 <- NA_real_
  pi1 <- NA_real_

  if (model_type == "symmetric") {
    pi <- params$pi %||% NA_real_
    if (!is.na(pi) && (pi < 0 || pi > 1))
      stop(sprintf("implied_ar2: pi out of [0,1]: %g", pi))
  } else if (model_type == "asymmetric") {
    pi0 <- params$pi0 %||% NA_real_
    pi1 <- params$pi1 %||% NA_real_
    if (!is.na(pi0) && (pi0 < 0 || pi0 > 1))
      stop(sprintf("implied_ar2: pi0 out of [0,1]: %g", pi0))
    if (!is.na(pi1) && (pi1 < 0 || pi1 > 1))
      stop(sprintf("implied_ar2: pi1 out of [0,1]: %g", pi1))
  }
  # model_type == "none": all pi fields remain NA.

  list(
    entry_rate       = entry_rate,
    exit_rate        = exit_rate,
    p_00             = p_00,
    p_10             = p_10,
    p_01             = p_01,
    p_11             = p_11,
    employment_rate  = employment_rate,
    pi               = pi,
    pi0              = pi0,
    pi1              = pi1
  )
}
