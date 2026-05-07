# ==============================================================================
# EM-baseline-ext: Implied probability transforms for extension models
# Created: 2026-05-06
#
# Companion to EM-baseline/R/implied_quantities.R. Handles the three extension
# model families:
#   - Covariate models  (probit link, individual-specific theta_i)
#   - FMM models        (two-type mixture)
#   - Inconsistency     (logistic-linked, individual/wave-specific pi)
#
# All functions follow the same contract: take a `params` list (as returned by
# the corresponding EM driver) plus any necessary data matrices, and return a
# named list of implied quantities. Bootstrap wrappers call these inside each
# replicate so SEs cover the full transform chain.
# ==============================================================================

# ---- Internal helpers -------------------------------------------------------

#' Compute individual-specific transition probabilities from probit parameters
#'
#' @param X       N × p numeric covariate matrix (with intercept).
#' @param beta0   p-vector: coefficients for job-finding probit.
#' @param beta1   p-vector: coefficients for employment-persistence probit.
#' @return A list with elements \code{theta0_i} and \code{theta1_i}, each an
#'   N-vector of probabilities in (0, 1).
#' @keywords internal
.probit_thetas <- function(X, beta0, beta1) {
  list(
    theta0_i = pnorm(as.vector(X %*% beta0)),
    theta1_i = pnorm(as.vector(X %*% beta1))
  )
}

#' Steady-state employment rate from individual-specific transition probabilities
#'
#' @param theta0_i N-vector of job-finding rates.
#' @param theta1_i N-vector of employment persistence rates.
#' @return N-vector of implied steady-state employment probabilities.
#' @keywords internal
.stationary_alpha <- function(theta0_i, theta1_i) {
  denom <- theta0_i + (1 - theta1_i)
  # Fast-path: avoid ifelse allocation when denominator is always positive
  if (any(denom <= 0)) ifelse(denom > 0, theta0_i / denom, NA_real_)
  else theta0_i / denom
}

# ---- Extension I: Covariate models ------------------------------------------

#' Compute implied probabilities for a covariate EM model
#'
#' Converts probit-link parameters \eqn{(\beta_0, \beta_1)} into
#' population-average implied probabilities and Average Marginal Effects (AMEs).
#'
#' The probit link is \eqn{\theta_{0i} = \Phi(x_i'\beta_0)}, with individual
#' AMEs \eqn{\partial \theta_{0i} / \partial x_{ik} = \phi(x_i'\beta_0)\beta_{0k}}.
#' The population-average AME averages these over all \eqn{i}.  For the exit
#' rate (\eqn{1 - \theta_{1i}}), the AME is the negative of the AME for
#' \eqn{\theta_{1i}}.
#'
#' @param params     Named list from \code{em_fit_covariates()}: must contain
#'   \code{beta0} (p-vector), \code{beta1} (p-vector), and optionally \code{pi}.
#' @param X          N × p numeric covariate matrix (with intercept in column 1).
#' @param model_type Character scalar: \code{"symmetric"} or \code{"none"}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{mean_entry_rate}{Population mean of \eqn{\Phi(x_i'\beta_0)}.}
#'     \item{mean_exit_rate}{Population mean of \eqn{1 - \Phi(x_i'\beta_1)}.}
#'     \item{mean_employment_rate}{Population mean of steady-state \eqn{\alpha_i}.}
#'     \item{pi}{Symmetric misclassification probability, or \code{NA}.}
#'     \item{ame_entry}{Named numeric vector: AME of each covariate on entry rate.}
#'     \item{ame_exit}{Named numeric vector: AME of each covariate on exit rate
#'       (positive = covariate raises exit rate).}
#'   }
#'
#' @examples
#' set.seed(1)
#' X    <- cbind(1, rnorm(100), rnorm(100))
#' beta <- c(qnorm(0.10), 0.5, -0.3)
#' p    <- list(beta0 = beta, beta1 = c(qnorm(0.90), -0.2, 0.1), pi = 0.05)
#' implied_covariates(p, X, "symmetric")
implied_covariates <- function(params, X, model_type = "symmetric") {
  stopifnot(
    is.list(params),
    is.numeric(params$beta0), !anyNA(params$beta0),
    is.numeric(params$beta1), !anyNA(params$beta1),
    is.matrix(X), is.numeric(X),
    ncol(X) == length(params$beta0),
    ncol(X) == length(params$beta1),
    model_type %in% c("symmetric", "none")
  )
  if (anyNA(X))
    stop("implied_covariates: X contains NA values — impute or filter before calling.")

  # Cache linear predictors to avoid recomputing X %*% beta for both pnorm and dnorm
  xb0 <- as.vector(X %*% params$beta0)
  xb1 <- as.vector(X %*% params$beta1)

  theta0_i <- pnorm(xb0)
  theta1_i <- pnorm(xb1)
  alpha_i  <- .stationary_alpha(theta0_i, theta1_i)

  # Population-average implied probabilities
  mean_entry_rate      <- mean(theta0_i,       na.rm = TRUE)
  mean_exit_rate       <- mean(1 - theta1_i,   na.rm = TRUE)
  mean_employment_rate <- mean(alpha_i,         na.rm = TRUE)

  # Average Marginal Effects (AME)
  # AME_k = (1/N) sum_i phi(x_i' beta) * beta_k = beta_k * mean(phi(x_i' beta))
  # For exit rate (1 - theta1): AME = -AME of theta1 (sign reversal from complement)
  phi0 <- dnorm(xb0)  # reuse cached linear predictor — no extra matrix-vector product
  phi1 <- dnorm(xb1)

  mean_phi0 <- mean(phi0)
  mean_phi1 <- mean(phi1)

  ame_entry <- params$beta0 * mean_phi0  # vectorized: one multiplication per coef
  ame_exit  <- -params$beta1 * mean_phi1 # negative: d(1-theta1)/dX = -d(theta1)/dX

  col_nms <- colnames(X)
  if (is.null(col_nms)) col_nms <- paste0("x", seq_len(ncol(X)))
  names(ame_entry) <- col_nms
  names(ame_exit)  <- col_nms

  pi <- if (model_type == "symmetric") params$pi %||% NA_real_ else NA_real_

  list(
    mean_entry_rate      = mean_entry_rate,
    mean_exit_rate       = mean_exit_rate,
    mean_employment_rate = mean_employment_rate,
    pi                   = pi,
    ame_entry            = ame_entry,
    ame_exit             = ame_exit
  )
}

# ---- Extension III: FMM models ----------------------------------------------

#' Compute implied probabilities for a two-type FMM EM model
#'
#' Converts type-specific Markov parameters and mixing weight \eqn{\phi} into
#' implied probabilities for each type separately and for the overall population
#' (weighted average by \eqn{\phi}).
#'
#' The weighted-average entry and exit rates represent the unconditional
#' population rates under the mixture: a randomly drawn individual faces entry
#' rate \eqn{\phi \theta_0^A + (1-\phi) \theta_0^B} and exit rate
#' \eqn{\phi (1-\theta_1^A) + (1-\phi)(1-\theta_1^B)}.
#'
#' @param params     Named list from \code{em_fit_fmm()}: must contain
#'   \code{theta0_A}, \code{theta1_A}, \code{theta0_B}, \code{theta1_B},
#'   \code{phi}, and optionally \code{pi}.
#' @param model_type Character scalar: \code{"symmetric"} or \code{"none"}.
#'
#' @return A named list:
#'   \describe{
#'     \item{entry_rate_A}{Job-finding rate for type A.}
#'     \item{exit_rate_A}{Separation rate for type A (\eqn{1 - \theta_1^A}).}
#'     \item{employment_rate_A}{Steady-state employment rate for type A.}
#'     \item{entry_rate_B}{Job-finding rate for type B.}
#'     \item{exit_rate_B}{Separation rate for type B.}
#'     \item{employment_rate_B}{Steady-state employment rate for type B.}
#'     \item{phi}{Mixing weight (probability of being type A).}
#'     \item{weighted_entry_rate}{Mixture-weighted average entry rate.}
#'     \item{weighted_exit_rate}{Mixture-weighted average exit rate.}
#'     \item{pi}{Symmetric misclassification probability, or \code{NA}.}
#'   }
#'
#' @examples
#' p <- list(theta0_A = 0.15, theta1_A = 0.95,
#'           theta0_B = 0.05, theta1_B = 0.70,
#'           phi = 0.4, pi = 0.05)
#' implied_fmm(p, "symmetric")
implied_fmm <- function(params, model_type = "symmetric") {
  required <- c("theta0_A", "theta1_A", "theta0_B", "theta1_B", "phi")
  missing_nms <- setdiff(required, names(params))
  if (length(missing_nms) > 0L)
    stop(sprintf("implied_fmm: missing params: %s",
                 paste(missing_nms, collapse = ", ")))
  stopifnot(model_type %in% c("symmetric", "none"))

  # Validate all probability parameters: must be numeric scalars in [0, 1], not NA
  for (nm in c("theta0_A", "theta1_A", "theta0_B", "theta1_B", "phi")) {
    v <- params[[nm]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1)
      stop(sprintf("implied_fmm: '%s' must be a numeric scalar in [0,1]; got: %s",
                   nm, deparse(v)))
  }

  theta0_A <- params$theta0_A
  theta1_A <- params$theta1_A
  theta0_B <- params$theta0_B
  theta1_B <- params$theta1_B
  phi      <- params$phi

  denom_A <- theta0_A + (1 - theta1_A)
  denom_B <- theta0_B + (1 - theta1_B)

  employment_rate_A <- if (denom_A > 0) theta0_A / denom_A else NA_real_
  employment_rate_B <- if (denom_B > 0) theta0_B / denom_B else NA_real_

  weighted_entry_rate <- phi * theta0_A + (1 - phi) * theta0_B
  weighted_exit_rate  <- phi * (1 - theta1_A) + (1 - phi) * (1 - theta1_B)

  pi <- if (model_type == "symmetric") params$pi %||% NA_real_ else NA_real_

  list(
    entry_rate_A         = theta0_A,
    exit_rate_A          = 1 - theta1_A,
    employment_rate_A    = employment_rate_A,
    entry_rate_B         = theta0_B,
    exit_rate_B          = 1 - theta1_B,
    employment_rate_B    = employment_rate_B,
    phi                  = phi,
    weighted_entry_rate  = weighted_entry_rate,
    weighted_exit_rate   = weighted_exit_rate,
    pi                   = pi
  )
}

# ---- Extension IV: Inconsistency models -------------------------------------

#' Compute implied probabilities for an inconsistency-augmented EM model
#'
#' The inconsistency model parameterises misclassification as individual- and
#' wave-specific via a logistic link:
#' \deqn{\pi_{it} = \frac{1}{2} \sigma(\delta_0 + \delta_1 Y_{it}^{\text{age}}
#'   + \delta_2 Y_{it}^{\text{edu}})}
#' where \eqn{\sigma} is the logistic function and \eqn{Y_{it}} are binary
#' inconsistency indicators.  The mean misclassification rate averages these
#' over all individual-wave observations.
#'
#' @param params      Named list from \code{em_fit_inconsistency()}: must
#'   contain \code{theta0}, \code{theta1}, \code{delta} (length-3 vector).
#' @param incons_mat  N × 6 integer matrix of inconsistency indicators.  Columns
#'   must be ordered: \code{Y_age_1, Y_age_2, Y_age_3, Y_edu_1, Y_edu_2,
#'   Y_edu_3}.
#'
#' @return A named list:
#'   \describe{
#'     \item{entry_rate}{Job-finding probability per period (\eqn{\theta_0}).}
#'     \item{exit_rate}{Separation probability (\eqn{1 - \theta_1}).}
#'     \item{employment_rate}{Steady-state employment rate.}
#'     \item{mean_pi}{Population-average misclassification rate, averaged over
#'       all N × 3 individual-wave observations.}
#'     \item{delta}{Length-3 delta coefficient vector (intercept, age, educ).}
#'     \item{pi_base}{Baseline misclassification rate (no inconsistency flags),
#'       \eqn{\frac{1}{2}\sigma(\delta_0)}.}
#'   }
#'
#' @examples
#' set.seed(1)
#' p     <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
#'               delta = c(-2.2, 0.5, 0.3))
#' imat  <- matrix(rbinom(300, 1, 0.1), ncol = 6)
#' implied_inconsistency(p, imat)
implied_inconsistency <- function(params, incons_mat) {
  stopifnot(
    is.list(params),
    is.numeric(params$theta0), length(params$theta0) == 1L, !is.na(params$theta0),
    params$theta0 >= 0, params$theta0 <= 1,
    is.numeric(params$theta1), length(params$theta1) == 1L, !is.na(params$theta1),
    params$theta1 >= 0, params$theta1 <= 1,
    is.numeric(params$delta),  length(params$delta)  == 3L,
    is.matrix(incons_mat), nrow(incons_mat) > 0L, ncol(incons_mat) == 6L
  )
  if (anyNA(incons_mat))
    stop("implied_inconsistency: incons_mat contains NA values.")
  if (!all(incons_mat %in% c(0L, 1L)))
    stop("implied_inconsistency: incons_mat values must be 0 or 1.")

  theta0 <- params$theta0
  theta1 <- params$theta1
  delta  <- params$delta

  denom <- theta0 + (1 - theta1)
  employment_rate <- if (denom > 0) theta0 / denom else NA_real_

  # For each individual-wave (i, t): pi_it = 0.5 * sigma(delta0 + delta1*Y_age_it + delta2*Y_edu_it)
  # Vectorised: compute all N*3 values at once using matrix column slices.
  # Columns 1-3 are age flags, columns 4-6 are edu flags.
  linpred_mat <- delta[1L] +
                 delta[2L] * incons_mat[, 1:3, drop = FALSE] +
                 delta[3L] * incons_mat[, 4:6, drop = FALSE]
  mean_pi <- mean(0.5 * plogis(linpred_mat))

  pi_base <- 0.5 * plogis(delta[1L])   # misclassification at zero flags

  list(
    entry_rate      = theta0,
    exit_rate       = 1 - theta1,
    employment_rate = employment_rate,
    mean_pi         = mean_pi,
    delta           = delta,
    pi_base         = pi_base
  )
}
