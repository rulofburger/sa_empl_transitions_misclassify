# ==============================================================================
# EM-tenure: Probability and CTMC transforms
# ==============================================================================
# Logit/inverse-logit and the continuous-time Markov chain (CTMC) link
# between discrete transition probabilities and exponential duration rates.
#
# CTMC link (two-state model, inter-wave spacing Delta = 0.25 years):
#
#   Employment (state 1):
#     theta1 = P(stay employed) = exp(-lambda_g * Delta)
#     lambda_g = -log(theta1) / Delta                    ... (Eq 3a)
#
#   Nonemployment (state 0):
#     theta0 = P(exit nonemployment) = 1 - exp(-lambda_d * Delta)
#     lambda_d = -log(1 - theta0) / Delta                ... (Eq 3b)
#
# All rates (lambda_g, lambda_d) are PER YEAR, matching emission durations
# which are measured in years.
# ==============================================================================

#' Logit transformation
#'
#' @param p Probability in (0, 1).
#' @return log(p / (1 - p)).
#' @export
logit <- function(p) log(p / (1 - p))

#' Inverse logit (sigmoid) transformation
#'
#' @param x Real number.
#' @return 1 / (1 + exp(-x)).
#' @export
inv_logit <- function(x) 1 / (1 + exp(-x))

#' Convert discrete-time transition probability to continuous-time rate
#'
#' For a transition (exit) probability: theta = 1 - exp(-lambda * delta),
#' so lambda = -log(1 - theta) / delta.
#'
#' Use this for theta0 (job-finding probability = exit from nonemployment).
#'
#' @param theta Discrete transition probability in (0, 1).
#' @param delta Spacing between waves in years (default 0.25 = one quarter).
#' @return Continuous-time rate lambda (per year).
#' @export
ctmc_lambda_from_transition <- function(theta, delta = .QUARTER_YEARS) {
  -log(1 - .bound01(theta)) / delta
}

#' Convert discrete-time persistence probability to continuous-time rate
#'
#' For a persistence (staying) probability: theta = exp(-lambda * delta),
#' so lambda = -log(theta) / delta.
#'
#' Use this for theta1 (employment persistence probability).
#'
#' @param theta Discrete persistence probability in (0, 1).
#' @param delta Spacing between waves in years (default 0.25 = one quarter).
#' @return Continuous-time rate lambda (per year).
#' @export
ctmc_lambda_from_persistence <- function(theta, delta = .QUARTER_YEARS) {
  -log(.bound01(theta)) / delta
}

#' Convert continuous-time rate to discrete-time transition probability
#'
#' theta = 1 - exp(-lambda * delta).
#'
#' Inverse of \code{ctmc_lambda_from_transition}.
#'
#' @param lambda Continuous-time rate (> 0, per year).
#' @param delta Spacing between waves in years (default 0.25).
#' @return Discrete transition probability in (0, 1).
#' @export
ctmc_transition_from_lambda <- function(lambda, delta = .QUARTER_YEARS) {
  1 - exp(-.pos(lambda) * delta)
}

#' Convert continuous-time rate to discrete-time persistence probability
#'
#' theta = exp(-lambda * delta).
#'
#' Inverse of \code{ctmc_lambda_from_persistence}.
#'
#' @param lambda Continuous-time rate (> 0, per year).
#' @param delta Spacing between waves in years (default 0.25).
#' @return Discrete persistence probability in (0, 1).
#' @export
ctmc_persistence_from_lambda <- function(lambda, delta = .QUARTER_YEARS) {
  exp(-.pos(lambda) * delta)
}

# --- Deprecated wrappers (will be removed in a future version) ---------------
# These existed before the theta1/theta0 distinction was corrected.
# They treated ALL thetas as transition probabilities with delta in months.
# Callers should migrate to ctmc_lambda_from_persistence (theta1) or
# ctmc_lambda_from_transition (theta0).

#' @rdname ctmc_lambda_from_transition
#' @usage NULL
#' @export
ctmc_lambda_from_theta <- function(theta, delta = .QUARTER_YEARS) {
  .Deprecated("ctmc_lambda_from_transition or ctmc_lambda_from_persistence")
  ctmc_lambda_from_transition(theta, delta)
}

#' @rdname ctmc_transition_from_lambda
#' @usage NULL
#' @export
ctmc_theta_from_lambda <- function(lambda, delta = .QUARTER_YEARS) {
  .Deprecated("ctmc_transition_from_lambda or ctmc_persistence_from_lambda")
  ctmc_transition_from_lambda(lambda, delta)
}
