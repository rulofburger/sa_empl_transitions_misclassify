# ==============================================================================
# EM-tenure: Probability and CTMC transforms
# ==============================================================================
# Logit/inverse-logit and the continuous-time Markov chain (CTMC) link
# between discrete transition probabilities and exponential duration rates.
#
# TeX reference: Eq (3) -- lambda_g = -log(1-theta_1)/3, lambda_d = -log(1-theta_0)/3
# where Delta = 3 months is the inter-wave spacing.
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
#' Under a two-state CTMC with inter-wave spacing Delta (in months),
#' the discrete transition probability theta and continuous rate lambda
#' satisfy: theta = 1 - exp(-lambda * Delta).
#'
#' Inverting: lambda = -log(1 - theta) / Delta.
#'
#' TeX ref: Eq (3), lambda_g = -log(1-theta_1)/3
#'
#' @param theta Discrete transition probability in (0, 1).
#' @param delta Spacing between waves in months (default 3).
#' @return Continuous-time rate lambda (per month).
#' @export
ctmc_lambda_from_theta <- function(theta, delta = 3) {
  -log(1 - .bound01(theta)) / delta
}

#' Convert continuous-time rate to discrete-time transition probability
#'
#' theta = 1 - exp(-lambda * delta).
#'
#' @param lambda Continuous-time rate (> 0).
#' @param delta Spacing between waves in months (default 3).
#' @return Discrete transition probability in (0, 1).
#' @export
ctmc_theta_from_lambda <- function(lambda, delta = 3) {
  1 - exp(-.pos(lambda) * delta)
}
