# ==============================================================================
# EM-baseline: Probability transforms
# Created: 2026-05-05
#
# Logit and inverse-logit transforms for parameter space mapping.
# No CTMC link is needed for the baseline model (no duration data).
# ==============================================================================

#' Logit transformation
#'
#' Maps a probability in (0, 1) to the real line. Values of exactly 0 or 1
#' produce -Inf or Inf respectively; values outside [0, 1] produce NaN.
#'
#' @param p Probability in (0, 1).
#' @return log(p / (1 - p)).
#' @examples
#' logit(0.5)  # 0
#' logit(0.9)  # approx 2.197
#' @export
logit <- function(p) log(p / (1 - p))

#' Inverse logit (sigmoid) transformation
#'
#' Maps a real number to a probability in (0, 1).
#'
#' @param x Real number.
#' @return 1 / (1 + exp(-x)).
#' @examples
#' inv_logit(0)  # 0.5
#' inv_logit(2)  # approx 0.88
#' @export
inv_logit <- function(x) 1 / (1 + exp(-x))
