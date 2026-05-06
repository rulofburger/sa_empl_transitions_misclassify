# ==============================================================================
# EM-baseline: Latent histories and Markov prior
# Created: 2026-05-05
#
# The model has 8 latent employment histories h = (h1, h2, h3) in {0,1}^3.
# The prior over histories follows a first-order Markov chain.
#
# TeX reference: EM baseline.tex, Eq (6) for the prior.
# ==============================================================================

#' Generate all 8 latent employment histories
#'
#' Each row is a history h = (h1, h2, h3) where h_t in {0, 1} indicates
#' true employment at wave t (1 = employed, 0 = nonemployed).
#'
#' @return An 8 x 3 integer matrix. Row ordering follows
#'   \code{expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1)}.
#' @examples
#' hmat <- latent_histories()
#' nrow(hmat)  # 8
#' @export
latent_histories <- function() {
  as.matrix(expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1))[, 1:3]
}

#' Compute Markov prior over latent histories
#'
#' P(h | alpha, theta0, theta1) = P(h1) * P(h2|h1) * P(h3|h2), where:
#' \itemize{
#'   \item P(h1 = 1) = alpha
#'   \item P(h_t = 1 | h_{t-1} = 1) = theta1  (persistence / stay-employed)
#'   \item P(h_t = 1 | h_{t-1} = 0) = theta0  (job-finding rate)
#' }
#'
#' TeX ref: Eq (6)
#'
#' @param hmat 8 x 3 integer matrix of latent histories (from
#'   \code{\link{latent_histories}}).
#' @param theta1 Employment persistence probability P(1->1) in (0, 1).
#' @param theta0 Job-finding probability P(0->1) in (0, 1).
#' @param alpha Initial employment probability P(h1 = 1) in (0, 1).
#' @return Length-8 numeric vector of prior probabilities summing to 1.
#' @examples
#' hmat <- latent_histories()
#' p <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 0.5)
#' sum(p)  # 1
#' @export
prior_over_histories <- function(hmat, theta1, theta0, alpha) {
  theta1 <- .bound01(theta1)
  theta0 <- .bound01(theta0)
  alpha  <- .bound01(alpha)

  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]

  # Initial state probability
  p_init <- ifelse(h1 == 1, alpha, 1 - alpha)

  # Transition probability: wave 1 -> wave 2
  # Vectorised indexing avoids nested ifelse() for readability
  p12 <- numeric(8L)
  p12[h1 == 1 & h2 == 1] <- theta1
  p12[h1 == 1 & h2 == 0] <- 1 - theta1
  p12[h1 == 0 & h2 == 1] <- theta0
  p12[h1 == 0 & h2 == 0] <- 1 - theta0

  # Transition probability: wave 2 -> wave 3
  p23 <- numeric(8L)
  p23[h2 == 1 & h3 == 1] <- theta1
  p23[h2 == 1 & h3 == 0] <- 1 - theta1
  p23[h2 == 0 & h3 == 1] <- theta0
  p23[h2 == 0 & h3 == 0] <- 1 - theta0

  p <- p_init * p12 * p23
  p / sum(p)
}

#' Compute stationary initial employment probability
#'
#' Under the stationarity restriction, the Markov chain is in its ergodic
#' steady state at wave 1, so alpha = theta0 / (theta0 + 1 - theta1).
#'
#' TeX ref: EM baseline.tex, Eq (4)
#'
#' @param theta0 Job-finding probability in (0, 1).
#' @param theta1 Employment persistence in (0, 1).
#' @return Scalar alpha in (0, 1).
#' @examples
#' stationary_alpha(theta0 = 0.1, theta1 = 0.9)  # 0.5
#' @export
stationary_alpha <- function(theta0, theta1) {
  theta0 <- .bound01(theta0)
  theta1 <- .bound01(theta1)
  .bound01(theta0 / (theta0 + 1 - theta1))
}
