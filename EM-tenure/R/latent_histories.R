# ==============================================================================
# EM-tenure: Latent histories, deterministic clocks, and Markov prior
# ==============================================================================
# The model has 8 latent employment histories h = (h1, h2, h3) in {0,1}^3.
# Given a history, tenure and nonemployment clocks evolve deterministically.
# The prior over histories follows a first-order Markov chain.
#
# TeX reference: Eq (1) for prior, Eq (4) for deterministic clocks.
# ==============================================================================

#' Generate all 8 latent employment histories
#'
#' Each row is a history h = (h1, h2, h3) where h_t in {0, 1}
#' indicates true employment at wave t.
#'
#' @return An 8 x 3 integer matrix. Row ordering follows
#'   expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1).
#' @export
latent_histories <- function() {
  as.matrix(expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1))[, 1:3]
}

#' Compute deterministic duration clocks from latent histories
#'
#' Given h, the true tenure g*(h) and nonemployment duration d*(h)
#' evolve on a quarterly grid:
#'   - Continuation (h_{t-1} = h_t = 1): g*_t = g*_{t-1} + increment
#'   - Start (h_{t-1} = 0, h_t = 1): g*_t = increment (clock resets)
#'   - Otherwise: g*_t = 0
#' Analogously for nonemployment.
#'
#' Wave 1: g*_1 = increment if h_1=1, else 0; d*_1 = increment if h_1=0, else 0.
#' (Wave 1 true duration is a random variable, not deterministic; see EMG.)
#'
#' TeX ref: Eq (4)
#'
#' @param hmat 8 x 3 matrix of latent histories.
#' @param increment Quarterly time increment in years (default .QUARTER_YEARS).
#' @return List with two 8 x 3 matrices: Gstar (tenure), Dstar (nonemployment).
#' @export
clocks_from_histories <- function(hmat, increment = .QUARTER_YEARS) {
  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]

  # Wave 1
  g1 <- h1 * increment
  d1 <- (1 - h1) * increment

  # Wave 2: continuation adds increment, start resets to increment, else 0
  g2 <- ifelse(h2 == 1, ifelse(h1 == 1, g1 + increment, increment), 0)
  d2 <- ifelse(h2 == 0, ifelse(h1 == 0, d1 + increment, increment), 0)

  # Wave 3
  g3 <- ifelse(h3 == 1, ifelse(h2 == 1, g2 + increment, increment), 0)
  d3 <- ifelse(h3 == 0, ifelse(h2 == 0, d2 + increment, increment), 0)

  list(
    Gstar = cbind(g1, g2, g3),
    Dstar = cbind(d1, d2, d3)
  )
}

#' Compute Markov prior over latent histories
#'
#' P(h | alpha, theta0, theta1) = P(h1) * P(h2|h1) * P(h3|h2)
#' where:
#'   P(h1 = 1) = alpha
#'   P(h_t = 1 | h_{t-1} = 1) = theta1
#'   P(h_t = 1 | h_{t-1} = 0) = theta0
#'
#' TeX ref: Eq (1)
#'
#' @param hmat 8 x 3 matrix of latent histories.
#' @param theta1 Employment persistence P(1->1).
#' @param theta0 Job finding P(0->1).
#' @param alpha Initial employment probability P(h1 = 1).
#' @return Length-8 vector of probabilities summing to 1.
#' @export
prior_over_histories <- function(hmat, theta1, theta0, alpha) {
  theta1 <- .bound01(theta1)
  theta0 <- .bound01(theta0)
  alpha  <- .bound01(alpha)

  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]

  p_init <- ifelse(h1 == 1, alpha, 1 - alpha)

  p12 <- ifelse(h1 == 1,
                ifelse(h2 == 1, theta1, 1 - theta1),
                ifelse(h2 == 1, theta0, 1 - theta0))

  p23 <- ifelse(h2 == 1,
                ifelse(h3 == 1, theta1, 1 - theta1),
                ifelse(h3 == 1, theta0, 1 - theta0))

  p <- p_init * p12 * p23
  p / sum(p)
}
