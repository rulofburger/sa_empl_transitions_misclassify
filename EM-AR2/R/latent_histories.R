# ==============================================================================
# EM-AR2: Latent histories, stationary distribution, and Markov prior
# Created: 2026-05-05
# ==============================================================================
# The AR(2) model has 2^4 = 16 latent employment histories h = (h1,h2,h3,h4)
# where h_t in {0,1} indicates true employment at wave t.
#
# The AR(2) transition probabilities are (TeX: Eq 2):
#   P(h_t=1 | h_{t-2}=j, h_{t-1}=k):
#     (0,0)->1: theta0
#     (1,0)->1: theta0 + theta01   [duration-dependence: had job 2 waves ago]
#     (1,1)->0: theta1
#     (0,1)->0: theta1 + theta10   [duration-dependence: just entered employment]
#
# The stationary distribution alpha(h1,h2) satisfies the balance equation
# with normalisation constant Phi = theta0*(theta10-1) + theta1*(theta01-1).
# (TeX: Eqs 3-4)
# ==============================================================================

#' Generate all 16 latent employment histories for the AR(2) 4-wave model
#'
#' Each row is a history h = (h1, h2, h3, h4) where h_t in {0,1}
#' indicates true employment at wave t.
#'
#' @return A 16 x 4 integer matrix. Row ordering follows
#'   expand.grid(h1=0:1, h2=0:1, h3=0:1, h4=0:1).
#' @export
latent_histories_ar2 <- function() {
  as.matrix(expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1, h4 = 0:1))
}

#' AR(2) transition probability P(h_t=1 | h_{t-2}=j, h_{t-1}=k)
#'
#' Returns the probability of being employed at t, given states at t-2 and t-1.
#' The AR(2) parameterisation is:
#'   p(1|0,0) = theta0
#'   p(1|1,0) = theta0 + theta01
#'   p(1|0,1) = 1 - theta1 - theta10
#'   p(1|1,1) = 1 - theta1
#'
#' @param h_prev2 Integer vector: states at t-2 (values 0 or 1).
#' @param h_prev1 Integer vector: states at t-1 (values 0 or 1). Same length as h_prev2.
#' @param theta0 Baseline job-finding probability.
#' @param theta01 Duration-dependence increment for job-finding.
#' @param theta1 Baseline separation probability.
#' @param theta10 Duration-dependence increment for separation.
#' @return Numeric vector: P(h_t=1 | h_{t-2}, h_{t-1}).
#' @export
ar2_trans_prob <- function(h_prev2, h_prev1, theta0, theta01, theta1, theta10) {
  # Vectorised: all four cases
  # p(1|0,0) = theta0
  # p(1|1,0) = theta0 + theta01
  # p(1|0,1) = 1 - theta1 - theta10
  # p(1|1,1) = 1 - theta1
  ifelse(
    h_prev1 == 0,
    ifelse(h_prev2 == 0, theta0, theta0 + theta01),
    ifelse(h_prev2 == 0, 1 - theta1 - theta10, 1 - theta1)
  )
}

#' Compute the AR(2) stationary joint distribution alpha(h1, h2)
#'
#' The stationary distribution over (h1, h2) pairs satisfies the balance
#' equations for the second-order Markov chain. With normalisation constant
#' Phi = theta0*(theta10-1) + theta1*(theta01-1), the four joint probabilities
#' are (TeX: Eqs 3-4):
#'   alpha(0,0) = theta1*(theta0+theta01-1) / Phi
#'   alpha(1,0) = -theta1*theta0 / Phi
#'   alpha(0,1) = -theta0*theta1 / Phi   [= alpha(1,0) by symmetry]
#'   alpha(1,1) = theta0*(theta1+theta10-1) / Phi
#'
#' @param theta0 Baseline job-finding probability in (0,1).
#' @param theta01 Duration-dependence increment for job-finding.
#' @param theta1 Baseline separation probability in (0,1).
#' @param theta10 Duration-dependence increment for separation.
#' @return Named numeric vector of length 4 with names "00", "10", "01", "11".
#'   Sums to 1.
#' @export
stationary_ar2 <- function(theta0, theta01, theta1, theta10) {
  theta0  <- .bound01(theta0)
  theta01 <- .bound01(theta01)
  theta1  <- .bound01(theta1)
  theta10 <- .bound01(theta10)

  Phi <- theta0 * (theta10 - 1) + theta1 * (theta01 - 1)

  if (abs(Phi) < 1e-12) {
    stop("Stationary distribution undefined: Phi is near zero. ",
         "Check that transition probabilities are not degenerate.")
  }

  alpha_00 <- theta1 * (theta0 + theta01 - 1) / Phi
  alpha_10 <- -theta1 * theta0 / Phi
  alpha_01 <- -theta0 * theta1 / Phi   # always equals alpha_10 for this AR(2) parameterisation
  alpha_11 <- theta0 * (theta1 + theta10 - 1) / Phi

  alpha <- c(`00` = alpha_00, `10` = alpha_10, `01` = alpha_01, `11` = alpha_11)

  # Guard against negative components — indicates violated model constraints
  if (any(alpha < -1e-9)) {
    stop("stationary_ar2: negative stationary probabilities. ",
         "Check that theta0+theta01 < 1 and theta1+theta10 < 1.")
  }

  # Normalise to sum to 1 (guards against floating-point accumulation)
  alpha / sum(alpha)
}

#' Compute the AR(2) Markov prior over all 16 latent histories
#'
#' P(H=h) = alpha(h1,h2) * P(h3|h1,h2) * P(h4|h2,h3)
#'
#' During EM estimation, `alpha` is passed as a free parameter (a named
#' length-4 vector with names "00","10","01","11") that is updated each
#' M-step. This is the Baum-Welch approach and guarantees non-decreasing LL.
#'
#' If `alpha` is NULL, the stationary distribution is computed from the
#' transition parameters via `stationary_ar2()`. Use this for initialization
#' only (before the EM starts) or for inference after convergence.
#'
#' @param hmat 16 x 4 integer matrix of latent histories (from latent_histories_ar2()).
#' @param theta0 Baseline job-finding probability.
#' @param theta01 Duration-dependence increment for job-finding.
#' @param theta1 Baseline separation probability.
#' @param theta10 Duration-dependence increment for separation.
#' @param alpha Optional named length-4 numeric vector: free initial pair
#'   distribution with names "00","10","01","11". If NULL (default), the
#'   stationary distribution is used.
#' @return Length-16 numeric vector of prior probabilities summing to 1.
#' @export
prior_over_histories_ar2 <- function(hmat, theta0, theta01, theta1, theta10,
                                      alpha = NULL) {
  if (is.null(alpha)) {
    alpha <- stationary_ar2(theta0, theta01, theta1, theta10)
  } else {
    expected_names <- c("00", "10", "01", "11")
    if (!setequal(names(alpha), expected_names)) {
      stop("prior_over_histories_ar2: alpha must be a named vector with names ",
           "'00','10','01','11'. Got: ", paste(names(alpha), collapse=","))
    }
    if (any(alpha < 0)) {
      stop("prior_over_histories_ar2: alpha contains negative values.")
    }
  }

  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]
  h4 <- hmat[, 4]

  # Prior over first pair: alpha(h1, h2)
  pair12 <- paste0(h1, h2)
  p_init <- alpha[pair12]

  # Transition wave 2->3 depends on (h1, h2)
  p_h3_1 <- ar2_trans_prob(h1, h2, theta0, theta01, theta1, theta10)
  p_trans23 <- ifelse(h3 == 1, p_h3_1, 1 - p_h3_1)

  # Transition wave 3->4 depends on (h2, h3)
  p_h4_1 <- ar2_trans_prob(h2, h3, theta0, theta01, theta1, theta10)
  p_trans34 <- ifelse(h4 == 1, p_h4_1, 1 - p_h4_1)

  prior <- p_init * p_trans23 * p_trans34
  prior / sum(prior)
}
