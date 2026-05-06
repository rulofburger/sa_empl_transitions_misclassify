# ==============================================================================
# EM-AR2: E-step — responsibilities and sufficient statistics
# Created: 2026-05-05
# ==============================================================================
# Computes gamma_{ih} = P(H=h | y_i, Theta) for each observation i and
# latent history h (16 histories for a 4-wave panel), plus all weighted
# sufficient statistics needed by the M-step.
#
# MODEL:
#   Emission: P(y_t=k | h_t, pi) = (1-pi) if k==h_t, else pi   [symmetric]
#   Asymmetric: P(y_t=0|h_t=1) = pi0,  P(y_t=1|h_t=0) = pi1
#
# SUFFICIENT STATISTICS:
#   D_{jk}   = sum_i sum_h gamma_{ih} * sum_{t in {3,4}} 1(h_{t-2}=j, h_{t-1}=k)
#   T_{jk1}  = sum_i sum_h gamma_{ih} * sum_{t in {3,4}} 1(h_{t-2}=j, h_{t-1}=k, h_t=1)
#   M        = sum_i sum_h gamma_{ih} * sum_t 1(y_{it} != h_t)       [symmetric]
#   M0, M1   = separate false-positive/false-negative counts          [asymmetric]
#   W        = sum_i w_i  (total weight)
#
# VECTORISATION: Operations are fully vectorised over observations (N). The
# only explicit loop is over H=16 latent histories (fixed size).
#
# TeX reference: Section on E-step and M-step sufficient statistics.
# ==============================================================================

#' Compute emission log-probabilities for all observations and histories
#'
#' For each observation i and history h, computes
#' log P(y_i | H=h, pi) = sum_t log P(y_{it} | h_t, pi).
#'
#' @param ymat N x 4 integer matrix of observed employment (0/1).
#' @param hmat 16 x 4 integer matrix of latent histories.
#' @param pi Scalar symmetric misclassification probability.
#' @param pi0 Scalar false non-employment probability P(y=0|h=1). Ignored if
#'   asymmetric=FALSE.
#' @param pi1 Scalar false employment probability P(y=1|h=0). Ignored if
#'   asymmetric=FALSE.
#' @param asymmetric Logical; if TRUE use pi0/pi1 instead of pi.
#' @return N x 16 matrix of emission log-probabilities.
#' @keywords internal
.emission_logprob <- function(ymat, hmat, pi = 0.05,
                               pi0 = NULL, pi1 = NULL,
                               asymmetric = FALSE) {
  N <- nrow(ymat)
  H <- nrow(hmat)  # 16

  # Log emission probabilities:
  # Correct classification: log(1 - pi)
  # Misclassification:      log(pi)
  if (!asymmetric) {
    log_correct <- log(1 - .bound01(pi))
    log_wrong   <- log(.bound01(pi))
  } else {
    # pi0 = P(y=0|h=1), pi1 = P(y=1|h=0)
    log_correct_emp   <- log(1 - .bound01(pi0))  # y=1, h=1
    log_wrong_emp     <- log(.bound01(pi0))       # y=0, h=1
    log_correct_nonemp <- log(1 - .bound01(pi1))  # y=0, h=0
    log_wrong_nonemp   <- log(.bound01(pi1))       # y=1, h=0
  }

  # Build N x H matrix of total emission log-densities
  log_emit <- matrix(0, nrow = N, ncol = H)

  for (j in seq_len(H)) {
    h <- hmat[j, ]  # length-4 binary vector

    # For each wave t, add log P(y_{it} | h_t, params)
    for (t in 1:4) {
      h_t <- h[t]
      y_t <- ymat[, t]  # length-N vector

      if (!asymmetric) {
        # match = y_t == h_t
        log_emit[, j] <- log_emit[, j] +
          ifelse(y_t == h_t, log_correct, log_wrong)
      } else {
        # Separate false-positive and false-negative
        if (h_t == 1) {
          log_emit[, j] <- log_emit[, j] +
            ifelse(y_t == 1, log_correct_emp, log_wrong_emp)
        } else {
          log_emit[, j] <- log_emit[, j] +
            ifelse(y_t == 0, log_correct_nonemp, log_wrong_nonemp)
        }
      }
    }
  }

  log_emit
}


#' E-step: compute responsibilities and sufficient statistics
#'
#' Computes gamma_{ih} = P(H=h | y_i, Theta) for all i,h and the weighted
#' sufficient statistics D_{jk}, T_{jk1}, M (and M0, M1 for asymmetric).
#'
#' @param ymat N x 4 integer matrix of observed employment (y1-y4).
#' @param w Length-N numeric weight vector.
#' @param hmat 16 x 4 integer matrix from latent_histories_ar2().
#' @param params Named list with:
#'   - theta0, theta01, theta1, theta10: AR(2) transition parameters
#'   - pi: misclassification probability (symmetric)
#'   - pi0, pi1: (optional) for asymmetric model
#'   - asymmetric: logical flag (default FALSE)
#' @return Named list:
#'   - gamma: N x 16 matrix of responsibilities (rows sum to 1)
#'   - loglik: scalar weighted observed-data log-likelihood
#'   - sufficient_stats: list with D (4x2x2 or named), T1 (named), M, M0, M1, W
#' @export
e_step_ar2 <- function(ymat, w, hmat, params) {
  N <- nrow(ymat)
  H <- nrow(hmat)  # 16

  asymmetric <- isTRUE(params$asymmetric)
  pi   <- params$pi   %||% 0.05
  pi0  <- params$pi0  %||% pi
  pi1  <- params$pi1  %||% pi

  # Prior over histories: length-16 vector
  # Use params$alpha as the free initial-pair distribution if provided.
  # This is the Baum-Welch approach; alpha is updated each M-step.
  prior <- prior_over_histories_ar2(hmat,
    theta0  = params$theta0,
    theta01 = params$theta01,
    theta1  = params$theta1,
    theta10 = params$theta10,
    alpha   = params$alpha   # NULL on first call → uses stationary distribution
  )
  log_prior <- log(prior)  # length-16; -Inf for zero-probability histories

  # Emission log-probabilities: N x 16
  log_emit <- .emission_logprob(ymat, hmat,
    pi = pi, pi0 = pi0, pi1 = pi1,
    asymmetric = asymmetric
  )

  # Joint log-probabilities: N x 16  (broadcast prior over rows)
  log_joint <- sweep(log_emit, 2, log_prior, "+")

  # Normalise per observation via log-sum-exp
  log_norm <- .row_logsumexp(log_joint)   # length-N
  log_gamma <- log_joint - log_norm       # N x 16
  gamma <- exp(log_gamma)                 # N x 16 responsibilities

  # Weighted log-likelihood
  loglik <- sum(w * log_norm)

  # --- Compute sufficient statistics ---
  # Weighted responsibilities: N x 16
  wgamma <- gamma * w  # each row scaled by weight_i

  # History columns
  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]
  h4 <- hmat[, 4]

  # Pre-compute column sums (sum over observations) once
  col_wgamma <- colSums(wgamma)  # length-16

  # S_{jk} = sufficient stats for free initial pair distribution alpha(h1,h2)
  # S_{jk} = sum_i sum_h gamma_{ih} * w_i * 1(h1=j, h2=k)
  S <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("0","1"), c("0","1")))
  for (j in 0:1) {
    for (k in 0:1) {
      ind_init <- (h1 == j) & (h2 == k)  # length-16
      S[j+1, k+1] <- sum(col_wgamma * ind_init)
    }
  }

  # D_{jk} = weighted exposure to AR(2) transitions at waves t=3 and t=4
  # t=3: conditioning on (h1, h2)
  # t=4: conditioning on (h2, h3)
  # For each (j,k) pair: D_{jk} = sum_i sum_h wgamma_{ih} * (1{h1=j,h2=k} + 1{h2=j,h3=k})
  D <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("0","1"), c("0","1")))
  T1 <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("0","1"), c("0","1")))

  for (j in 0:1) {
    for (k in 0:1) {
      # Indicators for (h_{t-2}=j, h_{t-1}=k) at t=3 and t=4
      ind_t3 <- (h1 == j) & (h2 == k)  # length-16
      ind_t4 <- (h2 == j) & (h3 == k)  # length-16

      # T1: transition to employment (h_t=1)
      ind_t3_to1 <- ind_t3 & (h3 == 1)
      ind_t4_to1 <- ind_t4 & (h4 == 1)

      D[j+1, k+1]  <- sum(col_wgamma * (ind_t3 + ind_t4))
      T1[j+1, k+1] <- sum(col_wgamma * (ind_t3_to1 + ind_t4_to1))
    }
  }

  W <- sum(w)

  # Misclassification count
  if (!asymmetric) {
    # M = sum_i sum_h gamma_{ih} * w_i * sum_t 1(y_{it} != h_t)
    M <- 0
    for (t in 1:4) {
      h_t <- hmat[, t]   # length-16
      y_t <- ymat[, t]   # length-N
      # indicator matrix N x 16: 1 if y_{it} != h_t
      mismatch <- outer(y_t, h_t, "!=") * 1L
      M <- M + sum(wgamma * mismatch)
    }
    M0 <- NA_real_
    M1 <- NA_real_
  } else {
    # M0 = false non-employment: y=0, h=1
    # M1 = false employment:     y=1, h=0
    M  <- NA_real_
    M0 <- 0
    M1 <- 0
    # Exposures for asymmetric denominators
    E0 <- 0  # sum_i sum_h gamma*w * sum_t 1(h_t=1)  [exposure for pi0]
    E1 <- 0  # sum_i sum_h gamma*w * sum_t 1(h_t=0)  [exposure for pi1]
    for (t in 1:4) {
      h_t <- hmat[, t]
      y_t <- ymat[, t]
      # False non-employment: h=1, y=0
      fn <- outer(y_t == 0, h_t == 1, "&") * 1L
      M0 <- M0 + sum(wgamma * fn)
      E0 <- E0 + sum(wgamma * outer(rep(1, N), h_t == 1, "*"))
      # False employment: h=0, y=1
      fe <- outer(y_t == 1, h_t == 0, "&") * 1L
      M1 <- M1 + sum(wgamma * fe)
      E1 <- E1 + sum(wgamma * outer(rep(1, N), h_t == 0, "*"))
    }
    attr(M0, "exposure") <- E0
    attr(M1, "exposure") <- E1
  }

  list(
    gamma  = gamma,
    loglik = loglik,
    sufficient_stats = list(
      S  = S,   # 2x2 matrix: S[j+1, k+1] = S_{jk} for alpha update
      D  = D,   # 2x2 matrix: D[j+1, k+1] = D_{jk}
      T1 = T1,  # 2x2 matrix: T1[j+1, k+1] = T_{jk->1}
      M  = M,
      M0 = M0,
      M1 = M1,
      W  = W
    )
  )
}
