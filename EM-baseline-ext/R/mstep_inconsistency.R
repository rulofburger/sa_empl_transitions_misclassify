# ==============================================================================
# EM-baseline-ext: M-step for Extension IV (inconsistency-augmented model)
# Created: 2026-05-06
#
# Theta updates: closed form (identical to baseline).
# Delta update:  one Fisher-scoring Newton-Raphson step on Q_miscl(delta).
#
# Q_miscl(delta) = sum_{i,t} w_i * [p_it * log(pi_it) + (1-p_it)*log(1-pi_it)]
# where pi_it = 0.5 * sigma(delta_0 + delta_1*Y_age_it + delta_2*Y_edu_it)
# and p_it = P(s_it != h_t | y_i) are held fixed (from E-step).
#
# Gradient: g_j = sum_{i,t} w_i * D_it * (p_it/pi_it - (1-p_it)/(1-pi_it)) * f_it_j
#   where D_it = 0.5 * sigma_it * (1-sigma_it),  f_it = (1, Y_age_it, Y_edu_it)
# Fisher information: I_jk = sum_{i,t} w_i * D_it^2 * (p_it/pi_it^2 + (1-p_it)/(1-pi_it)^2) * f_it_j * f_it_k
# NR step (Fisher scoring): delta_new = delta_old + (I + lambda*I3)^{-1} * g
#
# Armijo backtracking ensures Q_miscl does not decrease after the update.
#
# TeX ref: EM baseline.tex Section 8, Eqs (31)--(35).
# ==============================================================================

.INCNS_RIDGE <- 1e-6   # ridge for (I + lambda * I_3)^{-1} solve

#' Evaluate Q_miscl(delta) given fixed posterior mismatch probs
#'
#' Internal helper used for Armijo backtracking.
#'
#' @param delta Length-3 vector.
#' @param p_mat N x 3 matrix of posterior mismatch probabilities.
#' @param Y_age N x 3 matrix of age inconsistency indicators.
#' @param Y_edu N x 3 matrix of education inconsistency indicators.
#' @param w N-vector of observation weights.
#' @return Scalar Q_miscl value.
.eval_q_miscl <- function(delta, p_mat, Y_age, Y_edu, w) {
  p_mat   <- pmin(pmax(p_mat, 0), 1)  # clamp before log
  eta_mat <- delta[1L] + delta[2L] * Y_age + delta[3L] * Y_edu
  pi_mat  <- 0.5 * plogis(eta_mat)
  pi_mat  <- pmax(pmin(pi_mat, 1 - 1e-8), 1e-8)
  sum(w * rowSums(p_mat * log(pi_mat) + (1 - p_mat) * log(1 - pi_mat)))
}


#' Compute gradient and Fisher information for Q_miscl
#'
#' Internal helper.
#'
#' @param delta Length-3 vector.
#' @param p_mat N x 3 matrix of posterior mismatch probabilities.
#' @param Y_age,Y_edu N x 3 indicator matrices.
#' @param w N-vector of weights.
#' @return List with \code{grad} (length 3) and \code{fisher} (3 x 3 matrix).
.grad_fisher_q_miscl <- function(delta, p_mat, Y_age, Y_edu, w) {
  N       <- nrow(p_mat)
  p_mat   <- pmin(pmax(p_mat, 0), 1)  # clamp before division
  eta_mat <- delta[1L] + delta[2L] * Y_age + delta[3L] * Y_edu  # N x 3
  sig_mat <- plogis(eta_mat)
  pi_mat  <- 0.5 * sig_mat
  pi_mat  <- pmax(pmin(pi_mat, 1 - 1e-8), 1e-8)
  sig_mat <- pmin(pmax(sig_mat, 1e-8), 1 - 1e-8)

  D_mat    <- 0.5 * sig_mat * (1 - sig_mat)  # d pi_it / d eta_it  (N x 3)
  score_it <- D_mat * (p_mat / pi_mat - (1 - p_mat) / (1 - pi_mat))  # N x 3
  info_it  <- D_mat^2 * (p_mat / pi_mat^2 + (1 - p_mat) / (1 - pi_mat)^2)  # N x 3

  # f_it = (1, Y_age_it, Y_edu_it) -- 3-vector per (i,t)
  # Features: [N x 3] x 3 = effectively N*3 observations with 3 features
  # Gradient: g_j = sum_{i,t} w_i * score_it * f_it_j
  # Fisher:   I_jk = sum_{i,t} w_i * info_it * f_it_j * f_it_k
  # Expand weights: W_it = w_i (broadcast to N x 3)
  W3  <- matrix(w, nrow = N, ncol = 3L)  # N x 3
  WS  <- W3 * score_it
  WI  <- W3 * info_it

  # Gradient and Fisher information, treating intercept (j=1) separately
  # to avoid allocating matrix(1, N, 3L):
  grad   <- numeric(3L)
  fisher <- matrix(0, 3L, 3L)

  grad[1L]       <- sum(WS)
  grad[2L]       <- sum(WS * Y_age)
  grad[3L]       <- sum(WS * Y_edu)

  fisher[1L, 1L] <- sum(WI)
  fisher[1L, 2L] <- sum(WI * Y_age);         fisher[2L, 1L] <- fisher[1L, 2L]
  fisher[1L, 3L] <- sum(WI * Y_edu);         fisher[3L, 1L] <- fisher[1L, 3L]
  # Y_age and Y_edu are binary, so Y_age^2 == Y_age (squaring is a no-op)
  fisher[2L, 2L] <- sum(WI * Y_age)
  fisher[2L, 3L] <- sum(WI * Y_age * Y_edu); fisher[3L, 2L] <- fisher[2L, 3L]
  fisher[3L, 3L] <- sum(WI * Y_edu)

  list(grad = grad, fisher = fisher)
}


#' M-step for the inconsistency-augmented model
#'
#' Performs closed-form theta updates and one Newton-Raphson (Fisher scoring)
#' step for the misclassification coefficient vector \eqn{\delta}. Armijo
#' backtracking ensures \eqn{Q_{\text{miscl}}} does not decrease.
#'
#' TeX ref: EM baseline.tex Section 8, Eqs (31)--(35).
#'
#' @param suff Sufficient statistics list from \code{\link{e_step_inconsistency}}.
#'   Must include \code{p_mat} (N \times 3) and \code{weights} (N-vector).
#' @param incons_mat N × 6 matrix from \code{\link{compute_inconsistencies}}.
#' @param params_old Named list with current parameters (used for theta and
#'   for initialising the NR step).
#' @param stationary Logical (default \code{TRUE}).
#' @param theta_cap Upper bound for transition probabilities (default 0.999).
#' @param nr_max_backtrack Maximum Armijo backtracking steps (default 10L).
#' @return Named list with updated \code{theta0}, \code{theta1}, \code{alpha},
#'   \code{delta}.
#' @examples
#' \dontrun{
#'   df      <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                         y3 = rbinom(50,1,.6), weight = rep(1,50),
#'                         age1=25:74, age2=26:75, age3=27:76,
#'                         educ1=rep(3L,50), educ2=rep(3L,50), educ3=rep(3L,50))
#'   inc_mat <- as.matrix(compute_inconsistencies(df)[,
#'     c("Y_age_1","Y_age_2","Y_age_3","Y_edu_1","Y_edu_2","Y_edu_3")])
#'   out_e   <- e_step_inconsistency(df, inc_mat, init_params_inconsistency())
#'   m_step_inconsistency(out_e$suff, inc_mat, init_params_inconsistency())
#' }
#' @export
m_step_inconsistency <- function(suff, incons_mat, params_old,
                                  stationary        = TRUE,
                                  theta_cap         = 0.999,
                                  nr_max_backtrack  = 10L) {
  eps <- 1e-10

  # ---- Closed-form theta updates ------------------------------------------
  theta1 <- min(.bound01(suff$T11 / max(suff$D1, eps)), theta_cap)
  theta0 <- min(.bound01(suff$T01 / max(suff$D0, eps)), theta_cap)
  alpha  <- if (stationary) {
    stationary_alpha(theta0, theta1)
  } else {
    .bound01(suff$C1 / max(suff$C1 + suff$C0, eps))
  }

  # ---- NR step for delta ---------------------------------------------------
  N      <- length(suff$weights)
  w      <- suff$weights
  delta0 <- params_old$delta
  p_mat  <- suff$p_mat
  Y_age  <- incons_mat[, 1:3, drop = FALSE]
  Y_edu  <- incons_mat[, 4:6, drop = FALSE]

  gf     <- .grad_fisher_q_miscl(delta0, p_mat, Y_age, Y_edu, w)
  g      <- gf$grad
  I_mat  <- gf$fisher

  # Add ridge to Fisher information for numerical stability
  I_ridge <- I_mat + .INCNS_RIDGE * diag(3L)

  # Fisher scoring direction
  step <- tryCatch(
    solve(I_ridge, g),
    error = function(e) {
      warning("m_step_inconsistency: Fisher matrix singular, using gradient only")
      g
    }
  )

  # ---- Armijo backtracking -------------------------------------------------
  q0        <- .eval_q_miscl(delta0, p_mat, Y_age, Y_edu, w)
  alpha_arm <- 1.0  # step size
  delta_new <- delta0
  for (.bt in seq_len(nr_max_backtrack)) {
    delta_try <- delta0 + alpha_arm * step
    q_try     <- .eval_q_miscl(delta_try, p_mat, Y_age, Y_edu, w)
    if (q_try >= q0 - 1e-10) {
      delta_new <- delta_try
      break
    }
    alpha_arm <- alpha_arm / 2
  }

  if (identical(delta_new, delta0) && nr_max_backtrack > 0L) {
    # No improvement found — keep old delta (GEM guarantee: Q cannot decrease
    # if we stay put, which is valid since the zero step is always feasible)
    delta_new <- delta0
  }

  list(
    theta0 = theta0,
    theta1 = theta1,
    alpha  = alpha,
    delta  = delta_new
  )
}
