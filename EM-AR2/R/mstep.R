# ==============================================================================
# EM-AR2: M-step — closed-form parameter updates
# Created: 2026-05-05
# ==============================================================================
# All updates are closed-form because the AR(2) model has no duration
# emissions — only binary employment sequences. The key updates are:
#
# Transition probability updates (TeX: M-step section):
#   p_hat_{jk->1} = T_{jk->1} / D_{jk}    for (j,k) in {0,1}^2
#
# Recover AR(2) parameters from the four p_{jk->1} estimates:
#   theta0   = p_{00->1}
#   theta01  = p_{10->1} - p_{00->1}
#   theta1   = 1 - p_{11->1}
#   theta10  = p_{11->1} - p_{01->1}
#
# Misclassification update:
#   pi_hat = M / (4 * W)          [symmetric: M misclassified obs over 4*N total]
#   pi0_hat = M0 / E0              [asymmetric: false non-employment / exposure]
#   pi1_hat = M1 / E1              [asymmetric: false employment / exposure]
#
# IDENTITIES:
#   p_{00->1} = theta0                      (job-finding from (0,0))
#   p_{10->1} = theta0 + theta01            (job-finding from (1,0))
#   p_{01->1} = 1 - theta1 - theta10        (staying employed from (0,1))
#   p_{11->1} = 1 - theta1                  (staying employed from (1,1))
# ==============================================================================

#' M-step: update AR(2) parameters from sufficient statistics
#'
#' Closed-form parameter updates given weighted sufficient statistics from the
#' E-step. Implements the standard Baum-Welch approach:
#'   1. The free initial pair distribution alpha(h1,h2) is updated directly
#'      from the S_{jk} counts: alpha_hat(j,k) = S_{jk} / sum(S).
#'   2. The transition probabilities p_{jk->1} are updated from T/D counts.
#'   3. Structural parameters theta are recovered from the p's.
#'
#' Using free alpha (rather than stationarity-constrained alpha) guarantees
#' non-decreasing log-likelihood. After convergence, the stationary
#' distribution can be computed from theta for inference purposes.
#'
#' @param ss Sufficient statistics list from e_step_ar2()$sufficient_stats.
#'   Must contain S (2x2), D (2x2), T1 (2x2), M, M0, M1, W.
#' @param estimate_pi Logical; if TRUE estimate pi from data. If FALSE,
#'   use fixed_pi (default TRUE).
#' @param fixed_pi Scalar; value to use for pi when estimate_pi=FALSE
#'   (default 0 = no misclassification).
#' @param asymmetric Logical; if TRUE estimate separate pi0 and pi1 (default FALSE).
#'   Overrides estimate_pi when TRUE.
#' @param pi_cap Maximum allowed misclassification probability (default 0.49).
#' @param theta_eps Minimum bound for transition probabilities (default 1e-6).
#' @return Named list: alpha, theta0, theta01, theta1, theta10, pi [, pi0, pi1].
#'   Also includes the intermediate p_hat values for diagnostics.
#' @export
m_step_ar2 <- function(ss,
                        estimate_pi  = TRUE,
                        fixed_pi     = 0,
                        asymmetric   = FALSE,
                        pi_cap       = 0.49,
                        theta_eps    = 1e-6) {
  S  <- ss$S    # 2x2 matrix: S_{jk} for initial pair distribution
  D  <- ss$D    # 2x2 matrix, rows = j (prev2), cols = k (prev1)
  T1 <- ss$T1   # 2x2 matrix: T_{jk->1}
  M  <- ss$M
  M0 <- ss$M0
  M1 <- ss$M1
  W  <- ss$W

  # --- Update free initial pair distribution alpha(h1,h2) ---
  # alpha_hat(j,k) = S_{jk} / sum_{j',k'} S_{j'k'}
  S_total <- sum(S)
  alpha_new <- if (S_total > 1e-12) {
    a <- as.vector(S / S_total)
    names(a) <- c("00", "10", "01", "11")
    a
  } else {
    c(`00`=0.25, `10`=0.25, `01`=0.25, `11`=0.25)
  }
  alpha_new <- pmax(alpha_new, 1e-12)
  alpha_new <- alpha_new / sum(alpha_new)  # re-normalise for floating point

  # --- Update transition probabilities p_{jk->1} = T_{jk->1} / D_{jk} ---
  # Guard: if D_{jk} is near zero, use a neutral value of 0.5
  safe_div <- function(num, den) {
    ifelse(den > 1e-12, num / den, 0.5)
  }

  p_00_1 <- safe_div(T1[1, 1], D[1, 1])  # P(h_t=1 | h_{t-2}=0, h_{t-1}=0)
  p_10_1 <- safe_div(T1[2, 1], D[2, 1])  # P(h_t=1 | h_{t-2}=1, h_{t-1}=0)
  p_01_1 <- safe_div(T1[1, 2], D[1, 2])  # P(h_t=1 | h_{t-2}=0, h_{t-1}=1)
  p_11_1 <- safe_div(T1[2, 2], D[2, 2])  # P(h_t=1 | h_{t-2}=1, h_{t-1}=1)

  # Clamp to (0,1)
  p_00_1 <- .bound01(p_00_1, eps = theta_eps)
  p_10_1 <- .bound01(p_10_1, eps = theta_eps)
  p_01_1 <- .bound01(p_01_1, eps = theta_eps)
  p_11_1 <- .bound01(p_11_1, eps = theta_eps)

  # --- Recover AR(2) parameters from p_{jk->1} ---
  # theta0  = p_{00->1}
  # theta01 = p_{10->1} - p_{00->1}
  # theta1  = 1 - p_{11->1}
  # theta10 = p_{11->1} - p_{01->1}
  theta0  <- p_00_1
  theta01 <- p_10_1 - p_00_1
  theta1  <- 1 - p_11_1
  theta10 <- p_11_1 - p_01_1

  # Validate recovered parameters
  # theta01 and theta10 can in principle be negative if the transition matrix
  # estimate is inconsistent; clamp with a small negative floor and warn.
  if (theta01 < -0.01) {
    warning("M-step: theta01 = ", round(theta01, 4), " < 0. ",
            "Clamping to theta_eps. Consider checking initialisation.")
  }
  if (theta10 < -0.01) {
    warning("M-step: theta10 = ", round(theta10, 4), " < 0. ",
            "Clamping to theta_eps. Consider checking initialisation.")
  }

  theta0  <- .bound01(theta0,  eps = theta_eps)
  theta01 <- max(theta01,  theta_eps)
  theta1  <- .bound01(theta1,  eps = theta_eps)
  theta10 <- max(theta10, theta_eps)

  # --- Misclassification update ---
  if (asymmetric) {
    E0 <- attr(M0, "exposure") %||% (4 * W / 2)  # fallback
    E1 <- attr(M1, "exposure") %||% (4 * W / 2)

    pi0 <- .bound01(
      if (E0 > 1e-12) M0 / E0 else 0.05,
      eps = 1e-6
    )
    pi1 <- .bound01(
      if (E1 > 1e-12) M1 / E1 else 0.05,
      eps = 1e-6
    )

    pi0 <- min(pi0, pi_cap)
    pi1 <- min(pi1, pi_cap)

    return(list(
      alpha   = alpha_new,
      theta0  = theta0,
      theta01 = theta01,
      theta1  = theta1,
      theta10 = theta10,
      pi0 = pi0,
      pi1 = pi1,
      # Diagnostics
      .p_00_1 = p_00_1,
      .p_10_1 = p_10_1,
      .p_01_1 = p_01_1,
      .p_11_1 = p_11_1
    ))
  }

  if (!estimate_pi) {
    pi <- fixed_pi
  } else {
    # Symmetric: pi_hat = M / (4*W)
    pi <- if (W > 1e-12) M / (4 * W) else 0.05
    pi <- .bound01(pi, eps = 1e-6)
    pi <- min(pi, pi_cap)
  }

  list(
    alpha   = alpha_new,
    theta0  = theta0,
    theta01 = theta01,
    theta1  = theta1,
    theta10 = theta10,
    pi      = pi,
    # Diagnostics: intermediate p_hat values
    .p_00_1 = p_00_1,
    .p_10_1 = p_10_1,
    .p_01_1 = p_01_1,
    .p_11_1 = p_11_1
  )
}
