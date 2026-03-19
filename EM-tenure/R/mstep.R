# ==============================================================================
# EM-tenure: M-step — parameter updates
# ==============================================================================
# sigma2_g, sigma2_d, alpha, pi are closed-form.
# theta1, theta0 are found via Brent's method: the joint FOC combines the
# Markov prior gradient with the EMG emission gradient (TeX Eqs. foc_theta1
# and foc_theta0).  sigma2 is computed first (closed-form) so the Brent step
# uses the current-iteration sigma2 when evaluating log_emg_grad_lambda.
#
# Newton FOC for theta1 (TeX Eq. foc_theta1):
#   T11/theta1 - (D1-T11)/(1-theta1) + E_g * d_lambda_g/d_theta1 = 0
# where E_g = sum_i w_i * log_emg_grad_lambda(emg_g_x, lambda_g(theta1), sigma2_g)
# and d_lambda_g/d_theta1 = 1 / (3*(1-theta1))  [CTMC link, TeX Eq. 3].
# theta0 satisfies the symmetric equation with (T01, D0, emg_d_x/w, sigma2_d).
#
# The pooled variance estimator is from TeX Eq (21):
#   sigma_g^2 = (S_g + 2*S_g_start) / (2*(N_g + N_g_start))
#
# TeX reference: Section 2.7 (M-step), Eqs (19)-(22)
# ==============================================================================

# --- Internal: FOC evaluator for joint theta update --------------------------
#
# f(theta) = T_stay/theta - (D_from - T_stay)/(1-theta)
#              + E_emg(theta) * d_lambda/d_theta
# where d_lambda/d_theta = 1/(3*(1-theta)) [TeX Eq. 3: lambda = -log(1-theta)/3].
# E_emg = sum(emg_w * log_emg_grad_lambda(emg_x, lambda(theta), sigma2)).
#
# @keywords internal
.theta_foc <- function(theta, T_stay, D_from, emg_x, emg_w, sigma2) {
  lambda    <- ctmc_lambda_from_theta(theta)
  # d_lambda/d_theta = 1/(3*(1-theta)) for lambda = -log(1-theta)/3
  dl_dtheta <- 1 / (3 * (1 - theta))

  E_emg <- if (length(emg_x) > 0) {
    sum(emg_w * log_emg_grad_lambda(emg_x, lambda, sigma2), na.rm = TRUE)
  } else {
    0
  }

  T_stay / theta - (D_from - T_stay) / (1 - theta) + E_emg * dl_dtheta
}


# --- Internal: Brent's method solver for joint theta update ------------------
#
# Finds the root of .theta_foc() in (eps, theta_cap) using Brent's method.
# Brent's method combines bisection, secant, and inverse quadratic
# interpolation; it is bracketed and guaranteed to converge.
#
# The FOC is continuous and monotone decreasing in theta (the Markov terms
# dominate for theta away from boundaries), so a unique root exists in the
# feasible interval.
#
# @param T_stay   Scalar: weighted count of same-state transitions (T11 or T01).
# @param D_from   Scalar: weighted count of periods starting in state (D1 or D0).
# @param emg_x    Numeric vector: duration observations contributing EMG emission.
# @param emg_w    Numeric vector: corresponding weights (same length as emg_x).
# @param sigma2   Scalar: measurement variance for this state.
# @param theta_seed Scalar: starting value (used to guide bracket search).
# @param theta_cap  Scalar: upper bound for theta.
# @param eps      Scalar: lower bound for theta (default 1e-6).
# @param tol      Scalar: convergence tolerance (default 1e-10).
# @param max_iter Integer: maximum Brent iterations (default 100).
# @return Scalar theta in (eps, theta_cap).
# @keywords internal
.m_step_theta_newton <- function(T_stay, D_from, emg_x, emg_w, sigma2,
                                 theta_seed, theta_cap,
                                 eps = 1e-6, tol = 1e-10, max_iter = 100L) {
  lo <- eps
  hi <- theta_cap - eps

  f_lo <- .theta_foc(lo, T_stay, D_from, emg_x, emg_w, sigma2)
  f_hi <- .theta_foc(hi, T_stay, D_from, emg_x, emg_w, sigma2)

  # If FOC doesn't change sign, fall back to the prior-only MLE (seed)
  if (!is.finite(f_lo) || !is.finite(f_hi) || (f_lo * f_hi > 0)) {
    return(max(eps, min(theta_seed, theta_cap - eps)))
  }

  # Brent's method
  a <- lo; b <- hi
  fa <- f_lo; fb <- f_hi

  # Ensure |f(b)| <= |f(a)|
  if (abs(fa) < abs(fb)) {
    a <- hi; b <- lo
    fa <- f_hi; fb <- f_lo
  }

  c_pt  <- a; fc <- fa
  mflag <- TRUE
  s     <- b; fs <- fb
  d_pt  <- 0

  for (k in seq_len(max_iter)) {
    if (abs(fb) < tol || abs(b - a) < tol) break

    if (fa != fc && fb != fc) {
      # Inverse quadratic interpolation
      s <- a * fb * fc / ((fa - fb) * (fa - fc)) +
           b * fa * fc / ((fb - fa) * (fb - fc)) +
           c_pt * fa * fb / ((fc - fa) * (fc - fb))
    } else {
      # Secant
      s <- b - fb * (b - a) / (fb - fa)
    }

    # Conditions to fall back to bisection
    cond1 <- !(s > (3 * a + b) / 4 && s < b) &&
              !(s > b && s < (3 * a + b) / 4)
    cond2 <- mflag  && abs(s - b) >= abs(b - c_pt) / 2
    cond3 <- !mflag && abs(s - b) >= abs(c_pt - d_pt) / 2
    cond4 <- mflag  && abs(b - c_pt) < tol
    cond5 <- !mflag && abs(c_pt - d_pt) < tol

    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s <- (a + b) / 2
      mflag <- TRUE
    } else {
      mflag <- FALSE
    }

    s  <- max(eps, min(s, theta_cap - eps))
    fs <- .theta_foc(s, T_stay, D_from, emg_x, emg_w, sigma2)
    if (!is.finite(fs)) { s <- (a + b) / 2; fs <- .theta_foc(s, T_stay, D_from, emg_x, emg_w, sigma2) }

    d_pt <- c_pt
    c_pt <- b; fc <- fb

    if (fa * fs < 0) {
      b <- s; fb <- fs
    } else {
      a <- s; fa <- fs
    }

    if (abs(fa) < abs(fb)) {
      tmp <- a; a <- b; b <- tmp
      tmp <- fa; fa <- fb; fb <- tmp
    }
  }

  max(eps, min(b, theta_cap - eps))
}


#' M-step: update parameters from sufficient statistics
#'
#' @param suff Named list of sufficient statistics from \code{e_step()}.
#' @param total_weight Sum of survey weights (sum(w_i)).
#' @param misclassification Logical; if FALSE, fix pi=0 (no misclassification).
#' @param stationary Logical; if TRUE, impose alpha = theta0/(theta0 + 1 - theta1).
#' @param sigma_floor Minimum value for sigma2_g and sigma2_d (default 1e-8).
#' @param theta_cap Maximum value for theta1 and theta0 (default 0.999).
#' @param pi_cap Maximum value for pi (default 0.49).
#' @return Named list of updated parameters: alpha, theta0, theta1, pi,
#'   sigma2_g, sigma2_d, lambda_g, lambda_d.
#' @export
m_step <- function(suff, total_weight,
                   misclassification = TRUE,
                   stationary = FALSE,
                   sigma_floor = 1e-8,
                   theta_cap = 0.999,
                   pi_cap = 0.49) {
  # --- Misclassification (Eq 20) ---
  if (misclassification) {
    pi_hat <- suff$M / (3 * total_weight)
    pi_hat <- .bound01(pi_hat, eps = 0)
    pi_hat <- min(pi_hat, pi_cap)
  } else {
    pi_hat <- 0
  }

  # --- Measurement variances (Eq 21): pool continuations + starts ---
  # Computed BEFORE the theta Newton step so the M-step uses the current-
  # iteration sigma2 when evaluating log_emg_grad_lambda (TeX Section 2.7).
  # sigma_g^2 = (S_g + 2*S_g_start) / (2*(N_g + N_g_start))
  denom_g <- 2 * (suff$Ng + suff$Ng_start)
  sigma2_g <- if (denom_g > 0) {
    (suff$Sg + 2 * suff$Sg_start) / denom_g
  } else {
    sigma_floor
  }
  sigma2_g <- max(sigma2_g, sigma_floor)

  denom_d <- 2 * (suff$Nd + suff$Nd_start)
  sigma2_d <- if (denom_d > 0) {
    (suff$Sd + 2 * suff$Sd_start) / denom_d
  } else {
    sigma_floor
  }
  sigma2_d <- max(sigma2_d, sigma_floor)

  # --- Markov block (Eqs 19, foc_theta1, foc_theta0) ---
  # Seed from the prior-only closed-form MLE; Newton joint FOC refines it.
  # The upper cap theta_cap (default 0.999) keeps CTMC rates finite.
  theta1_seed <- .bound01(suff$T11 / suff$D1, eps = 1e-6)
  theta1 <- .m_step_theta_newton(
    T_stay = suff$T11, D_from = suff$D1,
    emg_x = suff$emg_g_x, emg_w = suff$emg_g_w,
    sigma2 = sigma2_g,
    theta_seed = theta1_seed, theta_cap = theta_cap
  )

  theta0_seed <- .bound01(suff$T01 / suff$D0, eps = 1e-6)
  theta0 <- .m_step_theta_newton(
    T_stay = suff$T01, D_from = suff$D0,
    emg_x = suff$emg_d_x, emg_w = suff$emg_d_w,
    sigma2 = sigma2_d,
    theta_seed = theta0_seed, theta_cap = theta_cap
  )

  if (stationary) {
    alpha <- theta0 / (theta0 + 1 - theta1)
  } else {
    alpha <- .bound01(suff$C1 / (suff$C1 + suff$C0), eps = 1e-6)
  }

  # --- Exponential rates (Eq 22): CTMC link ---
  lambda_g <- ctmc_lambda_from_theta(theta1)
  lambda_d <- ctmc_lambda_from_theta(theta0)

  list(
    alpha    = alpha,
    theta0   = theta0,
    theta1   = theta1,
    pi       = pi_hat,
    sigma2_g = sigma2_g,
    sigma2_d = sigma2_d,
    lambda_g = lambda_g,
    lambda_d = lambda_d
  )
}
