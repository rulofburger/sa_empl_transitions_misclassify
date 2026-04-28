# ==============================================================================
# EM-tenure: M-step — parameter updates
# ==============================================================================
# sigma2_g, alpha, pi are closed-form in both free and linked specifications.
#
# FREE specification (linked = FALSE, default):
#   theta1 = T11/D1 (closed form from Markov prior counts only)
#   theta0 = T01/D0 (closed form)
#   lambda_g via Brent on EMG score equation (emissions only)
#   lambda_d via Brent on interval-censored Exp score (emissions only)
#   The four blocks (theta, sigma, lambda_g, lambda_d) are independent.
#
# LINKED specification (linked = TRUE, legacy/testable restriction):
#   theta1 found via Brent on joint FOC combining Markov prior + EMG gradient,
#     because lambda_g = -log(1-theta1)/3 couples the two blocks.
#   theta0 found via Brent on joint FOC (discrete interval-censored variant).
#   lambda_g, lambda_d set deterministically from theta via CTMC link.
#
# TeX reference: Section 2.7 (M-step), Eqs mstep_alpha through foc_lambda_d
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


# --- Internal: discrete FOC evaluator for theta0 update (discrete mode) -----
#
# f(theta0) = T01/theta0 - (D0-T01)/(1-theta0)
#               + (d_lambda_d/d_theta0) * [E_marginal + E_trans]
# where:
#   E_marginal = sum(w_k * interval_grad_lambda_d(cat_k, lambda_d))
#   E_trans    = sum(w_k * transition_grad_lambda_d(curr_k, prev_k, lambda_d))
#   d_lambda_d/d_theta0 = 1/(3*(1-theta0))
#
# @param theta0    Scalar candidate value.
# @param T01       Scalar: weighted count of nonemployment-to-employment trans.
# @param D0        Scalar: weighted count of periods starting in state 0.
# @param cat_marg  Integer vector: category codes for marginal emissions.
# @param w_marg    Numeric vector: weights for marginal emissions.
# @param cat_curr  Integer vector: current-wave categories for Case 3 transitions.
# @param cat_prev  Integer vector: previous-wave categories for Case 3 transitions.
# @param w_trans   Numeric vector: weights for transition emissions.
# @keywords internal
.theta0_foc_discrete <- function(theta0, T01, D0,
                                  cat_marg, w_marg,
                                  cat_curr, cat_prev, w_trans) {
  lambda_d  <- ctmc_lambda_from_theta(theta0)
  dl_dtheta <- 1 / (3 * (1 - theta0))

  E_marg <- if (length(cat_marg) > 0) {
    sum(w_marg * interval_grad_lambda_d(cat_marg, lambda_d), na.rm = TRUE)
  } else {
    0
  }
  E_trans <- if (length(cat_curr) > 0) {
    sum(w_trans * transition_grad_lambda_d(cat_curr, cat_prev, lambda_d), na.rm = TRUE)
  } else {
    0
  }

  T01 / theta0 - (D0 - T01) / (1 - theta0) + (E_marg + E_trans) * dl_dtheta
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


# --- Internal: Brent's method solver for discrete theta0 update --------------
#
# Finds the root of .theta0_foc_discrete() in (eps, theta_cap).
# Signature mirrors .m_step_theta_newton for symmetry.
#
# @param T01,D0      Scalars: sufficient stats from E-step.
# @param cat_marg,w_marg   Marginal emission categories and weights.
# @param cat_curr,cat_prev,w_trans   Transition emission data.
# @param theta_seed  Starting value for fallback if bracket fails.
# @param theta_cap   Upper bound.
# @param eps,tol,max_iter  Numerics.
# @return Scalar theta0 in (eps, theta_cap).
# @keywords internal
.m_step_theta0_brent_discrete <- function(T01, D0,
                                           cat_marg, w_marg,
                                           cat_curr, cat_prev, w_trans,
                                           theta_seed, theta_cap,
                                           eps = 1e-6, tol = 1e-10,
                                           max_iter = 100L) {
  foc <- function(th) {
    .theta0_foc_discrete(th, T01, D0, cat_marg, w_marg, cat_curr, cat_prev, w_trans)
  }

  lo <- eps;  hi <- theta_cap - eps
  f_lo <- foc(lo);  f_hi <- foc(hi)

  if (!is.finite(f_lo) || !is.finite(f_hi) || (f_lo * f_hi > 0)) {
    return(max(eps, min(theta_seed, theta_cap - eps)))
  }

  a <- lo; b <- hi
  fa <- f_lo; fb <- f_hi

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
      s <- a * fb * fc / ((fa - fb) * (fa - fc)) +
           b * fa * fc / ((fb - fa) * (fb - fc)) +
           c_pt * fa * fb / ((fc - fa) * (fc - fb))
    } else {
      s <- b - fb * (b - a) / (fb - fa)
    }

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
    fs <- foc(s)
    if (!is.finite(fs)) { s <- (a + b) / 2; fs <- foc(s) }

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


# --- Internal: Brent solver for lambda_g (free specification) ----------------
#
# Finds the root of the EMG score equation (TeX Eq foc_lambda_g):
#   sum_i w_i * d/d(lambda_g) log f_EMG(x_i; lambda_g, sigma2_g) = 0
#
# This is a one-dimensional concave problem in lambda_g (given sigma2_g),
# so the root is unique when it exists.
#
# @param emg_x    Numeric vector: EMG-observed durations (in years).
# @param emg_w    Numeric vector: corresponding weights.
# @param sigma2   Scalar: measurement variance sigma2_g (already estimated).
# @param lambda_lo  Scalar: lower bound for lambda_g (default 1e-6).
# @param lambda_hi  Scalar: upper bound for lambda_g (default 50).
# @param tol      Scalar: convergence tolerance (default 1e-10).
# @param max_iter Integer: maximum Brent iterations (default 100).
# @return Scalar lambda_g in (lambda_lo, lambda_hi).
# @keywords internal
.m_step_lambda_g_brent <- function(emg_x, emg_w, sigma2,
                                   lambda_lo = 1e-6, lambda_hi = 50,
                                   tol = 1e-10, max_iter = 100L) {
  if (length(emg_x) == 0) return(0.01)  # fallback: no EMG data

  foc <- function(lam) {
    sum(emg_w * log_emg_grad_lambda(emg_x, lam, sigma2), na.rm = TRUE)
  }

  lo <- lambda_lo; hi <- lambda_hi
  f_lo <- foc(lo); f_hi <- foc(hi)

  # If bracket fails, fall back to moment estimator: lambda ~ 1/mean(x)
  if (!is.finite(f_lo) || !is.finite(f_hi) || (f_lo * f_hi > 0)) {
    wm <- sum(emg_w * emg_x, na.rm = TRUE) / sum(emg_w, na.rm = TRUE)
    return(max(lambda_lo, min(1 / max(wm, 0.01), lambda_hi)))
  }

  # Standard Brent's method (same implementation as .m_step_theta_newton)
  a <- lo; b <- hi; fa <- f_lo; fb <- f_hi
  if (abs(fa) < abs(fb)) { a <- hi; b <- lo; fa <- f_hi; fb <- f_lo }
  c_pt <- a; fc <- fa; mflag <- TRUE; d_pt <- 0

  for (k in seq_len(max_iter)) {
    if (abs(fb) < tol || abs(b - a) < tol) break
    if (fa != fc && fb != fc) {
      s <- a * fb * fc / ((fa - fb) * (fa - fc)) +
           b * fa * fc / ((fb - fa) * (fb - fc)) +
           c_pt * fa * fb / ((fc - fa) * (fc - fb))
    } else {
      s <- b - fb * (b - a) / (fb - fa)
    }
    cond1 <- !(s > (3 * a + b) / 4 && s < b) &&
              !(s > b && s < (3 * a + b) / 4)
    cond2 <- mflag  && abs(s - b) >= abs(b - c_pt) / 2
    cond3 <- !mflag && abs(s - b) >= abs(c_pt - d_pt) / 2
    cond4 <- mflag  && abs(b - c_pt) < tol
    cond5 <- !mflag && abs(c_pt - d_pt) < tol
    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s <- (a + b) / 2; mflag <- TRUE
    } else {
      mflag <- FALSE
    }
    s <- max(lambda_lo, min(s, lambda_hi))
    fs <- foc(s)
    if (!is.finite(fs)) { s <- (a + b) / 2; fs <- foc(s) }
    d_pt <- c_pt; c_pt <- b; fc <- fb
    if (fa * fs < 0) { b <- s; fb <- fs } else { a <- s; fa <- fs }
    if (abs(fa) < abs(fb)) {
      tmp <- a; a <- b; b <- tmp; tmp <- fa; fa <- fb; fb <- tmp
    }
  }
  max(lambda_lo, min(b, lambda_hi))
}


# --- Internal: Brent solver for lambda_d (free, discrete timegap) -----------
#
# Finds the root of the interval-censored exponential score (TeX Eq foc_lambda_d):
#   sum_k w_k * d/d(lambda_d) log P(D in [a_k, b_k)) = 0
# combining marginal emissions and transition emissions.
#
# @param cat_marg  Integer vector: marginal category codes (1:7).
# @param w_marg    Numeric vector: weights for marginal observations.
# @param cat_curr  Integer vector: current-wave categories for transitions.
# @param cat_prev  Integer vector: previous-wave categories for transitions.
# @param w_trans   Numeric vector: weights for transition observations.
# @param lambda_lo Scalar: lower bound (default 1e-6).
# @param lambda_hi Scalar: upper bound (default 50).
# @param tol       Scalar: convergence tolerance (default 1e-10).
# @param max_iter  Integer: maximum iterations (default 100).
# @return Scalar lambda_d in (lambda_lo, lambda_hi).
# @keywords internal
.m_step_lambda_d_brent <- function(cat_marg, w_marg,
                                   cat_curr, cat_prev, w_trans,
                                   lambda_lo = 1e-6, lambda_hi = 50,
                                   tol = 1e-10, max_iter = 100L) {
  # No data at all → fallback
  if (length(cat_marg) == 0 && length(cat_curr) == 0) return(0.01)

  foc <- function(lam) {
    E_marg <- if (length(cat_marg) > 0) {
      sum(w_marg * interval_grad_lambda_d(cat_marg, lam), na.rm = TRUE)
    } else { 0 }
    E_trans <- if (length(cat_curr) > 0) {
      sum(w_trans * transition_grad_lambda_d(cat_curr, cat_prev, lam), na.rm = TRUE)
    } else { 0 }
    E_marg + E_trans
  }

  lo <- lambda_lo; hi <- lambda_hi
  f_lo <- foc(lo); f_hi <- foc(hi)

  if (!is.finite(f_lo) || !is.finite(f_hi) || (f_lo * f_hi > 0)) {
    return(0.01)
  }

  a <- lo; b <- hi; fa <- f_lo; fb <- f_hi
  if (abs(fa) < abs(fb)) { a <- hi; b <- lo; fa <- f_hi; fb <- f_lo }
  c_pt <- a; fc <- fa; mflag <- TRUE; d_pt <- 0

  for (k in seq_len(max_iter)) {
    if (abs(fb) < tol || abs(b - a) < tol) break
    if (fa != fc && fb != fc) {
      s <- a * fb * fc / ((fa - fb) * (fa - fc)) +
           b * fa * fc / ((fb - fa) * (fb - fc)) +
           c_pt * fa * fb / ((fc - fa) * (fc - fb))
    } else {
      s <- b - fb * (b - a) / (fb - fa)
    }
    cond1 <- !(s > (3 * a + b) / 4 && s < b) &&
              !(s > b && s < (3 * a + b) / 4)
    cond2 <- mflag  && abs(s - b) >= abs(b - c_pt) / 2
    cond3 <- !mflag && abs(s - b) >= abs(c_pt - d_pt) / 2
    cond4 <- mflag  && abs(b - c_pt) < tol
    cond5 <- !mflag && abs(c_pt - d_pt) < tol
    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s <- (a + b) / 2; mflag <- TRUE
    } else {
      mflag <- FALSE
    }
    s <- max(lambda_lo, min(s, lambda_hi))
    fs <- foc(s)
    if (!is.finite(fs)) { s <- (a + b) / 2; fs <- foc(s) }
    d_pt <- c_pt; c_pt <- b; fc <- fb
    if (fa * fs < 0) { b <- s; fb <- fs } else { a <- s; fa <- fs }
    if (abs(fa) < abs(fb)) {
      tmp <- a; a <- b; b <- tmp; tmp <- fa; fa <- fb; fb <- tmp
    }
  }
  max(lambda_lo, min(b, lambda_hi))
}


#' M-step: update parameters from sufficient statistics
#'
#' @param suff Named list of sufficient statistics from \code{e_step()}.
#' @param total_weight Sum of survey weights (sum(w_i)).
#' @param stationary Logical; if TRUE, impose alpha = theta0/(theta0 + 1 - theta1).
#' @param linked Logical (default FALSE). If TRUE, use CTMC link (lambda = f(theta)),
#'   solving a joint FOC for theta. If FALSE (free specification), estimate theta
#'   from Markov counts only and lambda from emissions only.
#' @param discrete_timegap Logical (default TRUE). If TRUE, use interval-censored
#'   discrete FOC for theta0/lambda_d and skip sigma2_d estimation (returns NA).
#'   If FALSE, use legacy continuous EMG FOC and estimate sigma2_d.
#' @param sigma_floor Minimum value for sigma2_g and sigma2_d (default 1e-8).
#' @param theta_cap Maximum value for theta1 and theta0 (default 0.999).
#' @param pi_cap Maximum value for pi (default 0.49).
#' @return Named list of updated parameters: alpha, theta0, theta1, pi,
#'   sigma2_g, lambda_g, lambda_d. sigma2_d is included only when
#'   discrete_timegap = FALSE.
#' @export
m_step <- function(suff, total_weight,
                   stationary = FALSE,
                   linked = FALSE,
                   discrete_timegap = TRUE,
                   sigma_floor = 1e-8,
                   theta_cap = 0.999,
                   pi_cap = 0.49) {
  # --- Misclassification (Eq 20) ---
  pi_hat <- suff$M / (3 * total_weight)
  pi_hat <- .bound01(pi_hat, eps = 0)
  pi_hat <- min(pi_hat, pi_cap)

  # --- Measurement variance for tenure (Eq 21): pool continuations + starts ---
  # sigma_g^2 = (S_g + 2*S_g_start) / (2*(N_g + N_g_start))
  denom_g <- 2 * (suff$Ng + suff$Ng_start)
  sigma2_g <- if (denom_g > 0) {
    (suff$Sg + 2 * suff$Sg_start) / denom_g
  } else {
    sigma_floor
  }
  sigma2_g <- max(sigma2_g, sigma_floor)

  # sigma2_d only estimated in continuous mode
  if (!discrete_timegap) {
    denom_d <- 2 * (suff$Nd + suff$Nd_start)
    sigma2_d <- if (denom_d > 0) {
      (suff$Sd + 2 * suff$Sd_start) / denom_d
    } else {
      sigma_floor
    }
    sigma2_d <- max(sigma2_d, sigma_floor)
  }

  # --- theta1 ---
  theta1_seed <- .bound01(suff$T11 / suff$D1, eps = 1e-6)

  if (linked) {
    # LINKED: joint FOC combining Markov prior + EMG emission (TeX Eq foc_theta1_linked)
    theta1 <- .m_step_theta_newton(
      T_stay = suff$T11, D_from = suff$D1,
      emg_x = suff$emg_g_x, emg_w = suff$emg_g_w,
      sigma2 = sigma2_g,
      theta_seed = theta1_seed, theta_cap = theta_cap
    )
  } else {
    # FREE: closed form from Markov prior counts only (TeX Eq mstep_theta_free)
    theta1 <- max(1e-6, min(theta1_seed, theta_cap))
  }

  # --- theta0 ---
  theta0_seed <- .bound01(suff$T01 / suff$D0, eps = 1e-6)

  if (linked) {
    # LINKED: joint FOC (discrete or continuous variant)
    if (discrete_timegap) {
      theta0 <- .m_step_theta0_brent_discrete(
        T01 = suff$T01, D0 = suff$D0,
        cat_marg = suff$cat_d_marginal_c,
        w_marg   = suff$cat_d_marginal_w,
        cat_curr = suff$cat_d_trans_curr,
        cat_prev = suff$cat_d_trans_prev,
        w_trans  = suff$cat_d_trans_w,
        theta_seed = theta0_seed, theta_cap = theta_cap
      )
    } else {
      theta0 <- .m_step_theta_newton(
        T_stay = suff$T01, D_from = suff$D0,
        emg_x = suff$emg_d_x, emg_w = suff$emg_d_w,
        sigma2 = sigma2_d,
        theta_seed = theta0_seed, theta_cap = theta_cap
      )
    }
  } else {
    # FREE: closed form from Markov prior counts only
    theta0 <- max(1e-6, min(theta0_seed, theta_cap))
  }

  if (stationary) {
    alpha <- theta0 / (theta0 + 1 - theta1)
  } else {
    alpha <- .bound01(suff$C1 / (suff$C1 + suff$C0), eps = 1e-6)
  }

  # --- Exponential rates ---
  if (linked) {
    # CTMC link: deterministic from theta
    lambda_g <- ctmc_lambda_from_theta(theta1)
    lambda_d <- ctmc_lambda_from_theta(theta0)
  } else {
    # FREE: Brent on emission-only score equations (TeX Eqs foc_lambda_g, foc_lambda_d)
    lambda_g <- .m_step_lambda_g_brent(
      emg_x = suff$emg_g_x, emg_w = suff$emg_g_w,
      sigma2 = sigma2_g
    )
    if (discrete_timegap) {
      lambda_d <- .m_step_lambda_d_brent(
        cat_marg = suff$cat_d_marginal_c,
        w_marg   = suff$cat_d_marginal_w,
        cat_curr = suff$cat_d_trans_curr,
        cat_prev = suff$cat_d_trans_prev,
        w_trans  = suff$cat_d_trans_w
      )
    } else {
      # Continuous free: Brent on EMG score for d side
      lambda_d <- .m_step_lambda_g_brent(
        emg_x = suff$emg_d_x, emg_w = suff$emg_d_w,
        sigma2 = sigma2_d
      )
    }
  }

  out <- list(
    alpha    = alpha,
    theta0   = theta0,
    theta1   = theta1,
    pi       = pi_hat,
    sigma2_g = sigma2_g,
    lambda_g = lambda_g,
    lambda_d = lambda_d
  )
  if (!discrete_timegap) {
    out$sigma2_d <- sigma2_d
  }
  out
}
