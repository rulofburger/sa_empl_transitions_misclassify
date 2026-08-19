# ==============================================================================
# EM-tenure: M-step for the eps (Spec I) model
# ==============================================================================
# Created: 2026-04-30
#
# Internal helper
.clip_param <- function(x, lo, hi) max(lo, min(x, hi))

# Joint stationary theta block. Once alpha is constrained to the ergodic
# distribution, the initial-state part of Q depends on both transition
# probabilities. They therefore cannot be updated separately from transition
# counts. Under the CTMC link, the duration-emission terms also belong to this
# same two-dimensional block.
.q_theta_stationary_eps <- function(theta0, theta1, suff, linked = FALSE) {
  alpha <- theta0 / (theta0 + 1 - theta1)
  q <- suff$C1 * log(alpha) + suff$C0 * log1p(-alpha) +
    suff$T01 * log(theta0) + (suff$D0 - suff$T01) * log1p(-theta0) +
    suff$T11 * log(theta1) + (suff$D1 - suff$T11) * log1p(-theta1)
  if (!linked) return(q)

  lambda_g <- ctmc_lambda_from_persistence(theta1)
  lambda_d <- ctmc_lambda_from_transition(theta0)
  q <- q + suff$Lg_count * log(lambda_g) - lambda_g * suff$Lg_xsum
  if (length(suff$cat_d_marginal_c)) {
    use <- suff$cat_d_marginal_w > 0
    q <- q + sum(suff$cat_d_marginal_w[use] *
      log_emission_interval_d(suff$cat_d_marginal_c[use], lambda_d))
  }
  if (length(suff$cat_d_trans_curr)) {
    use <- suff$cat_d_trans_w > 0
    q <- q + sum(suff$cat_d_trans_w[use] * log_emission_transition_d(
      suff$cat_d_trans_curr[use], suff$cat_d_trans_prev[use], lambda_d))
  }
  q
}

.m_step_theta_stationary_eps <- function(suff, linked = FALSE,
                                          theta_cap = 0.999,
                                          theta_floor = 1e-6) {
  width <- theta_cap - 2 * theta_floor
  unpack <- function(z) theta_floor + width * plogis(z)
  pack <- function(x) qlogis((x - theta_floor) / width)
  seed0 <- .clip_param(suff$T01 / suff$D0, theta_floor * 2,
                       theta_cap - theta_floor * 2)
  seed1 <- .clip_param(suff$T11 / suff$D1, theta_floor * 2,
                       theta_cap - theta_floor * 2)
  scale <- max(suff$C0 + suff$C1, 1)
  objective <- function(z) {
    th <- unpack(z)
    -.q_theta_stationary_eps(th[1], th[2], suff, linked) / scale
  }
  starts <- list(c(seed0, seed1), c(0.05, 0.95), c(0.10, 0.90))
  fits <- lapply(starts, function(s) tryCatch(
    optim(pack(s), objective, method = "BFGS",
          control = list(maxit = 500L, reltol = 1e-12)),
    error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("stationary epsilon theta M-step failed from all starts")
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  theta <- unpack(best$par)
  list(theta0 = unname(theta[1]), theta1 = unname(theta[2]),
       optimizer_code = best$convergence,
       q_value = .q_theta_stationary_eps(theta[1], theta[2], suff, linked))
}
#
# Closed-form updates for: alpha, theta_0, theta_1 (free), pi, eps, lambda_g.
# Brent solver for: theta_1 / theta_0 (linked via CTMC), lambda_d (mirroring base).
#
# eps M-step (full 2^K per-wave contamination enumeration, per-wave Bernoulli MLE):
#   eps_hat = Eps_num / Eps_den
#           = sum(w * gamma * tau_sum) / sum(w * gamma * K)
# where tau_sum = E[# contaminated waves | g] in [0, K] from emissions_eps.R,
# and K is the number of observed tenures in the spell (K >= 2 only).
# The per-wave Bernoulli MLE is the correct fixed-point for the full 2^K emission
# (emission model and M-step formula share the same Q function).
# See EM tenure epsilon.tex, Section 4.3 and brainstorm 2026-05-01.
#
# lambda_g M-step (Exponential MLE; closed form):
#   lambda_g_hat = Lg_count / Lg_xsum
# where Lg_count and Lg_xsum aggregate per-spell Exp-emission counts and
# x-sums under the clean / contaminated decomposition (see emissions_eps.R).
#
# Companion: documents/EM tenure epsilon.tex (Section 4.3).
# ==============================================================================

#' Joint first-order condition for theta_1 under the CTMC link (eps model)
#'
#' Score: T11/theta - (D1-T11)/(1-theta) + (Lg_count/lambda - Lg_xsum) * dl/dtheta
#' where lambda_g = -log(theta)/Delta and dl/dtheta = -1/(theta * Delta).
#'
#' @param theta1 Persistence parameter in (0, 1).
#' @param T11 Sufficient statistic: E->E transition count.
#' @param D1 Sufficient statistic: employed-wave denominator.
#' @param Lg_count Exponential emission count for lambda_g M-step.
#' @param Lg_xsum Exponential emission x-sum for lambda_g M-step.
#' @return Scalar score (root gives the M-step optimum).
#' @keywords internal
.theta1_foc_eps <- function(theta1, T11, D1, Lg_count, Lg_xsum) {
  lambda    <- ctmc_lambda_from_persistence(theta1)
  dl_dtheta <- -1 / (theta1 * .QUARTER_YEARS)
  exp_score <- if (lambda > 0) Lg_count / lambda - Lg_xsum else 0
  T11 / theta1 - (D1 - T11) / (1 - theta1) + exp_score * dl_dtheta
}


#' Brent solver for theta_1 M-step under the CTMC link (eps model)
#'
#' Applies Brent's method to \code{.theta1_foc_eps()}. Falls back to the
#' prior-only MLE T11/D1 if the bracket endpoints have the same sign.
#'
#' @param T11,D1,Lg_count,Lg_xsum Sufficient statistics; see
#'   \code{.theta1_foc_eps}.
#' @param theta_seed Initial estimate for theta_1.
#' @param theta_cap Upper bound for theta_1.
#' @param eps_lo Lower bound for theta_1 (default 1e-6).
#' @param tol Convergence tolerance on the FOC (default 1e-10).
#' @param max_iter Maximum Brent iterations (default 100).
#' @return Scalar theta_1, constrained to [eps_lo, theta_cap - eps_lo].
#' @keywords internal
.m_step_theta1_eps_brent <- function(T11, D1, Lg_count, Lg_xsum,
                                     theta_seed, theta_cap,
                                     eps_lo = 1e-6, tol = 1e-10,
                                     max_iter = 100L) {
  lo <- eps_lo
  hi <- theta_cap - eps_lo

  f_lo <- .theta1_foc_eps(lo, T11, D1, Lg_count, Lg_xsum)
  f_hi <- .theta1_foc_eps(hi, T11, D1, Lg_count, Lg_xsum)

  if (!is.finite(f_lo) || !is.finite(f_hi) || (f_lo * f_hi > 0)) {
    return(max(eps_lo, min(theta_seed, theta_cap - eps_lo)))
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
    s  <- max(eps_lo, min(s, theta_cap - eps_lo))
    fs <- .theta1_foc_eps(s, T11, D1, Lg_count, Lg_xsum)
    if (!is.finite(fs)) {
      s  <- (a + b) / 2
      fs <- .theta1_foc_eps(s, T11, D1, Lg_count, Lg_xsum)
    }
    d_pt <- c_pt; c_pt <- b; fc <- fb
    if (fa * fs < 0) { b <- s; fb <- fs } else { a <- s; fa <- fs }
    if (abs(fa) < abs(fb)) {
      tmp <- a; a <- b; b <- tmp; tmp <- fa; fa <- fb; fb <- tmp
    }
  }
  max(eps_lo, min(b, theta_cap - eps_lo))
}


#' M-step for the eps (Spec I) model
#'
#' @param suff Named list of sufficient statistics from \code{e_step_eps()}.
#'   Required entries: C0, C1, D0, D1, T01, T11, M, Eps_num, Eps_den,
#'   Lg_count, Lg_xsum, cat_d_marginal_c, cat_d_marginal_w,
#'   cat_d_trans_curr, cat_d_trans_prev, cat_d_trans_w.
#' @param total_weight Sum of survey weights.
#' @param stationary Logical; if TRUE, impose stationarity on alpha.
#' @param linked Logical; if TRUE, use CTMC link (lambda = f(theta)).
#' @param theta_cap Maximum value for theta1 / theta0 (default 0.999).
#' @param pi_cap Maximum value for pi (default 0.49).
#' @param eps_cap Maximum value for eps (default 0.95).
#' @param eps_floor Minimum value for eps (default 1e-4).
#' @return Named list: alpha, theta0, theta1, pi, eps, lambda_g, lambda_d.
#'   (No sigma2_g.)
#' @references TeX: \emph{EM tenure epsilon.tex}, Section 4.3.
#' @examples
#' \dontrun{
#' # Typically called inside em_fit_tenure_eps(); direct use:
#' suff <- e_step_eps(df, init_params_eps(df))$suff
#' m_step_eps(suff, total_weight = nrow(df), stationary = FALSE)
#' }
#' @export
m_step_eps <- function(suff, total_weight,
                       stationary = FALSE,
                       linked = FALSE,
                       theta_cap = 0.999,
                       pi_cap = 0.49,
                       eps_cap = 0.95,
                       eps_floor = 1e-4) {

  # --- Misclassification (unchanged) ---
  pi_hat <- suff$M / (3 * total_weight)
  pi_hat <- .bound01(pi_hat, eps = 0)
  pi_hat <- min(pi_hat, pi_cap)

  # --- eps update (closed form) ---
  eps_hat <- if (suff$Eps_den > 0) suff$Eps_num / suff$Eps_den else eps_floor
  eps_hat <- .clip_param(eps_hat, eps_floor, eps_cap)

  # --- theta transition block ---
  theta1_seed <- if (suff$D1 > 0) {
    .bound01(suff$T11 / suff$D1, eps = 1e-6)
  } else {
    0.9
  }
  stationary_theta <- NULL
  if (stationary) {
    stationary_theta <- .m_step_theta_stationary_eps(
      suff, linked = linked, theta_cap = theta_cap)
    theta0 <- stationary_theta$theta0
    theta1 <- stationary_theta$theta1
  } else if (linked) {
    theta1 <- .m_step_theta1_eps_brent(
      T11 = suff$T11, D1 = suff$D1,
      Lg_count = suff$Lg_count, Lg_xsum = suff$Lg_xsum,
      theta_seed = theta1_seed, theta_cap = theta_cap
    )
  } else {
    theta1 <- .clip_param(theta1_seed, 1e-6, theta_cap)
  }

  # --- theta0 (unless jointly updated under stationarity) ---
  theta0_seed <- if (suff$D0 > 0) {
    .bound01(suff$T01 / suff$D0, eps = 1e-6)
  } else {
    0.1
  }
  if (stationary) {
    # Already updated jointly above.
  } else if (linked) {
    theta0 <- .m_step_theta0_brent_discrete(
      T01      = suff$T01, D0 = suff$D0,
      cat_marg = suff$cat_d_marginal_c,
      w_marg   = suff$cat_d_marginal_w,
      cat_curr = suff$cat_d_trans_curr,
      cat_prev = suff$cat_d_trans_prev,
      w_trans  = suff$cat_d_trans_w,
      theta_seed = theta0_seed, theta_cap = theta_cap
    )
  } else {
    theta0 <- .clip_param(theta0_seed, 1e-6, theta_cap)
  }

  # --- alpha ---
  if (stationary) {
    alpha <- theta0 / (theta0 + 1 - theta1)
  } else {
    alpha <- .bound01(suff$C1 / (suff$C1 + suff$C0), eps = 1e-6)
  }

  # --- Exponential rates ---
  if (linked) {
    lambda_g <- ctmc_lambda_from_persistence(theta1)
    lambda_d <- ctmc_lambda_from_transition(theta0)
  } else {
    # lambda_g closed form: Exponential MLE inverse-mean
    lambda_g <- if (suff$Lg_xsum > 0 && suff$Lg_count > 0) {
      .clip_param(suff$Lg_count / suff$Lg_xsum, 1e-6, 50)
    } else {
      0.01
    }

    lambda_d <- .m_step_lambda_d_brent(
      cat_marg = suff$cat_d_marginal_c,
      w_marg   = suff$cat_d_marginal_w,
      cat_curr = suff$cat_d_trans_curr,
      cat_prev = suff$cat_d_trans_prev,
      w_trans  = suff$cat_d_trans_w
    )
  }

  out <- list(
    alpha    = alpha,
    theta0   = theta0,
    theta1   = theta1,
    pi       = pi_hat,
    eps      = eps_hat,
    lambda_g = lambda_g,
    lambda_d = lambda_d
  )
  attr(out, "stationary_theta_diagnostics") <- stationary_theta
  out
}
