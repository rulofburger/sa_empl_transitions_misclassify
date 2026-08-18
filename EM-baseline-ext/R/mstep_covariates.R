# ==============================================================================
# EM-baseline-ext: monotone M/GEM step for observable heterogeneity
# ==============================================================================

# The stationary initial probability depends jointly on the two transition
# equations. Separate IRLS updates therefore need not increase the correct
# complete Q-function. This module maximizes the full coefficient block jointly
# and applies a backtracking guard.

.covariate_q <- function(par, suff, X, entry_active, stationary,
                         gradient = FALSE) {
  X_transition <- .as_transition_design(X)
  p <- ncol(X_transition$X12)
  p0 <- sum(entry_active)
  beta0 <- numeric(p)
  beta0[entry_active] <- par[seq_len(p0)]
  beta1 <- par[p0 + seq_len(p)]

  q <- 0
  grad0 <- numeric(p0)
  grad1 <- numeric(p)
  first <- NULL
  for (nm in c("X12", "X23")) {
    Xt <- X_transition[[nm]]
    st <- suff$transition[[nm]]
    eta0 <- as.vector(Xt %*% beta0)
    eta1 <- as.vector(Xt %*% beta1)
    theta0 <- pmin(pmax(pnorm(eta0), 1e-10), 1 - 1e-10)
    theta1 <- pmin(pmax(pnorm(eta1), 1e-10), 1 - 1e-10)
    fail0 <- pmax(st$eff_w_0 - st$eff_wy_0, 0)
    fail1 <- pmax(st$eff_w_1 - st$eff_wy_1, 0)
    q <- q + sum(st$eff_wy_0 * log(theta0) + fail0 * log1p(-theta0)) +
      sum(st$eff_wy_1 * log(theta1) + fail1 * log1p(-theta1))
    score_eta0 <- dnorm(eta0) *
      (st$eff_wy_0 / theta0 - fail0 / (1 - theta0))
    score_eta1 <- dnorm(eta1) *
      (st$eff_wy_1 / theta1 - fail1 / (1 - theta1))
    grad0 <- grad0 + as.vector(crossprod(Xt[, entry_active, drop = FALSE], score_eta0))
    grad1 <- grad1 + as.vector(crossprod(Xt, score_eta1))
    if (nm == "X12") first <- list(eta0 = eta0, eta1 = eta1,
                                      theta0 = theta0, theta1 = theta1)
  }

  if (stationary) {
    denom <- first$theta0 + 1 - first$theta1
    alpha <- pmin(pmax(first$theta0 / denom, 1e-10), 1 - 1e-10)
    q <- q + sum(suff$init_w1 * log(alpha) +
                 suff$init_w0 * log1p(-alpha))

    score_alpha <- suff$init_w1 / alpha - suff$init_w0 / (1 - alpha)
    score_eta0_init <- score_alpha * ((1 - first$theta1) / denom^2) *
      dnorm(first$eta0)
    score_eta1_init <- score_alpha * (first$theta0 / denom^2) *
      dnorm(first$eta1)
    grad0 <- grad0 + as.vector(crossprod(
      X_transition$X12[, entry_active, drop = FALSE], score_eta0_init))
    grad1 <- grad1 + as.vector(crossprod(X_transition$X12, score_eta1_init))
  }

  # Scaling leaves the maximizer unchanged and improves conditioning with
  # large survey weights.
  q_scale <- max(suff$total_weight, 1)
  if (!gradient) return(q / q_scale)
  c(grad0, grad1) / q_scale
}

.interpolate_cov_params <- function(old, proposed, fraction, entry_active) {
  out <- old
  out$beta0 <- old$beta0 + fraction * (proposed$beta0 - old$beta0)
  out$beta0[!entry_active] <- 0
  out$beta1 <- old$beta1 + fraction * (proposed$beta1 - old$beta1)
  if (!is.null(proposed$pi))
    out$pi <- old$pi + fraction * (proposed$pi - old$pi)
  if (!is.null(proposed$alpha))
    out$alpha <- old$alpha + fraction * (proposed$alpha - old$alpha)
  out
}

#' Monotone M/GEM step for the covariate extension
#'
#' Maximizes the transition-coefficient block jointly. Under stationarity the
#' objective includes the individual initial-state contribution, which couples
#' beta0 and beta1. A backtracking line search rejects coefficient proposals
#' that do not increase the complete Q-function.
#'
#' @param suff Sufficient statistics returned by e_step_covariates().
#' @param X Covariate matrix.
#' @param params_old Current parameter list.
#' @param model_type "symmetric" or "none".
#' @param stationary Whether initial employment is restricted to stationarity.
#' @param pi_cap Upper bound for the symmetric error probability.
#' @param entry_active Logical vector identifying coefficients active in the
#'   entry equation. Defaults to the entry_active attribute on X.
#' @return Updated parameter list, with diagnostics in attribute "mstep".
m_step_covariates <- function(suff, X, params_old, model_type = "symmetric",
                               stationary = TRUE, pi_cap = 0.49,
                               entry_active = attr(X, "entry_active")) {
  if (!model_type %in% c("symmetric", "none"))
    stop("m_step_covariates: model_type must be 'symmetric' or 'none'")

  X_transition <- .as_transition_design(X)
  p <- ncol(X_transition$X12)
  if (is.null(entry_active)) entry_active <- rep(TRUE, p)
  if (!is.logical(entry_active) || length(entry_active) != p || !any(entry_active))
    stop("m_step_covariates: entry_active must be a logical vector of length ncol(X)")

  old <- params_old
  old$beta0[!entry_active] <- 0
  par_old <- c(old$beta0[entry_active], old$beta1)

  opt <- optim(
    par = par_old,
    fn = function(z) -.covariate_q(z, suff, X, entry_active, stationary),
    gr = function(z) -.covariate_q(z, suff, X, entry_active, stationary,
                                    gradient = TRUE),
    method = "BFGS",
    # A handful of joint BFGS iterations is sufficient for a GEM step. The
    # outer EM loop and ascent safeguards finish the optimization without an
    # expensive fully-converged inner solve at every iteration.
    control = list(maxit = 12L, reltol = 1e-8)
  )

  proposed <- old
  p0 <- sum(entry_active)
  proposed$beta0[entry_active] <- opt$par[seq_len(p0)]
  proposed$beta0[!entry_active] <- 0
  proposed$beta1 <- opt$par[p0 + seq_len(p)]

  if (model_type == "symmetric") {
    pi_new <- suff$M / max(3 * suff$total_weight, 1e-10)
    proposed$pi <- min(max(pi_new, 0), pi_cap)
  }
  if (!stationary) {
    alpha_new <- suff$C1 / max(suff$C1 + suff$C0, 1e-10)
    proposed$alpha <- pmin(pmax(alpha_new, 1e-10), 1 - 1e-10)
  }

  q_old <- .covariate_q(par_old, suff, X, entry_active, stationary)
  par_new <- c(proposed$beta0[entry_active], proposed$beta1)
  q_new <- .covariate_q(par_new, suff, X, entry_active, stationary)

  proposal_full <- proposed
  fraction <- 1
  while ((!is.finite(q_new) || q_new < q_old - 1e-12) && fraction > 2^-20) {
    fraction <- fraction / 2
    proposed <- .interpolate_cov_params(old, proposal_full, fraction, entry_active)
    par_new <- c(proposed$beta0[entry_active], proposed$beta1)
    q_new <- .covariate_q(par_new, suff, X, entry_active, stationary)
  }
  if (!is.finite(q_new) || q_new < q_old - 1e-12) {
    proposed <- old
    q_new <- q_old
    fraction <- 0
  }

  attr(proposed, "mstep") <- list(
    q_old = q_old, q_new = q_new, step_fraction = fraction,
    optim_convergence = opt$convergence, optim_message = opt$message
  )
  proposed
}
