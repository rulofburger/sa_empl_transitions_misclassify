# ==============================================================================
# EM-baseline-ext: EM/GEM driver for Extension I (observable heterogeneity)
# Created: 2026-05-06
#
# Top-level entry point for fitting the covariate extension. Handles
# parameter initialisation, the GEM iteration loop, convergence checking,
# and result packaging.
#
# TeX ref: EM baseline.tex Section 5.
# ==============================================================================

#' Initialise parameters for the covariate extension
#'
#' Returns zero-vector starting values for \eqn{\beta_0} and \eqn{\beta_1},
#' corresponding to intercept-only probit (all individuals start with the
#' same transition probability). The intercept term is set so that the
#' implied probabilities match sensible baseline priors (θ₀ ≈ 0.10,
#' θ₁ ≈ 0.90).
#'
#' @param p Integer: number of columns in the design matrix X (including intercept).
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @return Named list with \code{beta0} (length-p), \code{beta1} (length-p),
#'   and optionally \code{pi}.
#' @examples
#' init_params_covariates(p = 4L, model_type = "symmetric")
#' @export
init_params_covariates <- function(p, model_type = "symmetric") {
  if (!model_type %in% c("symmetric", "none"))
    stop("init_params_covariates: model_type must be 'symmetric' or 'none'")
  if (!is.numeric(p) || length(p) != 1L || p < 1L)
    stop("init_params_covariates: p must be a positive integer")

  # Start: intercept set to match baseline defaults; other coefficients zero
  beta0 <- rep(0, p)
  beta1 <- rep(0, p)
  beta0[1L] <- qnorm(0.10)  # intercept -> theta0_i ≈ 0.10 for average individual
  beta1[1L] <- qnorm(0.90)  # intercept -> theta1_i ≈ 0.90

  params <- list(beta0 = beta0, beta1 = beta1)
  if (model_type == "symmetric") params$pi <- 0.05

  params
}


#' Fit the covariate extension EM/GEM model
#'
#' Iterates the GEM algorithm (E-step + one IRLS step per M-step) for the
#' observable heterogeneity extension until the relative change in the
#' observed-data log-likelihood drops below \code{tol}.
#'
#' TeX ref: EM baseline.tex Section 5.
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric).
#' @param X N × p design matrix from \code{\link{prepare_covariate_matrix}}.
#' @param model_type Character: \code{"symmetric"} (default) or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}). If \code{TRUE}, individual
#'   \eqn{\alpha_i} is derived from the stationarity restriction.
#' @param params0 Optional named list of starting parameters
#'   (\code{beta0}, \code{beta1}, and optionally \code{pi}). If \code{NULL},
#'   \code{\link{init_params_covariates}} is used.
#' @param max_iter Maximum number of GEM iterations (default 500L).
#' @param tol Convergence tolerance on relative change in log-likelihood
#'   (default 1e-6).
#' @param pi_cap Upper bound for \eqn{\pi} (default 0.49).
#' @param verbose Integer verbosity: 0 = silent, 1 = convergence summary
#'   (default 1L).
#' @return List with:
#'   \describe{
#'     \item{params}{Named list of final parameter estimates.}
#'     \item{loglik}{Final observed-data log-likelihood.}
#'     \item{history}{Data frame: iter, loglik, and one column per scalar
#'       parameter (intercepts of beta0/beta1 plus pi if applicable).}
#'     \item{converged}{Logical.}
#'     \item{iterations}{Integer: number of completed iterations.}
#'     \item{gamma}{N × 8 matrix of final responsibilities.}
#'     \item{model_type}{Character.}
#'     \item{stationary}{Logical.}
#'   }
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   X      <- matrix(1, nrow = 50, ncol = 1)
#'   em_fit_covariates(df, X, verbose = 0L)
#' }
#' @export
em_fit_covariates <- function(df,
                               X,
                               model_type = "symmetric",
                               stationary = TRUE,
                               params0    = NULL,
                               max_iter   = 500L,
                               tol        = 1e-6,
                               pi_cap     = 0.49,
                               verbose    = 1L) {
  if (!model_type %in% c("symmetric", "none"))
    stop("em_fit_covariates: model_type must be 'symmetric' or 'none'")
  .validate_panel_df(df)
  X_transition <- .as_transition_design(X, nrow(df))
  if (stationary && !.transition_design_is_time_invariant(X_transition))
    stop("em_fit_covariates: use stationary=FALSE with time-varying covariates")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1L)
    stop("em_fit_covariates: max_iter must be a positive integer")

  p      <- ncol(X_transition$X12)
  params <- params0 %||% init_params_covariates(p, model_type = model_type)
  coef_names <- colnames(X_transition$X12)
  if (is.null(coef_names)) coef_names <- paste0("x", seq_len(p))
  names(params$beta0) <- coef_names
  names(params$beta1) <- coef_names
  entry_active <- attr(X_transition$X12, "entry_active")
  if (is.null(entry_active)) entry_active <- rep(TRUE, p)
  persistence_active <- attr(X_transition$X12, "persistence_active")
  if (is.null(persistence_active)) persistence_active <- rep(TRUE, p)

  # Validate required fields
  for (nm in c("beta0", "beta1")) {
    if (is.null(params[[nm]]) || length(params[[nm]]) != p)
      stop(sprintf(
        "em_fit_covariates: params$%s must be a length-%d numeric vector", nm, p
      ))
  }
  if (model_type == "symmetric" && !is.null(params$pi)) {
    if (!is.numeric(params$pi) || length(params$pi) != 1L ||
        params$pi < 0 || params$pi >= 0.5)
      stop("em_fit_covariates: params$pi must be a scalar in [0, 0.5) for model_type='symmetric'")
  }
  params$beta0[!entry_active] <- 0
  params$beta1[!persistence_active] <- 0
  if (!stationary && is.null(params$alpha)) params$alpha <- 0.5

  # ---- Validate data quality once ------------------------------------------
  for (.col in c("y1", "y2", "y3")) {
    if (!all(df[[.col]] %in% c(0, 1)))
      stop(sprintf("em_fit_covariates: column '%s' must be binary (0/1 only)", .col))
  }
  if (anyNA(df$weight) || any(df$weight <= 0))
    stop("em_fit_covariates: weight must be non-NA and strictly positive")

  # Normalize weights for numerical conditioning. This does not change the
  # estimates or responsibilities; returned likelihoods retain the input scale.
  weight_scale <- mean(df$weight)
  df_fit <- df
  df_fit$weight <- df$weight / weight_scale
  converged <- FALSE

  # Include the starting point so monotonicity is directly testable.
  history_rows <- vector("list", max_iter + 1L)

  # Precompute static mismatch indicator matrix (N x 8) for symmetric model.
  # outer(s_t, h_t, "!=") is constant across iterations; precomputing avoids
  # recomputing 3 N x 8 outer products per E-step.
  hmat_static <- .hmat_cache()
  h1_s <- hmat_static[, 1L]; h2_s <- hmat_static[, 2L]; h3_s <- hmat_static[, 3L]
  mm_static <- if (model_type == "symmetric") {
    outer(as.integer(df_fit$y1), h1_s, "!=") +
    outer(as.integer(df_fit$y2), h2_s, "!=") +
    outer(as.integer(df_fit$y3), h3_s, "!=")
  } else {
    NULL
  }

  current <- e_step_covariates(
    df_fit, X, params, model_type, validate = FALSE,
    mm_precomp = mm_static, stationary = stationary
  )
  ll_current <- current$loglik
  history_rows[[1L]] <- c(
    iter = 0, loglik = ll_current * weight_scale,
    beta0_intercept = params$beta0[1L], beta1_intercept = params$beta1[1L],
    pi = if (model_type == "symmetric") params$pi else NA_real_,
    step_fraction = NA_real_
  )

  for (iter in seq_len(max_iter)) {
    candidate <- m_step_covariates(
      suff       = current$suff,
      X          = X,
      params_old = params,
      model_type = model_type,
      stationary = stationary,
      pi_cap     = pi_cap,
      entry_active = entry_active,
      persistence_active = persistence_active
    )
    mdiag <- attr(candidate, "mstep")
    attr(candidate, "mstep") <- NULL

    candidate_estep <- e_step_covariates(
      df_fit, X, candidate, model_type, validate = FALSE,
      mm_precomp = mm_static, stationary = stationary
    )
    ll_candidate <- candidate_estep$loglik

    # Observed-likelihood safeguard. A valid GEM Q-step should already pass;
    # this protects against optimizer and floating-point failures.
    candidate_full <- candidate
    obs_fraction <- 1
    while (ll_candidate < ll_current - 1e-10 && obs_fraction > 2^-20) {
      obs_fraction <- obs_fraction / 2
      candidate <- .interpolate_cov_params(params, candidate_full, obs_fraction,
                                            entry_active, persistence_active)
      candidate_estep <- e_step_covariates(
        df_fit, X, candidate, model_type, validate = FALSE,
        mm_precomp = mm_static, stationary = stationary
      )
      ll_candidate <- candidate_estep$loglik
    }
    if (ll_candidate < ll_current - 1e-10) {
      candidate <- params
      candidate_estep <- current
      ll_candidate <- ll_current
      obs_fraction <- 0
    }

    old_vec <- c(params$beta0[entry_active], params$beta1[persistence_active],
                 params$pi %||% numeric(0), params$alpha %||% numeric(0))
    new_vec <- c(candidate$beta0[entry_active], candidate$beta1[persistence_active],
                 candidate$pi %||% numeric(0), candidate$alpha %||% numeric(0))
    param_change <- max(abs(new_vec - old_vec))
    improvement <- ll_candidate - ll_current
    rel_improvement <- improvement / (abs(ll_current) + 1e-16)

    params <- candidate
    current <- candidate_estep
    ll_current <- ll_candidate
    history_rows[[iter + 1L]] <- c(
      iter = iter, loglik = ll_current * weight_scale,
      beta0_intercept = params$beta0[1L], beta1_intercept = params$beta1[1L],
      pi = if (model_type == "symmetric") params$pi else NA_real_,
      step_fraction = obs_fraction * (mdiag$step_fraction %||% 1)
    )

    if (improvement >= -1e-10 && rel_improvement < tol &&
        param_change < max(sqrt(tol), 1e-6)) {
      converged <- TRUE
      break
    }
  }

  if (verbose >= 1L) {
    status <- if (converged) "converged" else "did NOT converge"
    cat(sprintf(
      "em_fit_covariates [%s, stationary=%s]: %s in %d iters. loglik = %.4f\n",
      model_type, stationary, status, iter, ll_current * weight_scale
    ))
    cat(sprintf(
      "  beta0 intercept = %.4f (theta0 ~ %.4f)\n",
      params$beta0[1L], pnorm(params$beta0[1L])
    ))
    cat(sprintf(
      "  beta1 intercept = %.4f (theta1 ~ %.4f)\n",
      params$beta1[1L], pnorm(params$beta1[1L])
    ))
    if (model_type == "symmetric")
      cat(sprintf("  pi = %.4f\n", params$pi))
  }

  if (!converged && verbose >= 1L)
    warning(sprintf(
      "em_fit_covariates [%s, stationary=%s]: did not converge in %d iterations.",
      model_type, stationary, max_iter
    ))

  history_df <- as.data.frame(do.call(rbind, history_rows[seq_len(iter + 1L)]))

  list(
    params     = params,
    loglik     = current$loglik * weight_scale,
    history    = history_df,
    converged  = converged,
    iterations = iter,
    gamma      = current$gamma,
    model_type = model_type,
    stationary = stationary
  )
}


#' Evaluate the observed-data log-likelihood for the covariate extension
#'
#' @param df Data frame with columns y1, y2, y3, weight.
#' @param X N × p design matrix.
#' @param params Named list with beta0, beta1, pi (optional).
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @return Scalar observed-data log-likelihood.
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   X      <- matrix(1, nrow = 50, ncol = 1)
#'   params <- init_params_covariates(p = 1L)
#'   compute_observed_loglik_covariates(df, X, params)
#' }
#' @export
compute_observed_loglik_covariates <- function(df, X, params,
                                               model_type = "symmetric",
                                               stationary = TRUE) {
  e_step_covariates(df, X, params, model_type, validate = TRUE,
                    stationary = stationary)$loglik
}
