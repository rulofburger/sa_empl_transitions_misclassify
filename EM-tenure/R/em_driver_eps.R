# ==============================================================================
# EM-tenure: EM driver for the eps (Spec I) model
# ==============================================================================
# Created: 2026-04-30
#
# Top-level driver em_fit_tenure_eps() that orchestrates the eps E-step /
# M-step iteration with monotonicity guards, convergence monitoring, and
# parameter history logging. Mirrors em_fit_tenure_rho() in structure.
#
# Only supports discrete_timegap = TRUE.
#
# Companion: documents/EM tenure epsilon.tex.
# ==============================================================================

#' Fit the EM model with spell-pair contamination (eps / Spec I model)
#'
#' Estimates the seven-parameter eps model
#' \eqn{(\alpha, \theta_0, \theta_1, \pi, \varepsilon, \lambda_g, \lambda_d)}
#' via standard EM. Tenure emissions use the spell-pair point-mass + Exp
#' contamination mixture; timegap emissions use the interval-censored
#' Exp(\eqn{\lambda_d}) machinery from the base model.
#'
#' @param df Data frame with columns: y1-y3, tenure1-tenure3,
#'   timegap_cat1-timegap_cat3, weight.
#' @param params0 Optional named list of starting parameters (must include
#'   eps; see \code{init_params_eps()}).
#' @param stationary Logical; if TRUE, impose stationarity on alpha.
#' @param linked Logical; if TRUE, impose the CTMC link constraints
#'   \eqn{\lambda_g = -\log(\theta_1)/\Delta} and
#'   \eqn{\lambda_d = -\log(1-\theta_0)/\Delta}, where \eqn{\Delta = 0.25}
#'   years (quarterly survey spacing).
#' @param max_iter Maximum EM iterations (default 500).
#' @param tol Convergence tolerance on relative change in log-likelihood
#'   (default 1e-8).
#' @param theta_cap Maximum theta (default 0.999).
#' @param pi_cap Maximum pi (default 0.49).
#' @param eps_cap Maximum eps (default 0.95).
#' @param eps_floor Minimum eps (default 1e-4).
#' @param verbose Integer; 0 = silent, 1 = final summary, 2 = every iteration.
#' @return List with: \code{params}, \code{loglik}, \code{history},
#'   \code{converged}, \code{iterations}, \code{gamma}.
#' @references TeX: \emph{EM tenure epsilon.tex}.
#' @examples
#' \dontrun{
#' fit <- em_fit_tenure_eps(df_qlfs, max_iter = 200L, verbose = 1L)
#' print(fit$params)
#' print(fit$converged)
#' }
#' @export
em_fit_tenure_eps <- function(df,
                              params0 = NULL,
                              stationary = FALSE,
                              linked = FALSE,
                              max_iter = 500L,
                              tol = 1e-8,
                              theta_cap = 0.999,
                              pi_cap = 0.49,
                              eps_cap = 0.95,
                              eps_floor = 1e-4,
                              verbose = 1L) {
  # --- Validate inputs once (before the main loop) ---
  # validate_df_eps provides comprehensive column/type/NA/bounds checks.
  # Called here so e_step_eps() can skip df checks (check_df = FALSE)
  # on every iteration, saving O(N) traversals * max_iter per run.
  validate_df_eps(df)
  total_weight <- sum(df$weight)

  # --- Initialise parameters ---
  params <- if (!is.null(params0)) params0 else init_params_eps(df, linked = linked)

  if (is.null(params$eps)) {
    params$eps <- 0.20
  }
  # Drop sigma2_g if present in the warm-start (eps model has no sigma).
  params$sigma2_g <- NULL
  params$sigma2_d <- NULL

  if (linked) {
    params$lambda_g <- ctmc_lambda_from_persistence(params$theta1)
    params$lambda_d <- ctmc_lambda_from_transition(params$theta0)
  }

  # --- History storage ---
  history   <- vector("list", max_iter)
  ll_vec    <- numeric(max_iter)
  converged <- FALSE

  # cached_estep: holds the guard E-step result from the previous iteration so
  # that the next iteration reuses it instead of re-running the E-step.
  # This avoids 2x E-step cost per iteration when stationary || linked.
  cached_estep <- NULL

  for (iter in seq_len(max_iter)) {
    # --- E-step (reuse cached guard result when available) ---
    if (!is.null(cached_estep)) {
      estep_out    <- cached_estep
      cached_estep <- NULL
    } else {
      estep_out <- e_step_eps(df, params, check_df = FALSE)
    }
    ll_new       <- estep_out$loglik
    ll_vec[iter] <- ll_new

    # --- Convergence check ---
    if (iter > 1) {
      ll_change <- ll_new - ll_vec[iter - 1]
      if (ll_change < -1e-8) {
        warning(sprintf(
          "EM-eps iter %d: LL decreased by %.6e (from %.6e to %.6e)",
          iter, ll_change, ll_vec[iter - 1], ll_new
        ))
      }
      rel_change <- abs(ll_change) / (abs(ll_vec[iter - 1]) + 1e-16)
      if (rel_change < tol) {
        converged <- TRUE
        if (verbose >= 1) {
          message(sprintf(
            "EM-eps converged at iteration %d (rel_change = %.2e)",
            iter, rel_change
          ))
        }
        history[[iter]] <- unlist(params)
        break
      }
    }

    # --- M-step ---
    params_prev <- params
    params <- m_step_eps(
      suff         = estep_out$suff,
      total_weight = total_weight,
      stationary   = stationary,
      linked       = linked,
      theta_cap    = theta_cap,
      pi_cap       = pi_cap,
      eps_cap      = eps_cap,
      eps_floor    = eps_floor
    )

    # --- Monotonicity guard (all variants) ---
    # Runs for every iteration regardless of constraints. For constrained
    # models (stationary || linked), the result is cached and reused as the
    # next iteration's E-step, eliminating the duplicate E-step cost.
    estep_check <- e_step_eps(df, params, check_df = FALSE)
    if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
      if (verbose >= 2) {
        message(sprintf(
          "EM-eps iter %d: M-step decreased LL (%.6e -> %.6e); reverting.",
          iter, ll_new, estep_check$loglik
        ))
      }
      params    <- params_prev
      converged <- TRUE   # guard fired: best-known params retained; algorithm
                          # terminated with a valid (non-decreasing) parameter set
      history[[iter]] <- unlist(params)
      break
    }
    # Reuse the guard E-step as the next iteration's starting point (saves one
    # E-step call per iteration for all variants, not just constrained ones).
    cached_estep <- estep_check

    history[[iter]] <- unlist(params)

    if (verbose >= 2) {
      message(sprintf(
        paste0("EM-eps iter %3d | ll = %14.4f | alpha=%.4f theta1=%.4f ",
               "theta0=%.4f pi=%.4f eps=%.4f lam_g=%.4f lam_d=%.4f"),
        iter, ll_new, params$alpha, params$theta1, params$theta0,
        params$pi, params$eps, params$lambda_g, params$lambda_d
      ))
    }
  }

  if (!converged && verbose >= 1) {
    message(sprintf("EM-eps did not converge in %d iterations", max_iter))
  }

  # --- Trim history ---
  n_iter      <- min(iter, max_iter)
  history     <- history[seq_len(n_iter)]
  history_mat <- do.call(rbind, history)
  history_df  <- as.data.frame(history_mat)
  history_df$iteration <- seq_len(n_iter)
  history_df$loglik    <- as.numeric(ll_vec[seq_len(n_iter)])

  if (!converged) {
    # Use cached guard E-step if available; otherwise re-evaluate.
    estep_out <- if (!is.null(cached_estep)) cached_estep else e_step_eps(df, params, check_df = FALSE)
  }

  list(
    params     = params,
    loglik     = estep_out$loglik,
    history    = history_df,
    converged  = converged,
    iterations = n_iter,
    gamma      = estep_out$gamma
  )
}
