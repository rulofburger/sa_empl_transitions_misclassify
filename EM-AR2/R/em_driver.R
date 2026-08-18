# ==============================================================================
# EM-AR2: EM driver — parameter initialisation and main fitting function
# Created: 2026-05-05
# ==============================================================================
# Orchestrates E-step / M-step iteration with convergence monitoring and
# parameter history logging.
#
# TeX reference: EM algorithm convergence and initialisation.
# ==============================================================================

#' Initialise AR(2) EM parameters from moment-based estimates
#'
#' Provides crude moment-based starting values derived from the observed data.
#' Uses observed transition frequencies to seed theta0 and theta1, and starts
#' theta01 and theta10 near zero (AR(1) baseline). pi is initialised to a
#' small positive value to allow misclassification in the first E-step.
#'
#' The initial pair distribution alpha(h1,h2) is set to the stationary
#' distribution of the initial theta estimates. During EM iterations, alpha
#' is updated as a free parameter (Baum-Welch approach).
#'
#' @param df Data frame with columns y1, y2, y3, y4 (binary employment).
#' @param pi_init Starting misclassification probability (default 0.05).
#' @return Named list: alpha, theta0, theta01, theta1, theta10, pi.
#' @examples
#' \dontrun{
#' df <- data.frame(y1=c(0,1,0), y2=c(0,0,1), y3=c(1,0,0), y4=c(1,1,0), weight=c(1,1,1))
#' params <- init_params_ar2(df)
#' }
#' @export
init_params_ar2 <- function(df, pi_init = 0.05) {
  y <- as.matrix(df[, c("y1", "y2", "y3", "y4")])
  w <- df$weight %||% rep(1, nrow(df))
  transition_w <- rep(w, 3L)

  # Observed transition frequencies (crude moments, ignoring misclassification)
  # Job-finding: fraction of non-employed who become employed in the next wave
  from0 <- as.vector(y[, -4, drop = FALSE] == 0)
  to1 <- as.vector(y[, -1, drop = FALSE] == 1)
  theta0_init <- if (any(from0))
    weighted.mean(to1[from0], transition_w[from0]) else 0.1

  # Separation: fraction of employed who become non-employed in the next wave
  from1 <- !from0
  to0 <- !to1
  theta1_init <- if (any(from1))
    weighted.mean(to0[from1], transition_w[from1]) else 0.1

  theta0  <- .bound01(theta0_init, eps = 0.01)
  theta01 <- 0.05   # Start near zero: no duration-dependence by default
  theta1  <- .bound01(theta1_init, eps = 0.01)
  theta10 <- 0.05   # Start near zero: no duration-dependence by default

  # Clamp theta01/theta10 so joint constraints are satisfied
  theta01 <- min(theta01, (1 - theta0) - 1e-6)
  theta10 <- min(theta10, (1 - theta1) - 1e-6)

  # Initialize alpha from stationary distribution of initial theta estimates
  alpha <- tryCatch(
    stationary_ar2(theta0, theta01, theta1, theta10),
    error = function(e) c(`00`=0.25, `10`=0.25, `01`=0.25, `11`=0.25)
  )

  list(
    alpha   = alpha,
    theta0  = theta0,
    theta01 = theta01,
    theta1  = theta1,
    theta10 = theta10,
    pi      = .bound01(pi_init, eps = 1e-4)
  )
}

#' Initialise AR(2) EM parameters from existing MLE estimates
#'
#' Warm-starts the EM from a list of parameter estimates (e.g., from a prior
#' run or from the MLE). Useful for exact recovery comparisons.
#'
#' @param theta0 Baseline job-finding probability.
#' @param theta01 Duration-dependence increment for job-finding.
#' @param theta1 Baseline separation probability.
#' @param theta10 Duration-dependence increment for separation.
#' @param pi Misclassification probability (default 0.05).
#' @param pi0 False non-employment probability (asymmetric model).
#' @param pi1 False employment probability (asymmetric model).
#' @param alpha Optional named length-4 initial pair distribution.
#'   If NULL, the stationary distribution is computed from theta.
#' @return Named list suitable for passing to em_fit_ar2() as params0.
#' @export
params_from_values <- function(theta0, theta01, theta1, theta10,
                                pi = 0.05, pi0 = NULL, pi1 = NULL,
                                alpha = NULL) {
  theta0  <- .bound01(theta0)
  theta01 <- max(theta01, 0)
  theta1  <- .bound01(theta1)
  theta10 <- max(theta10, 0)

  if (theta0 + theta01 >= 1) {
    stop("params_from_values: theta0 + theta01 must be < 1. ",
         sprintf("Got %.4f + %.4f = %.4f.", theta0, theta01, theta0 + theta01))
  }
  if (theta1 + theta10 >= 1) {
    stop("params_from_values: theta1 + theta10 must be < 1. ",
         sprintf("Got %.4f + %.4f = %.4f.", theta1, theta10, theta1 + theta10))
  }

  if (is.null(alpha)) {
    alpha <- tryCatch(
      stationary_ar2(theta0, theta01, theta1, theta10),
      error = function(e) c(`00`=0.25, `10`=0.25, `01`=0.25, `11`=0.25)
    )
  }

  out <- list(
    alpha   = alpha,
    theta0  = theta0,
    theta01 = theta01,
    theta1  = theta1,
    theta10 = theta10,
    pi      = .bound01(pi)
  )
  if (!is.null(pi0)) out$pi0 <- .bound01(pi0)
  if (!is.null(pi1)) out$pi1 <- .bound01(pi1)
  out
}

#' Fit the AR(2) EM model for misclassified employment transitions
#'
#' Estimates the AR(2) Markov model for binary employment transitions observed
#' with misclassification error. Iterates E-step and M-step until the relative
#' change in log-likelihood is below tol or max_iter is reached.
#'
#' @param df Data frame with columns y1, y2, y3, y4 (binary employment)
#'   and weight (positive survey weights).
#' @param params0 Optional named list of starting parameters. If NULL, uses
#'   init_params_ar2(df).
#' @param estimate_pi Logical; if TRUE estimate pi from data (default TRUE).
#' @param fixed_pi Scalar; pi value when estimate_pi=FALSE (default 0).
#' @param asymmetric Logical; if TRUE use separate pi0, pi1 (default FALSE).
#' @param max_iter Maximum EM iterations (default 500).
#' @param tol Convergence tolerance on relative log-likelihood change
#'   (default 1e-8).
#' @param param_tol Convergence tolerance on the largest EM parameter update.
#'   Requiring both criteria avoids premature convergence on a flat likelihood.
#' @param pi_cap Maximum misclassification probability allowed (default 0.49).
#' @param verbose Integer; 0=silent, 1=final summary, 2=every 10 iterations,
#'   3=every iteration (default 1).
#' @param collapse_cells Logical; collapse identical observed histories before
#'   EM. This is likelihood-equivalent and strongly recommended for large data.
#' @return Named list:
#'   - params: final parameter estimates (named list)
#'   - loglik: final observed-data log-likelihood
#'   - history: data.frame of parameter values at each iteration
#'   - converged: logical
#'   - iterations: integer
#'   - gamma: N x 16 matrix of final responsibilities
#' @examples
#' \dontrun{
#' df <- data.frame(y1=rbinom(200,1,.3), y2=rbinom(200,1,.3),
#'                  y3=rbinom(200,1,.3), y4=rbinom(200,1,.3), weight=rep(1,200))
#' fit <- em_fit_ar2(df, estimate_pi = FALSE, fixed_pi = 0, max_iter = 50L, verbose = 0L)
#' fit$params
#' fit$loglik
#' }
#' @export
em_fit_ar2 <- function(df,
                        params0      = NULL,
                        estimate_pi  = TRUE,
                        fixed_pi     = 0,
                        asymmetric   = FALSE,
                        max_iter     = 500L,
                        tol          = 1e-8,
                        param_tol    = 1e-8,
                        pi_cap       = 0.49,
                        verbose      = 1L,
                        collapse_cells = TRUE) {
  # --- Validate inputs ---
  required_cols <- c("y1", "y2", "y3", "y4", "weight")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("em_fit_ar2: missing columns in df: ", paste(missing_cols, collapse = ", "))
  }

  y_original <- as.matrix(df[, c("y1", "y2", "y3", "y4")])
  if (any(is.na(y_original))) stop("em_fit_ar2: NA values in y1-y4 are not allowed.")
  if (!all(y_original %in% c(0L, 1L)))
    stop("em_fit_ar2: ymat must contain only 0/1 values.")
  if (any(is.na(df$weight))) stop("em_fit_ar2: NA values in weights are not allowed.")
  if (any(!is.finite(df$weight)) || any(df$weight <= 0))
    stop("em_fit_ar2: all weights must be finite and positive.")

  N_original <- nrow(df)
  df_em <- if (isTRUE(collapse_cells)) collapse_ar2_cells(df) else df

  # Prepare matrices
  ymat <- as.matrix(df_em[, c("y1", "y2", "y3", "y4")])
  w    <- df_em$weight
  N    <- nrow(ymat)

  if (any(is.na(ymat)))             stop("em_fit_ar2: NA values in y1-y4 are not allowed.")
  if (!all(ymat %in% c(0L, 1L)))   stop("em_fit_ar2: ymat must contain only 0/1 values.")
  if (any(is.na(w)))               stop("em_fit_ar2: NA values in weights are not allowed.")
  if (any(w <= 0))                 stop("em_fit_ar2: all weights must be positive.")

  # Generate latent histories
  hmat <- latent_histories_ar2()

  # Initialise parameters
  if (is.null(params0)) {
    params <- init_params_ar2(df_em)
  } else {
    params <- params0
  }

  # For no-ME variant: fix pi=0
  if (!estimate_pi && !asymmetric) {
    params$pi <- fixed_pi
  }

  if (asymmetric) {
    params$pi0 <- params$pi0 %||% params$pi %||% 0.05
    params$pi1 <- params$pi1 %||% params$pi %||% 0.05
    params$pi <- NULL
  }

  # Mark asymmetric flag inside params so e_step / m_step know which path
  params$asymmetric <- asymmetric

  # --- EM iteration ---
  loglik_prev <- -Inf
  history_rows <- vector("list", max_iter)
  gamma <- NULL
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # E-step
    estep_out <- e_step_ar2(ymat, w, hmat, params)
    loglik     <- estep_out$loglik
    gamma      <- estep_out$gamma
    ss         <- estep_out$sufficient_stats

    # Form the next EM update before checking convergence. A flat likelihood
    # alone is insufficient: weakly identified parameters may still be moving.
    new_params <- m_step_ar2(ss,
      estimate_pi = estimate_pi,
      fixed_pi    = fixed_pi,
      asymmetric  = asymmetric,
      pi_cap      = pi_cap
    )
    new_params$asymmetric <- asymmetric
    tracked <- intersect(c("theta0", "theta01", "theta1", "theta10",
                           "pi", "pi0", "pi1"), names(new_params))
    param_change <- max(c(
      abs(unlist(new_params[tracked]) - unlist(params[tracked])),
      abs(new_params$alpha - params$alpha)
    ))

    # Convergence check (after first iteration)
    if (iter > 1) {
      rel_change <- abs(loglik - loglik_prev) /
        (abs(loglik_prev) + .Machine$double.eps)

      if (verbose >= 3L || (verbose >= 2L && iter %% 10L == 0L)) {
        cat(sprintf("Iter %d: LL = %.4f  rel_change = %.2e\n",
                    iter, loglik, rel_change))
      }

      if (loglik < loglik_prev - 1e-6) {
        warning(sprintf(
          "EM-AR2: log-likelihood decreased at iteration %d (%.4f < %.4f). ",
          iter, loglik, loglik_prev),
          "This may indicate a numerical issue.")
      }

      if (rel_change < tol && param_change < param_tol) {
        if (verbose >= 1L) {
          cat(sprintf("EM-AR2: converged at iteration %d. LL = %.4f\n",
                      iter, loglik))
        }
        # Record final state and break
        history_rows[[iter]] <- .params_to_row(params, loglik, iter)
        converged <- TRUE
        break
      }
    }

    loglik_prev <- loglik

    # Record parameter history
    history_rows[[iter]] <- .params_to_row(params, loglik, iter)

    params <- new_params

    if (iter == max_iter) {
      warning(sprintf(
        "EM-AR2: did not converge after %d iterations. Last LL = %.4f",
        max_iter, loglik))
    }
  }

  iterations <- iter

  # Use gamma and loglik from the final EM loop iteration (already saved)
  loglik_final <- loglik

  if (verbose >= 1L && !converged) {
    message(sprintf("EM-AR2: stopped at max_iter=%d. Final LL = %.4f",
                    max_iter, loglik_final))
  }

  # Build history data frame
  history_rows <- Filter(Negate(is.null), history_rows)
  history <- if (length(history_rows) > 0) {
    do.call(rbind, lapply(history_rows, as.data.frame))
  } else {
    data.frame()
  }

  # Remove diagnostic .p_* fields from final params
  params_clean <- params[!grepl("^\\.p_", names(params))]

  list(
    params     = params_clean,
    loglik     = loglik_final,
    history    = history,
    converged  = converged,
    iterations = iterations,
    gamma      = gamma,
    cells      = df_em,
    n_obs      = N_original,
    estimator  = "collapsed_cell_em",
    diagnostics = list(fixed_point_error = param_change)
  )
}

# --- Internal helper: convert params list to a single-row data frame --------
#' @keywords internal
.params_to_row <- function(params, loglik, iter) {
  row <- list(
    iter    = iter,
    loglik  = loglik,
    theta0  = params$theta0,
    theta01 = params$theta01,
    theta1  = params$theta1,
    theta10 = params$theta10,
    alpha_00 = if (!is.null(params$alpha)) params$alpha["00"] else NA_real_,
    alpha_11 = if (!is.null(params$alpha)) params$alpha["11"] else NA_real_
  )
  if (!isTRUE(params$asymmetric)) {
    row$pi <- params$pi %||% NA_real_
  } else {
    row$pi0 <- params$pi0 %||% NA_real_
    row$pi1 <- params$pi1 %||% NA_real_
  }
  row
}
