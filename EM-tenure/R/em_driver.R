# ==============================================================================
# EM-tenure: EM driver — main fitting function
# ==============================================================================
# Orchestrates E-step / M-step iteration with convergence monitoring,
# monotonicity checks, and parameter history logging.
#
# TeX reference: Section 2.7 (Iteration and convergence)
# ==============================================================================

#' Initialise parameters for the EM algorithm
#'
#' @param df Data frame with y1, y2, y3 columns (and tenure1-3 / timegap_cat1-3
#'   or timegap1-3 depending on discrete_timegap).
#' @param misclassification Logical; if FALSE, set pi=0.
#' @param discrete_timegap Logical (default TRUE). If TRUE, sigma2_d is not
#'   initialised (discrete model has no sigma2_d). If FALSE, sigma2_d is
#'   initialised from observed timegap increments.
#' @return Named list of starting parameters.
#' @export
init_params <- function(df, misclassification = TRUE, discrete_timegap = TRUE) {
  # Crude moment-based initialisers
  emp_rate <- mean(c(df$y1, df$y2, df$y3), na.rm = TRUE)
  alpha  <- .bound01(emp_rate, eps = 0.01)
  theta1 <- .bound01(0.90, eps = 0.01)
  theta0 <- .bound01(0.10, eps = 0.01)

  pi_init <- if (misclassification) 0.05 else 0

  # Variance from observed increments (rough)
  dg2 <- df$tenure2 - df$tenure1
  dg3 <- df$tenure3 - df$tenure2

  # Use positive increments only (employed continuations)
  emp2 <- (df$y1 == 1) & (df$y2 == 1)
  emp3 <- (df$y2 == 1) & (df$y3 == 1)
  dg_vals <- c(dg2[emp2], dg3[emp3]) - 0.25
  sigma2_g <- if (length(dg_vals[!is.na(dg_vals)]) > 10) {
    max(var(dg_vals, na.rm = TRUE) / 2, 1e-4)
  } else {
    0.01
  }

  lambda_g <- ctmc_lambda_from_theta(theta1)
  lambda_d <- ctmc_lambda_from_theta(theta0)

  out <- list(
    alpha    = alpha,
    theta0   = theta0,
    theta1   = theta1,
    pi       = pi_init,
    sigma2_g = sigma2_g,
    lambda_g = lambda_g,
    lambda_d = lambda_d
  )

  if (!discrete_timegap) {
    # Initialise sigma2_d from observed timegap increments
    dd2 <- df$timegap2 - df$timegap1
    dd3 <- df$timegap3 - df$timegap2
    nonemp2 <- (df$y1 == 0) & (df$y2 == 0)
    nonemp3 <- (df$y2 == 0) & (df$y3 == 0)
    dd_vals <- c(dd2[nonemp2], dd3[nonemp3]) - 0.25
    sigma2_d <- if (length(dd_vals[!is.na(dd_vals)]) > 10) {
      max(var(dd_vals, na.rm = TRUE) / 2, 1e-4)
    } else {
      0.01
    }
    out$sigma2_d <- sigma2_d
  }

  out
}


#' Fit the EM model for employment transitions with tenure durations
#'
#' @param df Data frame with columns: y1-y3, tenure1-tenure3, weight.
#'   When discrete_timegap = TRUE: timegap_cat1-timegap_cat3 (integer 1-7).
#'   When discrete_timegap = FALSE: timegap1-timegap3 (continuous, years).
#' @param params0 Optional named list of starting parameters. If NULL, uses
#'   \code{init_params()}.
#' @param misclassification Logical; if FALSE, constrain pi=0.
#' @param stationary Logical; if TRUE, impose stationarity on alpha.
#' @param discrete_timegap Logical (default TRUE). If TRUE, use interval-censored
#'   discrete Exp(lambda_d) model for nonemployment durations. The timegap_cat1-3
#'   columns must be present and contain integers in 1:7. sigma2_d is not
#'   estimated. If FALSE, use legacy continuous EMG model (timegap1-3 required).
#' @param max_iter Maximum number of EM iterations (default 500).
#' @param tol Convergence tolerance on relative change in log-likelihood
#'   (default 1e-8).
#' @param sigma_floor Minimum variance (default 1e-8).
#' @param theta_cap Maximum transition probability (default 0.999).
#' @param pi_cap Maximum misclassification probability (default 0.49).
#' @param verbose Integer; 0=silent, 1=final summary, 2=every iteration.
#' @return List with:
#'   - params: final parameter estimates
#'   - loglik: final observed-data log-likelihood
#'   - history: data.frame of parameter values at each iteration
#'   - converged: logical
#'   - iterations: number of iterations performed
#'   - gamma: N x 8 matrix of final responsibilities
#' @export
em_fit_tenure <- function(df,
                          params0 = NULL,
                          misclassification = TRUE,
                          stationary = FALSE,
                          discrete_timegap = TRUE,
                          max_iter = 500L,
                          tol = 1e-8,
                          sigma_floor = 1e-8,
                          theta_cap = 0.999,
                          pi_cap = 0.49,
                          verbose = 1L) {
  # Validate required columns (common to both modes)
  required_base <- c("y1", "y2", "y3",
                     "tenure1", "tenure2", "tenure3",
                     "weight")
  missing_base <- setdiff(required_base, names(df))
  if (length(missing_base) > 0) {
    stop("Missing columns: ", paste(missing_base, collapse = ", "))
  }

  if (discrete_timegap) {
    required_cat <- c("timegap_cat1", "timegap_cat2", "timegap_cat3")
    missing_cat <- setdiff(required_cat, names(df))
    if (length(missing_cat) > 0) {
      stop("discrete_timegap=TRUE requires columns: ",
           paste(missing_cat, collapse = ", "))
    }
    if (!all(df$timegap_cat1 %in% 1:7, na.rm = TRUE) ||
        !all(df$timegap_cat2 %in% 1:7, na.rm = TRUE) ||
        !all(df$timegap_cat3 %in% 1:7, na.rm = TRUE)) {
      stop("timegap_cat1/2/3 must be integers in 1:7 (no NA allowed).")
    }
  } else {
    required_cont <- c("timegap1", "timegap2", "timegap3")
    missing_cont <- setdiff(required_cont, names(df))
    if (length(missing_cont) > 0) {
      stop("discrete_timegap=FALSE requires columns: ",
           paste(missing_cont, collapse = ", "))
    }
  }

  total_weight <- sum(df$weight)

  # Initialise parameters
  params <- if (!is.null(params0)) {
    params0
  } else {
    init_params(df, misclassification, discrete_timegap = discrete_timegap)
  }

  # Ensure lambdas are consistent with thetas at start
  params$lambda_g <- ctmc_lambda_from_theta(params$theta1)
  params$lambda_d <- ctmc_lambda_from_theta(params$theta0)

  # History storage
  history <- vector("list", max_iter)
  ll_vec  <- numeric(max_iter)

  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # --- E-step ---
    estep_out <- e_step(df, params, discrete_timegap = discrete_timegap)
    ll_new    <- estep_out$loglik
    ll_vec[iter] <- ll_new

    # --- Monotonicity check ---
    if (iter > 1) {
      ll_change <- ll_new - ll_vec[iter - 1]
      if (ll_change < -1e-8) {
        warning(sprintf(
          "EM iter %d: log-likelihood decreased by %.6e (from %.6e to %.6e)",
          iter, ll_change, ll_vec[iter - 1], ll_new
        ))
      }

      # Convergence check (relative change)
      rel_change <- abs(ll_change) / (abs(ll_vec[iter - 1]) + 1e-16)
      if (rel_change < tol) {
        converged <- TRUE
        if (verbose >= 1) {
          message(sprintf("EM converged at iteration %d (rel_change = %.2e)", iter, rel_change))
        }
        history[[iter]] <- unlist(params)
        break
      }
    }

    # --- M-step ---
    params_prev <- params
    params <- m_step(
      suff              = estep_out$suff,
      total_weight      = total_weight,
      misclassification = misclassification,
      stationary        = stationary,
      discrete_timegap  = discrete_timegap,
      sigma_floor       = sigma_floor,
      theta_cap         = theta_cap,
      pi_cap            = pi_cap
    )

    # --- Monotonicity guard ---
    # When the M-step uses constrained updates (e.g., stationary alpha), the
    # parameter update can violate the EM ascent property. We check by
    # evaluating LL at the new params. If LL decreased, revert to previous
    # params and declare convergence (the algorithm is at a constrained
    # stationary point).
    estep_check <- e_step(df, params, discrete_timegap = discrete_timegap)
    if (estep_check$loglik < ll_new - 1e-8 * abs(ll_new)) {
      if (verbose >= 2) {
        message(sprintf(
          "EM iter %d: M-step decreased LL (%.6e -> %.6e); reverting to previous params.",
          iter, ll_new, estep_check$loglik
        ))
      }
      params <- params_prev
      converged <- TRUE
      history[[iter]] <- unlist(params)
      break
    }

    history[[iter]] <- unlist(params)

    if (verbose >= 2) {
      if (discrete_timegap) {
        message(sprintf(
          "EM iter %3d | ll = %14.4f | alpha=%.4f theta1=%.4f theta0=%.4f pi=%.4f sig2g=%.6f",
          iter, ll_new, params$alpha, params$theta1, params$theta0,
          params$pi, params$sigma2_g
        ))
      } else {
        message(sprintf(
          "EM iter %3d | ll = %14.4f | alpha=%.4f theta1=%.4f theta0=%.4f pi=%.4f sig2g=%.6f sig2d=%.6f",
          iter, ll_new, params$alpha, params$theta1, params$theta0,
          params$pi, params$sigma2_g, params$sigma2_d
        ))
      }
    }
  }

  if (!converged && verbose >= 1) {
    message(sprintf("EM did not converge in %d iterations", max_iter))
  }

  # Trim history to actual iterations
  n_iter <- min(iter, max_iter)
  history <- history[seq_len(n_iter)]
  history_mat <- do.call(rbind, history)
  history_df <- as.data.frame(history_mat)
  history_df$iteration <- seq_len(n_iter)
  history_df$loglik <- as.numeric(ll_vec[seq_len(n_iter)])

  # Reuse last E-step if available; only recompute if M-step ran after it
  # (i.e., the algorithm did not converge on this iteration).
  if (!converged) {
    estep_out <- e_step(df, params, discrete_timegap = discrete_timegap)
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
