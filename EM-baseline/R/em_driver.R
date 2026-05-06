# ==============================================================================
# EM-baseline: EM driver — initialisation, iteration loop, and convergence
# Created: 2026-05-05
#
# Top-level entry point for fitting the baseline EM model. Handles parameter
# initialisation, the E-M iteration loop, convergence checking, and result
# packaging.
#
# TeX reference: EM baseline.tex, Section 2.4 (Iteration and convergence).
# ==============================================================================

#' Initialise parameters for the baseline EM model
#'
#' Returns sensible starting values for the EM algorithm. Starting values are
#' chosen to be in the interior of the parameter space and away from
#' degenerate boundaries.
#'
#' @param model_type Character: \code{"symmetric"} (default),
#'   \code{"asymmetric"}, or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}). If \code{TRUE}, alpha is
#'   derived from the stationarity restriction; if \code{FALSE}, alpha is set
#'   as an additional free parameter.
#' @return Named list of initial parameters: \code{alpha}, \code{theta0},
#'   \code{theta1}, and (if applicable) \code{pi}, \code{pi0}, \code{pi1}.
#' @examples
#' init_params("symmetric")
#' init_params("asymmetric", stationary = FALSE)
#' @export
init_params <- function(model_type = "symmetric", stationary = TRUE) {
  theta0 <- 0.10
  theta1 <- 0.90

  alpha <- if (stationary) {
    stationary_alpha(theta0, theta1)
  } else {
    theta0 / (theta0 + 1 - theta1)  # use stationary as warm-start for free alpha
  }

  params <- list(alpha = alpha, theta0 = theta0, theta1 = theta1)

  if (model_type == "symmetric") {
    params$pi <- 0.05
  } else if (model_type == "asymmetric") {
    params$pi0 <- 0.05
    params$pi1 <- 0.05
  }
  # model_type == "none": no misclassification params

  params
}

#' Evaluate the observed-data log-likelihood at given parameters
#'
#' Convenience function that calls \code{\link{e_step}} and returns only the
#' log-likelihood. Useful for LR tests and comparison across models.
#'
#' @param df Data frame with columns y1, y2, y3, weight.
#' @param params Named list of parameters (alpha, theta0, theta1, pi/pi0/pi1).
#' @param model_type Character: \code{"symmetric"}, \code{"asymmetric"}, or
#'   \code{"none"}.
#' @return Scalar: observed-data weighted log-likelihood.
#' @examples
#' df <- data.frame(
#'   y1 = c(1L, 0L), y2 = c(1L, 0L), y3 = c(1L, 0L), weight = c(1, 1)
#' )
#' params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
#' compute_observed_loglik(df, params, model_type = "none")
#' @export
compute_observed_loglik <- function(df, params, model_type = "symmetric") {
  e_step(df, params, model_type = model_type)$loglik
}

#' Fit the baseline EM model for employment transitions with misclassification
#'
#' Implements the EM algorithm described in \emph{EM baseline.tex}. Iterates
#' E-step (responsibilities) and M-step (closed-form parameter updates) until
#' the relative change in log-likelihood is below \code{tol}.
#'
#' TeX ref: EM baseline.tex, Section 2.4 (Eq 20 for convergence criterion).
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer, 0/1) and \code{weight} (positive numeric).
#' @param model_type Character: \code{"symmetric"} (default),
#'   \code{"asymmetric"}, or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}). If \code{TRUE}, alpha is
#'   derived from the stationarity restriction at every M-step.
#' @param params0 Optional named list of starting parameters. If \code{NULL},
#'   \code{\link{init_params}} is used.
#' @param max_iter Maximum number of EM iterations (default 1000L).
#' @param tol Convergence tolerance on relative change in log-likelihood
#'   (default 1e-8).
#' @param theta_cap Upper bound for theta0 and theta1 (default 0.999).
#' @param pi_cap Upper bound for pi, pi0, pi1 (default 0.49).
#' @param verbose Integer verbosity level. 0 = silent, 1 = convergence
#'   summary only, 2 = per-iteration output (default 1L).
#' @return List with:
#'   \describe{
#'     \item{params}{Named list of final parameter estimates.}
#'     \item{loglik}{Scalar: final observed-data log-likelihood.}
#'     \item{history}{Data frame with one row per iteration: iter, loglik,
#'       and one column per parameter.}
#'     \item{converged}{Logical: whether convergence was achieved.}
#'     \item{iterations}{Integer: number of completed iterations.}
#'     \item{gamma}{N x 8 matrix of final responsibilities.}
#'     \item{model_type}{Character: the model_type argument.}
#'     \item{stationary}{Logical: the stationary argument.}
#'   }
#' @details
#' \strong{Validation strategy}: data quality (factor check, NA, binary,
#' positive weights) is validated \emph{once} before entering the EM loop.
#' \code{\link{e_step}} is then called with \code{validate = FALSE} on
#' every iteration and on the final E-step, because the input data frame is
#' immutable between initialization and convergence. This avoids redundant
#' per-iteration scans without sacrificing safety.
#'
#' TeX ref: EM baseline.tex, Section 2.4 (Eq 20 for convergence criterion).
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y1 = sample(0:1, 200, TRUE), y2 = sample(0:1, 200, TRUE),
#'   y3 = sample(0:1, 200, TRUE), weight = rep(1, 200)
#' )
#' fit <- em_fit_baseline(df, model_type = "symmetric", verbose = 0L)
#' fit$loglik
#' }
#' @export
em_fit_baseline <- function(df,
                            model_type = "symmetric",
                            stationary = TRUE,
                            params0    = NULL,
                            max_iter   = 1000L,
                            tol        = 1e-8,
                            theta_cap  = 0.999,
                            pi_cap     = 0.49,
                            verbose    = 1L) {
  # --- Validate inputs -------------------------------------------------------
  .validate_panel_df(df)
  .validate_model_type(model_type)

  # Validate params0 before assigning (fail-fast: validate first, then assign)
  if (!is.null(params0)) {
    required_params <- c("alpha", "theta0", "theta1")
    if (model_type == "symmetric")  required_params <- c(required_params, "pi")
    if (model_type == "asymmetric") required_params <- c(required_params, "pi0", "pi1")
    missing_params <- setdiff(required_params, names(params0))
    if (length(missing_params) > 0)
      stop(sprintf("params0 is missing fields required for model_type='%s': %s",
                   model_type, paste(missing_params, collapse = ", ")))
    # Validate parameter values: must be finite numeric scalars in (0, 1)
    for (pname in required_params) {
      pval <- params0[[pname]]
      if (!is.numeric(pval) || length(pval) != 1L || !is.finite(pval) || pval <= 0 || pval >= 1)
        stop(sprintf(
          "params0$%s must be a single finite numeric in (0, 1); got: %s",
          pname, if (is.numeric(pval) && length(pval) == 1L && !is.na(pval)) format(pval) else class(pval)
        ))
    }
  }
  params <- params0 %||% init_params(model_type = model_type, stationary = stationary)

  # --- Validate data quality once before the loop ----------------------------
  # (e_step is called with validate = FALSE in the loop to avoid repeating on
  # immutable data at every iteration)
  for (.col in c("y1", "y2", "y3")) {
    if (is.factor(df[[.col]]))
      stop(sprintf(
        "Column '%s' is a factor; coerce to integer 0/1 before calling em_fit_baseline()",
        .col
      ))
  }
  if (anyNA(df$y1) || anyNA(df$y2) || anyNA(df$y3))
    stop("em_fit_baseline: y1/y2/y3 must not contain NA values")
  if (anyNA(df$weight) || any(df$weight <= 0))
    stop("em_fit_baseline: weight must be non-NA and strictly positive")
  for (.col in c("y1", "y2", "y3")) {
    vals <- df[[.col]]
    if (!all(vals %in% c(0, 1)))
      stop(sprintf("em_fit_baseline: column '%s' must be binary (0/1 only)", .col))
  }

  # --- EM iteration loop -----------------------------------------------------
  ll_prev   <- -Inf
  converged <- FALSE

  # Pre-allocate history matrix: O(1) per iteration (vs O(n²) for list-rbind).
  # Capture param_names upfront from initial params to ensure consistent column
  # order across iterations (m_step returns fields in a fixed order, but
  # selecting by name is robust to any future reordering).
  param_names <- names(params)
  col_names   <- c("iter", "loglik", param_names)
  history_mat <- matrix(NA_real_, nrow = max_iter, ncol = length(col_names),
                        dimnames = list(NULL, col_names))

  for (iter in seq_len(max_iter)) {

    # E-step (validate = FALSE: data quality checked once before the loop)
    estep_out <- e_step(df, params, model_type = model_type, validate = FALSE)
    ll_new    <- estep_out$loglik

    # Record history (O(1) assignment)
    history_mat[iter, ] <- c(iter, ll_new, unlist(params[param_names]))

    # Verbose output
    if (verbose >= 2L) {
      cat(sprintf("Iter %4d | loglik = %.6f | rel_change = %.2e\n",
                  iter, ll_new,
                  abs(ll_new - ll_prev) / (abs(ll_prev) + 1e-16)))
    }

    # Convergence check (after first E-step gives a valid ll)
    if (iter > 1L) {
      rel_change <- abs(ll_new - ll_prev) / (abs(ll_prev) + 1e-16)
      if (rel_change < tol) {
        converged <- TRUE
        if (verbose >= 1L) {
          cat(sprintf(
            "EM [%s, stationary=%s] converged in %d iterations. loglik = %.4f\n",
            model_type, stationary, iter, ll_new
          ))
        }
        break
      }
      # Monotonicity guard: warn if LL decreases (should not happen)
      if (ll_new < ll_prev - 1e-6) {
        warning(sprintf(
          "EM iteration %d: log-likelihood decreased by %.2e (numerical noise).",
          iter, ll_prev - ll_new
        ))
      }
    }

    ll_prev <- ll_new

    # M-step
    m_out  <- m_step(estep_out$suff,
                     model_type = model_type,
                     stationary = stationary,
                     theta_cap  = theta_cap,
                     pi_cap     = pi_cap)
    # Guard: catch any m_step contract violation (wrong/missing/extra fields)
    stopifnot(identical(sort(param_names), sort(names(m_out))))
    params <- m_out
  }

  if (!converged && verbose >= 1L) {
    warning(sprintf(
      "EM [%s, stationary=%s] did not converge in %d iterations. loglik = %.4f",
      model_type, stationary, max_iter, ll_new
    ))
  }

  # --- Final E-step for gamma ------------------------------------------------
  final_estep <- e_step(df, params, model_type = model_type, validate = FALSE)

  # --- Build history data frame (O(1) conversion from pre-allocated matrix) --
  history_df <- as.data.frame(history_mat[seq_len(iter), , drop = FALSE])

  list(
    params     = params,
    loglik     = final_estep$loglik,
    history    = history_df,
    converged  = converged,
    iterations = iter,
    gamma      = final_estep$gamma,
    model_type = model_type,
    stationary = stationary
  )
}
