# ==============================================================================
# EM-tenure: Bootstrap utilities for the contamination (epsilon) model
# Created: 2026-05-07
#
# Provides:
#   bootstrap_one_eps()  — single bootstrap replicate runner for the eps model
#
# Design notes:
#   - Reuses bootstrap_resample() and summarise_bootstrap() from
#     EM-baseline/R/bootstrap_utils.R (must be sourced before this file).
#   - Reuses implied_tenure_contamination() from
#     EM-tenure/R/implied_quantities_tenure_contamination.R.
#   - Warm-starts from the point-estimate parameter list (no multi-start).
#   - Quality flags: same .LL_THRESHOLD logic as baseline (50 nats).
#   - The resampler preserves all required eps-model columns:
#     y1–y3, tenure1–tenure3, timegap_cat1–timegap_cat3, weight.
#
# Companion: documents/EM tenure epsilon.tex
# Used by: EM-tenure/bootstrap_pipeline_tenure_contamination.R
# ==============================================================================

#' Run one bootstrap replicate for the eps contamination model
#'
#' Resamples the data and re-estimates the eps EM model, then computes
#' implied tenure-contamination quantities. Returns a compact list suitable
#' for collecting into a B-element list by \code{mclapply}.
#'
#' Internally calls \code{bootstrap_resample()} to resample rows, then
#' \code{em_fit_tenure_eps()} as the estimator (warm-started from
#' \code{params_start}), and finally \code{implied_tenure_contamination()}
#' to convert the fitted params to interpretable quantities.
#'
#' The function requires the following functions to already be loaded in the
#' session (via the respective source_all.R files):
#' \code{bootstrap_resample()}, \code{em_fit_tenure_eps()},
#' \code{implied_tenure_contamination()}.
#'
#' @param df           Data frame with columns: y1–y3, tenure1–tenure3,
#'   timegap_cat1–timegap_cat3, weight.
#' @param seed         Integer seed for this replicate (determines the
#'   resampled row indices and is passed to \code{bootstrap_resample()}).
#' @param stationary   Logical; passed to \code{em_fit_tenure_eps()}.
#' @param linked       Logical; passed to \code{em_fit_tenure_eps()}.
#' @param params_start Named list of point-estimate parameters used as
#'   warm start (must include \code{eps}).
#' @param point_loglik Scalar: log-likelihood of the point estimate; used to
#'   flag numerically failed reps via the .LL_THRESHOLD criterion.
#' @param theta_cap    Upper bound on theta1 / theta0. Default 0.999.
#' @param pi_cap       Upper bound on pi. Default 0.49.
#' @param eps_cap      Upper bound on eps. Default 0.95.
#' @param eps_floor    Lower bound on eps. Default 1e-4.
#' @param max_iter     Maximum EM iterations. Default 500L (matches
#'   \code{em_fit_tenure_eps()} default and \code{estimate_pipeline.R}).
#'
#' @return A named list:
#'   \describe{
#'     \item{params}{Named parameter list, or \code{NULL} on failure.}
#'     \item{implied}{Implied quantities list from
#'       \code{implied_tenure_contamination()}, or \code{NULL} on failure.}
#'     \item{loglik}{Scalar log-likelihood of the replicate fit.}
#'     \item{converged}{Logical: did the EM converge?}
#'     \item{flag}{Character: \code{"ok"}, \code{"no_converge"},
#'       \code{"low_loglik"}, or \code{"error"}.}
#'   }
#'
#' @seealso \code{bootstrap_resample()}, \code{em_fit_tenure_eps()},
#'   \code{implied_tenure_contamination()}, \code{summarise_bootstrap()}
#' @export
bootstrap_one_eps <- function(df, seed,
                              stationary, linked,
                              params_start, point_loglik,
                              theta_cap  = 0.999,
                              pi_cap     = 0.49,
                              eps_cap    = 0.95,
                              eps_floor  = 1e-4,
                              max_iter   = 500L) {
  # Resample rows (preserves all columns including tenure + timegap)
  boot <- bootstrap_resample(df, seed)

  # Run one EM replicate with warm-start
  fit <- tryCatch(
    em_fit_tenure_eps(
      df         = boot$df,
      params0    = params_start,
      stationary = stationary,
      linked     = linked,
      max_iter   = max_iter,
      tol        = 1e-8,
      theta_cap  = theta_cap,
      pi_cap     = pi_cap,
      eps_cap    = eps_cap,
      eps_floor  = eps_floor,
      verbose    = 0L
    ),
    error = function(e) list(.error = conditionMessage(e))
  )

  # Handle hard errors (e.g., data validation failure on a degenerate resample)
  if (!is.null(fit$.error)) {
    return(list(
      params    = NULL,
      implied   = NULL,
      loglik    = NA_real_,
      converged = FALSE,
      flag      = "error"
    ))
  }

  # Assign quality flag using the same criterion as the baseline bootstrap
  flag <- .flag_fit_eps(fit, point_loglik)

  # Compute implied quantities (safe: tryCatch in case params are degenerate)
  implied <- tryCatch(
    implied_tenure_contamination(fit$params),
    error = function(e) NULL
  )

  list(
    params    = fit$params,
    implied   = implied,
    loglik    = fit$loglik,
    converged = fit$converged,
    flag      = flag
  )
}

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Threshold for flagging bootstrap reps with anomalously low log-likelihood.
# Mirrors .LL_THRESHOLD in EM-baseline/R/bootstrap_utils.R (50 nats).
# Defined here independently so this file does not depend on the baseline
# constant being defined in the namespace at load time.
.LL_THRESHOLD_EPS <- 50

#' Assess quality of a bootstrap eps fit and assign a flag
#'
#' @param fit        Named list with at least \code{$converged} (logical) and
#'   \code{$loglik} (scalar numeric).
#' @param point_loglik Scalar numeric point-estimate log-likelihood.
#' @return Character scalar: \code{"ok"}, \code{"no_converge"},
#'   \code{"low_loglik"}, or \code{"error"}.
#' @keywords internal
.flag_fit_eps <- function(fit, point_loglik) {
  if (!isTRUE(fit$converged))    return("no_converge")
  if (is.na(fit$loglik))         return("error")
  if (!is.na(point_loglik) &&
      fit$loglik < point_loglik - .LL_THRESHOLD_EPS)
    return("low_loglik")
  "ok"
}
