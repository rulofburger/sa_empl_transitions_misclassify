# ==============================================================================
# EM-AR2: Bootstrap utilities
# Created: 2026-05-07
#
# Provides:
#   .find_latest_fit()       — find the most recent timestamped .rds fit file
#   bootstrap_one_ar2()      — single-rep estimator for AR(2) models
#
# Reuses from EM-baseline/R/bootstrap_utils.R:
#   bootstrap_resample()     — row-level nonparametric resampler
#   summarise_bootstrap()    — collapse B replicate lists to SE table
#
# The AR(2) data has columns y1, y2, y3, y4, weight.
# bootstrap_resample() validates y1, y2, y3, weight — the y4 column is
# preserved in the resampled data frame because row-index resampling copies
# all columns. This dependency is documented here in case bootstrap_resample()
# ever adds a y4 requirement.
#
# Quality flag thresholds mirror the baseline (.LL_THRESHOLD = 50 nats).
# ==============================================================================

# Minimum acceptable LL drop before flagging a rep as numerically failed.
# Mirrors .LL_THRESHOLD from EM-baseline/R/bootstrap_utils.R.
.LL_THRESHOLD_AR2 <- 50

#' Find the most recent timestamped AR(2) fit file for a model type
#'
#' Globs files matching \code{em_ar2_{model_key}_*.rds} in \code{results_dir},
#' sorts by timestamp suffix (YYYYMMDD_HHMMSS) in descending order, and
#' returns the path to the latest file.
#'
#' @param model_key Character scalar: one of \code{"nopi"}, \code{"sym"},
#'   \code{"asym"}.
#' @param results_dir Character scalar: path to the directory containing
#'   the fit \code{.rds} files.
#'
#' @return Character scalar: absolute path to the latest \code{.rds} file.
#'   Stops with an informative error if no matching file is found.
#' @examples
#' \dontrun{
#' path <- .find_latest_fit("sym", "EM-AR2/output/results")
#' fit  <- readRDS(path)
#' }
#' @keywords internal
.find_latest_fit <- function(model_key, results_dir) {
  stopifnot(
    is.character(model_key),   length(model_key)   == 1L,
    is.character(results_dir), length(results_dir) == 1L
  )
  pattern <- sprintf("em_ar2_%s_*.rds", model_key)
  matches <- Sys.glob(file.path(results_dir, pattern))
  if (length(matches) == 0L) {
    stop(sprintf(
      ".find_latest_fit: no files matching '%s' in '%s'.\n",
      pattern, results_dir
    ), "Run estimate_pipeline.R first to produce point estimates.")
  }
  # Timestamped names are the audit trail; modification times may change when
  # results are copied or restored and therefore cannot identify the latest run.
  bases <- basename(matches)
  stamp <- sub(sprintf("^em_ar2_%s_([0-9]{8}_[0-9]{6}).*$", model_key),
               "\\1", bases)
  ordered <- matches[order(stamp, bases, decreasing = TRUE)]
  ordered[[1L]]
}

#' Assess quality of an AR(2) bootstrap EM fit and assign a flag
#'
#' @param fit Named list with at least \code{$converged} and \code{$loglik}.
#' @param point_loglik Scalar point-estimate log-likelihood (or \code{NA} to
#'   skip LL threshold check).
#' @return Character scalar: \code{"ok"}, \code{"no_converge"},
#'   \code{"low_loglik"}, or \code{"error"}.
#' @keywords internal
.flag_fit_ar2 <- function(fit, point_loglik) {
  if (!isTRUE(fit$converged))  return("no_converge")
  if (is.na(fit$loglik))       return("error")
  if (!is.na(point_loglik) && fit$loglik < point_loglik - .LL_THRESHOLD_AR2)
    return("low_loglik")
  "ok"
}

#' Map AR(2) model type string to em_fit_ar2() arguments
#'
#' @param model_type Character scalar: \code{"none"}, \code{"symmetric"},
#'   or \code{"asymmetric"}.
#' @return Named list of arguments to pass to \code{em_fit_ar2()}.
#' @keywords internal
.ar2_fit_args <- function(model_type) {
  switch(model_type,
    none       = list(estimate_pi = FALSE, fixed_pi = 0,    asymmetric = FALSE),
    symmetric  = list(estimate_pi = TRUE,  fixed_pi = NULL, asymmetric = FALSE),
    asymmetric = list(estimate_pi = TRUE,  fixed_pi = NULL, asymmetric = TRUE),
    stop(sprintf(".ar2_fit_args: unknown model_type '%s'.", model_type))
  )
}

#' Run one bootstrap replicate for the AR(2) model
#'
#' Resamples the data (using \code{bootstrap_resample()} from EM-baseline),
#' re-estimates the AR(2) EM model warm-started from the point estimate, and
#' computes implied quantities. Returns a compact list suitable for collecting
#' into a B-element list by \code{parallel::mclapply}.
#'
#' @param df           Data frame with columns y1, y2, y3, y4, weight.
#' @param seed         Integer seed for this replicate (reproducibility).
#' @param model_type   Character scalar: \code{"none"}, \code{"symmetric"},
#'   or \code{"asymmetric"}.
#' @param params_start Named list of point-estimate parameters (warm start).
#' @param point_loglik Scalar: point-estimate log-likelihood for quality
#'   flagging. Pass \code{NA} to skip LL threshold check.
#' @param pi_cap       Maximum misclassification probability allowed in
#'   \code{em_fit_ar2()}. Default \code{0.49}.
#'
#' @return A named list:
#'   \describe{
#'     \item{params}{Named parameter list, or \code{NULL} on failure.}
#'     \item{implied}{Implied quantities list from \code{implied_ar2()},
#'       or \code{NULL} on failure.}
#'     \item{loglik}{Scalar log-likelihood.}
#'     \item{converged}{Logical.}
#'     \item{flag}{Character: \code{"ok"}, \code{"no_converge"},
#'       \code{"low_loglik"}, or \code{"error"}.}
#'   }
#' @export
bootstrap_one_ar2 <- function(df, seed, model_type, params_start,
                               point_loglik, pi_cap = 0.49) {
  # Guard: y4 is required by em_fit_ar2() but not checked by bootstrap_resample()
  if (!"y4" %in% names(df))
    stop("bootstrap_one_ar2: df missing required column 'y4'.")
  boot <- bootstrap_resample(df, seed)  # from EM-baseline/R/bootstrap_utils.R

  fit_args <- .ar2_fit_args(model_type)

  fit <- tryCatch(
    do.call(em_fit_ar2, c(
      list(
        df       = boot$df,
        params0  = params_start,
        max_iter = 1000L,
        tol      = 1e-8,
        pi_cap   = pi_cap,
        verbose  = 0L
      ),
      fit_args
    )),
    error = function(e) list(.error = conditionMessage(e))
  )

  if (!is.null(fit$.error)) {
    return(list(params   = NULL,
                implied  = NULL,
                loglik   = NA_real_,
                converged = FALSE,
                flag     = "error"))
  }

  flag    <- .flag_fit_ar2(fit, point_loglik)
  implied <- tryCatch(
    implied_ar2(fit$params, model_type),
    error = function(e) NULL
  )
  # Downgrade flag if implied quantities could not be computed: a rep flagged
  # "ok" with implied=NULL would inflate n_ok and corrupt summarise_bootstrap().
  if (is.null(implied) && flag == "ok") flag <- "error"

  # Strip pi fields from params to avoid duplicate quantity names when
  # summarise_bootstrap() concatenates params + implied (both contain pi/pi0/pi1).
  params_out <- fit$params
  params_out$pi  <- NULL
  params_out$pi0 <- NULL
  params_out$pi1 <- NULL

  list(
    params    = params_out,
    implied   = implied,
    loglik    = fit$loglik,
    converged = fit$converged,
    flag      = flag
  )
}
