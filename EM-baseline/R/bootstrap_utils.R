# ==============================================================================
# EM-baseline: Bootstrap utilities
# Created: 2026-05-06
#
# Provides:
#   bootstrap_resample()            — row-level nonparametric resampler
#   bootstrap_one_baseline()        — single-rep estimator for baseline models
#   summarise_bootstrap()           — collapse B replicate lists to SE table
#
# Extension-model wrappers (bootstrap_one_covariates, bootstrap_one_fmm,
# bootstrap_one_inconsistency, summarise_bootstrap_ame) live in
# EM-baseline-ext/R/bootstrap_utils_ext.R because they depend on EM-baseline-ext
# functions. This file is standalone: it only requires EM-baseline symbols.
#
# Design principles:
#  - Each bootstrap_one_*() warm-starts from the point estimate (params_start).
#    No multi-start: one EM run per replicate is sufficient because the bootstrap
#    DGP is a small perturbation of the original likelihood surface.
#  - Quality flags: reps where EM failed to converge OR where loglik is more
#    than .LL_THRESHOLD nats below the point-estimate loglik are flagged.
#    Flagged reps are stored but excluded from SE computation.
#  - Seeds: each rep receives a deterministic seed vector drawn at the top of
#    the pipeline, ensuring reproducibility regardless of parallelism order.
# ==============================================================================

# Minimum acceptable log-likelihood drop: reps with loglik < point_loglik - 50
# nats are treated as numerical failures and excluded from SE computation.
# 50 nats is a conservative threshold chosen to flag severe optimizer failures
# (e.g., convergence to a degenerate mode) while tolerating normal Monte Carlo
# variation across bootstrap samples.
.LL_THRESHOLD <- 50

#' Nonparametric bootstrap resample of a data frame
#'
#' Resamples rows of \code{df} with replacement (individual-level resampling,
#' preserving all columns including survey weights and all three waves). The
#' same row index is applied to \code{X} and/or \code{incons_mat} when provided,
#' ensuring that auxiliary matrices remain aligned with the resampled rows.
#'
#' @param df        Data frame of individual observations to resample.
#' @param seed      Integer seed for reproducibility.
#' @param X         Optional N × p numeric covariate matrix. Resampled using the
#'   same row index as \code{df}. Pass \code{NULL} to skip.
#' @param incons_mat Optional N × 6 inconsistency indicator matrix. Resampled
#'   using the same row index. Pass \code{NULL} to skip.
#'
#' @return A named list:
#'   \describe{
#'     \item{df}{Resampled data frame.}
#'     \item{X}{Resampled covariate matrix, or \code{NULL}.}
#'     \item{incons_mat}{Resampled inconsistency matrix, or \code{NULL}.}
#'   }
bootstrap_resample <- function(df, seed, X = NULL, incons_mat = NULL) {
  stopifnot(is.data.frame(df), is.numeric(seed), length(seed) == 1L)
  required_cols <- c("y1", "y2", "y3", "weight")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L)
    stop(sprintf("bootstrap_resample: df missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  set.seed(seed)
  idx <- sample.int(nrow(df), replace = TRUE)
  list(
    df         = df[idx, , drop = FALSE],
    X          = if (!is.null(X)) X[idx, , drop = FALSE] else NULL,
    incons_mat = if (!is.null(incons_mat)) incons_mat[idx, , drop = FALSE] else NULL
  )
}

#' Run one bootstrap replicate for the baseline model
#'
#' Resamples the data and re-estimates the baseline EM model, then computes
#' implied probabilities. Returns a compact list suitable for collecting into
#' a B-element list by \code{mclapply}.
#'
#' @param df           Data frame for the bootstrap replicate (already resampled
#'   OR original — this function calls \code{bootstrap_resample()} internally).
#' @param seed         Integer seed for this replicate.
#' @param model_type   Passed to \code{em_fit_baseline()}.
#' @param stationary   Passed to \code{em_fit_baseline()}.
#' @param params_start Named list of point-estimate parameters used as warm start.
#' @param point_loglik Scalar: log-likelihood of the point estimate. Used to
#'   flag numerically failed reps.
#' @param theta_cap    Passed to \code{em_fit_baseline()}.
#' @param pi_cap       Passed to \code{em_fit_baseline()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{params}{Named parameter list, or \code{NULL} on failure.}
#'     \item{implied}{Implied quantities list, or \code{NULL} on failure.}
#'     \item{loglik}{Scalar log-likelihood.}
#'     \item{converged}{Logical.}
#'     \item{flag}{Character: \code{"ok"}, \code{"no_converge"}, \code{"low_loglik"},
#'       or \code{"error"}.}
#'   }
bootstrap_one_baseline <- function(df, seed, model_type, stationary,
                                   params_start, point_loglik,
                                   theta_cap = 0.999, pi_cap = 0.49) {
  boot  <- bootstrap_resample(df, seed)
  .bootstrap_run_baseline(boot$df, model_type, stationary,
                          params_start, point_loglik, theta_cap, pi_cap)
}

#' Assess quality of a bootstrap EM fit and assign a flag
#'
#' Compares the replicate log-likelihood against the point estimate.  A rep is
#' flagged \code{"no_converge"} if the EM did not converge, \code{"low_loglik"}
#' if convergence was reported but the log-likelihood dropped more than
#' \code{.LL_THRESHOLD} nats below the point estimate (indicating the optimizer
#' reached a poor mode on an outlier bootstrap sample), and \code{"ok"}
#' otherwise.
#'
#' @param fit        Named list with at least \code{$converged} (logical) and
#'   \code{$loglik} (scalar numeric).
#' @param point_loglik Scalar numeric point-estimate log-likelihood, or
#'   \code{NA} (threshold check is skipped when \code{NA}).
#' @return Character scalar: \code{"ok"}, \code{"no_converge"},
#'   \code{"low_loglik"}, or \code{"error"} (loglik is NA despite convergence).
#' @keywords internal
.flag_fit <- function(fit, point_loglik) {
  if (!fit$converged)    return("no_converge")
  if (is.na(fit$loglik)) return("error")  # optimizer returned NA loglik
  if (!is.na(point_loglik) && fit$loglik < point_loglik - .LL_THRESHOLD)
    return("low_loglik")
  "ok"
}

#' Run baseline EM on already-resampled data and compute implied quantities
#'
#' Internal workhorse called by \code{bootstrap_one_baseline()}.
#'
#' @param df_boot      Resampled data frame (output of \code{bootstrap_resample()}).
#' @param model_type   \code{"symmetric"}, \code{"asymmetric"}, or \code{"none"}.
#' @param stationary   Logical: enforce stationarity constraint.
#' @param params_start Named list of warm-start parameters.
#' @param point_loglik Scalar point-estimate log-likelihood for quality flagging.
#' @param theta_cap    Upper bound on transition probabilities.
#' @param pi_cap       Upper bound on misclassification probability.
#' @return Named list: \code{params}, \code{implied}, \code{loglik},
#'   \code{converged}, \code{flag}.
#' @keywords internal
.bootstrap_run_baseline <- function(df_boot, model_type, stationary,
                                    params_start, point_loglik, theta_cap, pi_cap) {
  fit <- tryCatch(
    em_fit_baseline(
      df         = df_boot,
      model_type = model_type,
      stationary = stationary,
      params0    = params_start,
      max_iter   = 1000L,
      tol        = 1e-8,
      theta_cap  = theta_cap,
      pi_cap     = pi_cap,
      verbose    = 0L
    ),
    error = function(e) list(.error = conditionMessage(e))
  )

  if (!is.null(fit$.error)) {
    return(list(params = NULL, implied = NULL, loglik = NA_real_,
                converged = FALSE, flag = "error"))
  }

  flag    <- .flag_fit(fit, point_loglik)
  implied <- tryCatch(
    implied_baseline(fit$params, model_type),
    error = function(e) NULL
  )

  list(params    = fit$params,
       implied   = implied,
       loglik    = fit$loglik,
       converged = fit$converged,
       flag      = flag)
}

#' Summarise B bootstrap replicate results into a SE table
#'
#' Accepts the list of B results returned by \code{mclapply} (each element is
#' a named list with \code{params}, \code{implied}, \code{flag}).  Extracts
#' scalar quantities from \code{params} and \code{implied}, computes the
#' bootstrap SD across the \code{"ok"} reps (i.e., converged and not flagged),
#' and returns a tidy summary data frame.
#'
#' Only scalar quantities are summarised.  Non-scalar elements (e.g., the
#' \code{ame_entry} and \code{ame_exit} named vectors) are handled separately
#' by \code{summarise_bootstrap_ame()} in \code{EM-baseline-ext}.
#'
#' When the same quantity name appears in both \code{params} and \code{implied}
#' (e.g., \code{pi} in the symmetric baseline model), the \code{params}-side
#' value is kept and the duplicate is dropped to avoid double rows in the
#' output table.
#'
#' @param boot_results List of B named lists, each from a \code{bootstrap_one_*()}
#'   call.
#' @param point_params Named list of point-estimate parameters.
#' @param point_implied Named list of point-estimate implied quantities.
#' @param point_loglik  Scalar point-estimate log-likelihood.
#'
#' @return A data frame with columns: \code{quantity}, \code{estimate},
#'   \code{se}, \code{ci_lower}, \code{ci_upper}, \code{n_ok}, \code{n_reps}.
#'
#' @examples
#' pp     <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
#' pi_pt  <- list(entry_rate = 0.10, exit_rate = 0.10, employment_rate = 0.5,
#'                pi = 0.05, pi0 = NA_real_, pi1 = NA_real_)
#' boots  <- replicate(5, list(params = pp, implied = pi_pt,
#'                             loglik = -100, converged = TRUE, flag = "ok"),
#'                     simplify = FALSE)
#' summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
summarise_bootstrap <- function(boot_results, point_params, point_implied,
                                 point_loglik) {
  n_reps <- length(boot_results)
  flags  <- vapply(boot_results, function(r) r$flag, character(1L))
  ok_idx <- which(flags == "ok")
  n_ok   <- length(ok_idx)

  # Flatten a named list to a scalar numeric vector, dropping non-scalars.
  .flatten_scalar <- function(lst) {
    if (is.null(lst)) return(NULL)
    scalars <- Filter(function(x) is.numeric(x) && length(x) == 1L, lst)
    unlist(scalars)
  }

  boot_mat <- lapply(ok_idx, function(i) {
    r   <- boot_results[[i]]
    raw <- c(.flatten_scalar(r$params), .flatten_scalar(r$implied))
    # Deduplicate: keep params-side value when a name appears in both.
    raw[!duplicated(names(raw))]
  })

  # Reference names from point estimates (deduped the same way).
  point_raw  <- c(.flatten_scalar(point_params), .flatten_scalar(point_implied))
  point_flat <- point_raw[!duplicated(names(point_raw))]
  all_nms    <- names(point_flat)

  if (length(boot_mat) == 0L || length(all_nms) == 0L) {
    return(data.frame(
      quantity  = character(0),
      estimate  = numeric(0),
      se        = numeric(0),
      ci_lower  = numeric(0),
      ci_upper  = numeric(0),
      n_ok      = integer(0),
      n_reps    = integer(0)
    ))
  }

  # Build B × K matrix (only ok reps)
  K    <- length(all_nms)
  bmat <- matrix(NA_real_, nrow = length(boot_mat), ncol = K,
                 dimnames = list(NULL, all_nms))
  for (i in seq_along(boot_mat)) {
    row_vals <- boot_mat[[i]]
    shared   <- intersect(all_nms, names(row_vals))
    bmat[i, shared] <- row_vals[shared]
  }

  ses      <- collapse::fsd(bmat, na.rm = TRUE)  # single vectorised C pass
  ci_q     <- apply(bmat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  ci_lower <- ci_q[1L, ]
  ci_upper <- ci_q[2L, ]

  data.frame(
    quantity  = all_nms,
    estimate  = point_flat[all_nms],
    se        = ses,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    n_ok      = n_ok,
    n_reps    = n_reps,
    row.names = NULL
  )
}

