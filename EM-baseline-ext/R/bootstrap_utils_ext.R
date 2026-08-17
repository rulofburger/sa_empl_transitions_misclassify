# ==============================================================================
# EM-baseline-ext: Bootstrap wrappers for extension model families
# Created: 2026-05-06
#
# Provides model-specific bootstrap estimators and the AME summariser for the
# three extension model families (covariates, FMM, inconsistency).
#
# These functions live here — not in EM-baseline/R/bootstrap_utils.R — because
# they depend on em_fit_covariates, em_fit_fmm, em_fit_inconsistency,
# implied_covariates, implied_fmm, and implied_inconsistency, all of which are
# defined in EM-baseline-ext. Keeping them here preserves the one-way
# dependency: EM-baseline is standalone; EM-baseline-ext extends it.
#
# Requires:
#   EM-baseline/R/bootstrap_utils.R  (bootstrap_resample, .LL_THRESHOLD)
#   EM-baseline-ext/R/implied_quantities_ext.R
# ==============================================================================

#' Run one bootstrap replicate for a covariate model
#'
#' Resamples individuals with replacement (applying the same index to both
#' \code{df} and \code{X}), re-estimates the covariate EM model from a
#' warm start, and computes implied probabilities and AMEs.
#'
#' @param df           Original (pre-resample) data frame.
#' @param X            Original N × p covariate matrix.
#' @param seed         Integer seed for this replicate.
#' @param model_type   Passed to \code{em_fit_covariates()}: \code{"symmetric"}
#'   or \code{"none"}.
#' @param stationary   Logical. Passed to \code{em_fit_covariates()}.
#' @param params_start Named parameter list used as warm start.
#' @param point_loglik Scalar point-estimate log-likelihood for quality flagging.
#' @param pi_cap       Passed to \code{em_fit_covariates()}.
#'
#' @return A named list with elements \code{params}, \code{implied},
#'   \code{loglik}, \code{converged}, \code{flag}. See
#'   \code{bootstrap_one_baseline()} for flag values.
bootstrap_one_covariates <- function(df, X, seed, model_type, stationary,
                                     params_start, point_loglik, pi_cap = 0.49) {
  boot <- bootstrap_resample(df, seed, X = X)
  attr(boot$X, "entry_active") <- attr(X, "entry_active")
  fit  <- tryCatch(
    em_fit_covariates(
      df         = boot$df,
      X          = boot$X,
      model_type = model_type,
      stationary = stationary,
      params0    = params_start,
      max_iter   = 500L,
      tol        = 1e-6,
      pi_cap     = pi_cap,
      verbose    = 0L
    ),
    error = function(e) list(.error = conditionMessage(e))
  )

  if (!is.null(fit$.error)) {
    return(list(params = NULL, implied = NULL, loglik = NA_real_,
                converged = FALSE, flag = "error"))
  }

  # For symmetric models, also fit from the nested no-error solution. This
  # supplies a likelihood floor in every bootstrap sample rather than relying
  # on a single interior warm start.
  if (model_type == "symmetric") {
    fit_nested <- tryCatch({
      p_none <- params_start
      p_none$pi <- NULL
      restricted <- em_fit_covariates(
        df = boot$df, X = boot$X, model_type = "none",
        stationary = stationary, params0 = p_none,
        max_iter = 500L, tol = 1e-6, pi_cap = pi_cap, verbose = 0L
      )
      p_nested <- restricted$params
      p_nested$pi <- 1e-8
      candidate <- em_fit_covariates(
        df = boot$df, X = boot$X, model_type = "symmetric",
        stationary = stationary, params0 = p_nested,
        max_iter = 500L, tol = 1e-6, pi_cap = pi_cap, verbose = 0L
      )
      if (candidate$loglik < restricted$loglik - 1e-5)
        stop("symmetric bootstrap fit failed no-error nesting check")
      candidate
    }, error = function(e) NULL)
    if (!is.null(fit_nested) && fit_nested$loglik > fit$loglik)
      fit <- fit_nested
  }

  flag <- .flag_fit(fit, point_loglik)

  implied <- tryCatch(
    implied_covariates(fit$params, boot$X, model_type),
    error = function(e) NULL
  )

  list(params    = fit$params,
       implied   = implied,
       loglik    = fit$loglik,
       converged = fit$converged,
       flag      = flag)
}

#' Run one bootstrap replicate for an FMM model
#'
#' @param df           Original data frame.
#' @param seed         Integer seed.
#' @param model_type   Passed to \code{em_fit_fmm()}: \code{"symmetric"} or
#'   \code{"none"}.
#' @param stationary   Logical. Passed to \code{em_fit_fmm()}.
#' @param params_start Named parameter list for warm start.
#' @param point_loglik Scalar point-estimate log-likelihood.
#' @param theta_cap    Passed to \code{em_fit_fmm()}.
#' @param pi_cap       Passed to \code{em_fit_fmm()}.
#'
#' @return Same structure as \code{bootstrap_one_baseline()}.
bootstrap_one_fmm <- function(df, seed, model_type, stationary,
                               params_start, point_loglik,
                               theta_cap = 0.999, pi_cap = 0.49) {
  boot <- bootstrap_resample(df, seed)
  fit  <- tryCatch(
    em_fit_fmm(
      df         = boot$df,
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

  flag <- .flag_fit(fit, point_loglik)

  implied <- tryCatch(
    implied_fmm(fit$params, model_type),
    error = function(e) NULL
  )

  list(params    = fit$params,
       implied   = implied,
       loglik    = fit$loglik,
       converged = fit$converged,
       flag      = flag)
}

#' Run one bootstrap replicate for an inconsistency model
#'
#' @param df           Original data frame.
#' @param incons_mat   Original N × 6 inconsistency indicator matrix.
#' @param seed         Integer seed.
#' @param stationary   Logical. Passed to \code{em_fit_inconsistency()}.
#' @param params_start Named parameter list for warm start.
#' @param point_loglik Scalar point-estimate log-likelihood.
#' @param theta_cap    Passed to \code{em_fit_inconsistency()}.
#'
#' @return Same structure as \code{bootstrap_one_baseline()}.
bootstrap_one_inconsistency <- function(df, incons_mat, seed, stationary,
                                         params_start, point_loglik,
                                         theta_cap = 0.999) {
  boot <- bootstrap_resample(df, seed, incons_mat = incons_mat)
  fit  <- tryCatch(
    em_fit_inconsistency(
      df         = boot$df,
      incons_mat = boot$incons_mat,
      stationary = stationary,
      params0    = params_start,
      max_iter   = 1000L,
      tol        = 1e-8,
      theta_cap  = theta_cap,
      verbose    = 0L
    ),
    error = function(e) list(.error = conditionMessage(e))
  )

  if (!is.null(fit$.error)) {
    return(list(params = NULL, implied = NULL, loglik = NA_real_,
                converged = FALSE, flag = "error"))
  }

  flag <- .flag_fit(fit, point_loglik)

  implied <- tryCatch(
    implied_inconsistency(fit$params, boot$incons_mat),
    error = function(e) NULL
  )

  list(params    = fit$params,
       implied   = implied,
       loglik    = fit$loglik,
       converged = fit$converged,
       flag      = flag)
}

#' Summarise bootstrap AMEs (named vectors)
#'
#' Like \code{summarise_bootstrap()} but specifically for the \code{ame_entry}
#' and \code{ame_exit} named vectors returned by \code{implied_covariates()}.
#'
#' @param boot_results List of B results from \code{bootstrap_one_covariates()}.
#' @param point_ame_entry Named numeric vector: point-estimate AME on entry rate.
#' @param point_ame_exit  Named numeric vector: point-estimate AME on exit rate.
#'
#' @return A data frame with columns: \code{covariate}, \code{outcome}
#'   (\code{"entry"} or \code{"exit"}), \code{estimate}, \code{se},
#'   \code{ci_lower}, \code{ci_upper}, \code{n_ok}, \code{n_reps}.
#'
#' @examples
#' # Typically called inside .run_and_save() in bootstrap_pipeline.R
#' # after running B reps of bootstrap_one_covariates().
summarise_bootstrap_ame <- function(boot_results, point_ame_entry, point_ame_exit) {
  n_reps <- length(boot_results)
  flags  <- vapply(boot_results, function(r) r$flag, character(1L))
  ok_idx <- which(flags == "ok")
  n_ok   <- length(ok_idx)

  cov_nms <- names(point_ame_entry)

  .extract_ame <- function(type) {
    lapply(ok_idx, function(i) {
      r <- boot_results[[i]]
      if (is.null(r$implied)) return(rep(NA_real_, length(cov_nms)))
      r$implied[[paste0("ame_", type)]]
    })
  }

  .ame_summary <- function(point_ame, type) {
    vecs <- .extract_ame(type)
    K    <- length(cov_nms)
    bmat <- matrix(NA_real_, nrow = n_ok, ncol = K,
                   dimnames = list(NULL, cov_nms))
    for (i in seq_along(vecs)) {
      v <- vecs[[i]]
      if (!is.null(v) && length(v) == K) bmat[i, ] <- v
    }
    ses   <- collapse::fsd(bmat, na.rm = TRUE)  # single vectorised C pass
    ci_q  <- apply(bmat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    ci_lo <- ci_q[1L, ]
    ci_hi <- ci_q[2L, ]
    data.frame(
      covariate = cov_nms,
      outcome   = type,
      estimate  = point_ame[cov_nms],
      se        = ses,
      ci_lower  = ci_lo,
      ci_upper  = ci_hi,
      n_ok      = n_ok,
      n_reps    = n_reps,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  rbind(.ame_summary(point_ame_entry, "entry"),
        .ame_summary(point_ame_exit,  "exit"))
}
