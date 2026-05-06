# ==============================================================================
# EM-baseline-ext: EM driver for Extension III (2-type FMM)
# Created: 2026-05-06
#
# TeX ref: EM baseline.tex Section 7.
# ==============================================================================

#' Initialise parameters for the 2-type FMM extension
#'
#' Returns starting values with type A having higher employment persistence
#' and type B having lower persistence, separated enough to avoid immediate
#' label switching.
#'
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}).
#' @return Named list with all FMM parameters.
#' @examples
#' init_params_fmm("symmetric")
#' @export
init_params_fmm <- function(model_type = "symmetric", stationary = TRUE) {
  if (!model_type %in% c("symmetric", "none"))
    stop("init_params_fmm: model_type must be 'symmetric' or 'none'")

  theta0_A <- 0.12; theta1_A <- 0.92
  theta0_B <- 0.08; theta1_B <- 0.80

  alpha_A <- if (stationary) stationary_alpha(theta0_A, theta1_A) else 0.60
  alpha_B <- if (stationary) stationary_alpha(theta0_B, theta1_B) else 0.30

  params <- list(
    theta0_A = theta0_A, theta1_A = theta1_A, alpha_A = alpha_A,
    theta0_B = theta0_B, theta1_B = theta1_B, alpha_B = alpha_B,
    phi      = 0.5
  )
  if (model_type == "symmetric") params$pi <- 0.05

  params
}


#' Fit the 2-type FMM extension EM model
#'
#' Iterates the EM algorithm for the 2-type finite mixture model until
#' the relative change in observed-data log-likelihood drops below \code{tol}.
#' All M-step updates are in closed form.
#'
#' TeX ref: EM baseline.tex Section 7.
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric).
#' @param model_type Character: \code{"symmetric"} (default) or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}).
#' @param params0 Optional named list of starting parameters. If \code{NULL},
#'   \code{\link{init_params_fmm}} is used.
#' @param max_iter Maximum iterations (default 1000L).
#' @param tol Convergence tolerance on relative LL change (default 1e-8).
#' @param theta_cap Upper bound for transition probabilities (default 0.999).
#' @param pi_cap Upper bound for \eqn{\pi} (default 0.49).
#' @param verbose Integer: 0 = silent, 1 = summary (default 1L).
#' @return List with \code{params}, \code{loglik}, \code{history},
#'   \code{converged}, \code{iterations}, \code{gamma_A}, \code{gamma_B},
#'   \code{model_type}, \code{stationary}.
#' @examples
#' \dontrun{
#'   df  <- data.frame(y1 = rbinom(100,1,.6), y2 = rbinom(100,1,.6),
#'                     y3 = rbinom(100,1,.6), weight = rep(1,100))
#'   em_fit_fmm(df, verbose = 0L)
#' }
#' @export
em_fit_fmm <- function(df,
                        model_type = "symmetric",
                        stationary = TRUE,
                        params0    = NULL,
                        max_iter   = 1000L,
                        tol        = 1e-8,
                        theta_cap  = 0.999,
                        pi_cap     = 0.49,
                        verbose    = 1L) {
  if (!model_type %in% c("symmetric", "none"))
    stop("em_fit_fmm: model_type must be 'symmetric' or 'none'")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1L)
    stop("em_fit_fmm: max_iter must be a positive integer")
  .validate_panel_df(df)

  params <- params0 %||% init_params_fmm(model_type, stationary)

  if (model_type == "symmetric" && !is.null(params$pi)) {
    if (!is.numeric(params$pi) || length(params$pi) != 1L ||
        params$pi < 0 || params$pi >= 0.5)
      stop("em_fit_fmm: params$pi must be a scalar in [0, 0.5) for model_type='symmetric'")
  }

  # Validate data quality once
  for (.col in c("y1", "y2", "y3")) {
    if (!all(df[[.col]] %in% c(0, 1)))
      stop(sprintf("em_fit_fmm: column '%s' must be binary (0/1 only)", .col))
  }
  if (anyNA(df$weight) || any(df$weight <= 0))
    stop("em_fit_fmm: weight must be non-NA and strictly positive")

  ll_prev   <- -Inf
  converged <- FALSE
  history_rows <- vector("list", max_iter)

  # Precompute static mismatch indicator matrix (N x 8) for symmetric model.
  hmat_static <- .hmat_cache()
  h1_s <- hmat_static[, 1L]; h2_s <- hmat_static[, 2L]; h3_s <- hmat_static[, 3L]
  mm_static <- if (model_type == "symmetric") {
    outer(as.integer(df$y1), h1_s, "!=") +
    outer(as.integer(df$y2), h2_s, "!=") +
    outer(as.integer(df$y3), h3_s, "!=")
  } else {
    NULL
  }

  for (iter in seq_len(max_iter)) {
    estep_out <- e_step_fmm(df, params, model_type, validate = FALSE,
                             mm_precomp = mm_static)
    ll_new    <- estep_out$loglik

    history_rows[[iter]] <- c(
      iter      = iter,
      loglik    = ll_new,
      theta0_A  = params$theta0_A,
      theta1_A  = params$theta1_A,
      theta0_B  = params$theta0_B,
      theta1_B  = params$theta1_B,
      phi       = params$phi,
      pi        = if (model_type == "symmetric") params$pi else NA_real_
    )

    if (iter > 1L) {
      rel_change <- abs(ll_new - ll_prev) / (abs(ll_prev) + 1e-16)
      if (rel_change < tol) {
        converged <- TRUE
        break
      }
      if (ll_new < ll_prev - 1e-6)
        warning(sprintf("em_fit_fmm iter %d: LL decreased by %.2e.", iter, ll_prev - ll_new))
    }
    ll_prev <- ll_new

    params <- m_step_fmm(estep_out$suff, model_type, stationary,
                          theta_cap, pi_cap)
  }

  # Resolve label switching after convergence
  # NOTE: history columns reflect the labelling BEFORE this swap. If the swap
  # occurs, theta0_A/theta1_A in history rows correspond to what becomes type B
  # post-convergence. This is expected behaviour: the history tracks the
  # optimisation trajectory, not the final convention.
  params <- resolve_label_switching_fmm(params)

  if (verbose >= 1L) {
    status <- if (converged) "converged" else "did NOT converge"
    cat(sprintf(
      "em_fit_fmm [%s, stationary=%s]: %s in %d iters. loglik = %.4f\n",
      model_type, stationary, status, iter, ll_new
    ))
    cat(sprintf(
      "  Type A: theta0=%.4f, theta1=%.4f, alpha=%.4f\n",
      params$theta0_A, params$theta1_A, params$alpha_A
    ))
    cat(sprintf(
      "  Type B: theta0=%.4f, theta1=%.4f, alpha=%.4f\n",
      params$theta0_B, params$theta1_B, params$alpha_B
    ))
    cat(sprintf("  phi (type A share) = %.4f\n", params$phi))
    if (model_type == "symmetric")
      cat(sprintf("  pi = %.4f\n", params$pi))
  }

  if (!converged && verbose >= 1L)
    warning(sprintf(
      "em_fit_fmm [%s, stationary=%s]: did not converge in %d iterations.",
      model_type, stationary, max_iter
    ))

  final_estep <- e_step_fmm(df, params, model_type, validate = FALSE)
  history_df  <- as.data.frame(do.call(rbind, history_rows[seq_len(iter)]))

  list(
    params     = params,
    loglik     = final_estep$loglik,
    history    = history_df,
    converged  = converged,
    iterations = iter,
    gamma_A    = final_estep$gamma_A,
    gamma_B    = final_estep$gamma_B,
    model_type = model_type,
    stationary = stationary
  )
}


#' Evaluate observed-data log-likelihood for the FMM extension
#'
#' @param df Data frame with columns y1, y2, y3, weight.
#' @param params Named list of FMM parameters.
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @return Scalar observed-data log-likelihood.
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   compute_observed_loglik_fmm(df, init_params_fmm())
#' }
#' @export
compute_observed_loglik_fmm <- function(df, params,
                                         model_type = "symmetric") {
  e_step_fmm(df, params, model_type, validate = TRUE)$loglik
}
