# ==============================================================================
# EM-baseline-ext: EM driver for Extension IV (inconsistency-augmented model)
# Created: 2026-05-06
#
# TeX ref: EM baseline.tex Section 8.
# ==============================================================================

#' Initialise parameters for the inconsistency-augmented model
#'
#' Sets the intercept \eqn{\delta_0} so that the baseline (no-inconsistency)
#' misclassification probability equals 0.05. Slope coefficients
#' \eqn{\delta_1, \delta_2} start at zero (no inconsistency effect).
#'
#' @param stationary Logical (default \code{TRUE}).
#' @return Named list with \code{theta0}, \code{theta1}, \code{alpha},
#'   \code{delta} (length-3 vector).
#' @examples
#' init_params_inconsistency()
#' @export
init_params_inconsistency <- function(stationary = TRUE) {
  theta0 <- 0.10
  theta1 <- 0.90
  alpha  <- if (stationary) stationary_alpha(theta0, theta1) else 0.5

  # delta_0 = qlogis(2 * pi_base) so that pi_base = 0.5 * plogis(delta_0)
  # With pi_base = 0.05: delta_0 = qlogis(0.10) ~= -2.197
  delta <- c(qlogis(2 * 0.05), 0, 0)

  list(theta0 = theta0, theta1 = theta1, alpha = alpha, delta = delta)
}


#' Fit the inconsistency-augmented EM model
#'
#' Iterates GEM until the relative change in observed-data log-likelihood
#' drops below \code{tol}. Theta parameters are updated in closed form;
#' \eqn{\delta} is updated by one Fisher-scoring Newton-Raphson step per
#' iteration (GEM guarantee).
#'
#' TeX ref: EM baseline.tex Section 8.
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric).
#' @param incons_mat N × 6 numeric matrix from \code{\link{compute_inconsistencies}}.
#' @param stationary Logical (default \code{TRUE}).
#' @param params0 Optional named list of starting parameters. If \code{NULL},
#'   \code{\link{init_params_inconsistency}} is used.
#' @param max_iter Maximum iterations (default 1000L).
#' @param tol Convergence tolerance on relative LL change (default 1e-8).
#' @param theta_cap Upper bound for transition probabilities (default 0.999).
#' @param verbose Integer: 0 = silent, 1 = convergence summary (default 1L).
#' @return List with \code{params}, \code{loglik}, \code{history},
#'   \code{converged}, \code{iterations}, \code{gamma},
#'   \code{stationary}.
#' @examples
#' \dontrun{
#'   df      <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                         y3 = rbinom(50,1,.6), weight = rep(1,50),
#'                         age1=25:74, age2=26:75, age3=27:76,
#'                         educ1=rep(3L,50), educ2=rep(3L,50), educ3=rep(3L,50))
#'   inc_mat <- as.matrix(compute_inconsistencies(df)[,
#'     c("Y_age_1","Y_age_2","Y_age_3","Y_edu_1","Y_edu_2","Y_edu_3")])
#'   em_fit_inconsistency(df, inc_mat, verbose = 0L)
#' }
#' @export
em_fit_inconsistency <- function(df,
                                  incons_mat,
                                  stationary = TRUE,
                                  params0    = NULL,
                                  max_iter   = 1000L,
                                  tol        = 1e-8,
                                  theta_cap  = 0.999,
                                  verbose    = 1L) {
  .validate_panel_df(df)
  if (!is.matrix(incons_mat) || ncol(incons_mat) != 6L || nrow(incons_mat) != nrow(df))
    stop("em_fit_inconsistency: incons_mat must be an N x 6 matrix with same rows as df")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1L)
    stop("em_fit_inconsistency: max_iter must be a positive integer")
  # Validate incons_mat column order up-front (wrong order silently corrupts eta_mat)
  if (!is.null(colnames(incons_mat))) {
    expected_cols <- c("Y_age_1", "Y_age_2", "Y_age_3",
                       "Y_edu_1", "Y_edu_2", "Y_edu_3")
    if (!identical(colnames(incons_mat), expected_cols))
      stop(paste0(
        "em_fit_inconsistency: incons_mat columns out of order. ",
        "Expected: ", paste(expected_cols, collapse = ", ")
      ))
  }

  params <- params0 %||% init_params_inconsistency(stationary)

  if (length(params$delta) != 3L)
    stop("em_fit_inconsistency: params$delta must have length 3")

  for (.col in c("y1", "y2", "y3")) {
    if (!all(df[[.col]] %in% c(0, 1)))
      stop(sprintf("em_fit_inconsistency: column '%s' must be binary", .col))
  }
  if (anyNA(df$weight) || any(df$weight <= 0))
    stop("em_fit_inconsistency: weight must be non-NA and strictly positive")

  ll_prev   <- -Inf
  converged <- FALSE
  history_rows <- vector("list", max_iter)

  for (iter in seq_len(max_iter)) {
    estep_out <- e_step_inconsistency(df, incons_mat, params, validate = FALSE)
    ll_new    <- estep_out$loglik

    history_rows[[iter]] <- c(
      iter    = iter,
      loglik  = ll_new,
      theta0  = params$theta0,
      theta1  = params$theta1,
      alpha   = params$alpha,
      delta_0 = params$delta[1L],
      delta_1 = params$delta[2L],
      delta_2 = params$delta[3L],
      pi_base = 0.5 * plogis(params$delta[1L])
    )

    if (iter > 1L) {
      rel_change <- abs(ll_new - ll_prev) / (abs(ll_prev) + 1e-16)
      if (rel_change < tol) {
        converged <- TRUE
        break
      }
      if (ll_new < ll_prev - 1e-6)
        warning(sprintf(
          "em_fit_inconsistency iter %d: LL decreased by %.2e.", iter, ll_prev - ll_new
        ))
    }
    ll_prev <- ll_new

    params <- m_step_inconsistency(
      suff       = estep_out$suff,
      incons_mat = incons_mat,
      params_old = params,
      stationary = stationary,
      theta_cap  = theta_cap
    )
  }

  if (verbose >= 1L) {
    status  <- if (converged) "converged" else "did NOT converge"
    pi_base <- 0.5 * plogis(params$delta[1L])
    cat(sprintf(
      "em_fit_inconsistency [stationary=%s]: %s in %d iters. loglik = %.4f\n",
      stationary, status, iter, ll_new
    ))
    cat(sprintf(
      "  theta0=%.4f, theta1=%.4f, alpha=%.4f\n",
      params$theta0, params$theta1, params$alpha
    ))
    cat(sprintf(
      "  delta = (%.4f, %.4f, %.4f)  [pi_base = %.4f]\n",
      params$delta[1L], params$delta[2L], params$delta[3L], pi_base
    ))
  }

  if (!converged && verbose >= 1L)
    warning(sprintf(
      "em_fit_inconsistency [stationary=%s]: did not converge in %d iterations.",
      stationary, max_iter
    ))

  final_estep <- e_step_inconsistency(df, incons_mat, params, validate = FALSE)
  history_df  <- as.data.frame(do.call(rbind, history_rows[seq_len(iter)]))

  list(
    params     = params,
    loglik     = final_estep$loglik,
    history    = history_df,
    converged  = converged,
    iterations = iter,
    gamma      = final_estep$gamma,
    stationary = stationary
  )
}


#' Evaluate observed-data log-likelihood for the inconsistency model
#'
#' @param df Data frame with columns y1, y2, y3, weight.
#' @param incons_mat N × 6 matrix from \code{\link{compute_inconsistencies}}.
#' @param params Named list with theta0, theta1, alpha, delta.
#' @return Scalar observed-data log-likelihood.
#' @examples
#' \dontrun{
#'   df      <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                         y3 = rbinom(50,1,.6), weight = rep(1,50),
#'                         age1=25:74, age2=26:75, age3=27:76,
#'                         educ1=rep(3L,50), educ2=rep(3L,50), educ3=rep(3L,50))
#'   inc_mat <- as.matrix(compute_inconsistencies(df)[,
#'     c("Y_age_1","Y_age_2","Y_age_3","Y_edu_1","Y_edu_2","Y_edu_3")])
#'   compute_observed_loglik_inconsistency(df, inc_mat, init_params_inconsistency())
#' }
#' @export
compute_observed_loglik_inconsistency <- function(df, incons_mat, params) {
  e_step_inconsistency(df, incons_mat, params, validate = TRUE)$loglik
}
