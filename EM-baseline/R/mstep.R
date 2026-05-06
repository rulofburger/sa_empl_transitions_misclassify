# ==============================================================================
# EM-baseline: M-step — closed-form parameter updates
# Created: 2026-05-05
#
# All M-step updates are in closed form (no numerical optimisation needed).
# The complete-data log-likelihood decomposes additively into two blocks:
# (1) Markov prior block -> updates theta0, theta1 (and alpha)
# (2) Misclassification block -> updates pi (symmetric) or pi0, pi1 (asymmetric)
#
# TeX reference: EM baseline.tex, Eqs (18)-(19) for symmetric; Eq (23)
# for asymmetric; Eq (4) for stationarity restriction.
# ==============================================================================

#' M-step: compute closed-form parameter updates
#'
#' Given the sufficient statistics from \code{\link{e_step}}, computes the
#' M-step parameter updates. All updates are in closed form.
#'
#' TeX ref: Eqs (18)-(19) for symmetric; Eq (23) for asymmetric.
#'
#' @param suff Named list of sufficient statistics returned by
#'   \code{\link{e_step}}: C1, C0, D1, D0, T11, T01, M, M0, M1, H0, H1,
#'   total_weight.
#' @param model_type Character: \code{"symmetric"} (default),
#'   \code{"asymmetric"}, or \code{"none"}.
#' @param stationary Logical (default \code{TRUE}). If \code{TRUE}, alpha is
#'   derived from the stationarity restriction
#'   alpha = theta0 / (theta0 + 1 - theta1). If \code{FALSE}, alpha is
#'   estimated freely as C1 / (C1 + C0).
#' @param theta_cap Upper bound for theta0 and theta1 (default 0.999).
#' @param pi_cap Upper bound for pi, pi0, pi1 (default 0.49).
#' @return Named list of updated parameter estimates. Fields:
#'   \code{alpha}, \code{theta0}, \code{theta1}, and (depending on
#'   \code{model_type}): \code{pi}, \code{pi0} + \code{pi1}, or neither.
#' @examples
#' suff <- list(
#'   C1 = 30, C0 = 20, D1 = 80, D0 = 60, T11 = 72, T01 = 6,
#'   M = 5, M0 = 2, M1 = 3, H0 = 60, H1 = 90, total_weight = 50
#' )
#' m_step(suff, model_type = "symmetric")
#' @export
m_step <- function(suff, model_type = "symmetric", stationary = TRUE,
                   theta_cap = 0.999, pi_cap = 0.49) {
  .validate_model_type(model_type)  # fail-fast before any computation
  eps <- 1e-10

  # --- Transition probabilities (TeX Eq 18) ----------------------------------
  # theta1 = T11 / D1  (fraction of from-employment periods ending in employment)
  # theta0 = T01 / D0  (fraction of from-nonemployment periods ending in employment)
  if (suff$D1 < 1e-6)
    warning(sprintf("m_step: D1 = %.2e \u2248 0 \u2014 theta1 is unidentified; estimate unreliable", suff$D1))
  if (suff$D0 < 1e-6)
    warning(sprintf("m_step: D0 = %.2e \u2248 0 \u2014 theta0 is unidentified; estimate unreliable", suff$D0))
  theta1 <- .bound01(suff$T11 / max(suff$D1, eps), eps = eps)
  theta1 <- min(theta1, theta_cap)

  theta0 <- .bound01(suff$T01 / max(suff$D0, eps), eps = eps)
  theta0 <- min(theta0, theta_cap)

  # --- Initial employment probability (TeX Eq 4 vs free) ---------------------
  if (stationary) {
    # Stationarity restriction: alpha = theta0 / (theta0 + 1 - theta1)
    alpha <- stationary_alpha(theta0, theta1)
  } else {
    # Free alpha: MLE from initial-state mass
    alpha <- .bound01(suff$C1 / max(suff$C1 + suff$C0, eps))
  }

  params_out <- list(alpha = alpha, theta0 = theta0, theta1 = theta1)

  # --- Misclassification parameters ------------------------------------------
  if (model_type == "symmetric") {
    # pi = M / (3 * W)  (TeX Eq 19)
    pi <- suff$M / max(3 * suff$total_weight, eps)
    pi <- min(.bound01(pi, eps = 0), pi_cap)  # enforce [0, pi_cap]
    params_out$pi <- pi

  } else if (model_type == "asymmetric") {
    # pi0 = M0 / H0  (false-positive rate, TeX Eq 23)
    # pi1 = M1 / H1  (false-negative rate, TeX Eq 23)
    pi0 <- suff$M0 / max(suff$H0, eps)
    pi0 <- min(.bound01(pi0, eps = 0), pi_cap)
    pi1 <- suff$M1 / max(suff$H1, eps)
    pi1 <- min(.bound01(pi1, eps = 0), pi_cap)
    params_out$pi0 <- pi0
    params_out$pi1 <- pi1

  } else if (model_type == "none") {
    # No misclassification parameters — nothing to add

  } else {
    # Unreachable: .validate_model_type() above already stopped for unknown values.
    stop(sprintf("Unknown model_type: '%s'. Use 'symmetric', 'asymmetric', or 'none'.", model_type))  # nocov
  }

  params_out
}
