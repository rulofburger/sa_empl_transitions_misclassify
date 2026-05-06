# ==============================================================================
# EM-baseline-ext: M-step for Extension III (2-type FMM)
# Created: 2026-05-06
#
# All M-step updates are in closed form for the FMM. Each type has its own
# transition probability parameters updated from type-specific sufficient
# statistics. The mixing weight phi and misclassification probability pi are
# updated from pooled statistics.
#
# TeX ref: EM baseline.tex Section 7, Eqs (23)--(27).
# ==============================================================================

#' M-step for the 2-type FMM extension
#'
#' Computes closed-form parameter updates for the 2-type finite mixture model.
#' Type-specific Markov parameters are updated independently from their
#' respective sufficient statistics. The mixing weight \eqn{\phi} and
#' misclassification probability \eqn{\pi} are updated from pooled statistics.
#'
#' TeX ref: EM baseline.tex Section 7, Eqs (23)--(27).
#'
#' @param suff Sufficient statistics list from \code{\link{e_step_fmm}}.
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @param stationary Logical. If \code{TRUE}, type-specific \eqn{\alpha^k} is
#'   derived from the stationarity restriction.
#' @param theta_cap Upper bound for transition probabilities (default 0.999).
#' @param pi_cap Upper bound for \eqn{\pi} (default 0.49).
#' @param phi_cap Upper bound for \eqn{\phi} (default 0.999, prevents
#'   degenerate single-type solutions).
#' @return Named list with updated \code{theta0_A}, \code{theta1_A},
#'   \code{alpha_A}, \code{theta0_B}, \code{theta1_B}, \code{alpha_B},
#'   \code{phi}, and optionally \code{pi}.
#' @examples
#' \dontrun{
#'   df   <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                      y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   suff <- e_step_fmm(df, init_params_fmm())$suff
#'   m_step_fmm(suff)
#' }
#' @export
m_step_fmm <- function(suff, model_type = "symmetric", stationary = TRUE,
                        theta_cap = 0.999, pi_cap = 0.49, phi_cap = 0.999) {
  if (!model_type %in% c("symmetric", "none"))
    stop("m_step_fmm: model_type must be 'symmetric' or 'none'")

  eps <- 1e-10

  # ---- Helper: update one type's transition parameters ---------------------
  .update_type <- function(T11, T01, D1, D0, C1, C0) {
    theta1 <- .bound01(T11 / max(D1, eps), eps = eps)
    theta1 <- min(theta1, theta_cap)
    theta0 <- .bound01(T01 / max(D0, eps), eps = eps)
    theta0 <- min(theta0, theta_cap)
    if (stationary) {
      alpha <- stationary_alpha(theta0, theta1)
    } else {
      alpha <- .bound01(C1 / max(C1 + C0, eps))
    }
    list(theta0 = theta0, theta1 = theta1, alpha = alpha)
  }

  pA <- .update_type(suff$T11_A, suff$T01_A, suff$D1_A, suff$D0_A,
                     suff$C1_A, suff$C0_A)
  pB <- .update_type(suff$T11_B, suff$T01_B, suff$D1_B, suff$D0_B,
                     suff$C1_B, suff$C0_B)

  # ---- Mixing weight -------------------------------------------------------
  phi <- suff$phi_num / max(suff$total_weight, eps)
  # Clamp symmetrically: phi in [1-phi_cap, phi_cap] so neither type dominates
  phi_floor <- 1 - phi_cap
  phi <- min(max(.bound01(phi, eps = eps), phi_floor), phi_cap)

  # ---- Misclassification (shared, pooled) ----------------------------------
  params_out <- list(
    theta0_A = pA$theta0, theta1_A = pA$theta1, alpha_A = pA$alpha,
    theta0_B = pB$theta0, theta1_B = pB$theta1, alpha_B = pB$alpha,
    phi      = phi
  )

  if (model_type == "symmetric") {
    pi <- suff$M / max(3 * suff$total_weight, eps)
    pi <- min(max(pi, 0), pi_cap)
    params_out$pi <- pi
  }

  params_out
}


#' Resolve FMM label switching
#'
#' Ensures that type A is the more stable type (higher employment persistence)
#' by swapping the two types if \eqn{\theta_1^B > \theta_1^A}.
#'
#' @param params Named list of FMM parameters (as returned by
#'   \code{\link{m_step_fmm}}).
#' @return Possibly relabelled params list.
#' @export
resolve_label_switching_fmm <- function(params) {
  if (params$theta1_B > params$theta1_A) {
    # Swap type labels and adjust phi
    params_new <- params
    params_new$theta0_A <- params$theta0_B
    params_new$theta1_A <- params$theta1_B
    params_new$alpha_A  <- params$alpha_B
    params_new$theta0_B <- params$theta0_A
    params_new$theta1_B <- params$theta1_A
    params_new$alpha_B  <- params$alpha_A
    params_new$phi      <- 1 - params$phi
    return(params_new)
  }
  params
}
