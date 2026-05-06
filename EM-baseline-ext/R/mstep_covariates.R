# ==============================================================================
# EM-baseline-ext: M-step (GEM) for Extension I (observable heterogeneity)
# Created: 2026-05-06
#
# Performs a single IRLS (iteratively reweighted least squares) step for the
# weighted probit regression that updates beta0 and beta1. This is a GEM
# (Generalized EM) step: it strictly increases the Q function without finding
# the exact maximiser, guaranteeing monotone ascent of the observed-data LL.
#
# TeX ref: EM baseline.tex Section 5, Eqs (12)--(14).
# ==============================================================================

# Ridge constant for numerical stability of X'WX inversion
.IRLS_RIDGE <- 1e-6

#' M-step (GEM) for the covariate extension
#'
#' Given sufficient statistics from \code{\link{e_step_covariates}}, performs
#' one IRLS step for each of the two probit regressions (\eqn{\beta_0} and
#' \eqn{\beta_1}) and updates \eqn{\pi} (if applicable) in closed form.
#'
#' \strong{IRLS for weighted probit (GEM step)}:
#' Given current coefficient \eqn{\beta}, compute working quantities:
#' \deqn{\mu_i = \Phi(\eta_i), \quad \eta_i = x_i^\top \beta}
#' \deqn{z_i = \eta_i + (y_i - \mu_i)/\phi(\eta_i)}
#' \deqn{v_i = \phi(\eta_i)^2 / (\mu_i(1-\mu_i))}
#' then solve the weighted least squares:
#' \deqn{\beta_{\mathrm{new}} = (X^\top D X + \lambda I)^{-1} X^\top D z}
#' where \eqn{D = \mathrm{diag}(w_i \cdot v_i)} and \eqn{\lambda} is a
#' small ridge constant for numerical stability.
#'
#' TeX ref: EM baseline.tex Section 5, Eqs (12)--(14).
#'
#' @param suff Sufficient statistics list from \code{\link{e_step_covariates}}.
#' @param X N × p design matrix (same as passed to \code{e_step_covariates}).
#' @param params_old Named list with current \code{beta0}, \code{beta1},
#'   and optionally \code{pi}.
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @param stationary Logical. If \code{TRUE}, \eqn{\alpha} is not a free
#'   parameter (derived from individual-level stationarity). If \code{FALSE},
#'   returns a scalar \eqn{\alpha = C_1 / (C_1 + C_0)}.
#' @param pi_cap Upper bound for \eqn{\pi} (default 0.49).
#' @return Named list with updated \code{beta0}, \code{beta1},
#'   and optionally \code{pi} and \code{alpha} (free only).
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   X      <- matrix(1, nrow = 50, ncol = 1)
#'   params <- init_params_covariates(p = 1L)
#'   suff   <- e_step_covariates(df, X, params)$suff
#'   m_step_covariates(suff, X, params)
#' }
#' @export
m_step_covariates <- function(suff, X, params_old, model_type = "symmetric",
                               stationary = TRUE, pi_cap = 0.49) {
  if (!model_type %in% c("symmetric", "none"))
    stop("m_step_covariates: model_type must be 'symmetric' or 'none'")

  # ---- GEM: one IRLS step for beta1 (from state 1 -> persistence) ----------
  beta1_new <- .irls_probit_step(
    X          = X,
    eff_w      = suff$eff_w_1,
    eff_wy     = suff$eff_wy_1,
    beta_old   = params_old$beta1
  )

  # ---- GEM: one IRLS step for beta0 (from state 0 -> job finding) ----------
  beta0_new <- .irls_probit_step(
    X          = X,
    eff_w      = suff$eff_w_0,
    eff_wy     = suff$eff_wy_0,
    beta_old   = params_old$beta0
  )

  params_out <- list(beta0 = beta0_new, beta1 = beta1_new)

  # ---- pi update (closed form, TeX Eq 19) ----------------------------------
  if (model_type == "symmetric") {
    eps <- 1e-10
    pi  <- suff$M / max(3 * suff$total_weight, eps)
    pi  <- min(max(pi, 0), pi_cap)
    params_out$pi <- pi
  }

  # ---- Free alpha (scalar approximation) -----------------------------------
  if (!stationary) {
    eps   <- 1e-10
    alpha <- suff$C1 / max(suff$C1 + suff$C0, eps)
    alpha <- max(min(alpha, 1 - eps), eps)
    params_out$alpha <- alpha
  }

  params_out
}


# ------------------------------------------------------------------------------
# Internal: single IRLS step for probit
# ------------------------------------------------------------------------------

#' One IRLS step for weighted probit regression (GEM sub-routine)
#'
#' Performs a single iteration of IRLS for a weighted probit model where the
#' response is fractional (pseudo-observations from the E-step).
#'
#' @param X N × p design matrix.
#' @param eff_w N-vector: effective weights (responsibility-weighted transition
#'   counts at each individual, summed over both transitions).
#' @param eff_wy N-vector: effective weighted outcome numerators (responsibility
#'   × 1{h_t = 1} summed over transitions).
#' @param beta_old Current p-vector of probit coefficients.
#' @return Updated p-vector of probit coefficients.
#' @keywords internal
.irls_probit_step <- function(X, eff_w, eff_wy, beta_old) {
  N <- nrow(X)
  p <- ncol(X)

  # Mask: individuals with effective weight near zero contribute nothing
  active <- eff_w > 1e-10
  if (!any(active)) return(beta_old)  # no contributing observations

  # Fractional outcome at each individual (weighted average)
  # y_frac = eff_wy / eff_w for active rows; clamp to (0,1) for stability
  y_frac     <- rep(0.5, N)
  y_frac[active] <- eff_wy[active] / eff_w[active]
  y_frac     <- pmin(pmax(y_frac, 1e-6), 1 - 1e-6)

  # Linear predictor from current beta
  eta <- as.vector(X %*% beta_old)

  # IRLS working quantities (probit)
  mu      <- pnorm(eta)
  mu      <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
  phi_eta <- dnorm(eta)
  phi_eta <- pmax(phi_eta, 1e-10)

  # Working weights: v_i = phi(eta)^2 / (mu * (1 - mu))
  v <- (phi_eta^2) / (mu * (1 - mu))

  # Working response: z_i = eta_i + (y_i - mu_i) / phi(eta_i)
  z <- eta + (y_frac - mu) / phi_eta

  # Combined weight: observation weight * IRLS weight
  d <- eff_w * v  # N-vector

  # Weighted least squares: beta_new = (X' D X + ridge I)^{-1} X' D z
  # crossprod avoids materialising the p×N transposed intermediate.
  Xd   <- X * d               # N × p scaled X
  XtDX <- crossprod(Xd, X)    # p × p, no p×N temp
  XtDz <- crossprod(X, d * z) # p × 1

  # Ridge regularisation for numerical stability
  XtDX_reg <- XtDX + diag(.IRLS_RIDGE, p, p)

  beta_new <- tryCatch(
    as.vector(solve(XtDX_reg, XtDz)),
    error = function(e) {
      warning(paste0(
        ".irls_probit_step: solve() failed (", e$message,
        "). Returning current beta unchanged."
      ))
      beta_old
    }
  )

  beta_new
}
