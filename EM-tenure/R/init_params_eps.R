# ==============================================================================
# EM-tenure: Initial parameters for the eps (Spec I) model
# ==============================================================================
# Created: 2026-04-30
#
# init_params_eps() builds a starting parameter vector for the eps model.
# Derived from init_params() by:
#   - Removing sigma2_g (no Gaussian component in Spec I).
#   - Adding eps (per-wave contamination probability).
#   - Initialising lambda_g from the inverse-mean of observed tenures
#     (the Exp MLE), or via the CTMC link if linked = TRUE.
# ==============================================================================

#' Initialise parameters for the eps (Spec I) EM algorithm
#'
#' @param df Data frame with y1-y3, tenure1-tenure3, timegap_cat1-3.
#' @param linked Logical (default FALSE). If TRUE, initialise lambda_g and
#'   lambda_d from the CTMC link (lambda = f(theta)).
#' @param eps_init Starting value for eps in (0, 1). Default 0.20, motivated
#'   by the empirical 63.2% clock-consistent E-spell pairs (K = 2 implies
#'   eps = 1 - sqrt(0.632) ~ 0.205).
#' @return Named list: alpha, theta0, theta1, pi, eps, lambda_g, lambda_d.
#' @references TeX: \emph{EM tenure epsilon.tex}; design memo
#'   \code{EM-tenure/feedback/2026-04-30-epsilon-spec-design.md}.
#' @examples
#' \dontrun{
#' p0 <- init_params_eps(df_qlfs, linked = FALSE, eps_init = 0.20)
#' str(p0)
#' }
#' @export
init_params_eps <- function(df, linked = FALSE, eps_init = 0.20) {
  if (!is.finite(eps_init) || eps_init <= 0 || eps_init >= 1) {
    stop(sprintf("eps_init must be in (0, 1); got %.4g", eps_init))
  }

  req <- c("y1", "y2", "y3", "tenure1", "tenure2", "tenure3")
  missing_req <- setdiff(req, names(df))
  if (length(missing_req) > 0) {
    stop("init_params_eps: missing columns: ", paste(missing_req, collapse = ", "))
  }

  # Delegate to the shared initializer; strip Gaussian parameters (no sigma
  # component in Spec I). alpha = empirical employment rate; theta1/theta0 =
  # hardcoded defaults (0.90/0.10), not moment estimates.
  base <- init_params(df, discrete_timegap = TRUE, linked = linked)
  # Return in canonical eps-model order; eps is inserted at position 5.
  list(
    alpha    = base$alpha,
    theta0   = base$theta0,
    theta1   = base$theta1,
    pi       = base$pi,
    eps      = eps_init,
    lambda_g = base$lambda_g,
    lambda_d = base$lambda_d
  )
}
