# ==============================================================================
# EM-tenure: Implied quantity transforms for the contamination (epsilon) model
# Created: 2026-05-07
#
# Function `implied_tenure_contamination()` maps raw EM parameter lists for
# the eps (Spec I) model to economically and empirically interpretable quantities:
#
#   - Labour-market flow rates (entry, exit, employment)
#   - Misclassification and contamination rates (pi, eps)
#   - Spell-length interpretation (mean/median durations from lambda_g, lambda_d)
#   - Tenure measurement quality (P(clock-consistent), P(pair clean), etc.)
#
# These are computed for both point estimates and each bootstrap replicate, so
# that summarise_bootstrap() can directly produce SEs for all interpretable
# quantities (delta-method-free: each quantity is a 1-to-1 transform of the
# params, so bootstrap SDs are exact).
#
# Companion: documents/EM tenure epsilon.tex, Section 3.
# Used by: bootstrap_utils_tenure_contamination.R, build_tables_tenure_contamination.R
# ==============================================================================

#' Compute interpretable implied quantities from eps model parameter list
#'
#' Converts raw EM parameter estimates for the eps (Spec I) contamination model
#' into a set of economically and empirically interpretable quantities. The
#' function covers both the standard labour-market flow rates (shared with the
#' baseline model) and the tenure-specific outputs unique to the contamination
#' specification.
#'
#' All derived quantities are scalar; the function returns a flat named list
#' compatible with \code{summarise_bootstrap()}.
#'
#' @param params Named list of EM parameters. Must contain:
#'   \code{theta0}, \code{theta1}, \code{pi}, \code{eps},
#'   \code{lambda_g}, \code{lambda_d}.
#'   \code{alpha} is optional (if present it is returned as-is; if absent it is
#'   derived from the stationarity formula).
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{entry_rate}{\eqn{\theta_0}: job-finding probability per quarter.}
#'     \item{exit_rate}{\eqn{1 - \theta_1}: job-separation probability per quarter.}
#'     \item{employment_rate}{Steady-state employment rate
#'       \eqn{\theta_0 / (\theta_0 + 1 - \theta_1)}.}
#'     \item{pi}{Symmetric misclassification probability \eqn{\pi}.}
#'     \item{eps}{Per-wave tenure contamination probability \eqn{\varepsilon}.}
#'     \item{mean_spell_g_years}{Mean employment spell length in years:
#'       \eqn{1 / \lambda_g}.}
#'     \item{mean_spell_g_months}{Mean employment spell length in months:
#'       \eqn{12 / \lambda_g}.}
#'     \item{mean_spell_d_years}{Mean non-employment spell length in years:
#'       \eqn{1 / \lambda_d}.}
#'     \item{mean_spell_d_months}{Mean non-employment spell length in months:
#'       \eqn{12 / \lambda_d}.}
#'     \item{median_spell_g_months}{Median employment spell length in months:
#'       \eqn{12 \log(2) / \lambda_g}.}
#'     \item{p_clock_consistent}{Probability a single tenure report is
#'       clock-consistent (clean): \eqn{1 - \varepsilon}.}
#'     \item{p_pair_clean}{Probability both flanks of a K=2 spell are clean:
#'       \eqn{(1 - \varepsilon)^2}.}
#'     \item{p_triple_clean}{Probability all three waves of a K=3 spell are
#'       clean: \eqn{(1 - \varepsilon)^3}.}
#'     \item{contaminated_per_1000}{Expected number of contaminated tenure
#'       reports per 1,000 observations: \eqn{\varepsilon \times 1000}.}
#'   }
#'
#' @examples
#' p <- list(alpha = 0.53, theta0 = 0.10, theta1 = 0.92,
#'           pi = 0.088, eps = 0.171, lambda_g = 0.158, lambda_d = 0.5)
#' implied_tenure_contamination(p)
#'
#' @seealso \code{em_fit_tenure_eps()}, \code{bootstrap_one_eps()}
#' @export
implied_tenure_contamination <- function(params) {
  # --- Input validation ---
  required <- c("theta0", "theta1", "pi", "eps", "lambda_g", "lambda_d")
  missing_params <- setdiff(required, names(params))
  if (length(missing_params) > 0L) {
    stop(sprintf(
      "implied_tenure_contamination: missing required parameters: %s",
      paste(missing_params, collapse = ", ")
    ))
  }

  theta0   <- params$theta0
  theta1   <- params$theta1
  pi_val   <- params$pi
  eps_val  <- params$eps
  lambda_g <- params$lambda_g
  lambda_d <- params$lambda_d

  # Scalar and finite checks for core parameters
  for (nm in required) {
    val <- params[[nm]]
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop(sprintf(
        "implied_tenure_contamination: '%s' must be a finite scalar numeric; got: %s",
        nm, deparse(val)
      ))
    }
  }
  if (theta0 < 0 || theta0 > 1)
    stop(sprintf("implied_tenure_contamination: theta0 out of [0,1]: %g", theta0))
  if (theta1 < 0 || theta1 > 1)
    stop(sprintf("implied_tenure_contamination: theta1 out of [0,1]: %g", theta1))
  if (pi_val < 0 || pi_val > 1)
    stop(sprintf("implied_tenure_contamination: pi out of [0,1]: %g", pi_val))
  if (eps_val < 0 || eps_val > 1)
    stop(sprintf("implied_tenure_contamination: eps out of [0,1]: %g", eps_val))
  if (lambda_g < 0)
    stop(sprintf("implied_tenure_contamination: lambda_g must be >= 0; got: %g", lambda_g))
  if (lambda_d < 0)
    stop(sprintf("implied_tenure_contamination: lambda_d must be >= 0; got: %g", lambda_d))

  # --- Labour-market flow rates ---
  entry_rate      <- theta0
  exit_rate       <- 1 - theta1
  denom           <- theta0 + (1 - theta1)
  employment_rate <- if (denom > 0) theta0 / denom else NA_real_

  # --- Spell duration quantities (Exponential distribution) ---
  # Mean = 1/lambda; Median = log(2)/lambda. Return Inf when lambda = 0.
  mean_spell_g_years   <- if (lambda_g > 0) 1 / lambda_g else Inf
  mean_spell_g_months  <- if (lambda_g > 0) 12 / lambda_g else Inf
  mean_spell_d_years   <- if (lambda_d > 0) 1 / lambda_d else Inf
  mean_spell_d_months  <- if (lambda_d > 0) 12 / lambda_d else Inf
  median_spell_g_months <- if (lambda_g > 0) 12 * log(2) / lambda_g else Inf

  # --- Tenure measurement quality (contamination model quantities) ---
  # These are direct functions of eps, so bootstrap SDs are exact (no delta method).
  p_clean_1        <- 1 - eps_val          # P(single wave clean)
  p_pair_clean     <- (1 - eps_val)^2      # P(both flanks clean, K=2 spell)
  p_triple_clean   <- (1 - eps_val)^3      # P(all three waves clean, K=3 spell)
  contam_per_1000  <- eps_val * 1000       # expected contaminated per 1,000 obs

  list(
    entry_rate             = entry_rate,
    exit_rate              = exit_rate,
    employment_rate        = employment_rate,
    pi                     = pi_val,
    eps                    = eps_val,
    mean_spell_g_years     = mean_spell_g_years,
    mean_spell_g_months    = mean_spell_g_months,
    mean_spell_d_years     = mean_spell_d_years,
    mean_spell_d_months    = mean_spell_d_months,
    median_spell_g_months  = median_spell_g_months,
    p_clock_consistent     = p_clean_1,
    p_pair_clean           = p_pair_clean,
    p_triple_clean         = p_triple_clean,
    contaminated_per_1000  = contam_per_1000
  )
}
