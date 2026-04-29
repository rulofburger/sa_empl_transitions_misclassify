# ==============================================================================
# EM-tenure: Internal utility functions
# ==============================================================================
# These are low-level helpers used across all modules. They are not exported
# for direct use but are essential for numerical stability and convenience.
# ==============================================================================

#' Quarterly increment in years. All modules derive the inter-wave
#' time step from this single constant (TeX: Delta = 3 months = 0.25 years).
#' @keywords internal
.QUARTER_YEARS <- 0.25

# ==============================================================================
# Timegap category constants (discrete/censored timegap model)
# ==============================================================================
# The QLFS timegap variable is categorical with 7 observed categories
# (codes 1–7). Code 0 is "never worked" (NA after ingest). Codes 8 and 99
# are missing/unspecified (also NA).
#
# Category boundaries (in years):
#   k=1: [0,     0.25)   =  [0,  3) months
#   k=2: [0.25,  0.50)   =  [3,  6) months
#   k=3: [0.50,  0.75)   =  [6,  9) months
#   k=4: [0.75,  1.00)   = [9,  12) months
#   k=5: [1.00,  3.00)   = [12, 36) months
#   k=6: [3.00,  5.00)   = [36, 60) months
#   k=7: [5.00, Inf)     = [60, Inf) months
#
# TeX reference: DIAGNOSIS.md, Issue 3 Remedy (category interval table)
# ==============================================================================

#' Left boundaries of timegap categories (plus Inf sentinel).
#'
#' Length-8 vector: the 7 left endpoints followed by Inf (right endpoint of
#' the last category). So .TIMEGAP_BOUNDS_YEARS[k] = a_k and
#' .TIMEGAP_BOUNDS_YEARS[k+1] = b_k for category k in 1:7.
#' @keywords internal
.TIMEGAP_BOUNDS_YEARS <- c(0, 0.25, 0.5, 0.75, 1.0, 3.0, 5.0, Inf)

#' Number of timegap categories.
#' @keywords internal
.N_TIMEGAP_CATS <- 7L

#' Midpoints of timegap categories in months (for descriptive use and
#' backward-compatible continuous timegap midpoint mapping).
#' These midpoints correspond to the raw QLFS codes 1–7.
#' @keywords internal
.TIMEGAP_MIDPOINTS_MONTHS <- c(1.5, 4.5, 7.5, 10.5, 24.0, 48.0, 90.0)

#' Floor for wrong-state duration imputation (nearest-non-zero, all-same-state
#' fallback). Equals 0.25 months / 12 ≈ 0.021 years ≈ 1 week.
#'
#' Used when an individual has no valid donor wave for nearest-non-zero
#' imputation (e.g., employed at all three waves — no nonemployment duration
#' available at any wave). The very small floor keeps the EMG density finite
#' (so misclassification is not structurally ruled out) but assigns low
#' likelihood (so the EM does not over-attribute misclassification).
#' @keywords internal
.DURATION_FLOOR <- 0.25 / 12  # ≈ 0.021 years

#' Return the interval [a_k, b_k) for timegap category k.
#'
#' @param k Integer category code in 1:7.
#' @return Numeric vector of length 2: c(a_k, b_k). b_k = Inf for k = 7.
#' @keywords internal
.timegap_interval <- function(k) {
  stopifnot(length(k) == 1L, k >= 1L, k <= .N_TIMEGAP_CATS)
  c(.TIMEGAP_BOUNDS_YEARS[k], .TIMEGAP_BOUNDS_YEARS[k + 1L])
}

#' Map a continuous nonemployment duration (years) to a timegap category code.
#'
#' Uses the left-closed, right-open interval convention: category k covers
#' [.TIMEGAP_BOUNDS_YEARS[k], .TIMEGAP_BOUNDS_YEARS[k+1]).
#'
#' @param d_years Numeric vector of durations in years (non-negative).
#' @return Integer vector of category codes in 1:7. Values below 0 or NA
#'   return NA.
#' @keywords internal
.continuous_to_cat <- function(d_years) {
  # findInterval with left.open = FALSE gives the index of the largest
  # boundary <= d, which equals the category number.
  # Boundaries used: the first 7 entries (left endpoints of cats 1–7).
  # Categories are 1-indexed so no offset needed.
  cats <- findInterval(d_years, .TIMEGAP_BOUNDS_YEARS[seq_len(.N_TIMEGAP_CATS)],
                       left.open = FALSE)
  # Cap at 7 (any d >= 5 gets category 7)
  cats <- pmin(cats, .N_TIMEGAP_CATS)
  # Values below 0 or NA stay as 0 → set to NA
  cats[cats == 0L | is.na(d_years)] <- NA_integer_
  as.integer(cats)
}

#' Null coalescing operator
#' Returns `a` if non-NULL, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

#' Clip probabilities to (eps, 1-eps)
#'
#' Prevents log(0) and log(1) in likelihood computations.
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive bound (default 1e-10).
#' @return Clipped probabilities in (eps, 1-eps).
#' @keywords internal
.bound01 <- function(p, eps = 1e-10) {
  pmin(1 - eps, pmax(eps, p))
}

#' Enforce strict positivity with a floor
#'
#' Used for variances, rates, and denominators to avoid division by zero.
#'
#' @param x Numeric vector.
#' @param eps Small positive floor (default 1e-12).
#' @return Values >= eps.
#' @keywords internal
.pos <- function(x, eps = 1e-12) {
  pmax(eps, x)
}

#' Safe logit transform
#'
#' Returns NA at boundaries 0 and 1 (where logit is undefined).
#'
#' @param p Probability in (0, 1).
#' @return log(p/(1-p)), or NA if p is on boundary.
#' @keywords internal
.safe_logit <- function(p) {
  ifelse(p <= 0 | p >= 1, NA_real_, log(p / (1 - p)))
}

#' Numerically stable log-sum-exp
#'
#' Computes log(sum(exp(x))) without overflow/underflow.
#' Shifts by max(x) before exponentiating.
#'
#' @param x Numeric vector (may contain -Inf).
#' @return Scalar: log(sum(exp(x))).
#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)  # all -Inf
  m + log(sum(exp(x - m)))
}

#' Complementary error function
#'
#' Computes erfc(x) = 1 - erf(x) = 2 * Phi(-sqrt(2) * x),
#' where Phi is the standard Normal CDF. Numerically stable for all finite x.
#'
#' @param x Numeric vector.
#' @return erfc(x).
#' @keywords internal
erfc <- function(x) {
  2 * pnorm(-sqrt(2) * x)
}

# ==============================================================================
# Warm-start utilities
# ==============================================================================

#' Load per-model warm-start parameters from the most recent saved results
#'
#' Scans \code{results_dir} for the most recent \code{.rds} file for each of
#' the 8 model variants and builds a named list of initial-parameter objects
#' ready to pass as \code{params0} to \code{em_fit_tenure()} or
#' \code{em_fit_tenure_rho()}. When no prior result exists for a model, the
#' corresponding list element is \code{NULL} (caller should fall back to
#' \code{init_params()} / \code{init_params_rho()} defaults).
#'
#' File-name convention (set by the pipeline's \code{saveRDS()} calls):
#' \preformatted{
#'   fit_miscl_YYYYMMDD_HHMMSS.rds
#'   fit_miscl_stationary_YYYYMMDD_HHMMSS.rds
#'   fit_miscl_linked_YYYYMMDD_HHMMSS.rds
#'   fit_miscl_stationary_linked_YYYYMMDD_HHMMSS.rds
#'   fit_rho_YYYYMMDD_HHMMSS.rds
#'   fit_rho_stationary_YYYYMMDD_HHMMSS.rds
#'   fit_rho_linked_YYYYMMDD_HHMMSS.rds
#'   fit_rho_stationary_linked_YYYYMMDD_HHMMSS.rds
#' }
#'
#' @param results_dir Character. Path to the directory containing saved
#'   \code{.rds} result files (relative to the working directory or absolute).
#' @param df Data frame. Passed to \code{init_params()} / \code{init_params_rho()}
#'   to build the baseline parameter structure before overwriting with loaded
#'   values. Must contain the standard model columns.
#' @param verbose Logical (default \code{TRUE}). If \code{TRUE}, prints a
#'   one-line status message per model.
#' @return Named list with elements:
#'   \describe{
#'     \item{miscl}{Params for free non-stationary model, or \code{NULL}.}
#'     \item{miscl_stationary}{Params for free stationary model, or \code{NULL}.}
#'     \item{miscl_linked}{Params for CTMC-linked non-stationary model, or \code{NULL}.}
#'     \item{miscl_stationary_linked}{Params for CTMC-linked stationary model, or \code{NULL}.}
#'     \item{rho}{Params for rho free non-stationary model, or \code{NULL}.}
#'     \item{rho_stationary}{Params for rho free stationary model, or \code{NULL}.}
#'     \item{rho_linked}{Params for rho CTMC-linked non-stationary model, or \code{NULL}.}
#'     \item{rho_stationary_linked}{Params for rho CTMC-linked stationary model, or \code{NULL}.}
#'   }
#' @export
load_warm_starts <- function(results_dir, df, verbose = TRUE) {
  # Map: model label → list(prefix, is_linked, is_rho)
  # Prefixes are ordered from most specific to least specific so that, e.g.,
  # "fit_miscl_stationary_linked_" does not match "fit_miscl_" accidentally.
  .model_specs <- list(
    miscl_stationary_linked  = list(prefix = "fit_miscl_stationary_linked_",  linked = TRUE,  rho = FALSE),
    miscl_linked             = list(prefix = "fit_miscl_linked_",             linked = TRUE,  rho = FALSE),
    miscl_stationary         = list(prefix = "fit_miscl_stationary_",         linked = FALSE, rho = FALSE),
    miscl                    = list(prefix = "fit_miscl_",                    linked = FALSE, rho = FALSE),
    rho_stationary_linked    = list(prefix = "fit_rho_stationary_linked_",    linked = TRUE,  rho = TRUE),
    rho_linked               = list(prefix = "fit_rho_linked_",               linked = TRUE,  rho = TRUE),
    rho_stationary           = list(prefix = "fit_rho_stationary_",           linked = FALSE, rho = TRUE),
    rho                      = list(prefix = "fit_rho_",                      linked = FALSE, rho = TRUE)
  )

  out <- vector("list", length(.model_specs))
  names(out) <- names(.model_specs)

  for (.label in names(.model_specs)) {
    spec    <- .model_specs[[.label]]
    pattern <- paste0("^", spec$prefix, "\\d{8}_\\d{6}\\.rds$")
    files   <- sort(list.files(results_dir, pattern = pattern, full.names = TRUE))

    if (length(files) == 0L) {
      if (verbose) {
        message(sprintf("  [warm-start] %-30s — no prior run found (using defaults)", .label))
      }
      out[[.label]] <- NULL
      next
    }

    .latest <- tail(files, 1L)
    .loaded <- tryCatch(
      readRDS(.latest),
      error = function(e) {
        warning(sprintf("load_warm_starts: could not read '%s': %s", .latest, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(.loaded)) {
      if (verbose) {
        message(sprintf("  [warm-start] %-30s — file unreadable (using defaults)", .label))
      }
      out[[.label]] <- NULL
      next
    }

    # Build a type-appropriate baseline, then overwrite with loaded keys.
    .base <- if (spec$rho) {
      init_params_rho(df, linked = spec$linked)
    } else {
      init_params(df, discrete_timegap = TRUE, linked = spec$linked)
    }
    for (.k in intersect(names(.loaded$params), names(.base))) {
      .base[[.k]] <- .loaded$params[[.k]]
    }

    if (verbose) {
      message(sprintf(
        "  [warm-start] %-30s — %s  [loglik = %.4f, iter = %d, converged = %s]",
        .label,
        basename(.latest),
        .loaded$loglik,
        .loaded$iterations,
        if (isTRUE(.loaded$converged)) "YES" else "NO"
      ))
    }

    out[[.label]] <- .base
    rm(.latest, .loaded, .base, .k)
  }

  out
}
