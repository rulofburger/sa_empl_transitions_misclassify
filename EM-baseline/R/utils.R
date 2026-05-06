# ==============================================================================
# EM-baseline: Shared utilities and private helpers
# Created: 2026-05-05
#
# Private helpers used throughout the EM-baseline module. All internal
# functions are prefixed with '.' and tagged @keywords internal.
# ==============================================================================

#' Clamp a probability to a safe open interval
#'
#' Prevents log(0) and division-by-zero by keeping probabilities strictly
#' inside (eps, 1 - eps).
#'
#' @param p Numeric vector of probabilities.
#' @param eps Boundary tolerance (default 1e-10).
#' @return Numeric vector with values in [eps, 1 - eps].
#' @keywords internal
.bound01 <- function(p, eps = 1e-10) {
  pmin(1 - eps, pmax(eps, p))
}

#' Numerically stable log-sum-exp
#'
#' Computes log(sum(exp(x))) without overflow or underflow by shifting by
#' max(x) before exponentiating.
#'
#' @param x Numeric vector.
#' @return Scalar: log(sum(exp(x))).
#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Null-coalescing operator
#'
#' Returns \code{a} if it is not NULL, otherwise returns \code{b}.
#'
#' @param a First value.
#' @param b Fallback value.
#' @return \code{a} if \code{!is.null(a)}, else \code{b}.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Validate that a data frame has the required panel columns
#'
#' Stops with an informative error if any of y1, y2, y3, weight are absent.
#'
#' @param df Data frame to validate.
#' @return Invisible TRUE.
#' @keywords internal
.validate_panel_df <- function(df) {
  required_cols <- c("y1", "y2", "y3", "weight")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("df is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  invisible(TRUE)
}

#' Validate model_type argument
#'
#' Stops with an informative error if model_type is not one of the three
#' recognised values.
#'
#' @param model_type Character string to validate.
#' @return Invisible TRUE.
#' @keywords internal
.validate_model_type <- function(model_type) {
  valid_types <- c("symmetric", "asymmetric", "none")
  if (!model_type %in% valid_types) {
    stop(sprintf("model_type must be one of: %s",
                 paste(valid_types, collapse = ", ")))
  }
  invisible(TRUE)
}
