# ==============================================================================
# EM-AR2: Internal utility functions
# Created: 2026-05-05
# ==============================================================================
# Low-level helpers used across all EM-AR2 modules. Not exported for direct
# use. Essential for numerical stability and boundary enforcement.
# ==============================================================================

#' Clamp a value to the open interval (eps, 1-eps)
#'
#' Prevents probabilities from hitting exactly 0 or 1, which would cause
#' log(0) = -Inf or division by zero in the EM updates.
#'
#' @param x Numeric vector.
#' @param eps Boundary epsilon (default 1e-12).
#' @return Numeric vector with values clamped to (eps, 1-eps).
#' @keywords internal
.bound01 <- function(x, eps = 1e-12) {
  pmin(pmax(x, eps), 1 - eps)
}

#' Numerically stable log-sum-exp
#'
#' Computes log(sum(exp(x))) in a numerically stable way by subtracting the
#' maximum before exponentiating. Used in the E-step to normalise
#' log-responsibilities without overflow/underflow.
#'
#' @param x Numeric vector of log-values.
#' @return Scalar: log(sum(exp(x))).
#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}

#' Row-wise log-sum-exp for a matrix
#'
#' Applies .logsumexp() to each row of a matrix. Used in the E-step to
#' compute normalisation constants for each observation.
#'
#' @param mat Numeric matrix (N x H).
#' @return Length-N vector: log(sum(exp(mat[i,]))) for each row i.
#' @keywords internal
.row_logsumexp <- function(mat) {
  # Per-row maximum for numerical stability
  row_max <- apply(mat, 1, max)
  row_max + log(rowSums(exp(mat - row_max)))
}

#' Null-coalescing operator
#'
#' Returns lhs if non-NULL, otherwise rhs.
#'
#' @param lhs Left-hand side.
#' @param rhs Right-hand side (default value).
#' @return lhs if !is.null(lhs), else rhs.
#' @keywords internal
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs
