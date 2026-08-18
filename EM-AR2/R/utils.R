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

#' Collapse a four-wave panel to its 16 observed histories
#'
#' The likelihood and EM responsibilities are identical for observations with
#' the same four binary outcomes.  Collapsing therefore makes estimation on
#' the full QLFS panel fast while retaining squared weights and counts needed
#' for individual-level sandwich standard errors.
collapse_ar2_cells <- function(df) {
  required <- c("y1", "y2", "y3", "y4", "weight")
  missing <- setdiff(required, names(df))
  if (length(missing))
    stop("collapse_ar2_cells: missing columns: ", paste(missing, collapse = ", "))
  y <- as.matrix(df[, paste0("y", 1:4), drop = FALSE])
  if (anyNA(y) || !all(y %in% 0:1))
    stop("collapse_ar2_cells: y1-y4 must be non-missing binary values")
  if (anyNA(df$weight) || any(!is.finite(df$weight)) || any(df$weight <= 0))
    stop("collapse_ar2_cells: weights must be finite and strictly positive")

  index <- 1L + y[, 1L] + 2L * y[, 2L] + 4L * y[, 3L] + 8L * y[, 4L]
  histories <- as.data.frame(expand.grid(y1 = 0:1, y2 = 0:1,
                                         y3 = 0:1, y4 = 0:1))
  weight <- weight_sq <- numeric(16L)
  count <- tabulate(index, nbins = 16L)
  sw <- rowsum(df$weight, index, reorder = FALSE)
  sw2 <- rowsum(df$weight^2, index, reorder = FALSE)
  weight[as.integer(rownames(sw))] <- sw[, 1L]
  weight_sq[as.integer(rownames(sw2))] <- sw2[, 1L]
  cells <- cbind(histories, weight = weight,
                 weight_sq = weight_sq, count = count)
  cells <- cells[cells$count > 0L, , drop = FALSE]
  attr(cells, "n_obs") <- nrow(df)
  attr(cells, "weight_sum") <- sum(df$weight)
  cells
}
