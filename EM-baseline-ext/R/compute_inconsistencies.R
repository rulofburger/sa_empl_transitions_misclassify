# ==============================================================================
# EM-baseline-ext: Compute wave-attributed inconsistency indicators
# Created: 2026-05-06
#
# Computes wave-specific inconsistency dummies for age and education
# following the wave-attribution rules in EM baseline.tex Section 8
# (Extension IV), Eq (27).
#
# Age rule: between consecutive quarterly waves, age should remain constant
# or increase by exactly 1 year.
#
# Education rule: education level should be non-decreasing across waves and
# cannot jump by more than one category per quarter.
#
# Wave-attribution rule (for both variables):
#   V_12 = 1{ pairwise gap 1-2 is inconsistent }
#   V_23 = 1{ pairwise gap 2-3 is inconsistent }
#   Y_1  = V_12 * (1 - V_23)   <- only gap 1-2 inconsistent: error at wave 1
#   Y_2  = V_12 * V_23          <- both gaps inconsistent: error at wave 2
#   Y_3  = (1 - V_12) * V_23   <- only gap 2-3 inconsistent: error at wave 3
#
# If the original value at either relevant wave is NA, the pairwise indicator
# is set to 0 (conservative: we cannot confirm an inconsistency).
# ==============================================================================

#' Compute wave-attributed age and education inconsistency indicators
#'
#' Given a panel data frame with columns \code{age1}--\code{age3} and
#' \code{educ1}--\code{educ3}, computes six wave-attributed inconsistency
#' dummies following TeX Eq~(27).
#'
#' \strong{Age rule}: between consecutive quarterly waves, age must remain
#' constant or increase by exactly 1 (\code{age_{t+1} - age_t ∈ {0, 1}}).
#'
#' \strong{Education rule}: education must be non-decreasing and cannot jump
#' by more than one category (\code{0 ≤ educ_{t+1} - educ_t ≤ 1}).
#'
#' \strong{Wave-attribution}: a pairwise inconsistency in gap \eqn{t,t+1}
#' but not \eqn{t+1,t+2} is attributed to wave \eqn{t}; a simultaneous
#' inconsistency in both gaps is attributed to wave \eqn{t+1} (single
#' erroneous wave generates both pairwise discrepancies); an inconsistency
#' only in gap \eqn{t+1,t+2} is attributed to wave \eqn{t+2}.
#'
#' @param df Data frame with numeric columns \code{age1}, \code{age2},
#'   \code{age3}, \code{educ1}, \code{educ2}, \code{educ3}. Other columns
#'   are passed through unchanged.
#' @return The input data frame with six additional integer (0/1) columns:
#'   \code{Y_age_1}, \code{Y_age_2}, \code{Y_age_3},
#'   \code{Y_edu_1}, \code{Y_edu_2}, \code{Y_edu_3}.
#'   NA values in the original age/education columns produce a 0 indicator
#'   (conservative: inconsistency not confirmed).
#' @examples
#' df <- data.frame(
#'   age1 = c(25, 30, 40), age2 = c(26, 32, 41), age3 = c(27, 33, 42),
#'   educ1 = c(3L, 2L, 4L), educ2 = c(3L, 2L, 4L), educ3 = c(3L, 2L, 4L)
#' )
#' out <- compute_inconsistencies(df)
#' # Age gap 30->32 (jump of 2) is inconsistent at waves 1-2.
#' # Age gap 32->33 (normal) is consistent at waves 2-3.
#' # Attribution: Y_age_1 = 1 for that row.
#' @export
compute_inconsistencies <- function(df) {
  required_cols <- c("age1", "age2", "age3", "educ1", "educ2", "educ3")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L)
    stop(sprintf(
      "compute_inconsistencies: missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))

  # Type check: arithmetic on ordered factors throws cryptic errors
  for (col in required_cols) {
    if (is.logical(df[[col]]) && all(is.na(df[[col]])))
      df[[col]] <- as.numeric(df[[col]])
    if (!is.numeric(df[[col]]) && !is.integer(df[[col]]))
      stop(sprintf(
        "compute_inconsistencies: column '%s' must be numeric or integer (got %s). Convert factors before calling.",
        col, class(df[[col]])[1L]
      ))
  }

  # ---- Age inconsistency ------------------------------------------------
  # Pairwise: age should change by 0 or +1 between consecutive waves
  # NA in either wave → treat pairwise gap as consistent (conservative)
  d12_age <- df$age2 - df$age1
  d23_age <- df$age3 - df$age2

  # Use tolerance-based check rather than %in% c(0L,1L) to handle floating-point
  # age differences gracefully (e.g. 1.0000001 from double arithmetic).
  .age_ok <- function(d) abs(d - 0) < 0.01 | abs(d - 1) < 0.01
  V12_age <- as.integer(!is.na(d12_age) & !.age_ok(d12_age))
  V23_age <- as.integer(!is.na(d23_age) & !.age_ok(d23_age))

  df$Y_age_1 <- V12_age * (1L - V23_age)
  df$Y_age_2 <- V12_age * V23_age
  df$Y_age_3 <- (1L - V12_age) * V23_age

  # ---- Education inconsistency ------------------------------------------
  # Pairwise: educ must be non-decreasing and cannot jump by more than 1
  d12_edu <- df$educ2 - df$educ1
  d23_edu <- df$educ3 - df$educ2

  # Education: non-decreasing with at most +1 per wave
  .edu_ok <- function(d) abs(d - 0) < 0.01 | abs(d - 1) < 0.01
  V12_edu <- as.integer(!is.na(d12_edu) & !.edu_ok(d12_edu))
  V23_edu <- as.integer(!is.na(d23_edu) & !.edu_ok(d23_edu))

  df$Y_edu_1 <- V12_edu * (1L - V23_edu)
  df$Y_edu_2 <- V12_edu * V23_edu
  df$Y_edu_3 <- (1L - V12_edu) * V23_edu

  df
}

#' Compute the magnitude of wave-attributed inconsistencies
#'
#' Magnitude is the distance of a consecutive-wave change from the admissible
#' interval [0, 1].  When both adjacent gaps are inconsistent, their magnitudes
#' are summed and attributed to the middle wave, matching the indicator rule.
#' Missing gaps contribute zero.
#' @export
compute_inconsistency_extent <- function(df) {
  df <- compute_inconsistencies(df)
  distance <- function(d) ifelse(is.na(d), 0, ifelse(d < 0, -d,
                                                     ifelse(d > 1, d - 1, 0)))
  attribute_extent <- function(v12, v23, s12, s23) cbind(
    s12 * v12 * (1L - v23),
    (s12 + s23) * v12 * v23,
    s23 * (1L - v12) * v23
  )
  d12_age <- df$age2 - df$age1; d23_age <- df$age3 - df$age2
  d12_edu <- df$educ2 - df$educ1; d23_edu <- df$educ3 - df$educ2
  va12 <- as.integer(df$Y_age_1 == 1L | df$Y_age_2 == 1L)
  va23 <- as.integer(df$Y_age_2 == 1L | df$Y_age_3 == 1L)
  ve12 <- as.integer(df$Y_edu_1 == 1L | df$Y_edu_2 == 1L)
  ve23 <- as.integer(df$Y_edu_2 == 1L | df$Y_edu_3 == 1L)
  age_extent <- attribute_extent(va12, va23, distance(d12_age), distance(d23_age))
  edu_extent <- attribute_extent(ve12, ve23, distance(d12_edu), distance(d23_edu))
  df[paste0("extent_age_", 1:3)] <- age_extent
  df[paste0("extent_edu_", 1:3)] <- edu_extent
  df
}
