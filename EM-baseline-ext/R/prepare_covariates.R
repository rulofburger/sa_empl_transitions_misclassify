# ==============================================================================
# EM-baseline-ext: Build covariate design matrix for Extension I
# Created: 2026-05-06
#
# Constructs one of three fixed covariate specifications for the observable
# heterogeneity extension (EM baseline.tex Section 5). Uses wave-1 values
# only (time-invariant covariates assumption). Continuous variables are
# standardised to mean 0, sd 1 for numerical stability of the probit IRLS.
#
# Covariate sets:
#   Set 1 (parsimonious): intercept, age, age^2, education
#   Set 2 (demographics): Set 1 + race dummies + female
#   Set 3 (full):         Set 2 + contract type dummies
# ==============================================================================

#' Build the probit covariate design matrix
#'
#' Constructs an N × p design matrix for the covariate extension of the
#' baseline EM model (TeX Section~5). Wave-1 values are used for all
#' variables (time-invariant covariates assumption). Continuous predictors
#' are standardised to mean 0 / sd 1 on the provided sample.
#'
#' @param df Data frame containing the following columns (from
#'   \code{ingest_data_3waves_SA.R}):
#'   \code{age1}, \code{educ1} (required for all sets);
#'   \code{race1}, \code{female1} (required for sets 2 and 3);
#'   \code{contracttype1} (required for set 3).
#' @param covariate_set Integer scalar: \code{1L} (parsimonious),
#'   \code{2L} (demographics), or \code{3L} (full). Default \code{1L}.
#' @return A named list with:
#'   \describe{
#'     \item{X}{N × p numeric matrix with an intercept column as the first
#'       column. Column names describe each covariate.}
#'     \item{col_names}{Character vector of column names of \code{X}.}
#'     \item{center}{Named numeric vector: centering constants applied to
#'       each continuous column (0 for binary/intercept columns).}
#'     \item{scale}{Named numeric vector: scaling constants applied to each
#'       continuous column (1 for binary/intercept columns). Use
#'       \code{(x - center) / scale} to reproduce standardisation.}
#'   }
#' @details
#' \strong{Variable transformations}:
#' \itemize{
#'   \item \code{age1}: standardised to mean 0 / sd 1.
#'   \item \code{age1^2}: computed from raw (unstandardised) \code{age1},
#'     then standardised separately.
#'   \item \code{educ1}: standardised to mean 0 / sd 1.
#'   \item \code{race1}: one-hot encoded (first category dropped as reference).
#'   \item \code{female1}: included as-is (already binary 0/1).
#'   \item \code{contracttype1}: one-hot encoded (first category dropped; wave-1 value).
#' }
#' @examples
#' \dontrun{
#' source("scripts/ingest_data_3waves_SA.R")
#' result <- prepare_covariate_matrix(df_qlfs, covariate_set = 2L)
#' dim(result$X)
#' }
#' @export
prepare_covariate_matrix <- function(df, covariate_set = 1L) {
  if (!covariate_set %in% c(1L, 2L, 3L))
    stop(sprintf(
      "prepare_covariate_matrix: covariate_set must be 1, 2, or 3. Got: %s",
      covariate_set
    ))

  # Required columns per set
  required <- c("age1", "educ1")
  if (covariate_set >= 2L) required <- c(required, "race1", "female1")
  if (covariate_set >= 3L) required <- c(required, "contracttype1")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0L)
    stop(sprintf(
      "prepare_covariate_matrix: missing required columns for set %d: %s",
      covariate_set, paste(missing_cols, collapse = ", ")
    ))

  # Guard against NA in required columns
  for (col in required) {
    if (anyNA(df[[col]]))
      stop(sprintf(
        "prepare_covariate_matrix: column '%s' contains NA values. Drop incomplete rows before calling this function.",
        col
      ))
  }

  N <- nrow(df)

  # ---- Tracking vectors for centering / scaling ----------------------------
  center_vals <- numeric(0)
  scale_vals  <- numeric(0)

  # ---- Helper: standardise a numeric vector --------------------------------
  # Returns a named list: mu (center), sig (scale), val (standardised values)
  .std <- function(x, name) {
    mu  <- mean(x, na.rm = TRUE)
    sig <- sd(x, na.rm = TRUE)
    if (sig < 1e-10) {
      warning(sprintf(
        "prepare_covariate_matrix: column '%s' has near-zero variance (sd=%.2e); not standardised.",
        name, sig
      ))
      sig <- 1
    }
    list(mu = mu, sig = sig, val = (x - mu) / sig)
  }

  # ---- Helper: one-hot encode a categorical variable -----------------------
  # Drops the first (lowest-valued) category as reference.
  # Returns a named list: mat (N x k matrix), centers (named numeric), scales (named numeric)
  .onehot <- function(x, prefix) {
    lvls <- sort(unique(x[!is.na(x)]))
    if (length(lvls) < 2L) {
      warning(sprintf(
        "prepare_covariate_matrix: '%s' has fewer than 2 unique values; no dummy created.",
        prefix
      ))
      return(list(mat = matrix(nrow = length(x), ncol = 0L),
                  centers = numeric(0), scales = numeric(0)))
    }
    # Drop first level as reference category
    dummy_lvls <- lvls[-1L]
    mat <- matrix(0L, nrow = length(x), ncol = length(dummy_lvls))
    col_nms <- paste0(prefix, "_", dummy_lvls)
    colnames(mat) <- col_nms
    for (j in seq_along(dummy_lvls)) {
      mat[, j] <- as.integer(x == dummy_lvls[j])
    }
    # Dummy columns: center=0, scale=1
    centers <- setNames(rep(0, length(col_nms)), col_nms)
    scales  <- setNames(rep(1, length(col_nms)), col_nms)
    list(mat = mat, centers = centers, scales = scales)
  }

  # ---- Build column list ---------------------------------------------------
  col_parts <- list()

  # Intercept (always first)
  intercept_mat <- matrix(1, nrow = N, ncol = 1L, dimnames = list(NULL, "intercept"))
  center_vals["intercept"] <- 0
  scale_vals["intercept"]  <- 1
  col_parts[["intercept"]] <- intercept_mat

  # age (standardised)
  age_raw    <- df$age1
  age_res    <- .std(age_raw, "age")
  age_std    <- matrix(age_res$val, ncol = 1L, dimnames = list(NULL, "age"))
  center_vals["age"] <- age_res$mu
  scale_vals["age"]  <- age_res$sig
  col_parts[["age"]] <- age_std

  # age^2 (compute from raw, then standardise)
  age2_raw   <- age_raw^2
  age2_res   <- .std(age2_raw, "age_sq")
  age2_std   <- matrix(age2_res$val, ncol = 1L, dimnames = list(NULL, "age_sq"))
  center_vals["age_sq"] <- age2_res$mu
  scale_vals["age_sq"]  <- age2_res$sig
  col_parts[["age_sq"]] <- age2_std

  # education (standardised)
  educ_res   <- .std(df$educ1, "educ")
  educ_std   <- matrix(educ_res$val, ncol = 1L, dimnames = list(NULL, "educ"))
  center_vals["educ"] <- educ_res$mu
  scale_vals["educ"]  <- educ_res$sig
  col_parts[["educ"]] <- educ_std

  if (covariate_set >= 2L) {
    # race dummies (one-hot, reference = first race category)
    race_res <- .onehot(df$race1, "race")
    if (ncol(race_res$mat) > 0L) {
      col_parts[["race"]] <- race_res$mat
      center_vals         <- c(center_vals, race_res$centers)
      scale_vals          <- c(scale_vals,  race_res$scales)
    }

    # female (binary, no transformation)
    # Guard: factors give {1, 2} not {0, 1} — catch before silent bias
    if (is.factor(df$female1)) {
      stop("prepare_covariate_matrix: 'female1' is a factor. Convert to integer 0/1 before calling.")
    } else if (!all(df$female1[!is.na(df$female1)] %in% c(0, 1))) {
      stop("prepare_covariate_matrix: 'female1' must be binary (0/1 only).")
    }
    center_vals["female"] <- 0
    scale_vals["female"]  <- 1
    female_mat <- matrix(as.integer(df$female1), ncol = 1L,
                         dimnames = list(NULL, "female"))
    col_parts[["female"]] <- female_mat
  }

  if (covariate_set >= 3L) {
    # contract type dummies (one-hot, reference = first category)
    ct_res <- .onehot(df$contracttype1, "contracttype")
    if (ncol(ct_res$mat) > 0L) {
      col_parts[["contracttype"]] <- ct_res$mat
      center_vals                  <- c(center_vals, ct_res$centers)
      scale_vals                   <- c(scale_vals,  ct_res$scales)
    }
  }

  # ---- Assemble final matrix -----------------------------------------------
  X <- do.call(cbind, col_parts)

  list(
    X         = X,
    col_names = colnames(X),
    center    = center_vals,
    scale     = scale_vals
  )
}
