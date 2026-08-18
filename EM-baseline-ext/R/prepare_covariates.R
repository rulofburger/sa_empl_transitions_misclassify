# ==============================================================================
# EM-baseline-ext: Build covariate design matrix for Extension I
# Created: 2026-05-06
#
# Constructs one of three fixed covariate specifications for the observable
# heterogeneity extension (EM baseline.tex Section 5). Uses wave-1 values
# for demographics. Contract type and sector use the origin-wave value for
# each transition. Continuous variables are
# standardised to mean 0, sd 1 for numerical stability of the probit IRLS.
#
# Covariate sets:
#   Set 1 (parsimonious): intercept, age, age^2, education
#   Set 2 (demographics): Set 1 + race dummies + female
#   Set 3 (exit-rich):    Set 2 + contract type dummies and an informal-sector
#                         indicator in the persistence equation only. Both are
#                         defined only in employment and excluded from entry.
# ==============================================================================

#' Attach transition-origin informal-sector indicators to a three-wave panel
#'
#' The estimation panel was built with person-matching method A, while sector
#' is retained only in the upstream long-format file. This helper performs a
#' validated many-to-one merge using household ID, method-A person ID, and
#' survey wave. Informal sector uses sector2 == 2, the QLFS definition that
#' allocates agriculture between formal and informal employment.
#'
#' @param df Three-wave panel containing hhnr, pnr, period1/period2, and y1/y2.
#' @param sector_long Upstream long-format QLFS data containing hhnr,
#'   pnr_methodA, wave, and sector2.
#' @return df with sector2_1/sector2_2 and informal_sector1/informal_sector2.
attach_transition_informal_sector <- function(df, sector_long) {
  panel_required <- c("hhnr", "pnr", "period1", "period2", "y1", "y2")
  source_required <- c("hhnr", "pnr_methodA", "wave", "sector2")
  missing_panel <- setdiff(panel_required, names(df))
  missing_source <- setdiff(source_required, names(sector_long))
  if (length(missing_panel))
    stop("attach_wave1_informal_sector: panel is missing: ",
         paste(missing_panel, collapse = ", "))
  if (length(missing_source))
    stop("attach_wave1_informal_sector: sector source is missing: ",
         paste(missing_source, collapse = ", "))

  relevant <- sector_long$wave %in% unique(c(df$period1, df$period2)) &
    !is.na(sector_long$pnr_methodA)
  src <- sector_long[relevant, source_required, drop = FALSE]
  source_key <- paste(src$hhnr, src$pnr_methodA, src$wave, sep = "|")
  if (anyDuplicated(source_key))
    stop("attach_wave1_informal_sector: upstream merge key is not unique")

  for (t in 1:2) {
    panel_key <- paste(df$hhnr, df$pnr, df[[paste0("period", t)]], sep = "|")
    idx <- match(panel_key, source_key)
    if (anyNA(idx))
      stop(sprintf(
        "attach_transition_informal_sector: wave %d: %d of %d rows did not match",
        t, sum(is.na(idx)), length(idx)
      ))
    sector2 <- as.numeric(unclass(src$sector2[idx]))
    employed <- df[[paste0("y", t)]] == 1L
    if (anyNA(sector2[employed]))
      stop(sprintf("attach_transition_informal_sector: sector2_%d missing for employed respondents", t))
    if (!all(sector2[employed] %in% c(1, 2, 4)))
      stop(sprintf("attach_transition_informal_sector: unexpected sector2_%d code", t))
    df[[paste0("sector2_", t)]] <- sector2
    df[[paste0("informal_sector", t)]] <-
      as.integer(employed & !is.na(sector2) & sector2 == 2)
  }
  df
}

# Backward-compatible wave-1 helper used by older scripts/tests.
attach_wave1_informal_sector <- function(df, sector_long) {
  if (!"period2" %in% names(df)) df$period2 <- df$period1
  if (!"y2" %in% names(df)) df$y2 <- df$y1
  attach_transition_informal_sector(df, sector_long)
}

#' Build the probit covariate design matrix
#'
#' Constructs one or two N × p design matrices for the covariate extension of
#' the baseline EM model (TeX Section~5). Wave-1 values are used for fixed
#' demographics; contract type and sector vary by transition. Continuous predictors
#' are standardised to mean 0 / sd 1 on the provided sample.
#'
#' @param df Data frame containing the following columns (from
#'   \code{ingest_data_3waves_SA.R}):
#'   \code{age1}, \code{educ1} (required for all sets);
#'   \code{race1}, \code{female1} (required for sets 2 and 3);
#'   \code{contracttype1}, \code{contracttype2}, \code{informal_sector1},
#'   \code{informal_sector2} (required for set 3).
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
#'   \item Contract type and informal sector use their transition-origin
#'     values: wave 1 for transition 1--2 and wave 2 for transition 2--3.
#'     They enter the persistence equation only.
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
  if (covariate_set >= 3L)
    required <- c(required, "contracttype1", "contracttype2",
                  "informal_sector1", "informal_sector2")
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
  contract_t2_mat <- NULL

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
    # Use the union of origin-wave categories so both transition matrices have
    # an identical coefficient layout.
    ct_res <- .onehot(c(df$contracttype1, df$contracttype2), "contracttype")
    if (ncol(ct_res$mat) > 0L) {
      col_parts[["contracttype"]] <- ct_res$mat[seq_len(N), , drop = FALSE]
      contract_t2_mat <- ct_res$mat[N + seq_len(N), , drop = FALSE]
      center_vals                  <- c(center_vals, ct_res$centers)
      scale_vals                   <- c(scale_vals,  ct_res$scales)
    }

    if (!all(df$informal_sector1 %in% c(0, 1)) ||
        !all(df$informal_sector2 %in% c(0, 1)))
      stop("prepare_covariate_matrix: informal-sector indicators must be binary (0/1 only).")
    informal_mat <- matrix(as.integer(df$informal_sector1), ncol = 1L,
                           dimnames = list(NULL, "informal_sector"))
    col_parts[["informal_sector"]] <- informal_mat
    center_vals["informal_sector"] <- 0
    scale_vals["informal_sector"] <- 1
  }

  # ---- Assemble final matrix -----------------------------------------------
  X <- do.call(cbind, col_parts)
  X_transition <- list(X12 = X, X23 = X)
  if (covariate_set >= 3L) {
    contract_cols <- grepl("^contracttype_", colnames(X))
    if (any(contract_cols))
      X_transition$X23[, contract_cols] <- contract_t2_mat
    X_transition$X23[, "informal_sector"] <- as.integer(df$informal_sector2)
  }

  # Contract type and informal sector are defined only for employed
  # respondents. Keep a common coefficient layout, but mark both inactive in
  # the entry equation. The EM driver fixes their beta0 coefficients at zero.
  entry_active <- !grepl("^(contracttype_|informal_sector$)", colnames(X))
  attr(X, "entry_active") <- entry_active
  attr(X_transition$X12, "entry_active") <- entry_active
  attr(X_transition$X23, "entry_active") <- entry_active

  list(
    X         = X,
    col_names = colnames(X),
    center    = center_vals,
    scale     = scale_vals,
    entry_active = entry_active,
    X_transition = X_transition
  )
}
