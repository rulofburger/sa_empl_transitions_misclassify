# ==============================================================================
# EM-baseline-ext: Tests for prepare_covariate_matrix()
# Created: 2026-05-06
# ==============================================================================

# Helper: synthetic data frame mimicking df_qlfs columns
.make_cov_df <- function(n = 20L) {
  set.seed(42L)
  data.frame(
    age1         = sample(20:50, n, replace = TRUE),
    educ1        = sample(1:5,   n, replace = TRUE),
    race1        = sample(1:4,   n, replace = TRUE),
    female1      = sample(0:1,   n, replace = TRUE),
    contracttype1 = sample(1:3,  n, replace = TRUE)
  )
}

# ---- Column count per set --------------------------------------------------

test_that("set 1 returns intercept + age + age_sq + educ (p=4)", {
  df  <- .make_cov_df()
  out <- prepare_covariate_matrix(df, covariate_set = 1L)
  expect_equal(colnames(out$X)[1L], "intercept")
  expect_true(all(c("age", "age_sq", "educ") %in% colnames(out$X)))
  expect_equal(ncol(out$X), 4L)
})

test_that("set 2 includes race dummies and female on top of set 1", {
  df  <- .make_cov_df()
  out1 <- prepare_covariate_matrix(df, covariate_set = 1L)
  out2 <- prepare_covariate_matrix(df, covariate_set = 2L)
  # Set 2 must have at least set 1 columns + female + race dummies (3 levels - 1 = 3 dummies)
  expect_true(ncol(out2$X) > ncol(out1$X))
  expect_true("female" %in% colnames(out2$X))
  expect_true(any(grepl("^race_", colnames(out2$X))))
})

test_that("set 3 includes contract type dummies on top of set 2", {
  df  <- .make_cov_df()
  out2 <- prepare_covariate_matrix(df, covariate_set = 2L)
  out3 <- prepare_covariate_matrix(df, covariate_set = 3L)
  expect_true(ncol(out3$X) > ncol(out2$X))
  expect_true(any(grepl("^contracttype_", colnames(out3$X))))
  contract_cols <- grepl("^contracttype_", colnames(out3$X))
  expect_true(all(!out3$entry_active[contract_cols]))
  expect_true(all(out3$entry_active[!contract_cols]))
})

test_that("set 3 contract coefficients remain zero in the entry equation", {
  set.seed(12L)
  n <- 250L
  cov_df <- .make_cov_df(n)
  prep <- prepare_covariate_matrix(cov_df, covariate_set = 3L)
  panel <- data.frame(
    y1 = rbinom(n, 1L, 0.5),
    y2 = rbinom(n, 1L, 0.5),
    y3 = rbinom(n, 1L, 0.5),
    weight = rep(1, n)
  )
  fit <- em_fit_covariates(panel, prep$X, model_type = "symmetric",
                           max_iter = 100L, verbose = 0L)
  contract_cols <- grepl("^contracttype_", colnames(prep$X))
  expect_equal(fit$params$beta0[contract_cols], rep(0, sum(contract_cols)))
})

# ---- Intercept -------------------------------------------------------------

test_that("first column is all-ones intercept", {
  df  <- .make_cov_df()
  out <- prepare_covariate_matrix(df, covariate_set = 1L)
  expect_true(all(out$X[, "intercept"] == 1))
})

# ---- Standardisation -------------------------------------------------------

test_that("continuous columns (age, age_sq, educ) are standardised (mean ~0, sd ~1)", {
  df  <- .make_cov_df(n = 200L)
  out <- prepare_covariate_matrix(df, covariate_set = 1L)
  for (col in c("age", "age_sq", "educ")) {
    expect_equal(mean(out$X[, col]), 0, tolerance = 1e-10,
                 label = paste(col, "mean should be 0"))
    expect_equal(sd(out$X[, col]),   1, tolerance = 1e-10,
                 label = paste(col, "sd should be 1"))
  }
})

test_that("intercept column has center=0 and scale=1", {
  df  <- .make_cov_df()
  out <- prepare_covariate_matrix(df, covariate_set = 1L)
  expect_equal(out$center[["intercept"]], 0)
  expect_equal(out$scale[["intercept"]],  1)
})

test_that("center and scale vectors have same names as X columns", {
  df  <- .make_cov_df()
  out <- prepare_covariate_matrix(df, covariate_set = 2L)
  expect_setequal(names(out$center), colnames(out$X))
  expect_setequal(names(out$scale),  colnames(out$X))
})

# ---- Dummy encoding --------------------------------------------------------

test_that("race dummies are binary 0/1", {
  df   <- .make_cov_df()
  out  <- prepare_covariate_matrix(df, covariate_set = 2L)
  race_cols <- grep("^race_", colnames(out$X), value = TRUE)
  for (col in race_cols) {
    expect_true(all(out$X[, col] %in% c(0L, 1L)),
                label = paste(col, "must be 0 or 1"))
  }
})

test_that("contract type dummies are binary 0/1", {
  df   <- .make_cov_df()
  out  <- prepare_covariate_matrix(df, covariate_set = 3L)
  ct_cols <- grep("^contracttype_", colnames(out$X), value = TRUE)
  for (col in ct_cols) {
    expect_true(all(out$X[, col] %in% c(0L, 1L)),
                label = paste(col, "must be 0 or 1"))
  }
})

test_that("race dummies: each row has at most one active dummy (one-hot)", {
  df  <- .make_cov_df(n = 100L)
  out <- prepare_covariate_matrix(df, covariate_set = 2L)
  race_cols <- grep("^race_", colnames(out$X), value = TRUE)
  row_sums  <- rowSums(out$X[, race_cols, drop = FALSE])
  # Rows with reference category have sum 0; all others have sum 1
  expect_true(all(row_sums %in% c(0, 1)))
})

# ---- Dimensions and nrow ---------------------------------------------------

test_that("number of rows in X equals nrow(df)", {
  df  <- .make_cov_df(n = 30L)
  out <- prepare_covariate_matrix(df, covariate_set = 3L)
  expect_equal(nrow(out$X), 30L)
})

# ---- Return value structure ------------------------------------------------

test_that("return value has components X, col_names, center, scale", {
  df  <- .make_cov_df()
  out <- prepare_covariate_matrix(df, covariate_set = 1L)
  expect_true(is.matrix(out$X))
  expect_type(out$col_names, "character")
  expect_type(out$center, "double")
  expect_type(out$scale, "double")
  expect_equal(out$col_names, colnames(out$X))
})

# ---- Error handling --------------------------------------------------------

test_that("invalid covariate_set raises an error", {
  df <- .make_cov_df()
  expect_error(prepare_covariate_matrix(df, covariate_set = 0L), "covariate_set must be")
  expect_error(prepare_covariate_matrix(df, covariate_set = 4L), "covariate_set must be")
})

test_that("missing required columns raise an error", {
  df_bad <- data.frame(age1 = 1:5)   # missing educ1 and others
  expect_error(prepare_covariate_matrix(df_bad, covariate_set = 1L), "missing required columns")
})

test_that("NA in required column raises an error", {
  df      <- .make_cov_df()
  df$age1[3] <- NA
  expect_error(prepare_covariate_matrix(df, covariate_set = 1L), "contains NA")
})

test_that("near-zero variance column triggers warning", {
  df_const        <- .make_cov_df()
  df_const$age1   <- rep(30L, nrow(df_const))  # constant age -> zero std dev
  expect_warning(prepare_covariate_matrix(df_const, covariate_set = 1L),
                 "near-zero variance")
})

test_that("race dummy mapping: race=1 is reference (all dummies zero)", {
  # Verify that the first race category is the reference level (no dummy column)
  # and that each other category maps to exactly one active dummy.
  df          <- .make_cov_df(n = 60L)
  df$race1    <- rep(c(1L, 2L, 3L, 4L), 15L)   # balanced
  out         <- prepare_covariate_matrix(df, covariate_set = 2L)
  race_cols   <- grep("^race_", colnames(out$X), value = TRUE)
  # Rows where race1 == 1 should have all race dummies = 0 (reference)
  ref_rows    <- which(df$race1 == 1L)
  expect_true(all(out$X[ref_rows, race_cols] == 0))
  # Each non-reference category should activate exactly one dummy
  for (r in c(2L, 3L, 4L)) {
    rows_r <- which(df$race1 == r)
    active  <- rowSums(out$X[rows_r, race_cols, drop = FALSE])
    expect_true(all(active == 1L),
                label = sprintf("race=%d rows must have exactly 1 active race dummy", r))
  }
})
