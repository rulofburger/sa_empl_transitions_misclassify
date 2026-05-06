# ==============================================================================
# EM-baseline-ext: Tests for compute_inconsistencies()
# Created: 2026-05-06
# ==============================================================================

# Helper: build a minimal 1-row data frame with given age/educ values
.make_incons_df <- function(age1, age2, age3, educ1, educ2, educ3) {
  data.frame(
    age1  = age1,  age2  = age2,  age3  = age3,
    educ1 = educ1, educ2 = educ2, educ3 = educ3
  )
}

# ---- Output structure -------------------------------------------------------

test_that("compute_inconsistencies returns a data frame with 6 new columns", {
  df  <- .make_incons_df(25, 26, 27, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  new_cols <- c("Y_age_1", "Y_age_2", "Y_age_3", "Y_edu_1", "Y_edu_2", "Y_edu_3")
  expect_true(all(new_cols %in% names(out)))
  expect_equal(nrow(out), nrow(df))
})

test_that("compute_inconsistencies returns integer 0/1 columns", {
  df  <- .make_incons_df(25, 26, 27, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  for (col in c("Y_age_1", "Y_age_2", "Y_age_3",
                "Y_edu_1", "Y_edu_2", "Y_edu_3")) {
    expect_true(all(out[[col]] %in% c(0L, 1L)),
                label = paste(col, "must be 0 or 1"))
  }
})

test_that("compute_inconsistencies preserves original columns", {
  df  <- .make_incons_df(25, 26, 27, 3L, 3L, 3L)
  df$extra <- "keep"
  out <- compute_inconsistencies(df)
  expect_true("extra" %in% names(out))
  expect_equal(out$extra, df$extra)
})

# ---- Clean data: all indicators zero ---------------------------------------

test_that("clean data produces all-zero indicators (age constant)", {
  df  <- .make_incons_df(30, 30, 30, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 0L)
})

test_that("clean data produces all-zero indicators (age +1 each wave)", {
  df  <- .make_incons_df(30, 31, 32, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 0L)
})

test_that("clean data produces all-zero indicators (education constant)", {
  df  <- .make_incons_df(30, 31, 32, 4L, 4L, 4L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 0L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 0L)
})

test_that("clean data produces all-zero indicators (education +1 once)", {
  df  <- .make_incons_df(30, 31, 32, 3L, 4L, 4L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 0L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 0L)
})

# ---- Age inconsistency attribution -----------------------------------------

test_that("age inconsistency in gap 1-2 only: attributed to wave 1", {
  # age jump of 2 in gap 1-2, normal in gap 2-3 -> V12=1, V23=0 -> Y_age_1=1
  df  <- .make_incons_df(30, 32, 33, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 1L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 0L)
})

test_that("age inconsistency in gap 2-3 only: attributed to wave 3", {
  # normal gap 1-2, age decrease in gap 2-3 -> V12=0, V23=1 -> Y_age_3=1
  df  <- .make_incons_df(30, 31, 29, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 1L)
})

test_that("age inconsistency in both gaps: attributed to wave 2", {
  # V12=1, V23=1 -> Y_age_2=1 (single bad wave 2 causes both gaps)
  df  <- .make_incons_df(30, 35, 31, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 1L)
  expect_equal(out$Y_age_3, 0L)
})

test_that("age decrease in gap 1-2 only: attributed to wave 1", {
  # age drops from 31 to 29 in gap 1-2, then normal -> Y_age_1=1
  df  <- .make_incons_df(31, 29, 30, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 1L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 0L)
})

# ---- Education inconsistency attribution -----------------------------------

test_that("education decrease in gap 2-3 only: attributed to wave 3", {
  # educ drops from 4 to 3 in gap 2-3 -> V12=0, V23=1 -> Y_edu_3=1
  df  <- .make_incons_df(30, 31, 32, 3L, 4L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 0L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 1L)
})

test_that("education jump > 1 in gap 1-2 only: attributed to wave 1", {
  # educ jumps from 2 to 5 -> V12=1, V23=0 -> Y_edu_1=1
  df  <- .make_incons_df(30, 31, 32, 2L, 5L, 5L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 1L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 0L)
})

test_that("education inconsistency in both gaps: attributed to wave 2", {
  # educ: 3 -> 6 -> 3 (jump up then jump down) -> V12=1, V23=1 -> Y_edu_2=1
  df  <- .make_incons_df(30, 31, 32, 3L, 6L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 0L)
  expect_equal(out$Y_edu_2, 1L)
  expect_equal(out$Y_edu_3, 0L)
})

# ---- NA handling ------------------------------------------------------------

test_that("NA in age2 results in 0 indicators for age (conservative)", {
  df  <- .make_incons_df(30, NA, 31, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  # Both V12 and V23 involve age2: both pairwise gaps undefined -> all 0
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 0L)
})

test_that("NA in age1 results in 0 for V12, so Y_age_1 and Y_age_2 are 0", {
  # age1=NA: V12 undefined -> V12=0; gap 2-3 inconsistent -> V23=1 -> Y_age_3=1
  df  <- .make_incons_df(NA, 30, 28, 3L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, 0L)
  expect_equal(out$Y_age_2, 0L)
  expect_equal(out$Y_age_3, 1L)
})

test_that("NA in educ column results in 0 education indicators", {
  df  <- .make_incons_df(30, 31, 32, NA, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 0L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 0L)
})

# ---- Error handling ---------------------------------------------------------

test_that("compute_inconsistencies errors on missing required columns", {
  df_bad <- data.frame(age1 = 30, age2 = 31, age3 = 32)  # missing educ cols
  expect_error(compute_inconsistencies(df_bad), "missing required columns")
})

# ---- Multiple rows ----------------------------------------------------------

test_that("compute_inconsistencies works correctly on multiple rows", {
  df <- data.frame(
    age1  = c(25, 30, 40),
    age2  = c(26, 32, 41),  # row 2: age jump of 2 -> Y_age_1=1
    age3  = c(27, 33, 42),
    educ1 = c(3L, 2L, 4L),
    educ2 = c(3L, 2L, 4L),
    educ3 = c(3L, 2L, 4L)
  )
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_age_1, c(0L, 1L, 0L))
  expect_equal(out$Y_age_2, c(0L, 0L, 0L))
  expect_equal(out$Y_age_3, c(0L, 0L, 0L))
  # Education all clean
  expect_equal(out$Y_edu_1, c(0L, 0L, 0L))
})

test_that("education decrease in gap 1-2 only: attributed to wave 1", {
  # educ drops from 4 to 3 (decrease) -> V12_edu = 1; stays at 3 -> V23_edu = 0
  # Attribution: Y_edu_1 = 1 (error at wave 1)
  df  <- .make_incons_df(30, 31, 32, 4L, 3L, 3L)
  out <- compute_inconsistencies(df)
  expect_equal(out$Y_edu_1, 1L)
  expect_equal(out$Y_edu_2, 0L)
  expect_equal(out$Y_edu_3, 0L)
})

# ---- Integration: compute_inconsistencies -> e_step_inconsistency ----------

test_that("compute_inconsistencies output flows into e_step_inconsistency", {
  set.seed(55L)
  n    <- 50L
  age1 <- sample(25:45, n, replace = TRUE)
  edu1 <- sample(1L:5L, n, replace = TRUE)
  df_panel <- data.frame(
    y1 = rbinom(n, 1L, 0.6),
    y2 = rbinom(n, 1L, 0.6),
    y3 = rbinom(n, 1L, 0.6),
    weight = rep(1, n),
    age1   = age1,
    age2   = age1 + sample(0:1, n, replace = TRUE),
    age3   = age1 + sample(1:2, n, replace = TRUE),
    educ1  = edu1,
    educ2  = pmin(edu1 + sample(0:1, n, replace = TRUE), 5L),
    educ3  = pmin(edu1 + sample(0:2, n, replace = TRUE), 5L)
  )
  inc_df  <- compute_inconsistencies(df_panel)
  inc_mat <- as.matrix(inc_df[, c("Y_age_1", "Y_age_2", "Y_age_3",
                                   "Y_edu_1", "Y_edu_2", "Y_edu_3")])
  params  <- init_params_inconsistency()
  # Should run without error and return valid output
  out     <- e_step_inconsistency(df_panel, inc_mat, params)
  expect_true(is.finite(out$loglik))
  expect_equal(nrow(out$gamma), n)
  expect_equal(ncol(out$gamma), 8L)
  expect_equal(rowSums(out$gamma), rep(1, n), tolerance = 1e-10)
})
