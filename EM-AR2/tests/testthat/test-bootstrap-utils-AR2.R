# ==============================================================================
# EM-AR2 tests: bootstrap_utils_AR2.R
# Created: 2026-05-07
#
# Tests:
#   .find_latest_fit()     — file discovery by timestamp suffix
#   .flag_fit_ar2()        — quality flag assignment
#   .ar2_fit_args()        — model type → em_fit_ar2 argument mapping
#
# bootstrap_one_ar2() is an integration test that requires the EM algorithm
# to run; it is skipped unless actual data and model fits are available.
# ==============================================================================

# ---------------------------------------------------------------------------
# .find_latest_fit() tests
# ---------------------------------------------------------------------------

test_that(".find_latest_fit returns the most recent file", {
  tmp <- withr::local_tempdir()
  # Create fake files with different timestamp suffixes
  file.create(file.path(tmp, "em_ar2_sym_20260501_100000.rds"))
  file.create(file.path(tmp, "em_ar2_sym_20260505_232625.rds"))
  file.create(file.path(tmp, "em_ar2_sym_20260430_090000.rds"))

  result <- .find_latest_fit("sym", tmp)
  expect_equal(basename(result), "em_ar2_sym_20260505_232625.rds")
})

test_that(".find_latest_fit errors when no matching file exists", {
  tmp <- withr::local_tempdir()
  expect_error(.find_latest_fit("nopi", tmp), regexp = "no files matching")
})

test_that(".find_latest_fit errors when directory does not exist", {
  expect_error(.find_latest_fit("sym", "/no/such/dir"), regexp = "no files matching|cannot open")
})

test_that(".find_latest_fit errors when directory is empty", {
  tmp <- withr::local_tempdir()
  # Directory exists but has no matching .rds files
  expect_error(.find_latest_fit("sym", tmp), regexp = "no files matching")
})

test_that(".find_latest_fit returns consistent result for same-basename files (tiebreaking)", {
  tmp <- withr::local_tempdir()
  # Same timestamp basename: tiebreaking via file.mtime should be consistent
  f1 <- file.path(tmp, "em_ar2_sym_20260505_232625.rds")
  f2 <- file.path(tmp, "em_ar2_sym_20260505_232625_v2.rds")
  file.create(f1)
  file.create(f2)
  result1 <- .find_latest_fit("sym", tmp)
  result2 <- .find_latest_fit("sym", tmp)
  # Must be deterministic: same file returned on repeated calls
  expect_equal(result1, result2)
})

test_that(".find_latest_fit only returns the key-specific file", {
  tmp <- withr::local_tempdir()
  file.create(file.path(tmp, "em_ar2_sym_20260505_232625.rds"))
  file.create(file.path(tmp, "em_ar2_asym_20260506_100000.rds"))

  result_sym <- .find_latest_fit("sym", tmp)
  expect_match(basename(result_sym), "^em_ar2_sym_")

  result_asym <- .find_latest_fit("asym", tmp)
  expect_match(basename(result_asym), "^em_ar2_asym_")
})

# ---------------------------------------------------------------------------
# .flag_fit_ar2() tests
# ---------------------------------------------------------------------------

test_that(".flag_fit_ar2 returns 'ok' for converged fit near point LL", {
  fit <- list(converged = TRUE, loglik = -1000)
  expect_equal(.flag_fit_ar2(fit, -1000), "ok")
})

test_that(".flag_fit_ar2 returns 'no_converge' when converged is FALSE", {
  fit <- list(converged = FALSE, loglik = -1000)
  expect_equal(.flag_fit_ar2(fit, -1000), "no_converge")
})

test_that(".flag_fit_ar2 returns 'error' when loglik is NA", {
  fit <- list(converged = TRUE, loglik = NA_real_)
  expect_equal(.flag_fit_ar2(fit, -1000), "error")
})

test_that(".flag_fit_ar2 returns 'low_loglik' when LL drops too much", {
  fit <- list(converged = TRUE, loglik = -1100)
  # point LL = -1000, threshold = 50 nats → flag if loglik < -1050
  expect_equal(.flag_fit_ar2(fit, -1000), "low_loglik")
})

test_that(".flag_fit_ar2 skips LL check when point_loglik is NA", {
  fit <- list(converged = TRUE, loglik = -9999)
  expect_equal(.flag_fit_ar2(fit, NA_real_), "ok")
})

test_that(".flag_fit_ar2 returns 'ok' when LL drop is exactly at threshold", {
  fit <- list(converged = TRUE, loglik = -1050)
  # -1050 is exactly at the boundary; not strictly less than -1050
  expect_equal(.flag_fit_ar2(fit, -1000), "ok")
})

test_that(".flag_fit_ar2 returns 'low_loglik' when LL drop strictly exceeds threshold", {
  fit <- list(converged = TRUE, loglik = -1050.1)
  # strictly below -1050 => flagged
  expect_equal(.flag_fit_ar2(fit, -1000), "low_loglik")
})

test_that(".flag_fit_ar2 returns 'ok' when LL drop is just inside threshold", {
  fit <- list(converged = TRUE, loglik = -1049.9)
  expect_equal(.flag_fit_ar2(fit, -1000), "ok")
})

# ---------------------------------------------------------------------------
# .ar2_fit_args() tests
# ---------------------------------------------------------------------------

test_that(".ar2_fit_args returns correct args for 'none'", {
  args <- .ar2_fit_args("none")
  expect_false(args$estimate_pi)
  expect_equal(args$fixed_pi, 0)
  expect_false(args$asymmetric)
})

test_that(".ar2_fit_args returns correct args for 'symmetric'", {
  args <- .ar2_fit_args("symmetric")
  expect_true(args$estimate_pi)
  expect_false(args$asymmetric)
})

test_that(".ar2_fit_args returns correct args for 'asymmetric'", {
  args <- .ar2_fit_args("asymmetric")
  expect_true(args$estimate_pi)
  expect_true(args$asymmetric)
})

test_that(".ar2_fit_args errors on unknown model_type", {
  expect_error(.ar2_fit_args("unknown"), regexp = "unknown model_type")
})
