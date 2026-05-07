# ==============================================================================
# EM-baseline: Tests for implied_quantities.R
# Created: 2026-05-06
# ==============================================================================

# ---------------------------------------------------------------------------
# implied_baseline()
# ---------------------------------------------------------------------------

test_that("implied_baseline returns correct names for symmetric model", {
  p   <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
  out <- implied_baseline(p, "symmetric")
  expect_named(out, c("entry_rate", "exit_rate", "employment_rate",
                       "pi", "pi0", "pi1"))
})

test_that("implied_baseline entry_rate equals theta0", {
  p   <- list(theta0 = 0.12, theta1 = 0.85, pi = 0.05)
  out <- implied_baseline(p, "symmetric")
  expect_equal(out$entry_rate, 0.12)
})

test_that("implied_baseline exit_rate equals 1 - theta1", {
  p   <- list(theta0 = 0.12, theta1 = 0.85, pi = 0.05)
  out <- implied_baseline(p, "symmetric")
  expect_equal(out$exit_rate, 0.15)
})

test_that("implied_baseline employment_rate uses stationarity formula", {
  theta0 <- 0.10
  theta1 <- 0.90
  expected_alpha <- theta0 / (theta0 + 1 - theta1)
  p   <- list(theta0 = theta0, theta1 = theta1, pi = 0.05)
  out <- implied_baseline(p, "symmetric")
  expect_equal(out$employment_rate, expected_alpha, tolerance = 1e-12)
})

test_that("implied_baseline pi is correct for symmetric model", {
  p   <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.07)
  out <- implied_baseline(p, "symmetric")
  expect_equal(out$pi, 0.07)
  expect_true(is.na(out$pi0))
  expect_true(is.na(out$pi1))
})

test_that("implied_baseline pi0 and pi1 correct for asymmetric model", {
  p   <- list(theta0 = 0.10, theta1 = 0.90, pi0 = 0.03, pi1 = 0.07)
  out <- implied_baseline(p, "asymmetric")
  expect_equal(out$pi0, 0.03)
  expect_equal(out$pi1, 0.07)
  expect_true(is.na(out$pi))
})

test_that("implied_baseline all error rates NA for no-error model", {
  p   <- list(theta0 = 0.10, theta1 = 0.90)
  out <- implied_baseline(p, "none")
  expect_true(is.na(out$pi))
  expect_true(is.na(out$pi0))
  expect_true(is.na(out$pi1))
})

test_that("implied_baseline employment_rate in (0, 1) for valid params", {
  p   <- list(theta0 = 0.15, theta1 = 0.80, pi = 0.04)
  out <- implied_baseline(p, "symmetric")
  expect_gt(out$employment_rate, 0)
  expect_lt(out$employment_rate, 1)
})

test_that("implied_baseline uses stationarity formula even when alpha present", {
  # When params also contains a free alpha (from stationary = FALSE models),
  # implied_baseline should still compute employment_rate from theta0/theta1.
  theta0 <- 0.10
  theta1 <- 0.90
  alpha_free <- 0.60   # differs from stationarity value
  p <- list(theta0 = theta0, theta1 = theta1, alpha = alpha_free, pi = 0.05)
  out <- implied_baseline(p, "symmetric")
  expected_alpha <- theta0 / (theta0 + 1 - theta1)
  expect_equal(out$employment_rate, expected_alpha, tolerance = 1e-12)
})

test_that("implied_baseline errors on missing theta0", {
  p <- list(theta1 = 0.90, pi = 0.05)
  expect_error(implied_baseline(p, "symmetric"))
})

test_that("implied_baseline errors on invalid model_type", {
  p <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
  expect_error(implied_baseline(p, "wrong"))
})
