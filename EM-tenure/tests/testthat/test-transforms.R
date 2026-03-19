# ==============================================================================
# Tests for EM-tenure: transforms.R
# ==============================================================================

test_that("logit and inv_logit are inverses", {
  p_vals <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (p in p_vals) {
    expect_equal(inv_logit(logit(p)), p, tolerance = 1e-12)
  }
})

test_that("logit handles boundaries", {
  expect_equal(logit(0.5), 0)
  expect_true(logit(0.01) < 0)
  expect_true(logit(0.99) > 0)
})

test_that("CTMC link round-trips correctly", {
  theta_vals <- c(0.05, 0.1, 0.5, 0.9, 0.95)
  for (th in theta_vals) {
    lam <- ctmc_lambda_from_theta(th)
    th_back <- ctmc_theta_from_lambda(lam)
    expect_equal(th_back, th, tolerance = 1e-12)
  }
})

test_that("CTMC lambda is positive for theta in (0,1)", {
  expect_true(ctmc_lambda_from_theta(0.5) > 0)
  expect_true(ctmc_lambda_from_theta(0.9) > 0)
  expect_true(ctmc_lambda_from_theta(0.01) > 0)
})

test_that("CTMC link matches TeX formula: lambda = -log(1-theta)/3", {
  theta <- 0.85
  expected <- -log(1 - theta) / 3
  expect_equal(ctmc_lambda_from_theta(theta), expected, tolerance = 1e-12)
})

test_that("CTMC link handles boundary theta via .bound01 clamping", {
  # theta=0 and theta=1 are clamped by .bound01, so lambda stays finite
  expect_true(is.finite(ctmc_lambda_from_theta(0)))
  expect_true(is.finite(ctmc_lambda_from_theta(1)))
  # Near boundaries still produce sensible values
  expect_true(is.finite(ctmc_lambda_from_theta(0.9999)))
  # theta=0 (clamped to eps) should produce same result as eps directly
  expect_equal(ctmc_lambda_from_theta(0), ctmc_lambda_from_theta(1e-10))
})
