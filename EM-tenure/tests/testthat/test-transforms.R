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

test_that("CTMC transition link round-trips correctly", {
  # theta0-style: transition probability (exit from state 0)
  theta_vals <- c(0.05, 0.1, 0.5, 0.9, 0.95)
  for (th in theta_vals) {
    lam <- ctmc_lambda_from_transition(th)
    th_back <- ctmc_transition_from_lambda(lam)
    expect_equal(th_back, th, tolerance = 1e-12)
  }
})

test_that("CTMC persistence link round-trips correctly", {
  # theta1-style: persistence probability (stay in state 1)
  theta_vals <- c(0.05, 0.1, 0.5, 0.9, 0.95)
  for (th in theta_vals) {
    lam <- ctmc_lambda_from_persistence(th)
    th_back <- ctmc_persistence_from_lambda(lam)
    expect_equal(th_back, th, tolerance = 1e-12)
  }
})

test_that("CTMC lambdas are positive for theta in (0,1)", {
  expect_true(ctmc_lambda_from_transition(0.5) > 0)
  expect_true(ctmc_lambda_from_transition(0.9) > 0)
  expect_true(ctmc_lambda_from_transition(0.01) > 0)
  expect_true(ctmc_lambda_from_persistence(0.5) > 0)
  expect_true(ctmc_lambda_from_persistence(0.9) > 0)
  expect_true(ctmc_lambda_from_persistence(0.01) > 0)
})

test_that("CTMC transition matches formula: lambda = -log(1-theta)/0.25", {
  theta <- 0.10
  expected <- -log(1 - theta) / 0.25
  expect_equal(ctmc_lambda_from_transition(theta), expected, tolerance = 1e-12)
})

test_that("CTMC persistence matches formula: lambda = -log(theta)/0.25", {
  theta <- 0.96
  expected <- -log(theta) / 0.25
  expect_equal(ctmc_lambda_from_persistence(theta), expected, tolerance = 1e-12)
})

test_that("CTMC gives correct expected tenure for theta1=0.96", {
  # theta1 = 0.96 → quarterly exit rate 4% → E[tenure] ≈ 6.1 years
  lambda_g <- ctmc_lambda_from_persistence(0.96)
  expected_tenure <- 1 / lambda_g
  expect_equal(expected_tenure, 1 / (-log(0.96) / 0.25), tolerance = 1e-6)
  expect_true(expected_tenure > 5.5 && expected_tenure < 7.0)
})

test_that("CTMC persistence != CTMC transition for same theta", {

  # For theta = 0.96:
  #   persistence: lambda = -log(0.96)/0.25 ≈ 0.163
  #   transition:  lambda = -log(1-0.96)/0.25 = -log(0.04)/0.25 ≈ 12.88
  # These should be very different
  lam_pers <- ctmc_lambda_from_persistence(0.96)
  lam_tran <- ctmc_lambda_from_transition(0.96)
  expect_true(lam_tran > 10 * lam_pers)
})

test_that("CTMC link handles boundary theta via .bound01 clamping", {
  # theta=0 and theta=1 are clamped by .bound01, so lambda stays finite
  expect_true(is.finite(ctmc_lambda_from_transition(0)))
  expect_true(is.finite(ctmc_lambda_from_transition(1)))
  expect_true(is.finite(ctmc_lambda_from_persistence(0)))
  expect_true(is.finite(ctmc_lambda_from_persistence(1)))
  # Near boundaries still produce sensible values
  expect_true(is.finite(ctmc_lambda_from_transition(0.9999)))
  expect_true(is.finite(ctmc_lambda_from_persistence(0.0001)))
})
