# ==============================================================================
# EM-baseline: Tests for transforms.R
# Created: 2026-05-05
# ==============================================================================

test_that("logit and inv_logit are exact inverses", {
  p_vals <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  expect_equal(inv_logit(logit(p_vals)), p_vals, tolerance = 1e-12)
})

test_that("inv_logit and logit are exact inverses (reverse direction)", {
  x_vals <- c(-5, -1, 0, 1, 5)
  expect_equal(logit(inv_logit(x_vals)), x_vals, tolerance = 1e-12)
})

test_that("logit(0.5) == 0", {
  expect_equal(logit(0.5), 0)
})

test_that("inv_logit(0) == 0.5", {
  expect_equal(inv_logit(0), 0.5)
})

test_that("logit is monotone increasing", {
  p <- seq(0.1, 0.9, by = 0.1)
  expect_true(all(diff(logit(p)) > 0))
})

test_that("logit boundary: p=0 -> -Inf, p=1 -> +Inf", {
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1),  Inf)
})

test_that("logit invalid: out-of-range inputs produce NaN or Inf", {
  # p < 0: log of negative -> NaN
  expect_true(is.nan(logit(-0.1)))
  # p > 1: log of negative (p/(1-p) < 0) -> NaN
  expect_true(is.nan(logit(1.1)))
})

test_that("inv_logit large positive -> near 1, large negative -> near 0", {
  expect_equal(inv_logit(100),  1, tolerance = 1e-10)
  expect_equal(inv_logit(-100), 0, tolerance = 1e-10)
})
