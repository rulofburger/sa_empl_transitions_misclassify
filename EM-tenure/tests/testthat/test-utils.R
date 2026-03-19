# ==============================================================================
# Tests for EM-tenure: utils.R
# ==============================================================================

test_that(".bound01 clamps to (eps, 1-eps)", {
  expect_equal(.bound01(0), 1e-10)
  expect_equal(.bound01(1), 1 - 1e-10)
  expect_equal(.bound01(0.5), 0.5)
  expect_equal(.bound01(-1), 1e-10)
  expect_equal(.bound01(2), 1 - 1e-10)
})

test_that(".pos enforces positivity", {
  expect_equal(.pos(0), 1e-12)
  expect_equal(.pos(-5), 1e-12)
  expect_equal(.pos(3), 3)
})

test_that(".logsumexp handles extreme values", {
  # Known result: log(exp(1000) + exp(1001)) = 1001 + log(1 + exp(-1))
  result <- .logsumexp(c(1000, 1001))
  expected <- 1001 + log1p(exp(-1))
  expect_equal(result, expected, tolerance = 1e-10)

  # Single element

  expect_equal(.logsumexp(5), 5)

  # All equal: log(n * exp(x)) = x + log(n)
  expect_equal(.logsumexp(rep(3, 4)), 3 + log(4), tolerance = 1e-10)
})

test_that("erfc matches 2*pnorm(-x*sqrt(2))", {
  x_vals <- c(-2, -1, 0, 0.5, 1, 2, 5)
  for (x in x_vals) {
    expect_equal(erfc(x), 2 * pnorm(-x * sqrt(2)), tolerance = 1e-10)
  }
})

test_that("%||% provides default for NULL", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
})
