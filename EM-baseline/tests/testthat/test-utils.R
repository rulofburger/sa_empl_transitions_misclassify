# ==============================================================================
# EM-baseline: Tests for utils.R
# Created: 2026-05-05
# ==============================================================================

test_that(".bound01 clamps values to (eps, 1-eps)", {
  eps <- 1e-10
  expect_equal(.bound01(0.5), 0.5)
  expect_equal(.bound01(0),   eps)
  expect_equal(.bound01(1),   1 - eps)
  expect_equal(.bound01(-1),  eps)
  expect_equal(.bound01(2),   1 - eps)
  # Vector input
  result <- .bound01(c(-1, 0.3, 2))
  expect_equal(result, c(eps, 0.3, 1 - eps))
})

test_that(".bound01 respects custom eps", {
  expect_equal(.bound01(0, eps = 1e-4), 1e-4)
  expect_equal(.bound01(1, eps = 1e-4), 1 - 1e-4)
})

test_that(".logsumexp matches naive computation for small vectors", {
  x <- c(1, 2, 3)
  expect_equal(.logsumexp(x), log(sum(exp(x))), tolerance = 1e-12)
})

test_that(".logsumexp handles large values without overflow", {
  x <- c(700, 701, 702)
  # naive exp(700) would overflow; .logsumexp should not
  result <- .logsumexp(x)
  expect_true(is.finite(result))
  expect_equal(result, log(exp(0) + exp(1) + exp(2)) + 700, tolerance = 1e-10)
})

test_that(".logsumexp handles single element", {
  expect_equal(.logsumexp(5), 5)
})

test_that("%||% returns first value when not NULL", {
  expect_equal("a" %||% "b", "a")
  expect_equal(42L %||% 0L, 42L)
})

test_that("%||% returns second value when first is NULL", {
  expect_equal(NULL %||% "b", "b")
  expect_equal(NULL %||% 99, 99)
})
