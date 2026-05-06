# ==============================================================================
# EM-AR2 tests: utils.R
# ==============================================================================

test_that(".bound01 clamps values to (eps, 1-eps)", {
  expect_equal(.bound01(0), 1e-12)
  expect_equal(.bound01(1), 1 - 1e-12)
  expect_equal(.bound01(0.5), 0.5)
  expect_equal(.bound01(-1), 1e-12)
  expect_equal(.bound01(2), 1 - 1e-12)
})

test_that(".bound01 works with custom eps", {
  expect_equal(.bound01(0, eps = 0.01), 0.01)
  expect_equal(.bound01(1, eps = 0.01), 0.99)
  expect_equal(.bound01(0.5, eps = 0.01), 0.5)
})

test_that(".bound01 is vectorised", {
  x <- c(-1, 0, 0.3, 0.7, 1, 2)
  result <- .bound01(x)
  expect_length(result, 6)
  expect_true(all(result > 0))
  expect_true(all(result < 1))
})

test_that(".logsumexp is accurate for typical values", {
  x <- c(1, 2, 3)
  expected <- log(exp(1) + exp(2) + exp(3))
  expect_equal(.logsumexp(x), expected, tolerance = 1e-10)
})

test_that(".logsumexp handles large values without overflow", {
  x <- c(800, 801, 802)
  result <- .logsumexp(x)
  expect_true(is.finite(result))
  expect_true(result > 802)
})

test_that(".logsumexp handles a single value", {
  expect_equal(.logsumexp(3.14), 3.14)
})

test_that(".logsumexp returns -Inf for -Inf input", {
  # All -Inf: sum of exp(-Inf) = 0 -> log(0) = -Inf
  # After max subtraction: max = -Inf is not finite -> returns -Inf directly
  expect_equal(.logsumexp(c(-Inf, -Inf)), -Inf)
})

test_that(".row_logsumexp matches row-wise .logsumexp", {
  set.seed(42)
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  row_result <- .row_logsumexp(mat)
  expected <- apply(mat, 1, .logsumexp)
  expect_equal(row_result, expected, tolerance = 1e-10)
})

test_that("%||% returns lhs when non-NULL", {
  expect_equal(5L %||% 10L, 5L)
  expect_equal("hello" %||% "world", "hello")
})

test_that("%||% returns rhs when lhs is NULL", {
  expect_equal(NULL %||% 10L, 10L)
  x <- NULL
  expect_equal(x %||% "default", "default")
})
