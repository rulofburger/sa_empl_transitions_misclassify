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

# ==============================================================================
# Timegap category constants (discrete/censored timegap model)
# ==============================================================================

test_that(".TIMEGAP_BOUNDS_YEARS has length 8 with correct endpoints", {
  expect_length(.TIMEGAP_BOUNDS_YEARS, 8)
  expect_equal(.TIMEGAP_BOUNDS_YEARS[1], 0)
  expect_equal(.TIMEGAP_BOUNDS_YEARS[8], Inf)
  # Strictly increasing (finite part)
  finite_bounds <- .TIMEGAP_BOUNDS_YEARS[is.finite(.TIMEGAP_BOUNDS_YEARS)]
  expect_true(all(diff(finite_bounds) > 0))
})

test_that(".N_TIMEGAP_CATS equals 7", {
  expect_equal(.N_TIMEGAP_CATS, 7L)
})

test_that(".TIMEGAP_MIDPOINTS_MONTHS has length 7 and is in months", {
  expect_length(.TIMEGAP_MIDPOINTS_MONTHS, 7)
  # Each midpoint should lie within its category interval (converted to years)
  for (k in seq_len(.N_TIMEGAP_CATS)) {
    m_yr <- .TIMEGAP_MIDPOINTS_MONTHS[k] / 12
    iv <- .timegap_interval(k)
    expect_true(m_yr >= iv[1] && (is.infinite(iv[2]) || m_yr < iv[2]),
                label = paste0("midpoint of cat ", k, " in interval"))
  }
})

test_that(".DURATION_FLOOR is positive and very small", {
  expect_gt(.DURATION_FLOOR, 0)
  expect_lt(.DURATION_FLOOR, 0.1)  # less than 1.2 months
})

test_that(".timegap_interval returns correct bounds for all categories", {
  expected <- list(
    c(0,    0.25),
    c(0.25, 0.50),
    c(0.50, 0.75),
    c(0.75, 1.00),
    c(1.00, 3.00),
    c(3.00, 5.00),
    c(5.00, Inf)
  )
  for (k in seq_len(.N_TIMEGAP_CATS)) {
    expect_equal(.timegap_interval(k), expected[[k]],
                 label = paste0(".timegap_interval(", k, ")"))
  }
})

test_that(".timegap_interval errors on out-of-range k", {
  expect_error(.timegap_interval(0))
  expect_error(.timegap_interval(8))
})

test_that(".continuous_to_cat maps durations to correct categories", {
  # Representative values within each category
  expect_equal(.continuous_to_cat(0.0),   1L)   # left boundary of cat 1
  expect_equal(.continuous_to_cat(0.1),   1L)
  expect_equal(.continuous_to_cat(0.25),  2L)   # left boundary of cat 2
  expect_equal(.continuous_to_cat(0.3),   2L)
  expect_equal(.continuous_to_cat(0.5),   3L)   # left boundary of cat 3
  expect_equal(.continuous_to_cat(0.75),  4L)   # left boundary of cat 4
  expect_equal(.continuous_to_cat(1.0),   5L)   # left boundary of cat 5
  expect_equal(.continuous_to_cat(2.5),   5L)
  expect_equal(.continuous_to_cat(3.0),   6L)   # left boundary of cat 6
  expect_equal(.continuous_to_cat(5.0),   7L)   # left boundary of cat 7
  expect_equal(.continuous_to_cat(100.0), 7L)   # large value → cat 7
})

test_that(".continuous_to_cat returns NA for NA input", {
  expect_true(is.na(.continuous_to_cat(NA_real_)))
})

test_that(".continuous_to_cat is vectorised", {
  d <- c(0.1, 0.3, 0.6, 0.9, 1.5, 4.0, 6.0)
  result <- .continuous_to_cat(d)
  expect_equal(result, c(1L, 2L, 3L, 4L, 5L, 6L, 7L))
})
