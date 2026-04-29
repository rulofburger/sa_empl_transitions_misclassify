# ==============================================================================
# Tests for EM-tenure: discrete interval-censored timegap emissions
# ==============================================================================
# Tests the five new emission functions added for Issue 3 (discrete timegap):
#   log_emission_interval_d()
#   log_emission_transition_d()
#   log_emission_start_d_cat()
#   interval_grad_lambda_d()
#   transition_grad_lambda_d()
# ==============================================================================

# ==============================================================================
# log_emission_interval_d
# ==============================================================================

test_that("log_emission_interval_d: probs sum to 1 over all categories", {
  for (lambda_d in c(0.1, 0.5, 1.0, 2.0, 5.0)) {
    log_probs <- log_emission_interval_d(1:.N_TIMEGAP_CATS, lambda_d)
    total <- sum(exp(log_probs))
    expect_equal(total, 1.0, tolerance = 1e-8,
                 label = paste0("lambda_d=", lambda_d))
  }
})

test_that("log_emission_interval_d: returns -Inf for invalid categories", {
  lambda_d <- 0.5
  expect_equal(log_emission_interval_d(0L,  lambda_d), -Inf)
  expect_equal(log_emission_interval_d(8L,  lambda_d), -Inf)
  expect_equal(log_emission_interval_d(-1L, lambda_d), -Inf)
  expect_true(is.infinite(log_emission_interval_d(NA_integer_, lambda_d)))
})

test_that("log_emission_interval_d: matches numerical integration of dexp", {
  lambda_d <- 0.8
  for (k in 1:.N_TIMEGAP_CATS) {
    iv <- .timegap_interval(k)
    expected <- integrate(function(x) dexp(x, rate = lambda_d),
                          lower = iv[1],
                          upper = if (is.infinite(iv[2])) 1e4 else iv[2])$value
    # For k=7, upper=Inf is approximated by 1e4; exact = exp(-lambda*5)
    if (k == .N_TIMEGAP_CATS) expected <- exp(-lambda_d * iv[1])
    got <- exp(log_emission_interval_d(k, lambda_d))
    expect_equal(got, expected, tolerance = 1e-6,
                 label = paste0("cat=", k, " lambda_d=", lambda_d))
  }
})

test_that("log_emission_interval_d: vectorised over cat", {
  lambda_d <- 0.5
  cats <- c(1L, 3L, 5L, 7L)
  result <- log_emission_interval_d(cats, lambda_d)
  expect_length(result, 4)
  expect_true(all(result <= 0))   # log probabilities <= 0
  expect_true(all(is.finite(result)))
})

test_that("log_emission_interval_d: larger lambda concentrates mass in lower cats", {
  # At high lambda, Exp(lambda) concentrates near 0 -> most mass in cat 1
  p_large_lambda <- exp(log_emission_interval_d(1L, lambda_d = 10))
  p_small_lambda <- exp(log_emission_interval_d(1L, lambda_d = 0.1))
  expect_gt(p_large_lambda, p_small_lambda)
})

# ==============================================================================
# log_emission_transition_d
# ==============================================================================

test_that("log_emission_transition_d: row sums to 1 for each source category", {
  lambda_d <- 0.5
  for (j in 1:.N_TIMEGAP_CATS) {
    log_probs <- log_emission_transition_d(1:.N_TIMEGAP_CATS, rep(j, .N_TIMEGAP_CATS), lambda_d)
    # Only finite entries; should sum to 1
    finite_probs <- exp(log_probs[is.finite(log_probs)])
    expect_equal(sum(finite_probs), 1.0, tolerance = 1e-6,
                 label = paste0("source cat=", j, " lambda_d=", lambda_d))
  }
})

test_that("log_emission_transition_d: categories 1-4 are deterministic", {
  # Categories 1-4 have narrow 3-month intervals; adding 0.25 years shifts
  # the entire interval exactly into the next category.
  lambda_d <- 0.5
  deterministic_pairs <- list(
    c(2L, 1L),  # cat 1 -> cat 2
    c(3L, 2L),  # cat 2 -> cat 3
    c(4L, 3L),  # cat 3 -> cat 4
    c(5L, 4L)   # cat 4 -> cat 5
  )
  for (pair in deterministic_pairs) {
    k <- pair[1]; j <- pair[2]
    lp <- log_emission_transition_d(k, j, lambda_d)
    expect_equal(exp(lp), 1.0, tolerance = 1e-8,
                 label = paste0("cat_prev=", j, " -> cat_curr=", k))
  }
})

test_that("log_emission_transition_d: unreachable transitions return -Inf", {
  lambda_d <- 0.5
  # Cat 1 -> Cat 1 is unreachable: [0, 0.25) + 0.25 = [0.25, 0.5) = cat 2
  expect_equal(log_emission_transition_d(1L, 1L, lambda_d), -Inf)
  # Cat 1 -> Cat 3 is unreachable
  expect_equal(log_emission_transition_d(3L, 1L, lambda_d), -Inf)
  # Cat 7 -> Cat 1 is unreachable
  expect_equal(log_emission_transition_d(1L, 7L, lambda_d), -Inf)
})

test_that("log_emission_transition_d: cat 7 is absorbing (stays in 7)", {
  lambda_d <- 0.5
  lp <- log_emission_transition_d(7L, 7L, lambda_d)
  expect_equal(exp(lp), 1.0, tolerance = 1e-8)
})

test_that("log_emission_transition_d: vectorised over paired cats", {
  lambda_d <- 0.5
  cat_prev <- c(1L, 2L, 3L, 4L, 7L)
  cat_curr <- c(2L, 3L, 4L, 5L, 7L)
  result <- log_emission_transition_d(cat_curr, cat_prev, lambda_d)
  expect_length(result, 5)
  # All deterministic transitions should give prob = 1
  expect_true(all(abs(exp(result) - 1.0) < 1e-6))
})

test_that("log_emission_transition_d: invalid categories return -Inf", {
  lambda_d <- 0.5
  expect_equal(log_emission_transition_d(0L, 1L, lambda_d), -Inf)
  expect_equal(log_emission_transition_d(2L, 8L, lambda_d), -Inf)
  expect_true(is.infinite(log_emission_transition_d(NA_integer_, 1L, lambda_d)))
})

# ==============================================================================
# log_emission_start_d_cat
# ==============================================================================

test_that("log_emission_start_d_cat: returns 0 for cat=1", {
  expect_equal(log_emission_start_d_cat(1L), 0)
})

test_that("log_emission_start_d_cat: returns -Inf for cat > 1", {
  for (k in 2:.N_TIMEGAP_CATS) {
    expect_equal(log_emission_start_d_cat(k), -Inf,
                 label = paste0("cat=", k))
  }
})

test_that("log_emission_start_d_cat: returns -Inf for NA", {
  expect_equal(log_emission_start_d_cat(NA_integer_), -Inf)
})

test_that("log_emission_start_d_cat: vectorised", {
  result <- log_emission_start_d_cat(c(1L, 2L, 3L, 1L))
  expect_equal(result, c(0, -Inf, -Inf, 0))
})

# ==============================================================================
# interval_grad_lambda_d
# ==============================================================================

test_that("interval_grad_lambda_d: matches finite difference for all categories", {
  lambda_d <- 0.6
  eps <- 1e-7
  for (k in 1:.N_TIMEGAP_CATS) {
    analytical <- interval_grad_lambda_d(k, lambda_d)
    numerical  <- (log_emission_interval_d(k, lambda_d + eps) -
                   log_emission_interval_d(k, lambda_d - eps)) / (2 * eps)
    expect_equal(analytical, numerical, tolerance = 1e-4,
                 label = paste0("cat=", k))
  }
})

test_that("interval_grad_lambda_d: returns NA for invalid categories", {
  expect_true(is.na(interval_grad_lambda_d(0L, 0.5)))
  expect_true(is.na(interval_grad_lambda_d(8L, 0.5)))
  expect_true(is.na(interval_grad_lambda_d(NA_integer_, 0.5)))
})

test_that("interval_grad_lambda_d: vectorised", {
  result <- interval_grad_lambda_d(1:.N_TIMEGAP_CATS, lambda_d = 0.4)
  expect_length(result, .N_TIMEGAP_CATS)
  expect_true(all(!is.na(result)))
  # Gradient should be negative (increasing lambda decreases mean duration,
  # shifting mass toward lower categories) — especially for cat 7
  # (rightmost). Gradient = -a_k for cat 7 = -5.0.
  expect_equal(result[.N_TIMEGAP_CATS], -5.0, tolerance = 1e-10)
})

test_that("interval_grad_lambda_d: stable at large lambda", {
  result <- interval_grad_lambda_d(3L, lambda_d = 10)
  expect_true(is.finite(result))
})

# ==============================================================================
# transition_grad_lambda_d
# ==============================================================================

test_that("transition_grad_lambda_d: matches finite difference for valid pairs", {
  lambda_d <- 0.5
  eps <- 1e-7
  # Use deterministic pairs (categories 1-4 -> one step forward)
  pairs <- list(c(2L, 1L), c(3L, 2L), c(5L, 4L), c(7L, 7L))
  for (pair in pairs) {
    k <- pair[1]; j <- pair[2]
    analytical <- transition_grad_lambda_d(k, j, lambda_d)
    numerical  <- (log_emission_transition_d(k, j, lambda_d + eps) -
                   log_emission_transition_d(k, j, lambda_d - eps)) / (2 * eps)
    # Both may be NA for degenerate cases; skip those
    if (!is.na(analytical) && is.finite(numerical)) {
      expect_equal(analytical, numerical, tolerance = 1e-4,
                   label = paste0("cat_prev=", j, " cat_curr=", k))
    }
  }
})

test_that("transition_grad_lambda_d: returns NA for unreachable transitions", {
  # cat 1 -> cat 1 is unreachable (log_emission_transition_d = -Inf)
  result <- transition_grad_lambda_d(1L, 1L, lambda_d = 0.5)
  expect_true(is.na(result))
})

test_that("transition_grad_lambda_d: vectorised", {
  cat_prev <- c(1L, 2L, 3L, 4L)
  cat_curr <- c(2L, 3L, 4L, 5L)
  result <- transition_grad_lambda_d(cat_curr, cat_prev, lambda_d = 0.5)
  expect_length(result, 4)
})
