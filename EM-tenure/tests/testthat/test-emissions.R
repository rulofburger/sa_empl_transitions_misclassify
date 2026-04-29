# ==============================================================================
# Tests for EM-tenure: emissions.R
# ==============================================================================

test_that("log_misclass_prob returns correct length", {
  hmat <- latent_histories()
  lq <- log_misclass_prob(c(1, 1, 1), hmat, pi = 0.05)
  expect_length(lq, 8)
})

test_that("log_misclass_prob: no misclassification when pi=0", {
  hmat <- latent_histories()
  s <- c(1, 0, 1)
  lq <- log_misclass_prob(s, hmat, pi = 0)
  # Only history (1,0,1) should have log-prob = 0
  idx_match <- which(hmat[, 1] == 1 & hmat[, 2] == 0 & hmat[, 3] == 1)
  expect_equal(lq[idx_match], 0)
  # All others should be -Inf
  expect_true(all(lq[-idx_match] == -Inf))
})

test_that("log_misclass_prob: pi=0.5 gives equal weight to all histories", {
  hmat <- latent_histories()
  s <- c(1, 0, 1)
  lq <- log_misclass_prob(s, hmat, pi = 0.5)
  # All should be equal: 3 * log(0.5)
  expect_true(all(abs(lq - 3 * log(0.5)) < 1e-10))
})

test_that("log_emg returns finite values for positive inputs", {
  x_vals <- c(0.1, 0.5, 1, 2, 5)
  for (x in x_vals) {
    result <- log_emg(x, lambda = 0.5, sigma2 = 0.01)
    expect_true(is.finite(result))
  }
})

test_that("log_emg: EMG integrates to 1 (numerical check)", {
  lambda <- 0.3
  sigma2 <- 0.02
  # Integrate exp(log_emg) over a wide range
  f <- function(x) exp(log_emg(x, lambda, sigma2))
  integral <- integrate(f, lower = -10, upper = 100, subdivisions = 2000)
  expect_equal(integral$value, 1, tolerance = 0.05)
})

test_that("log_emission_increment_g returns correct Normal(0, 2*sigma2) density", {
  sigma2 <- 0.05
  delta <- 0.1
  # N(0, 2*sigma2) log-density
  expected <- dnorm(delta, mean = 0, sd = sqrt(2 * sigma2), log = TRUE)
  result <- log_emission_increment_g(delta, sigma2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("log_emission_increment_d matches g version with its own sigma", {
  sigma2 <- 0.03
  delta <- -0.05
  expected <- dnorm(delta, mean = 0, sd = sqrt(2 * sigma2), log = TRUE)
  result <- log_emission_increment_d(delta, sigma2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("log_emission_start_g returns N(.QUARTER_YEARS, sigma2_g) density", {
  sigma2 <- 0.04
  g_val <- 0.3
  expected <- dnorm(g_val, mean = .QUARTER_YEARS, sd = sqrt(sigma2), log = TRUE)
  result <- log_emission_start_g(g_val, sigma2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("log_emission_start_d returns N(.QUARTER_YEARS, sigma2_d) density", {
  sigma2 <- 0.02
  d_val <- 0.2
  expected <- dnorm(d_val, mean = .QUARTER_YEARS, sd = sqrt(sigma2), log = TRUE)
  result <- log_emission_start_d(d_val, sigma2)
  expect_equal(result, expected, tolerance = 1e-10)
})

# --- log_emg_grad_lambda tests -----------------------------------------------

test_that("log_emg_grad_lambda matches numerical finite difference", {
  lambda <- 0.4
  sigma2 <- 0.02
  x_vals <- c(0.1, 0.5, 1.0, 2.0)
  eps    <- 1e-6

  for (x in x_vals) {
    analytical <- log_emg_grad_lambda(x, lambda, sigma2)
    numerical  <- (log_emg(x, lambda + eps, sigma2) - log_emg(x, lambda - eps, sigma2)) / (2 * eps)
    expect_equal(analytical, numerical, tolerance = 1e-5,
                 label = paste0("grad at x=", x))
  }
})

test_that("log_emg_grad_lambda is vectorised and returns correct length", {
  x_vec  <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  result <- log_emg_grad_lambda(x_vec, lambda = 0.3, sigma2 = 0.01)
  expect_length(result, length(x_vec))
  expect_true(all(is.finite(result)))
})

test_that("log_emg_grad_lambda returns NA for x <= 0", {
  result_zero <- log_emg_grad_lambda(0, lambda = 0.3, sigma2 = 0.01)
  result_neg  <- log_emg_grad_lambda(-1, lambda = 0.3, sigma2 = 0.01)
  expect_true(is.na(result_zero))
  expect_true(is.na(result_neg))
})

test_that("log_emg_grad_lambda is stable for large lambda*sigma2 (large neg u)", {
  # When lambda*sigma2 >> x, u = (x - lambda*sigma2)/sigma is very negative
  # Phi(u) -> 0, and the log-space computation avoids underflow
  result <- log_emg_grad_lambda(x = 0.1, lambda = 10, sigma2 = 1)
  expect_true(is.finite(result))
})

test_that("log_emg_grad_lambda: different lambda/sigma2 combinations", {
  # Cross-check: two different parameterisations against numerical gradient
  cases <- list(
    list(lambda = 0.1, sigma2 = 0.001, x = 0.5),
    list(lambda = 2.0, sigma2 = 0.05,  x = 1.0),
    list(lambda = 0.5, sigma2 = 0.1,   x = 0.3)
  )
  eps <- 1e-7
  for (case in cases) {
    analytical <- log_emg_grad_lambda(case$x, case$lambda, case$sigma2)
    numerical  <- (log_emg(case$x, case$lambda + eps, case$sigma2) -
                   log_emg(case$x, case$lambda - eps, case$sigma2)) / (2 * eps)
    expect_equal(analytical, numerical, tolerance = 1e-4,
                 label = paste0("lambda=", case$lambda, " sigma2=", case$sigma2))
  }
})
