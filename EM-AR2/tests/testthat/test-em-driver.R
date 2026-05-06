# ==============================================================================
# EM-AR2 tests: em_driver.R
# ==============================================================================

# Helper: simulate data from AR(2) model with misclassification
.simulate_ar2 <- function(N, theta0, theta01, theta1, theta10, pi, seed = 1) {
  set.seed(seed)

  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat, theta0, theta01, theta1, theta10)

  # Sample latent histories
  h_idx <- sample.int(16, size = N, replace = TRUE, prob = prior)
  h     <- hmat[h_idx, ]

  # Apply misclassification
  y <- h
  for (t in 1:4) {
    flip <- runif(N) < pi
    y[, t] <- ifelse(flip, 1L - h[, t], h[, t])
  }

  data.frame(
    y1 = y[, 1], y2 = y[, 2], y3 = y[, 3], y4 = y[, 4],
    weight = rep(1, N)
  )
}

test_that("em_fit_ar2 errors on missing columns", {
  bad_df <- data.frame(y1 = 1, y2 = 1, y3 = 1)  # missing y4 and weight
  expect_error(em_fit_ar2(bad_df), regexp = "missing columns")
})

test_that("em_fit_ar2 errors on non-binary y values", {
  bad_df <- data.frame(
    y1 = c(0L, 1L, 2L), y2 = c(0L, 1L, 0L),
    y3 = c(1L, 0L, 1L), y4 = c(0L, 1L, 0L),
    weight = c(1, 1, 1)
  )
  expect_error(em_fit_ar2(bad_df), regexp = "must contain only 0/1 values")
})

test_that("em_fit_ar2 errors on NA weights", {
  bad_df <- data.frame(
    y1 = c(0L, 1L), y2 = c(1L, 0L), y3 = c(0L, 1L), y4 = c(1L, 0L),
    weight = c(1, NA_real_)
  )
  expect_error(em_fit_ar2(bad_df), regexp = "NA values in weights")
})

test_that("params_from_values errors when theta0 + theta01 >= 1", {
  expect_error(
    params_from_values(theta0 = 0.6, theta01 = 0.45, theta1 = 0.1, theta10 = 0.05),
    regexp = "theta0 \\+ theta01 must be < 1"
  )
})

test_that("params_from_values errors when theta1 + theta10 >= 1", {
  expect_error(
    params_from_values(theta0 = 0.1, theta01 = 0.05, theta1 = 0.6, theta10 = 0.45),
    regexp = "theta1 \\+ theta10 must be < 1"
  )
})

test_that("em_fit_ar2 returns correct output structure", {
  df <- .simulate_ar2(200, 0.10, 0.15, 0.08, 0.12, pi = 0.05)

  fit <- em_fit_ar2(df, max_iter = 20L, verbose = 0L)

  expect_named(fit, c("params", "loglik", "history", "converged", "iterations", "gamma"))
  expect_true(is.list(fit$params))
  expect_true(is.numeric(fit$loglik))
  expect_true(is.data.frame(fit$history))
  expect_equal(dim(fit$gamma), c(200L, 16L))
})

test_that("em_fit_ar2 log-likelihood is non-decreasing across iterations", {
  df  <- .simulate_ar2(300, 0.10, 0.15, 0.08, 0.12, pi = 0.05, seed = 7)
  fit <- em_fit_ar2(df, max_iter = 50L, verbose = 0L)

  ll_history <- fit$history$loglik
  diffs <- diff(ll_history)

  # Allow tiny numerical decrease (floating point)
  expect_true(all(diffs > -1e-4))
})

test_that("em_fit_ar2 no-ME variant converges", {
  df  <- .simulate_ar2(300, 0.10, 0.15, 0.08, 0.12, pi = 0, seed = 99)
  fit <- em_fit_ar2(df, estimate_pi = FALSE, fixed_pi = 0,
                    max_iter = 100L, verbose = 0L)

  expect_equal(fit$params$pi, 0)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_ar2 asymmetric variant runs and returns pi0/pi1", {
  df  <- .simulate_ar2(300, 0.10, 0.15, 0.08, 0.12, pi = 0.05, seed = 11)
  fit <- em_fit_ar2(df, asymmetric = TRUE, max_iter = 30L, verbose = 0L)

  expect_true("pi0" %in% names(fit$params))
  expect_true("pi1" %in% names(fit$params))
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_ar2 converges to near-true params on N=2000 data", {
  # Large sample to ensure EM converges close to truth
  true_theta0 <- 0.10; true_theta01 <- 0.15
  true_theta1 <- 0.08; true_theta10 <- 0.12
  pi <- 0.00  # No misclassification for simplicity

  df  <- .simulate_ar2(2000, true_theta0, true_theta01, true_theta1, true_theta10,
                        pi = pi, seed = 123)
  fit <- em_fit_ar2(df, estimate_pi = FALSE, fixed_pi = 0,
                    max_iter = 200L, tol = 1e-8, verbose = 0L)

  expect_equal(fit$params$theta0,  true_theta0,  tolerance = 0.02)
  expect_equal(fit$params$theta1,  true_theta1,  tolerance = 0.02)
  expect_equal(fit$params$theta01, true_theta01, tolerance = 0.04)
  expect_equal(fit$params$theta10, true_theta10, tolerance = 0.04)
})

test_that("init_params_ar2 returns all required parameter names including alpha", {
  df     <- .simulate_ar2(100, 0.1, 0.15, 0.08, 0.12, pi=0.05)
  params <- init_params_ar2(df)
  expected_names <- c("alpha", "theta0", "theta01", "theta1", "theta10", "pi")
  expect_true(all(expected_names %in% names(params)))
  expect_equal(sum(params$alpha), 1, tolerance = 1e-10)
})

test_that("init_params_ar2 returns bounded values", {
  df     <- .simulate_ar2(100, 0.1, 0.15, 0.08, 0.12, pi=0.05)
  params <- init_params_ar2(df)
  # alpha sums to 1 and all components are positive
  expect_equal(sum(params$alpha), 1, tolerance = 1e-10)
  expect_true(all(params$alpha > 0))
  # theta parameters and pi must be in (0, 1)
  scalar_params <- c(params$theta0, params$theta01, params$theta1, params$theta10, params$pi)
  expect_true(all(scalar_params > 0))
  expect_true(all(scalar_params < 1))
  # Joint constraints
  expect_true(params$theta0 + params$theta01 < 1)
  expect_true(params$theta1 + params$theta10 < 1)
})

test_that("params_from_values wraps scalar inputs correctly and includes alpha", {
  params <- params_from_values(0.1, 0.2, 0.08, 0.15, pi = 0.05)
  expect_equal(params$theta0,  0.1)
  expect_equal(params$theta01, 0.2)
  expect_equal(params$theta1,  0.08)
  expect_equal(params$theta10, 0.15)
  expect_equal(params$pi,      0.05)
  expect_true("alpha" %in% names(params))
  expect_equal(sum(params$alpha), 1, tolerance = 1e-10)
})
