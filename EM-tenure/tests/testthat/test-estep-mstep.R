# ==============================================================================
# Tests for EM-tenure: E-step and M-step
# ==============================================================================

# --- E-step tests ---

test_that("e_step returns correct structure", {
  df <- simulate_panel(n = 50, seed = 42)
  params <- init_params(df)
  params$lambda_g <- ctmc_lambda_from_theta(params$theta1)
  params$lambda_d <- ctmc_lambda_from_theta(params$theta0)

  result <- e_step(df, params)

  expect_true(is.list(result))
  expect_true("gamma" %in% names(result))
  expect_true("loglik" %in% names(result))
  expect_true("suff" %in% names(result))

  # gamma should be N x 8

  expect_equal(dim(result$gamma), c(50, 8))

  # responsibilities sum to 1 per row
  row_sums <- rowSums(result$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # loglik should be finite
  expect_true(is.finite(result$loglik))
})

test_that("e_step: responsibilities sum correctly", {
  df <- simulate_panel(n = 100, seed = 123)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  result <- e_step(df, params)

  # All gamma values should be in [0, 1]
  expect_true(all(result$gamma >= -1e-15))
  expect_true(all(result$gamma <= 1 + 1e-15))
})

test_that("e_step: sufficient stats are non-negative", {
  df <- simulate_panel(n = 100, seed = 456)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  result <- e_step(df, params)

  expect_true(result$suff$C1 >= 0)
  expect_true(result$suff$C0 >= 0)
  expect_true(result$suff$D1 >= 0)
  expect_true(result$suff$D0 >= 0)
  expect_true(result$suff$Ng >= 0)
  expect_true(result$suff$Nd >= 0)
})

# --- M-step tests ---

test_that("m_step returns all required parameters", {
  df <- simulate_panel(n = 100, seed = 789)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params)
  new_params <- m_step(estep_out$suff, total_weight = sum(df$weight))

  expect_true(all(c("alpha", "theta0", "theta1", "pi",
                     "sigma2_g", "sigma2_d", "lambda_g", "lambda_d") %in%
                   names(new_params)))

  # All in valid ranges
  expect_true(new_params$alpha > 0 && new_params$alpha < 1)
  expect_true(new_params$theta1 > 0 && new_params$theta1 < 1)
  expect_true(new_params$theta0 > 0 && new_params$theta0 < 1)
  expect_true(new_params$pi >= 0 && new_params$pi < 0.5)
  expect_true(new_params$sigma2_g > 0)
  expect_true(new_params$sigma2_d > 0)
  expect_true(new_params$lambda_g > 0)
  expect_true(new_params$lambda_d > 0)
})

test_that("m_step without misclassification sets pi=0", {
  df <- simulate_panel(n = 100, seed = 101)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params)
  new_params <- m_step(estep_out$suff, sum(df$weight), misclassification = FALSE)
  expect_equal(new_params$pi, 0)
})

test_that("m_step with stationarity enforces alpha = theta0/(theta0+1-theta1)", {
  df <- simulate_panel(n = 200, seed = 202)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params)
  new_params <- m_step(estep_out$suff, sum(df$weight), stationary = TRUE)
  expected_alpha <- new_params$theta0 / (new_params$theta0 + 1 - new_params$theta1)
  expect_equal(new_params$alpha, expected_alpha, tolerance = 1e-10)
})

test_that("m_step lambda linked to theta", {
  df <- simulate_panel(n = 100, seed = 303)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params)
  new_params <- m_step(estep_out$suff, sum(df$weight))
  expect_equal(new_params$lambda_g, ctmc_lambda_from_theta(new_params$theta1), tolerance = 1e-12)
  expect_equal(new_params$lambda_d, ctmc_lambda_from_theta(new_params$theta0), tolerance = 1e-12)
})
