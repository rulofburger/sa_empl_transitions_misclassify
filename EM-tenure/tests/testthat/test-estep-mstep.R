# ==============================================================================
# Tests for EM-tenure: E-step and M-step
# ==============================================================================

# --- E-step tests ---

test_that("e_step returns correct structure", {
  df <- simulate_panel(n = 50, seed = 42, discrete_timegap = FALSE)
  params <- init_params(df, discrete_timegap = FALSE)
  params$lambda_g <- ctmc_lambda_from_theta(params$theta1)
  params$lambda_d <- ctmc_lambda_from_theta(params$theta0)

  result <- e_step(df, params, discrete_timegap = FALSE)

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
  df <- simulate_panel(n = 100, seed = 123, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  result <- e_step(df, params, discrete_timegap = FALSE)

  # All gamma values should be in [0, 1]
  expect_true(all(result$gamma >= -1e-15))
  expect_true(all(result$gamma <= 1 + 1e-15))
})

test_that("e_step: sufficient stats are non-negative", {
  df <- simulate_panel(n = 100, seed = 456, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  result <- e_step(df, params, discrete_timegap = FALSE)

  expect_true(result$suff$C1 >= 0)
  expect_true(result$suff$C0 >= 0)
  expect_true(result$suff$D1 >= 0)
  expect_true(result$suff$D0 >= 0)
  expect_true(result$suff$Ng >= 0)
  expect_true(result$suff$Nd >= 0)
})

# --- M-step tests ---

test_that("m_step returns all required parameters", {
  df <- simulate_panel(n = 100, seed = 789, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params, discrete_timegap = FALSE)
  new_params <- m_step(estep_out$suff, total_weight = sum(df$weight),
                       discrete_timegap = FALSE)

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
  df <- simulate_panel(n = 100, seed = 101, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params, discrete_timegap = FALSE)
  new_params <- m_step(estep_out$suff, sum(df$weight),
                       misclassification = FALSE, discrete_timegap = FALSE)
  expect_equal(new_params$pi, 0)
})

test_that("m_step with stationarity enforces alpha = theta0/(theta0+1-theta1)", {
  df <- simulate_panel(n = 200, seed = 202, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params, discrete_timegap = FALSE)
  new_params <- m_step(estep_out$suff, sum(df$weight),
                       stationary = TRUE, discrete_timegap = FALSE)
  expected_alpha <- new_params$theta0 / (new_params$theta0 + 1 - new_params$theta1)
  expect_equal(new_params$alpha, expected_alpha, tolerance = 1e-10)
})

test_that("m_step lambda linked to theta", {
  df <- simulate_panel(n = 100, seed = 303, discrete_timegap = FALSE)
  params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01,
    lambda_g = ctmc_lambda_from_theta(0.9),
    lambda_d = ctmc_lambda_from_theta(0.1)
  )
  estep_out <- e_step(df, params, discrete_timegap = FALSE)
  new_params <- m_step(estep_out$suff, sum(df$weight), discrete_timegap = FALSE)
  expect_equal(new_params$lambda_g, ctmc_lambda_from_theta(new_params$theta1), tolerance = 1e-12)
  expect_equal(new_params$lambda_d, ctmc_lambda_from_theta(new_params$theta0), tolerance = 1e-12)
})

# ==============================================================================
# e_step with discrete_timegap = TRUE
# ==============================================================================

test_that("e_step discrete: returns correct structure", {
  df <- simulate_panel(n = 80, seed = 501, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  out <- e_step(df, params, discrete_timegap = TRUE)

  expect_true(is.list(out))
  expect_true(all(c("gamma", "loglik", "suff") %in% names(out)))
  # gamma: n x 8

  expect_equal(nrow(out$gamma), nrow(df))
  expect_equal(ncol(out$gamma), 8)
  # responsibilities sum to 1
  row_sums <- rowSums(out$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))
  expect_true(is.finite(out$loglik))
})

test_that("e_step discrete: suff stats contain discrete category fields", {
  df <- simulate_panel(n = 80, seed = 502, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  out <- e_step(df, params, discrete_timegap = TRUE)

  suff <- out$suff
  # Discrete mode should have categorical marginals and transitions
  expect_true("cat_d_marginal_c" %in% names(suff))
  expect_true("cat_d_marginal_w" %in% names(suff))
  # Should NOT have continuous sufficient stats Sd, Nd
  expect_null(suff$Sd)
  expect_null(suff$Nd)
})

test_that("e_step discrete: perfect-classification data gives pi near 0", {
  df <- simulate_panel(n = 200, pi = 0, seed = 503, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  params$pi <- 0.01  # small starting pi
  out <- e_step(df, params, discrete_timegap = TRUE)

  # The diagonal histories (no misclassification) should dominate
  # gamma columns 1 (000→000 identity) and 8 (111→111 identity) should carry most weight
  diag_weight <- sum(out$gamma[, 1]) + sum(out$gamma[, 8])
  expect_true(diag_weight / nrow(df) > 0.8)
})

test_that("e_step discrete: rejects data with NA tenure", {
  df <- simulate_panel(n = 50, seed = 504, discrete_timegap = TRUE)
  df$tenure1[1] <- NA
  params <- init_params(simulate_panel(n = 50, seed = 505, discrete_timegap = TRUE),
                        discrete_timegap = TRUE)
  expect_error(e_step(df, params, discrete_timegap = TRUE), "NA")
})

# ==============================================================================
# m_step with discrete_timegap = TRUE
# ==============================================================================

test_that("m_step discrete: returns all required parameters (no sigma2_d)", {
  df <- simulate_panel(n = 100, seed = 601, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  estep_out <- e_step(df, params, discrete_timegap = TRUE)
  new_params <- m_step(estep_out$suff, sum(df$weight), discrete_timegap = TRUE)

  # Must have: alpha, theta1, theta0, pi, sigma2_g, lambda_g
  expect_true(all(c("alpha", "theta1", "theta0", "pi", "sigma2_g", "lambda_g") %in%
                  names(new_params)))
  # Must NOT have sigma2_d (only for continuous mode)
  expect_null(new_params$sigma2_d)
  # Parameters in valid ranges
  expect_true(new_params$alpha > 0 && new_params$alpha < 1)
  expect_true(new_params$theta1 > 0 && new_params$theta1 < 1)
  expect_true(new_params$theta0 > 0 && new_params$theta0 < 1)
  expect_true(new_params$pi >= 0 && new_params$pi < 0.5)
  expect_true(new_params$sigma2_g > 0)
})

test_that("m_step discrete: without misclassification sets pi=0", {
  df <- simulate_panel(n = 80, seed = 602, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  estep_out <- e_step(df, params, discrete_timegap = TRUE)
  new_params <- m_step(estep_out$suff, sum(df$weight),
                       misclassification = FALSE, discrete_timegap = TRUE)
  expect_equal(new_params$pi, 0)
})

test_that("m_step discrete: stationarity enforces alpha constraint", {
  df <- simulate_panel(n = 100, seed = 603, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  estep_out <- e_step(df, params, discrete_timegap = TRUE)
  new_params <- m_step(estep_out$suff, sum(df$weight),
                       stationary = TRUE, discrete_timegap = TRUE)
  expected_alpha <- new_params$theta0 / (new_params$theta0 + 1 - new_params$theta1)
  expect_equal(new_params$alpha, expected_alpha, tolerance = 1e-10)
})

test_that("m_step discrete: lambda linked to theta via FOC", {
  df <- simulate_panel(n = 100, seed = 604, discrete_timegap = TRUE)
  params <- init_params(df, discrete_timegap = TRUE)
  estep_out <- e_step(df, params, discrete_timegap = TRUE)
  new_params <- m_step(estep_out$suff, sum(df$weight), discrete_timegap = TRUE)

  # lambda_g should be consistent with theta1
  expect_equal(new_params$lambda_g, ctmc_lambda_from_theta(new_params$theta1),
               tolerance = 1e-12)
  # theta0 is solved via Brent; lambda_d should match
  expect_equal(new_params$lambda_d, ctmc_lambda_from_theta(new_params$theta0),
               tolerance = 1e-12)
})
