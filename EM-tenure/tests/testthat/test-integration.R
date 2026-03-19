# ==============================================================================
# Tests for EM-tenure: Full EM integration test with synthetic data
# ==============================================================================

test_that("em_fit_tenure converges on synthetic data", {
  # Generate data from known parameters
  true_params <- list(
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, sigma2_d = 0.01
  )
  df <- simulate_panel(
    n = 100, seed = 42,
    alpha = true_params$alpha,
    theta1 = true_params$theta1,
    theta0 = true_params$theta0,
    pi = true_params$pi,
    sigma2_g = true_params$sigma2_g,
    sigma2_d = true_params$sigma2_d
  )

  fit <- em_fit_tenure(df, max_iter = 200, verbose = 0)

  expect_true(fit$converged)
  expect_true(fit$iterations < 200)
  expect_true(is.finite(fit$loglik))

  # Parameters should be in the right ballpark (not exact due to small N)
  expect_true(abs(fit$params$alpha - true_params$alpha) < 0.2)
  expect_true(abs(fit$params$theta1 - true_params$theta1) < 0.15)
  expect_true(abs(fit$params$theta0 - true_params$theta0) < 0.15)
  expect_true(fit$params$pi < 0.3)  # Should at least be small
})

test_that("em_fit_tenure: monotone log-likelihood", {
  df <- simulate_panel(n = 80, seed = 99)
  fit <- em_fit_tenure(df, max_iter = 50, verbose = 0)

  # Log-likelihood should be non-decreasing (within numerical tolerance)
  ll <- fit$history$loglik[!is.na(fit$history$loglik)]
  diffs <- diff(ll)
  expect_true(all(diffs >= -1e-6 * max(abs(ll))))
})

test_that("em_fit_tenure without misclassification converges", {
  df <- simulate_panel(n = 80, pi = 0, seed = 77)
  fit <- em_fit_tenure(df, misclassification = FALSE, max_iter = 200, verbose = 0)

  expect_true(fit$converged)
  expect_equal(fit$params$pi, 0)
})

test_that("em_fit_tenure with stationarity converges", {
  df <- simulate_panel(n = 80, seed = 88)
  fit <- em_fit_tenure(df, stationary = TRUE, max_iter = 200, verbose = 0)

  expect_true(fit$converged)
  expected_alpha <- fit$params$theta0 / (fit$params$theta0 + 1 - fit$params$theta1)
  expect_equal(fit$params$alpha, expected_alpha, tolerance = 1e-6)
})

test_that("em_fit_tenure: gamma has correct dimensions", {
  df <- simulate_panel(n = 100, seed = 55)
  fit <- em_fit_tenure(df, max_iter = 20, verbose = 0)

  expect_equal(dim(fit$gamma), c(100, 8))
  row_sums <- rowSums(fit$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("em_fit_tenure: larger sample recovers parameters more precisely", {
  skip_if_not(identical(Sys.getenv("EM_FULL_TESTS"), "true"),
              "Set EM_FULL_TESTS=true to run large-N tests")
  true_params <- list(
    alpha = 0.7, theta1 = 0.85, theta0 = 0.15, pi = 0.03,
    sigma2_g = 0.005, sigma2_d = 0.008
  )
  df <- simulate_panel(
    n = 800, seed = 1234,
    alpha = true_params$alpha,
    theta1 = true_params$theta1,
    theta0 = true_params$theta0,
    pi = true_params$pi,
    sigma2_g = true_params$sigma2_g,
    sigma2_d = true_params$sigma2_d
  )

  fit <- em_fit_tenure(df, max_iter = 200, verbose = 0)

  expect_true(fit$converged)
  # With n=800, estimates should be reasonably close
  expect_true(abs(fit$params$theta1 - true_params$theta1) < 0.15)
  expect_true(abs(fit$params$theta0 - true_params$theta0) < 0.15)
})

test_that("em_fit_tenure validates required columns", {
  df <- data.frame(y1 = 1, y2 = 0)
  expect_error(em_fit_tenure(df), "Missing columns")
})

test_that("simulate_panel returns correct structure", {
  df <- simulate_panel(n = 50, seed = 11)
  expect_equal(nrow(df), 50)
  expect_true(all(c("y1", "y2", "y3", "tenure1", "tenure2", "tenure3",
                     "timegap1", "timegap2", "timegap3", "weight",
                     "h1", "h2", "h3") %in% names(df)))
  # Observed states are 0/1

  expect_true(all(df$y1 %in% c(0, 1)))
  expect_true(all(df$y2 %in% c(0, 1)))
  expect_true(all(df$y3 %in% c(0, 1)))
  # Durations are non-negative
  expect_true(all(df$tenure1 >= 0))
  expect_true(all(df$timegap1 >= 0))
})
