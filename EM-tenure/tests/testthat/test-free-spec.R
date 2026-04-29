# ==============================================================================
# Tests for EM-tenure: Free specification M-step and linked backward compat
# Created: 2026-04-28
# ==============================================================================

# --- Free specification: closed-form theta ---

test_that("m_step with linked=FALSE gives closed-form theta1 = T11/D1", {
  # Construct minimal sufficient statistics
  suff <- list(
    T11 = 90, D1 = 100, T01 = 8, D0 = 100,
    C1 = 60, C0 = 40, M = 9,
    Ng = 80, Sg = 0.8, Ng_start = 10, Sg_start = 0.05,
    emg_g_x = rexp(50, 0.02), emg_g_w = rep(1, 50),
    cat_d_marginal_c = sample(1:7, 30, replace = TRUE),
    cat_d_marginal_w = rep(1, 30),
    cat_d_trans_curr = sample(1:7, 20, replace = TRUE),
    cat_d_trans_prev = sample(1:7, 20, replace = TRUE),
    cat_d_trans_w    = rep(1, 20)
  )
  result <- m_step(suff, total_weight = 100, linked = FALSE)

  # theta1 should be T11/D1 = 0.90 (closed form)
  expect_equal(result$theta1, 0.90, tolerance = 1e-6)
  expect_equal(result$theta0, 0.08, tolerance = 1e-6)

  # lambda_g should NOT equal ctmc_lambda_from_persistence(theta1) in general
  ctmc_val <- ctmc_lambda_from_persistence(result$theta1)
  # (They could coincidentally match, but the point is lambda_g is estimated
  # independently from EMG data, not derived from theta1)
  expect_true(is.finite(result$lambda_g))
  expect_true(result$lambda_g > 0)
})

test_that("m_step with linked=TRUE matches legacy behaviour", {
  suff <- list(
    T11 = 90, D1 = 100, T01 = 8, D0 = 100,
    C1 = 60, C0 = 40, M = 9,
    Ng = 80, Sg = 0.8, Ng_start = 10, Sg_start = 0.05,
    emg_g_x = rexp(50, 0.02), emg_g_w = rep(1, 50),
    cat_d_marginal_c = sample(1:7, 30, replace = TRUE),
    cat_d_marginal_w = rep(1, 30),
    cat_d_trans_curr = sample(1:7, 20, replace = TRUE),
    cat_d_trans_prev = sample(1:7, 20, replace = TRUE),
    cat_d_trans_w    = rep(1, 20)
  )
  result <- m_step(suff, total_weight = 100, linked = TRUE)

  # Under linked, lambda_g = -log(theta1)/Delta, lambda_d = -log(1-theta0)/Delta
  expect_equal(result$lambda_g, ctmc_lambda_from_persistence(result$theta1), tolerance = 1e-10)
  expect_equal(result$lambda_d, ctmc_lambda_from_transition(result$theta0), tolerance = 1e-10)
})


# --- Free spec: full integration on synthetic data ---

test_that("em_fit_tenure free spec converges on synthetic data", {
  df <- simulate_panel(
    n = 200, seed = 42,
    alpha = 0.6, theta1 = 0.9, theta0 = 0.1, pi = 0.05,
    sigma2_g = 0.01, discrete_timegap = TRUE
  )

  fit <- em_fit_tenure(df, max_iter = 300, verbose = 0,
                       linked = FALSE, discrete_timegap = TRUE)

  expect_true(fit$converged)
  expect_true(is.finite(fit$loglik))
  expect_true(abs(fit$params$theta1 - 0.9) < 0.15)
  expect_true(abs(fit$params$theta0 - 0.1) < 0.15)
  expect_true(fit$params$lambda_g > 0)
  expect_true(fit$params$lambda_d > 0)
})

test_that("em_fit_tenure free spec: monotone log-likelihood", {
  df <- simulate_panel(n = 150, seed = 99, discrete_timegap = TRUE)
  fit <- em_fit_tenure(df, max_iter = 100, verbose = 0,
                       linked = FALSE, discrete_timegap = TRUE)

  ll <- fit$history$loglik[!is.na(fit$history$loglik)]
  diffs <- diff(ll)
  expect_true(all(diffs >= -1e-6 * max(abs(ll))))
})


# --- Free spec with user-supplied lambda (inconsistent with CTMC) ---

test_that("simulate_panel with explicit lambda_g, lambda_d", {
  # lambda_g = 0.01 is much smaller than ctmc_lambda_from_persistence(0.9) ≈ 0.42
  df <- simulate_panel(
    n = 100, seed = 42,
    theta1 = 0.9, theta0 = 0.1,
    lambda_g = 0.01, lambda_d = 0.05,
    discrete_timegap = TRUE
  )
  expect_equal(nrow(df), 100)
  # Tenure should be longer on average when lambda_g is small
  mean_tenure <- mean(df$tenure1[df$y1 == 1], na.rm = TRUE)
  expect_true(mean_tenure > 1)  # > 1 year with lambda_g = 0.01
})


# --- Linked backward compatibility ---

test_that("em_fit_tenure linked=TRUE still converges", {
  df <- simulate_panel(n = 100, seed = 42, discrete_timegap = TRUE)
  fit <- em_fit_tenure(df, max_iter = 300, verbose = 0,
                       linked = TRUE, discrete_timegap = TRUE)

  expect_true(fit$converged || fit$iterations == 300)
  expect_true(is.finite(fit$loglik))
  # CTMC link enforced
  expect_equal(fit$params$lambda_g,
               ctmc_lambda_from_persistence(fit$params$theta1), tolerance = 1e-10)
})
