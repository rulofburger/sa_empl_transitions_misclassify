# ==============================================================================
# EM-baseline: Tests for em_driver.R (init_params, compute_observed_loglik,
#              em_fit_baseline)
# Created: 2026-05-05
# ==============================================================================

# ---- shared synthetic data generator ----------------------------------------

.make_panel <- function(n = 200, theta0 = 0.10, theta1 = 0.90,
                        pi = 0.05, seed = 99) {
  set.seed(seed)
  alpha <- theta0 / (theta0 + 1 - theta1)

  h1 <- rbinom(n, 1, alpha)
  h2 <- ifelse(h1 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  h3 <- ifelse(h2 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))

  flip <- function(h) ifelse(rbinom(length(h), 1, pi) == 1, 1 - h, h)
  data.frame(y1 = flip(h1), y2 = flip(h2), y3 = flip(h3), weight = rep(1, n))
}

# ---- tests: init_params -----------------------------------------------------

test_that("init_params symmetric returns required fields", {
  p <- init_params("symmetric", stationary = TRUE)
  expect_true(all(c("alpha", "theta0", "theta1", "pi") %in% names(p)))
})

test_that("init_params asymmetric returns pi0 and pi1, not pi", {
  p <- init_params("asymmetric", stationary = TRUE)
  expect_true(all(c("alpha", "theta0", "theta1", "pi0", "pi1") %in% names(p)))
  expect_false("pi" %in% names(p))
})

test_that("init_params none has no pi fields", {
  p <- init_params("none", stationary = TRUE)
  expect_false("pi"  %in% names(p))
  expect_false("pi0" %in% names(p))
  expect_false("pi1" %in% names(p))
})

test_that("init_params stationary: alpha = theta0/(theta0+1-theta1)", {
  p <- init_params("symmetric", stationary = TRUE)
  expected <- p$theta0 / (p$theta0 + 1 - p$theta1)
  expect_equal(p$alpha, expected, tolerance = 1e-12)
})

test_that("init_params returns all finite values in (0,1)", {
  for (mt in c("symmetric", "asymmetric", "none")) {
    p <- init_params(mt, stationary = TRUE)
    probs <- unlist(p)
    expect_true(all(is.finite(probs)))
    expect_true(all(probs > 0 & probs < 1))
  }
})

# ---- tests: compute_observed_loglik -----------------------------------------

test_that("compute_observed_loglik returns finite negative scalar", {
  df     <- .make_panel(n = 50)
  params <- init_params("symmetric")
  ll     <- compute_observed_loglik(df, params, model_type = "symmetric")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
  expect_true(ll <= 0)
})

test_that("compute_observed_loglik matches e_step loglik", {
  df     <- .make_panel(n = 50)
  params <- init_params("symmetric")
  ll1    <- compute_observed_loglik(df, params, "symmetric")
  ll2    <- e_step(df, params, "symmetric")$loglik
  expect_equal(ll1, ll2, tolerance = 1e-12)
})

# ---- tests: em_fit_baseline — structure -------------------------------------

test_that("em_fit_baseline returns correct structure", {
  df     <- .make_panel(n = 100)
  result <- em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                            max_iter = 20L, verbose = 0L)
  expect_named(result,
    c("params", "loglik", "history", "converged", "iterations", "gamma",
      "model_type", "stationary"))
  expect_equal(dim(result$gamma), c(100L, 8L))
  expect_true(is.data.frame(result$history))
  expect_true(is.finite(result$loglik))
})

# ---- tests: em_fit_baseline — convergence and monotonicity ------------------

test_that("EM log-likelihood is monotonically non-decreasing (free-alpha)", {
  # The free-alpha model has no stationarity approximation, so EM is
  # exactly monotone. The stationarity model uses the short-panel
  # approximation (update theta from transition counts, then derive alpha),
  # which can introduce tiny LL decreases (see EM baseline.tex, Sec 2.3).
  df     <- .make_panel(n = 300, seed = 7)
  result <- em_fit_baseline(df, model_type = "symmetric", stationary = FALSE,
                            max_iter = 200L, tol = 1e-10, verbose = 0L)
  ll_seq <- result$history$loglik
  diffs  <- diff(ll_seq)
  expect_true(all(diffs >= -1e-8))
})

test_that("EM converges under stationarity approximation", {
  # The stationarity approximation (EM baseline.tex, Sec 2.3) trades exact
  # EM monotonicity for closed-form updates. Strict monotonicity is only
  # guaranteed for the free-alpha model. For the stationary model we verify
  # that convergence is achieved and the final loglik is finite.
  df     <- .make_panel(n = 300, seed = 7)
  result <- em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                            max_iter = 500L, tol = 1e-8, verbose = 0L)
  expect_true(result$converged)
  expect_true(is.finite(result$loglik))
  expect_true(result$loglik <= 0)
})

test_that("EM converges on moderately-sized synthetic data", {
  df     <- .make_panel(n = 500, seed = 42)
  result <- em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                            max_iter = 500L, tol = 1e-8, verbose = 0L)
  expect_true(result$converged)
})

# ---- tests: em_fit_baseline — parameter recovery ----------------------------

test_that("no-error EM recovers theta0 and theta1 on clean Markov data", {
  # With no misclassification, EM should recover the true params accurately
  n       <- 2000
  theta0  <- 0.12
  theta1  <- 0.88
  alpha   <- theta0 / (theta0 + 1 - theta1)
  set.seed(123)
  h1 <- rbinom(n, 1, alpha)
  h2 <- ifelse(h1 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  h3 <- ifelse(h2 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  df <- data.frame(y1 = h1, y2 = h2, y3 = h3, weight = rep(1, n))

  result <- em_fit_baseline(df, model_type = "none", stationary = TRUE,
                            max_iter = 500L, verbose = 0L)
  expect_equal(result$params$theta0, theta0, tolerance = 0.02)
  expect_equal(result$params$theta1, theta1, tolerance = 0.02)
})

# ---- tests: model variants --------------------------------------------------

test_that("asymmetric model converges and returns pi0 and pi1", {
  df     <- .make_panel(n = 300, seed = 77)
  result <- em_fit_baseline(df, model_type = "asymmetric", stationary = TRUE,
                            max_iter = 300L, verbose = 0L)
  expect_true("pi0" %in% names(result$params))
  expect_true("pi1" %in% names(result$params))
  expect_true(result$params$pi0 >= 0 & result$params$pi0 < 0.49)
  expect_true(result$params$pi1 >= 0 & result$params$pi1 < 0.49)
})

test_that("free-alpha model has loglik >= stationary model loglik", {
  # Free model is less restricted, so LL should be weakly higher
  df       <- .make_panel(n = 400, seed = 55)
  fit_stat <- em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                              max_iter = 300L, verbose = 0L)
  fit_free <- em_fit_baseline(df, model_type = "symmetric", stationary = FALSE,
                              max_iter = 300L, verbose = 0L)
  expect_true(fit_free$loglik >= fit_stat$loglik - 1e-4)
})

test_that("compute_observed_loglik matches em_fit_baseline final loglik", {
  df     <- .make_panel(n = 100, seed = 11)
  result <- em_fit_baseline(df, model_type = "symmetric", stationary = TRUE,
                            max_iter = 100L, verbose = 0L)
  ll_check <- compute_observed_loglik(df, result$params, model_type = "symmetric")
  expect_equal(result$loglik, ll_check, tolerance = 1e-8)
})

# ---- tests: error handling --------------------------------------------------

test_that("em_fit_baseline errors on missing columns", {
  bad_df <- data.frame(y1 = 1, y2 = 0, weight = 1)  # missing y3
  expect_error(em_fit_baseline(bad_df), regexp = "missing required columns")
})

test_that("em_fit_baseline errors on unknown model_type", {
  df <- .make_panel(n = 10)
  expect_error(em_fit_baseline(df, model_type = "fmm"), regexp = "model_type must be")
})

# ---- tests: em_fit_baseline — params0 validation ----------------------------

test_that("em_fit_baseline errors on params0 with NA value", {
  df      <- .make_panel(n = 30)
  params0 <- list(alpha = NA_real_, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  expect_error(em_fit_baseline(df, params0 = params0, model_type = "symmetric",
                               verbose = 0L),
               regexp = "finite numeric")
})

test_that("em_fit_baseline errors on params0 with Inf value", {
  df      <- .make_panel(n = 30)
  params0 <- list(alpha = Inf, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  expect_error(em_fit_baseline(df, params0 = params0, model_type = "symmetric",
                               verbose = 0L),
               regexp = "finite numeric")
})

test_that("em_fit_baseline errors on params0 with value > 1", {
  df      <- .make_panel(n = 30)
  params0 <- list(alpha = 1.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  expect_error(em_fit_baseline(df, params0 = params0, model_type = "symmetric",
                               verbose = 0L),
               regexp = "\\(0, 1\\)")
})

test_that("em_fit_baseline errors on params0 with value < 0", {
  df      <- .make_panel(n = 30)
  params0 <- list(alpha = -0.1, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  expect_error(em_fit_baseline(df, params0 = params0, model_type = "symmetric",
                               verbose = 0L),
               regexp = "\\(0, 1\\)")
})

test_that("em_fit_baseline errors on params0 missing required fields", {
  df      <- .make_panel(n = 30)
  params0 <- list(alpha = 0.5, theta1 = 0.9, pi = 0.05)  # missing theta0
  expect_error(em_fit_baseline(df, params0 = params0, model_type = "symmetric",
                               verbose = 0L),
               regexp = "missing fields")
})
