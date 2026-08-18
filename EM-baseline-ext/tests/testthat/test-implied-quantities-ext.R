# ==============================================================================
# EM-baseline-ext: Tests for implied_quantities_ext.R
# Created: 2026-05-06
# ==============================================================================

# ---------------------------------------------------------------------------
# implied_covariates()
# ---------------------------------------------------------------------------

test_that("implied_covariates returns correct names", {
  X    <- cbind(1, c(-1, 0, 1))
  p    <- list(beta0 = c(qnorm(0.1), 0.3), beta1 = c(qnorm(0.9), -0.2), pi = 0.05)
  out  <- implied_covariates(p, X, "symmetric")
  expect_named(out, c(
    "mean_entry_rate", "mean_exit_rate", "mean_employment_rate",
    "entry_flow", "exit_flow", "total_churn_flow",
    "entry_p10", "entry_median", "entry_p90",
    "exit_p10", "exit_median", "exit_p90",
    "contract_exit_effect", "informal_exit_effect", "alpha", "pi",
    "ame_entry", "ame_exit"
  ))
})

test_that("implied_covariates mean_entry_rate is mean(pnorm(Xbeta0))", {
  set.seed(1)
  X    <- cbind(1, rnorm(50))
  beta0 <- c(-1.0, 0.5)
  p    <- list(beta0 = beta0, beta1 = c(1.5, -0.3), pi = 0.05)
  out  <- implied_covariates(p, X, "symmetric")
  expected <- mean(pnorm(X %*% beta0))
  expect_equal(out$mean_entry_rate, expected, tolerance = 1e-12)
})

test_that("implied_covariates mean_exit_rate is mean(1 - pnorm(Xbeta1))", {
  set.seed(2)
  X    <- cbind(1, rnorm(50))
  beta1 <- c(1.5, -0.2)
  p    <- list(beta0 = c(-1.0, 0.5), beta1 = beta1, pi = 0.05)
  out  <- implied_covariates(p, X, "symmetric")
  expected <- mean(1 - pnorm(X %*% beta1))
  expect_equal(out$mean_exit_rate, expected, tolerance = 1e-12)
})

test_that("implied_covariates AME intercept sign correct for entry rate", {
  # Positive beta0[1] -> positive AME on entry (use index, not name)
  X    <- cbind(intercept = 1, x2 = rnorm(50, sd = 0))
  p    <- list(beta0 = c(1.0, 0.0), beta1 = c(1.5, 0.0), pi = 0.05)
  out  <- implied_covariates(p, X, "symmetric")
  expect_gt(out$ame_entry[["intercept"]], 0)
})

test_that("implied_covariates AME on exit rate is negative AME on theta1", {
  set.seed(3)
  X    <- cbind(1, rnorm(30))
  p    <- list(beta0 = c(-1.0, 0.5), beta1 = c(1.5, -0.3), pi = 0.05)
  out  <- implied_covariates(p, X, "symmetric")
  # AME exit = -AME theta1; AME of x2 on theta1 should be negative (beta1[2] < 0)
  # so AME exit for x2 should be positive
  expect_gt(out$ame_exit[2], 0)
})

test_that("implied_covariates AME names match colnames(X)", {
  X        <- cbind(intercept = 1, age = rnorm(20), educ = rnorm(20))
  p        <- list(beta0 = c(-1, 0.2, 0.1), beta1 = c(1.5, -0.1, 0.05), pi = 0.05)
  out      <- implied_covariates(p, X, "symmetric")
  expect_named(out$ame_entry, c("intercept", "age", "educ"))
  expect_named(out$ame_exit,  c("intercept", "age", "educ"))
})

test_that("implied_covariates pi is NA for no-error model", {
  X   <- cbind(1, rnorm(20))
  p   <- list(beta0 = c(-1, 0.2), beta1 = c(1.5, -0.1))
  out <- implied_covariates(p, X, "none")
  expect_true(is.na(out$pi))
})

test_that("implied_covariates errors on mismatched X and beta0 dimensions", {
  X <- cbind(1, rnorm(20))
  p <- list(beta0 = c(-1, 0.2, 0.05), beta1 = c(1.5, -0.1, 0.0), pi = 0.05)
  # X has 2 cols but betas have 3 elements
  expect_error(implied_covariates(p, X, "symmetric"))
})

test_that("implied_covariates mean rates are in (0,1)", {
  set.seed(4)
  X   <- cbind(1, rnorm(200), rnorm(200))
  p   <- list(beta0 = c(-1.0, 0.4, -0.2), beta1 = c(1.5, -0.3, 0.1), pi = 0.06)
  out <- implied_covariates(p, X, "symmetric")
  expect_gt(out$mean_entry_rate, 0); expect_lt(out$mean_entry_rate, 1)
  expect_gt(out$mean_exit_rate,  0); expect_lt(out$mean_exit_rate,  1)
  expect_gt(out$mean_employment_rate, 0); expect_lt(out$mean_employment_rate, 1)
})

# ---------------------------------------------------------------------------
# implied_fmm()
# ---------------------------------------------------------------------------

test_that("implied_fmm returns correct names", {
  p   <- list(theta0_A = 0.15, theta1_A = 0.95,
              theta0_B = 0.05, theta1_B = 0.70, phi = 0.4, pi = 0.05)
  out <- implied_fmm(p, "symmetric")
  expect_named(out, c("entry_rate_A", "exit_rate_A", "employment_rate_A",
                       "entry_rate_B", "exit_rate_B", "employment_rate_B",
                       "phi", "weighted_entry_rate", "weighted_exit_rate", "pi"))
})

test_that("implied_fmm exit rates are 1 - theta1", {
  p   <- list(theta0_A = 0.15, theta1_A = 0.95,
              theta0_B = 0.05, theta1_B = 0.70, phi = 0.4, pi = 0.05)
  out <- implied_fmm(p, "symmetric")
  expect_equal(out$exit_rate_A, 0.05)
  expect_equal(out$exit_rate_B, 0.30)
})

test_that("implied_fmm weighted_entry_rate is phi-weighted average", {
  phi <- 0.4
  p   <- list(theta0_A = 0.15, theta1_A = 0.95,
              theta0_B = 0.05, theta1_B = 0.70, phi = phi, pi = 0.05)
  out <- implied_fmm(p, "symmetric")
  expected <- phi * 0.15 + (1 - phi) * 0.05
  expect_equal(out$weighted_entry_rate, expected, tolerance = 1e-12)
})

test_that("implied_fmm weighted_exit_rate is phi-weighted average", {
  phi <- 0.4
  p   <- list(theta0_A = 0.15, theta1_A = 0.95,
              theta0_B = 0.05, theta1_B = 0.70, phi = phi, pi = 0.05)
  out <- implied_fmm(p, "symmetric")
  expected <- phi * (1 - 0.95) + (1 - phi) * (1 - 0.70)
  expect_equal(out$weighted_exit_rate, expected, tolerance = 1e-12)
})

test_that("implied_fmm employment rates use stationarity formula", {
  p   <- list(theta0_A = 0.10, theta1_A = 0.90,
              theta0_B = 0.20, theta1_B = 0.80, phi = 0.5, pi = 0.05)
  out <- implied_fmm(p, "symmetric")
  expect_equal(out$employment_rate_A, 0.10 / (0.10 + 0.10), tolerance = 1e-12)
  expect_equal(out$employment_rate_B, 0.20 / (0.20 + 0.20), tolerance = 1e-12)
})

test_that("implied_fmm pi is NA for no-error model", {
  p   <- list(theta0_A = 0.15, theta1_A = 0.95,
              theta0_B = 0.05, theta1_B = 0.70, phi = 0.4)
  out <- implied_fmm(p, "none")
  expect_true(is.na(out$pi))
})

test_that("implied_fmm errors on missing required params", {
  p <- list(theta0_A = 0.15, theta1_A = 0.95)  # missing B-type and phi
  expect_error(implied_fmm(p, "symmetric"))
})

# ---------------------------------------------------------------------------
# implied_inconsistency()
# ---------------------------------------------------------------------------

test_that("implied_inconsistency returns correct names", {
  set.seed(5)
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(rbinom(60, 1, 0.1), ncol = 6)
  out  <- implied_inconsistency(p, imat)
  expect_named(out, c("entry_rate", "exit_rate", "employment_rate",
                       "mean_pi", "delta", "pi_base",
                       "pi_age_additional", "pi_edu_additional"))
})

test_that("implied_inconsistency entry/exit rates match theta0/theta1", {
  set.seed(6)
  p    <- list(theta0 = 0.12, theta1 = 0.88, alpha = 0.5,
               delta = c(-2.0, 0.3, 0.2))
  imat <- matrix(rbinom(120, 1, 0.1), ncol = 6)
  out  <- implied_inconsistency(p, imat)
  expect_equal(out$entry_rate, 0.12)
  expect_equal(out$exit_rate,  0.12)
})

test_that("implied_inconsistency mean_pi in (0, 0.5)", {
  set.seed(7)
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(rbinom(300, 1, 0.2), ncol = 6)
  out  <- implied_inconsistency(p, imat)
  expect_gt(out$mean_pi, 0)
  expect_lt(out$mean_pi, 0.5)
})

test_that("implied_inconsistency pi_base is 0.5 * plogis(delta[1])", {
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(0L, nrow = 5, ncol = 6)   # all zeros: only intercept matters
  out  <- implied_inconsistency(p, imat)
  expected_pi_base <- 0.5 * plogis(-2.2)
  expect_equal(out$pi_base,  expected_pi_base, tolerance = 1e-12)
  # With all zeros, mean_pi should equal pi_base
  expect_equal(out$mean_pi, expected_pi_base, tolerance = 1e-12)
})

test_that("implied_inconsistency errors on wrong ncol of incons_mat", {
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(0L, nrow = 10, ncol = 3)  # wrong: need 6 cols
  expect_error(implied_inconsistency(p, imat))
})

test_that("implied_inconsistency errors on wrong delta length", {
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5))  # length 2, not 3
  imat <- matrix(0L, nrow = 5, ncol = 6)
  expect_error(implied_inconsistency(p, imat))
})

test_that("implied_inconsistency errors on NA values in incons_mat", {
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(0L, nrow = 5, ncol = 6)
  imat[1, 1] <- NA_integer_
  expect_error(implied_inconsistency(p, imat), "NA values")
})

test_that("implied_inconsistency errors on non-binary values in incons_mat", {
  p    <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  imat <- matrix(0L, nrow = 5, ncol = 6)
  imat[2, 3] <- 2L   # out-of-range
  expect_error(implied_inconsistency(p, imat), "0 or 1")
})

# ---------------------------------------------------------------------------
# implied_fmm() boundary tests (phi = 0 or 1)
# ---------------------------------------------------------------------------

test_that("implied_fmm with phi=1: weighted rates equal type-A rates", {
  p <- list(theta0_A = 0.20, theta1_A = 0.90,
            theta0_B = 0.05, theta1_B = 0.70, phi = 1.0, pi = NA_real_)
  out <- implied_fmm(p, "none")
  expect_equal(out$weighted_entry_rate, out$entry_rate_A)
  expect_equal(out$weighted_exit_rate,  out$exit_rate_A)
})

test_that("implied_fmm with phi=0: weighted rates equal type-B rates", {
  p <- list(theta0_A = 0.20, theta1_A = 0.90,
            theta0_B = 0.05, theta1_B = 0.70, phi = 0.0, pi = NA_real_)
  out <- implied_fmm(p, "none")
  expect_equal(out$weighted_entry_rate, out$entry_rate_B)
  expect_equal(out$weighted_exit_rate,  out$exit_rate_B)
})

# ---------------------------------------------------------------------------
# implied_covariates() anyNA guard
# ---------------------------------------------------------------------------

test_that("implied_covariates errors when X contains NA", {
  X     <- cbind(1, rnorm(10))
  X[3, 2] <- NA_real_
  params <- list(beta0 = c(qnorm(0.1), 0.5), beta1 = c(qnorm(0.9), -0.2), pi = 0.05)
  expect_error(implied_covariates(params, X, "symmetric"), "NA values")
})
