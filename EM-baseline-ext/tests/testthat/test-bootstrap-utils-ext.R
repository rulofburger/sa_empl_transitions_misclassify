# ==============================================================================
# EM-baseline-ext: Tests for bootstrap_utils_ext.R
# Created: 2026-05-06
# ==============================================================================

# ---------------------------------------------------------------------------
# summarise_bootstrap_ame()
# ---------------------------------------------------------------------------

test_that("summarise_bootstrap_ame returns correct column names", {
  set.seed(1L)
  ame_e <- c(intercept = 0.02, age = -0.01, educ = 0.005)
  ame_x <- c(intercept = -0.01, age = 0.005, educ = -0.002)
  make_rep <- function(i) {
    set.seed(i)
    imp <- list(ame_entry = ame_e * (1 + runif(1, -0.01, 0.01)),
                ame_exit  = ame_x * (1 + runif(1, -0.01, 0.01)))
    list(params = list(beta0 = c(1, 2, 3), beta1 = c(1, 2, 3)),
         implied = imp, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- lapply(1:10, make_rep)
  out   <- summarise_bootstrap_ame(boots, ame_e, ame_x)
  expect_true(all(c("covariate", "outcome", "estimate", "se",
                    "ci_lower", "ci_upper", "n_ok", "n_reps") %in% names(out)))
})

test_that("summarise_bootstrap_ame has one row per covariate per outcome", {
  ame_e <- c(x1 = 0.02, x2 = -0.01)
  ame_x <- c(x1 = -0.01, x2 = 0.005)
  make_rep <- function() {
    imp <- list(ame_entry = ame_e, ame_exit = ame_x)
    list(params = list(), implied = imp, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- replicate(5, make_rep(), simplify = FALSE)
  out   <- summarise_bootstrap_ame(boots, ame_e, ame_x)
  # 2 covariates × 2 outcomes = 4 rows
  expect_equal(nrow(out), 4L)
  expect_equal(sort(unique(out$outcome)), c("entry", "exit"))
  expect_equal(sort(unique(out$covariate)), c("x1", "x2"))
})

test_that("summarise_bootstrap_ame SE is 0 for constant replicates", {
  ame_e <- c(x1 = 0.02, x2 = -0.01)
  ame_x <- c(x1 = -0.01, x2 = 0.005)
  make_rep <- function() {
    imp <- list(ame_entry = ame_e, ame_exit = ame_x)
    list(params = list(), implied = imp, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- replicate(10, make_rep(), simplify = FALSE)
  out   <- summarise_bootstrap_ame(boots, ame_e, ame_x)
  expect_true(all(out$se == 0))
})

test_that("summarise_bootstrap_ame excludes non-ok reps from SE", {
  set.seed(2L)
  ame_e <- c(x1 = 0.02)
  ame_x <- c(x1 = -0.01)
  make_rep <- function(i, flag) {
    set.seed(i)
    imp <- list(ame_entry = ame_e * (1 + runif(1, -0.01, 0.01)),
                ame_exit  = ame_x)
    list(params = list(), implied = imp, loglik = -100, converged = TRUE, flag = flag)
  }
  boots <- c(lapply(1:8,  function(i) make_rep(i, "ok")),
             lapply(9:10, function(i) make_rep(i, "no_converge")))
  out   <- summarise_bootstrap_ame(boots, ame_e, ame_x)
  expect_equal(unique(out$n_ok),   8L)
  expect_equal(unique(out$n_reps), 10L)
})

test_that("summarise_bootstrap_ame handles all-failed boots", {
  ame_e <- c(x1 = 0.02)
  ame_x <- c(x1 = -0.01)
  boots <- replicate(5, list(params = NULL, implied = NULL, loglik = NA,
                             converged = FALSE, flag = "error"),
                     simplify = FALSE)
  out   <- summarise_bootstrap_ame(boots, ame_e, ame_x)
  expect_true(is.data.frame(out))
  expect_true("covariate" %in% names(out))
  # 1 covariate × 2 outcomes = 2 rows; all SEs NA; n_ok = 0
  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$n_ok), 0L)
  expect_true(all(is.na(out$se)))
})

# ---------------------------------------------------------------------------
# bootstrap_one_covariates() — synthetic data smoke test
# ---------------------------------------------------------------------------

test_that("bootstrap_one_covariates returns expected structure", {
  skip_if_not(exists("em_fit_covariates"),
              "em_fit_covariates not loaded — skip in isolated test runs")
  set.seed(7L)
  N  <- 150L
  df <- data.frame(y1 = rbinom(N, 1, 0.5), y2 = rbinom(N, 1, 0.5),
                   y3 = rbinom(N, 1, 0.5), weight = rep(1, N))
  X  <- cbind(1, rnorm(N))
  p0 <- list(beta0 = c(qnorm(0.1), 0.1), beta1 = c(qnorm(0.9), -0.1), pi = 0.05)
  res <- bootstrap_one_covariates(df, X, seed = 1L, model_type = "symmetric",
                                   stationary = TRUE, params_start = p0,
                                   point_loglik = -1e6)
  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "error"))
  expect_false(res$flag == "low_loglik")  # point_loglik=-1e6 is very low
})

test_that("bootstrap_one_covariates is reproducible with same seed", {
  skip_if_not(exists("em_fit_covariates"),
              "em_fit_covariates not loaded — skip in isolated test runs")
  set.seed(7L)
  N  <- 150L
  df <- data.frame(y1 = rbinom(N, 1, 0.5), y2 = rbinom(N, 1, 0.5),
                   y3 = rbinom(N, 1, 0.5), weight = rep(1, N))
  X  <- cbind(1, rnorm(N))
  p0 <- list(beta0 = c(qnorm(0.1), 0.1), beta1 = c(qnorm(0.9), -0.1), pi = 0.05)
  r1 <- bootstrap_one_covariates(df, X, seed = 42L, model_type = "symmetric",
                                  stationary = TRUE, params_start = p0,
                                  point_loglik = -1e6)
  r2 <- bootstrap_one_covariates(df, X, seed = 42L, model_type = "symmetric",
                                  stationary = TRUE, params_start = p0,
                                  point_loglik = -1e6)
  expect_identical(r1$loglik, r2$loglik)
})

test_that("bootstrap_one_covariates returns flag='error' on bad df", {
  # Passing a df without required columns causes bootstrap_resample to error;
  # bootstrap_one_* does not swallow this (it is a programming error, not data).
  bad_df <- data.frame(a = 1:5, b = 1:5)
  X      <- cbind(1, 1:5)
  p0     <- list(beta0 = c(0, 0), beta1 = c(0, 0))
  expect_error(
    bootstrap_one_covariates(bad_df, X, seed = 1L, model_type = "symmetric",
                              stationary = TRUE, params_start = p0,
                              point_loglik = -100),
    "missing required columns"
  )
})

# ---------------------------------------------------------------------------
# bootstrap_one_fmm() — synthetic data smoke test
# ---------------------------------------------------------------------------

test_that("bootstrap_one_fmm returns expected structure", {
  skip_if_not(exists("em_fit_fmm"),
              "em_fit_fmm not loaded — skip in isolated test runs")
  set.seed(7L)
  N  <- 150L
  df <- data.frame(y1 = rbinom(N, 1, 0.5), y2 = rbinom(N, 1, 0.5),
                   y3 = rbinom(N, 1, 0.5), weight = rep(1, N))
  p0 <- list(theta0_A = 0.20, theta1_A = 0.85,
             theta0_B = 0.05, theta1_B = 0.70, phi = 0.4, pi = 0.05)
  res <- bootstrap_one_fmm(df, seed = 2L, model_type = "symmetric",
                             stationary = TRUE, params_start = p0,
                             point_loglik = -1e6)
  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "error"))
})

test_that("bootstrap_one_fmm returns flag='error' on bad df", {
  bad_df <- data.frame(a = 1:5)
  p0 <- list(theta0_A = 0.2, theta1_A = 0.9, theta0_B = 0.05, theta1_B = 0.7,
             phi = 0.4, pi = 0.05)
  expect_error(
    bootstrap_one_fmm(bad_df, seed = 1L, model_type = "symmetric",
                       stationary = TRUE, params_start = p0,
                       point_loglik = -100),
    "missing required columns"
  )
})

# ---------------------------------------------------------------------------
# bootstrap_one_inconsistency() — synthetic data smoke test
# ---------------------------------------------------------------------------

test_that("bootstrap_one_inconsistency returns expected structure", {
  skip_if_not(exists("em_fit_inconsistency"),
              "em_fit_inconsistency not loaded — skip in isolated test runs")
  set.seed(7L)
  N  <- 150L
  df <- data.frame(y1 = rbinom(N, 1, 0.5), y2 = rbinom(N, 1, 0.5),
                   y3 = rbinom(N, 1, 0.5), weight = rep(1, N))
  imat <- matrix(rbinom(N * 6, 1, 0.1), nrow = N, ncol = 6L)
  p0   <- list(theta0 = 0.10, theta1 = 0.90, alpha = 0.5,
               delta = c(-2.2, 0.5, 0.3))
  res  <- bootstrap_one_inconsistency(df, imat, seed = 3L, stationary = TRUE,
                                       params_start = p0, point_loglik = -1e6)
  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "error"))
})

test_that("bootstrap_one_inconsistency returns flag='error' on bad df", {
  bad_df <- data.frame(a = 1:5)
  imat   <- matrix(0L, nrow = 5, ncol = 6)
  p0     <- list(theta0 = 0.10, theta1 = 0.90, delta = c(-2, 0.5, 0.3))
  expect_error(
    bootstrap_one_inconsistency(bad_df, imat, seed = 1L, stationary = TRUE,
                                 params_start = p0, point_loglik = -100),
    "missing required columns"
  )
})

# ---------------------------------------------------------------------------
# .flag_fit() — quality flag logic (accessible via bootstrap_one_baseline env)
# ---------------------------------------------------------------------------

test_that(".flag_fit returns 'ok' when converged and loglik acceptable", {
  fit <- list(converged = TRUE, loglik = -100)
  expect_equal(.flag_fit(fit, point_loglik = -100), "ok")
})

test_that(".flag_fit returns 'no_converge' when converged=FALSE", {
  fit <- list(converged = FALSE, loglik = -100)
  expect_equal(.flag_fit(fit, point_loglik = -100), "no_converge")
})

test_that(".flag_fit returns 'low_loglik' when loglik drop > threshold", {
  fit <- list(converged = TRUE, loglik = -160)   # -100 - 51 < -100 - 50
  expect_equal(.flag_fit(fit, point_loglik = -100), "low_loglik")
})

test_that(".flag_fit returns 'ok' when loglik drop exactly at threshold", {
  fit <- list(converged = TRUE, loglik = -150)   # -100 - 50 = exactly threshold
  expect_equal(.flag_fit(fit, point_loglik = -100), "ok")
})

test_that(".flag_fit skips loglik check when point_loglik is NA", {
  fit <- list(converged = TRUE, loglik = -1e6)
  expect_equal(.flag_fit(fit, point_loglik = NA_real_), "ok")
})
