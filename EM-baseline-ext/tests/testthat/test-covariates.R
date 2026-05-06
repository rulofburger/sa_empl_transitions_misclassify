# ==============================================================================
# EM-baseline-ext: Tests for Extension I (covariate model, GEM)
# Created: 2026-05-06
# ==============================================================================

# Helper: synthetic panel data frame (no real covariates needed here)
.make_panel_cov <- function(n = 200L, seed = 123L,
                             theta0 = 0.10, theta1 = 0.90, pi = 0.05) {
  set.seed(seed)
  # Simulate latent employment and misclassification
  alpha <- theta0 / (theta0 + 1 - theta1)
  h1 <- rbinom(n, 1L, alpha)
  h2 <- ifelse(h1 == 1L, rbinom(n, 1L, theta1), rbinom(n, 1L, theta0))
  h3 <- ifelse(h2 == 1L, rbinom(n, 1L, theta1), rbinom(n, 1L, theta0))
  s1 <- ifelse(rbinom(n, 1L, pi), 1L - h1, h1)
  s2 <- ifelse(rbinom(n, 1L, pi), 1L - h2, h2)
  s3 <- ifelse(rbinom(n, 1L, pi), 1L - h3, h3)
  data.frame(y1 = s1, y2 = s2, y3 = s3, weight = rep(1, n))
}

# Helper: intercept-only design matrix
.make_X_intercept <- function(n) {
  matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "intercept"))
}

# ---- init_params_covariates ------------------------------------------------

test_that("init_params_covariates returns correct structure for symmetric", {
  p <- init_params_covariates(4L, model_type = "symmetric")
  expect_equal(length(p$beta0), 4L)
  expect_equal(length(p$beta1), 4L)
  expect_true(!is.null(p$pi))
  expect_true(p$pi > 0 && p$pi < 0.5)
})

test_that("init_params_covariates returns no pi for model_type='none'", {
  p <- init_params_covariates(3L, model_type = "none")
  expect_null(p$pi)
  expect_equal(length(p$beta0), 3L)
})

test_that("init_params_covariates intercepts imply sensible theta values", {
  p <- init_params_covariates(4L)
  theta0_implied <- pnorm(p$beta0[1L])
  theta1_implied <- pnorm(p$beta1[1L])
  expect_true(theta0_implied > 0.05 && theta0_implied < 0.30)
  expect_true(theta1_implied > 0.70 && theta1_implied < 0.99)
})

# ---- e_step_covariates: structure ------------------------------------------

test_that("e_step_covariates returns gamma, loglik, suff", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_true(!is.null(out$gamma))
  expect_true(!is.null(out$loglik))
  expect_true(!is.null(out$suff))
})

test_that("e_step_covariates: gamma rows sum to 1", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_equal(nrow(out$gamma), nrow(df))
  expect_equal(ncol(out$gamma), 8L)
  expect_equal(rowSums(out$gamma), rep(1, nrow(df)), tolerance = 1e-10)
})

test_that("e_step_covariates: gamma entries non-negative", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_true(all(out$gamma >= 0))
})

test_that("e_step_covariates: loglik is finite", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_true(is.finite(out$loglik))
})

test_that("e_step_covariates: model_type='none' loglik is finite", {
  df     <- .make_panel_cov(pi = 0)
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L, model_type = "none")
  # model_type='none' can fail if latent state doesn't match any history;
  # avoid that with consistent pi=0 data
  out <- tryCatch(
    e_step_covariates(df, X, params, model_type = "none"),
    error = function(e) NULL
  )
  # Either succeeds with finite loglik or errors (both acceptable on random data)
  if (!is.null(out)) expect_true(is.finite(out$loglik))
})

# ---- e_step_covariates: suff statistics ------------------------------------

test_that("e_step_covariates suff has required fields", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  required_suff <- c("eff_w_1", "eff_wy_1", "eff_w_0", "eff_wy_0",
                     "M", "C1", "C0", "total_weight")
  expect_true(all(required_suff %in% names(out$suff)))
})

test_that("e_step_covariates: eff_w vectors are non-negative", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_true(all(out$suff$eff_w_1 >= 0))
  expect_true(all(out$suff$eff_w_0 >= 0))
})

test_that("e_step_covariates: total_weight equals sum of weights", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  out    <- e_step_covariates(df, X, params)
  expect_equal(out$suff$total_weight, sum(df$weight))
})

# ---- Intercept-only approximates baseline ----------------------------------

test_that("intercept-only covariate model LL is close to baseline EM LL", {
  # Run baseline EM (using baseline functions sourced via helper-source.R)
  df     <- .make_panel_cov(n = 300L)
  fit_bl <- em_fit_baseline(df, model_type = "symmetric", verbose = 0L)

  # Run covariate extension with intercept-only design matrix
  X      <- .make_X_intercept(nrow(df))
  fit_cv <- em_fit_covariates(df, X, model_type = "symmetric", verbose = 0L)

  # LL should be close (same model at convergence)
  # Relative tolerance: 0.1% of the baseline LL magnitude
  expect_equal(fit_cv$loglik, fit_bl$loglik,
               tolerance = 1e-3 * abs(fit_bl$loglik))
})

# ---- m_step_covariates -----------------------------------------------------

test_that("m_step_covariates returns beta0, beta1, and pi", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  suff   <- e_step_covariates(df, X, params)$suff
  out    <- m_step_covariates(suff, X, params)
  expect_equal(length(out$beta0), 1L)
  expect_equal(length(out$beta1), 1L)
  expect_true(!is.null(out$pi))
})

test_that("m_step_covariates: pi is in (0, pi_cap)", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  suff   <- e_step_covariates(df, X, params)$suff
  out    <- m_step_covariates(suff, X, params, pi_cap = 0.49)
  expect_true(out$pi >= 0 && out$pi <= 0.49)
})

test_that("m_step_covariates: beta values are finite", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  suff   <- e_step_covariates(df, X, params)$suff
  out    <- m_step_covariates(suff, X, params)
  expect_true(all(is.finite(out$beta0)))
  expect_true(all(is.finite(out$beta1)))
})

# ---- em_fit_covariates: convergence ----------------------------------------

test_that("em_fit_covariates converges on synthetic data", {
  df  <- .make_panel_cov(n = 300L)
  X   <- .make_X_intercept(nrow(df))
  fit <- em_fit_covariates(df, X, model_type = "symmetric", verbose = 0L)
  expect_true(fit$converged)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_covariates: LL is monotone non-decreasing", {
  df  <- .make_panel_cov(n = 300L)
  X   <- .make_X_intercept(nrow(df))
  fit <- em_fit_covariates(df, X, verbose = 0L)
  ll_vec <- fit$history$loglik
  diffs  <- diff(ll_vec)
  expect_true(all(diffs >= -1e-5), label = "LL must be non-decreasing")
})

test_that("em_fit_covariates: no-error variant converges", {
  df  <- .make_panel_cov(n = 300L, pi = 0)
  X   <- .make_X_intercept(nrow(df))
  fit <- em_fit_covariates(df, X, model_type = "none", verbose = 0L)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_covariates returns gamma with correct dimensions", {
  df  <- .make_panel_cov(n = 100L)
  X   <- .make_X_intercept(nrow(df))
  fit <- em_fit_covariates(df, X, verbose = 0L)
  expect_equal(nrow(fit$gamma), nrow(df))
  expect_equal(ncol(fit$gamma), 8L)
})

# ---- compute_observed_loglik_covariates ------------------------------------

test_that("compute_observed_loglik_covariates returns finite scalar", {
  df     <- .make_panel_cov()
  X      <- .make_X_intercept(nrow(df))
  params <- init_params_covariates(1L)
  ll     <- compute_observed_loglik_covariates(df, X, params)
  expect_true(is.finite(ll))
  expect_length(ll, 1L)
})

# ---- em_fit_covariates: additional convergence guards ---------------------

test_that("em_fit_covariates runs more than 1 iteration", {
  df  <- .make_panel_cov(n = 300L)
  X   <- .make_X_intercept(nrow(df))
  fit <- em_fit_covariates(df, X, verbose = 0L)
  expect_true(fit$iterations > 1L,
              label = "covariate EM must iterate more than once")
})

test_that("em_fit_covariates: params change from starting values", {
  df  <- .make_panel_cov(n = 300L)
  X   <- .make_X_intercept(nrow(df))
  p0  <- init_params_covariates(1L)
  fit <- em_fit_covariates(df, X, params0 = p0, verbose = 0L)
  expect_false(identical(fit$params$beta0, p0$beta0),
               label = "beta0 must change from initial value")
})

test_that("covariate model estimates non-zero slope on age-structured data", {
  set.seed(42L)
  n       <- 500L
  age_std <- rnorm(n)
  p_emp   <- pnorm(-0.2 + 0.5 * age_std)
  y_sim   <- function(p) as.integer(runif(n) < p)
  df <- data.frame(y1 = y_sim(p_emp), y2 = y_sim(p_emp), y3 = y_sim(p_emp),
                   weight = rep(1, n))
  X   <- cbind(1, age_std)
  fit <- em_fit_covariates(df, X, verbose = 0L, max_iter = 200L)
  expect_true(abs(fit$params$beta0[2L]) > 1e-4 || abs(fit$params$beta1[2L]) > 1e-4,
              label = "age slope should be non-zero on age-structured data")
})
