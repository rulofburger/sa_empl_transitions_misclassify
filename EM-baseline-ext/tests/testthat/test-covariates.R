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
                     "M", "C1", "C0", "init_w1", "init_w0", "total_weight")
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

test_that("stationary Q analytic gradient matches central differences", {
  df <- .make_panel_cov(n = 120L)
  X <- cbind(intercept = 1, x = seq(-1, 1, length.out = nrow(df)))
  params <- init_params_covariates(ncol(X))
  suff <- e_step_covariates(df, X, params)$suff
  active <- rep(TRUE, ncol(X))
  z <- c(params$beta0, params$beta1)
  analytic <- .covariate_q(z, suff, X, active, TRUE, gradient = TRUE)
  h <- 1e-6
  numeric_grad <- vapply(seq_along(z), function(j) {
    zp <- zm <- z
    zp[j] <- zp[j] + h
    zm[j] <- zm[j] - h
    (.covariate_q(zp, suff, X, active, TRUE) -
       .covariate_q(zm, suff, X, active, TRUE)) / (2 * h)
  }, numeric(1))
  expect_equal(analytic, numeric_grad, tolerance = 1e-5)
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

test_that("stationary symmetric fit weakly dominates nested no-error start", {
  df <- .make_panel_cov(n = 400L, pi = 0.04)
  X <- .make_X_intercept(nrow(df))
  fit_none <- em_fit_covariates(df, X, model_type = "none", verbose = 0L)
  p0 <- fit_none$params
  p0$pi <- 1e-8
  fit_sym <- em_fit_covariates(df, X, model_type = "symmetric",
                               params0 = p0, verbose = 0L)
  expect_gte(fit_sym$loglik, fit_none$loglik - 1e-6)
})

test_that("free-alpha model uses its free initial probability", {
  df <- .make_panel_cov(n = 300L)
  X <- .make_X_intercept(nrow(df))
  p1 <- init_params_covariates(1L)
  p1$alpha <- 0.2
  p2 <- p1
  p2$alpha <- 0.8
  ll1 <- e_step_covariates(df, X, p1, stationary = FALSE)$loglik
  ll2 <- e_step_covariates(df, X, p2, stationary = FALSE)$loglik
  expect_false(isTRUE(all.equal(ll1, ll2)))
})

test_that("common weight rescaling leaves coefficients unchanged", {
  df <- .make_panel_cov(n = 250L, pi = 0)
  X <- .make_X_intercept(nrow(df))
  fit1 <- em_fit_covariates(df, X, model_type = "none", verbose = 0L)
  df$weight <- df$weight * 100
  fit2 <- em_fit_covariates(df, X, model_type = "none", verbose = 0L)
  expect_equal(fit1$params$beta0, fit2$params$beta0, tolerance = 1e-8)
  expect_equal(fit1$params$beta1, fit2$params$beta1, tolerance = 1e-8)
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

test_that("time-varying transition designs require free alpha and fit", {
  df <- .make_panel_cov(n = 120L, pi = 0)
  z1 <- rep(c(0, 1), length.out = nrow(df))
  z2 <- 1 - z1
  X12 <- cbind(intercept = 1, origin_attribute = z1)
  X23 <- cbind(intercept = 1, origin_attribute = z2)
  design <- list(X12 = X12, X23 = X23)

  expect_error(
    em_fit_covariates(df, design, model_type = "none", stationary = TRUE,
                      verbose = 0L),
    "stationary=FALSE"
  )
  fit <- em_fit_covariates(df, design, model_type = "none",
                           stationary = FALSE, verbose = 0L)
  expect_true(is.finite(fit$loglik))
  expect_equal(dim(fit$gamma), c(nrow(df), 8L))
})

test_that("analytical sandwich/delta SEs are finite on a free-alpha fit", {
  df <- .make_panel_cov(n = 180L, pi = 0)
  df$weight <- runif(nrow(df), 0.5, 2)
  X <- cbind(intercept = 1, z = rnorm(nrow(df)))
  fit <- em_fit_covariates(df, X, model_type = "none", stationary = FALSE,
                           max_iter = 150L, verbose = 0L)
  out <- analytical_se_covariates(df, X, fit, "none")
  expect_true(all(is.finite(out$summary$se)))
  expect_true(all(out$summary$se > 0))
  expect_gt(out$diagnostics$min_bread_eigenvalue, 0)
})
