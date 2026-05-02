# test-rho-model.R
# Tests for rho-augmented duration contamination model.
#
# Created: 2026-04-29
# TeX ref: "EM tenure rho.tex"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.make_rho_data <- function(n = 200L, seed = 42L) {
  set.seed(seed)
  s <- sample(0:1, n, replace = TRUE, prob = c(0.35, 0.65))
  list(
    y1          = s,
    y2          = s,
    y3          = s,
    tenure1     = pmax(0.26, s * rexp(n, 2)),
    tenure2     = pmax(0.26, s * rexp(n, 2) + 0.25),
    tenure3     = pmax(0.26, s * rexp(n, 2) + 0.50),
    timegap_cat1 = sample(1:7, n, replace = TRUE),
    timegap_cat2 = sample(1:7, n, replace = TRUE),
    timegap_cat3 = sample(1:7, n, replace = TRUE),
    weight      = rep(1, n)
  )
}

.make_rho_params <- function(rho = 0.1) {
  list(
    alpha    = 0.6,
    theta1   = 0.85,
    theta0   = 0.15,
    pi       = 0.05,
    rho      = rho,
    sigma2_g = 0.01,
    lambda_g = 2.0,
    lambda_d = 1.5
  )
}

# ---------------------------------------------------------------------------
# 1. .log_mix_rho: numerical stability at extremes
# ---------------------------------------------------------------------------

test_that(".log_mix_rho is numerically stable at extremes", {
  # large positive log_clock >> log_pop -> mixture ~ log_clock
  log_clock <- 1000
  log_pop   <- -1000
  res <- .log_mix_rho(log_clock, log_pop, rho = 0.1)
  expect_true(is.finite(res))

  # large positive log_pop >> log_clock -> mixture ~ log_pop
  log_clock2 <- -1000
  log_pop2   <- 1000
  res2 <- .log_mix_rho(log_clock2, log_pop2, rho = 0.9)
  expect_true(is.finite(res2))

  # equal inputs: mixture log[(1-rho)*exp(a) + rho*exp(a)] = a exactly
  a <- -2; b <- -2
  res3 <- .log_mix_rho(a, b, rho = 0.5)
  expect_equal(res3, a, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. .omega_rho: boundary behaviour
# ---------------------------------------------------------------------------

test_that(".omega_rho is in [0, 1] and approaches limits correctly", {
  # when clock >> pop, omega (contamination prob) -> 0
  omega_lo <- .omega_rho(log_clock = 100, log_pop = -100, rho = 0.1)
  expect_lt(omega_lo, 0.01)

  # when pop >> clock, omega -> 1
  omega_hi <- .omega_rho(log_clock = -100, log_pop = 100, rho = 0.1)
  expect_gt(omega_hi, 0.99)

  # always in [0, 1]
  set.seed(1)
  lc <- rnorm(50, mean = -2, sd = 2)
  lp <- rnorm(50, mean = -2, sd = 2)
  omegas <- .omega_rho(lc, lp, rho = 0.2)
  expect_true(all(omegas >= 0) && all(omegas <= 1))
})

# ---------------------------------------------------------------------------
# 3. e_step_rho: output structure
# ---------------------------------------------------------------------------

test_that("e_step_rho returns correct structure", {
  df     <- as.data.frame(.make_rho_data())
  params <- .make_rho_params()

  out <- e_step_rho(df, params)

  expect_named(out, c("gamma", "loglik", "suff"))
  expect_equal(dim(out$gamma), c(200L, 8L))
  expect_true(is.finite(out$loglik))
  expect_true(!is.null(out$suff$Omega))
  expect_true(is.finite(out$suff$Omega))
  expect_true(out$suff$Omega >= 0)
})

# ---------------------------------------------------------------------------
# 4. e_step_rho: row sums of gamma are 1
# ---------------------------------------------------------------------------

test_that("e_step_rho gamma rows sum to 1", {
  df     <- as.data.frame(.make_rho_data())
  params <- .make_rho_params()

  out <- e_step_rho(df, params)
  row_sums <- rowSums(out$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

# ---------------------------------------------------------------------------
# 5. e_step_rho: input validation — bad rho
# ---------------------------------------------------------------------------

test_that("e_step_rho rejects out-of-range rho", {
  df     <- as.data.frame(.make_rho_data())

  expect_error(e_step_rho(df, .make_rho_params(rho = 0)),    "params\\$rho")
  expect_error(e_step_rho(df, .make_rho_params(rho = 1)),    "params\\$rho")
  expect_error(e_step_rho(df, .make_rho_params(rho = -0.1)), "params\\$rho")
  expect_error(e_step_rho(df, .make_rho_params(rho = 1.5)),  "params\\$rho")
})

# ---------------------------------------------------------------------------
# 6. e_step_rho: input validation — bad timegap_cat
# ---------------------------------------------------------------------------

test_that("e_step_rho rejects out-of-range timegap_cat", {
  df <- as.data.frame(.make_rho_data())
  df$timegap_cat1[1] <- 99L
  expect_error(e_step_rho(df, .make_rho_params()), "timegap_cat")
})

# ---------------------------------------------------------------------------
# 7. m_step_rho: rho update closed-form check
# ---------------------------------------------------------------------------

test_that("m_step_rho rho update uses Omega / (3 * total_weight)", {
  suff <- list(
    C1 = 50, C0 = 50,
    D1 = 80, D0 = 80,
    T11 = 60, T01 = 20,
    M  = 5,
    Omega = 90,
    Sg = 0.5, Ng = 100,
    Sg_start = 0.1, Ng_start = 20,
    emg_g_x          = numeric(0),
    emg_g_w          = numeric(0),
    cat_d_marginal_c  = integer(0),
    cat_d_marginal_w  = numeric(0),
    cat_d_trans_curr  = integer(0),
    cat_d_trans_prev  = integer(0),
    cat_d_trans_w     = numeric(0)
  )
  total_weight <- 150

  out <- m_step_rho(suff, total_weight, stationary = FALSE, linked = FALSE,
                    rho_cap = 0.49)

  expected_rho <- min(max(90 / (3 * 150), 1e-4), 0.49)
  expect_equal(out$rho, expected_rho, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 8. LL monotonicity: em_fit_tenure_rho never decreases log-likelihood
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_rho log-likelihood is non-decreasing", {
  df <- as.data.frame(.make_rho_data(n = 150, seed = 7))

  fit <- em_fit_tenure_rho(df, max_iter = 30L, tol = 1e-6, verbose = 0L)

  ll_hist <- fit$history$loglik
  ll_hist <- ll_hist[is.finite(ll_hist)]
  diffs   <- diff(ll_hist)
  # Allow tiny floating-point drops but no real decreases
  expect_true(all(diffs >= -1e-6),
              info = paste("LL decreased at iter:", which(diffs < -1e-6)))
})

# ---------------------------------------------------------------------------
# 9. Convergence: em_fit_tenure_rho converges on clean data
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_rho converges in reasonable iterations", {
  df  <- as.data.frame(.make_rho_data(n = 200, seed = 1))
  fit <- em_fit_tenure_rho(df, max_iter = 200L, tol = 1e-6, verbose = 0L)

  expect_true(fit$converged)
  expect_lt(fit$iterations, 200L)
  expect_true(is.finite(fit$loglik))
})

# ---------------------------------------------------------------------------
# 10. init_params_rho: rho_init out of range is rejected
# ---------------------------------------------------------------------------

test_that("init_params_rho rejects invalid rho_init", {
  df <- as.data.frame(.make_rho_data())
  expect_error(init_params_rho(df, rho_init = 0),    "rho_init")
  expect_error(init_params_rho(df, rho_init = 1),    "rho_init")
  expect_error(init_params_rho(df, rho_init = -0.1), "rho_init")
})

# ---------------------------------------------------------------------------
# 11. stationary=TRUE: alpha respects stationarity formula
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_rho stationary=TRUE converges and alpha ~ theta0/(theta0 + 1 - theta1)", {
  df  <- as.data.frame(.make_rho_data(n = 200, seed = 3))
  fit <- em_fit_tenure_rho(df, max_iter = 100L, tol = 1e-5, stationary = TRUE,
                           verbose = 0L)

  p <- fit$params
  alpha_implied <- p$theta0 / (p$theta0 + 1 - p$theta1)
  expect_equal(p$alpha, alpha_implied, tolerance = 1e-8)
  expect_true(is.finite(fit$loglik))
})

# ---------------------------------------------------------------------------
# 12. linked=TRUE: lambda_g and lambda_d are CTMC-derived from theta
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_rho linked=TRUE lambda_g and lambda_d match CTMC link", {
  df  <- as.data.frame(.make_rho_data(n = 200, seed = 5))
  fit <- em_fit_tenure_rho(df, max_iter = 100L, tol = 1e-5, linked = TRUE,
                           verbose = 0L)

  p <- fit$params
  # CTMC link: lambda = -log(theta) / (1/4 year interval)
  lambda_g_expected <- ctmc_lambda_from_persistence(p$theta1)
  lambda_d_expected <- ctmc_lambda_from_transition(p$theta0)
  expect_equal(p$lambda_g, lambda_g_expected, tolerance = 1e-8)
  expect_equal(p$lambda_d, lambda_d_expected, tolerance = 1e-8)
  expect_true(is.finite(fit$loglik))
})
