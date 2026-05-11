# test-implied-tenure-contamination.R
# Tests for implied_tenure_contamination() and bootstrap_one_eps()
#
# Created: 2026-05-07
# Plan ref: .cg-docs/plans/2026-05-07-tenure-contamination-bootstrap-tables.md

# Note: helper-source.R (sourced automatically by testthat) loads all eps
# model symbols. implied_quantities_tenure_contamination.R and
# bootstrap_utils_tenure_contamination.R are sourced explicitly below
# because they may not yet appear in source_all.R.

if (requireNamespace("here", quietly = TRUE)) {
  .r_dir <- here::here("EM-tenure", "R")
} else {
  .r_dir <- normalizePath("EM-tenure/R", mustWork = FALSE)
}

source(file.path(.r_dir, "implied_quantities_tenure_contamination.R"))

# bootstrap_utils_tenure_contamination requires bootstrap_resample from baseline
.baseline_r_dir <- if (requireNamespace("here", quietly = TRUE)) {
  here::here("EM-baseline", "R")
} else {
  normalizePath("EM-baseline/R", mustWork = FALSE)
}
source(file.path(.baseline_r_dir, "bootstrap_utils.R"))
source(file.path(.r_dir, "bootstrap_utils_tenure_contamination.R"))

rm(.r_dir, .baseline_r_dir)

# ---------------------------------------------------------------------------
# Helper: canonical parameter set used across tests
# ---------------------------------------------------------------------------

.p_full <- list(
  alpha    = 0.53,
  theta0   = 0.10,
  theta1   = 0.92,
  pi       = 0.088,
  eps      = 0.17,
  lambda_g = 0.16,
  lambda_d = 0.50
)

# ---------------------------------------------------------------------------
# 1. implied_tenure_contamination: basic correctness
# ---------------------------------------------------------------------------

test_that("implied_tenure_contamination: entry/exit/employment_rate correct", {
  out <- implied_tenure_contamination(.p_full)

  expect_equal(out$entry_rate, .p_full$theta0)
  expect_equal(out$exit_rate,  1 - .p_full$theta1)

  denom <- .p_full$theta0 + (1 - .p_full$theta1)
  expect_equal(out$employment_rate, .p_full$theta0 / denom, tolerance = 1e-12)
})

test_that("implied_tenure_contamination: pi and eps passed through", {
  out <- implied_tenure_contamination(.p_full)
  expect_equal(out$pi,  .p_full$pi)
  expect_equal(out$eps, .p_full$eps)
})

test_that("implied_tenure_contamination: spell duration quantities correct", {
  out <- implied_tenure_contamination(.p_full)

  expect_equal(out$mean_spell_g_years,   1  / .p_full$lambda_g, tolerance = 1e-12)
  expect_equal(out$mean_spell_g_months,  12 / .p_full$lambda_g, tolerance = 1e-12)
  expect_equal(out$mean_spell_d_years,   1  / .p_full$lambda_d, tolerance = 1e-12)
  expect_equal(out$mean_spell_d_months,  12 / .p_full$lambda_d, tolerance = 1e-12)
  expect_equal(out$median_spell_g_months,
               12 * log(2) / .p_full$lambda_g, tolerance = 1e-12)
})

test_that("implied_tenure_contamination: contamination probability quantities", {
  out <- implied_tenure_contamination(.p_full)
  eps <- .p_full$eps

  expect_equal(out$p_clock_consistent,   1    - eps,       tolerance = 1e-12)
  expect_equal(out$p_pair_clean,         (1   - eps)^2,    tolerance = 1e-12)
  expect_equal(out$p_triple_clean,       (1   - eps)^3,    tolerance = 1e-12)
  expect_equal(out$contaminated_per_1000, eps * 1000,      tolerance = 1e-12)
})

test_that("implied_tenure_contamination: returns 14 named scalar elements", {
  out <- implied_tenure_contamination(.p_full)

  expected_names <- c(
    "entry_rate", "exit_rate", "employment_rate",
    "pi", "eps",
    "mean_spell_g_years", "mean_spell_g_months",
    "mean_spell_d_years", "mean_spell_d_months",
    "median_spell_g_months",
    "p_clock_consistent", "p_pair_clean", "p_triple_clean",
    "contaminated_per_1000"
  )
  expect_setequal(names(out), expected_names)
  for (nm in expected_names) {
    expect_true(is.numeric(out[[nm]]), info = nm)
    expect_equal(length(out[[nm]]), 1L, info = nm)
  }
})

# ---------------------------------------------------------------------------
# 2. implied_tenure_contamination: edge cases
# ---------------------------------------------------------------------------

test_that("implied_tenure_contamination: eps=0 → all clean, 0 contaminated", {
  p0 <- .p_full
  p0$eps <- 0
  out <- implied_tenure_contamination(p0)

  expect_equal(out$p_clock_consistent,   1,    tolerance = 1e-12)
  expect_equal(out$p_pair_clean,         1,    tolerance = 1e-12)
  expect_equal(out$p_triple_clean,       1,    tolerance = 1e-12)
  expect_equal(out$contaminated_per_1000, 0,   tolerance = 1e-12)
})

test_that("implied_tenure_contamination: eps=1 → all contaminated", {
  p1 <- .p_full
  p1$eps <- 1
  out <- implied_tenure_contamination(p1)

  expect_equal(out$p_clock_consistent,   0,    tolerance = 1e-12)
  expect_equal(out$p_pair_clean,         0,    tolerance = 1e-12)
  expect_equal(out$p_triple_clean,       0,    tolerance = 1e-12)
  expect_equal(out$contaminated_per_1000, 1000, tolerance = 1e-12)
})

test_that("implied_tenure_contamination: lambda_g=0 → Inf spell durations", {
  p_zero_lg <- .p_full
  p_zero_lg$lambda_g <- 0
  out <- implied_tenure_contamination(p_zero_lg)

  expect_true(is.infinite(out$mean_spell_g_years))
  expect_true(is.infinite(out$mean_spell_g_months))
  expect_true(is.infinite(out$median_spell_g_months))
  # lambda_d unchanged
  expect_equal(out$mean_spell_d_years, 1 / .p_full$lambda_d, tolerance = 1e-12)
})

test_that("implied_tenure_contamination: stationarity formula verified", {
  # employment_rate = theta0 / (theta0 + 1 - theta1)
  p <- list(alpha = 0.99, theta0 = 0.20, theta1 = 0.80,
            pi = 0.05, eps = 0.10, lambda_g = 1.0, lambda_d = 0.5)
  out <- implied_tenure_contamination(p)
  # 0.20 / (0.20 + 0.20) = 0.50
  expect_equal(out$employment_rate, 0.50, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. implied_tenure_contamination: input validation
# ---------------------------------------------------------------------------

test_that("implied_tenure_contamination: missing param throws error", {
  p_bad <- .p_full
  p_bad$eps <- NULL
  expect_error(implied_tenure_contamination(p_bad), "missing required parameters")
})

test_that("implied_tenure_contamination: pi out of [0,1] throws error", {
  p_bad <- .p_full
  p_bad$pi <- 1.5
  expect_error(implied_tenure_contamination(p_bad), "out of \\[0,1\\]")
})

test_that("implied_tenure_contamination: eps out of [0,1] throws error", {
  p_bad <- .p_full
  p_bad$eps <- -0.01
  expect_error(implied_tenure_contamination(p_bad), "out of \\[0,1\\]")
})

test_that("implied_tenure_contamination: lambda_g < 0 throws error", {
  p_bad <- .p_full
  p_bad$lambda_g <- -1
  expect_error(implied_tenure_contamination(p_bad), "lambda_g must be >= 0")
})

test_that("implied_tenure_contamination: NA param throws error", {
  p_bad <- .p_full
  p_bad$theta0 <- NA_real_
  expect_error(implied_tenure_contamination(p_bad))
})

# ---------------------------------------------------------------------------
# 4. bootstrap_one_eps: structural output tests (small synthetic data)
# ---------------------------------------------------------------------------

.make_eps_boot_data <- function(n = 100L, seed = 123L) {
  set.seed(seed)
  s  <- sample(0:1, n, replace = TRUE, prob = c(0.4, 0.6))
  g1 <- pmax(0.26, rexp(n, 2))
  g2 <- ifelse(runif(n) < 0.63, g1 + 0.25, pmax(0.26, rexp(n, 2)))
  g3 <- ifelse(runif(n) < 0.63, g2 + 0.25, pmax(0.26, rexp(n, 2)))
  as.data.frame(list(
    y1           = s, y2 = s, y3 = s,
    tenure1      = g1, tenure2 = g2, tenure3 = g3,
    timegap_cat1 = sample(1:7, n, replace = TRUE),
    timegap_cat2 = sample(1:7, n, replace = TRUE),
    timegap_cat3 = sample(1:7, n, replace = TRUE),
    weight       = rep(1, n)
  ))
}

test_that("bootstrap_one_eps: returns correct output structure", {
  df   <- .make_eps_boot_data()
  params_start <- list(
    alpha = 0.6, theta0 = 0.10, theta1 = 0.90,
    pi = 0.05, eps = 0.20, lambda_g = 2.0, lambda_d = 1.5
  )

  res <- bootstrap_one_eps(
    df            = df,
    seed          = 42L,
    stationary    = FALSE,
    linked        = FALSE,
    params_start  = params_start,
    point_loglik  = -1e6
  )

  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "low_loglik", "error"))
  expect_true(is.numeric(res$loglik))
})

test_that("bootstrap_one_eps: implied has correct names when flag='ok'", {
  df <- .make_eps_boot_data(n = 200L, seed = 7L)
  params_start <- list(
    alpha = 0.6, theta0 = 0.10, theta1 = 0.90,
    pi = 0.05, eps = 0.20, lambda_g = 2.0, lambda_d = 1.5
  )

  res <- bootstrap_one_eps(
    df           = df,
    seed         = 99L,
    stationary   = FALSE,
    linked       = FALSE,
    params_start = params_start,
    point_loglik = NA_real_   # skip LL threshold check
  )

  if (res$flag == "ok") {
    expect_false(is.null(res$implied))
    expect_true("p_clock_consistent" %in% names(res$implied))
    expect_true("mean_spell_g_months" %in% names(res$implied))
    expect_true("contaminated_per_1000" %in% names(res$implied))
  } else {
    skip("bootstrap rep did not converge on synthetic data — structural test skipped")
  }
})

test_that(".flag_fit_eps: flags low loglik correctly", {
  fit_ok     <- list(converged = TRUE, loglik = -1000)
  fit_low    <- list(converged = TRUE, loglik = -2000)
  fit_nc     <- list(converged = FALSE, loglik = -1000)
  fit_na     <- list(converged = TRUE, loglik = NA_real_)

  expect_equal(.flag_fit_eps(fit_ok,  point_loglik = -1000), "ok")
  expect_equal(.flag_fit_eps(fit_low, point_loglik = -1000), "low_loglik")
  expect_equal(.flag_fit_eps(fit_nc,  point_loglik = -1000), "no_converge")
  expect_equal(.flag_fit_eps(fit_na,  point_loglik = -1000), "error")
})

test_that(".flag_fit_eps: point_loglik=NA skips low_loglik check → ok", {
  # When point_loglik is NA, the LL-threshold check is skipped entirely.
  # A rep far below any reasonable loglik still gets "ok" (no reference to compare).
  fit_low <- list(converged = TRUE, loglik = -2000)
  expect_equal(.flag_fit_eps(fit_low, point_loglik = NA_real_), "ok")
})

# ---------------------------------------------------------------------------
# 5. implied_tenure_contamination: additional edge cases
# ---------------------------------------------------------------------------

test_that("implied_tenure_contamination: lambda_d=0 → Inf non-employment spell", {
  p_zero_ld <- .p_full
  p_zero_ld$lambda_d <- 0
  out <- implied_tenure_contamination(p_zero_ld)

  expect_true(is.infinite(out$mean_spell_d_years))
  expect_true(is.infinite(out$mean_spell_d_months))
  # lambda_g unchanged
  expect_equal(out$mean_spell_g_years, 1 / .p_full$lambda_g, tolerance = 1e-12)
})

test_that("implied_tenure_contamination: employment_rate=NA when denom=0", {
  p_zero_denom <- .p_full
  p_zero_denom$theta0 <- 0
  p_zero_denom$theta1 <- 1
  out <- implied_tenure_contamination(p_zero_denom)

  expect_true(is.na(out$employment_rate))
  expect_equal(out$entry_rate, 0)
  expect_equal(out$exit_rate,  0)
})

test_that("implied_tenure_contamination: Inf lambda_g raises error (not finite)", {
  p_inf <- .p_full
  p_inf$lambda_g <- Inf
  expect_error(implied_tenure_contamination(p_inf), "finite scalar numeric")
})

# ---------------------------------------------------------------------------
# 6. bootstrap_one_eps: variant coverage
# ---------------------------------------------------------------------------

test_that("bootstrap_one_eps: output structure with stationary=TRUE", {
  df <- .make_eps_boot_data(n = 200L, seed = 11L)
  params_start <- list(
    alpha = 0.6, theta0 = 0.10, theta1 = 0.90,
    pi = 0.05, eps = 0.20, lambda_g = 2.0, lambda_d = 1.5
  )
  res <- bootstrap_one_eps(
    df           = df,
    seed         = 42L,
    stationary   = TRUE,
    linked       = FALSE,
    params_start = params_start,
    point_loglik = NA_real_
  )
  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "low_loglik", "error"))
})

test_that("bootstrap_one_eps: output structure with linked=TRUE", {
  df <- .make_eps_boot_data(n = 200L, seed = 12L)
  params_start <- list(
    alpha = 0.6, theta0 = 0.10, theta1 = 0.90,
    pi = 0.05, eps = 0.20, lambda_g = 2.0, lambda_d = 1.5
  )
  res <- bootstrap_one_eps(
    df           = df,
    seed         = 43L,
    stationary   = FALSE,
    linked       = TRUE,
    params_start = params_start,
    point_loglik = NA_real_
  )
  expect_true(is.list(res))
  expect_true(all(c("params", "implied", "loglik", "converged", "flag") %in% names(res)))
  expect_true(res$flag %in% c("ok", "no_converge", "low_loglik", "error"))
})
