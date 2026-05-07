# ==============================================================================
# EM-baseline: Tests for bootstrap_utils.R
# Created: 2026-05-06
# ==============================================================================

# ---------------------------------------------------------------------------
# bootstrap_resample()
# ---------------------------------------------------------------------------

test_that("bootstrap_resample preserves dimensions", {
  df  <- data.frame(y1 = c(1,0,1), y2 = c(0,1,1), y3 = c(1,1,0), weight = c(1,1,1))
  res <- bootstrap_resample(df, seed = 1L)
  expect_equal(nrow(res$df), nrow(df))
  expect_equal(ncol(res$df), ncol(df))
  expect_null(res$X)
  expect_null(res$incons_mat)
})

test_that("bootstrap_resample preserves column names", {
  df  <- data.frame(y1 = c(1,0,1), y2 = c(0,1,1), y3 = c(1,1,0), weight = c(1,1,1))
  res <- bootstrap_resample(df, seed = 1L)
  expect_equal(names(res$df), names(df))
})

test_that("bootstrap_resample is reproducible given same seed", {
  df   <- data.frame(y1 = 0:9, y2 = 0:9, y3 = 0:9, weight = rep(1, 10))
  res1 <- bootstrap_resample(df, seed = 42L)
  res2 <- bootstrap_resample(df, seed = 42L)
  expect_identical(res1$df, res2$df)
})

test_that("bootstrap_resample produces different samples for different seeds", {
  df   <- data.frame(y1 = 1:50, y2 = 1:50, y3 = 1:50, weight = rep(1, 50))
  res1 <- bootstrap_resample(df, seed = 1L)
  res2 <- bootstrap_resample(df, seed = 2L)
  expect_false(identical(res1$df, res2$df))
})

test_that("bootstrap_resample resamples X with same row index as df", {
  df  <- data.frame(y1 = 1:5, y2 = 1:5, y3 = 1:5, weight = rep(1, 5))
  X   <- cbind(1, c(10, 20, 30, 40, 50))
  res <- bootstrap_resample(df, seed = 7L, X = X)
  # The resampled X rows should correspond to the resampled df rows
  # i.e., X[i, 2] == 10 * df[i, "y1"] since y1 == 1..5 and X[,2] == y1*10
  expect_equal(res$X[, 2], res$df$y1 * 10)
})

test_that("bootstrap_resample resamples incons_mat with same row index", {
  df       <- data.frame(y1 = 1:5, y2 = 1:5, y3 = 1:5, weight = rep(1, 5))
  imat     <- matrix(seq_len(30), nrow = 5, ncol = 6)
  res      <- bootstrap_resample(df, seed = 3L, incons_mat = imat)
  # Row 1 of imat should match corresponding df row
  df_y1    <- res$df$y1         # resampled y1 values (1-based row ids)
  expect_equal(res$incons_mat[, 1], imat[df_y1, 1])
})

# ---------------------------------------------------------------------------
# summarise_bootstrap()
# ---------------------------------------------------------------------------

test_that("summarise_bootstrap returns correct column names", {
  pp     <- list(theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  pi_pt  <- list(entry_rate = 0.1, exit_rate = 0.1, employment_rate = 0.5, pi = 0.05,
                 pi0 = NA_real_, pi1 = NA_real_)
  # Construct fake boot_results with flag = "ok"
  make_rep <- function(t0, t1, pi) {
    p <- list(theta0 = t0, theta1 = t1, pi = pi)
    i <- list(entry_rate = t0, exit_rate = 1-t1, employment_rate = t0/(t0+1-t1),
              pi = pi, pi0 = NA_real_, pi1 = NA_real_)
    list(params = p, implied = i, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- lapply(1:10, function(i) make_rep(0.1 + i*0.001, 0.9 - i*0.001, 0.05))
  out   <- summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
  expect_true(all(c("quantity","estimate","se","ci_lower","ci_upper","n_ok","n_reps")
                  %in% names(out)))
})

test_that("summarise_bootstrap SE is 0 for constant bootstrap samples", {
  pp    <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
  pi_pt <- list(entry_rate = 0.10, exit_rate = 0.10, employment_rate = 0.5,
                pi = 0.05, pi0 = NA_real_, pi1 = NA_real_)
  make_const_rep <- function() {
    list(params = pp, implied = pi_pt, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- replicate(20, make_const_rep(), simplify = FALSE)
  out   <- summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
  # Only check rows where the point estimate is non-NA (NA quantities produce NA SE)
  scalar_rows <- out[!is.na(out$estimate) & !is.na(out$se), ]
  expect_true(all(scalar_rows$se == 0))
})

test_that("summarise_bootstrap n_ok excludes flagged reps", {
  pp    <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
  pi_pt <- list(entry_rate = 0.10, exit_rate = 0.10, employment_rate = 0.5,
                pi = 0.05, pi0 = NA_real_, pi1 = NA_real_)
  make_rep <- function(flag) {
    list(params = pp, implied = pi_pt, loglik = -100, converged = TRUE, flag = flag)
  }
  boots <- c(replicate(8, make_rep("ok"), simplify = FALSE),
             replicate(2, make_rep("no_converge"), simplify = FALSE))
  out   <- summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
  expect_equal(unique(out$n_ok),   8L)
  expect_equal(unique(out$n_reps), 10L)
})

test_that("summarise_bootstrap handles all-failed boots gracefully", {
  pp    <- list(theta0 = 0.10, theta1 = 0.90)
  pi_pt <- list(entry_rate = 0.10, exit_rate = 0.10, employment_rate = 0.5,
                pi = NA_real_, pi0 = NA_real_, pi1 = NA_real_)
  boots <- replicate(5, list(params = NULL, implied = NULL, loglik = NA,
                             converged = FALSE, flag = "error"), simplify = FALSE)
  out   <- summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
  expect_equal(nrow(out), 0L)
})

# ---------------------------------------------------------------------------
# Quality flag logic (tested via bootstrap_one_baseline on synthetic data)
# ---------------------------------------------------------------------------

test_that("bootstrap_one_baseline returns flag='ok' for well-conditioned data", {
  set.seed(99L)
  N  <- 200L
  y1 <- rbinom(N, 1, 0.5); y2 <- rbinom(N, 1, 0.5); y3 <- rbinom(N, 1, 0.5)
  df <- data.frame(y1 = y1, y2 = y2, y3 = y3, weight = rep(1, N))
  p0 <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05, alpha = 0.5)
  res <- bootstrap_one_baseline(df, seed = 1L, model_type = "symmetric",
                                stationary = TRUE, params_start = p0,
                                point_loglik = -1e6)  # very low threshold -> never low_loglik
  expect_equal(res$flag, "ok")
})

test_that("summarise_bootstrap output covers all unique quantities exactly once", {
  # Symmetric baseline: 'pi' appears in both params and implied — dedup must keep it once.
  # The output should contain every unique scalar quantity name from params + implied,
  # each appearing exactly once.
  pp    <- list(theta0 = 0.10, theta1 = 0.90, pi = 0.05)
  pi_pt <- list(entry_rate = 0.10, exit_rate = 0.10, employment_rate = 0.5,
                pi = 0.05, pi0 = NA_real_, pi1 = NA_real_)
  make_rep <- function() {
    list(params = pp, implied = pi_pt, loglik = -100, converged = TRUE, flag = "ok")
  }
  boots <- replicate(5, make_rep(), simplify = FALSE)
  out   <- summarise_bootstrap(boots, pp, pi_pt, point_loglik = -100)
  # All unique scalar names from params + implied should appear exactly once
  all_scalars <- unique(c(
    names(Filter(function(x) is.numeric(x) && length(x) == 1L, pp)),
    names(Filter(function(x) is.numeric(x) && length(x) == 1L, pi_pt))
  ))
  expect_equal(sort(out$quantity), sort(all_scalars))
})

test_that("bootstrap_resample errors when df is missing required columns", {
  bad_df <- data.frame(a = 1:5, b = 1:5)
  expect_error(bootstrap_resample(bad_df, seed = 1L), "missing required columns")
})
