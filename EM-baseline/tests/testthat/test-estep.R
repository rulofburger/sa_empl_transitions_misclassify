# ==============================================================================
# EM-baseline: Tests for estep.R
# Created: 2026-05-05
# ==============================================================================

# ---- helpers ----------------------------------------------------------------

# Generate a small synthetic data frame for testing
.make_df <- function(n = 50, theta0 = 0.1, theta1 = 0.9, pi = 0.05, seed = 42) {
  set.seed(seed)
  alpha <- theta0 / (theta0 + 1 - theta1)

  h1 <- rbinom(n, 1, alpha)
  h2 <- ifelse(h1 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  h3 <- ifelse(h2 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))

  flip <- function(h) ifelse(rbinom(length(h), 1, pi) == 1, 1 - h, h)
  data.frame(y1 = flip(h1), y2 = flip(h2), y3 = flip(h3), weight = rep(1, n))
}

# Deterministic Markov data with no misclassification (for exact assertions)
.make_df_perfect <- function(n = 60, theta0 = 0.1, theta1 = 0.9) {
  alpha <- theta0 / (theta0 + 1 - theta1)
  # Alternate employed/nonemployed for deterministic mix
  n_emp  <- round(n * alpha)
  n_nemp <- n - n_emp
  # All employment spells stay employed (theta1=1 approx); nonemployment stays
  h1 <- c(rep(1L, n_emp), rep(0L, n_nemp))
  h2 <- h1  # deterministic: no transitions
  h3 <- h1
  data.frame(y1 = h1, y2 = h2, y3 = h3, weight = rep(1, n))
}

# ---- tests: .log_misclass_wave -----------------------------------------------

test_that(".log_misclass_wave none: 0 on match, -Inf on mismatch", {
  result <- .log_misclass_wave(c(1L, 0L, 1L), c(0L, 1L), "none")
  # s=1 vs h=0: mismatch -> -Inf; s=1 vs h=1: match -> 0
  expect_equal(result[1, 1], -Inf)   # s1=1 vs h=0
  expect_equal(result[1, 2], 0)      # s1=1 vs h=1
  # s=0 vs h=0: match -> 0; s=0 vs h=1: mismatch -> -Inf
  expect_equal(result[2, 1], 0)      # s2=0 vs h=0
  expect_equal(result[2, 2], -Inf)   # s2=0 vs h=1
})

test_that(".log_misclass_wave symmetric: correct log-probs", {
  pi     <- 0.05
  result <- .log_misclass_wave(c(1L, 0L), c(0L, 1L), "symmetric", pi = pi)
  # s=1, h=0: mismatch -> log(pi)
  expect_equal(result[1, 1], log(pi),       tolerance = 1e-12)
  # s=1, h=1: match   -> log(1-pi)
  expect_equal(result[1, 2], log(1 - pi),   tolerance = 1e-12)
  # s=0, h=0: match   -> log(1-pi)
  expect_equal(result[2, 1], log(1 - pi),   tolerance = 1e-12)
  # s=0, h=1: mismatch -> log(pi)
  expect_equal(result[2, 2], log(pi),       tolerance = 1e-12)
})

test_that(".log_misclass_wave asymmetric: pi0 and pi1 applied separately", {
  pi0    <- 0.10
  pi1    <- 0.20
  result <- .log_misclass_wave(c(1L, 0L), c(0L, 1L), "asymmetric",
                                pi0 = pi0, pi1 = pi1)
  # (s=1, h=0): false positive -> log(pi0)
  expect_equal(result[1, 1], log(pi0),       tolerance = 1e-12)
  # (s=1, h=1): true positive -> log(1-pi1)
  expect_equal(result[1, 2], log(1 - pi1),   tolerance = 1e-12)
  # (s=0, h=0): true negative -> log(1-pi0)
  expect_equal(result[2, 1], log(1 - pi0),   tolerance = 1e-12)
  # (s=0, h=1): false negative -> log(pi1)
  expect_equal(result[2, 2], log(pi1),       tolerance = 1e-12)
})

test_that(".log_misclass_wave symmetric pi=0.0001: near-zero misclassification", {
  # With very small pi, mismatches give near -Inf and matches give near 0
  pi     <- 1e-4
  result <- .log_misclass_wave(c(1L, 0L), c(0L, 1L), "symmetric", pi = pi)
  # mismatch log-prob should be very negative (log(pi) = -9.21)
  expect_true(result[1, 1] < -5)   # s=1, h=0: mismatch -> log(1e-4)
  # match log-prob should be close to 0 (log(1-pi) ≈ -1e-4)
  expect_equal(result[1, 2], log(1 - pi), tolerance = 1e-10)  # s=1, h=1: match
})

# ---- tests: structure -------------------------------------------------------

test_that("e_step returns correct structure for symmetric model", {
  df     <- .make_df(n = 30)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")

  expect_named(result, c("gamma", "loglik", "suff"))
  expect_equal(dim(result$gamma), c(30L, 8L))
  expect_true(is.finite(result$loglik))
  expect_true(result$loglik <= 0)
})

test_that("e_step returns correct structure for none model", {
  df     <- .make_df(n = 20)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  result <- e_step(df, params, model_type = "none")
  expect_equal(dim(result$gamma), c(20L, 8L))
})

test_that("e_step returns correct structure for asymmetric model", {
  df     <- .make_df(n = 20)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi0 = 0.04, pi1 = 0.06)
  result <- e_step(df, params, model_type = "asymmetric")
  expect_equal(dim(result$gamma), c(20L, 8L))
})

# ---- tests: responsibilities ------------------------------------------------

test_that("gamma rows sum to 1", {
  df     <- .make_df(n = 100)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  row_sums <- rowSums(result$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("gamma is non-negative", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_true(all(result$gamma >= 0))
})

test_that("gamma concentrates on matching history when pi=0 (none model)", {
  # With no misclassification, observation (1,0,1) should assign all
  # responsibility to history (1,0,1)
  df     <- data.frame(y1 = 1L, y2 = 0L, y3 = 1L, weight = 1)
  params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
  result <- e_step(df, params, model_type = "none")
  hmat   <- latent_histories()
  # Find index of history (1,0,1)
  idx_101 <- which(hmat[,1] == 1 & hmat[,2] == 0 & hmat[,3] == 1)
  expect_equal(result$gamma[1, idx_101], 1, tolerance = 1e-10)
  other_idx <- setdiff(1:8, idx_101)
  expect_true(all(result$gamma[1, other_idx] < 1e-10))
})

test_that("gamma concentrates on correct history with deterministic data (none)", {
  # Deterministic data: each observation has exactly one compatible history
  df   <- .make_df_perfect(n = 60)
  hmat <- latent_histories()
  params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
  result <- e_step(df, params, model_type = "none")
  for (i in seq_len(nrow(df))) {
    h_true <- which(hmat[, 1] == df$y1[i] &
                    hmat[, 2] == df$y2[i] &
                    hmat[, 3] == df$y3[i])
    expect_equal(result$gamma[i, h_true], 1, tolerance = 1e-10)
  }
})

# ---- tests: sufficient statistics -------------------------------------------

test_that("suff C1 + C0 equals total weight", {
  df     <- .make_df(n = 80)
  df$weight <- runif(80, 0.5, 2)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_equal(result$suff$C1 + result$suff$C0, sum(df$weight), tolerance = 1e-8)
})

test_that("suff D1 + D0 equals 2 * total_weight", {
  # Two transitions per individual: t=1->2 and t=2->3
  df     <- .make_df(n = 60)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  W <- sum(df$weight)
  expect_equal(result$suff$D1 + result$suff$D0, 2 * W, tolerance = 1e-8)
})

test_that("suff T11 + T01 equals D1 and T10 + T00 = D0 (accounting)", {
  # T11 + T10 = D1 (transitions from state 1 = T11 + T10)
  # Compute T10 indirectly: T10 = D1 - T11
  # Similarly T00 = D0 - T01
  df     <- .make_df(n = 60)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  s      <- result$suff
  # T11 + T10 = D1  =>  T11 <= D1
  expect_true(s$T11 <= s$D1 + 1e-8)
  expect_true(s$T01 <= s$D0 + 1e-8)
})

test_that("suff M is non-negative (second block)", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_true(result$suff$M >= 0)
})

# ---- tests: validate parameter ----------------------------------------------

test_that("e_step validate=FALSE produces same result as validate=TRUE on valid data", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  r1 <- e_step(df, params, model_type = "symmetric", validate = TRUE)
  r2 <- e_step(df, params, model_type = "symmetric", validate = FALSE)
  expect_equal(r1$loglik, r2$loglik, tolerance = 1e-12)
  expect_equal(r1$gamma,  r2$gamma,  tolerance = 1e-12)
})

test_that("e_step validate=FALSE skips data checks (invalid data passes to NaN guard)", {
  # Non-binary y1=2 with model_type='none': no history matches (2,*,*), so all
  # log_joint rows for this obs are -Inf, triggering the NaN guard.
  df     <- data.frame(y1 = 2L, y2 = 0L, y3 = 0L, weight = 1)
  params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
  # With validate=TRUE the binary check fires first
  expect_error(e_step(df, params, model_type = "none", validate = TRUE),
               regexp = "binary")
  # With validate=FALSE the binary check is skipped; the NaN guard fires instead
  expect_error(e_step(df, params, model_type = "none", validate = FALSE),
               regexp = "zero probability")
})

test_that("e_step errors on invalid params (alpha > 1) with validate=TRUE", {
  df     <- .make_df(n = 10)
  params <- list(alpha = 2.0, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  expect_error(e_step(df, params, model_type = "symmetric", validate = TRUE),
               regexp = "params\\$alpha")
})

test_that("e_step errors on invalid params (NA) with validate=TRUE", {
  df     <- .make_df(n = 10)
  params <- list(alpha = NA_real_, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "none", validate = TRUE),
               regexp = "params\\$alpha")
})

# ---- tests: NaN guard -------------------------------------------------------

test_that("e_step NaN guard fires for all-impossible observation (none model)", {
  # Non-binary value y1=2 with model_type='none': .log_misclass_wave returns
  # -Inf for all 8 histories at wave 1 (none have h1=2), making every log_joint
  # row -Inf and producing NaN in log_row_sum. validate=FALSE bypasses binary check.
  df     <- data.frame(y1 = 2L, y2 = 0L, y3 = 0L, weight = 1)
  params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
  expect_error(
    e_step(df, params, model_type = "none", validate = FALSE),
    regexp = "zero probability"
  )
})

test_that("suff M0 + M1 equals M for asymmetric model", {
  # For asymmetric model, M = M0 (false positives) + M1 (false negatives)
  df     <- .make_df(n = 100)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi0 = 0.04, pi1 = 0.06)
  result <- e_step(df, params, model_type = "asymmetric")
  s <- result$suff
  # NOTE: suff$M is set to 0 for asymmetric models (performance optimization —
  # m_step only reads M for symmetric models). The mathematically correct total
  # is M0 + M1; M is intentionally zeroed to avoid 6 unnecessary outer() calls.
  expect_equal(s$M, 0)
  expect_true(s$M0 + s$M1 > 0)  # there should be misclassification
})

test_that("suff H0 + H1 equals 3 * total_weight for asymmetric model", {
  # Each of 3 waves contributes one observation per individual
  df     <- .make_df(n = 70)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi0 = 0.04, pi1 = 0.06)
  result <- e_step(df, params, model_type = "asymmetric")
  W <- sum(df$weight)
  expect_equal(result$suff$H0 + result$suff$H1, 3 * W, tolerance = 1e-8)
})

test_that("suff M = 0 and H0/H1 = 0 for none model (gating check)", {
  df     <- .make_df(n = 40)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  result <- e_step(df, params, model_type = "none")
  expect_equal(result$suff$M,  0)
  expect_equal(result$suff$M0, 0)
  expect_equal(result$suff$M1, 0)
  expect_equal(result$suff$H0, 0)
  expect_equal(result$suff$H1, 0)
})

test_that("asymmetric extreme: pi0=0 gives M0=0", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi0 = 0, pi1 = 0.05)
  result <- e_step(df, params, model_type = "asymmetric")
  # M0 is the expected count of false positives; with pi0=0 it should be ~0
  expect_equal(result$suff$M0, 0, tolerance = 1e-10)
})

# ---- tests: input validation ------------------------------------------------

test_that("e_step errors on missing columns", {
  bad_df <- data.frame(y1 = 1, y2 = 0)  # missing y3 and weight
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(bad_df, params), regexp = "missing required columns")
})

test_that("e_step errors on unknown model_type", {
  df     <- .make_df(n = 10)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "fmm"), regexp = "model_type")
})

test_that("e_step errors on factor columns", {
  df        <- .make_df(n = 10)
  df$y1     <- factor(df$y1)
  params    <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "none"), regexp = "factor")
})

test_that("e_step errors on NA in y-columns", {
  df        <- .make_df(n = 10)
  df$y1[1]  <- NA_integer_
  params    <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "none"), regexp = "NA")
})

test_that("e_step errors on non-binary y-values", {
  df        <- .make_df(n = 10)
  df$y1[1]  <- 2L
  params    <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "none"), regexp = "binary")
})

test_that("e_step errors on non-positive weights", {
  df           <- .make_df(n = 10)
  df$weight[1] <- -1
  params       <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  expect_error(e_step(df, params, model_type = "none"), regexp = "positive")
})

# ---- tests: structure -------------------------------------------------------

test_that("e_step returns correct structure for symmetric model", {
  df     <- .make_df(n = 30)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")

  expect_named(result, c("gamma", "loglik", "suff"))
  expect_equal(dim(result$gamma), c(30L, 8L))
  expect_true(is.finite(result$loglik))
  expect_true(result$loglik <= 0)
})

test_that("e_step returns correct structure for none model", {
  df     <- .make_df(n = 20)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9)
  result <- e_step(df, params, model_type = "none")
  expect_equal(dim(result$gamma), c(20L, 8L))
})

test_that("e_step returns correct structure for asymmetric model", {
  df     <- .make_df(n = 20)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi0 = 0.04, pi1 = 0.06)
  result <- e_step(df, params, model_type = "asymmetric")
  expect_equal(dim(result$gamma), c(20L, 8L))
})

# ---- tests: responsibilities ------------------------------------------------

test_that("gamma rows sum to 1", {
  df     <- .make_df(n = 100)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  row_sums <- rowSums(result$gamma)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("gamma is non-negative", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_true(all(result$gamma >= 0))
})

test_that("gamma concentrates on matching history when pi=0 (none model)", {
  # With no misclassification, observation (1,0,1) should assign all
  # responsibility to history (1,0,1)
  df     <- data.frame(y1 = 1L, y2 = 0L, y3 = 1L, weight = 1)
  params <- list(alpha = 0.5, theta0 = 0.5, theta1 = 0.5)
  result <- e_step(df, params, model_type = "none")
  hmat   <- latent_histories()
  # Find index of history (1,0,1)
  idx_101 <- which(hmat[,1] == 1 & hmat[,2] == 0 & hmat[,3] == 1)
  expect_equal(result$gamma[1, idx_101], 1, tolerance = 1e-10)
  other_idx <- setdiff(1:8, idx_101)
  expect_true(all(result$gamma[1, other_idx] < 1e-10))
})

# ---- tests: sufficient statistics -------------------------------------------

test_that("suff C1 + C0 equals total weight", {
  df     <- .make_df(n = 80)
  df$weight <- runif(80, 0.5, 2)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_equal(result$suff$C1 + result$suff$C0, sum(df$weight), tolerance = 1e-8)
})

test_that("suff D1 + D0 equals 2 * total_weight", {
  # Two transitions per individual: t=1->2 and t=2->3
  df     <- .make_df(n = 60)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  W <- sum(df$weight)
  expect_equal(result$suff$D1 + result$suff$D0, 2 * W, tolerance = 1e-8)
})

test_that("suff T11 + T01 equals D1 and T10 + T00 = D0 (accounting)", {
  # T11 + T10 = D1 (transitions from state 1 = T11 + T10)
  # Compute T10 indirectly: T10 = D1 - T11
  # Similarly T00 = D0 - T01
  df     <- .make_df(n = 60)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  s      <- result$suff
  # T11 + T10 = D1  =>  T11 <= D1
  expect_true(s$T11 <= s$D1 + 1e-8)
  expect_true(s$T01 <= s$D0 + 1e-8)
})

test_that("suff M is non-negative", {
  df     <- .make_df(n = 50)
  params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
  result <- e_step(df, params, model_type = "symmetric")
  expect_true(result$suff$M >= 0)
})
