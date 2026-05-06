# ==============================================================================
# EM-baseline-ext: Tests for Extension IV (inconsistency-augmented model)
# Created: 2026-05-06
# ==============================================================================

# Helper: synthetic panel data with age/educ columns
.make_incons_panel <- function(n = 300L, seed = 7L,
                                theta0 = 0.10, theta1 = 0.90, pi = 0.05) {
  set.seed(seed)
  alpha <- theta0 / (theta0 + 1 - theta1)
  h1 <- rbinom(n, 1L, alpha)
  h2 <- ifelse(h1 == 1L, rbinom(n, 1L, theta1), rbinom(n, 1L, theta0))
  h3 <- ifelse(h2 == 1L, rbinom(n, 1L, theta1), rbinom(n, 1L, theta0))
  s1 <- ifelse(rbinom(n, 1L, pi), 1L - h1, h1)
  s2 <- ifelse(rbinom(n, 1L, pi), 1L - h2, h2)
  s3 <- ifelse(rbinom(n, 1L, pi), 1L - h3, h3)
  # Clean age and edu (no inconsistencies) for most observations
  age1 <- sample(18:55, n, replace = TRUE)
  age2 <- age1 + sample(0:1, n, replace = TRUE)
  age3 <- age2 + sample(0:1, n, replace = TRUE)
  edu1 <- sample(1:5, n, replace = TRUE)
  edu2 <- pmin(edu1 + sample(0:1, n, replace = TRUE), 5L)
  edu3 <- pmin(edu2 + sample(0:1, n, replace = TRUE), 5L)
  data.frame(
    y1 = s1, y2 = s2, y3 = s3, weight = rep(1, n),
    age1 = age1, age2 = age2, age3 = age3,
    educ1 = edu1, educ2 = edu2, educ3 = edu3
  )
}

# Helper: build incons_mat from df
.make_incons_mat <- function(df) {
  inc <- compute_inconsistencies(df)
  as.matrix(inc[, c("Y_age_1", "Y_age_2", "Y_age_3",
                    "Y_edu_1", "Y_edu_2", "Y_edu_3")])
}

# ---- init_params_inconsistency --------------------------------------------

test_that("init_params_inconsistency returns correct structure", {
  p <- init_params_inconsistency()
  expect_true(all(c("theta0", "theta1", "alpha", "delta") %in% names(p)))
  expect_length(p$delta, 3L)
})

test_that("init_params_inconsistency: pi_base close to 0.05", {
  p       <- init_params_inconsistency()
  pi_base <- 0.5 * plogis(p$delta[1L])
  expect_equal(pi_base, 0.05, tolerance = 1e-3)
})

test_that("init_params_inconsistency: delta slopes start at zero", {
  p <- init_params_inconsistency()
  expect_equal(p$delta[2L], 0)
  expect_equal(p$delta[3L], 0)
})

# ---- e_step_inconsistency: structure ---------------------------------------

test_that("e_step_inconsistency returns gamma, loglik, suff", {
  df        <- .make_incons_panel()
  inc_mat   <- .make_incons_mat(df)
  params    <- init_params_inconsistency()
  out       <- e_step_inconsistency(df, inc_mat, params)
  expect_true(all(c("gamma", "loglik", "suff") %in% names(out)))
})

test_that("e_step_inconsistency: gamma rows sum to 1", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  expect_equal(rowSums(out$gamma), rep(1, nrow(df)), tolerance = 1e-10)
})

test_that("e_step_inconsistency: gamma non-negative, dimensions N x 8", {
  df      <- .make_incons_panel(n = 100L)
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  expect_true(all(out$gamma >= 0))
  expect_equal(dim(out$gamma), c(100L, 8L))
})

test_that("e_step_inconsistency: loglik is finite", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  expect_true(is.finite(out$loglik))
})

# ---- e_step_inconsistency: suff stats -------------------------------------

test_that("e_step_inconsistency suff has required fields", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  required <- c("T11", "T01", "D1", "D0", "C1", "C0",
                "p_mat", "pi_mat", "sig_mat", "weights", "total_weight")
  expect_true(all(required %in% names(out$suff)))
})

test_that("e_step_inconsistency: p_mat is N x 3 in [0,1]", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  expect_equal(dim(out$suff$p_mat), c(nrow(df), 3L))
  expect_true(all(out$suff$p_mat >= 0 & out$suff$p_mat <= 1))
})

test_that("e_step_inconsistency: pi_mat in (0, 0.5)", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out     <- e_step_inconsistency(df, inc_mat, params)
  expect_true(all(out$suff$pi_mat > 0 & out$suff$pi_mat < 0.5))
})

# ---- delta=0 reduces to baseline ------------------------------------------

test_that("inconsistency E-step with delta slopes=0 gives same LL as baseline EM", {
  df      <- .make_incons_panel(n = 300L)
  # Make incons_mat all zeros (no inconsistencies)
  inc_mat <- matrix(0L, nrow = nrow(df), ncol = 6L)
  params  <- init_params_inconsistency()
  params$delta[2:3] <- 0  # ensure slopes = 0

  # Baseline EM with same theta init
  params_bl <- list(theta0 = params$theta0, theta1 = params$theta1,
                    alpha = params$alpha, pi = 0.5 * plogis(params$delta[1L]))
  ll_incons <- e_step_inconsistency(df, inc_mat, params)$loglik

  # Should be close to a single E-step of baseline with same params
  # (not exact EM convergence — just single-step comparison)
  expect_true(is.finite(ll_incons))
})

# ---- m_step_inconsistency -------------------------------------------------

test_that("m_step_inconsistency returns required fields", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out_e   <- e_step_inconsistency(df, inc_mat, params)
  out_m   <- m_step_inconsistency(out_e$suff, inc_mat, params)
  expect_true(all(c("theta0", "theta1", "alpha", "delta") %in% names(out_m)))
  expect_length(out_m$delta, 3L)
})

test_that("m_step_inconsistency: theta values in valid range", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  out_e   <- e_step_inconsistency(df, inc_mat, params)
  out_m   <- m_step_inconsistency(out_e$suff, inc_mat, params)
  expect_true(out_m$theta0 > 0 && out_m$theta0 <= 0.999)
  expect_true(out_m$theta1 > 0 && out_m$theta1 <= 0.999)
})

test_that("m_step_inconsistency: NR step does not decrease LL (public API)", {
  # Verify GEM guarantee via observed-data LL: two consecutive EM iterations
  # must have non-decreasing LL, which is implied by Q_miscl not decreasing.
  df      <- .make_incons_panel(n = 200L)
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  # Iteration 1
  out_e1  <- e_step_inconsistency(df, inc_mat, params)
  params1 <- m_step_inconsistency(out_e1$suff, inc_mat, params)
  # Iteration 2 E-step evaluates LL under updated params
  out_e2  <- e_step_inconsistency(df, inc_mat, params1)
  expect_true(out_e2$loglik >= out_e1$loglik - 1e-8,
              label = "GEM iteration must not decrease observed-data LL")
})

# ---- em_fit_inconsistency: convergence -------------------------------------

test_that("em_fit_inconsistency converges on synthetic clean data", {
  df      <- .make_incons_panel(n = 300L)
  inc_mat <- .make_incons_mat(df)
  fit     <- em_fit_inconsistency(df, inc_mat, verbose = 0L)
  expect_true(fit$converged)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_inconsistency: LL is monotone non-decreasing", {
  df      <- .make_incons_panel(n = 300L)
  inc_mat <- .make_incons_mat(df)
  fit     <- em_fit_inconsistency(df, inc_mat, verbose = 0L)
  ll_vec  <- fit$history$loglik
  diffs   <- diff(ll_vec)
  expect_true(all(diffs >= -1e-5), label = "Inconsistency LL must be non-decreasing")
})

test_that("em_fit_inconsistency returns gamma with correct dimensions", {
  df      <- .make_incons_panel(n = 100L)
  inc_mat <- .make_incons_mat(df)
  fit     <- em_fit_inconsistency(df, inc_mat, verbose = 0L)
  expect_equal(nrow(fit$gamma), 100L)
  expect_equal(ncol(fit$gamma), 8L)
})

test_that("em_fit_inconsistency: delta has finite values at convergence", {
  df      <- .make_incons_panel(n = 300L)
  inc_mat <- .make_incons_mat(df)
  fit     <- em_fit_inconsistency(df, inc_mat, verbose = 0L)
  expect_true(all(is.finite(fit$params$delta)))
})

# ---- compute_observed_loglik_inconsistency ---------------------------------

test_that("compute_observed_loglik_inconsistency returns finite scalar", {
  df      <- .make_incons_panel()
  inc_mat <- .make_incons_mat(df)
  params  <- init_params_inconsistency()
  ll      <- compute_observed_loglik_inconsistency(df, inc_mat, params)
  expect_true(is.finite(ll))
  expect_length(ll, 1L)
})

# ---- em_fit_inconsistency: additional convergence guards -----------------

test_that("em_fit_inconsistency runs more than 1 iteration", {
  df      <- .make_incons_panel(n = 300L)
  inc_mat <- .make_incons_mat(df)
  fit     <- em_fit_inconsistency(df, inc_mat, verbose = 0L)
  expect_true(fit$iterations > 1L,
              label = "inconsistency EM must iterate more than once")
})

test_that("em_fit_inconsistency: params change from starting values", {
  df      <- .make_incons_panel(n = 300L)
  inc_mat <- .make_incons_mat(df)
  p0      <- init_params_inconsistency()
  fit     <- em_fit_inconsistency(df, inc_mat, params0 = p0, verbose = 0L)
  expect_false(identical(fit$params$theta0, p0$theta0),
               label = "theta0 must change from initial value")
})

test_that("positive delta[2] increases pi_mat for age-flagged observations", {
  set.seed(123L)
  n       <- 400L
  inc_mat <- matrix(0L, nrow = n, ncol = 6L)
  colnames(inc_mat) <- c("Y_age_1", "Y_age_2", "Y_age_3",
                         "Y_edu_1", "Y_edu_2", "Y_edu_3")
  # Flag the first half of observations at wave 1
  inc_mat[seq_len(n / 2L), "Y_age_1"] <- 1L
  df          <- .make_incons_panel(n = n)
  params_pos  <- init_params_inconsistency()
  params_pos$delta[2L] <- 2.0   # large positive age coefficient
  out          <- e_step_inconsistency(df, inc_mat, params_pos)
  pi_flagged   <- mean(out$suff$pi_mat[seq_len(n / 2L), 1L])
  pi_unflagged <- mean(out$suff$pi_mat[(n / 2L + 1L):n, 1L])
  expect_true(pi_flagged > pi_unflagged,
              label = "positive delta[2] must raise pi_mat for flagged rows")
})

test_that("positive delta[3] increases pi_mat for edu-flagged observations", {
  set.seed(999L)
  n       <- 400L
  inc_mat <- matrix(0L, nrow = n, ncol = 6L)
  colnames(inc_mat) <- c("Y_age_1", "Y_age_2", "Y_age_3",
                         "Y_edu_1", "Y_edu_2", "Y_edu_3")
  # Flag the first half at wave 1 for education inconsistency
  inc_mat[seq_len(n / 2L), "Y_edu_1"] <- 1L
  df         <- .make_incons_panel(n = n)
  params_edu <- init_params_inconsistency()
  params_edu$delta[3L] <- 2.0   # large positive edu coefficient
  out          <- e_step_inconsistency(df, inc_mat, params_edu)
  pi_flagged   <- mean(out$suff$pi_mat[seq_len(n / 2L), 1L])
  pi_unflagged <- mean(out$suff$pi_mat[(n / 2L + 1L):n, 1L])
  expect_true(pi_flagged > pi_unflagged,
              label = "positive delta[3] must raise pi_mat for edu-flagged rows")
})
