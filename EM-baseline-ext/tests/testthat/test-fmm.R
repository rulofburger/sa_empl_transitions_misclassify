# ==============================================================================
# EM-baseline-ext: Tests for Extension III (2-type FMM)
# Created: 2026-05-06
# ==============================================================================

# Helper: simulate 2-type mixture panel data
.make_fmm_panel <- function(n = 300L, seed = 42L,
                             theta0_A = 0.15, theta1_A = 0.92,
                             theta0_B = 0.05, theta1_B = 0.75,
                             phi = 0.4, pi = 0.05) {
  set.seed(seed)
  type    <- ifelse(runif(n) < phi, "A", "B")
  th0     <- ifelse(type == "A", theta0_A, theta0_B)
  th1     <- ifelse(type == "A", theta1_A, theta1_B)
  alpha_A <- theta0_A / (theta0_A + 1 - theta1_A)
  alpha_B <- theta0_B / (theta0_B + 1 - theta1_B)
  alpha   <- ifelse(type == "A", alpha_A, alpha_B)
  h1      <- rbinom(n, 1L, alpha)
  h2      <- ifelse(h1 == 1L, rbinom(n, 1L, th1), rbinom(n, 1L, th0))
  h3      <- ifelse(h2 == 1L, rbinom(n, 1L, th1), rbinom(n, 1L, th0))
  s1      <- ifelse(rbinom(n, 1L, pi), 1L - h1, h1)
  s2      <- ifelse(rbinom(n, 1L, pi), 1L - h2, h2)
  s3      <- ifelse(rbinom(n, 1L, pi), 1L - h3, h3)
  data.frame(y1 = s1, y2 = s2, y3 = s3, weight = rep(1, n))
}

# ---- init_params_fmm -------------------------------------------------------

test_that("init_params_fmm returns correct structure for symmetric", {
  p <- init_params_fmm("symmetric")
  required <- c("theta0_A", "theta1_A", "alpha_A",
                "theta0_B", "theta1_B", "alpha_B", "phi", "pi")
  expect_true(all(required %in% names(p)))
  expect_true(p$phi > 0 && p$phi < 1)
  expect_true(p$pi  > 0 && p$pi  < 0.5)
})

test_that("init_params_fmm: type A more stable than type B by default", {
  p <- init_params_fmm()
  expect_true(p$theta1_A >= p$theta1_B)
})

test_that("init_params_fmm: no pi for model_type='none'", {
  p <- init_params_fmm("none")
  expect_null(p$pi)
})

# ---- e_step_fmm: output structure -----------------------------------------

test_that("e_step_fmm returns gamma_A, gamma_B, loglik, suff", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  expect_true(all(c("gamma_A", "gamma_B", "loglik", "suff") %in% names(out)))
})

test_that("e_step_fmm: gamma_A + gamma_B rows sum to 1", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  total  <- out$gamma_A + out$gamma_B
  expect_equal(rowSums(total), rep(1, nrow(df)), tolerance = 1e-10)
})

test_that("e_step_fmm: gamma_A and gamma_B are non-negative", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  expect_true(all(out$gamma_A >= 0))
  expect_true(all(out$gamma_B >= 0))
})

test_that("e_step_fmm: loglik is finite", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  expect_true(is.finite(out$loglik))
})

test_that("e_step_fmm: gamma matrices are N x 8", {
  df     <- .make_fmm_panel(n = 100L)
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  expect_equal(nrow(out$gamma_A), 100L)
  expect_equal(ncol(out$gamma_A), 8L)
  expect_equal(dim(out$gamma_B), dim(out$gamma_A))
})

# ---- e_step_fmm: suff statistics ------------------------------------------

test_that("e_step_fmm suff has required fields", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  required_suff <- c("T11_A", "T01_A", "D1_A", "D0_A", "C1_A", "C0_A",
                     "T11_B", "T01_B", "D1_B", "D0_B", "C1_B", "C0_B",
                     "phi_num", "total_weight", "M")
  expect_true(all(required_suff %in% names(out$suff)))
})

test_that("e_step_fmm: phi_num / total_weight is in (0,1)", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  phi_implied <- out$suff$phi_num / out$suff$total_weight
  expect_true(phi_implied > 0 && phi_implied < 1)
})

test_that("e_step_fmm: type suff stats are non-negative", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  out    <- e_step_fmm(df, params)
  s <- out$suff
  expect_true(all(c(s$T11_A, s$T01_A, s$D1_A, s$D0_A) >= 0))
  expect_true(all(c(s$T11_B, s$T01_B, s$D1_B, s$D0_B) >= 0))
})

# ---- phi -> 1: approximates single-type baseline --------------------------

test_that("FMM with phi~1 approximates single-type baseline LL", {
  df     <- .make_fmm_panel(n = 300L, phi = 0.999)
  # Force phi very close to 1 (nearly single-type)
  params_fmm <- init_params_fmm()
  params_fmm$phi <- 0.99
  params_fmm$theta0_B <- params_fmm$theta0_A
  params_fmm$theta1_B <- params_fmm$theta1_A
  params_fmm$alpha_B  <- params_fmm$alpha_A
  ll_fmm <- compute_observed_loglik_fmm(df, params_fmm)
  expect_true(is.finite(ll_fmm))
})

# ---- m_step_fmm ------------------------------------------------------------

test_that("m_step_fmm returns all required fields", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  suff   <- e_step_fmm(df, params)$suff
  out    <- m_step_fmm(suff)
  required <- c("theta0_A", "theta1_A", "alpha_A",
                "theta0_B", "theta1_B", "alpha_B", "phi", "pi")
  expect_true(all(required %in% names(out)))
})

test_that("m_step_fmm: transition probabilities in (0, theta_cap)", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  suff   <- e_step_fmm(df, params)$suff
  out    <- m_step_fmm(suff, theta_cap = 0.999)
  for (nm in c("theta0_A", "theta1_A", "theta0_B", "theta1_B")) {
    expect_true(out[[nm]] > 0 && out[[nm]] <= 0.999,
                label = paste(nm, "must be in (0, 0.999]"))
  }
})

test_that("m_step_fmm: pi in (0, pi_cap)", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  suff   <- e_step_fmm(df, params)$suff
  out    <- m_step_fmm(suff, pi_cap = 0.49)
  expect_true(out$pi >= 0 && out$pi <= 0.49)
})

test_that("m_step_fmm: phi in (0, 1)", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  suff   <- e_step_fmm(df, params)$suff
  out    <- m_step_fmm(suff)
  expect_true(out$phi > 0 && out$phi < 1)
})

# ---- resolve_label_switching_fmm ------------------------------------------

test_that("resolve_label_switching swaps when theta1_B > theta1_A", {
  params <- init_params_fmm()
  # Force a flipped situation
  params$theta1_A <- 0.75
  params$theta1_B <- 0.92
  params$phi      <- 0.4
  resolved        <- resolve_label_switching_fmm(params)
  expect_true(resolved$theta1_A >= resolved$theta1_B)
  expect_equal(resolved$phi, 1 - params$phi)
})

test_that("resolve_label_switching is identity when already labelled correctly", {
  params  <- init_params_fmm()   # theta1_A > theta1_B by default
  resolved <- resolve_label_switching_fmm(params)
  expect_equal(resolved$theta1_A, params$theta1_A)
  expect_equal(resolved$phi,      params$phi)
})

# ---- em_fit_fmm: convergence -----------------------------------------------

test_that("em_fit_fmm converges on synthetic two-type data", {
  df  <- .make_fmm_panel(n = 300L)
  fit <- em_fit_fmm(df, model_type = "symmetric", verbose = 0L)
  expect_true(fit$converged)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_fmm: LL is monotone non-decreasing", {
  df  <- .make_fmm_panel(n = 300L)
  fit <- em_fit_fmm(df, verbose = 0L)
  ll_vec <- fit$history$loglik
  diffs  <- diff(ll_vec)
  expect_true(all(diffs >= -1e-5), label = "FMM LL must be non-decreasing")
})

test_that("em_fit_fmm: no-error variant converges", {
  df  <- .make_fmm_panel(n = 300L, pi = 0)
  fit <- em_fit_fmm(df, model_type = "none", verbose = 0L)
  expect_true(is.finite(fit$loglik))
})

test_that("em_fit_fmm: result has correct gamma dimensions", {
  df  <- .make_fmm_panel(n = 100L)
  fit <- em_fit_fmm(df, verbose = 0L)
  expect_equal(nrow(fit$gamma_A), 100L)
  expect_equal(ncol(fit$gamma_A), 8L)
  expect_equal(dim(fit$gamma_B), dim(fit$gamma_A))
})

test_that("em_fit_fmm: label switching resolved (type A more stable)", {
  df  <- .make_fmm_panel(n = 300L)
  fit <- em_fit_fmm(df, verbose = 0L)
  expect_true(fit$params$theta1_A >= fit$params$theta1_B)
})

# ---- compute_observed_loglik_fmm -------------------------------------------

test_that("compute_observed_loglik_fmm returns finite scalar", {
  df     <- .make_fmm_panel()
  params <- init_params_fmm()
  ll     <- compute_observed_loglik_fmm(df, params)
  expect_true(is.finite(ll))
  expect_length(ll, 1L)
})

# ---- em_fit_fmm: additional convergence guards ---------------------------

test_that("em_fit_fmm runs more than 1 iteration", {
  df  <- .make_fmm_panel(n = 300L)
  fit <- em_fit_fmm(df, verbose = 0L)
  expect_true(fit$iterations > 1L,
              label = "FMM EM must iterate more than once")
})

test_that("em_fit_fmm: params change from starting values", {
  df  <- .make_fmm_panel(n = 300L)
  p0  <- init_params_fmm()
  fit <- em_fit_fmm(df, params0 = p0, verbose = 0L)
  expect_false(identical(fit$params$phi, p0$phi),
               label = "phi must change from initial value")
})

test_that("em_fit_fmm: phi stays in (0, 1) throughout all iterations", {
  df  <- .make_fmm_panel(n = 300L)
  fit <- em_fit_fmm(df, verbose = 0L, max_iter = 200L)
  phi_history <- fit$history$phi
  expect_true(all(phi_history > 0 & phi_history < 1),
              label = "phi must remain in (0, 1) at every iteration")
})
