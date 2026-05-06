# ==============================================================================
# EM-AR2 tests: mstep.R
# ==============================================================================

# Helper: build artificial sufficient statistics from known parameters
# (as if the EM had perfect data and exact responsibilities)
.make_perfect_ss <- function(theta0, theta01, theta1, theta10, W = 1000) {
  # True transition probabilities
  p_00_1 <- theta0
  p_10_1 <- theta0 + theta01
  p_01_1 <- 1 - theta1 - theta10
  p_11_1 <- 1 - theta1

  # Artificial D_{jk}: distribute W evenly across the 4 states
  D <- matrix(W / 4, nrow = 2, ncol = 2,
              dimnames = list(c("0","1"), c("0","1")))

  # T1_{jk} = p_{jk->1} * D_{jk}
  T1 <- matrix(c(p_00_1, p_10_1, p_01_1, p_11_1) * W / 4,
               nrow = 2, ncol = 2,
               dimnames = list(c("0","1"), c("0","1")))

  # S: uniform initial pair distribution
  S <- matrix(W / 4, nrow = 2, ncol = 2,
              dimnames = list(c("0","1"), c("0","1")))

  # M = expected misclassification: use pi=0.05, M = 0.05 * 4 * W
  list(S=S, D=D, T1=T1, M=0.05 * 4 * W, M0=NA_real_, M1=NA_real_, W=W)
}

# .make_perfect_ss now includes S; .make_perfect_ss_with_S is an alias for compatibility
.make_perfect_ss_with_S <- .make_perfect_ss

test_that("m_step_ar2 recovers true parameters from perfect sufficient stats", {
  theta0  <- 0.10; theta01 <- 0.15
  theta1  <- 0.08; theta10 <- 0.12

  ss  <- .make_perfect_ss_with_S(theta0, theta01, theta1, theta10)
  out <- m_step_ar2(ss, estimate_pi = TRUE)

  expect_equal(out$theta0,  theta0,  tolerance = 1e-8)
  expect_equal(out$theta01, theta01, tolerance = 1e-8)
  expect_equal(out$theta1,  theta1,  tolerance = 1e-8)
  expect_equal(out$theta10, theta10, tolerance = 1e-8)
})

test_that("m_step_ar2 pi estimate = M/(4W)", {
  theta0 <- 0.10; theta01 <- 0.15; theta1 <- 0.08; theta10 <- 0.12
  W <- 500
  ss <- .make_perfect_ss_with_S(theta0, theta01, theta1, theta10, W = W)
  ss$M <- 0.07 * 4 * W   # set M to match pi=0.07

  out <- m_step_ar2(ss, estimate_pi = TRUE)
  expect_equal(out$pi, 0.07, tolerance = 1e-8)
})

test_that("m_step_ar2 with estimate_pi=FALSE returns fixed_pi", {
  ss  <- .make_perfect_ss_with_S(0.1, 0.15, 0.08, 0.12)
  out <- m_step_ar2(ss, estimate_pi = FALSE, fixed_pi = 0)
  expect_equal(out$pi, 0)
})

test_that("m_step_ar2 returns alpha", {
  ss  <- .make_perfect_ss_with_S(0.1, 0.15, 0.08, 0.12)
  out <- m_step_ar2(ss, estimate_pi = TRUE)
  expect_true("alpha" %in% names(out))
  expect_equal(sum(out$alpha), 1, tolerance = 1e-10)
  expect_true(all(out$alpha >= 0))
})

test_that("m_step_ar2 all returned probabilities are in (0,1)", {
  ss  <- .make_perfect_ss_with_S(0.1, 0.15, 0.08, 0.12)
  out <- m_step_ar2(ss, estimate_pi = TRUE)

  expect_gt(out$theta0,  0); expect_lt(out$theta0,  1)
  expect_gt(out$theta01, 0); expect_lt(out$theta01, 1)
  expect_gt(out$theta1,  0); expect_lt(out$theta1,  1)
  expect_gt(out$theta10, 0); expect_lt(out$theta10, 1)
  expect_gt(out$pi,      0); expect_lt(out$pi,      1)
})

test_that("m_step_ar2 theta01 identity: p10_1 - p00_1 = theta01", {
  ss  <- .make_perfect_ss_with_S(0.1, 0.18, 0.08, 0.12)
  out <- m_step_ar2(ss, estimate_pi = TRUE)

  expect_equal(out$.p_10_1 - out$.p_00_1, out$theta01, tolerance = 1e-8)
})

test_that("m_step_ar2 theta10 identity: p11_1 - p01_1 = theta10", {
  ss  <- .make_perfect_ss_with_S(0.1, 0.18, 0.08, 0.14)
  out <- m_step_ar2(ss, estimate_pi = TRUE)

  expect_equal(out$.p_11_1 - out$.p_01_1, out$theta10, tolerance = 1e-8)
})

test_that("m_step_ar2 pi is capped at pi_cap", {
  ss    <- .make_perfect_ss_with_S(0.1, 0.15, 0.08, 0.12)
  ss$M  <- 0.6 * 4 * ss$W   # pi would be 0.6 without cap
  out   <- m_step_ar2(ss, estimate_pi = TRUE, pi_cap = 0.49)
  expect_lte(out$pi, 0.49)
})

test_that("m_step_ar2 asymmetric mode returns pi0 and pi1", {
  ss <- .make_perfect_ss(0.1, 0.15, 0.08, 0.12, W = 1000)
  # Attach M0 and M1 with exposure attributes
  M0 <- 30.0; attr(M0, "exposure") <- 400
  M1 <- 20.0; attr(M1, "exposure") <- 600
  ss$M0 <- M0; ss$M1 <- M1

  out <- m_step_ar2(ss, asymmetric = TRUE)

  expect_true("pi0" %in% names(out))
  expect_true("pi1" %in% names(out))
  expect_equal(out$pi0, 30 / 400, tolerance = 1e-8)
  expect_equal(out$pi1, 20 / 600, tolerance = 1e-8)
})
