# ==============================================================================
# EM-baseline: Tests for mstep.R
# Created: 2026-05-05
# ==============================================================================

# ---- helper: build hand-crafted sufficient stats ----------------------------

.make_suff <- function(T11 = 60, D1 = 80, T01 = 15, D0 = 80,
                       C1 = 55, C0 = 45,
                       M = 6, M0 = 3, M1 = 3, H0 = 120, H1 = 120,
                       total_weight = 100) {
  list(T11 = T11, D1 = D1, T01 = T01, D0 = D0,
       C1 = C1, C0 = C0,
       M = M, M0 = M0, M1 = M1, H0 = H0, H1 = H1,
       total_weight = total_weight)
}

# ---- tests: transition probabilities ----------------------------------------

test_that("m_step computes theta1 = T11/D1 and theta0 = T01/D0", {
  suff   <- .make_suff(T11 = 72, D1 = 80, T01 = 8, D0 = 80)
  result <- m_step(suff, model_type = "none", stationary = FALSE)
  expect_equal(result$theta1, 72 / 80, tolerance = 1e-12)
  expect_equal(result$theta0, 8 / 80,  tolerance = 1e-12)
})

test_that("m_step clamps theta to theta_cap", {
  suff   <- .make_suff(T11 = 80, D1 = 80, T01 = 80, D0 = 80)  # would give 1.0
  result <- m_step(suff, model_type = "none", theta_cap = 0.999)
  expect_true(result$theta1 <= 0.999)
  expect_true(result$theta0 <= 0.999)
})

test_that("m_step with stationary=TRUE enforces alpha = theta0/(theta0+1-theta1)", {
  suff   <- .make_suff(T11 = 60, D1 = 80, T01 = 15, D0 = 80)
  result <- m_step(suff, model_type = "none", stationary = TRUE)
  expected_alpha <- result$theta0 / (result$theta0 + 1 - result$theta1)
  expect_equal(result$alpha, expected_alpha, tolerance = 1e-12)
})

test_that("m_step with stationary=FALSE uses alpha = C1/(C1+C0)", {
  suff   <- .make_suff(C1 = 60, C0 = 40)
  result <- m_step(suff, model_type = "none", stationary = FALSE)
  expect_equal(result$alpha, 60 / 100, tolerance = 1e-12)
})

# ---- tests: symmetric misclassification -------------------------------------

test_that("m_step symmetric computes pi = M/(3*W)", {
  W      <- 100
  M      <- 12
  suff   <- .make_suff(M = M, total_weight = W)
  result <- m_step(suff, model_type = "symmetric", stationary = TRUE)
  expect_equal(result$pi, M / (3 * W), tolerance = 1e-12)
})

test_that("m_step symmetric clamps pi to pi_cap", {
  suff   <- .make_suff(M = 200, total_weight = 100)  # would give > 0.49
  result <- m_step(suff, model_type = "symmetric", pi_cap = 0.49)
  expect_true(result$pi <= 0.49)
})

test_that("m_step symmetric has pi field but no pi0/pi1", {
  suff   <- .make_suff()
  result <- m_step(suff, model_type = "symmetric")
  expect_true("pi" %in% names(result))
  expect_false("pi0" %in% names(result))
  expect_false("pi1" %in% names(result))
})

# ---- tests: asymmetric misclassification ------------------------------------

test_that("m_step asymmetric computes pi0 = M0/H0 and pi1 = M1/H1", {
  suff   <- .make_suff(M0 = 5, H0 = 100, M1 = 8, H1 = 100)
  result <- m_step(suff, model_type = "asymmetric", stationary = TRUE)
  expect_equal(result$pi0, 5 / 100, tolerance = 1e-12)
  expect_equal(result$pi1, 8 / 100, tolerance = 1e-12)
})

test_that("m_step asymmetric has pi0 and pi1 fields but no pi", {
  suff   <- .make_suff()
  result <- m_step(suff, model_type = "asymmetric")
  expect_true("pi0" %in% names(result))
  expect_true("pi1" %in% names(result))
  expect_false("pi" %in% names(result))
})

# ---- tests: no misclassification --------------------------------------------

test_that("m_step none has no pi/pi0/pi1 fields", {
  suff   <- .make_suff()
  result <- m_step(suff, model_type = "none")
  expect_false("pi"  %in% names(result))
  expect_false("pi0" %in% names(result))
  expect_false("pi1" %in% names(result))
})

# ---- tests: numerical edge cases --------------------------------------------

test_that("m_step handles D1=0 without NaN", {
  suff   <- .make_suff(T11 = 0, D1 = 0)  # degenerate: no from-employed periods
  result <- m_step(suff, model_type = "symmetric")
  expect_true(is.finite(result$theta1))
  expect_true(result$theta1 >= 0)
})

test_that("m_step handles zero misclassification count (M=0)", {
  suff   <- .make_suff(M = 0, M0 = 0, M1 = 0)
  result <- m_step(suff, model_type = "symmetric")
  expect_equal(result$pi, 0, tolerance = 1e-12)
})

test_that("m_step errors on unknown model_type", {
  suff <- .make_suff()
  expect_error(m_step(suff, model_type = "fmm"), regexp = "model_type")
})
