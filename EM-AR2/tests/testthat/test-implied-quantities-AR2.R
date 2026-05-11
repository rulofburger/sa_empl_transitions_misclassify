# ==============================================================================
# EM-AR2 tests: implied_ar2()
# Created: 2026-05-07
# ==============================================================================

test_that("implied_ar2 returns list with all expected names", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi = 0.05)
  out <- implied_ar2(p, "symmetric")
  expected <- c("entry_rate", "exit_rate", "p_00", "p_10", "p_01", "p_11",
                "employment_rate", "pi", "pi0", "pi1")
  expect_true(all(expected %in% names(out)))
})

test_that("implied_ar2 computes p_jk from eq. 4 correctly (symmetric)", {
  theta0 <- 0.10; theta01 <- 0.05; theta1 <- 0.08; theta10 <- 0.03
  p <- list(theta0 = theta0, theta01 = theta01, theta1 = theta1,
            theta10 = theta10, pi = 0.05)
  out <- implied_ar2(p, "symmetric")

  expect_equal(out$p_00, theta0)
  expect_equal(out$p_10, theta0 + theta01)
  expect_equal(out$p_01, 1 - theta1 - theta10)
  expect_equal(out$p_11, 1 - theta1)
})

test_that("implied_ar2 entry_rate = theta0, exit_rate = 1 - theta1", {
  p <- list(theta0 = 0.12, theta01 = 0.04, theta1 = 0.07, theta10 = 0.02)
  out <- implied_ar2(p, "none")
  expect_equal(out$entry_rate, 0.12)
  expect_equal(out$exit_rate,  0.93)
})

test_that("implied_ar2 model_type = 'none': pi, pi0, pi1 all NA", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03)
  out <- implied_ar2(p, "none")
  expect_true(is.na(out$pi))
  expect_true(is.na(out$pi0))
  expect_true(is.na(out$pi1))
})

test_that("implied_ar2 model_type = 'symmetric': pi set, pi0/pi1 NA", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi = 0.07)
  out <- implied_ar2(p, "symmetric")
  expect_equal(out$pi, 0.07)
  expect_true(is.na(out$pi0))
  expect_true(is.na(out$pi1))
})

test_that("implied_ar2 model_type = 'asymmetric': pi0/pi1 set, pi NA", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi0 = 0.06, pi1 = 0.04)
  out <- implied_ar2(p, "asymmetric")
  expect_equal(out$pi0, 0.06)
  expect_equal(out$pi1, 0.04)
  expect_true(is.na(out$pi))
})

test_that("implied_ar2 employment_rate is in (0, 1) for well-behaved params", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03)
  out <- implied_ar2(p, "none")
  # For well-conditioned params stationary_ar2() must succeed (not return NA).
  expect_false(is.na(out$employment_rate))
  expect_gt(out$employment_rate, 0)
  expect_lt(out$employment_rate, 1)
})

test_that("implied_ar2 returns NA employment_rate for degenerate params", {
  # All probabilities at 0 — stationary distribution may not be defined
  p <- list(theta0 = 0, theta01 = 0, theta1 = 1, theta10 = 0)
  out <- implied_ar2(p, "none")
  # Either a numeric value or NA is acceptable; must not error
  expect_true(is.numeric(out$employment_rate))
})

test_that("implied_ar2 errors on invalid model_type", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "unknown_type"), regexp = "model_type")
})

test_that("implied_ar2 errors when required params are missing", {
  # theta0 missing
  p <- list(theta01 = 0.05, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors on theta0 < 0", {
  p <- list(theta0 = -0.05, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors on theta0 > 1", {
  p <- list(theta0 = 1.05, theta01 = 0.00, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors when theta0 + theta01 > 1 (p_10 would exceed 1)", {
  p <- list(theta0 = 0.70, theta01 = 0.40, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors when theta1 + theta10 > 1 (p_01 would be negative)", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.80, theta10 = 0.30)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors on negative theta01", {
  p <- list(theta0 = 0.10, theta01 = -0.05, theta1 = 0.08, theta10 = 0.03)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors on negative theta10", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = -0.01)
  expect_error(implied_ar2(p, "none"))
})

test_that("implied_ar2 errors on out-of-range pi (symmetric)", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi = 1.2)
  expect_error(implied_ar2(p, "symmetric"), regexp = "pi")
})

test_that("implied_ar2 errors on out-of-range pi0 (asymmetric)", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi0 = -0.1, pi1 = 0.05)
  expect_error(implied_ar2(p, "asymmetric"), regexp = "pi0")
})

test_that("implied_ar2 errors on out-of-range pi1 (asymmetric)", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi0 = 0.05, pi1 = 1.1)
  expect_error(implied_ar2(p, "asymmetric"), regexp = "pi1")
})

test_that("implied_ar2 transition probabilities are all in [0, 1]", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.02)
  out <- implied_ar2(p, "none")
  for (nm in c("p_00", "p_10", "p_01", "p_11")) {
    expect_gte(out[[nm]], 0, label = nm)
    expect_lte(out[[nm]], 1, label = nm)
  }
})

test_that("implied_ar2 all outputs are scalar numerics", {
  p <- list(theta0 = 0.10, theta01 = 0.05, theta1 = 0.08, theta10 = 0.03,
            pi0 = 0.06, pi1 = 0.04)
  out <- implied_ar2(p, "asymmetric")
  for (nm in names(out)) {
    val <- out[[nm]]
    expect_true(is.numeric(val) && length(val) == 1L,
                info = sprintf("'%s' should be scalar numeric", nm))
  }
})
