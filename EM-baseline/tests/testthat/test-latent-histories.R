# ==============================================================================
# EM-baseline: Tests for latent_histories.R
# Created: 2026-05-05
# ==============================================================================

test_that("latent_histories returns 8 x 3 integer matrix", {
  hmat <- latent_histories()
  expect_equal(dim(hmat), c(8L, 3L))
  expect_true(all(hmat %in% c(0L, 1L)))
})

test_that("latent_histories contains all binary combinations", {
  hmat <- latent_histories()
  # Each row should be unique
  expect_equal(nrow(unique(hmat)), 8L)
})

test_that("prior_over_histories sums to 1", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 0.5)
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("prior_over_histories returns non-negative probabilities", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 0.8, theta0 = 0.15, alpha = 0.4)
  expect_true(all(p >= 0))
})

test_that("prior_over_histories correct for simple case (theta1=1, theta0=0)", {
  hmat <- latent_histories()
  # theta1=1: once employed, always employed; theta0=0: never enter employment
  # alpha=0.5: so h1=1 -> must have h2=1, h3=1 (prob = 0.5)
  #            h1=0 -> must have h2=0, h3=0 (prob = 0.5)
  p <- prior_over_histories(hmat, theta1 = 1 - 1e-10, theta0 = 1e-10, alpha = 0.5)
  # History (1,1,1) should have ~0.5, history (0,0,0) ~0.5, rest ~0
  idx_111 <- which(hmat[,1]==1 & hmat[,2]==1 & hmat[,3]==1)
  idx_000 <- which(hmat[,1]==0 & hmat[,2]==0 & hmat[,3]==0)
  other   <- setdiff(1:8, c(idx_111, idx_000))
  expect_equal(p[idx_111] + p[idx_000], 1, tolerance = 1e-6)
  expect_true(all(p[other] < 1e-6))
})

test_that("stationary_alpha equals theta0/(theta0 + 1 - theta1)", {
  theta0 <- 0.12
  theta1 <- 0.88
  expected <- theta0 / (theta0 + 1 - theta1)
  expect_equal(stationary_alpha(theta0, theta1), expected, tolerance = 1e-12)
})

test_that("prior marginal at wave 1 equals alpha under stationarity", {
  theta0 <- 0.10
  theta1 <- 0.90
  alpha  <- stationary_alpha(theta0, theta1)
  hmat   <- latent_histories()
  p      <- prior_over_histories(hmat, theta1, theta0, alpha)
  # Marginal P(h1 = 1) = sum of prior over histories where h1 = 1
  p_empl_wave1 <- sum(p[hmat[, 1] == 1])
  expect_equal(p_empl_wave1, alpha, tolerance = 1e-10)
})

test_that("prior_over_histories sums to 1 when theta0 near 0", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 1e-8, alpha = 0.5)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_true(all(p >= 0))
})

test_that("prior_over_histories sums to 1 when theta1 near 1", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 1 - 1e-8, theta0 = 0.1, alpha = 0.5)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_true(all(p >= 0))
})

test_that("prior_over_histories sums to 1 when alpha near 0", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 1e-8)
  expect_equal(sum(p), 1, tolerance = 1e-10)
})

test_that("prior_over_histories sums to 1 when alpha near 1", {
  hmat <- latent_histories()
  p <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 1 - 1e-8)
  expect_equal(sum(p), 1, tolerance = 1e-10)
})
