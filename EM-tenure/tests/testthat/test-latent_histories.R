# ==============================================================================
# Tests for EM-tenure: latent_histories.R
# ==============================================================================

test_that("latent_histories returns 8x3 binary matrix", {
  hmat <- latent_histories()
  expect_equal(dim(hmat), c(8, 3))
  expect_true(all(hmat %in% c(0, 1)))
})

test_that("latent_histories has all 2^3 combinations", {
  hmat <- latent_histories()
  # Encode as integers (0-7)
  codes <- hmat[, 1] * 4 + hmat[, 2] * 2 + hmat[, 3]
  expect_equal(sort(codes), 0:7)
})

test_that("prior_over_histories sums to 1", {
  hmat <- latent_histories()
  prior <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 0.6)
  expect_equal(sum(prior), 1, tolerance = 1e-10)
})

test_that("prior_over_histories with alpha=1 gives zero weight to h1=0", {
  hmat <- latent_histories()
  prior <- prior_over_histories(hmat, theta1 = 0.9, theta0 = 0.1, alpha = 1 - 1e-10)
  expect_true(all(prior[hmat[, 1] == 0] < 1e-8))
})

test_that("clocks_from_histories are deterministic and correct", {
  hmat <- latent_histories()
  clocks <- clocks_from_histories(hmat)

  # History (1,1,1): tenure = 0.25, 0.50, 0.75
  idx <- which(hmat[, 1] == 1 & hmat[, 2] == 1 & hmat[, 3] == 1)
  expect_equal(as.numeric(clocks$Gstar[idx, ]), c(0.25, 0.50, 0.75), tolerance = 1e-12)
  expect_equal(as.numeric(clocks$Dstar[idx, ]), c(0, 0, 0))

  # History (0,0,0): timegap = 0.25, 0.50, 0.75
  idx <- which(hmat[, 1] == 0 & hmat[, 2] == 0 & hmat[, 3] == 0)
  expect_equal(as.numeric(clocks$Dstar[idx, ]), c(0.25, 0.50, 0.75), tolerance = 1e-12)
  expect_equal(as.numeric(clocks$Gstar[idx, ]), c(0, 0, 0))

  # History (1,0,1): tenure = 0.25, 0, 0.25; timegap = 0, 0.25, 0
  idx <- which(hmat[, 1] == 1 & hmat[, 2] == 0 & hmat[, 3] == 1)
  expect_equal(as.numeric(clocks$Gstar[idx, ]), c(0.25, 0, 0.25), tolerance = 1e-12)
  expect_equal(as.numeric(clocks$Dstar[idx, ]), c(0, 0.25, 0), tolerance = 1e-12)
})
