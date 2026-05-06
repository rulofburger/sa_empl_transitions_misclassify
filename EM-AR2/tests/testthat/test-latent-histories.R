# ==============================================================================
# EM-AR2 tests: latent_histories.R
# ==============================================================================

test_that("latent_histories_ar2 returns 16 x 4 binary matrix", {
  hmat <- latent_histories_ar2()
  expect_equal(dim(hmat), c(16L, 4L))
  expect_true(all(hmat %in% c(0L, 1L)))
})

test_that("latent_histories_ar2 contains all 16 unique histories", {
  hmat <- latent_histories_ar2()
  unique_rows <- unique(hmat)
  expect_equal(nrow(unique_rows), 16L)
})

test_that("ar2_trans_prob returns correct probabilities for known cases", {
  theta0 <- 0.1; theta01 <- 0.2; theta1 <- 0.08; theta10 <- 0.15

  # P(h_t=1 | 0,0) = theta0
  expect_equal(ar2_trans_prob(0, 0, theta0, theta01, theta1, theta10), theta0)

  # P(h_t=1 | 1,0) = theta0 + theta01
  expect_equal(ar2_trans_prob(1, 0, theta0, theta01, theta1, theta10), theta0 + theta01)

  # P(h_t=1 | 0,1) = 1 - theta1 - theta10
  expect_equal(ar2_trans_prob(0, 1, theta0, theta01, theta1, theta10), 1 - theta1 - theta10)

  # P(h_t=1 | 1,1) = 1 - theta1
  expect_equal(ar2_trans_prob(1, 1, theta0, theta01, theta1, theta10), 1 - theta1)
})

test_that("ar2_trans_prob is vectorised", {
  result <- ar2_trans_prob(c(0,1,0,1), c(0,0,1,1), 0.1, 0.2, 0.08, 0.15)
  expect_length(result, 4)
  expect_equal(result[1], 0.1)
  expect_equal(result[2], 0.3)
  expect_equal(result[3], 1 - 0.08 - 0.15)
  expect_equal(result[4], 1 - 0.08)
})

test_that("stationary_ar2 sums to 1 and is non-negative", {
  alpha <- stationary_ar2(theta0=0.1, theta01=0.2, theta1=0.08, theta10=0.15)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(all(alpha >= 0))
})

test_that("stationary_ar2 has correct names", {
  alpha <- stationary_ar2(0.1, 0.2, 0.08, 0.15)
  expect_setequal(names(alpha), c("00", "10", "01", "11"))
})

test_that("stationary_ar2 reduces to AR(1) ergodic when theta01=theta10=0", {
  # With theta01=theta10=0, the AR(2) collapses to AR(1).
  # AR(1) ergodic: mu = theta0 / (theta0 + theta1)
  # alpha(1,1) + alpha(0,1) = mu (share employed in steady state)
  theta0 <- 0.12; theta1 <- 0.08
  # Use small but non-zero theta01/theta10 to avoid Phi=0 singularity
  alpha <- stationary_ar2(theta0, 0, theta1, 0)
  mu_ar1 <- theta0 / (theta0 + theta1)
  mu_from_alpha <- alpha["01"] + alpha["11"]
  expect_equal(as.numeric(mu_from_alpha), mu_ar1, tolerance = 1e-8)
})

test_that("stationary_ar2 errors when Phi is near zero", {
  # Phi = theta0*(theta10-1) + theta1*(theta01-1) = 0 when both increments are 1
  # This should throw an informative error
  expect_error(stationary_ar2(0.1, 1, 0.08, 1), regexp = "Phi")
})

test_that("prior_over_histories_ar2 sums to 1", {
  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15)
  expect_equal(sum(prior), 1, tolerance = 1e-10)
})

test_that("prior_over_histories_ar2 is non-negative", {
  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15)
  expect_true(all(prior >= 0))
})

test_that("prior_over_histories_ar2 returns length-16 vector", {
  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15)
  expect_length(prior, 16L)
})

test_that("prior_over_histories_ar2 marginal h1,h2 matches alpha when supplied", {
  hmat  <- latent_histories_ar2()
  theta0 <- 0.1; theta01 <- 0.2; theta1 <- 0.08; theta10 <- 0.15
  alpha <- stationary_ar2(theta0, theta01, theta1, theta10)
  # Pass alpha explicitly
  prior <- prior_over_histories_ar2(hmat, theta0, theta01, theta1, theta10, alpha = alpha)

  # Marginalise prior over h3, h4 to get joint distribution of (h1,h2)
  h1  <- hmat[, 1]; h2 <- hmat[, 2]
  for (j in 0:1) {
    for (k in 0:1) {
      marginal_jk <- sum(prior[(h1 == j) & (h2 == k)])
      alpha_jk    <- alpha[paste0(j, k)]
      expect_equal(marginal_jk, as.numeric(alpha_jk), tolerance = 1e-8)
    }
  }
})

test_that("prior_over_histories_ar2 with NULL alpha uses stationary", {
  hmat  <- latent_histories_ar2()
  theta0 <- 0.1; theta01 <- 0.2; theta1 <- 0.08; theta10 <- 0.15
  # NULL alpha should give same result as passing stationary
  alpha <- stationary_ar2(theta0, theta01, theta1, theta10)
  prior_null <- prior_over_histories_ar2(hmat, theta0, theta01, theta1, theta10)
  prior_stat <- prior_over_histories_ar2(hmat, theta0, theta01, theta1, theta10, alpha = alpha)
  expect_equal(prior_null, prior_stat, tolerance = 1e-10)
})

test_that("prior_over_histories_ar2 errors on unnamed alpha", {
  hmat <- latent_histories_ar2()
  unnamed_alpha <- c(0.25, 0.25, 0.25, 0.25)  # no names
  expect_error(
    prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15, alpha = unnamed_alpha),
    regexp = "alpha must be a named vector"
  )
})

test_that("prior_over_histories_ar2 errors on wrong alpha names", {
  hmat <- latent_histories_ar2()
  bad_alpha <- c(a = 0.25, b = 0.25, c = 0.25, d = 0.25)  # wrong names
  expect_error(
    prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15, alpha = bad_alpha),
    regexp = "alpha must be a named vector with names"
  )
})

test_that("prior_over_histories_ar2 errors on negative alpha components", {
  hmat <- latent_histories_ar2()
  neg_alpha <- c(`00` = -0.1, `10` = 0.4, `01` = 0.4, `11` = 0.3)
  expect_error(
    prior_over_histories_ar2(hmat, 0.1, 0.2, 0.08, 0.15, alpha = neg_alpha),
    regexp = "alpha contains negative values"
  )
})
