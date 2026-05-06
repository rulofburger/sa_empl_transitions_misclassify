# ==============================================================================
# EM-AR2 tests: estep.R
# ==============================================================================

# Helper: create small synthetic dataset
.make_synthetic_df <- function(N = 50, seed = 42) {
  set.seed(seed)
  data.frame(
    y1 = sample(0:1, N, replace = TRUE),
    y2 = sample(0:1, N, replace = TRUE),
    y3 = sample(0:1, N, replace = TRUE),
    y4 = sample(0:1, N, replace = TRUE),
    weight = rep(1, N)
  )
}

# Helper: typical params
.default_params <- function(pi = 0.05) {
  list(theta0=0.1, theta01=0.2, theta1=0.08, theta10=0.15, pi=pi)
}

test_that("e_step_ar2 returns correct shapes", {
  df    <- .make_synthetic_df(30)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)

  expect_named(out, c("gamma", "loglik", "sufficient_stats"))
  expect_equal(dim(out$gamma), c(30L, 16L))
  expect_true(is.numeric(out$loglik))
  expect_true(is.finite(out$loglik))
})

test_that("e_step_ar2 responsibilities sum to 1 for each observation", {
  df    <- .make_synthetic_df(20)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  row_sums <- rowSums(out$gamma)
  expect_equal(row_sums, rep(1, 20), tolerance = 1e-10)
})

test_that("e_step_ar2 with pi=0: responsibility is 1 for matching history", {
  hmat <- latent_histories_ar2()

  # Use row 5 of hmat as the observed sequence (0,0,1,0)
  h_row <- 5L
  ymat  <- matrix(hmat[h_row, ], nrow = 1)
  w     <- 1
  params <- .default_params(pi = 0)
  params$pi <- 1e-12  # effectively zero

  out <- e_step_ar2(ymat, w, hmat, params)

  # The responsibility for the matching history should be essentially 1
  expect_gt(out$gamma[1, h_row], 0.999)
  expect_lt(sum(out$gamma[1, -h_row]), 0.001)
})

test_that("e_step_ar2 log-likelihood is finite and negative or zero", {
  df    <- .make_synthetic_df(100)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  expect_true(is.finite(out$loglik))
  expect_lte(out$loglik, 0)
})

test_that("e_step_ar2 sufficient stats D, T1, S are non-negative", {
  df    <- .make_synthetic_df(50)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  ss  <- out$sufficient_stats

  expect_true(all(ss$D  >= 0))
  expect_true(all(ss$T1 >= 0))
  expect_true(all(ss$S  >= 0))
})

test_that("e_step_ar2 T1[j,k] <= D[j,k] for all (j,k)", {
  df    <- .make_synthetic_df(50)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  ss  <- out$sufficient_stats

  expect_true(all(ss$T1 <= ss$D + 1e-10))
})

test_that("e_step_ar2 M is non-negative", {
  df    <- .make_synthetic_df(50)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  expect_gte(out$sufficient_stats$M, 0)
})

test_that("e_step_ar2 with uniform weights == equal weights", {
  df     <- .make_synthetic_df(20)
  hmat   <- latent_histories_ar2()
  ymat   <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out1 <- e_step_ar2(ymat, df$weight, hmat, params)           # weight=1
  out2 <- e_step_ar2(ymat, 2 * df$weight, hmat, params)       # weight=2

  # Responsibilities should be identical (weights don't change normalisation per obs)
  expect_equal(out1$gamma, out2$gamma, tolerance = 1e-10)
  # LL should scale by 2
  expect_equal(out2$loglik, 2 * out1$loglik, tolerance = 1e-8)
})

test_that("e_step_ar2 sufficient_stats S sums to total weight", {
  df    <- .make_synthetic_df(40)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- .default_params()

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  # S_{jk} sums over all (j,k) should equal total weight
  expect_equal(sum(out$sufficient_stats$S), sum(df$weight), tolerance = 1e-10)
})

test_that("e_step_ar2 asymmetric mode returns M0 and M1", {
  df    <- .make_synthetic_df(30)
  hmat  <- latent_histories_ar2()
  ymat  <- as.matrix(df[, c("y1","y2","y3","y4")])
  params <- list(theta0=0.1, theta01=0.2, theta1=0.08, theta10=0.15,
                 pi0=0.04, pi1=0.06, asymmetric=TRUE)

  out <- e_step_ar2(ymat, df$weight, hmat, params)
  ss  <- out$sufficient_stats

  expect_true(is.numeric(ss$M0))
  expect_true(is.numeric(ss$M1))
  expect_gte(ss$M0, 0)
  expect_gte(ss$M1, 0)
})
