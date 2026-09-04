make_hybrid_test_data <- function(n = 80L, seed = 882L) {
  set.seed(seed)
  y <- sample(0:1, n, replace = TRUE, prob = c(.35, .65))
  g1 <- pmax(.26, rexp(n, 2))
  g2 <- ifelse(runif(n) < .6, g1 + .25, pmax(.26, rexp(n, 2)))
  g3 <- ifelse(runif(n) < .6, g2 + .25, pmax(.26, rexp(n, 2)))
  data.frame(y1 = y, y2 = y, y3 = y, tenure1 = g1, tenure2 = g2,
    tenure3 = g3, timegap_cat1 = sample(1:7, n, TRUE),
    timegap_cat2 = sample(1:7, n, TRUE),
    timegap_cat3 = sample(1:7, n, TRUE), weight = 1)
}

make_hybrid_test_params <- function() list(alpha = .6, theta1 = .85,
  theta0 = .15, pi = .05, eps = .20, lambda_g = 2, lambda_d = 1.5)

test_that("hybrid tenure emission nests marginal contamination at zero local mass", {
  lambda <- c(.32, .30, .18, .15, .13)
  cases <- list(
    list(g = matrix(c(2, 2.25), nrow = 1), s = matrix(TRUE, 1, 2), d = 0:1),
    list(g = matrix(c(2, 4), nrow = 1), s = matrix(TRUE, 1, 2), d = 0:1),
    list(g = matrix(c(2, 2.25, 2.5), nrow = 1),
         s = matrix(TRUE, 1, 3), d = 0:2))
  for (x in cases) {
    old <- log_emission_spell_g(x$g, x$s, as.integer(x$d), lambda, .25)
    new <- log_emission_spell_g_local_gross(x$g, x$s, as.integer(x$d),
      lambda, eps_local = 0, eps_gross = .25)
    expect_equal(new$loglik, old$loglik, tolerance = 1e-10)
  }
})

test_that("local tenure branch rewards bounded month discrepancies", {
  lambda <- c(.32, .30, .18, .15, .13)
  s <- matrix(TRUE, 1, 2); d <- 0:1
  near <- matrix(c(2, 2.25 + 1/12), nrow = 1)
  far <- matrix(c(2, 2.25 + 18/12), nrow = 1)
  near_hybrid <- log_emission_spell_g_local_gross(near, s, d, lambda,
    eps_local = .12, eps_gross = .13)$loglik
  near_gross <- log_emission_spell_g_local_gross(near, s, d, lambda,
    eps_local = 0, eps_gross = .25)$loglik
  far_hybrid <- log_emission_spell_g_local_gross(far, s, d, lambda,
    eps_local = .12, eps_gross = .13)$loglik
  expect_gt(near_hybrid, near_gross)
  expect_true(is.finite(far_hybrid))
})

test_that("hybrid transform preserves total tenure contamination", {
  p <- list(alpha = .48, pi = .025, eps = .25, eps_local = .10,
    eps_d = .16, lambda_g = c(.32, .30, .18, .15, .13),
    lambda_d = c(.44, .21, .09, .10, .11))
  z <- .piecewise_hybrid_pack(p)
  out <- .piecewise_hybrid_unpack(z)
  expect_equal(out$eps_local, .10, tolerance = 1e-10)
  expect_equal(out$eps_gross, .15, tolerance = 1e-10)
  expect_equal(out$eps_local + out$eps_gross, out$eps, tolerance = 1e-12)
})

test_that("hybrid tenure E-step returns a valid observed likelihood", {
  df <- make_hybrid_test_data()
  p <- make_hybrid_test_params()
  p$eps_local <- .08
  p$eps_gross <- p$eps - p$eps_local
  p$tenure_contamination_model <- "local_gross"
  out <- e_step_eps(df, p, suff_stats = FALSE)
  expect_true(is.finite(out$loglik))
  expect_equal(rowSums(out$gamma), rep(1, nrow(df)), tolerance = 1e-10)
})
