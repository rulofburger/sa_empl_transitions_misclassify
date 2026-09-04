test_that("zero job-change probability exactly nests the existing emission", {
  g <- matrix(c(2, 2.25, 2.50,
                2, 0.10, 0.35,
                1, 1.25, 1.50), nrow = 3, byrow = TRUE)
  s <- matrix(c(TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE,
                TRUE, FALSE, TRUE), nrow = 3, byrow = TRUE)
  old <- log_emission_spell_g(g, s, 0:2, rep(.2, 5), .25)
  new <- log_emission_spell_g_job_change(
    g, s, 0:2, rep(.2, 5), .25, job_change_prob = 0)
  expect_identical(new$loglik, old$loglik)
  expect_equal(new$job_changes, rep(0, 3))
  expect_equal(new$job_change_opportunities, rep(2, 3))
})

test_that("a persistent short-tenure reset receives high job-change weight", {
  g <- matrix(c(2, .10, .35), nrow = 1)
  s <- matrix(TRUE, nrow = 1, ncol = 3)
  old <- log_emission_spell_g(g, s, 0:2, rep(.2, 5), .25)
  new <- log_emission_spell_g_job_change(
    g, s, 0:2, rep(.2, 5), .25, job_change_prob = .05)
  expect_true(is.finite(new$loglik))
  expect_gt(new$loglik, old$loglik)
  expect_gt(new$job_changes, .90)
})

test_that("unchanged tenure clocks do not spuriously favor a reset", {
  g <- matrix(c(2, 2.25, 2.50), nrow = 1)
  s <- matrix(TRUE, nrow = 1, ncol = 3)
  new <- log_emission_spell_g_job_change(
    g, s, 0:2, rep(.2, 5), .25, job_change_prob = .05)
  expect_true(is.finite(new$loglik))
  expect_lt(new$job_changes, .01)
})

test_that("the complete E-step is exactly nested at q zero", {
  df <- data.frame(
    y1 = c(1L, 1L, 0L, 1L), y2 = c(1L, 1L, 1L, 0L),
    y3 = c(1L, 0L, 1L, 0L),
    tenure1 = c(2, 1, NA, .5), tenure2 = c(2.25, .1, .2, NA),
    tenure3 = c(2.5, NA, .45, NA),
    timegap_cat1 = c(NA, NA, 4L, NA),
    timegap_cat2 = c(NA, NA, NA, 2L),
    timegap_cat3 = c(NA, 1L, NA, 3L), weight = rep(1, 4))
  p <- list(alpha = .6, theta0 = .08, theta1 = .95, pi = .02,
    eps = .25, eps_d = .15, lambda_g = rep(.2, 5), beta_g = 0,
    lambda_d = rep(.3, 5), beta_d = 0,
    timegap_contamination_model = "joint_marginal")
  old <- e_step_eps(df, p, suff_stats = FALSE)
  p$job_change_prob <- 0
  new <- e_step_eps(df, p, suff_stats = FALSE)
  expect_identical(new$loglik, old$loglik)
  expect_identical(new$gamma, old$gamma)
})

test_that("job-change parameter validation is strict", {
  g <- matrix(c(1, 1.25), nrow = 1)
  s <- matrix(TRUE, nrow = 1, ncol = 2)
  expect_error(log_emission_spell_g_job_change(
    g, s, 0:1, rep(.2, 5), .2, -0.01), "must be in")
  expect_error(log_emission_spell_g_job_change(
    g, s, 0:1, rep(.2, 5), .2, 1), "must be in")
})

test_that("the job-change probability is recovered in an emission simulation", {
  set.seed(17031)
  n <- 800L
  q_true <- .06
  eps <- .15
  lambda <- rep(.2, 5)
  g <- matrix(NA_real_, n, 3L)
  g[, 1L] <- rexp(n, .2)
  for (j in 2:3) {
    reset <- runif(n) < q_true
    g[, j] <- ifelse(reset, runif(n, 0, .25), g[, j - 1L] + .25)
  }
  contaminated <- matrix(runif(3L * n) < eps, n, 3L)
  gross <- matrix(rexp(3L * n, .2), n, 3L)
  g[contaminated] <- gross[contaminated]
  observed <- matrix(TRUE, n, 3L)
  criterion <- function(q) -sum(log_emission_spell_g_job_change(
    g, observed, 0:2, lambda, eps, q)$loglik)
  q_hat <- optimize(criterion, c(1e-5, .20), tol = 1e-5)$minimum
  expect_equal(q_hat, q_true, tolerance = .02)
})

test_that("the full-likelihood simulator returns valid observed data", {
  p <- list(alpha=.55, pi=.03, eps=.20, eps_d=.12,
    job_change_prob=.04, lambda_g=rep(.20, 5), lambda_d=rep(.30, 5))
  d <- simulate_eps_piecewise_job_change(2000, p, seed=991)
  expect_silent(validate_df_eps(d))
  expect_true(all(is.na(d$tenure1[d$y1 == 0L])))
  expect_true(all(is.na(d$timegap_cat1[d$y1 == 1L])))
  opportunities <- (d$h1 == 1L & d$h2 == 1L) +
    (d$h2 == 1L & d$h3 == 1L)
  expect_equal(sum(d$reset12) + sum(d$reset23),
    p$job_change_prob * sum(opportunities), tolerance=15)
})

test_that("fixed-q optimisation retains the requested reset probability", {
  p <- list(alpha=.55, pi=.03, eps=.20, eps_d=.12,
    job_change_prob=.04, lambda_g=rep(.20, 5), beta_g=0,
    lambda_d=rep(.30, 5), beta_d=0,
    timegap_contamination_model="joint_marginal")
  d <- simulate_eps_piecewise_job_change(300, p, seed=992)
  fit <- suppressWarnings(fit_eps_piecewise_job_change(
    d, p, q_fixed=.025, maxit=1))
  expect_identical(fit$params$job_change_prob, .025)
  expect_false("job_change" %in% names(fit$par_unconstrained))
})
