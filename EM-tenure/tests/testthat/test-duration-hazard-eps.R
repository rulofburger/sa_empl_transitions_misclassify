make_duration_eps_data <- function(n = 100L, seed = 42L) {
  set.seed(seed)
  s <- sample(0:1, n, replace = TRUE)
  g1 <- pmax(.26, rexp(n, 2))
  data.frame(y1 = s, y2 = s, y3 = s,
    tenure1 = g1, tenure2 = g1 + .25, tenure3 = g1 + .5,
    timegap_cat1 = sample(1:7, n, TRUE),
    timegap_cat2 = sample(1:7, n, TRUE),
    timegap_cat3 = sample(1:7, n, TRUE), weight = 1)
}

make_duration_eps_params <- function() list(alpha=.6,theta0=.15,theta1=.85,
  pi=.05,eps=.2,lambda_g=2,lambda_d=1.5)

test_that("beta zero nests the exponential duration distribution", {
  x <- c(.1, .5, 2, 8)
  lambda <- .23
  expect_equal(.log_duration_density(x, lambda, 0),
               log(lambda) - lambda * x, tolerance = 1e-12)
  for (k in 1:7) {
    expect_equal(log_emission_interval_d(k, lambda, 0),
                 log_emission_interval_d(k, lambda), tolerance = 1e-12)
  }
})

test_that("beta zero history prior equals the constant linked prior", {
  df <- make_duration_eps_data(n = 25L, seed = 812L)
  h <- latent_histories()
  alpha <- .57; lambda_g <- .18; lambda_d <- .22
  theta1 <- exp(-lambda_g * .QUARTER_YEARS)
  theta0 <- 1 - exp(-lambda_d * .QUARTER_YEARS)
  lp <- .log_duration_history_prior_eps(
    h, alpha, lapply(1:3, function(t) df[[paste0("y", t)]]),
    lapply(1:3, function(t) df[[paste0("tenure", t)]]),
    lapply(1:3, function(t) df[[paste0("timegap_cat", t)]]),
    lambda_g, 0, lambda_d, 0)
  target <- log(prior_over_histories(h, theta1, theta0, alpha))
  expect_equal(lp, matrix(target, nrow(df), 8, byrow = TRUE),
               tolerance = 1e-11)
})

test_that("duration-dependent E-step nests linked epsilon likelihood", {
  df <- make_duration_eps_data(n = 80L, seed = 925L)
  p <- make_duration_eps_params()
  p$lambda_g <- ctmc_lambda_from_persistence(p$theta1)
  p$lambda_d <- ctmc_lambda_from_transition(p$theta0)
  old <- e_step_eps(df, p, suff_stats = FALSE)
  p$beta_g <- 0; p$beta_d <- 0
  nested <- e_step_eps(df, p, suff_stats = FALSE)
  expect_equal(nested$loglik, old$loglik, tolerance = 1e-9)
  expect_equal(nested$gamma, old$gamma, tolerance = 1e-10)
})

test_that("negative beta produces declining quarterly transition rates", {
  d <- c(0, .5, 2, 5)
  p <- .duration_transition_probability(d, .2, -.5)
  expect_true(all(diff(p) < 0))
  expect_true(all(p > 0 & p < 1))
})

test_that("marginal missing-clock transition risk is finite and nests constant hazard", {
  lambda <- .27
  expect_equal(.duration_marginal_transition_probability(lambda, 0),
    .duration_transition_probability(0, lambda, 0), tolerance=1e-12)
  p <- .duration_marginal_transition_probability(lambda, -.75)
  expect_true(is.finite(p) && p > 0 && p < 1)
})

test_that("category-integrated transition risks are finite for heavy tails", {
  p <- .duration_category_transition_probability(1:7, .4, -.9)
  expect_true(all(is.finite(p)))
  expect_true(all(p > 0 & p < 1))
  expect_true(all(diff(p) < 0))
})
