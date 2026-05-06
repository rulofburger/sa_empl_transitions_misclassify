# ==============================================================================
# EM-AR2 tests: inference.R
# ==============================================================================

# Helper: known params for no-ME model (pi=0)
.known_params <- function(pi = 0) {
  list(theta0=0.10, theta01=0.15, theta1=0.08, theta10=0.12, pi=pi)
}

test_that("model_cell_probs_ar2 returns 16 rows summing to 1", {
  probs <- model_cell_probs_ar2(.known_params())
  expect_equal(nrow(probs), 16L)
  expect_equal(sum(probs$model_prob), 1, tolerance = 1e-10)
})

test_that("model_cell_probs_ar2 all probabilities are non-negative", {
  probs <- model_cell_probs_ar2(.known_params())
  expect_true(all(probs$model_prob >= 0))
})

test_that("model_cell_probs_ar2 with pi=0 matches marginalised latent prior", {
  # With pi=0, P(y=h) so P(y1,y2,y3,y4) = P(h1,h2,h3,h4)
  params <- .known_params(pi = 0)
  hmat   <- latent_histories_ar2()
  prior  <- prior_over_histories_ar2(hmat,
    params$theta0, params$theta01, params$theta1, params$theta10)

  probs  <- model_cell_probs_ar2(params)

  # With pi=0, each cell prob should equal the prior of the matching history
  for (i in seq_len(nrow(probs))) {
    y   <- as.integer(probs[i, 1:4])
    idx <- which(hmat[, 1] == y[1] & hmat[, 2] == y[2] &
                 hmat[, 3] == y[3] & hmat[, 4] == y[4])
    expect_equal(probs$model_prob[i], unname(prior[idx]), tolerance = 1e-10)
  }
})

test_that("implied_transitions_ar2 returns named vector of length 4", {
  trans <- implied_transitions_ar2(.known_params())
  expect_length(trans, 4L)
  expected_names <- c("p_emp_w4_from0", "p_ever_emp_from0",
                      "p_nonemp_w4_from1", "p_ever_nonemp_from1")
  expect_named(trans, expected_names)
})

test_that("implied_transitions_ar2 all values are in [0,1]", {
  trans <- implied_transitions_ar2(.known_params())
  expect_true(all(trans >= 0))
  expect_true(all(trans <= 1))
})

test_that("implied_transitions_ar2 ever-transition >= wave-4 transition", {
  # P(ever employed from 0) >= P(employed at wave 4 from 0)
  trans <- implied_transitions_ar2(.known_params())
  expect_gte(trans["p_ever_emp_from0"],   trans["p_emp_w4_from0"])
  expect_gte(trans["p_ever_nonemp_from1"], trans["p_nonemp_w4_from1"])
})

test_that("goodness_of_fit_ar2 returns correct structure", {
  # Small synthetic data
  set.seed(5)
  N  <- 100
  df <- data.frame(
    y1 = sample(0:1, N, replace=TRUE),
    y2 = sample(0:1, N, replace=TRUE),
    y3 = sample(0:1, N, replace=TRUE),
    y4 = sample(0:1, N, replace=TRUE),
    weight = rep(1, N)
  )
  gof <- goodness_of_fit_ar2(df, .known_params())

  expect_named(gof, c("table", "ssr"))
  expect_equal(nrow(gof$table), 16L)
  expect_true(is.numeric(gof$ssr))
  expect_gte(gof$ssr, 0)
})

test_that("goodness_of_fit_ar2 empirical probs sum to 1", {
  set.seed(6)
  N  <- 80
  df <- data.frame(
    y1 = sample(0:1, N, replace=TRUE), y2 = sample(0:1, N, replace=TRUE),
    y3 = sample(0:1, N, replace=TRUE), y4 = sample(0:1, N, replace=TRUE),
    weight = rep(1, N)
  )
  gof <- goodness_of_fit_ar2(df, .known_params())
  expect_equal(sum(gof$table$empirical), 1, tolerance = 1e-10)
})

test_that("goodness_of_fit_ar2 SSR is near zero when model fits well", {
  # Generate data from the exact model (large N for good empirical approximation)
  set.seed(42)
  theta0 <- 0.10; theta01 <- 0.15; theta1 <- 0.08; theta10 <- 0.12
  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat, theta0, theta01, theta1, theta10)

  N     <- 5000
  h_idx <- sample.int(16, size = N, replace = TRUE, prob = prior)
  h     <- hmat[h_idx, ]
  df    <- data.frame(y1=h[,1], y2=h[,2], y3=h[,3], y4=h[,4], weight=rep(1,N))

  params <- list(theta0=theta0, theta01=theta01, theta1=theta1, theta10=theta10, pi=0)
  gof    <- goodness_of_fit_ar2(df, params)

  # With N=5000 and pi=0, SSR should be small (< 0.01)
  expect_lt(gof$ssr, 0.01)
})

test_that("implied_transitions_ar2 errors on degenerate prior (no mass on h1=0)", {
  # Supply alpha with zero mass on all h1=0 states: "00" (h1=0,h2=0) and "01" (h1=0,h2=1)
  degenerate_params <- list(
    theta0  = 0.10, theta01 = 0.15,
    theta1  = 0.08, theta10 = 0.12,
    pi      = 0,
    alpha   = c(`00` = 0, `01` = 0, `10` = 0.5, `11` = 0.5)
  )
  expect_error(
    implied_transitions_ar2(degenerate_params),
    regexp = "no prior mass on h1=0"
  )
})

test_that("implied_transitions_ar2 errors on degenerate prior (no mass on h1=1)", {
  # Supply alpha with zero mass on all h1=1 states: "10" (h1=1,h2=0) and "11" (h1=1,h2=1)
  degenerate_params <- list(
    theta0  = 0.10, theta01 = 0.15,
    theta1  = 0.08, theta10 = 0.12,
    pi      = 0,
    alpha   = c(`00` = 0.5, `01` = 0.5, `10` = 0, `11` = 0)
  )
  expect_error(
    implied_transitions_ar2(degenerate_params),
    regexp = "no prior mass on h1=1"
  )
})
