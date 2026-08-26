# Exact baseline MLE and analytical inference

.make_mle_panel <- function(n = 200, theta0 = 0.10, theta1 = 0.90,
                            pi = 0.05, seed = 99) {
  set.seed(seed)
  alpha <- theta0 / (theta0 + 1 - theta1)
  h1 <- rbinom(n, 1, alpha)
  h2 <- ifelse(h1 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  h3 <- ifelse(h2 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  flip <- function(h) ifelse(rbinom(length(h), 1, pi) == 1, 1 - h, h)
  data.frame(y1 = flip(h1), y2 = flip(h2), y3 = flip(h3), weight = rep(1, n))
}

test_that("eight-cell likelihood matches the row-level likelihood", {
  df <- .make_mle_panel(n = 500, theta0 = 0.12, theta1 = 0.88, pi = 0.03, seed = 81)
  p <- list(alpha = 0.5, theta0 = 0.11, theta1 = 0.89, pi = 0.04)
  cells <- collapse_baseline_cells(df)
  ll_cell <- sum(cells$weight * log(baseline_cell_probabilities(p, "symmetric")))
  ll_row <- compute_observed_loglik(df, p, "symmetric")
  expect_equal(ll_cell, ll_row, tolerance = 1e-8)
})

test_that("exact MLE passes convergence and nesting checks", {
  df <- .make_mle_panel(n = 1500, theta0 = 0.10, theta1 = 0.90, pi = 0.03, seed = 92)
  none <- fit_baseline_mle(df, "none", TRUE, verbose = 0L)
  sym_start <- none$params; sym_start$pi <- 0.02
  sym <- fit_baseline_mle(df, "symmetric", TRUE,
                          starts = list(sym_start, init_params("symmetric", TRUE)),
                          verbose = 0L)
  expect_true(none$converged)
  expect_true(sym$converged)
  expect_gte(sym$loglik, none$loglik - 1e-5)
  expect_lt(sym$diagnostics$max_abs_score, 1e-6)
})

test_that("analytical covariance is finite and delta quantities agree", {
  df <- .make_mle_panel(n = 1800, theta0 = 0.12, theta1 = 0.91, pi = 0.025, seed = 103)
  fit <- fit_baseline_mle(df, "symmetric", FALSE,
                          starts = list(init_params("symmetric", FALSE)),
                          verbose = 0L)
  inf <- analytical_se_baseline(df, fit)
  expect_true(all(is.finite(inf$summary$se)))
  expect_true(all(inf$summary$se > 0))
  expect_equal(
    inf$summary$estimate[inf$summary$quantity == "exit_rate"],
    1 - fit$params$theta1,
    tolerance = 1e-10
  )
  expect_gt(inf$diagnostics$information_min_eigenvalue, 0)
})

test_that("analytical inference rejects a mismatched sample", {
  df <- .make_mle_panel(n = 600, seed = 104)
  fit <- fit_baseline_mle(df, "none", FALSE, verbose = 0L)
  expect_error(analytical_se_baseline(df[-1L, ], fit), "signatures differ")
})

test_that("exact MLE records explicit and attributed panel provenance", {
  df <- .make_mle_panel(n = 600, seed = 105)
  fit_explicit <- fit_baseline_mle(
    df, "none", FALSE, verbose = 0L, source_panel = "df_qlfs_B.rds"
  )
  expect_identical(fit_explicit$sample$source_panel, "df_qlfs_B.rds")

  attr(df, "source_panel") <- "df_qlfs_C.rds"
  fit_attributed <- fit_baseline_mle(df, "none", FALSE, verbose = 0L)
  expect_identical(fit_attributed$sample$source_panel, "df_qlfs_C.rds")
})
