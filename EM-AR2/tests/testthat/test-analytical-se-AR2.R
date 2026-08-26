test_that("collapse_ar2_cells preserves likelihood inputs", {
  set.seed(19)
  df <- data.frame(y1 = sample(0:1, 100, TRUE), y2 = sample(0:1, 100, TRUE),
                   y3 = sample(0:1, 100, TRUE), y4 = sample(0:1, 100, TRUE),
                   weight = runif(100, 0.5, 2))
  cells <- collapse_ar2_cells(df)
  expect_equal(sum(cells$count), 100)
  expect_equal(sum(cells$weight), sum(df$weight), tolerance = 1e-12)
  expect_equal(sum(cells$weight_sq), sum(df$weight^2), tolerance = 1e-12)
})

test_that("analytical_se_ar2 returns finite delta-method standard errors", {
  set.seed(31)
  histories <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(histories, 0.10, 0.15, 0.08, 0.12)
  latent <- histories[sample.int(16L, 5000L, TRUE, prior), ]
  observed <- abs(latent - matrix(runif(20000L) < 0.03, 5000L, 4L))
  df <- data.frame(y1 = observed[, 1], y2 = observed[, 2],
                   y3 = observed[, 3], y4 = observed[, 4], weight = 1)
  fit <- em_fit_ar2(df, max_iter = 1000L, tol = 1e-10,
                    param_tol = 1e-8, verbose = 0L)
  inf <- analytical_se_ar2(df, fit, "symmetric")
  expect_true(all(is.finite(inf$summary$se)))
  expect_true(all(inf$summary$se >= 0))
  expect_true(all(c("p_00", "p_10", "p_01", "p_11", "employment_rate") %in%
                  inf$summary$quantity))
  expect_equal(inf$n_obs, nrow(df))
})
