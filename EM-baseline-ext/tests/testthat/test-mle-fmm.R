test_that("FMM cell probabilities are valid and invariant to label switching", {
  p <- init_params_fmm("symmetric", stationary = FALSE)
  probs <- fmm_cell_probabilities(p, "symmetric")
  expect_equal(sum(probs), 1, tolerance = 1e-12)
  expect_true(all(probs > 0))

  swapped <- p
  for (stem in c("theta0", "theta1", "alpha")) {
    swapped[[paste0(stem, "_A")]] <- p[[paste0(stem, "_B")]]
    swapped[[paste0(stem, "_B")]] <- p[[paste0(stem, "_A")]]
  }
  swapped$phi <- 1 - p$phi
  expect_equal(fmm_cell_probabilities(swapped, "symmetric"), probs,
               tolerance = 1e-12)
})

test_that("parameter count exposes unrestricted symmetric FMM underidentification", {
  expect_equal(fmm_parameter_count("none", TRUE), 5L)
  expect_equal(fmm_parameter_count("symmetric", TRUE), 6L)
  expect_equal(fmm_parameter_count("none", FALSE), 7L)
  expect_equal(fmm_parameter_count("symmetric", FALSE), 8L)
  expect_gt(fmm_parameter_count("symmetric", FALSE), 7L)
})
