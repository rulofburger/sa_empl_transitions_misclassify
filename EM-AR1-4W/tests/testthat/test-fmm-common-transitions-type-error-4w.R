test_that("common-transition type-error probabilities are valid and label invariant", {
  p <- list(theta0=.04, theta1=.95, alpha=.48,
            phi=.7, pi_A=.01, pi_B=.07)
  swapped <- list(theta0=p$theta0, theta1=p$theta1, alpha=p$alpha,
                  phi=1-p$phi, pi_A=p$pi_B, pi_B=p$pi_A)
  probs <- fmm_ctte_4w_cell_probabilities(p)
  expect_equal(sum(probs), 1, tolerance=1e-12)
  expect_true(all(probs > 0))
  expect_equal(probs, fmm_ctte_4w_cell_probabilities(swapped), tolerance=1e-12)
  expect_equal(pack_fmm_ctte_4w(resolve_fmm_ctte_labels(p)),
               pack_fmm_ctte_4w(resolve_fmm_ctte_labels(swapped)), tolerance=1e-12)
})

test_that("common-transition type-error model has full generic local rank", {
  eta <- pack_fmm_ctte_4w(list(theta0=.04, theta1=.95,
    alpha=.48, phi=.72, pi_A=.008, pi_B=.065))
  J <- .ar1_4w_jacobian(function(z) setNames(
    fmm_ctte_4w_cell_probabilities(unpack_fmm_ctte_4w(z))[1:15],
    paste0("cell", 1:15)), eta)
  expect_equal(qr(J, tol=1e-8)$rank, length(eta))
})
