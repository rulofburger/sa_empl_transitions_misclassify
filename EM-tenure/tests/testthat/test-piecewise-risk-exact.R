test_that("exact piecewise risks agree with independent numerical integrals", {
  for (h in list(rep(.2, 5), c(.25,.28,.17,.15,.13), c(.45,.21,.09,.1,.06))) {
    for (delta in c(.25, .5)) for (category in 0:7) {
      interval <- if (category == 0) c(0, Inf) else .timegap_interval(category)
      lower <- interval[1]; upper <- interval[2]
      knots <- .DURATION_PIECEWISE_KNOTS
      cuts <- sort(unique(c(lower, upper, knots, knots - delta)))
      cuts <- cuts[cuts >= lower & cuts <= upper]
      ha <- .piecewise_duration_cumhaz(lower, h)
      integrand <- function(x) {
        hx <- .piecewise_duration_cumhaz(x, h)
        .piecewise_duration_hazard(x, h) * exp(-(hx - ha)) *
          (-expm1(-(.piecewise_duration_cumhaz(x + delta, h) - hx)))
      }
      numerator <- sum(vapply(seq_len(length(cuts)-1L), function(j)
        integrate(integrand, cuts[j], cuts[j+1], rel.tol = 1e-11,
          abs.tol = 1e-13)$value, numeric(1)))
      denom <- if (is.infinite(upper)) 1 else
        -expm1(-(.piecewise_duration_cumhaz(upper, h) - ha))
      expect_equal(.piecewise_mean_transition_risk_exact(h, lower, upper, delta),
        numerator / denom, tolerance = 1e-10)
    }
  }
})

test_that("constant-hazard and zero-gap risks have their exact limits", {
  expect_equal(.piecewise_mean_transition_risk_exact(rep(.2,5), 100, Inf),
    -expm1(-.2*.25), tolerance = 1e-12)
  expect_equal(.piecewise_mean_transition_risk_exact(rep(.2,5), delta = 0),
    0, tolerance = 1e-12)
})
