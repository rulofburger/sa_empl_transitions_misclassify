# Exact E[1-exp{-(H(D+delta)-H(D))} | lower <= D < upper]
# for a piecewise-constant duration hazard. Split at both the hazard knots
# and their delta-shifted counterparts, then integrate the exponential tails.
.piecewise_mean_transition_risk_exact <- function(hazards, lower = 0,
    upper = Inf, delta = .QUARTER_YEARS) {
  knots <- .DURATION_PIECEWISE_KNOTS
  if (length(hazards) != length(knots)-1L || any(!is.finite(hazards)) ||
      any(hazards <= 0) || lower < 0 || lower >= upper || delta < 0)
    stop("Invalid piecewise risk inputs")
  cuts <- sort(unique(c(lower, upper,
    knots[is.finite(knots) & knots > lower & knots < upper],
    (knots - delta)[is.finite(knots) & knots - delta > lower &
      knots - delta < upper])))
  origin <- cuts[-length(cuts)]
  end <- cuts[-1L]
  current <- .piecewise_duration_hazard(origin, hazards)
  future <- .piecewise_duration_hazard(origin + delta, hazards)
  ha <- .piecewise_duration_cumhaz(lower, hazards)
  # Scale by S(lower), avoiding underflow in old-duration tail categories.
  stay_mass <- sum(current / future *
    exp(-(.piecewise_duration_cumhaz(origin + delta, hazards) - ha)) *
    (-expm1(-future * (end - origin))))
  category_mass <- if (is.infinite(upper)) 1 else
    -expm1(-(.piecewise_duration_cumhaz(upper, hazards) - ha))
  pmin(1, pmax(0, 1 - stay_mass / category_mass))
}

.piecewise_category_transition_risk_exact <- function(category, hazards) {
  lookup <- vapply(1:7, function(k) {
    interval <- .timegap_interval(k)
    .piecewise_mean_transition_risk_exact(hazards, interval[1], interval[2])
  }, numeric(1))
  lookup[category]
}
