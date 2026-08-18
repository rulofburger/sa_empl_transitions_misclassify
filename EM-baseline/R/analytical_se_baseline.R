# ==============================================================================
# EM-baseline: survey-weighted sandwich covariance and delta-method SEs
# ==============================================================================

.baseline_numeric_jacobian <- function(fn, x, rel_step = 1e-5) {
  f0 <- fn(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(x),
              dimnames = list(names(f0), names(x)))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x; xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
    J[, j] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  J
}

.baseline_reported_quantities <- function(eta, model_type, stationary) {
  p <- unpack_baseline_eta(eta, model_type, stationary)
  implied <- implied_baseline(p, model_type)
  out <- c(
    alpha = p$alpha,
    theta0 = p$theta0,
    theta1 = p$theta1
  )
  if (model_type == "symmetric") out <- c(out, pi = p$pi)
  if (model_type == "asymmetric") out <- c(out, pi0 = p$pi0, pi1 = p$pi1)
  c(out,
    entry_rate = implied$entry_rate,
    exit_rate = implied$exit_rate,
    employment_rate = implied$employment_rate)
}

analytical_se_baseline <- function(df, fit, finite_sample = TRUE) {
  stopifnot(is.list(fit), identical(fit$estimator, "direct_eight_cell_mle"))
  cells <- collapse_baseline_cells(df)
  if (!identical(baseline_sample_signature(cells), fit$sample$signature))
    stop("analytical_se_baseline: fit and data sample signatures differ")

  eta <- fit$eta
  model_type <- fit$model_type
  stationary <- fit$stationary
  log_cell <- function(z) {
    p <- unpack_baseline_eta(z, model_type, stationary)
    setNames(log(pmax(baseline_cell_probabilities(p, model_type), 1e-300)),
             apply(cells$histories, 1L, paste0, collapse = ""))
  }
  cell_scores <- .baseline_numeric_jacobian(log_cell, eta)

  fn <- function(z) .baseline_average_nll(z, cells, model_type, stationary)
  bread <- optimHess(eta, fn)
  bread <- (bread + t(bread)) / 2
  meat <- matrix(0, nrow = length(eta), ncol = length(eta),
                 dimnames = list(names(eta), names(eta)))
  for (j in seq_len(8L)) {
    s <- cell_scores[j, ]
    meat <- meat + (cells$weight_sq[j] / cells$weight_sum^2) * tcrossprod(s)
  }
  bread_inv <- solve(bread)
  vcov_eta <- bread_inv %*% meat %*% bread_inv
  K <- length(eta)
  if (finite_sample && cells$n > K) vcov_eta <- vcov_eta * cells$n / (cells$n - K)
  dimnames(vcov_eta) <- list(names(eta), names(eta))

  qfun <- function(z) .baseline_reported_quantities(z, model_type, stationary)
  estimates <- qfun(eta)
  delta <- .baseline_numeric_jacobian(qfun, eta)
  vcov_q <- delta %*% vcov_eta %*% t(delta)
  ses <- sqrt(pmax(diag(vcov_q), 0))
  eig <- eigen(bread, symmetric = TRUE, only.values = TRUE)$values

  list(
    summary = data.frame(quantity = names(estimates), estimate = unname(estimates),
                         se = unname(ses), row.names = NULL),
    vcov_eta = vcov_eta,
    vcov_quantities = vcov_q,
    bread = bread,
    meat = meat,
    delta_jacobian = delta,
    diagnostics = list(
      information_min_eigenvalue = min(eig),
      information_condition = max(eig) / min(eig),
      max_abs_score = fit$diagnostics$max_abs_score
    ),
    label = fit$label %||% NA_character_,
    model_type = model_type,
    stationary = stationary,
    sample_signature = cells$signature %||% baseline_sample_signature(cells),
    method = "individual_survey_weighted_sandwich_delta"
  )
}
