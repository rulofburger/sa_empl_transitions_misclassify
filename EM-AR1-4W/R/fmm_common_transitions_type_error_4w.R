# ==============================================================================
# Two-type four-wave AR(1) FMM with common transitions and type-specific error
# ==============================================================================

.fmm_ctte_names <- c("theta0", "theta1", "alpha", "phi", "pi_A", "pi_B")

pack_fmm_ctte_4w <- function(params) {
  clamp <- function(x, eps = 1e-8) pmin(pmax(x, eps), 1 - eps)
  out <- setNames(numeric(length(.fmm_ctte_names)), .fmm_ctte_names)
  for (nm in c("theta0", "theta1", "alpha", "phi"))
    out[nm] <- qlogis(clamp(params[[nm]]))
  for (nm in c("pi_A", "pi_B"))
    out[nm] <- qlogis(clamp(params[[nm]] / .5))
  out
}

unpack_fmm_ctte_4w <- function(eta) {
  eta <- setNames(as.numeric(eta), .fmm_ctte_names)
  p <- as.list(plogis(eta))
  p$pi_A <- .5 * plogis(eta["pi_A"])
  p$pi_B <- .5 * plogis(eta["pi_B"])
  p
}

resolve_fmm_ctte_labels <- function(params) {
  if (params$pi_A <= params$pi_B) return(params)
  out <- params
  out$pi_A <- params$pi_B
  out$pi_B <- params$pi_A
  out$phi <- 1 - params$phi
  out
}

fmm_ctte_4w_cell_probabilities <- function(params) {
  common <- list(theta0 = params$theta0, theta1 = params$theta1)
  pA <- ar1_4w_cell_probabilities(c(common,
    list(alpha = params$alpha, pi = params$pi_A)), "symmetric")
  pB <- ar1_4w_cell_probabilities(c(common,
    list(alpha = params$alpha, pi = params$pi_B)), "symmetric")
  params$phi * pA + (1 - params$phi) * pB
}

initial_fmm_ctte_4w <- function() {
  list(theta0 = .045, theta1 = .955, alpha = .48,
       phi = .75, pi_A = .006, pi_B = .045)
}

.boundary_fmm_ctte_starts <- function() {
  lapply(c(1e-8, 1e-6, 1e-4), function(low_error)
    list(theta0 = .026, theta1 = .975, alpha = .479,
         phi = .80, pi_A = low_error, pi_B = .175))
}

.random_fmm_ctte_start <- function() {
  resolve_fmm_ctte_labels(list(
    theta0 = runif(1, .005, .18), theta1 = runif(1, .82, .995),
    alpha = runif(1, .25, .75),
    phi = runif(1, .08, .92), pi_A = runif(1, .0005, .04),
    pi_B = runif(1, .015, .18)))
}

fit_fmm_ctte_4w <- function(df, starts = NULL, random_starts = 120L,
                            seed = 20260828L, refine_top = 12L,
                            maxit = 5000L, reltol = 1e-13, verbose = 1L) {
  cells <- collapse_ar1_4w_cells(df)
  if (is.null(starts)) starts <- c(list(initial_fmm_ctte_4w()),
                                   .boundary_fmm_ctte_starts())
  set.seed(seed)
  if (random_starts > 0L) for (ii in seq_len(random_starts))
    starts[[length(starts) + 1L]] <- .random_fmm_ctte_start()
  fn <- function(z) {
    probs <- pmax(fmm_ctte_4w_cell_probabilities(unpack_fmm_ctte_4w(z)), 1e-300)
    -sum(cells$weight * log(probs)) / sum(cells$weight)
  }
  # A moderate first pass screens broadly; only the best distinct basins receive
  # the stringent likelihood polish used for inference.
  screened <- lapply(starts, function(start) tryCatch(
    optim(pack_fmm_ctte_4w(start), fn, method = "BFGS",
          control = list(maxit = 700L, reltol = 1e-9)), error = function(e) NULL))
  screened <- Filter(function(x) !is.null(x) && is.finite(x$value), screened)
  if (!length(screened)) stop("all common-transition FMM starts failed")
  screened <- screened[order(vapply(screened, `[[`, numeric(1), "value"))]
  screened <- screened[seq_len(min(length(screened), refine_top))]
  candidates <- lapply(screened, function(first) tryCatch(
    optim(first$par, fn, method = "BFGS",
          control = list(maxit = maxit, reltol = reltol)), error = function(e) NULL))
  candidates <- Filter(function(x) !is.null(x) && is.finite(x$value), candidates)
  if (!length(candidates)) stop("all common-transition FMM refinements failed")
  best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
  params <- resolve_fmm_ctte_labels(unpack_fmm_ctte_4w(best$par))
  eta <- pack_fmm_ctte_4w(params)
  # Re-polish after label normalization so the reported score is evaluated at
  # the exact parameter vector stored in the result.
  opt <- optim(eta, fn, method = "BFGS",
               control = list(maxit = maxit, reltol = reltol))
  params <- resolve_fmm_ctte_labels(unpack_fmm_ctte_4w(opt$par))
  eta <- pack_fmm_ctte_4w(params)
  probs <- fmm_ctte_4w_cell_probabilities(params)
  loglik <- sum(cells$weight * log(pmax(probs, 1e-300)))
  score <- .ar1_4w_gradient(fn, eta)
  H <- optimHess(eta, fn); H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  prob_jac <- .ar1_4w_jacobian(function(z) setNames(
    fmm_ctte_4w_cell_probabilities(unpack_fmm_ctte_4w(z))[1:15],
    paste0("cell", 1:15)), eta)
  rank <- qr(prob_jac, tol = 1e-8)$rank
  K <- length(eta)
  candidate_table <- do.call(rbind, lapply(candidates, function(x) {
    p <- resolve_fmm_ctte_labels(unpack_fmm_ctte_4w(x$par))
    data.frame(loglik = -x$value * sum(cells$weight),
      entry = p$theta0, exit = 1 - p$theta1, alpha = p$alpha,
      share_A = p$phi, pi_A = p$pi_A, pi_B = p$pi_B)
  }))
  candidate_table <- candidate_table[order(candidate_table$loglik, decreasing = TRUE), ]
  # A type-specific error probability may be estimated on the zero boundary.
  # In that case logit-scale curvature vanishes even though the optimizer score
  # and probability Jacobian are satisfactory; retain the curvature diagnostic
  # rather than incorrectly declaring optimizer failure.
  converged <- opt$convergence == 0L && max(abs(score)) < 1e-6 && min(eig) > 0
  if (verbose) cat(sprintf(
    "FMM-AR1-4W [common transitions, type error]: ll=%.3f score=%.2e rank=%d/%d minEig=%.2e starts=%d\n",
    loglik, max(abs(score)), rank, K, min(eig), nrow(candidate_table)))
  list(params = params, eta = eta, loglik = loglik, converged = converged,
       identified = rank == K, model_type = "common_transitions_type_error",
       stationary = FALSE, cell_probabilities = probs, candidates = candidate_table,
       diagnostics = list(optimizer_code = opt$convergence,
         max_abs_score = max(abs(score)), information_min_eigenvalue = min(eig),
         information_condition = max(eig) / min(eig),
         weak_logit_curvature = min(eig) < 1e-9,
         probability_jacobian_rank = rank, parameter_count = K,
         start_count = nrow(candidate_table)),
       sample = list(n = attr(cells, "n_obs"),
                     weight_sum = attr(cells, "weight_sum")),
       estimator = "direct_sixteen_cell_two_type_fmm_mle")
}

.fmm_ctte_quantities <- function(eta) {
  p <- resolve_fmm_ctte_labels(unpack_fmm_ctte_4w(eta))
  c(reliable_share = unname(p$phi), initial_employment = unname(p$alpha),
    entry_rate = unname(p$theta0),
    exit_rate = unname(1 - p$theta1), A_pi = unname(p$pi_A), B_pi = unname(p$pi_B),
    aggregate_initial_employment = unname(p$alpha),
    steady_state_employment = stationary_alpha_ar1(p$theta0, p$theta1),
    aggregate_pi = unname(p$phi * p$pi_A + (1 - p$phi) * p$pi_B))
}

analytical_se_fmm_ctte_4w <- function(df, fit, finite_sample = TRUE) {
  if (!fit$identified || !fit$converged)
    stop("common-transition FMM fit must be identified and converged")
  cells <- collapse_ar1_4w_cells(df); eta <- fit$eta
  logcells <- function(z) setNames(log(pmax(
    fmm_ctte_4w_cell_probabilities(unpack_fmm_ctte_4w(z)), 1e-300)),
    paste0("cell", 1:16))
  scores <- .ar1_4w_jacobian(logcells, eta)
  avg_nll <- function(z) -sum(cells$weight * logcells(z)) / sum(cells$weight)
  bread <- optimHess(eta, avg_nll); bread <- (bread + t(bread)) / 2
  meat <- matrix(0, length(eta), length(eta), dimnames = list(names(eta), names(eta)))
  for (jj in 1:16) meat <- meat + cells$weight_sq[jj] / sum(cells$weight)^2 *
    tcrossprod(scores[jj, ])
  inv <- solve(bread); vcov_eta <- inv %*% meat %*% inv
  n <- attr(cells, "n_obs"); K <- length(eta)
  if (finite_sample && n > K) vcov_eta <- vcov_eta * n / (n - K)
  estimates <- .fmm_ctte_quantities(eta)
  delta <- .ar1_4w_jacobian(.fmm_ctte_quantities, eta)
  vcov_q <- delta %*% vcov_eta %*% t(delta)
  list(summary = data.frame(quantity = names(estimates),
      estimate = unname(estimates), se = sqrt(pmax(diag(vcov_q), 0))),
    vcov_eta = vcov_eta, vcov_quantities = vcov_q, bread = bread,
    meat = meat, delta_jacobian = delta,
    method = "individual_survey_weighted_sandwich_delta")
}
