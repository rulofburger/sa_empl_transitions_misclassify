# ==============================================================================
# Exact observed-data MLE and identification diagnostics for the two-type FMM
# ==============================================================================

.fmm_eta_names <- function(model_type, stationary) {
  stopifnot(model_type %in% c("none", "symmetric"))
  out <- c("theta0_A", "theta1_A")
  if (!stationary) out <- c(out, "alpha_A")
  out <- c(out, "theta0_B", "theta1_B")
  if (!stationary) out <- c(out, "alpha_B")
  out <- c(out, "phi")
  if (model_type == "symmetric") out <- c(out, "pi")
  out
}

fmm_parameter_count <- function(model_type = "symmetric", stationary = TRUE) {
  length(.fmm_eta_names(model_type, stationary))
}

pack_fmm_params <- function(params, model_type = "symmetric", stationary = TRUE) {
  nms <- .fmm_eta_names(model_type, stationary)
  clamp <- function(x, lo = 1e-8, hi = 1 - 1e-8) pmin(pmax(x, lo), hi)
  out <- setNames(numeric(length(nms)), nms)
  for (nm in setdiff(nms, "pi")) out[nm] <- qlogis(clamp(params[[nm]]))
  if ("pi" %in% nms) out["pi"] <- qlogis(clamp(params$pi / 0.5))
  out
}

unpack_fmm_eta <- function(eta, model_type = "symmetric", stationary = TRUE) {
  eta <- setNames(as.numeric(eta), .fmm_eta_names(model_type, stationary))
  p <- as.list(plogis(eta))
  if (stationary) {
    p$alpha_A <- stationary_alpha(p$theta0_A, p$theta1_A)
    p$alpha_B <- stationary_alpha(p$theta0_B, p$theta1_B)
  }
  if (model_type == "symmetric") p$pi <- 0.5 * plogis(eta["pi"])
  p
}

.fmm_cell_cache <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      observed <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))[, 1:3]
      latent <- latent_histories()
      mismatch <- matrix(0L, nrow = 8L, ncol = 8L)
      for (tt in seq_len(3L)) mismatch <- mismatch + outer(observed[, tt], latent[, tt], "!=")
      cache <<- list(observed = observed, latent = latent, mismatch = mismatch)
    }
    cache
  }
})

.fmm_type_cell_probabilities <- function(theta0, theta1, alpha, pi = 0) {
  cache <- .fmm_cell_cache()
  latent <- cache$latent
  prior <- prior_over_histories(latent, theta1, theta0, alpha)
  emission <- (1 - pi)^(3L - cache$mismatch) * pi^cache$mismatch
  out <- as.vector(emission %*% prior)
  out / sum(out)
}

fmm_cell_probabilities <- function(params, model_type = "symmetric") {
  stopifnot(model_type %in% c("none", "symmetric"))
  pi <- if (model_type == "symmetric") params$pi else 0
  pA <- .fmm_type_cell_probabilities(
    params$theta0_A, params$theta1_A, params$alpha_A, pi
  )
  pB <- .fmm_type_cell_probabilities(
    params$theta0_B, params$theta1_B, params$alpha_B, pi
  )
  params$phi * pA + (1 - params$phi) * pB
}

.fmm_numeric_jacobian <- function(fn, x, rel_step = 1e-5) {
  f0 <- fn(x)
  out <- matrix(NA_real_, nrow = length(f0), ncol = length(x))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x; xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
    out[, j] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  dimnames(out) <- list(names(f0), names(x))
  out
}

.fmm_matrix_rank <- function(x, tolerance = NULL) {
  singular <- svd(x, nu = 0L, nv = 0L)$d
  if (!length(singular)) return(0L)
  if (is.null(tolerance))
    tolerance <- max(dim(x)) * max(singular) * sqrt(.Machine$double.eps)
  sum(singular > tolerance)
}

.resolve_fmm_labels <- function(params) {
  if (params$theta1_B <= params$theta1_A) return(params)
  out <- params
  for (stem in c("theta0", "theta1", "alpha")) {
    out[[paste0(stem, "_A")]] <- params[[paste0(stem, "_B")]]
    out[[paste0(stem, "_B")]] <- params[[paste0(stem, "_A")]]
  }
  out$phi <- 1 - params$phi
  out
}

fit_fmm_mle <- function(df, model_type = "symmetric", stationary = TRUE,
                        starts, maxit = 4000L, reltol = 1e-12,
                        verbose = 1L) {
  .validate_panel_df(df)
  if (!model_type %in% c("none", "symmetric"))
    stop("fit_fmm_mle: model_type must be 'none' or 'symmetric'")
  if (!is.list(starts) || !length(starts))
    stop("fit_fmm_mle: starts must be a non-empty list")

  # Reuse the validated baseline cell collapse when available.
  if (!exists("collapse_baseline_cells", mode = "function"))
    stop("fit_fmm_mle requires EM-baseline/R/mle_baseline.R to be sourced first")
  cells <- collapse_baseline_cells(df)
  fn <- function(eta) {
    p <- unpack_fmm_eta(eta, model_type, stationary)
    probs <- pmax(fmm_cell_probabilities(p, model_type), 1e-300)
    -sum(cells$weight * log(probs)) / cells$weight_sum
  }

  candidates <- lapply(starts, function(start) {
    eta0 <- pack_fmm_params(start, model_type, stationary)
    opt <- tryCatch(
      optim(eta0, fn, method = "BFGS",
            control = list(maxit = maxit, reltol = reltol)),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$value)) return(NULL)
    params <- .resolve_fmm_labels(unpack_fmm_eta(opt$par, model_type, stationary))
    eta <- pack_fmm_params(params, model_type, stationary)
    probs <- fmm_cell_probabilities(params, model_type)
    list(opt = opt, params = params, eta = eta,
         loglik = sum(cells$weight * log(pmax(probs, 1e-300))))
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("fit_fmm_mle: every optimization start failed")
  best <- candidates[[which.max(vapply(candidates, `[[`, numeric(1L), "loglik"))]]

  eta <- best$eta
  avg_nll <- function(z) fn(z)
  grad <- if (exists(".numeric_gradient", mode = "function")) {
    .numeric_gradient(avg_nll, eta)
  } else {
    as.numeric(.fmm_numeric_jacobian(function(z) avg_nll(z), eta))
  }
  hessian <- optimHess(eta, avg_nll)
  hessian <- (hessian + t(hessian)) / 2
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  # Seven cell probabilities are independent because the eighth is their sum
  # complement. Rank below K is direct evidence of local nonidentification.
  prob_fn <- function(z) {
    out <- fmm_cell_probabilities(unpack_fmm_eta(z, model_type, stationary), model_type)[1:7]
    setNames(out, sprintf("cell_%d", 1:7))
  }
  probability_jacobian <- .fmm_numeric_jacobian(prob_fn, eta)
  jacobian_rank <- .fmm_matrix_rank(probability_jacobian)
  K <- length(eta)
  identified <- jacobian_rank == K
  converged <- identical(best$opt$convergence, 0L) && max(abs(grad)) < 1e-6

  candidate_table <- do.call(rbind, lapply(candidates, function(x) {
    p <- x$params
    data.frame(loglik = x$loglik, theta0_A = p$theta0_A,
               exit_A = 1 - p$theta1_A, alpha_A = p$alpha_A,
               theta0_B = p$theta0_B, exit_B = 1 - p$theta1_B,
               alpha_B = p$alpha_B, phi = p$phi,
               pi = p$pi %||% NA_real_, stringsAsFactors = FALSE)
  }))
  candidate_table <- candidate_table[order(candidate_table$loglik, decreasing = TRUE), ]

  if (verbose >= 1L) {
    cat(sprintf(
      "FMM MLE [%s, stationary=%s]: ll=%.3f code=%d max|score|=%.2e rank=%d/%d minEig=%.2e identified=%s\n",
      model_type, stationary, best$loglik, best$opt$convergence,
      max(abs(grad)), jacobian_rank, K, min(eig), identified
    ))
  }

  list(params = best$params, eta = eta, loglik = best$loglik,
       converged = converged, identified = identified,
       model_type = model_type, stationary = stationary,
       cell_probabilities = fmm_cell_probabilities(best$params, model_type),
       candidates = candidate_table,
       diagnostics = list(optimizer_code = best$opt$convergence,
                          max_abs_score = max(abs(grad)),
                          information_min_eigenvalue = min(eig),
                          information_eigenvalues = eig,
                          probability_jacobian_rank = jacobian_rank,
                          parameter_count = K,
                          start_count = length(starts)),
       sample = list(n = cells$n, weight_sum = cells$weight_sum,
                     signature = baseline_sample_signature(cells)),
       estimator = "direct_eight_cell_fmm_mle")
}
