# ==============================================================================
# EM-baseline: exact observed-data MLE on the eight history cells
# ==============================================================================

.baseline_observed_histories <- function() {
  as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))[, 1:3]
}

collapse_baseline_cells <- function(df) {
  .validate_panel_df(df)
  for (nm in c("y1", "y2", "y3")) {
    if (anyNA(df[[nm]]) || !all(df[[nm]] %in% c(0, 1)))
      stop("collapse_baseline_cells: y1/y2/y3 must be non-missing binary values")
  }
  if (anyNA(df$weight) || any(!is.finite(df$weight)) || any(df$weight <= 0))
    stop("collapse_baseline_cells: weights must be finite and strictly positive")

  index <- 1L + as.integer(df$y1) + 2L * as.integer(df$y2) +
    4L * as.integer(df$y3)
  weight_by_cell <- numeric(8L)
  weight_sq_by_cell <- numeric(8L)
  weight_sum <- rowsum(df$weight, index, reorder = FALSE)
  weight_sq_sum <- rowsum(df$weight^2, index, reorder = FALSE)
  weight_by_cell[as.integer(rownames(weight_sum))] <- weight_sum[, 1L]
  weight_sq_by_cell[as.integer(rownames(weight_sq_sum))] <- weight_sq_sum[, 1L]
  list(
    histories = .baseline_observed_histories(),
    weight = weight_by_cell,
    weight_sq = weight_sq_by_cell,
    count = as.numeric(tabulate(index, nbins = 8L)),
    n = nrow(df),
    weight_sum = sum(df$weight)
  )
}

baseline_sample_signature <- function(cells) {
  paste(c(cells$n, format(cells$weight, digits = 17, scientific = TRUE)),
        collapse = "|")
}

.baseline_eta_names <- function(model_type, stationary) {
  out <- c("theta0", "theta1")
  if (!stationary) out <- c(out, "alpha")
  if (model_type == "symmetric") out <- c(out, "pi")
  if (model_type == "asymmetric") out <- c(out, "pi0", "pi1")
  out
}

pack_baseline_params <- function(params, model_type, stationary) {
  nms <- .baseline_eta_names(model_type, stationary)
  out <- numeric(length(nms)); names(out) <- nms
  clamp <- function(x, lo = 1e-8, hi = 1 - 1e-8) pmin(pmax(x, lo), hi)
  out["theta0"] <- qlogis(clamp(params$theta0))
  out["theta1"] <- qlogis(clamp(params$theta1))
  if (!stationary) out["alpha"] <- qlogis(clamp(params$alpha))
  if (model_type == "symmetric") out["pi"] <- qlogis(clamp(params$pi / 0.5))
  if (model_type == "asymmetric") {
    out["pi0"] <- qlogis(clamp(params$pi0 / 0.5))
    out["pi1"] <- qlogis(clamp(params$pi1 / 0.5))
  }
  out
}

unpack_baseline_eta <- function(eta, model_type, stationary) {
  eta <- setNames(as.numeric(eta), .baseline_eta_names(model_type, stationary))
  theta0 <- plogis(eta["theta0"])
  theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha(theta0, theta1) else plogis(eta["alpha"])
  out <- list(alpha = unname(alpha), theta0 = unname(theta0), theta1 = unname(theta1))
  if (model_type == "symmetric") out$pi <- unname(0.5 * plogis(eta["pi"]))
  if (model_type == "asymmetric") {
    out$pi0 <- unname(0.5 * plogis(eta["pi0"]))
    out$pi1 <- unname(0.5 * plogis(eta["pi1"]))
  }
  out
}

baseline_cell_probabilities <- function(params, model_type = "symmetric") {
  .validate_model_type(model_type)
  observed <- .baseline_observed_histories()
  latent <- latent_histories()
  prior <- prior_over_histories(latent, params$theta1, params$theta0, params$alpha)
  probs <- numeric(8L)

  pi <- params$pi %||% 0
  pi0 <- params$pi0 %||% 0
  pi1 <- params$pi1 %||% 0
  for (j in seq_len(8L)) {
    emission <- rep(1, 8L)
    for (tt in seq_len(3L)) {
      s <- observed[j, tt]
      h <- latent[, tt]
      if (model_type == "none") {
        emission <- emission * as.numeric(s == h)
      } else if (model_type == "symmetric") {
        emission <- emission * ifelse(s == h, 1 - pi, pi)
      } else {
        wave_prob <- ifelse(
          h == 0,
          ifelse(s == 1, pi0, 1 - pi0),
          ifelse(s == 0, pi1, 1 - pi1)
        )
        emission <- emission * wave_prob
      }
    }
    probs[j] <- sum(prior * emission)
  }
  probs / sum(probs)
}

.baseline_average_nll <- function(eta, cells, model_type, stationary) {
  params <- unpack_baseline_eta(eta, model_type, stationary)
  probs <- pmax(baseline_cell_probabilities(params, model_type), 1e-300)
  -sum(cells$weight * log(probs)) / cells$weight_sum
}

.numeric_gradient <- function(fn, x, rel_step = 1e-5) {
  out <- numeric(length(x))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x; xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
    out[j] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  out
}

fit_baseline_mle <- function(df, model_type = "symmetric", stationary = TRUE,
                             starts = NULL, maxit = 2000L, reltol = 1e-12,
                             compute_gamma = FALSE, verbose = 1L) {
  .validate_panel_df(df); .validate_model_type(model_type)
  cells <- collapse_baseline_cells(df)
  if (is.null(starts)) starts <- list(init_params(model_type, stationary))
  if (!is.list(starts) || !length(starts)) stop("fit_baseline_mle: starts must be a non-empty list")

  fn <- function(z) .baseline_average_nll(z, cells, model_type, stationary)
  candidates <- vector("list", length(starts))
  for (s in seq_along(starts)) {
    eta0 <- pack_baseline_params(starts[[s]], model_type, stationary)
    opt <- tryCatch(
      optim(eta0, fn, method = "BFGS", hessian = FALSE,
            control = list(maxit = maxit, reltol = reltol)),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value)) candidates[[s]] <- opt
  }
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("fit_baseline_mle: all starts failed")
  opt <- candidates[[which.min(vapply(candidates, `[[`, numeric(1L), "value"))]]

  eta_hat <- setNames(opt$par, .baseline_eta_names(model_type, stationary))
  params <- unpack_baseline_eta(eta_hat, model_type, stationary)
  probs <- baseline_cell_probabilities(params, model_type)
  loglik <- sum(cells$weight * log(pmax(probs, 1e-300)))
  grad <- .numeric_gradient(fn, eta_hat)
  hessian <- optimHess(eta_hat, fn)
  eig <- eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values
  converged <- identical(opt$convergence, 0L) && max(abs(grad)) < 1e-6 && min(eig) > 1e-9

  if (verbose >= 1L) {
    cat(sprintf("MLE [%s, stationary=%s]: ll=%.4f, code=%d, max|score|=%.2e, minEig=%.2e\n",
                model_type, stationary, loglik, opt$convergence,
                max(abs(grad)), min(eig)))
  }

  gamma <- if (compute_gamma)
    e_step(df, params, model_type = model_type, validate = FALSE)$gamma else NULL
  list(
    params = params,
    eta = eta_hat,
    loglik = loglik,
    converged = converged,
    iterations = unname(opt$counts[["function"]]),
    gamma = gamma,
    model_type = model_type,
    stationary = stationary,
    cell_probabilities = probs,
    diagnostics = list(
      optimizer_code = opt$convergence,
      optimizer_message = opt$message,
      max_abs_score = max(abs(grad)),
      information_min_eigenvalue = min(eig),
      information_condition = max(eig) / min(eig),
      start_count = length(starts)
    ),
    sample = list(
      n = cells$n,
      weight_sum = cells$weight_sum,
      cell_weights = cells$weight,
      signature = baseline_sample_signature(cells),
      source_panel = "df_qlfs_A.rds"
    ),
    estimator = "direct_eight_cell_mle"
  )
}

check_baseline_nesting <- function(fits, tolerance = 1e-4) {
  comparisons <- list(
    c("none_stat", "sym_stat"), c("sym_stat", "asym_stat"),
    c("none_free", "sym_free"), c("sym_free", "asym_free"),
    c("none_stat", "none_free"), c("sym_stat", "sym_free"),
    c("asym_stat", "asym_free")
  )
  failures <- vapply(comparisons, function(x) {
    fits[[x[2L]]]$loglik < fits[[x[1L]]]$loglik - tolerance
  }, logical(1L))
  if (any(failures)) {
    bad <- vapply(comparisons[failures], function(x) paste(x, collapse = " -> "), character(1L))
    stop("Baseline likelihood nesting check failed: ", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}
