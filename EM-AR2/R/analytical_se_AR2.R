# ==============================================================================
# EM-AR2: observed-information sandwich covariance and delta-method SEs
# ==============================================================================

.ar2_numeric_jacobian <- function(fn, x, rel_step = 1e-5) {
  f0 <- fn(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(x),
              dimnames = list(names(f0), names(x)))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  J
}

.ar2_model_type <- function(fit) {
  if (!is.null(fit$params$pi0) || !is.null(fit$params$pi1)) return("asymmetric")
  if (!is.null(fit$params$pi) && fit$params$pi > 0) return("symmetric")
  "none"
}

.ar2_eta_names <- function(model_type) {
  out <- c("alpha10", "alpha01", "alpha11",
           "p00", "gap10", "p01", "gap11")
  if (model_type == "symmetric") out <- c(out, "pi")
  if (model_type == "asymmetric") out <- c(out, "pi0", "pi1")
  out
}

pack_ar2_params <- function(params, model_type) {
  clamp <- function(x, eps = 1e-8) pmin(pmax(x, eps), 1 - eps)
  alpha <- params$alpha
  if (is.null(alpha))
    alpha <- stationary_ar2(params$theta0, params$theta01,
                            params$theta1, params$theta10)
  alpha <- alpha[c("00", "10", "01", "11")]
  p00 <- params$theta0
  p10 <- params$theta0 + params$theta01
  p01 <- 1 - params$theta1 - params$theta10
  p11 <- 1 - params$theta1
  out <- c(
    alpha10 = log(alpha["10"] / alpha["00"]),
    alpha01 = log(alpha["01"] / alpha["00"]),
    alpha11 = log(alpha["11"] / alpha["00"]),
    p00 = qlogis(clamp(p00)),
    gap10 = qlogis(clamp((p10 - p00) / (1 - p00))),
    p01 = qlogis(clamp(p01)),
    gap11 = qlogis(clamp((p11 - p01) / (1 - p01)))
  )
  names(out) <- .ar2_eta_names("none")
  if (model_type == "symmetric") out <- c(out, pi = qlogis(clamp(params$pi / 0.5)))
  if (model_type == "asymmetric") {
    out <- c(out,
             pi0 = qlogis(clamp(params$pi0 / 0.5)),
             pi1 = qlogis(clamp(params$pi1 / 0.5)))
  }
  setNames(as.numeric(out), .ar2_eta_names(model_type))
}

unpack_ar2_eta <- function(eta, model_type) {
  eta <- setNames(as.numeric(eta), .ar2_eta_names(model_type))
  alpha_raw <- c(1, exp(eta[c("alpha10", "alpha01", "alpha11")]))
  alpha <- alpha_raw / sum(alpha_raw)
  names(alpha) <- c("00", "10", "01", "11")

  p00 <- plogis(eta["p00"])
  p10 <- p00 + (1 - p00) * plogis(eta["gap10"])
  p01 <- plogis(eta["p01"])
  p11 <- p01 + (1 - p01) * plogis(eta["gap11"])
  out <- list(
    alpha = alpha,
    theta0 = unname(p00),
    theta01 = unname(p10 - p00),
    theta1 = unname(1 - p11),
    theta10 = unname(p11 - p01),
    asymmetric = identical(model_type, "asymmetric")
  )
  if (model_type == "none") out$pi <- 0
  if (model_type == "symmetric") out$pi <- unname(0.5 * plogis(eta["pi"]))
  if (model_type == "asymmetric") {
    out$pi0 <- unname(0.5 * plogis(eta["pi0"]))
    out$pi1 <- unname(0.5 * plogis(eta["pi1"]))
  }
  out
}

.ar2_cell_prob_vector <- function(eta, model_type) {
  model_cell_probs_ar2(unpack_ar2_eta(eta, model_type))$model_prob
}

.ar2_reported_quantities <- function(eta, model_type) {
  p <- unpack_ar2_eta(eta, model_type)
  implied <- implied_ar2(p, model_type)
  out <- c(
    alpha_00 = unname(p$alpha["00"]), alpha_10 = unname(p$alpha["10"]),
    alpha_01 = unname(p$alpha["01"]), alpha_11 = unname(p$alpha["11"]),
    theta0 = p$theta0, theta01 = p$theta01,
    theta1 = p$theta1, theta10 = p$theta10
  )
  if (model_type == "symmetric") out <- c(out, pi = p$pi)
  if (model_type == "asymmetric") out <- c(out, pi0 = p$pi0, pi1 = p$pi1)
  c(out,
    p_00 = implied$p_00, p_10 = implied$p_10,
    p_01 = implied$p_01, p_11 = implied$p_11,
    employment_rate = implied$employment_rate)
}

#' Analytical standard errors for an AR(2) EM fit
#'
#' Uses the observed-data Hessian as bread, individual survey-weighted scores
#' as meat, and a numerical delta method for reported transition quantities.
#' No bootstrap resampling is involved.
analytical_se_ar2 <- function(df, fit, model_type = NULL,
                              finite_sample = TRUE, rel_step = 1e-5) {
  if (is.null(model_type)) model_type <- .ar2_model_type(fit)
  if (!model_type %in% c("none", "symmetric", "asymmetric"))
    stop("analytical_se_ar2: invalid model_type")
  cells <- collapse_ar2_cells(df)
  eta <- pack_ar2_params(fit$params, model_type)
  probs <- pmax(.ar2_cell_prob_vector(eta, model_type), 1e-300)
  cell_index <- 1L + cells$y1 + 2L * cells$y2 + 4L * cells$y3 + 8L * cells$y4
  weight_sum <- sum(cells$weight)

  average_nll <- function(z) {
    pr <- pmax(.ar2_cell_prob_vector(z, model_type), 1e-300)
    -sum(cells$weight * log(pr[cell_index])) / weight_sum
  }
  log_cell <- function(z) {
    setNames(log(pmax(.ar2_cell_prob_vector(z, model_type), 1e-300)),
             apply(expand.grid(y1 = 0:1, y2 = 0:1,
                               y3 = 0:1, y4 = 0:1), 1L, paste0, collapse = ""))
  }
  scores_all <- .ar2_numeric_jacobian(log_cell, eta, rel_step)
  bread <- optimHess(eta, average_nll)
  bread <- (bread + t(bread)) / 2
  dimnames(bread) <- list(names(eta), names(eta))

  meat <- matrix(0, length(eta), length(eta),
                 dimnames = list(names(eta), names(eta)))
  for (j in seq_len(nrow(cells))) {
    s <- scores_all[cell_index[j], ]
    meat <- meat + cells$weight_sq[j] / weight_sum^2 * tcrossprod(s)
  }
  eig <- eigen(bread, symmetric = TRUE)
  tol <- max(eig$values) * 1e-10
  rank <- sum(eig$values > tol)
  if (rank < length(eta) || min(eig$values) <= 0)
    stop("analytical_se_ar2: observed information is not positive definite; ",
         "the model is locally unidentified or the EM has not reached an interior maximum")
  bread_inv <- solve(bread)
  vcov_eta <- bread_inv %*% meat %*% bread_inv
  K <- length(eta)
  n <- attr(cells, "n_obs")
  if (finite_sample && n > K) vcov_eta <- vcov_eta * n / (n - K)
  dimnames(vcov_eta) <- list(names(eta), names(eta))

  qfun <- function(z) .ar2_reported_quantities(z, model_type)
  estimates <- qfun(eta)
  delta <- .ar2_numeric_jacobian(qfun, eta, rel_step)
  vcov_q <- delta %*% vcov_eta %*% t(delta)
  se <- sqrt(pmax(diag(vcov_q), 0))
  score <- colSums(scores_all[cell_index, , drop = FALSE] * cells$weight) / weight_sum

  list(
    summary = data.frame(quantity = names(estimates),
                         estimate = unname(estimates), se = unname(se)),
    vcov_eta = vcov_eta,
    vcov_quantities = vcov_q,
    bread = bread,
    meat = meat,
    delta_jacobian = delta,
    diagnostics = list(
      max_abs_average_score = max(abs(score)),
      information_min_eigenvalue = min(eig$values),
      information_condition = max(eig$values) / min(eig$values),
      information_rank = rank
    ),
    model_type = model_type,
    n_obs = n,
    method = "individual_survey_weighted_sandwich_delta"
  )
}
