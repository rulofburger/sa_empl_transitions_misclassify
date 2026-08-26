# Combined three-wave model with covariate-dependent transitions and
# reliability-dependent symmetric misclassification.

.as_reliability_design <- function(Z, N = NULL) {
  if (is.array(Z) && length(dim(Z)) == 3L && dim(Z)[2L] == 3L) {
    out <- lapply(seq_len(3L), function(tt) Z[, tt, , drop = FALSE][, 1L, ])
  } else if (is.list(Z) && length(Z) == 3L) {
    out <- unname(Z)
  } else {
    stop("Z must be an N x 3 x q array or a list of three matrices")
  }
  out <- lapply(out, function(z) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    z
  })
  if (!is.null(N) && any(vapply(out, nrow, integer(1L)) != N))
    stop("Each reliability design must have one row per observation")
  if (length(unique(vapply(out, ncol, integer(1L)))) != 1L ||
      !all(vapply(out[-1L], function(z) identical(colnames(z), colnames(out[[1L]])),
                  logical(1L))))
    stop("Reliability designs must have identical columns")
  out
}

.covrel_layout <- function(X) {
  Xt <- .as_transition_design(X)
  active0 <- attr(Xt$X12, "entry_active") %||% rep(TRUE, ncol(Xt$X12))
  active1 <- attr(Xt$X12, "persistence_active") %||% rep(TRUE, ncol(Xt$X12))
  list(X = Xt, active0 = active0, active1 = active1,
       p0 = sum(active0), p1 = sum(active1))
}

pack_covariate_reliability <- function(params, X) {
  layout <- .covrel_layout(X)
  out <- c(params$beta0[layout$active0], params$beta1[layout$active1],
           initial_employment_logit = qlogis(params$alpha), params$delta)
  names(out) <- c(paste0("entry_", colnames(layout$X$X12)[layout$active0]),
                  paste0("persistence_", colnames(layout$X$X12)[layout$active1]),
                  "initial_employment_logit", names(params$delta))
  out
}

unpack_covariate_reliability <- function(eta, X, delta_names) {
  layout <- .covrel_layout(X)
  p <- ncol(layout$X$X12); q <- length(delta_names)
  if (length(eta) != layout$p0 + layout$p1 + 1L + q)
    stop("Parameter vector has the wrong length")
  beta0 <- beta1 <- numeric(p)
  names(beta0) <- names(beta1) <- colnames(layout$X$X12)
  beta0[layout$active0] <- eta[seq_len(layout$p0)]
  beta1[layout$active1] <- eta[layout$p0 + seq_len(layout$p1)]
  offset <- layout$p0 + layout$p1
  delta <- eta[offset + 1L + seq_len(q)]
  names(delta) <- delta_names
  list(beta0 = beta0, beta1 = beta1,
       alpha = plogis(eta[offset + 1L]), delta = delta)
}

.covrel_components <- function(eta, data, retain_scores = FALSE,
                               retain_posterior = FALSE) {
  df <- data$df; Xt <- .as_transition_design(data$X, nrow(df))
  Z <- .as_reliability_design(data$Z, nrow(df))
  params <- unpack_covariate_reliability(eta, Xt, colnames(Z[[1L]]))
  h <- .hmat_cache(); N <- nrow(df); H <- nrow(h)

  theta0 <- lapply(Xt, function(x)
    pmin(pmax(pnorm(as.vector(x %*% params$beta0)), 1e-9), 1 - 1e-9))
  theta1 <- lapply(Xt, function(x)
    pmin(pmax(pnorm(as.vector(x %*% params$beta1)), 1e-9), 1 - 1e-9))
  pi_wave <- lapply(Z, function(z)
    pmin(pmax(0.5 * plogis(as.vector(z %*% params$delta)), 1e-9), 0.5 - 1e-9))

  log_joint <- matrix(0, N, H)
  log_joint[, h[, 1L] == 1L] <- log(params$alpha)
  log_joint[, h[, 1L] == 0L] <- log1p(-params$alpha)
  log_joint <- log_joint +
    .log_markov_trans_indiv(h[, 1L], h[, 2L], theta1[[1L]], theta0[[1L]]) +
    .log_markov_trans_indiv(h[, 2L], h[, 3L], theta1[[2L]], theta0[[2L]])

  mismatch_masks <- vector("list", 3L)
  for (tt in seq_len(3L)) {
    mismatch <- outer(df[[paste0("y", tt)]], h[, tt], "!=")
    mismatch_masks[[tt]] <- mismatch
    emit <- matrix(log1p(-pi_wave[[tt]]), N, H)
    emit[mismatch] <- matrix(log(pi_wave[[tt]]), N, H)[mismatch]
    log_joint <- log_joint + emit
  }
  mx <- apply(log_joint, 1L, max)
  denom <- rowSums(exp(log_joint - mx))
  loglik_i <- mx + log(denom)
  posterior <- exp(log_joint - mx) / denom
  out <- list(loglik_i = loglik_i, params = params,
              theta0 = theta0, theta1 = theta1, pi_wave = pi_wave)
  if (!retain_scores && !retain_posterior) return(out)
  if (retain_posterior) out$posterior <- posterior
  if (!retain_scores) return(out)

  layout <- .covrel_layout(Xt)
  blocks0 <- vector("list", 2L); blocks1 <- vector("list", 2L)
  for (tt in seq_len(2L)) {
    hf <- h[, tt]; ht <- h[, tt + 1L]
    risk0 <- as.vector(posterior %*% as.integer(hf == 0L))
    risk1 <- as.vector(posterior %*% as.integer(hf == 1L))
    succ0 <- as.vector(posterior %*% as.integer(hf == 0L & ht == 1L))
    succ1 <- as.vector(posterior %*% as.integer(hf == 1L & ht == 1L))
    eta0 <- as.vector(Xt[[tt]] %*% params$beta0)
    eta1 <- as.vector(Xt[[tt]] %*% params$beta1)
    q0 <- theta0[[tt]]; q1 <- theta1[[tt]]
    s0 <- dnorm(eta0) * (succ0 / q0 - (risk0 - succ0) / (1 - q0))
    s1 <- dnorm(eta1) * (succ1 / q1 - (risk1 - succ1) / (1 - q1))
    blocks0[[tt]] <- Xt[[tt]][, layout$active0, drop = FALSE] * s0
    blocks1[[tt]] <- Xt[[tt]][, layout$active1, drop = FALSE] * s1
  }
  score <- cbind(blocks0[[1L]] + blocks0[[2L]],
                 blocks1[[1L]] + blocks1[[2L]],
                 initial_employment_logit =
                   as.vector(posterior %*% h[, 1L]) - params$alpha)
  delta_score <- matrix(0, N, length(params$delta),
                        dimnames = list(NULL, names(params$delta)))
  for (tt in seq_len(3L)) {
    expected_mismatch <- rowSums(posterior * mismatch_masks[[tt]])
    pi <- pi_wave[[tt]]
    scalar <- (expected_mismatch - pi) * (1 - 2 * pi) / (1 - pi)
    delta_score <- delta_score + Z[[tt]] * scalar
  }
  out$scores <- cbind(score, delta_score)
  out$posterior <- posterior
  out
}

collapse_covariate_reliability <- function(df, X, Z) {
  Xt <- .as_transition_design(X, nrow(df))
  Zt <- .as_reliability_design(Z, nrow(df))
  key_df <- data.frame(y1 = df$y1, y2 = df$y2, y3 = df$y3,
                       Xt$X12, Xt$X23, do.call(cbind, Zt), check.names = FALSE)
  key <- do.call(paste, c(key_df, sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  copy_attrs <- function(mat) {
    attr(mat, "entry_active") <- attr(Xt$X12, "entry_active")
    attr(mat, "persistence_active") <- attr(Xt$X12, "persistence_active")
    mat
  }
  out_X <- list(X12 = copy_attrs(Xt$X12[first, , drop = FALSE]),
                X23 = copy_attrs(Xt$X23[first, , drop = FALSE]))
  out_Z <- lapply(Zt, function(z) z[first, , drop = FALSE])
  list(df = data.frame(y1 = df$y1[first], y2 = df$y2[first],
                       y3 = df$y3[first],
                       weight = as.vector(rowsum(df$weight, group,
                                                 reorder = FALSE))),
       X = out_X, Z = out_Z,
       weight_sq = as.vector(rowsum(df$weight^2, group, reorder = FALSE)),
       count = as.vector(rowsum(rep(1, nrow(df)), group, reorder = FALSE)),
       n_original = nrow(df), weight_sum = sum(df$weight))
}

fit_covariate_reliability <- function(data, start, maxit = 2000L,
                                      reltol = 1e-10) {
  eta0 <- if (is.list(start)) pack_covariate_reliability(start, data$X) else start
  weight_total <- sum(data$df$weight)
  cache <- new.env(parent = emptyenv())
  cache$eta <- NULL
  evaluate <- function(z) {
    if (!is.null(cache$eta) && identical(as.numeric(z), cache$eta))
      return(cache$value)
    value <- .covrel_components(z, data, retain_scores = TRUE)
    cache$eta <- as.numeric(z); cache$value <- value
    value
  }
  objective <- function(z)
    -sum(data$df$weight * evaluate(z)$loglik_i) / weight_total
  gradient <- function(z) {
    detail <- evaluate(z)
    -colSums(detail$scores * data$df$weight) / weight_total
  }
  opt <- optim(eta0, objective, gradient, method = "BFGS",
               control = list(maxit = maxit, reltol = reltol,
                              trace = 1L, REPORT = 10L))
  detail <- .covrel_components(opt$par, data, retain_scores = TRUE,
                               retain_posterior = TRUE)
  score <- colSums(detail$scores * data$df$weight) / weight_total
  list(params = detail$params, eta = setNames(opt$par, names(eta0)),
       loglik = sum(data$df$weight * detail$loglik_i),
       converged = opt$convergence == 0L,
       optimizer_code = opt$convergence,
       iterations = unname(opt$counts[["function"]]),
       max_abs_score = max(abs(score)), scores = detail$scores,
       posterior = detail$posterior)
}

.covrel_quantities <- function(eta, data, components = NULL) {
  detail <- components %||% .covrel_components(eta, data, retain_posterior = TRUE)
  p <- detail$params; post <- detail$posterior; h <- .hmat_cache()
  w <- data$df$weight
  risk0 <- list(as.vector(post %*% as.integer(h[, 1L] == 0L)),
                as.vector(post %*% as.integer(h[, 2L] == 0L)))
  risk1 <- lapply(risk0, function(z) 1 - z)
  entry_num <- risk0[[1L]] * detail$theta0[[1L]] +
    risk0[[2L]] * detail$theta0[[2L]]
  exit_prob <- lapply(detail$theta1, function(z) 1 - z)
  exit_num <- risk1[[1L]] * exit_prob[[1L]] +
    risk1[[2L]] * exit_prob[[2L]]
  weighted_ratio <- function(a, d) sum(w * a) / sum(w * d)
  scenario <- function(...) {
    values <- c(...); z <- setNames(rep(0, length(p$delta)), names(p$delta))
    z[names(values)] <- values
    unname(0.5 * plogis(sum(z * p$delta)))
  }
  q <- c(
    mean_entry_rate = weighted_ratio(entry_num, risk0[[1L]] + risk0[[2L]]),
    mean_exit_rate = weighted_ratio(exit_num, risk1[[1L]] + risk1[[2L]]),
    mean_misclassification_rate =
      sum(w * Reduce(`+`, detail$pi_wave)) / (3 * sum(w)),
    mean_employment_rate = weighted_ratio((risk1[[1L]] + risk1[[2L]]) / 2,
                                           rep(1, length(w))),
    initial_employment = unname(p$alpha),
    pi_base = scenario(error_intercept = 1),
    pi_age_mild = scenario(error_intercept = 1, age_inconsistency = 1),
    pi_age_severe = scenario(error_intercept = 1, age_inconsistency = 1,
                             large_age_inconsistency = 1),
    pi_education_mild = scenario(error_intercept = 1,
                                 education_inconsistency = 1),
    pi_education_severe = scenario(error_intercept = 1,
      education_inconsistency = 1, large_education_inconsistency = 1),
    pi_both_mild = scenario(error_intercept = 1, age_inconsistency = 1,
                            education_inconsistency = 1),
    pi_both_both_severe = scenario(error_intercept = 1,
      age_inconsistency = 1, education_inconsistency = 1,
      large_age_inconsistency = 1, large_education_inconsistency = 1),
    pi_B_not_C = scenario(error_intercept = 1, panel_B_not_C = 1),
    pi_A_not_B = scenario(error_intercept = 1, panel_A_not_B = 1)
  )
  q
}

analytical_se_covariate_reliability <- function(data, fit,
                                                step_scale = 1e-5) {
  wt <- sum(data$df$weight); eta <- fit$eta; K <- length(eta)
  objective <- function(z)
    -sum(data$df$weight * .covrel_components(z, data)$loglik_i) / wt
  gradient <- function(z) {
    x <- .covrel_components(z, data, retain_scores = TRUE)
    -colSums(x$scores * data$df$weight) / wt
  }
  bread <- optimHess(eta, objective, gradient)
  bread <- (bread + t(bread)) / 2
  eig <- eigen(bread, symmetric = TRUE)
  threshold <- max(eig$values) * 1e-9
  if (min(eig$values) <= threshold)
    stop("Combined model information matrix is not positive definite")
  inv_bread <- solve(bread)
  meat <- crossprod(fit$scores * sqrt(data$weight_sq)) / wt^2
  vcov_eta <- inv_bread %*% meat %*% inv_bread
  if (data$n_original > K)
    vcov_eta <- vcov_eta * data$n_original / (data$n_original - K)
  dimnames(vcov_eta) <- list(names(eta), names(eta))
  average_score <- colSums(fit$scores * data$df$weight) / wt
  newton_correction <- as.vector(inv_bread %*% average_score)
  names(newton_correction) <- names(eta)

  q0 <- .covrel_quantities(eta, data)
  jac <- matrix(NA_real_, length(q0), K,
                dimnames = list(names(q0), names(eta)))
  for (j in seq_len(K)) {
    hh <- step_scale * max(1, abs(eta[j])); plus <- minus <- eta
    plus[j] <- plus[j] + hh; minus[j] <- minus[j] - hh
    jac[, j] <- (.covrel_quantities(plus, data) -
                  .covrel_quantities(minus, data)) / (2 * hh)
  }
  q_vcov <- jac %*% vcov_eta %*% t(jac)
  param_est <- c(fit$params$beta0[.covrel_layout(data$X)$active0],
                 fit$params$beta1[.covrel_layout(data$X)$active1],
                 initial_employment = fit$params$alpha, fit$params$delta)
  param_names <- c(names(eta)[seq_len(K - length(fit$params$delta) - 1L)],
                   "initial_employment", names(fit$params$delta))
  names(param_est) <- param_names
  deriv <- c(rep(1, K - length(fit$params$delta) - 1L),
             fit$params$alpha * (1 - fit$params$alpha),
             rep(1, length(fit$params$delta)))
  param_vcov <- diag(deriv) %*% vcov_eta %*% diag(deriv)
  summary <- rbind(
    data.frame(quantity = names(param_est), estimate = unname(param_est),
               std_error = sqrt(pmax(diag(param_vcov), 0))),
    data.frame(quantity = names(q0), estimate = unname(q0),
               std_error = sqrt(pmax(diag(q_vcov), 0))))
  list(summary = summary, vcov_eta = vcov_eta, vcov_quantities = q_vcov,
       delta_jacobian = jac, bread = bread, meat = meat,
       diagnostics = list(N = data$n_original, cells = nrow(data$df), K = K,
         information_rank = sum(eig$values > threshold),
         min_information_eigenvalue = min(eig$values),
         information_condition = max(eig$values) / min(eig$values),
         max_abs_score = fit$max_abs_score,
         max_abs_newton_correction = max(abs(newton_correction)),
         newton_correction = newton_correction))
}
