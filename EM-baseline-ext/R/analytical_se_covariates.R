# Approximate analytical inference for free-alpha covariate models.
#
# The covariance estimator is an individual-level, survey-weighted sandwich.
# The bread is a numerical derivative of the analytic observed-data score; the
# meat is the outer product of weighted individual scores. Nonlinear rates and
# flows use a central-difference delta method. This does not replace a full
# survey-design variance estimator when strata/PSU identifiers are available.

.pack_covariate_params <- function(params, X, model_type) {
  Xt <- .as_transition_design(X)
  active <- attr(Xt$X12, "entry_active") %||% rep(TRUE, ncol(Xt$X12))
  out <- c(params$beta0[active], params$beta1, alpha_logit = qlogis(params$alpha))
  names(out) <- c(paste0("beta0_", colnames(Xt$X12)[active]),
                  paste0("beta1_", colnames(Xt$X12)), "alpha_logit")
  if (model_type == "symmetric") out <- c(out, pi_logit = qlogis(params$pi))
  out
}

.unpack_covariate_params <- function(eta, template, X, model_type) {
  Xt <- .as_transition_design(X)
  active <- attr(Xt$X12, "entry_active") %||% rep(TRUE, ncol(Xt$X12))
  p <- ncol(Xt$X12); p0 <- sum(active)
  out <- template
  out$beta0[] <- 0
  out$beta0[active] <- eta[seq_len(p0)]
  out$beta1 <- eta[p0 + seq_len(p)]
  out$alpha <- plogis(eta[p0 + p + 1L])
  if (model_type == "symmetric") out$pi <- plogis(eta[p0 + p + 2L])
  names(out$beta0) <- names(template$beta0) %||% colnames(Xt$X12)
  names(out$beta1) <- names(template$beta1) %||% colnames(Xt$X12)
  out
}

.collapse_covariate_information_cells <- function(df, X) {
  Xt <- .as_transition_design(X, nrow(df))
  key_df <- data.frame(y1 = df$y1, y2 = df$y2, y3 = df$y3,
                       Xt$X12, Xt$X23, check.names = FALSE)
  key <- do.call(paste, c(key_df, sep = "\r"))
  first <- !duplicated(key)
  group <- match(key, key[first])
  out_X12 <- Xt$X12[first, , drop = FALSE]
  out_X23 <- Xt$X23[first, , drop = FALSE]
  attr(out_X12, "entry_active") <- attr(Xt$X12, "entry_active")
  attr(out_X23, "entry_active") <- attr(Xt$X12, "entry_active")
  list(
    df = data.frame(y1 = df$y1[first], y2 = df$y2[first], y3 = df$y3[first],
                    weight = as.vector(rowsum(df$weight, group, reorder = FALSE))),
    X = list(X12 = out_X12, X23 = out_X23)
  )
}

.individual_covariate_scores <- function(df, X, params, model_type, gamma) {
  Xt <- .as_transition_design(X, nrow(df))
  h <- .hmat_cache()
  active <- attr(Xt$X12, "entry_active") %||% rep(TRUE, ncol(Xt$X12))
  blocks0 <- list(); blocks1 <- list()
  for (tt in seq_len(2L)) {
    hf <- h[, tt]; ht <- h[, tt + 1L]
    risk0 <- as.vector(gamma %*% as.integer(hf == 0L))
    risk1 <- as.vector(gamma %*% as.integer(hf == 1L))
    succ0 <- as.vector(gamma %*% as.integer(hf == 0L & ht == 1L))
    succ1 <- as.vector(gamma %*% as.integer(hf == 1L & ht == 1L))
    eta0 <- as.vector(Xt[[tt]] %*% params$beta0)
    eta1 <- as.vector(Xt[[tt]] %*% params$beta1)
    th0 <- pmin(pmax(pnorm(eta0), 1e-9), 1 - 1e-9)
    th1 <- pmin(pmax(pnorm(eta1), 1e-9), 1 - 1e-9)
    s0 <- dnorm(eta0) * (succ0 / th0 - (risk0 - succ0) / (1 - th0))
    s1 <- dnorm(eta1) * (succ1 / th1 - (risk1 - succ1) / (1 - th1))
    blocks0[[tt]] <- sweep(Xt[[tt]][, active, drop = FALSE], 1L, s0, "*")
    blocks1[[tt]] <- sweep(Xt[[tt]], 1L, s1, "*")
  }
  score <- cbind(blocks0[[1L]] + blocks0[[2L]],
                 blocks1[[1L]] + blocks1[[2L]])
  score <- cbind(score,
                 alpha_logit = as.vector(gamma %*% h[, 1L]) - params$alpha)
  if (model_type == "symmetric") {
    mismatch <- outer(df$y1, h[, 1L], "!=") +
      outer(df$y2, h[, 2L], "!=") + outer(df$y3, h[, 3L], "!=")
    score <- cbind(score, pi_logit = rowSums(gamma * mismatch) - 3 * params$pi)
  }
  score
}

.individual_implied_components <- function(df, X, params, gamma) {
  Xt <- .as_transition_design(X, nrow(df)); h <- .hmat_cache()
  risk0 <- list(as.vector(gamma %*% as.integer(h[, 1L] == 0L)),
                as.vector(gamma %*% as.integer(h[, 2L] == 0L)))
  risk1 <- lapply(risk0, function(z) 1 - z)
  theta0 <- lapply(Xt, function(z) pnorm(as.vector(z %*% params$beta0)))
  exit <- lapply(Xt, function(z) 1 - pnorm(as.vector(z %*% params$beta1)))
  entry_num <- risk0[[1L]] * theta0[[1L]] + risk0[[2L]] * theta0[[2L]]
  exit_num <- risk1[[1L]] * exit[[1L]] + risk1[[2L]] * exit[[2L]]
  out <- list(
    mean_entry_rate = list(a = entry_num, d = risk0[[1L]] + risk0[[2L]]),
    mean_exit_rate = list(a = exit_num, d = risk1[[1L]] + risk1[[2L]]),
    mean_employment_rate = list(a = (risk1[[1L]] + risk1[[2L]]) / 2, d = rep(1, nrow(df))),
    entry_flow = list(a = entry_num / 2, d = rep(1, nrow(df))),
    exit_flow = list(a = exit_num / 2, d = rep(1, nrow(df))),
    total_churn_flow = list(a = (entry_num + exit_num) / 2, d = rep(1, nrow(df)))
  )
  effect_columns <- c(contract_exit_effect = "contracttype_1",
                      informal_exit_effect = "informal_sector")
  for (nm in names(effect_columns)) {
    column <- effect_columns[[nm]]
    j <- match(column, colnames(Xt$X12))
    if (is.na(j)) next
    effect_num <- numeric(nrow(df))
    for (tt in seq_len(2L)) {
      eta_base <- as.vector(Xt[[tt]] %*% params$beta1) -
        Xt[[tt]][, j] * params$beta1[j]
      delta <- (1 - pnorm(eta_base + params$beta1[j])) - (1 - pnorm(eta_base))
      effect_num <- effect_num + risk1[[tt]] * delta
    }
    out[[nm]] <- list(a = effect_num, d = risk1[[1L]] + risk1[[2L]])
  }
  out
}

#' Approximate sandwich and delta-method standard errors
analytical_se_covariates <- function(df, X, fit, model_type,
                                     step_scale = 1e-4) {
  Xt <- .as_transition_design(X, nrow(df))
  template <- fit$params
  coef_names <- colnames(Xt$X12)
  names(template$beta0) <- coef_names
  names(template$beta1) <- coef_names
  eta0 <- .pack_covariate_params(template, Xt, model_type)
  K <- length(eta0); N <- nrow(df); weight_scale <- mean(df$weight)
  collapsed <- .collapse_covariate_information_cells(df, Xt)

  q_names <- c("mean_entry_rate", "mean_exit_rate", "mean_employment_rate",
               "entry_flow", "exit_flow", "total_churn_flow",
               "contract_exit_effect", "informal_exit_effect")
  evaluate <- function(eta) {
    params <- .unpack_covariate_params(eta, template, collapsed$X, model_type)
    estep <- e_step_covariates(collapsed$df, collapsed$X, params, model_type,
                               stationary = FALSE)
    score <- .individual_covariate_scores(collapsed$df, collapsed$X, params,
                                           model_type, estep$gamma)
    total_score <- colSums(score * (collapsed$df$weight / weight_scale))
    implied <- implied_covariates(params, collapsed$X, model_type,
                                  df = collapsed$df, gamma = estep$gamma)
    q <- unlist(implied[q_names], use.names = TRUE)
    list(score = total_score, q = q)
  }

  base <- evaluate(eta0)
  jac_score <- matrix(NA_real_, K, K, dimnames = list(names(eta0), names(eta0)))
  q_keep <- names(base$q)[is.finite(base$q)]
  jac_q <- matrix(NA_real_, length(q_keep), K,
                  dimnames = list(q_keep, names(eta0)))
  for (j in seq_len(K)) {
    hh <- step_scale * (1 + abs(eta0[j]))
    plus <- minus <- eta0
    plus[j] <- plus[j] + hh; minus[j] <- minus[j] - hh
    ep <- evaluate(plus); em <- evaluate(minus)
    jac_score[, j] <- (ep$score - em$score) / (2 * hh)
    jac_q[, j] <- (ep$q[q_keep] - em$q[q_keep]) / (2 * hh)
  }
  bread <- -(jac_score + t(jac_score)) / 2

  gamma_full <- fit$gamma
  if (is.null(gamma_full) || nrow(gamma_full) != N)
    gamma_full <- e_step_covariates(df, Xt, template, model_type,
                                    stationary = FALSE)$gamma
  score_i <- .individual_covariate_scores(df, Xt, template, model_type, gamma_full)
  weighted_score_i <- score_i * (df$weight / weight_scale)
  meat <- crossprod(weighted_score_i) * N / max(N - K, 1L)

  eig <- eigen(bread, symmetric = TRUE)
  floor_value <- max(abs(eig$values)) * 1e-10
  inv_bread <- eig$vectors %*% diag(1 / pmax(eig$values, floor_value), K) %*%
    t(eig$vectors)
  vcov_eta <- inv_bread %*% meat %*% inv_bread
  dimnames(vcov_eta) <- list(names(eta0), names(eta0))

  natural_est <- c(template$beta0[active <- (attr(Xt$X12, "entry_active") %||%
                                                rep(TRUE, ncol(Xt$X12)))],
                   template$beta1, alpha = template$alpha)
  names(natural_est) <- c(paste0("beta0_", coef_names[active]),
                          paste0("beta1_", coef_names), "alpha")
  deriv <- c(rep(1, sum(active) + ncol(Xt$X12)),
             template$alpha * (1 - template$alpha))
  if (model_type == "symmetric") {
    natural_est <- c(natural_est, pi = template$pi)
    deriv <- c(deriv, template$pi * (1 - template$pi))
  }
  vcov_natural <- diag(deriv) %*% vcov_eta %*% diag(deriv)
  param_se <- sqrt(pmax(diag(vcov_natural), 0))
  # Full linearization: direct sampling variation in each weighted empirical
  # functional plus the parameter-estimation component from the sandwich IF.
  components <- .individual_implied_components(df, Xt, template, gamma_full)
  direct_if <- matrix(0, N, length(q_keep), dimnames = list(NULL, q_keep))
  w_norm <- df$weight / weight_scale
  for (nm in q_keep) {
    cc <- components[[nm]]
    denom <- sum(w_norm * cc$d)
    direct_if[, nm] <- w_norm * (cc$a - base$q[[nm]] * cc$d) / denom
  }
  eta_if <- weighted_score_i %*% t(inv_bread)
  parameter_if <- eta_if %*% t(jac_q)
  total_q_if <- direct_if + parameter_if
  q_vcov <- crossprod(total_q_if) * N / max(N - K, 1L)
  q_se <- sqrt(pmax(diag(q_vcov), 0))
  summary <- data.frame(
    quantity = c(names(natural_est), q_keep),
    estimate = c(natural_est, base$q[q_keep]),
    se = c(param_se, q_se), stringsAsFactors = FALSE
  )
  summary <- summary[!duplicated(summary$quantity), , drop = FALSE]
  newton_correction <- as.vector(inv_bread %*% base$score)
  names(newton_correction) <- names(eta0)
  list(summary = summary, vcov_eta = vcov_eta, bread = bread, meat = meat,
       delta_jacobian = jac_q, vcov_quantities = q_vcov,
       diagnostics = list(N = N, K = K, cells = nrow(collapsed$df),
                          min_bread_eigenvalue = min(eig$values),
                          condition_number = max(abs(eig$values)) /
                            max(min(abs(eig$values)), .Machine$double.eps),
                          max_abs_score = max(abs(base$score)),
                          total_score = base$score,
                          newton_correction = newton_correction,
                          max_abs_newton_correction = max(abs(newton_correction))))
}
