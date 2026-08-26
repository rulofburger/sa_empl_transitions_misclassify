# Four-wave, one-type AR(2) model with observed tenure and non-employment duration.

prepare_ar2_duration_4w <- function(panel_path = "data/raw/df_qlfs_A.rds",
                                    collapse = TRUE) {
  panel <- readRDS(panel_path)
  required <- c(paste0("employed", 1:4), paste0("age", 1:3),
                paste0("educ", 1:3), paste0("weight", 1:4),
                paste0("tenure", 1:3), paste0("timegap", 1:3),
                paste0("neverworked", 1:3))
  missing <- setdiff(required, names(panel))
  if (length(missing)) stop("Four-wave panel is missing: ", paste(missing, collapse = ", "))
  keep <- panel$age1 > 17 & panel$age1 < 56 &
    complete.cases(panel[paste0("employed", 1:4)])
  panel <- panel[keep, required, drop = FALSE]
  for (nm in names(panel)) panel[[nm]] <- as.numeric(unclass(panel[[nm]]))
  y <- as.matrix(panel[paste0("employed", 1:4)])
  weight <- with(panel, (weight1 * weight2 * weight3 * weight4) ^ .25)
  if (any(!is.finite(weight) | weight <= 0)) stop("Invalid survey weights")
  weight <- nrow(panel) * weight / sum(weight)

  # This is the same coding used in the three-wave covariate model and the
  # duration-expanded four-wave FMM. Durations are observed-state covariates:
  # tenure is active for reported employment and timegap for reported
  # non-employment. Missing indicators retain all complete employment histories.
  timegap_months <- c(`1` = 1.5, `2` = 4.5, `3` = 7.5, `4` = 10.5,
                      `5` = 24, `6` = 48, `7` = 90)
  tenure <- timegap <- never <- tenure_missing <- timegap_missing <- vector("list", 3)
  for (tt in 1:3) {
    reported <- y[, tt]
    ten <- panel[[paste0("tenure", tt)]]
    ten[ten < 0] <- NA_real_
    ten[reported == 0] <- 0
    gap_code <- panel[[paste0("timegap", tt)]]
    gap <- unname(timegap_months[as.character(gap_code)])
    gap[gap_code == 0] <- 0
    gap[reported == 1] <- 0
    nw <- panel[[paste0("neverworked", tt)]]
    use_proxy <- reported == 0 & nw == 1
    gap[use_proxy] <- pmax(panel[[paste0("age", tt)]][use_proxy] -
                            panel[[paste0("educ", tt)]][use_proxy] - 6, 0)
    tenure_missing[[tt]] <- as.integer(reported == 1 & is.na(ten))
    timegap_missing[[tt]] <- as.integer(reported == 0 & is.na(gap))
    tenure[[tt]] <- ifelse(reported == 1, ifelse(is.na(ten), 0, ten / 12), 0)
    timegap[[tt]] <- ifelse(reported == 0, ifelse(is.na(gap), 0, gap / 12), 0)
    never[[tt]] <- ifelse(reported == 0, ifelse(is.na(nw), 0, nw), 0)
  }
  standardise_across_waves <- function(values) {
    raw <- log1p(unlist(values, use.names = FALSE))
    center <- mean(raw); scale <- sd(raw)
    if (!is.finite(scale) || scale < 1e-10) stop("Duration covariate is constant")
    list(values = split((raw - center) / scale, rep(1:3, each = nrow(panel))),
         center = center, scale = scale)
  }
  tenure_z <- standardise_across_waves(tenure)
  timegap_z <- standardise_across_waves(timegap)
  # AR(2) transitions are 2->3 and 3->4, so only origin waves 2 and 3 enter.
  X_entry <- lapply(2:3, function(tt) cbind(
    intercept = 1, lag2 = 0,
    log_time_since_work = timegap_z$values[[tt]],
    never_worked = never[[tt]], timegap_missing = timegap_missing[[tt]]))
  X_persistence <- lapply(2:3, function(tt) cbind(
    intercept = 1, lag2 = 0,
    log_tenure = tenure_z$values[[tt]], tenure_missing = tenure_missing[[tt]]))
  out <- list(y = y, weight = weight, weight_sq = weight^2,
       X_entry = X_entry, X_persistence = X_persistence,
       scaling = list(log_tenure = tenure_z[c("center", "scale")],
                      log_time_since_work = timegap_z[c("center", "scale")]),
       n_original = nrow(panel))
  if (!collapse) return(out)
  key_matrix <- cbind(y, do.call(cbind, X_entry), do.call(cbind, X_persistence))
  key <- do.call(paste, c(as.data.frame(key_matrix), sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  collapse_sum <- function(x) as.vector(rowsum(x, group, reorder = FALSE))
  out$y <- out$y[first, , drop = FALSE]
  out$X_entry <- lapply(out$X_entry, function(x) x[first, , drop = FALSE])
  out$X_persistence <- lapply(out$X_persistence, function(x) x[first, , drop = FALSE])
  out$weight <- collapse_sum(out$weight)
  out$weight_sq <- collapse_sum(out$weight_sq)
  out$n_cells <- length(out$weight)
  out
}

.ar2_duration_unpack <- function(z, data) {
  pa <- 3L; p0 <- ncol(data$X_entry[[1]]); p1 <- ncol(data$X_persistence[[1]])
  alpha_raw <- c(1, exp(z[seq_len(pa)] - max(c(0, z[seq_len(pa)]))))
  # Recompute the baseline element on the same stabilised scale.
  shift <- max(c(0, z[seq_len(pa)]))
  alpha_raw <- c(exp(-shift), exp(z[seq_len(pa)] - shift))
  list(alpha = setNames(alpha_raw / sum(alpha_raw), c("00", "10", "01", "11")),
       beta0 = setNames(z[pa + seq_len(p0)], colnames(data$X_entry[[1]])),
       beta1 = setNames(z[pa + p0 + seq_len(p1)], colnames(data$X_persistence[[1]])),
       pi = .5 * plogis(z[pa + p0 + p1 + 1L]))
}

.ar2_duration_log_components <- function(z, data, retain_scores = FALSE) {
  p <- .ar2_duration_unpack(z, data)
  h <- latent_histories_ar2(); n <- nrow(data$y); H <- nrow(h)
  logc <- matrix(0, n, H)
  pair <- paste0(h[, 1], h[, 2])
  logc <- sweep(logc, 2, log(p$alpha[pair]), "+")
  for (tt in 1:2) {
    for (j in seq_len(H)) {
      X <- if (h[j, tt + 1] == 0) data$X_entry[[tt]] else data$X_persistence[[tt]]
      beta <- if (h[j, tt + 1] == 0) p$beta0 else p$beta1
      X[, "lag2"] <- h[j, tt]
      q <- pmin(pmax(pnorm(as.vector(X %*% beta)), 1e-10), 1 - 1e-10)
      logc[, j] <- logc[, j] + if (h[j, tt + 2] == 1) log(q) else log1p(-q)
    }
  }
  pi <- pmin(pmax(p$pi, 1e-10), .5 - 1e-10)
  for (j in seq_len(H)) {
    mismatch <- rowSums(data$y != matrix(h[j, ], n, 4, byrow = TRUE))
    logc[, j] <- logc[, j] + mismatch * log(pi) + (4 - mismatch) * log1p(-pi)
  }
  mx <- apply(logc, 1, max); denom <- rowSums(exp(logc - mx)); loglik_i <- mx + log(denom)
  if (!retain_scores) return(list(loglik_i = loglik_i, params = p))
  post <- exp(logc - mx) / denom
  K <- length(z); scores <- matrix(0, n, K)
  alpha_ind <- cbind(h[,1] == 1 & h[,2] == 0,
                     h[,1] == 0 & h[,2] == 1,
                     h[,1] == 1 & h[,2] == 1) * 1
  scores[, 1:3] <- post %*% alpha_ind - matrix(p$alpha[c("10","01","11")], n, 3, byrow=TRUE)
  offset0 <- 3L; offset1 <- 3L + length(p$beta0)
  for (tt in 1:2) for (j in seq_len(H)) {
    origin <- h[j, tt + 1]; target <- h[j, tt + 2]
    X <- if (origin == 0) data$X_entry[[tt]] else data$X_persistence[[tt]]
    beta <- if (origin == 0) p$beta0 else p$beta1
    X[, "lag2"] <- h[j, tt]
    eta <- as.vector(X %*% beta); q <- pmin(pmax(pnorm(eta), 1e-10), 1 - 1e-10)
    scalar <- dnorm(eta) * (target - q) / (q * (1 - q)) * post[, j]
    cols <- if (origin == 0) offset0 + seq_along(beta) else offset1 + seq_along(beta)
    scores[, cols] <- scores[, cols, drop=FALSE] + X * scalar
  }
  s <- 2 * pi; dpi <- .5 * s * (1 - s)
  mismatch_expectation <- numeric(n)
  for (j in seq_len(H))
    mismatch_expectation <- mismatch_expectation + post[,j] *
      rowSums(data$y != matrix(h[j,], n, 4, byrow=TRUE))
  scores[, K] <- dpi * (mismatch_expectation / pi - (4 - mismatch_expectation) / (1 - pi))
  list(loglik_i = loglik_i, params = p, scores = scores, posterior = post)
}

fit_ar2_duration_4w <- function(data, start = NULL, maxit = 1000L, reltol = 1e-10) {
  p0 <- ncol(data$X_entry[[1]]); p1 <- ncol(data$X_persistence[[1]])
  if (is.null(start)) {
    start <- c(0, 0, 0, qnorm(.04), rep(0, p0 - 1),
               qnorm(.96), rep(0, p1 - 1), qlogis(.04))
  }
  weight_total <- sum(data$weight)
  objective <- function(z) -sum(data$weight * .ar2_duration_log_components(z, data)$loglik_i) / weight_total
  gradient <- function(z) {
    x <- .ar2_duration_log_components(z, data, TRUE)
    -colSums(x$scores * data$weight) / weight_total
  }
  opt <- optim(start, objective, gradient, method = "BFGS",
               control = list(maxit = maxit, reltol = reltol))
  detail <- .ar2_duration_log_components(opt$par, data, TRUE)
  score <- colSums(detail$scores * data$weight) / weight_total
  list(params = detail$params, eta = opt$par,
       loglik = sum(data$weight * detail$loglik_i), converged = opt$convergence == 0,
       iterations = unname(opt$counts[["function"]]), optimizer_code = opt$convergence,
       max_abs_score = max(abs(score)), scores = detail$scores,
       posterior = detail$posterior)
}

analytical_se_ar2_duration_4w <- function(data, fit) {
  wt <- sum(data$weight)
  objective <- function(z) -sum(data$weight * .ar2_duration_log_components(z, data)$loglik_i) / wt
  gradient <- function(z) -colSums(.ar2_duration_log_components(z, data, TRUE)$scores * data$weight) / wt
  bread <- optimHess(fit$eta, objective, gradient); bread <- (bread + t(bread)) / 2
  eig <- eigen(bread, symmetric = TRUE)
  if (min(eig$values) <= max(eig$values) * 1e-10)
    stop("AR(2) duration information matrix is not positive definite")
  meat <- crossprod(fit$scores * sqrt(data$weight_sq)) / wt^2
  inv <- solve(bread); vcov <- inv %*% meat %*% inv
  K <- length(fit$eta); n <- data$n_original
  if (n > K) vcov <- vcov * n / (n - K)
  se_eta <- sqrt(pmax(diag(vcov), 0))
  p <- fit$params
  names_eta <- c(paste0("alpha_logit_", c("10","01","11")),
                 paste0("entry_", names(p$beta0)),
                 paste0("persistence_", names(p$beta1)), "misclassification_logit")
  names(se_eta) <- names_eta
  # Delta-method SEs for probabilities at reference covariate values.
  reported <- function(z) {
    pp <- .ar2_duration_unpack(z, data)
    risk_mean <- function(equation, lag2) {
      values <- weights <- numeric()
      for (tt in 1:2) {
        origin <- data$y[, tt + 1]
        use <- if (equation == "entry") origin == 0 else origin == 1
        X <- if (equation == "entry") data$X_entry[[tt]] else data$X_persistence[[tt]]
        beta <- if (equation == "entry") pp$beta0 else pp$beta1
        X[, "lag2"] <- lag2
        probability <- pnorm(as.vector(X %*% beta))
        if (equation == "persistence") probability <- 1 - probability
        values <- c(values, probability[use]); weights <- c(weights, data$weight[use])
      }
      weighted.mean(values, weights)
    }
    c(reference_entry_00 = unname(pnorm(pp$beta0["intercept"])),
      reference_entry_10 = unname(pnorm(pp$beta0["intercept"] + pp$beta0["lag2"])),
      reference_exit_01 = unname(1 - pnorm(pp$beta1["intercept"])),
      reference_exit_11 = unname(1 - pnorm(pp$beta1["intercept"] + pp$beta1["lag2"])),
      risk_weighted_entry_00 = risk_mean("entry", 0),
      risk_weighted_entry_10 = risk_mean("entry", 1),
      risk_weighted_exit_01 = risk_mean("persistence", 0),
      risk_weighted_exit_11 = risk_mean("persistence", 1),
      pi = unname(pp$pi))
  }
  jacobian <- .ar2_numeric_jacobian(reported, fit$eta)
  prob_se <- sqrt(pmax(diag(jacobian %*% vcov %*% t(jacobian)), 0))
  list(coefficient_summary = data.frame(term = names_eta, estimate = fit$eta, se = se_eta),
       probability_summary = data.frame(quantity = names(reported(fit$eta)),
         estimate = unname(reported(fit$eta)), se = unname(prob_se)),
       vcov_eta = vcov, bread = bread, meat = meat,
       diagnostics = list(information_rank = sum(eig$values > max(eig$values) * 1e-10),
         information_min_eigenvalue = min(eig$values),
         information_condition = max(eig$values) / min(eig$values)))
}
