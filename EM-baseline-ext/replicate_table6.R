# ==============================================================================
# Replicate and audit Table 6: inconsistency-dependent misclassification
#
# This script estimates the observed-data likelihood directly after collapsing
# identical outcome/inconsistency histories.  It avoids the non-monotone
# stationary GEM update in the legacy implementation and reports analytical
# survey-weighted sandwich standard errors with delta-method transformations.
# ==============================================================================

library(here)
library(dplyr)
library(tidyverse)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

.inc_names <- function(stationary) {
  out <- c("theta0", "theta1")
  if (!stationary) out <- c(out, "alpha")
  c(out, "delta0", "delta_age", "delta_education")
}

.unpack_inc <- function(eta, stationary) {
  eta <- setNames(as.numeric(eta), .inc_names(stationary))
  theta0 <- plogis(eta["theta0"])
  theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha(theta0, theta1) else plogis(eta["alpha"])
  list(
    theta0 = unname(theta0), theta1 = unname(theta1), alpha = unname(alpha),
    delta = unname(eta[c("delta0", "delta_age", "delta_education")])
  )
}

.pack_inc <- function(params, stationary) {
  clamp <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)
  out <- c(theta0 = qlogis(clamp(params$theta0)),
           theta1 = qlogis(clamp(params$theta1)))
  if (!stationary) out <- c(out, alpha = qlogis(clamp(params$alpha)))
  c(out, delta0 = params$delta[1L], delta_age = params$delta[2L],
    delta_education = params$delta[3L])
}

.collapse_inc_cells <- function(df, inc_mat) {
  key_df <- data.frame(y1 = df$y1, y2 = df$y2, y3 = df$y3,
                       inc_mat, check.names = FALSE)
  key <- do.call(paste, c(key_df, sep = "\r"))
  first <- !duplicated(key)
  group <- match(key, key[first])
  w <- as.vector(rowsum(df$weight, group, reorder = FALSE))
  w2 <- as.vector(rowsum(df$weight^2, group, reorder = FALSE))
  n <- as.vector(rowsum(rep(1, nrow(df)), group, reorder = FALSE))
  list(pattern = key_df[first, , drop = FALSE], weight = w, weight_sq = w2,
       count = n, n = nrow(df), weight_sum = sum(df$weight))
}

.inc_cell_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  prior <- prior_over_histories(latent, params$theta1, params$theta0, params$alpha)
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  age <- as.matrix(pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(pattern[, paste0("Y_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(params$delta[1L] + params$delta[2L] * age +
                           params$delta[3L] * edu)
  out <- numeric(nrow(pattern))
  for (h in seq_len(nrow(latent))) {
    emission <- matrix(1, nrow(pattern), 3L)
    for (tt in 1:3)
      emission[, tt] <- ifelse(y[, tt] == latent[h, tt],
                               1 - pi_mat[, tt], pi_mat[, tt])
    out <- out + prior[h] * apply(emission, 1L, prod)
  }
  out
}

.numeric_jacobian_inc <- function(fn, x, rel_step = 1e-5) {
  f0 <- fn(x)
  out <- matrix(NA_real_, length(f0), length(x))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
    out[, j] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  dimnames(out) <- list(names(f0), names(x))
  out
}

.matrix_rank_inc <- function(x) {
  d <- svd(x, nu = 0L, nv = 0L)$d
  sum(d > max(dim(x)) * max(d) * sqrt(.Machine$double.eps))
}

.fit_inc_mle <- function(cells, stationary, starts, verbose = TRUE) {
  fn <- function(eta) {
    p <- pmax(.inc_cell_probabilities(.unpack_inc(eta, stationary), cells$pattern),
              1e-300)
    -sum(cells$weight * log(p)) / cells$weight_sum
  }
  candidates <- lapply(starts, function(start) {
    eta0 <- .pack_inc(start, stationary)
    opt <- tryCatch(optim(eta0, fn, method = "BFGS",
                          control = list(maxit = 4000L, reltol = 1e-12)),
                    error = function(e) NULL)
    if (is.null(opt) || !is.finite(opt$value)) return(NULL)
    list(opt = opt, eta = setNames(opt$par, .inc_names(stationary)),
         loglik = -opt$value * cells$weight_sum)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("All inconsistency-model optimization starts failed")
  best <- candidates[[which.max(vapply(candidates, `[[`, numeric(1L), "loglik"))]]
  eta <- best$eta
  params <- .unpack_inc(eta, stationary)
  grad <- .numeric_gradient(fn, eta)
  hessian <- optimHess(eta, fn)
  hessian <- (hessian + t(hessian)) / 2
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  # Conditional probabilities provide seven independent outcome probabilities
  # for every distinct six-indicator history.  Stack their Jacobians to test
  # whether all local parameter directions change the observable distribution.
  x_names <- c(paste0("Y_age_", 1:3), paste0("Y_edu_", 1:3))
  x_unique <- unique(cells$pattern[, x_names, drop = FALSE])
  outcomes <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
  probability_map <- function(z) {
    blocks <- lapply(seq_len(nrow(x_unique)), function(i) {
      pp <- cbind(as.data.frame(outcomes),
                  x_unique[rep(i, 8L), , drop = FALSE])
      .inc_cell_probabilities(.unpack_inc(z, stationary), pp)[1:7]
    })
    unlist(blocks, use.names = FALSE)
  }
  probability_jacobian <- .numeric_jacobian_inc(probability_map, eta)
  rank <- .matrix_rank_inc(probability_jacobian)

  # Individual-level survey-weighted sandwich covariance, calculated exactly
  # from the collapsed pattern cells using sums of w_i and w_i^2.
  score <- .numeric_jacobian_inc(function(z) {
    log(pmax(.inc_cell_probabilities(.unpack_inc(z, stationary), cells$pattern),
             1e-300))
  }, eta)
  meat <- matrix(0, length(eta), length(eta))
  for (i in seq_len(nrow(score)))
    meat <- meat + cells$weight_sq[i] * tcrossprod(score[i, ])
  meat <- meat / cells$weight_sum^2
  bread_inv <- tryCatch(solve(hessian), error = function(e) MASS::ginv(hessian))
  vcov <- bread_inv %*% meat %*% bread_inv
  dimnames(vcov) <- list(names(eta), names(eta))

  if (verbose) cat(sprintf(
    "Table 6 direct MLE [stationary=%s]: ll=%.3f code=%d max|score|=%.2e minEig=%.2e rank=%d/%d\n",
    stationary, best$loglik, best$opt$convergence, max(abs(grad)), min(eig),
    rank, length(eta)))
  list(params = params, eta = eta, loglik = best$loglik, vcov = vcov,
       stationary = stationary,
       converged = best$opt$convergence == 0L && max(abs(grad)) < 1e-6 && min(eig) > 0,
       identified = rank == length(eta),
       diagnostics = list(code = best$opt$convergence, max_abs_score = max(abs(grad)),
                          min_information_eigenvalue = min(eig), rank = rank,
                          parameter_count = length(eta), starts = length(starts)))
}

.inc_quantities <- function(eta, stationary, cells) {
  p <- .unpack_inc(eta, stationary)
  pi0 <- 0.5 * plogis(p$delta[1L])
  pi_age <- 0.5 * plogis(p$delta[1L] + p$delta[2L])
  pi_edu <- 0.5 * plogis(p$delta[1L] + p$delta[3L])
  pi_both <- 0.5 * plogis(sum(p$delta))
  age <- as.matrix(cells$pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(cells$pattern[, paste0("Y_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(p$delta[1L] + p$delta[2L] * age + p$delta[3L] * edu)
  mean_pi_weighted <- sum(cells$weight * rowSums(pi_mat)) / (3 * cells$weight_sum)
  mean_pi_unweighted <- sum(cells$count * rowSums(pi_mat)) / (3 * cells$n)
  steady <- p$theta0 / (p$theta0 + 1 - p$theta1)
  c(entry_rate = p$theta0, exit_rate = 1 - p$theta1,
    initial_employment = p$alpha, steady_employment = steady,
    stationarity_gap = p$alpha - steady,
    delta0 = p$delta[1L], delta_age = p$delta[2L],
    delta_education = p$delta[3L], pi_base = pi0,
    pi_age = pi_age, age_effect = pi_age - pi0,
    pi_education = pi_edu, education_effect = pi_edu - pi0,
    pi_both = pi_both, mean_pi_unweighted = mean_pi_unweighted,
    mean_pi_survey_weighted = mean_pi_weighted)
}

.inc_inference <- function(fit, cells) {
  q <- function(z) .inc_quantities(z, fit$stationary, cells)
  estimate <- q(fit$eta)
  J <- .numeric_jacobian_inc(q, fit$eta)
  variance <- J %*% fit$vcov %*% t(J)
  data.frame(quantity = names(estimate), estimate = unname(estimate),
             std_error = sqrt(pmax(diag(variance), 0)), row.names = NULL)
}

.fit_custom_inc <- function(cells, eta_names, unpack, probabilities, starts,
                            label) {
  fn <- function(eta) {
    p <- pmax(probabilities(unpack(eta), cells$pattern), 1e-300)
    -sum(cells$weight * log(p)) / cells$weight_sum
  }
  candidates <- lapply(starts, function(eta0) {
    eta0 <- setNames(as.numeric(eta0), eta_names)
    opt <- tryCatch(optim(eta0, fn, method = "BFGS",
                          control = list(maxit = 5000L, reltol = 1e-12)),
                    error = function(e) NULL)
    if (is.null(opt) || !is.finite(opt$value)) return(NULL)
    list(opt = opt, eta = setNames(opt$par, eta_names),
         loglik = -opt$value * cells$weight_sum)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop(label, ": all optimization starts failed")
  best <- candidates[[which.max(vapply(candidates, `[[`, numeric(1L), "loglik"))]]
  eta <- best$eta
  grad <- .numeric_gradient(fn, eta)
  hessian <- optimHess(eta, fn); hessian <- (hessian + t(hessian)) / 2
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  design_names <- setdiff(names(cells$pattern), c("y1", "y2", "y3"))
  design_unique <- unique(cells$pattern[, design_names, drop = FALSE])
  outcomes <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
  probability_map <- function(z) unlist(lapply(seq_len(nrow(design_unique)), function(i) {
    pp <- cbind(as.data.frame(outcomes), design_unique[rep(i, 8L), , drop = FALSE])
    probabilities(unpack(z), pp)[1:7]
  }), use.names = FALSE)
  rank <- .matrix_rank_inc(.numeric_jacobian_inc(probability_map, eta))

  score <- .numeric_jacobian_inc(function(z)
    log(pmax(probabilities(unpack(z), cells$pattern), 1e-300)), eta)
  meat <- matrix(0, length(eta), length(eta))
  for (i in seq_len(nrow(score)))
    meat <- meat + cells$weight_sq[i] * tcrossprod(score[i, ])
  meat <- meat / cells$weight_sum^2
  bread_inv <- tryCatch(solve(hessian), error = function(e) MASS::ginv(hessian))
  vcov <- bread_inv %*% meat %*% bread_inv
  dimnames(vcov) <- list(eta_names, eta_names)
  cat(sprintf("%s: ll=%.3f code=%d max|score|=%.2e minEig=%.2e rank=%d/%d\n",
              label, best$loglik, best$opt$convergence, max(abs(grad)),
              min(eig), rank, length(eta)))
  list(params = unpack(eta), eta = eta, loglik = best$loglik, vcov = vcov,
       converged = best$opt$convergence == 0L && max(abs(grad)) < 1e-6 && min(eig) > 0,
       identified = rank == length(eta),
       diagnostics = list(code = best$opt$convergence, max_abs_score = max(abs(grad)),
                          min_information_eigenvalue = min(eig), rank = rank,
                          parameter_count = length(eta), starts = length(starts)))
}

.group_eta_names <- c("theta0_reliable", "theta1_reliable", "alpha_reliable",
                      "theta0_unreliable", "theta1_unreliable", "alpha_unreliable",
                      "delta0", "delta_age", "delta_education")

.unpack_group <- function(eta) {
  eta <- setNames(as.numeric(eta), .group_eta_names)
  list(theta0 = plogis(eta[c("theta0_reliable", "theta0_unreliable")]),
       theta1 = plogis(eta[c("theta1_reliable", "theta1_unreliable")]),
       alpha = plogis(eta[c("alpha_reliable", "alpha_unreliable")]),
       delta = unname(eta[c("delta0", "delta_age", "delta_education")]))
}

.group_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  age <- as.matrix(pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(pattern[, paste0("Y_edu_", 1:3)])
  unreliable <- as.integer(rowSums(cbind(age, edu)) > 0) + 1L
  pi_mat <- 0.5 * plogis(params$delta[1L] + params$delta[2L] * age +
                           params$delta[3L] * edu)
  out <- numeric(nrow(pattern))
  for (g in 1:2) {
    idx <- which(unreliable == g)
    if (!length(idx)) next
    prior <- prior_over_histories(latent, params$theta1[g], params$theta0[g],
                                  params$alpha[g])
    for (h in seq_len(nrow(latent))) {
      emission <- matrix(1, length(idx), 3L)
      for (tt in 1:3)
        emission[, tt] <- ifelse(y[idx, tt] == latent[h, tt],
                                 1 - pi_mat[idx, tt], pi_mat[idx, tt])
      out[idx] <- out[idx] + prior[h] * apply(emission, 1L, prod)
    }
  }
  out
}

.group_quantities <- function(eta, cells) {
  p <- .unpack_group(eta)
  pi0 <- 0.5 * plogis(p$delta[1L])
  pi_age <- 0.5 * plogis(p$delta[1L] + p$delta[2L])
  pi_edu <- 0.5 * plogis(p$delta[1L] + p$delta[3L])
  age <- as.matrix(cells$pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(cells$pattern[, paste0("Y_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(p$delta[1L] + p$delta[2L] * age + p$delta[3L] * edu)
  c(entry_reliable = unname(p$theta0[1L]),
    exit_reliable = unname(1 - p$theta1[1L]),
    initial_reliable = unname(p$alpha[1L]),
    entry_unreliable = unname(p$theta0[2L]),
    exit_unreliable = unname(1 - p$theta1[2L]),
    initial_unreliable = unname(p$alpha[2L]),
    delta0 = p$delta[1L], delta_age = p$delta[2L], delta_education = p$delta[3L],
    pi_base = pi0, pi_age = pi_age, age_effect = pi_age - pi0,
    pi_education = pi_edu, education_effect = pi_edu - pi0,
    mean_pi_survey_weighted = sum(cells$weight * rowSums(pi_mat)) /
      (3 * cells$weight_sum))
}

.extent_eta_names <- c("theta0", "theta1", "alpha", "delta0", "delta_age",
                       "delta_education", "delta_age_severe", "delta_education_severe")

.unpack_extent <- function(eta) {
  eta <- setNames(as.numeric(eta), .extent_eta_names)
  list(theta0 = plogis(eta["theta0"]), theta1 = plogis(eta["theta1"]),
       alpha = plogis(eta["alpha"]),
       delta = unname(eta[c("delta0", "delta_age", "delta_education",
                            "delta_age_severe", "delta_education_severe")]))
}

.extent_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  prior <- prior_over_histories(latent, params$theta1, params$theta0, params$alpha)
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  age <- as.matrix(pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(pattern[, paste0("Y_edu_", 1:3)])
  age_extent <- as.matrix(pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(pattern[, paste0("severe_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(params$delta[1L] + params$delta[2L] * age +
    params$delta[3L] * edu + params$delta[4L] * age_extent +
    params$delta[5L] * edu_extent)
  out <- numeric(nrow(pattern))
  for (h in seq_len(nrow(latent))) {
    emission <- matrix(1, nrow(pattern), 3L)
    for (tt in 1:3)
      emission[, tt] <- ifelse(y[, tt] == latent[h, tt],
                               1 - pi_mat[, tt], pi_mat[, tt])
    out <- out + prior[h] * apply(emission, 1L, prod)
  }
  out
}

.extent_quantities <- function(eta, cells) {
  p <- .unpack_extent(eta); d <- p$delta
  scenario_pi <- function(age, edu, age_extent, edu_extent)
    0.5 * plogis(d[1L] + d[2L] * age + d[3L] * edu +
                   d[4L] * age_extent + d[5L] * edu_extent)
  age <- as.matrix(cells$pattern[, paste0("Y_age_", 1:3)])
  edu <- as.matrix(cells$pattern[, paste0("Y_edu_", 1:3)])
  age_extent <- as.matrix(cells$pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(cells$pattern[, paste0("severe_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(d[1L] + d[2L] * age + d[3L] * edu +
                          d[4L] * age_extent + d[5L] * edu_extent)
  c(entry_rate = unname(p$theta0), exit_rate = unname(1 - p$theta1),
    initial_employment = unname(p$alpha),
    delta0 = d[1L], delta_age = d[2L], delta_education = d[3L],
    delta_age_severe = d[4L], delta_education_severe = d[5L],
    pi_base = scenario_pi(0, 0, 0, 0),
    pi_age_mild = scenario_pi(1, 0, 0, 0),
    pi_age_severe = scenario_pi(1, 0, 1, 0),
    age_severity_effect = scenario_pi(1, 0, 1, 0) - scenario_pi(1, 0, 0, 0),
    pi_education_mild = scenario_pi(0, 1, 0, 0),
    pi_education_severe = scenario_pi(0, 1, 0, 1),
    education_severity_effect = scenario_pi(0, 1, 0, 1) - scenario_pi(0, 1, 0, 0),
    mean_pi_survey_weighted = sum(cells$weight * rowSums(pi_mat)) /
      (3 * cells$weight_sum))
}

.custom_inference <- function(fit, quantities, cells) {
  estimate <- quantities(fit$eta, cells)
  J <- .numeric_jacobian_inc(function(z) quantities(z, cells), fit$eta)
  variance <- J %*% fit$vcov %*% t(J)
  data.frame(quantity = names(estimate), estimate = unname(estimate),
             std_error = sqrt(pmax(diag(variance), 0)), row.names = NULL)
}

# Match the legacy Table 6 analysis sample.
keep <- complete.cases(df_qlfs[, c("y1", "y2", "y3", "weight",
                                   "age1", "age2", "age3",
                                   "educ1", "educ2", "educ3")]) &
  df_qlfs$weight > 0
df <- as.data.frame(df_qlfs[keep, , drop = FALSE])
df$y1 <- as.integer(df$y1); df$y2 <- as.integer(df$y2); df$y3 <- as.integer(df$y3)
df$weight <- as.numeric(df$weight)
inc_df <- compute_inconsistency_extent(df)
inc_names <- c(paste0("Y_age_", 1:3), paste0("Y_edu_", 1:3))
inc_mat <- as.matrix(inc_df[, inc_names])
cells <- .collapse_inc_cells(df, inc_mat)
extent_names <- c(paste0("extent_age_", 1:3), paste0("extent_edu_", 1:3))
severe_names <- c(paste0("severe_age_", 1:3), paste0("severe_edu_", 1:3))
inc_df[severe_names] <- lapply(inc_df[extent_names], function(x) as.integer(x >= 2))
cells_extent <- .collapse_inc_cells(
  df, as.matrix(inc_df[, c(inc_names, severe_names)])
)

cat(sprintf("Table 6 sample: N=%s; %d collapsed likelihood cells\n",
            format(nrow(df), big.mark = ","), nrow(cells$pattern)))
cat(sprintf("Severity robustness: %d collapsed likelihood cells\n",
            nrow(cells_extent$pattern)))
flag_summary <- data.frame(
  indicator = c(inc_names, "any_age", "any_education", "any_inconsistency"),
  unweighted_percent = 100 * c(colMeans(inc_mat),
    mean(rowSums(inc_mat[, 1:3, drop = FALSE]) > 0),
    mean(rowSums(inc_mat[, 4:6, drop = FALSE]) > 0),
    mean(rowSums(inc_mat) > 0)),
  survey_weighted_percent = 100 * c(
    colSums(inc_mat * df$weight) / sum(df$weight),
    weighted.mean(rowSums(inc_mat[, 1:3, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(inc_mat[, 4:6, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(inc_mat) > 0, df$weight))
)
cat("\nInconsistency prevalence (%):\n")
print(flag_summary, row.names = FALSE, digits = 4)
severity_summary <- data.frame(
  variable = c("age", "education"),
  severe_share_of_flagged_percent = 100 * c(
    sum(df$weight * rowSums(as.matrix(inc_df[, paste0("severe_age_", 1:3)]))) /
      sum(df$weight * rowSums(inc_mat[, 1:3, drop = FALSE])),
    sum(df$weight * rowSums(as.matrix(inc_df[, paste0("severe_edu_", 1:3)]))) /
      sum(df$weight * rowSums(inc_mat[, 4:6, drop = FALSE]))
  )
)
cat("\nShare of attributed flags at least two units beyond the admissible range (%):\n")
print(severity_summary, row.names = FALSE, digits = 4)

# Obtain nested homogeneous estimates for stable starts, then add dispersed
# slope starts to guard against local optima.
set.seed(62026L)
make_starts <- function(stationary) {
  base <- fit_baseline_mle(df, "symmetric", stationary,
                           starts = list(init_params("symmetric", stationary)),
                           verbose = 0L)$params
  center <- list(theta0 = base$theta0, theta1 = base$theta1, alpha = base$alpha,
                 delta = c(qlogis(2 * base$pi), 0, 0))
  fixed <- list(center,
    within(center, delta <- c(delta[1L], 2.5, 1.0)),
    within(center, delta <- c(delta[1L], 1.0, 2.5)),
    within(center, delta <- c(delta[1L], -1.0, 1.0)))
  random <- lapply(1:8, function(i) {
    z <- center
    z$theta0 <- plogis(qlogis(z$theta0) + rnorm(1, 0, .35))
    z$theta1 <- plogis(qlogis(z$theta1) + rnorm(1, 0, .35))
    if (!stationary) z$alpha <- plogis(qlogis(z$alpha) + rnorm(1, 0, .25))
    z$delta <- z$delta + c(rnorm(1, 0, .35), rnorm(2, 0, 1.25))
    z
  })
  c(fixed, random)
}

fit_stat <- .fit_inc_mle(cells, TRUE, make_starts(TRUE))
fit_free <- .fit_inc_mle(cells, FALSE, make_starts(FALSE))
inf_stat <- .inc_inference(fit_stat, cells)
inf_free <- .inc_inference(fit_free, cells)

# Robustness 1: allow the true transition process and initial employment
# probability to differ between records with and without any inconsistency.
group_center <- c(
  theta0_reliable = unname(fit_free$eta["theta0"]),
  theta1_reliable = unname(fit_free$eta["theta1"]),
  alpha_reliable = unname(fit_free$eta["alpha"]),
  theta0_unreliable = unname(fit_free$eta["theta0"]),
  theta1_unreliable = unname(fit_free$eta["theta1"]),
  alpha_unreliable = unname(fit_free$eta["alpha"]),
  delta0 = unname(fit_free$eta["delta0"]),
  delta_age = unname(fit_free$eta["delta_age"]),
  delta_education = unname(fit_free$eta["delta_education"])
)
group_starts <- c(list(group_center), lapply(1:11, function(i)
  group_center + rnorm(length(group_center), 0,
                       c(rep(.35, 6), .25, .75, .75))))
fit_group <- .fit_custom_inc(cells, .group_eta_names, .unpack_group,
                             .group_probabilities, group_starts,
                             "Reliability-group transition robustness")
inf_group <- .custom_inference(fit_group, .group_quantities, cells)

# Robustness 2: add the distance from the admissible [0,1] age/education
# change to test whether more severe inconsistencies imply larger error rates.
extent_center <- c(fit_free$eta, delta_age_severe = 0,
                   delta_education_severe = 0)
extent_center <- extent_center[.extent_eta_names]
extent_starts <- c(list(extent_center), lapply(1:11, function(i)
  extent_center + rnorm(length(extent_center), 0,
                        c(rep(.30, 3), .25, .75, .75, .25, .25))))
fit_extent <- .fit_custom_inc(cells_extent, .extent_eta_names, .unpack_extent,
                              .extent_probabilities, extent_starts,
                              "Inconsistency-extent robustness")
inf_extent <- .custom_inference(fit_extent, .extent_quantities, cells_extent)

table6 <- merge(inf_stat, inf_free, by = "quantity", all = TRUE,
                suffixes = c("_stationary", "_free"), sort = FALSE)
order_q <- .inc_quantities(fit_stat$eta, TRUE, cells)
table6 <- table6[match(names(order_q), table6$quantity), ]
table6$estimate_stationary[table6$quantity == "initial_employment"] <- NA_real_
table6$std_error_stationary[table6$quantity == "initial_employment"] <- NA_real_

cat("\nReplicated Table 6 (probabilities are proportions; analytical SE):\n")
print(table6, row.names = FALSE, digits = 6)
cat(sprintf("\nLog likelihood: stationary %.3f; free alpha %.3f\n",
            fit_stat$loglik, fit_free$loglik))
cat(sprintf("LR test of stationarity: LR=%.3f, p=%.4g\n",
            2 * (fit_free$loglik - fit_stat$loglik),
            pchisq(2 * (fit_free$loglik - fit_stat$loglik), 1L, lower.tail = FALSE)))
gap_row <- inf_free[inf_free$quantity == "stationarity_gap", ]
cat(sprintf("Robust Wald test of stationarity: gap=%.4f, SE=%.4f, z=%.2f, p=%.4g\n",
            gap_row$estimate, gap_row$std_error,
            gap_row$estimate / gap_row$std_error,
            2 * pnorm(-abs(gap_row$estimate / gap_row$std_error))))
cat(sprintf("Identified/converged: stationary %s/%s; free %s/%s\n",
            fit_stat$identified, fit_stat$converged,
            fit_free$identified, fit_free$converged))
cat("\nRobustness: group-specific true transition rates:\n")
print(inf_group, row.names = FALSE, digits = 6)
cat(sprintf("Identified/converged: %s/%s\n", fit_group$identified, fit_group$converged))
cat("\nRobustness: inconsistency extent:\n")
print(inf_extent, row.names = FALSE, digits = 6)
cat(sprintf("Identified/converged: %s/%s\n", fit_extent$identified, fit_extent$converged))

out_dir <- here::here("EM-baseline-ext", "output", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(table6, file.path(out_dir, "table6_inconsistency_audit.csv"), row.names = FALSE)
write.csv(flag_summary, file.path(out_dir, "table6_inconsistency_prevalence.csv"), row.names = FALSE)
write.csv(severity_summary, file.path(out_dir, "table6_inconsistency_severity_prevalence.csv"), row.names = FALSE)
write.csv(inf_group, file.path(out_dir, "table6_reliability_group_robustness.csv"), row.names = FALSE)
write.csv(inf_extent, file.path(out_dir, "table6_inconsistency_extent.csv"), row.names = FALSE)
saveRDS(list(stationary = fit_stat, free = fit_free, table = table6,
             reliability_group = fit_group, reliability_group_table = inf_group,
             extent = fit_extent, extent_table = inf_extent,
             prevalence = flag_summary, severity_prevalence = severity_summary),
        file.path(out_dir, "fit_table6_inconsistency_audit.rds"))
