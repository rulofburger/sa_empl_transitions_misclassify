# ==============================================================================
# Replicate and audit current Table 3 (legacy script/output names retain table6)
#
# This script estimates the observed-data likelihood directly after collapsing
# identical outcome/inconsistency histories.  It avoids the non-monotone
# stationary GEM update in the legacy implementation and reports analytical
# survey-weighted sandwich standard errors with delta-method transformations.
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

.base_delta_names <- c("delta0", "delta_age", "delta_education",
                       "delta_race", "delta_gender",
                       "delta_two_inconsistencies",
                       "delta_three_inconsistencies",
                       "delta_four_inconsistencies")

.inc_names <- function(stationary) {
  out <- c("theta0", "theta1")
  if (!stationary) out <- c(out, "alpha")
  c(out, .base_delta_names)
}

.unpack_inc <- function(eta, stationary) {
  eta <- setNames(as.numeric(eta), .inc_names(stationary))
  theta0 <- plogis(eta["theta0"])
  theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha(theta0, theta1) else plogis(eta["alpha"])
  list(
    theta0 = unname(theta0), theta1 = unname(theta1), alpha = unname(alpha),
    delta = eta[.base_delta_names]
  )
}

.pack_inc <- function(params, stationary) {
  clamp <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)
  out <- c(theta0 = qlogis(clamp(params$theta0)),
           theta1 = qlogis(clamp(params$theta1)))
  if (!stationary) out <- c(out, alpha = qlogis(clamp(params$alpha)))
  c(out, setNames(unname(params$delta[.base_delta_names]), .base_delta_names))
}

.base_error_index <- function(delta, pattern) {
  delta <- setNames(as.numeric(delta), .base_delta_names)
  delta["delta0"] +
    delta["delta_age"] * as.matrix(pattern[, paste0("Y_age_", 1:3)]) +
    delta["delta_education"] * as.matrix(pattern[, paste0("Y_edu_", 1:3)]) +
    delta["delta_race"] * as.matrix(pattern[, paste0("Y_race_", 1:3)]) +
    delta["delta_gender"] * as.matrix(pattern[, paste0("Y_gender_", 1:3)]) +
    delta["delta_two_inconsistencies"] *
      as.matrix(pattern[, paste0("Y_exactly_2_", 1:3)]) +
    delta["delta_three_inconsistencies"] *
      as.matrix(pattern[, paste0("Y_exactly_3_", 1:3)]) +
    delta["delta_four_inconsistencies"] *
      as.matrix(pattern[, paste0("Y_exactly_4_", 1:3)])
}

.base_error_design <- function(pattern) {
  lapply(1:3, function(tt) {
    out <- cbind(
      delta0 = 1,
      delta_age = pattern[[paste0("Y_age_", tt)]],
      delta_education = pattern[[paste0("Y_edu_", tt)]],
      delta_race = pattern[[paste0("Y_race_", tt)]],
      delta_gender = pattern[[paste0("Y_gender_", tt)]],
      delta_two_inconsistencies = pattern[[paste0("Y_exactly_2_", tt)]],
      delta_three_inconsistencies = pattern[[paste0("Y_exactly_3_", tt)]],
      delta_four_inconsistencies = pattern[[paste0("Y_exactly_4_", tt)]])
    storage.mode(out) <- "double"
    out
  })
}

.inc_cell_detail <- function(params, pattern, stationary, design = NULL) {
  latent <- latent_histories(); H <- nrow(latent); n <- nrow(pattern)
  prior <- prior_over_histories(latent, params$theta1, params$theta0,
                                params$alpha)
  joint <- matrix(prior, n, H, byrow = TRUE)
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  if (is.null(design)) design <- .base_error_design(pattern)
  mismatch_masks <- vector("list", 3L)
  pi_wave <- vector("list", 3L)
  for (tt in 1:3) {
    pi <- pmin(pmax(0.5 * plogis(as.vector(
      design[[tt]] %*% params$delta)), 1e-9), 0.5 - 1e-9)
    mismatch <- outer(y[, tt], latent[, tt], "!=")
    emission <- matrix(1 - pi, n, H)
    emission[mismatch] <- matrix(pi, n, H)[mismatch]
    joint <- joint * emission
    mismatch_masks[[tt]] <- mismatch
    pi_wave[[tt]] <- pi
  }
  probability <- rowSums(joint)
  posterior <- joint / probability

  score0 <- score1 <- numeric(n)
  for (tt in 1:2) {
    hf <- latent[, tt]; ht <- latent[, tt + 1L]
    risk0 <- as.vector(posterior %*% as.integer(hf == 0L))
    risk1 <- as.vector(posterior %*% as.integer(hf == 1L))
    succ0 <- as.vector(posterior %*% as.integer(hf == 0L & ht == 1L))
    succ1 <- as.vector(posterior %*% as.integer(hf == 1L & ht == 1L))
    score0 <- score0 + succ0 - params$theta0 * risk0
    score1 <- score1 + succ1 - params$theta1 * risk1
  }
  initial_residual <- as.vector(posterior %*% latent[, 1L]) - params$alpha
  if (stationary) {
    score0 <- score0 + initial_residual * (1 - params$theta0)
    score1 <- score1 + initial_residual * params$theta1
    score <- cbind(theta0 = score0, theta1 = score1)
  } else {
    score <- cbind(theta0 = score0, theta1 = score1,
                   alpha = initial_residual)
  }
  delta_score <- matrix(0, n, length(params$delta),
                        dimnames = list(NULL, names(params$delta)))
  for (tt in 1:3) {
    expected_mismatch <- rowSums(posterior * mismatch_masks[[tt]])
    pi <- pi_wave[[tt]]
    scalar <- (expected_mismatch - pi) * (1 - 2 * pi) / (1 - pi)
    delta_score <- delta_score + design[[tt]] * scalar
  }
  list(probability = probability, scores = cbind(score, delta_score))
}

.conditional_pi_by_count <- function(pi_mat, cells) {
  components <- c("age", "edu", "race", "gender")
  count_mat <- Reduce(`+`, lapply(components, function(x)
    as.matrix(cells$pattern[, paste0("Y_", x, "_", 1:3)])))
  out <- vapply(0:4, function(kk) {
    selected <- count_mat == kk
    denominator <- sum(cells$weight * rowSums(selected))
    if (denominator <= 0) return(NA_real_)
    sum(cells$weight * rowSums(pi_mat * selected)) / denominator
  }, numeric(1L))
  names(out) <- paste0("pi_count_", 0:4)
  out
}

.pi_distribution_summary <- function(pi_mat, cells) {
  values <- as.vector(pi_mat)
  weights <- rep(cells$weight, times = ncol(pi_mat))
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]; weights <- weights[keep]
  ord <- order(values)
  values <- values[ord]; weights <- weights[ord]
  median_index <- which(cumsum(weights) >= sum(weights) / 2)[1L]
  c(min_pi_survey_weighted = values[1L],
    median_pi_survey_weighted = values[median_index],
    max_pi_survey_weighted = values[length(values)])
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

.inc_cell_probabilities <- function(params, pattern, stationary) {
  .inc_cell_detail(params, pattern, stationary)$probability
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
    p <- pmax(.inc_cell_detail(.unpack_inc(eta, stationary), cells$pattern,
                              stationary)$probability, 1e-300)
    -sum(cells$weight * log(p)) / cells$weight_sum
  }
  gr <- function(eta) {
    detail <- .inc_cell_detail(.unpack_inc(eta, stationary), cells$pattern,
                               stationary)
    -colSums(detail$scores * cells$weight) / cells$weight_sum
  }
  candidates <- lapply(seq_along(starts), function(i) {
    if (verbose) cat(sprintf("  start %d/%d [stationary=%s]\n",
                             i, length(starts), stationary))
    start <- starts[[i]]
    eta0 <- .pack_inc(start, stationary)
    opt <- tryCatch(optim(eta0, fn, gr, method = "BFGS",
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
  grad <- gr(eta)
  hessian <- optimHess(eta, fn, gr)
  hessian <- (hessian + t(hessian)) / 2
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  # Conditional probabilities provide seven independent outcome probabilities
  # for every distinct six-indicator history.  Stack their Jacobians to test
  # whether all local parameter directions change the observable distribution.
  x_names <- setdiff(names(cells$pattern), c("y1", "y2", "y3"))
  x_unique <- unique(cells$pattern[, x_names, drop = FALSE])
  outcomes <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
  probability_map <- function(z) {
    blocks <- lapply(seq_len(nrow(x_unique)), function(i) {
      pp <- cbind(as.data.frame(outcomes),
                  x_unique[rep(i, 8L), , drop = FALSE])
      .inc_cell_probabilities(.unpack_inc(z, stationary), pp, stationary)[1:7]
    })
    unlist(blocks, use.names = FALSE)
  }
  probability_jacobian <- .numeric_jacobian_inc(probability_map, eta)
  rank <- .matrix_rank_inc(probability_jacobian)

  # Individual-level survey-weighted sandwich covariance, calculated exactly
  # from the collapsed pattern cells using sums of w_i and w_i^2.
  score <- .inc_cell_detail(params, cells$pattern, stationary)$scores
  meat <- matrix(0, length(eta), length(eta))
  for (i in seq_len(nrow(score)))
    meat <- meat + cells$weight_sq[i] * tcrossprod(score[i, ])
  meat <- meat / cells$weight_sum^2
  bread_inv <- tryCatch(solve(hessian), error = function(e) MASS::ginv(hessian))
  vcov <- bread_inv %*% meat %*% bread_inv
  dimnames(vcov) <- list(names(eta), names(eta))

  if (verbose) cat(sprintf(
    "Table 3 direct MLE [stationary=%s]: ll=%.3f code=%d max|score|=%.2e minEig=%.2e rank=%d/%d\n",
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
  scenario <- function(...) {
    values <- c(...); z <- setNames(rep(0, length(.base_delta_names) - 1L),
      setdiff(.base_delta_names, "delta0"))
    z[names(values)] <- values
    unname(0.5 * plogis(p$delta["delta0"] + sum(p$delta[names(z)] * z)))
  }
  pi0 <- scenario()
  pi_age <- scenario(delta_age = 1)
  pi_edu <- scenario(delta_education = 1)
  pi_race <- scenario(delta_race = 1)
  pi_gender <- scenario(delta_gender = 1)
  pi_both <- scenario(delta_age = 1, delta_education = 1,
                      delta_two_inconsistencies = 1)
  pi_mat <- 0.5 * plogis(.base_error_index(p$delta, cells$pattern))
  mean_pi_weighted <- sum(cells$weight * rowSums(pi_mat)) / (3 * cells$weight_sum)
  mean_pi_unweighted <- sum(cells$count * rowSums(pi_mat)) / (3 * cells$n)
  steady <- p$theta0 / (p$theta0 + 1 - p$theta1)
  c(entry_rate = p$theta0, exit_rate = 1 - p$theta1,
    initial_employment = p$alpha, steady_employment = steady,
    stationarity_gap = p$alpha - steady,
    delta0 = unname(p$delta["delta0"]),
    delta_age = unname(p$delta["delta_age"]),
    delta_education = unname(p$delta["delta_education"]),
    delta_race = unname(p$delta["delta_race"]),
    delta_gender = unname(p$delta["delta_gender"]),
    delta_two_inconsistencies = unname(p$delta["delta_two_inconsistencies"]),
    delta_three_inconsistencies = unname(p$delta["delta_three_inconsistencies"]),
    delta_four_inconsistencies = unname(p$delta["delta_four_inconsistencies"]),
    pi_base = pi0,
    pi_age = pi_age, age_effect = pi_age - pi0,
    pi_education = pi_edu, education_effect = pi_edu - pi0,
    pi_race = pi_race, race_effect = pi_race - pi0,
    pi_gender = pi_gender, gender_effect = pi_gender - pi0,
    pi_both = pi_both, .conditional_pi_by_count(pi_mat, cells),
    .pi_distribution_summary(pi_mat, cells),
    mean_pi_unweighted = mean_pi_unweighted,
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
                            label, design_function, stationary) {
  design <- design_function(cells$pattern)
  fn <- function(eta) {
    p <- pmax(.inc_cell_detail(unpack(eta), cells$pattern, stationary,
                              design)$probability, 1e-300)
    -sum(cells$weight * log(p)) / cells$weight_sum
  }
  gr <- function(eta) {
    detail <- .inc_cell_detail(unpack(eta), cells$pattern, stationary, design)
    -colSums(detail$scores * cells$weight) / cells$weight_sum
  }
  candidates <- lapply(seq_along(starts), function(i) {
    cat(sprintf("  %s: start %d/%d\n", label, i, length(starts)))
    eta0 <- starts[[i]]
    eta0 <- setNames(as.numeric(eta0), eta_names)
    opt <- tryCatch(optim(eta0, fn, gr, method = "BFGS",
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
  grad <- gr(eta)
  hessian <- optimHess(eta, fn, gr); hessian <- (hessian + t(hessian)) / 2
  eig <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values

  design_names <- setdiff(names(cells$pattern), c("y1", "y2", "y3"))
  design_unique <- unique(cells$pattern[, design_names, drop = FALSE])
  outcomes <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
  probability_map <- function(z) unlist(lapply(seq_len(nrow(design_unique)), function(i) {
    pp <- cbind(as.data.frame(outcomes), design_unique[rep(i, 8L), , drop = FALSE])
    probabilities(unpack(z), pp)[1:7]
  }), use.names = FALSE)
  rank <- .matrix_rank_inc(.numeric_jacobian_inc(probability_map, eta))

  score <- .inc_cell_detail(unpack(eta), cells$pattern, stationary, design)$scores
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
                      .base_delta_names)

.unpack_group <- function(eta) {
  eta <- setNames(as.numeric(eta), .group_eta_names)
  list(theta0 = plogis(eta[c("theta0_reliable", "theta0_unreliable")]),
       theta1 = plogis(eta[c("theta1_reliable", "theta1_unreliable")]),
       alpha = plogis(eta[c("alpha_reliable", "alpha_unreliable")]),
       delta = eta[.base_delta_names])
}

.group_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  design_names <- unlist(lapply(c("age", "edu", "race", "gender"), function(x)
    paste0("Y_", x, "_", 1:3)))
  unreliable <- as.integer(rowSums(as.matrix(pattern[, design_names])) > 0) + 1L
  pi_mat <- 0.5 * plogis(.base_error_index(params$delta, pattern))
  out <- numeric(nrow(pattern))
  for (g in 1:2) {
    idx <- which(unreliable == g)
    if (!length(idx)) next
    prior <- prior_over_histories(latent, params$theta1[g], params$theta0[g],
                                  params$alpha[g])
    for (h in seq_len(nrow(latent))) {
      e1 <- ifelse(y[idx, 1L] == latent[h, 1L],
                   1 - pi_mat[idx, 1L], pi_mat[idx, 1L])
      e2 <- ifelse(y[idx, 2L] == latent[h, 2L],
                   1 - pi_mat[idx, 2L], pi_mat[idx, 2L])
      e3 <- ifelse(y[idx, 3L] == latent[h, 3L],
                   1 - pi_mat[idx, 3L], pi_mat[idx, 3L])
      out[idx] <- out[idx] + prior[h] * e1 * e2 * e3
    }
  }
  out
}

.group_quantities <- function(eta, cells) {
  p <- .unpack_group(eta)
  pi0 <- unname(0.5 * plogis(p$delta["delta0"]))
  pi_age <- unname(0.5 * plogis(p$delta["delta0"] + p$delta["delta_age"]))
  pi_edu <- unname(0.5 * plogis(p$delta["delta0"] + p$delta["delta_education"]))
  pi_race <- unname(0.5 * plogis(p$delta["delta0"] + p$delta["delta_race"]))
  pi_gender <- unname(0.5 * plogis(p$delta["delta0"] + p$delta["delta_gender"]))
  pi_mat <- 0.5 * plogis(.base_error_index(p$delta, cells$pattern))
  c(entry_reliable = unname(p$theta0[1L]),
    exit_reliable = unname(1 - p$theta1[1L]),
    initial_reliable = unname(p$alpha[1L]),
    entry_unreliable = unname(p$theta0[2L]),
    exit_unreliable = unname(1 - p$theta1[2L]),
    initial_unreliable = unname(p$alpha[2L]),
    delta0 = unname(p$delta["delta0"]),
    delta_age = unname(p$delta["delta_age"]),
    delta_education = unname(p$delta["delta_education"]),
    delta_race = unname(p$delta["delta_race"]),
    delta_gender = unname(p$delta["delta_gender"]),
    pi_base = pi0, pi_age = pi_age, age_effect = pi_age - pi0,
    pi_education = pi_edu, education_effect = pi_edu - pi0,
    pi_race = pi_race, race_effect = pi_race - pi0,
    pi_gender = pi_gender, gender_effect = pi_gender - pi0,
    mean_pi_survey_weighted = sum(cells$weight * rowSums(pi_mat)) /
      (3 * cells$weight_sum))
}

.extent_eta_names <- function(stationary) c("theta0", "theta1",
  if (!stationary) "alpha", .base_delta_names,
  "delta_age_severe", "delta_education_severe")

.extent_error_design <- function(pattern) {
  base <- .base_error_design(pattern)
  lapply(1:3, function(tt) cbind(base[[tt]],
    delta_age_severe = pattern[[paste0("severe_age_", tt)]],
    delta_education_severe = pattern[[paste0("severe_edu_", tt)]]))
}

.unpack_extent <- function(eta, stationary) {
  eta <- setNames(as.numeric(eta), .extent_eta_names(stationary))
  theta0 <- plogis(eta["theta0"]); theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha(theta0, theta1) else plogis(eta["alpha"])
  list(theta0 = theta0, theta1 = theta1, alpha = alpha,
       delta = eta[c(.base_delta_names,
                     "delta_age_severe", "delta_education_severe")])
}

.extent_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  prior <- prior_over_histories(latent, params$theta1, params$theta0, params$alpha)
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  age_extent <- as.matrix(pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(pattern[, paste0("severe_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(.base_error_index(params$delta, pattern) +
    params$delta["delta_age_severe"] * age_extent +
    params$delta["delta_education_severe"] * edu_extent)
  out <- numeric(nrow(pattern))
  for (h in seq_len(nrow(latent))) {
    e1 <- ifelse(y[, 1L] == latent[h, 1L], 1 - pi_mat[, 1L], pi_mat[, 1L])
    e2 <- ifelse(y[, 2L] == latent[h, 2L], 1 - pi_mat[, 2L], pi_mat[, 2L])
    e3 <- ifelse(y[, 3L] == latent[h, 3L], 1 - pi_mat[, 3L], pi_mat[, 3L])
    out <- out + prior[h] * e1 * e2 * e3
  }
  out
}

.extent_quantities <- function(eta, cells, stationary) {
  p <- .unpack_extent(eta, stationary); d <- p$delta
  scenario_pi <- function(age = 0, edu = 0, race = 0, gender = 0,
                          age_extent = 0, edu_extent = 0,
                          two = 0, three = 0, four = 0)
    unname(0.5 * plogis(d["delta0"] + d["delta_age"] * age +
      d["delta_education"] * edu + d["delta_race"] * race +
      d["delta_gender"] * gender +
      d["delta_two_inconsistencies"] * two +
      d["delta_three_inconsistencies"] * three +
      d["delta_four_inconsistencies"] * four +
      d["delta_age_severe"] * age_extent +
      d["delta_education_severe"] * edu_extent))
  age_extent <- as.matrix(cells$pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(cells$pattern[, paste0("severe_edu_", 1:3)])
  pi_mat <- 0.5 * plogis(.base_error_index(d, cells$pattern) +
    d["delta_age_severe"] * age_extent +
    d["delta_education_severe"] * edu_extent)
  steady <- stationary_alpha(p$theta0, p$theta1)
  c(entry_rate = unname(p$theta0), exit_rate = unname(1 - p$theta1),
    initial_employment = unname(p$alpha), steady_employment = unname(steady),
    stationarity_gap = unname(p$alpha - steady),
    delta0 = unname(d["delta0"]), delta_age = unname(d["delta_age"]),
    delta_education = unname(d["delta_education"]),
    delta_race = unname(d["delta_race"]),
    delta_gender = unname(d["delta_gender"]),
    delta_two_inconsistencies = unname(d["delta_two_inconsistencies"]),
    delta_three_inconsistencies = unname(d["delta_three_inconsistencies"]),
    delta_four_inconsistencies = unname(d["delta_four_inconsistencies"]),
    delta_age_severe = unname(d["delta_age_severe"]),
    delta_education_severe = unname(d["delta_education_severe"]),
    pi_base = scenario_pi(),
    pi_age_mild = scenario_pi(age = 1),
    pi_age_severe = scenario_pi(age = 1, age_extent = 1),
    age_severity_effect = scenario_pi(age = 1, age_extent = 1) - scenario_pi(age = 1),
    pi_education_mild = scenario_pi(edu = 1),
    pi_education_severe = scenario_pi(edu = 1, edu_extent = 1),
    education_severity_effect = scenario_pi(edu = 1, edu_extent = 1) - scenario_pi(edu = 1),
    pi_race = scenario_pi(race = 1),
    race_effect = scenario_pi(race = 1) - scenario_pi(),
    pi_gender = scenario_pi(gender = 1),
    gender_effect = scenario_pi(gender = 1) - scenario_pi(),
    pi_both_mild = scenario_pi(age = 1, edu = 1, two = 1),
    pi_both_age_severe = scenario_pi(age = 1, edu = 1, age_extent = 1, two = 1),
    pi_both_education_severe = scenario_pi(age = 1, edu = 1, edu_extent = 1,
                                           two = 1),
    pi_both_both_severe = scenario_pi(age = 1, edu = 1,
                                      age_extent = 1, edu_extent = 1, two = 1),
    .conditional_pi_by_count(pi_mat, cells),
    .pi_distribution_summary(pi_mat, cells),
    mean_pi_survey_weighted = sum(cells$weight * rowSums(pi_mat)) /
      (3 * cells$weight_sum))
}

.matching_eta_names <- function(stationary) c(
  .extent_eta_names(stationary), "delta_B_not_C", "delta_A_not_B")

.matching_error_design <- function(pattern) {
  extent <- .extent_error_design(pattern)
  lapply(1:3, function(tt) cbind(extent[[tt]],
    delta_B_not_C = pattern[[paste0("panel_B_not_C_", tt)]],
    delta_A_not_B = pattern[[paste0("panel_A_not_B_", tt)]]))
}

.unpack_matching <- function(eta, stationary) {
  eta <- setNames(as.numeric(eta), .matching_eta_names(stationary))
  theta0 <- plogis(eta["theta0"]); theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha(theta0, theta1) else plogis(eta["alpha"])
  list(theta0 = theta0, theta1 = theta1, alpha = alpha,
       delta = eta[c(.base_delta_names, "delta_age_severe",
                     "delta_education_severe", "delta_B_not_C",
                     "delta_A_not_B")])
}

.matching_probabilities <- function(params, pattern) {
  latent <- latent_histories()
  prior <- prior_over_histories(latent, params$theta1, params$theta0, params$alpha)
  y <- as.matrix(pattern[, c("y1", "y2", "y3")])
  age_extent <- as.matrix(pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(pattern[, paste0("severe_edu_", 1:3)])
  b_not_c <- as.matrix(pattern[, paste0("panel_B_not_C_", 1:3)])
  a_not_b <- as.matrix(pattern[, paste0("panel_A_not_B_", 1:3)])
  d <- params$delta
  pi_mat <- 0.5 * plogis(.base_error_index(d, pattern) +
    d["delta_age_severe"] * age_extent +
    d["delta_education_severe"] * edu_extent +
    d["delta_B_not_C"] * b_not_c + d["delta_A_not_B"] * a_not_b)
  out <- numeric(nrow(pattern))
  for (h in seq_len(nrow(latent))) {
    e1 <- ifelse(y[, 1L] == latent[h, 1L], 1 - pi_mat[, 1L], pi_mat[, 1L])
    e2 <- ifelse(y[, 2L] == latent[h, 2L], 1 - pi_mat[, 2L], pi_mat[, 2L])
    e3 <- ifelse(y[, 3L] == latent[h, 3L], 1 - pi_mat[, 3L], pi_mat[, 3L])
    out <- out + prior[h] * e1 * e2 * e3
  }
  out
}

.matching_quantities <- function(eta, cells, stationary) {
  p <- .unpack_matching(eta, stationary); d <- p$delta
  scenario_pi <- function(age = 0, edu = 0, race = 0, gender = 0,
                          age_extent = 0, edu_extent = 0,
                          b_not_c = 0, a_not_b = 0,
                          two = 0, three = 0, four = 0)
    unname(0.5 * plogis(d["delta0"] + d["delta_age"] * age +
      d["delta_education"] * edu + d["delta_race"] * race +
      d["delta_gender"] * gender +
      d["delta_two_inconsistencies"] * two +
      d["delta_three_inconsistencies"] * three +
      d["delta_four_inconsistencies"] * four +
      d["delta_age_severe"] * age_extent +
      d["delta_education_severe"] * edu_extent +
      d["delta_B_not_C"] * b_not_c + d["delta_A_not_B"] * a_not_b))
  age_extent <- as.matrix(cells$pattern[, paste0("severe_age_", 1:3)])
  edu_extent <- as.matrix(cells$pattern[, paste0("severe_edu_", 1:3)])
  b_not_c <- as.matrix(cells$pattern[, paste0("panel_B_not_C_", 1:3)])
  a_not_b <- as.matrix(cells$pattern[, paste0("panel_A_not_B_", 1:3)])
  pi_mat <- 0.5 * plogis(.base_error_index(d, cells$pattern) +
    d["delta_age_severe"] * age_extent +
    d["delta_education_severe"] * edu_extent +
    d["delta_B_not_C"] * b_not_c + d["delta_A_not_B"] * a_not_b)
  steady <- stationary_alpha(p$theta0, p$theta1)
  pi_base <- scenario_pi()
  pi_age_mild <- scenario_pi(age = 1)
  pi_age_severe <- scenario_pi(age = 1, age_extent = 1)
  pi_education_mild <- scenario_pi(edu = 1)
  pi_education_severe <- scenario_pi(edu = 1, edu_extent = 1)
  pi_race <- scenario_pi(race = 1)
  pi_gender <- scenario_pi(gender = 1)
  pi_both_mild <- scenario_pi(age = 1, edu = 1, two = 1)
  pi_both_both_severe <- scenario_pi(age = 1, edu = 1,
                                      age_extent = 1, edu_extent = 1, two = 1)
  pi_b_not_c <- scenario_pi(b_not_c = 1)
  pi_a_not_b <- scenario_pi(a_not_b = 1)
  c(entry_rate = unname(p$theta0), exit_rate = unname(1 - p$theta1),
    initial_employment = unname(p$alpha), steady_employment = unname(steady),
    stationarity_gap = unname(p$alpha - steady),
    delta0 = unname(d["delta0"]), delta_age = unname(d["delta_age"]),
    delta_education = unname(d["delta_education"]),
    delta_race = unname(d["delta_race"]),
    delta_gender = unname(d["delta_gender"]),
    delta_two_inconsistencies = unname(d["delta_two_inconsistencies"]),
    delta_three_inconsistencies = unname(d["delta_three_inconsistencies"]),
    delta_four_inconsistencies = unname(d["delta_four_inconsistencies"]),
    delta_age_severe = unname(d["delta_age_severe"]),
    delta_education_severe = unname(d["delta_education_severe"]),
    delta_B_not_C = unname(d["delta_B_not_C"]),
    delta_A_not_B = unname(d["delta_A_not_B"]),
    pi_base = pi_base,
    pi_age_mild = pi_age_mild,
    pi_age_severe = pi_age_severe,
    age_severity_effect = pi_age_severe - pi_age_mild,
    pi_education_mild = pi_education_mild,
    pi_education_severe = pi_education_severe,
    education_severity_effect = pi_education_severe - pi_education_mild,
    pi_race = pi_race, race_effect = pi_race - pi_base,
    pi_gender = pi_gender, gender_effect = pi_gender - pi_base,
    pi_both_mild = pi_both_mild,
    pi_both_both_severe = pi_both_both_severe,
    pi_B_not_C = pi_b_not_c,
    B_not_C_effect = pi_b_not_c - pi_base,
    pi_A_not_B = pi_a_not_b,
    A_not_B_effect = pi_a_not_b - pi_base,
    .conditional_pi_by_count(pi_mat, cells),
    .pi_distribution_summary(pi_mat, cells),
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

.table2_person_wave_keys <- function(panel) {
  raw <- readRDS(here::here("data", "raw", paste0("df_qlfs_", panel, ".rds")))
  num <- function(x) as.numeric(unclass(x))
  weight <- (num(raw$weight1) * num(raw$weight2) * num(raw$weight3))^(1 / 3)
  keep <- num(raw$age1) > 17 & num(raw$age1) < 56 &
    complete.cases(raw[paste0("employed", 1:3)]) &
    is.finite(weight) & weight > 0
  raw <- raw[keep, c("hhnr", "pnr", paste0("period", 1:3)), drop = FALSE]
  keys <- unique(unlist(lapply(1:3, function(tt)
    paste(raw$hhnr, raw$pnr, raw[[paste0("period", tt)]], sep = "|")),
    use.names = FALSE))
  rm(raw); gc(verbose = FALSE)
  keys
}

# Match the inconsistency-model analysis sample.
keep <- complete.cases(df_qlfs[, c("y1", "y2", "y3", "weight",
                                   "age1", "age2", "age3",
                                   "educ1", "educ2", "educ3")]) &
  df_qlfs$weight > 0
df <- as.data.frame(df_qlfs[keep, , drop = FALSE])
df$y1 <- as.integer(df$y1); df$y2 <- as.integer(df$y2); df$y3 <- as.integer(df$y3)
df$weight <- as.numeric(df$weight)
inc_df <- add_inconsistency_count_dummies(
  compute_demographic_inconsistencies(compute_inconsistency_extent(df)))
component_inc_names <- unlist(lapply(c("age", "edu", "race", "gender"), function(x)
  paste0("Y_", x, "_", 1:3)))
count_inc_names <- unlist(lapply(2:4, function(k) paste0("Y_exactly_", k, "_", 1:3)))
inc_names <- c(component_inc_names, count_inc_names)
inc_mat <- as.matrix(inc_df[, inc_names])
component_inc_mat <- as.matrix(inc_df[, component_inc_names])
cells <- .collapse_inc_cells(df, inc_mat)
extent_names <- c(paste0("extent_age_", 1:3), paste0("extent_edu_", 1:3))
severe_names <- c(paste0("severe_age_", 1:3), paste0("severe_edu_", 1:3))
inc_df[severe_names] <- lapply(inc_df[extent_names], function(x) as.integer(x >= 2))
cells_extent <- .collapse_inc_cells(
  df, as.matrix(inc_df[, c(inc_names, severe_names)])
)

# Construct person-wave set differences using exactly the complete-panel and
# weight filters used for Table 2. These are wave-varying measurement-quality
# indicators, not person-level sample restrictions.
panel_B_keys <- .table2_person_wave_keys("B")
panel_C_keys <- .table2_person_wave_keys("C")
panel_A_key_mat <- sapply(1:3, function(tt)
  paste(df$hhnr, df$pnr, df[[paste0("period", tt)]], sep = "|"))
in_B <- matrix(panel_A_key_mat %in% panel_B_keys, nrow(df), 3L)
in_C <- matrix(panel_A_key_mat %in% panel_C_keys, nrow(df), 3L)
panel_B_not_C <- 1L * (in_B & !in_C)
panel_A_not_B <- 1L * !in_B
colnames(panel_B_not_C) <- paste0("panel_B_not_C_", 1:3)
colnames(panel_A_not_B) <- paste0("panel_A_not_B_", 1:3)
stopifnot(sum(panel_B_not_C * panel_A_not_B) == 0L,
          sum(panel_B_not_C) > 0L, sum(panel_A_not_B) > 0L)
inc_df <- cbind(inc_df, as.data.frame(panel_B_not_C),
                as.data.frame(panel_A_not_B))
matching_names <- c(inc_names, severe_names,
                    colnames(panel_B_not_C), colnames(panel_A_not_B))
cells_matching <- .collapse_inc_cells(df, as.matrix(inc_df[, matching_names]))
membership_summary <- data.frame(
  indicator = c("Panel B but not C", "Panel A but not B"),
  unweighted_percent = 100 * c(mean(panel_B_not_C), mean(panel_A_not_B)),
  survey_weighted_percent = 100 * c(
    sum(panel_B_not_C * df$weight) / (3 * sum(df$weight)),
    sum(panel_A_not_B * df$weight) / (3 * sum(df$weight)))
)
rm(panel_B_keys, panel_C_keys, panel_A_key_mat, in_B, in_C)
gc(verbose = FALSE)

cat(sprintf("Table 3 sample: N=%s; %d collapsed likelihood cells\n",
            format(nrow(df), big.mark = ","), nrow(cells$pattern)))
cat(sprintf("Severity robustness: %d collapsed likelihood cells\n",
            nrow(cells_extent$pattern)))
cat(sprintf("Matching-quality extension: %d collapsed likelihood cells\n",
            nrow(cells_matching$pattern)))
flag_summary <- data.frame(
  indicator = c(inc_names, "any_age", "any_education", "any_race",
                "any_gender", "any_inconsistency"),
  unweighted_percent = 100 * c(colMeans(inc_mat),
    mean(rowSums(component_inc_mat[, 1:3, drop = FALSE]) > 0),
    mean(rowSums(component_inc_mat[, 4:6, drop = FALSE]) > 0),
    mean(rowSums(component_inc_mat[, 7:9, drop = FALSE]) > 0),
    mean(rowSums(component_inc_mat[, 10:12, drop = FALSE]) > 0),
    mean(rowSums(component_inc_mat) > 0)),
  survey_weighted_percent = 100 * c(
    colSums(inc_mat * df$weight) / sum(df$weight),
    weighted.mean(rowSums(component_inc_mat[, 1:3, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(component_inc_mat[, 4:6, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(component_inc_mat[, 7:9, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(component_inc_mat[, 10:12, drop = FALSE]) > 0, df$weight),
    weighted.mean(rowSums(component_inc_mat) > 0, df$weight))
)
cat("\nInconsistency prevalence (%):\n")
print(flag_summary, row.names = FALSE, digits = 4)
severity_summary <- data.frame(
  variable = c("age", "education"),
  severe_share_of_flagged_percent = 100 * c(
    sum(df$weight * rowSums(as.matrix(inc_df[, paste0("severe_age_", 1:3)]))) /
      sum(df$weight * rowSums(component_inc_mat[, 1:3, drop = FALSE])),
    sum(df$weight * rowSums(as.matrix(inc_df[, paste0("severe_edu_", 1:3)]))) /
      sum(df$weight * rowSums(component_inc_mat[, 4:6, drop = FALSE]))
  )
)
cat("\nShare of attributed flags at least two units beyond the admissible range (%):\n")
print(severity_summary, row.names = FALSE, digits = 4)
cat("\nTable 2 matching-set indicators (% of person-wave observations):\n")
print(membership_summary, row.names = FALSE, digits = 4)

# Obtain nested homogeneous estimates for stable starts, then add dispersed
# slope starts to guard against local optima.
set.seed(62026L)
save_retry <- identical(Sys.getenv("TABLE3_SAVE_RETRY"), "1")
prior_audit_path <- here::here("EM-baseline-ext", "output", "results",
                              "fit_table6_inconsistency_audit.rds")
if (!file.exists(prior_audit_path))
  stop("Existing Table 3 audit required for nested warm starts: ",
       prior_audit_path)
prior_audit <- readRDS(prior_audit_path)
reuse_audited_fits <- identical(Sys.getenv("TABLE3_REUSE_AUDITED_FITS"), "1")
make_starts <- function(stationary) {
  old <- if (stationary) prior_audit$stationary else prior_audit$free
  delta <- setNames(rep(0, length(.base_delta_names)), .base_delta_names)
  delta[names(old$params$delta)] <- old$params$delta
  center <- list(theta0 = old$params$theta0, theta1 = old$params$theta1,
                 alpha = old$params$alpha, delta = delta)
  fixed <- list(center,
    within(center, delta[c("delta_two_inconsistencies",
                           "delta_three_inconsistencies")] <- c(0.5, -0.5)),
    within(center, delta[c("delta_two_inconsistencies",
                           "delta_four_inconsistencies")] <- c(-0.5, 1.0)),
    within(center, delta[c("delta_three_inconsistencies",
                           "delta_four_inconsistencies")] <- c(0.5, -1.0)))
  if (save_retry) fixed[1L] else fixed
}

if (reuse_audited_fits) {
  fit_stat <- prior_audit$stationary
  fit_free <- prior_audit$free
} else {
  fit_stat <- .fit_inc_mle(cells, TRUE, make_starts(TRUE))
  fit_free <- .fit_inc_mle(cells, FALSE, make_starts(FALSE))
}
inf_stat <- .inc_inference(fit_stat, cells)
inf_free <- .inc_inference(fit_free, cells)

# The reliability-group transition robustness is not one of the three Table 3
# specifications. Retain its last audited result instead of needlessly
# re-estimating it whenever the Table 3 error equation is extended.
fit_group <- prior_audit$reliability_group
inf_group <- prior_audit$reliability_group_table

# Robustness 2: add the distance from the admissible [0,1] age/education
# change to test whether more severe inconsistencies imply larger error rates.
fit_extent_model <- function(stationary, base_fit, label) {
  center <- c(base_fit$eta, delta_age_severe = 0,
              delta_education_severe = 0)
  names_eta <- .extent_eta_names(stationary)
  center <- center[names_eta]
  if (save_retry) {
    prob <- if (stationary)
      c(theta0 = 0.0227593, theta1 = 1 - 0.0247272) else
      c(theta0 = 0.025703608, theta1 = 1 - 0.021409684,
        alpha = 0.476008308)
    center <- c(qlogis(prob),
      delta0 = if (stationary) -3.06737 else -3.066903156,
      delta_age = if (stationary) 1.13303 else 1.132267586,
      delta_education = if (stationary) 1.01634 else 1.016228431,
      delta_race = if (stationary) 0.00223739 else 0.013256203,
      delta_gender = if (stationary) 1.50392 else 1.506068283,
      delta_two_inconsistencies = 0,
      delta_three_inconsistencies = 0,
      delta_four_inconsistencies = 0,
      delta_age_severe = if (stationary) 1.23743 else 1.240929301,
      delta_education_severe = if (stationary) 0.114517 else 0.116561070)
    center <- center[names_eta]
  }
  perturb_sd <- ifelse(grepl("^theta|^alpha", names_eta), .30,
    ifelse(names_eta == "delta0" | grepl("severe", names_eta), .25, .75))
  starts <- c(list(center), lapply(1:3, function(i)
    center + rnorm(length(center), 0, perturb_sd)))
  if (save_retry) starts <- starts[1L]
  .fit_custom_inc(cells_extent, names_eta,
    function(z) .unpack_extent(z, stationary), .extent_probabilities,
    starts, label, .extent_error_design, stationary)
}
if (reuse_audited_fits) {
  fit_extent_stat <- prior_audit$extent_stationary
  fit_extent_free <- prior_audit$extent_free
} else {
  fit_extent_stat <- fit_extent_model(TRUE, fit_stat,
    "Stationary inconsistency-increment model")
  fit_extent_free <- fit_extent_model(FALSE, fit_free,
    "Free-alpha inconsistency-increment model")
}
inf_extent_stat <- .custom_inference(fit_extent_stat,
  function(z, cc) .extent_quantities(z, cc, TRUE), cells_extent)
inf_extent_free <- .custom_inference(fit_extent_free,
  function(z, cc) .extent_quantities(z, cc, FALSE), cells_extent)
extent_table <- merge(inf_extent_stat, inf_extent_free, by = "quantity",
  all = TRUE, suffixes = c("_stationary", "_free"), sort = FALSE)
extent_order <- names(.extent_quantities(fit_extent_stat$eta, cells_extent, TRUE))
extent_table <- extent_table[match(extent_order, extent_table$quantity), ]

# Main extension requested for Table 3: retain the age/education indicators and
# severity increments, and add both Table 2 person-wave set-difference flags to
# the misclassification equation.
fit_matching_model <- function(stationary, base_fit, label) {
  center <- c(base_fit$eta, delta_B_not_C = 0, delta_A_not_B = 0)
  names_eta <- .matching_eta_names(stationary)
  center <- center[names_eta]
  if (save_retry) {
    prob <- if (stationary)
      c(theta0 = 0.02239912, theta1 = 1 - 0.02433976) else
      c(theta0 = 0.02535111, theta1 = 1 - 0.02100869,
        alpha = 0.47595867)
    center <- c(qlogis(prob),
      delta0 = if (stationary) -3.17466244 else -3.17425895,
      delta_age = if (stationary) 0.93009547 else 0.92920944,
      delta_education = if (stationary) 0.82233147 else 0.82188567,
      delta_race = if (stationary) -0.17621279 else -0.16970159,
      delta_gender = if (stationary) 1.30379635 else 1.30369674,
      delta_two_inconsistencies = 0,
      delta_three_inconsistencies = 0,
      delta_four_inconsistencies = 0,
      delta_age_severe = if (stationary) 1.28652910 else 1.28990145,
      delta_education_severe = if (stationary) 0.09994649 else 0.10206079,
      delta_B_not_C = if (stationary) 0.39652703 else 0.39663685,
      delta_A_not_B = if (stationary) 0.39649365 else 0.39843502)
    center <- center[names_eta]
  }
  perturb_sd <- ifelse(grepl("^theta|^alpha", names_eta), .30,
    ifelse(names_eta == "delta0" | grepl("severe", names_eta), .25, .75))
  starts <- c(list(center), lapply(1:3, function(i)
    center + rnorm(length(center), 0, perturb_sd)))
  if (save_retry) starts <- starts[1L]
  .fit_custom_inc(cells_matching, names_eta,
    function(z) .unpack_matching(z, stationary), .matching_probabilities,
    starts, label, .matching_error_design, stationary)
}
if (reuse_audited_fits) {
  fit_matching_stat <- prior_audit$matching_stationary
  fit_matching_free <- prior_audit$matching_free
} else {
  fit_matching_stat <- fit_matching_model(TRUE, fit_extent_stat,
    "Stationary matching-quality inconsistency model")
  fit_matching_free <- fit_matching_model(FALSE, fit_extent_free,
    "Free-alpha matching-quality inconsistency model")
}
if (!fit_matching_stat$converged || !fit_matching_stat$identified ||
    !fit_matching_free$converged || !fit_matching_free$identified)
  stop("Matching-quality Table 3 specifications failed convergence or identification checks")
if (fit_matching_stat$loglik < fit_extent_stat$loglik - 1e-4 ||
    fit_matching_free$loglik < fit_extent_free$loglik - 1e-4)
  stop("Matching-quality Table 3 likelihood does not dominate its nested model")
inf_matching_stat <- .custom_inference(fit_matching_stat,
  function(z, cc) .matching_quantities(z, cc, TRUE), cells_matching)
inf_matching_free <- .custom_inference(fit_matching_free,
  function(z, cc) .matching_quantities(z, cc, FALSE), cells_matching)
matching_table <- merge(inf_matching_stat, inf_matching_free, by = "quantity",
  all = TRUE, suffixes = c("_stationary", "_free"), sort = FALSE)
matching_order <- names(.matching_quantities(
  fit_matching_stat$eta, cells_matching, TRUE))
matching_table <- matching_table[match(matching_order, matching_table$quantity), ]

table6 <- merge(inf_stat, inf_free, by = "quantity", all = TRUE,
                suffixes = c("_stationary", "_free"), sort = FALSE)
order_q <- .inc_quantities(fit_stat$eta, TRUE, cells)
table6 <- table6[match(names(order_q), table6$quantity), ]
table6$estimate_stationary[table6$quantity == "initial_employment"] <- NA_real_
table6$std_error_stationary[table6$quantity == "initial_employment"] <- NA_real_

cat("\nReplicated Table 3 (probabilities are proportions; analytical SE):\n")
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
print(extent_table, row.names = FALSE, digits = 6)
cat(sprintf("Identified/converged: stationary %s/%s; free %s/%s\n",
            fit_extent_stat$identified, fit_extent_stat$converged,
            fit_extent_free$identified, fit_extent_free$converged))
cat("\nTable 3 extension: matching-quality indicators in the error equation:\n")
print(matching_table, row.names = FALSE, digits = 6)
cat(sprintf("Identified/converged: stationary %s/%s; free %s/%s\n",
            fit_matching_stat$identified, fit_matching_stat$converged,
            fit_matching_free$identified, fit_matching_free$converged))

out_dir <- here::here("EM-baseline-ext", "output", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(table6, file.path(out_dir, "table6_inconsistency_audit.csv"), row.names = FALSE)
write.csv(flag_summary, file.path(out_dir, "table6_inconsistency_prevalence.csv"), row.names = FALSE)
write.csv(severity_summary, file.path(out_dir, "table6_inconsistency_severity_prevalence.csv"), row.names = FALSE)
write.csv(inf_group, file.path(out_dir, "table6_reliability_group_robustness.csv"), row.names = FALSE)
write.csv(inf_extent_free, file.path(out_dir, "table6_inconsistency_extent.csv"), row.names = FALSE)
write.csv(extent_table, file.path(out_dir, "table6_inconsistency_extent_audit.csv"),
          row.names = FALSE)
write.csv(membership_summary,
          file.path(out_dir, "table6_matching_membership_prevalence.csv"),
          row.names = FALSE)
write.csv(matching_table,
          file.path(out_dir, "table6_matching_quality_audit.csv"),
          row.names = FALSE)
saveRDS(list(stationary = fit_stat, free = fit_free, table = table6,
             reliability_group = fit_group, reliability_group_table = inf_group,
             extent = fit_extent_free, extent_table = inf_extent_free,
             extent_stationary = fit_extent_stat, extent_free = fit_extent_free,
             extent_comparison_table = extent_table,
             matching_stationary = fit_matching_stat,
             matching_free = fit_matching_free,
             matching_comparison_table = matching_table,
             matching_membership_prevalence = membership_summary,
             prevalence = flag_summary, severity_prevalence = severity_summary),
        file.path(out_dir, "fit_table6_inconsistency_audit.rds"))
