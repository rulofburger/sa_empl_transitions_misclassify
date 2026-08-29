# Exact four-wave, one-type AR(2) likelihood with the full Table 4 Set 2
# transition covariates and reliability-dependent symmetric misclassification.

.ar2r_clamp <- function(x, lo = 1e-10, hi = 1 - 1e-10)
  pmin(pmax(x, lo), hi)

compute_inconsistency_extent_4w <- function(df) {
  required <- c(paste0("age", 1:4), paste0("educ", 1:4),
                paste0("race", 1:4), paste0("female", 1:4))
  if (length(missing <- setdiff(required, names(df))))
    stop("Missing four-wave reliability variables: ", paste(missing, collapse = ", "))
  edge <- function(prefix) {
    x <- as.matrix(df[paste0(prefix, 1:4)])
    d <- x[, 2:4, drop = FALSE] - x[, 1:3, drop = FALSE]
    bad <- !is.na(d) & !(abs(d) < .01 | abs(d - 1) < .01)
    distance <- ifelse(is.na(d), 0, ifelse(d < 0, -d, ifelse(d > 1, d - 1, 0)))
    list(bad = bad, distance = distance)
  }
  attribute_wave <- function(bad, distance) {
    # This is the four-wave rule already used by the FMM module: adjacent bad
    # edges are attributed to their shared report; an isolated first/last edge
    # is attributed to the outside endpoint.
    v1 <- bad[,1]; v2 <- bad[,2]; v3 <- bad[,3]
    d1 <- distance[,1]; d2 <- distance[,2]; d3 <- distance[,3]
    flag <- cbind(v1 & !v2, v1 & v2,
                  (!v1 & v2) | (v2 & v3), !v2 & v3) * 1L
    extent <- cbind(d1 * v1 * !v2,
                    (d1 + d2) * v1 * v2,
                    d2 * !v1 * v2 + (d2 + d3) * v2 * v3,
                    d3 * !v2 * v3)
    list(flag = flag, extent = extent)
  }
  age <- do.call(attribute_wave, edge("age"))
  educ <- do.call(attribute_wave, edge("educ"))
  categorical_edge <- function(prefix) {
    x <- as.matrix(df[paste0(prefix, 1:4)])
    !is.na(x[, 2:4, drop = FALSE]) & !is.na(x[, 1:3, drop = FALSE]) &
      x[, 2:4, drop = FALSE] != x[, 1:3, drop = FALSE]
  }
  categorical_wave <- function(bad)
    cbind(bad[,1] & !bad[,2], bad[,1] & bad[,2],
          (!bad[,1] & bad[,2]) | (bad[,2] & bad[,3]),
          !bad[,2] & bad[,3]) * 1L
  race <- categorical_wave(categorical_edge("race"))
  gender <- categorical_wave(categorical_edge("female"))
  out <- list()
  for (tt in 1:4) {
    out[[paste0("Y_age_", tt)]] <- age$flag[, tt]
    out[[paste0("Y_edu_", tt)]] <- educ$flag[, tt]
    out[[paste0("Y_race_", tt)]] <- race[, tt]
    out[[paste0("Y_gender_", tt)]] <- gender[, tt]
    out[[paste0("extent_age_", tt)]] <- age$extent[, tt]
    out[[paste0("extent_edu_", tt)]] <- educ$extent[, tt]
    count <- age$flag[, tt] + educ$flag[, tt] + race[, tt] + gender[, tt]
    for (kk in 2:4)
      out[[paste0("Y_exactly_", kk, "_", tt)]] <- as.integer(count == kk)
  }
  as.data.frame(out)
}

prepare_ar2_set4_reliability_4w <- function(
    panel_path = "data/raw/df_qlfs_A.rds",
    sector_path = "data/raw/QLFSmerged_mapped.rds", collapse = TRUE) {
  if (!exists("prepare_fmm_covariates_inconsistency_4w", mode = "function"))
    stop("Source EM-AR1-4W/R/source_all.R before preparing the Set 4 AR(2) data")
  base <- prepare_fmm_covariates_inconsistency_4w(
    panel_path, sector_path, collapse = FALSE)
  raw <- readRDS(panel_path)
  required <- c("hhnr", "pnr", paste0("period", 1:4), paste0("employed", 1:4),
                paste0("age", 1:4), paste0("educ", 1:4),
                paste0("race", 1:4), paste0("female", 1:4),
                paste0("weight", 1:4))
  filter_required <- c(paste0("employed", 1:4), paste0("age", 1:4),
                       paste0("educ", 1:4), "race1", "female1",
                       paste0("weight", 1:4))
  keep <- raw$age1 > 17 & raw$age1 < 56 &
    complete.cases(raw[filter_required])
  raw <- raw[keep, required, drop = FALSE]
  for (nm in c(paste0("age", 1:4), paste0("educ", 1:4),
               paste0("race", 1:4), paste0("female", 1:4)))
    raw[[nm]] <- as.numeric(unclass(raw[[nm]]))
  if (nrow(raw) != base$n_original) stop("Set 2 and reliability samples do not align")
  inc <- compute_inconsistency_extent_4w(raw)
  Z <- lapply(1:4, function(tt) cbind(
    error_intercept = 1,
    age_inconsistency = inc[[paste0("Y_age_", tt)]],
    education_inconsistency = inc[[paste0("Y_edu_", tt)]],
    race_inconsistency = inc[[paste0("Y_race_", tt)]],
    gender_inconsistency = inc[[paste0("Y_gender_", tt)]],
    two_inconsistencies = inc[[paste0("Y_exactly_2_", tt)]],
    three_inconsistencies = inc[[paste0("Y_exactly_3_", tt)]],
    four_inconsistencies = inc[[paste0("Y_exactly_4_", tt)]]))
  out <- list(y = base$y, weight = base$weight, weight_sq = base$weight_sq,
    X = base$X[2:3], Z = Z, entry_active = base$entry_active,
    persistence_active = base$persistence_active,
    covariate_names = base$covariate_names, scaling = base$scaling,
    n_original = base$n_original,
    reliability_prevalence = setNames(vapply(2:ncol(Z[[1]]), function(j)
      mean(vapply(Z, function(z) mean(z[, j]), numeric(1))), numeric(1)),
      colnames(Z[[1]])[-1]))
  if (!collapse) return(out)
  key_parts <- cbind(out$y, do.call(cbind, out$X),
                     do.call(cbind, lapply(Z, function(z) z[, -1, drop = FALSE])))
  key <- do.call(paste, c(as.data.frame(key_parts), sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  collapse_sum <- function(x) as.vector(rowsum(x, group, reorder = FALSE))
  out$y <- out$y[first, , drop = FALSE]
  out$X <- lapply(out$X, function(x) x[first, , drop = FALSE])
  out$Z <- lapply(out$Z, function(x) x[first, , drop = FALSE])
  out$weight <- collapse_sum(out$weight); out$weight_sq <- collapse_sum(out$weight_sq)
  out$n_cells <- length(out$weight)
  out
}

.ar2_piecewise_duration_bins <- function(months, prefix) {
  months <- as.numeric(months)
  months[!is.finite(months) | months < 0] <- 0
  lower <- c(3, 6, 12, 24, 60)
  upper <- c(6, 12, 24, 60, Inf)
  labels <- c("4_6m", "7_12m", "13_24m", "25_60m", "over_60m")
  out <- vapply(seq_along(lower), function(j)
    as.numeric(months > lower[j] & months <= upper[j]), numeric(length(months)))
  if (is.null(dim(out))) out <- matrix(out, ncol = length(lower))
  colnames(out) <- paste0(prefix, "_", labels)
  out
}

# Prepare the Table 5 column (4) model with piecewise-constant rather than
# log-linear duration hazards. The estimation sample and all non-duration
# covariates are identical to the Set 2 specification.
prepare_ar2_set2_piecewise_4w <- function(
    panel_path = "data/raw/df_qlfs_A.rds",
    sector_path = "data/raw/QLFSmerged_mapped.rds", collapse = TRUE) {
  if (!exists("prepare_fmm_covariates_inconsistency_4w", mode = "function"))
    stop("Source EM-AR1-4W/R/source_all.R before preparing piecewise AR(2) data")
  base <- prepare_fmm_covariates_inconsistency_4w(
    panel_path, sector_path, collapse = FALSE)
  panel <- readRDS(panel_path)
  filter_required <- c(paste0("employed", 1:4), paste0("age", 1:4),
    paste0("educ", 1:4), "race1", "female1", paste0("weight", 1:4))
  duration_required <- c(paste0("tenure", 1:3), paste0("timegap", 1:3),
                         paste0("neverworked", 1:3))
  missing <- setdiff(c(filter_required, duration_required), names(panel))
  if (length(missing))
    stop("Four-wave panel is missing: ", paste(missing, collapse = ", "))
  keep <- panel$age1 > 17 & panel$age1 < 56 &
    complete.cases(panel[filter_required])
  panel <- panel[keep, c(filter_required, duration_required), drop = FALSE]
  for (nm in names(panel)) panel[[nm]] <- as.numeric(unclass(panel[[nm]]))
  if (nrow(panel) != base$n_original)
    stop("Piecewise-duration and Set 2 samples do not align")

  timegap_months <- c(`1` = 1.5, `2` = 4.5, `3` = 7.5, `4` = 10.5,
                      `5` = 24, `6` = 48, `7` = 90)
  tenure_bins <- timegap_bins <- vector("list", 3L)
  for (tt in 1:3) {
    reported <- panel[[paste0("employed", tt)]]
    tenure <- panel[[paste0("tenure", tt)]]
    tenure[tenure < 0] <- NA_real_
    tenure <- ifelse(reported == 1 & !is.na(tenure), tenure, 0)

    gap_code <- panel[[paste0("timegap", tt)]]
    gap <- unname(timegap_months[as.character(gap_code)])
    gap[gap_code == 0] <- 0
    never <- panel[[paste0("neverworked", tt)]]
    # Time since work is undefined for never-worked respondents; their
    # separate indicator remains in the entry equation.
    gap[reported == 1 | (reported == 0 & never == 1) | is.na(gap)] <- 0
    tenure_bins[[tt]] <- .ar2_piecewise_duration_bins(tenure, "tenure")
    timegap_bins[[tt]] <- .ar2_piecewise_duration_bins(gap, "timegap")
  }

  X <- lapply(1:3, function(tt) {
    retained <- base$X[[tt]][, !colnames(base$X[[tt]]) %in%
      c("log_tenure", "log_time_since_work"), drop = FALSE]
    cbind(retained, tenure_bins[[tt]], timegap_bins[[tt]])
  })
  xnames <- colnames(X[[1]])
  tenure_names <- colnames(tenure_bins[[1]])
  timegap_names <- colnames(timegap_bins[[1]])
  entry_active <- !xnames %in% c(tenure_names, "tenure_missing",
    "permanent_contract", "informal_sector")
  persistence_active <- !xnames %in% c(timegap_names, "never_worked",
    "timegap_missing")
  Z <- lapply(1:4, function(tt)
    matrix(1, nrow(panel), 1L, dimnames = list(NULL, "error_intercept")))
  out <- list(y = base$y, weight = base$weight, weight_sq = base$weight_sq,
    X = X[2:3], Z = Z, entry_active = entry_active,
    persistence_active = persistence_active, covariate_names = xnames,
    scaling = c(base$scaling, list(duration_reference = "0--3 months",
      duration_bins = c("4--6", "7--12", "13--24", "25--60", ">60"))),
    n_original = base$n_original, stage = "piecewise_duration")
  if (!collapse) return(out)
  key_parts <- cbind(out$y, do.call(cbind, out$X))
  key <- do.call(paste, c(as.data.frame(key_parts), sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  collapse_sum <- function(x) as.vector(rowsum(x, group, reorder = FALSE))
  out$y <- out$y[first, , drop = FALSE]
  out$X <- lapply(out$X, function(x) x[first, , drop = FALSE])
  out$Z <- lapply(out$Z, function(x) x[first, , drop = FALSE])
  out$weight <- collapse_sum(out$weight)
  out$weight_sq <- collapse_sum(out$weight_sq)
  out$n_cells <- length(out$weight)
  out
}

subset_ar2_reliability_stage <- function(data,
    stage = c("constant", "table3_column1", "inconsistency", "reliability")) {
  stage <- match.arg(stage)
  keep <- switch(stage, constant = 1L,
                 table3_column1 = seq_len(ncol(data$Z[[1]])),
                 inconsistency = seq_len(min(5L, ncol(data$Z[[1]]))),
                 reliability = seq_len(ncol(data$Z[[1]])))
  data$Z <- lapply(data$Z, function(z) z[, keep, drop = FALSE])
  # The full preparation is collapsed on all reliability variables. Re-collapse
  # after selecting a stage so omitted reliability fields do not create
  # redundant likelihood cells in the earlier nested models.
  key_parts <- cbind(data$y, do.call(cbind, data$X),
    do.call(cbind, lapply(data$Z, function(z) z[, -1, drop = FALSE])))
  key <- do.call(paste, c(as.data.frame(key_parts), sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  if (sum(first) < length(first)) {
    collapse_sum <- function(x) as.vector(rowsum(x, group, reorder = FALSE))
    data$y <- data$y[first, , drop = FALSE]
    data$X <- lapply(data$X, function(x) x[first, , drop = FALSE])
    data$Z <- lapply(data$Z, function(x) x[first, , drop = FALSE])
    data$weight <- collapse_sum(data$weight)
    data$weight_sq <- collapse_sum(data$weight_sq)
  }
  data$n_cells <- length(data$weight)
  data$stage <- stage
  data
}

.ar2r_layout <- function(data) {
  active0 <- data$entry_active
  active1 <- data$persistence_active
  list(p0 = 2L + sum(active0), p1 = 2L + sum(active1), pd = ncol(data$Z[[1]]),
       active0 = active0, active1 = active1)
}

ar2_reliability_names <- function(data) {
  L <- .ar2r_layout(data)
  c(paste0("alpha_logit_", c("10", "01", "11")),
    paste0("entry_", c("intercept", "lag2", data$covariate_names[L$active0])),
    paste0("persistence_", c("intercept", "lag2", data$covariate_names[L$active1])),
    colnames(data$Z[[1]]))
}

.ar2r_unpack <- function(eta, data) {
  L <- .ar2r_layout(data); eta <- setNames(as.numeric(eta), ar2_reliability_names(data))
  shift <- max(c(0, eta[1:3])); a <- c(exp(-shift), exp(eta[1:3] - shift))
  j0 <- 3L + seq_len(L$p0); j1 <- 3L + L$p0 + seq_len(L$p1)
  jd <- 3L + L$p0 + L$p1 + seq_len(L$pd)
  list(alpha = setNames(a / sum(a), c("00", "10", "01", "11")),
    beta0 = eta[j0], beta1 = eta[j1], delta = eta[jd])
}

.ar2r_components <- function(eta, data, retain_scores = FALSE,
                             retain_posterior = FALSE) {
  p <- .ar2r_unpack(eta, data); h <- latent_histories_ar2()
  n <- nrow(data$y); H <- nrow(h); logc <- matrix(0, n, H)
  pair <- paste0(h[, 1], h[, 2]); logc <- sweep(logc, 2, log(p$alpha[pair]), "+")
  q0 <- q1 <- vector("list", 2)
  for (tt in 1:2) {
    x0 <- cbind(intercept = 1, lag2 = 0,
      data$X[[tt]][, .ar2r_layout(data)$active0, drop = FALSE])
    x1 <- cbind(intercept = 1, lag2 = 0,
      data$X[[tt]][, .ar2r_layout(data)$active1, drop = FALSE])
    q0[[tt]] <- q1[[tt]] <- vector("list", 2)
    for (lag in 0:1) {
      x0[, "lag2"] <- lag; x1[, "lag2"] <- lag
      q0[[tt]][[lag + 1L]] <- .ar2r_clamp(pnorm(as.vector(x0 %*% p$beta0)))
      q1[[tt]][[lag + 1L]] <- .ar2r_clamp(pnorm(as.vector(x1 %*% p$beta1)))
    }
    for (j in seq_len(H)) {
      q <- if (h[j, tt + 1L] == 0) q0[[tt]][[h[j, tt] + 1L]] else
        q1[[tt]][[h[j, tt] + 1L]]
      logc[, j] <- logc[, j] + if (h[j, tt + 2L]) log(q) else log1p(-q)
    }
  }
  pi_wave <- lapply(data$Z, function(z)
    .ar2r_clamp(.5 * plogis(as.vector(z %*% p$delta)), hi = .5 - 1e-10))
  for (j in seq_len(H)) for (tt in 1:4) {
    mismatch <- data$y[, tt] != h[j, tt]; pi <- pi_wave[[tt]]
    logc[, j] <- logc[, j] + ifelse(mismatch, log(pi), log1p(-pi))
  }
  mx <- apply(logc, 1, max); denom <- rowSums(exp(logc - mx))
  loglik_i <- mx + log(denom)
  if (!retain_scores && !retain_posterior)
    return(list(loglik_i = loglik_i, params = p, q0 = q0, q1 = q1,
                pi_wave = pi_wave))
  post <- exp(logc - mx) / denom
  if (!retain_scores) return(list(loglik_i = loglik_i, params = p,
    posterior = post, q0 = q0, q1 = q1, pi_wave = pi_wave))
  K <- length(eta); scores <- matrix(0, n, K)
  alpha_ind <- cbind(h[,1] == 1 & h[,2] == 0, h[,1] == 0 & h[,2] == 1,
                     h[,1] == 1 & h[,2] == 1) * 1
  scores[, 1:3] <- post %*% alpha_ind -
    matrix(p$alpha[c("10", "01", "11")], n, 3, byrow = TRUE)
  L <- .ar2r_layout(data); off0 <- 3L; off1 <- 3L + L$p0
  for (tt in 1:2) for (origin in 0:1) for (lag in 0:1) {
    active <- if (!origin) L$active0 else L$active1
    X <- cbind(intercept = 1, lag2 = lag, data$X[[tt]][, active, drop = FALSE])
    beta <- if (!origin) p$beta0 else p$beta1
    eta_t <- as.vector(X %*% beta); q <- .ar2r_clamp(pnorm(eta_t))
    # Histories that differ only outside this transition have the same score.
    # Summing their posterior mass first avoids calculating the same n-vector
    # twice for each target state in the four-wave likelihood.
    for (target in 0:1) {
      mask <- h[, tt] == lag & h[, tt + 1L] == origin &
        h[, tt + 2L] == target
      posterior_mass <- as.vector(post %*% as.numeric(mask))
      scalar <- dnorm(eta_t) * (target - q) / (q * (1 - q)) * posterior_mass
      cols <- if (!origin) off0 + seq_along(beta) else off1 + seq_along(beta)
      scores[, cols] <- scores[, cols, drop = FALSE] + X * scalar
    }
  }
  delta_cols <- 3L + L$p0 + L$p1 + seq_len(L$pd)
  for (tt in 1:4) {
    mismatch_prob <- rowSums(post * (data$y[,tt] !=
      matrix(h[,tt], n, H, byrow = TRUE)))
    pi <- pi_wave[[tt]]; s <- 2 * pi; dpi <- .5 * s * (1 - s)
    scalar <- dpi * (mismatch_prob / pi - (1 - mismatch_prob) / (1 - pi))
    scores[, delta_cols] <- scores[, delta_cols, drop = FALSE] + data$Z[[tt]] * scalar
  }
  list(loglik_i = loglik_i, params = p, scores = scores, posterior = post,
       q0 = q0, q1 = q1, pi_wave = pi_wave)
}

fit_ar2_set4_reliability <- function(data, start, maxit = 1500L, reltol = 1e-10) {
  wt <- sum(data$weight); cache_eta <- cache <- NULL
  evaluate <- function(z, scores = FALSE) {
    if (is.null(cache_eta) || !identical(as.numeric(z), cache_eta) ||
        (scores && is.null(cache$scores))) {
      cache_eta <<- as.numeric(z)
      cache <<- .ar2r_components(z, data, retain_scores = scores,
                                 retain_posterior = scores)
    }
    cache
  }
  objective <- function(z) -sum(data$weight * evaluate(z)$loglik_i) / wt
  gradient <- function(z) -colSums(evaluate(z, TRUE)$scores * data$weight) / wt
  opt <- optim(start, objective, gradient, method = "BFGS",
    control = list(maxit = maxit, reltol = reltol))
  detail <- .ar2r_components(opt$par, data, TRUE, TRUE)
  score <- colSums(detail$scores * data$weight) / wt
  list(params = detail$params, eta = setNames(opt$par, ar2_reliability_names(data)),
    loglik = sum(data$weight * detail$loglik_i), converged = opt$convergence == 0L,
    optimizer_code = opt$convergence, iterations = unname(opt$counts[["function"]]),
    max_abs_score = max(abs(score)), scores = detail$scores,
    posterior = detail$posterior, n_obs = data$n_original, n_cells = data$n_cells,
    stage = data$stage)
}

ar2_reliability_quantities <- function(eta, data, detail = NULL) {
  d <- detail %||% .ar2r_components(eta, data, retain_posterior = TRUE)
  h <- latent_histories_ar2(); w <- data$weight; post <- d$posterior
  rate <- function(origin, lag, exit = FALSE) {
    num <- den <- numeric(nrow(data$y))
    for (tt in 1:2) {
      mask <- h[,tt] == lag & h[,tt+1L] == origin
      risk <- as.vector(post %*% as.numeric(mask))
      q <- if (!origin) d$q0[[tt]][[lag + 1L]] else d$q1[[tt]][[lag + 1L]]
      if (exit) q <- 1 - q
      num <- num + risk * q; den <- den + risk
    }
    sum(w * num) / sum(w * den)
  }
  c(entry_00 = rate(0, 0), entry_10 = rate(0, 1),
    exit_01 = rate(1, 0, TRUE), exit_11 = rate(1, 1, TRUE),
    mean_misclassification = sum(w * Reduce(`+`, d$pi_wave)) / (4 * sum(w)),
    initial_wave2_employment = d$params$alpha["01"] + d$params$alpha["11"])
}

#' Decompose an apparent two-wave transition under the AR(2) model
#'
#' Conditions only on the two reported states that form an apparent transition.
#' The later report is integrated out, so the result is a prospective implication
#' of the fitted model rather than a retrospective classification using the full
#' observed history. The three reported events are mutually exclusive: no latent
#' transition, a latent transition that reverses next period, and a latent
#' transition that persists next period.
#'
#' @param data Collapsed reliability-stage data returned by
#'   [subset_ar2_reliability_stage()].
#' @param fit Fitted object returned by [fit_ar2_set4_reliability()].
#' @return Data frame with pooled, apparent-entry, and apparent-exit shares.
prospective_apparent_transition_ar2 <- function(data, fit) {
  h <- latent_histories_ar2()
  detail <- .ar2r_components(fit$eta, data)
  n <- nrow(data$y); H <- nrow(h)

  latent_log <- matrix(0, n, H)
  initial_pair <- paste0(h[, 1], h[, 2])
  latent_log <- sweep(latent_log, 2,
                      log(fit$params$alpha[initial_pair]), "+")
  for (tt in 1:2) for (j in seq_len(H)) {
    origin <- h[j, tt + 1L]; lag <- h[j, tt]
    q <- if (origin == 0L) detail$q0[[tt]][[lag + 1L]] else
      detail$q1[[tt]][[lag + 1L]]
    latent_log[, j] <- latent_log[, j] +
      if (h[j, tt + 2L] == 1L) log(q) else log1p(-q)
  }

  directions <- c("All apparent transitions", "Apparent entries",
                  "Apparent exits")
  numerator <- matrix(0, length(directions), 3L,
    dimnames = list(directions,
      c("classification_only", "true_reversal", "true_persistent")))
  denominator <- setNames(numeric(length(directions)), directions)

  for (tt in 1:2) {
    log_joint <- latent_log
    for (ss in tt:(tt + 1L)) for (j in seq_len(H)) {
      pi <- detail$pi_wave[[ss]]
      log_joint[, j] <- log_joint[, j] +
        ifelse(data$y[, ss] != h[j, ss], log(pi), log1p(-pi))
    }
    row_max <- Reduce(pmax, lapply(seq_len(H), function(j) log_joint[, j]))
    posterior <- exp(log_joint - row_max)
    posterior <- posterior / rowSums(posterior)
    event_masks <- list(
      classification_only = h[, tt] == h[, tt + 1L],
      true_reversal = h[, tt] != h[, tt + 1L] &
        h[, tt + 2L] == h[, tt],
      true_persistent = h[, tt] != h[, tt + 1L] &
        h[, tt + 2L] == h[, tt + 1L]
    )
    event_probability <- sapply(event_masks, function(mask)
      as.vector(posterior %*% mask))
    selections <- list(
      `All apparent transitions` = data$y[, tt] != data$y[, tt + 1L],
      `Apparent entries` = data$y[, tt] == 0L & data$y[, tt + 1L] == 1L,
      `Apparent exits` = data$y[, tt] == 1L & data$y[, tt + 1L] == 0L)
    for (direction in directions) {
      w <- data$weight * selections[[direction]]
      numerator[direction, ] <- numerator[direction, ] +
        colSums(event_probability * w)
      denominator[direction] <- denominator[direction] + sum(w)
    }
  }
  shares <- numerator / denominator
  data.frame(direction = directions, shares, row.names = NULL,
             check.names = FALSE)
}

analytical_se_ar2_set4_reliability <- function(data, fit, step = 1e-5) {
  wt <- sum(data$weight)
  objective <- function(z) -sum(data$weight * .ar2r_components(z, data)$loglik_i) / wt
  gradient <- function(z) -colSums(.ar2r_components(z, data, TRUE)$scores * data$weight) / wt
  bread <- optimHess(fit$eta, objective, gradient); bread <- (bread + t(bread)) / 2
  eig <- eigen(bread, symmetric = TRUE); threshold <- max(eig$values) * 1e-9
  inv <- solve(bread); meat <- crossprod(fit$scores * sqrt(data$weight_sq)) / wt^2
  vcov <- inv %*% meat %*% inv
  if (data$n_original > length(fit$eta))
    vcov <- vcov * data$n_original / (data$n_original - length(fit$eta))
  dimnames(vcov) <- list(names(fit$eta), names(fit$eta))
  q0 <- ar2_reliability_quantities(fit$eta, data)
  J <- matrix(NA_real_, length(q0), length(fit$eta),
              dimnames = list(names(q0), names(fit$eta)))
  for (j in seq_along(fit$eta)) {
    hh <- step * max(1, abs(fit$eta[j])); zp <- zm <- fit$eta
    zp[j] <- zp[j] + hh; zm[j] <- zm[j] - hh
    J[,j] <- (ar2_reliability_quantities(zp, data) -
                ar2_reliability_quantities(zm, data)) / (2 * hh)
  }
  qvcov <- J %*% vcov %*% t(J)
  score <- colSums(fit$scores * data$weight) / wt
  correction <- as.vector(solve(bread, score)); names(correction) <- names(fit$eta)
  list(coefficient_summary = data.frame(term = names(fit$eta),
      estimate = unname(fit$eta), se = sqrt(pmax(diag(vcov), 0))),
    probability_summary = data.frame(quantity = names(q0), estimate = unname(q0),
      se = sqrt(pmax(diag(qvcov), 0))), vcov_eta = vcov, vcov_quantities = qvcov,
    diagnostics = list(information_rank = sum(eig$values > threshold),
      parameter_count = length(fit$eta), min_information_eigenvalue = min(eig$values),
      information_condition = max(eig$values) / min(eig$values),
      max_abs_score = fit$max_abs_score,
      max_abs_newton_correction = max(abs(correction)), newton_correction = correction))
}
