# Decomposition of apparent two-wave employment transitions for Table 8.
# Except for the explicitly labelled Table 7/9 posterior decompositions, each row
# conditions only on the two status reports forming the apparent transition and
# integrates the subsequent report out of the fitted model.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
})
source("EM-baseline/R/source_all.R")
source("EM-baseline-ext/R/source_all.R")
source("EM-AR1-4W/R/source_all.R")
source("EM-AR2/R/source_all.R")
source("EM-tenure/R/source_all.R")

.decomposition_template <- function() {
  directions <- c("All apparent transitions", "Apparent entries",
                  "Apparent exits")
  list(
    numerator = matrix(0, length(directions), 3L,
      dimnames = list(directions,
        c("classification_only", "true_reversal", "true_persistent"))),
    denominator = setNames(numeric(length(directions)), directions))
}

.add_decomposition_pair <- function(acc, event_probability, y0, y1, weight) {
  selected <- list(
    `All apparent transitions` = y0 != y1,
    `Apparent entries` = y0 == 0L & y1 == 1L,
    `Apparent exits` = y0 == 1L & y1 == 0L)
  for (direction in names(selected)) {
    w <- weight * selected[[direction]]
    acc$numerator[direction, ] <- acc$numerator[direction, ] +
      colSums(event_probability * w)
    acc$denominator[direction] <- acc$denominator[direction] + sum(w)
  }
  acc
}

.finish_decomposition <- function(acc) {
  shares <- acc$numerator / acc$denominator
  data.frame(direction = rownames(shares), shares, row.names = NULL,
             check.names = FALSE)
}

# Generic prospective calculation. `log_prior` is the model-implied latent-path
# distribution before seeing reports. `pi_at(t)` may return either a common N
# vector or an N x H matrix when error varies across latent types.
.prospective_from_log_prior <- function(h, log_prior, y, weight, pair_times,
                                        pi_at) {
  n <- nrow(y); H <- nrow(h)
  if (nrow(log_prior) != n || ncol(log_prior) != H)
    stop("Latent prior has the wrong dimensions")
  acc <- .decomposition_template()
  for (tt in pair_times) {
    log_joint <- log_prior
    for (ss in tt:(tt + 1L)) {
      pi <- pi_at(ss)
      if (is.vector(pi)) pi <- matrix(pi, n, H)
      if (!identical(dim(pi), c(n, H))) stop("pi_at returned the wrong dimensions")
      mismatch <- outer(y[, ss], h[, ss], "!=")
      log_joint <- log_joint + ifelse(mismatch, log(pi), log1p(-pi))
    }
    row_max <- apply(log_joint, 1L, max)
    posterior <- exp(log_joint - row_max)
    posterior <- posterior / rowSums(posterior)
    masks <- list(
      classification_only = h[, tt] == h[, tt + 1L],
      true_reversal = h[, tt] != h[, tt + 1L] & h[, tt + 2L] == h[, tt],
      true_persistent = h[, tt] != h[, tt + 1L] &
        h[, tt + 2L] == h[, tt + 1L])
    event_probability <- sapply(masks, function(mask)
      as.vector(posterior %*% mask))
    acc <- .add_decomposition_pair(acc, event_probability,
      y[, tt], y[, tt + 1L], weight)
  }
  .finish_decomposition(acc)
}

prospective_apparent_transition_baseline <- function(fit) {
  h <- latent_histories()
  prior <- prior_over_histories(h, fit$params$theta1,
                                fit$params$theta0, fit$params$alpha)
  observed <- .baseline_observed_histories()
  pair_index <- 1L + observed[, 1] + 2L * observed[, 2]
  cell_weight <- fit$sample$cell_weights
  pair_weight <- numeric(4L)
  sums <- rowsum(cell_weight, pair_index, reorder = FALSE)
  pair_weight[as.integer(rownames(sums))] <- sums[, 1L]
  pairs <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1))
  acc <- .decomposition_template()
  for (i in seq_len(nrow(pairs))) {
    if (pairs[i, 1] == pairs[i, 2]) next
    emission <- rep(1, nrow(h))
    for (tt in 1:2) {
      if (fit$model_type == "none") emission <- emission *
        (h[, tt] == pairs[i, tt]) else emission <- emission *
        ifelse(h[, tt] == pairs[i, tt], 1 - fit$params$pi, fit$params$pi)
    }
    posterior <- prior * emission / sum(prior * emission)
    event_probability <- matrix(c(
      sum(posterior[h[, 1] == h[, 2]]),
      sum(posterior[h[, 1] != h[, 2] & h[, 3] == h[, 1]]),
      sum(posterior[h[, 1] != h[, 2] & h[, 3] == h[, 2]])), nrow = 1L)
    acc <- .add_decomposition_pair(acc, event_probability,
      pairs[i, 1], pairs[i, 2], pair_weight[i])
  }
  .finish_decomposition(acc)
}

.markov_covariate_log_prior <- function(h, X, params) {
  Xt <- .as_transition_design(X)
  n <- nrow(Xt[[1L]]); H <- nrow(h)
  q0 <- lapply(Xt, function(x) pmin(pmax(
    pnorm(as.vector(x %*% params$beta0)), 1e-9), 1 - 1e-9))
  q1 <- lapply(Xt, function(x) pmin(pmax(
    pnorm(as.vector(x %*% params$beta1)), 1e-9), 1 - 1e-9))
  out <- matrix(0, n, H)
  for (j in seq_len(H)) {
    out[, j] <- if (h[j, 1L]) log(params$alpha) else log1p(-params$alpha)
    for (tt in 1:2) {
      q <- if (h[j, tt] == 0L) q0[[tt]] else q1[[tt]]
      out[, j] <- out[, j] + if (h[j, tt + 1L]) log(q) else log1p(-q)
    }
  }
  out
}

prospective_apparent_transition_covariates <- function(df, X, params, pi_wave) {
  h <- latent_histories(); y <- as.matrix(df[paste0("y", 1:3)])
  .prospective_from_log_prior(h, .markov_covariate_log_prior(h, X, params),
    y, df$weight, 1L, function(tt) pi_wave[[tt]])
}

prospective_apparent_transition_ar2_constant <- function(df, fit) {
  h <- latent_histories_ar2(); p <- fit$params
  prior <- prior_over_histories_ar2(h, p$theta0, p$theta01,
    p$theta1, p$theta10, alpha = p$alpha)
  y <- as.matrix(df[paste0("y", 1:4)])
  log_prior <- matrix(log(prior), nrow(y), nrow(h), byrow = TRUE)
  .prospective_from_log_prior(h, log_prior, y, df$weight, 1:2,
    function(tt) rep(p$pi, nrow(y)))
}

.fmm_constant_path_prior <- function(h, p, type) {
  suffix <- c("A", "B")[type]
  q0 <- p[[paste0("theta0_", suffix)]]
  q1 <- p[[paste0("theta1_", suffix)]]
  alpha <- p[[paste0("alpha_", suffix)]]
  ans <- numeric(nrow(h))
  for (j in seq_len(nrow(h))) {
    val <- if (h[j, 1L]) alpha else 1 - alpha
    for (tt in 1:3) {
      q <- if (h[j, tt] == 0L) q0 else q1
      val <- val * if (h[j, tt + 1L]) q else 1 - q
    }
    ans[j] <- val
  }
  ans
}

prospective_apparent_transition_fmm_constant <- function(data, fit) {
  h0 <- latent_histories_ar2(); h <- rbind(h0, h0); p <- fit$params
  prior <- c(p$phi * .fmm_constant_path_prior(h0, p, 1L),
    (1 - p$phi) * .fmm_constant_path_prior(h0, p, 2L))
  .prospective_from_log_prior(h,
    matrix(log(prior), nrow(data$y), nrow(h), byrow = TRUE),
    data$y, data$weight, 1:2,
    function(tt) rep(unname(p$pi), nrow(data$y)))
}

.fmm_controlled_log_prior <- function(data, p) {
  h0 <- latent_histories_ar2(); h <- rbind(h0, h0)
  n <- nrow(data$y); out <- matrix(NA_real_, n, nrow(h))
  type_mass <- c(p$phi, 1 - p$phi)
  active1 <- data$persistence_active %||% rep(TRUE, length(data$entry_active))
  for (k in 1:2) {
    q0 <- lapply(1:3, function(tt) .fmm_covinc_transition(
      data$X[[tt]], p$beta0, k, data$entry_active))
    q1 <- lapply(1:3, function(tt) .fmm_covinc_transition(
      data$X[[tt]], p$beta1, k, active1))
    for (j in seq_len(nrow(h0))) {
      col <- (k - 1L) * nrow(h0) + j
      out[, col] <- log(type_mass[k]) +
        if (h0[j, 1L]) log(p$alpha[k]) else log1p(-p$alpha[k])
      for (tt in 1:3) {
        q <- if (h0[j, tt] == 0L) q0[[tt]] else q1[[tt]]
        out[, col] <- out[, col] +
          if (h0[j, tt + 1L]) log(q) else log1p(-q)
      }
    }
  }
  list(h = h, log_prior = out)
}

prospective_apparent_transition_fmm_controlled <- function(data, fit) {
  p <- fit$params; latent <- .fmm_controlled_log_prior(data, p)
  n <- nrow(data$y); H0 <- nrow(latent$h) / 2L
  .prospective_from_log_prior(latent$h, latent$log_prior,
    data$y, data$weight, 1:2, function(tt) {
      a <- .fmm_covinc_pi(data$Z[[tt]], p$delta, 1L)
      b <- .fmm_covinc_pi(data$Z[[tt]], p$delta, 2L)
      cbind(matrix(a, n, H0), matrix(b, n, H0))
    })
}

# Duration reports are measurement equations in Table 7 and impose consistency
# restrictions on latent paths. This row therefore uses the fitted full-record
# posterior rather than pretending it is a status-pair-only calculation.
posterior_apparent_transition_tenure <- function(df, fit) {
  h <- latent_histories(); gamma <- fit$gamma
  masks <- list(
    classification_only = h[, 1] == h[, 2],
    true_reversal = h[, 1] != h[, 2] & h[, 3] == h[, 1],
    true_persistent = h[, 1] != h[, 2] & h[, 3] == h[, 2])
  event_probability <- sapply(masks, function(mask)
    as.vector(gamma %*% mask))
  acc <- .decomposition_template()
  acc <- .add_decomposition_pair(acc, event_probability,
    df$y1, df$y2, df$weight)
  .finish_decomposition(acc)
}

# Table 1.
message("Decomposing Table 1 models")
fit_raw <- readRDS("EM-baseline/output/results/fit_none_free.rds")
fit_adjusted <- readRDS("EM-baseline/output/results/fit_sym_free.rds")
result_list <- list(
  transform(prospective_apparent_transition_baseline(fit_raw),
            model = "Table 1 col 1: raw/no error", conditioning = "pair_only"),
  transform(prospective_apparent_transition_baseline(fit_adjusted),
            model = "Table 1 col 2: symmetric error", conditioning = "pair_only"))
rm(fit_raw, fit_adjusted); gc(verbose = FALSE)

# Tables 4 and 7 use the common three-wave ingestion.
message("Preparing the Table 4 and Table 7 samples")
source("scripts/ingest_data_3waves_SA.R")
sector_source <- as.data.frame(readRDS("data/raw/QLFSmerged_mapped.rds"))
df_qlfs <- attach_transition_informal_sector(df_qlfs, sector_source)
keep <- complete.cases(df_qlfs[, c("y1", "y2", "y3", "weight",
  "age1", "age2", "age3", "educ1", "educ2", "educ3")]) &
  df_qlfs$weight > 0
df4 <- as.data.frame(df_qlfs[keep, , drop = FALSE])
for (nm in c("y1", "y2", "y3")) df4[[nm]] <- as.integer(df4[[nm]])
df4$weight <- as.numeric(df4$weight)
for (nm in c("contracttype1", "contracttype2"))
  df4[[nm]] <- ifelse(is.na(df4[[nm]]), 0L, as.integer(df4[[nm]]))
X4 <- prepare_covariate_matrix(df4, covariate_set = 4L)$X_transition
inc <- add_inconsistency_count_dummies(
  compute_demographic_inconsistencies(compute_inconsistency_extent(df4)))
Z4 <- lapply(1:3, function(tt) cbind(
  error_intercept = 1,
  age_inconsistency = inc[[paste0("Y_age_", tt)]],
  education_inconsistency = inc[[paste0("Y_edu_", tt)]],
  race_inconsistency = inc[[paste0("Y_race_", tt)]],
  gender_inconsistency = inc[[paste0("Y_gender_", tt)]],
  two_inconsistencies = inc[[paste0("Y_exactly_2_", tt)]],
  three_inconsistencies = inc[[paste0("Y_exactly_3_", tt)]],
  four_inconsistencies = inc[[paste0("Y_exactly_4_", tt)]]))
fit_t4_constant <- readRDS(
  "EM-baseline-ext/output/results/fit_cov_s4_sym_free.rds")
fit_t4_reliability <- readRDS(
  "EM-baseline-ext/output/results/fit_cov_s4_reliability_free.rds")
pi4_constant <- lapply(1:3, function(tt)
  rep(fit_t4_constant$params$pi, nrow(df4)))
pi4_reliability <- lapply(Z4, function(z) pmin(pmax(
  .5 * plogis(as.vector(z %*% fit_t4_reliability$params$delta)), 1e-9),
  .5 - 1e-9))

message("Decomposing Table 4 columns 4 and 5")
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_covariates(
    df4, X4, fit_t4_constant$params, pi4_constant),
  model = "Table 4 col 4: Set 2, constant error", conditioning = "pair_only")
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_covariates(
    df4, X4, fit_t4_reliability$params, pi4_reliability),
  model = "Table 4 col 5: Set 2, reliability error", conditioning = "pair_only")
rm(df4, X4, Z4, fit_t4_constant, fit_t4_reliability,
   pi4_constant, pi4_reliability, sector_source, inc); gc(verbose = FALSE)

message("Decomposing Table 7 column 4 using its full duration posterior")
df_tenure <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
fit_t7 <- readRDS(
  "EM-tenure/output/results/timegap_contamination_robustness/fits_latest.rds")$marginal
result_list[[length(result_list) + 1L]] <- transform(
  posterior_apparent_transition_tenure(df_tenure, fit_t7),
  model = "Table 7 col 4: independent duration error",
  conditioning = "full_duration_posterior")
rm(df_qlfs, df_tenure, fit_t7); gc(verbose = FALSE)

# Table 5 column (2), on its own complete four-wave estimation sample.
message("Decomposing Table 5 column 2")
source("scripts/ingest_data_4waves_SA.R")
df_ar2 <- df_qlfs[, c(paste0("y", 1:4), "weight")]
df_ar2$weight <- nrow(df_ar2) * df_ar2$weight / sum(df_ar2$weight)
ar2_files <- list.files("EM-AR2/output/results",
  pattern = "^em_ar2_sym_[0-9_]+[.]rds$", full.names = TRUE)
fit_t5_ar2 <- readRDS(tail(sort(ar2_files), 1L))
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_ar2_constant(df_ar2, fit_t5_ar2),
  model = "Table 5 col 2: AR(2)", conditioning = "pair_only")
rm(df_qlfs, df_ar2, fit_t5_ar2, ar2_files); gc(verbose = FALSE)

# Table 5 column (6).
message("Decomposing Table 5 column 6")
data_t5_final <- subset_ar2_reliability_stage(
  prepare_ar2_set4_reliability_4w(), "table3_column1")
fit_t5_final <- readRDS(paste0(
  "EM-AR2/output/results/set4_reliability/",
  "fit_ar2_set4_table3_column1_latest.rds"))
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_ar2(data_t5_final, fit_t5_final),
  model = "Table 5 col 6: AR(2), reliability error", conditioning = "pair_only")
rm(data_t5_final, fit_t5_final); gc(verbose = FALSE)

# Table 6 columns (3)/(4) and (7)/(8), both on the common controlled sample.
message("Decomposing Table 6 columns 3/4 and 7/8")
data_t6 <- prepare_fmm_covariates_inconsistency_4w(
  error_design = "intercept_only")
fit_t6_common <- readRDS(paste0(
  "EM-AR1-4W/output/results/",
  "fmm_ar1_4w_sym_free_same_sample_comparator.rds"))
fit_t6_type_intercepts <- readRDS(paste0(
  "EM-AR1-4W/output/results/",
  "fmm_type_error_intercept_only_4w_latest.rds"))
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_fmm_constant(data_t6, fit_t6_common),
  model = "Table 6 cols 3/4: FMM, common error", conditioning = "pair_only")
result_list[[length(result_list) + 1L]] <- transform(
  prospective_apparent_transition_fmm_controlled(data_t6, fit_t6_type_intercepts),
  model = "Table 6 cols 7/8: FMM, type error", conditioning = "pair_only")
source("EM-tenure/R/four_wave_duration_implications.R")
result_list[[length(result_list) + 1L]] <- build_duration_implications_4w()
results <- do.call(rbind, result_list)
results <- results[, c("direction", "model", "conditioning",
  "classification_only", "true_reversal", "true_persistent")]
stopifnot(
  nrow(results) == 33L,
  nrow(unique(results[c("direction", "model")])) == 33L,
  setequal(unique(results$direction),
    c("All apparent transitions", "Apparent entries", "Apparent exits")),
  max(abs(rowSums(results[, 4:6]) - 1)) < 1e-10,
  all(unlist(results[, 4:6]) >= 0),
  all(unlist(results[, 4:6]) <= 1))

output_path <- paste0(
  "EM-AR2/output/results/set4_reliability/",
  "apparent_transition_decomposition.csv")
write.csv(results, output_path, row.names = FALSE)
print(transform(results,
  classification_only = 100 * classification_only,
  true_reversal = 100 * true_reversal,
  true_persistent = 100 * true_persistent), digits = 5, row.names = FALSE)
