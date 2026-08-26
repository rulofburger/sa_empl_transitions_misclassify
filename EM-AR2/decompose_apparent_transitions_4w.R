# Prospective decomposition of apparent two-wave employment transitions.
# Produces Appendix Table A4C inputs by comparing the Table 1 nonstationary
# no-error and symmetric-error models with the final Table 5 AR(2) model.

source("EM-baseline/R/source_all.R")
source("EM-AR1-4W/R/source_all.R")
source("EM-AR2/R/source_all.R")

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
  directions <- c("All apparent transitions", "Apparent entries",
                  "Apparent exits")
  numerator <- matrix(0, length(directions), 3L,
    dimnames = list(directions,
      c("classification_only", "true_reversal", "true_persistent")))
  denominator <- setNames(numeric(length(directions)), directions)

  for (i in seq_len(nrow(pairs))) {
    if (pairs[i, 1] == pairs[i, 2]) next
    emission <- rep(1, nrow(h))
    for (tt in 1:2) {
      if (fit$model_type == "none") {
        emission <- emission * (h[, tt] == pairs[i, tt])
      } else if (fit$model_type == "symmetric") {
        emission <- emission * ifelse(h[, tt] == pairs[i, tt],
                                      1 - fit$params$pi, fit$params$pi)
      } else stop("The appendix comparison expects none or symmetric error")
    }
    posterior <- prior * emission
    posterior <- posterior / sum(posterior)
    event_probability <- c(
      classification_only = sum(posterior[h[, 1] == h[, 2]]),
      true_reversal = sum(posterior[h[, 1] != h[, 2] & h[, 3] == h[, 1]]),
      true_persistent = sum(posterior[h[, 1] != h[, 2] & h[, 3] == h[, 2]]))
    selected <- c(`All apparent transitions` = TRUE,
      `Apparent entries` = unname(pairs[i, 1] == 0L),
      `Apparent exits` = unname(pairs[i, 1] == 1L))
    for (direction in directions) if (selected[[direction]]) {
      numerator[direction, ] <- numerator[direction, ] +
        pair_weight[i] * event_probability
      denominator[direction] <- denominator[direction] + pair_weight[i]
    }
  }
  shares <- numerator / denominator
  data.frame(direction = directions, shares, row.names = NULL,
             check.names = FALSE)
}

fit_raw <- readRDS("EM-baseline/output/results/fit_none_free.rds")
fit_adjusted <- readRDS("EM-baseline/output/results/fit_sym_free.rds")
full_data <- prepare_ar2_set4_reliability_4w()
data_final <- subset_ar2_reliability_stage(full_data, "reliability")
fit_final <- readRDS(paste0(
  "EM-AR2/output/results/set4_reliability/",
  "fit_ar2_set4_reliability_latest.rds"))

results <- rbind(
  transform(prospective_apparent_transition_baseline(fit_raw),
            model = "Table 1: raw/no error"),
  transform(prospective_apparent_transition_baseline(fit_adjusted),
            model = "Table 1: symmetric error"),
  transform(prospective_apparent_transition_ar2(data_final, fit_final),
            model = "Table 5: final AR(2)"))
results <- results[, c("direction", "model", "classification_only",
                       "true_reversal", "true_persistent")]
stopifnot(max(abs(rowSums(results[, 3:5]) - 1)) < 1e-10)

output_path <- paste0(
  "EM-AR2/output/results/set4_reliability/",
  "apparent_transition_decomposition.csv")
write.csv(results, output_path, row.names = FALSE)
print(transform(results,
  classification_only = 100 * classification_only,
  true_reversal = 100 * true_reversal,
  true_persistent = 100 * true_persistent), digits = 5, row.names = FALSE)
