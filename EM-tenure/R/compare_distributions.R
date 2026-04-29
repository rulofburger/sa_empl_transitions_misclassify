# ==============================================================================
# EM-tenure: Distributional comparison helpers
# ==============================================================================
# Functions for comparing real vs simulated data across multiple dimensions.
# Used by EM-tenure/simulation_diagnostic.R.
#
# All functions accept a data frame with at minimum:
#   y1, y2, y3, tenure1-3, timegap_cat1-3, weight
# and return data.tables or lists suitable for printing/plotting.
# ==============================================================================

library(data.table)

# --- (a) Employment rates ----------------------------------------------------

#' Compute employment rates at each wave
#'
#' @param df Data frame with y1, y2, y3 and weight columns.
#' @param label Character label for this dataset.
#' @return data.table with one row per wave.
compute_employment_rates <- function(df, label = "data") {
  dt <- as.data.table(df)
  data.table(
    dataset = label,
    wave    = 1:3,
    emp_rate = c(
      weighted.mean(dt$y1, dt$weight, na.rm = TRUE),
      weighted.mean(dt$y2, dt$weight, na.rm = TRUE),
      weighted.mean(dt$y3, dt$weight, na.rm = TRUE)
    ),
    n = nrow(dt)
  )
}


# --- (b) Transition matrices -------------------------------------------------

#' Compute weighted 2x2 transition matrix between consecutive waves
#'
#' @param df Data frame with y columns and weight.
#' @param from,to Column names (e.g., "y1", "y2").
#' @param label Character label.
#' @return data.table with columns: dataset, from, to, n_weighted, row_pct.
compute_transitions <- function(df, from = "y1", to = "y2", label = "data") {
  dt <- as.data.table(df)
  trans <- dt[, .(n_weighted = sum(weight)), by = c(from, to)]
  setnames(trans, c(from, to), c("from_state", "to_state"))
  trans[, row_total := sum(n_weighted), by = from_state]
  trans[, row_pct := round(n_weighted / row_total * 100, 2)]
  trans[, dataset := label]
  trans[, waves := paste0(from, "->", to)]
  setcolorder(trans, c("dataset", "waves", "from_state", "to_state",
                        "n_weighted", "row_total", "row_pct"))
  trans[order(from_state, to_state)]
}


# --- (c) Tenure distribution --------------------------------------------------

#' Compute weighted tenure quantiles for employed observations
#'
#' @param df Data frame.
#' @param label Character label.
#' @return data.table with summary statistics.
compute_tenure_stats <- function(df, label = "data") {
  dt <- as.data.table(df)
  # Pool tenure from all waves where y == 1
  tenure_vals <- c(
    dt[y1 == 1, tenure1],
    dt[y2 == 1, tenure2],
    dt[y3 == 1, tenure3]
  )
  tenure_wts <- c(
    dt[y1 == 1, weight],
    dt[y2 == 1, weight],
    dt[y3 == 1, weight]
  )

  # Weighted quantiles via Hmisc-free approach (sort + cumulative weight)
  ord <- order(tenure_vals)
  tv  <- tenure_vals[ord]
  tw  <- tenure_wts[ord]
  cw  <- cumsum(tw) / sum(tw)
  wq  <- function(p) tv[which.max(cw >= p)]

  data.table(
    dataset  = label,
    n_obs    = length(tenure_vals),
    mean     = round(weighted.mean(tenure_vals, tenure_wts), 3),
    median   = round(wq(0.50), 3),
    q25      = round(wq(0.25), 3),
    q75      = round(wq(0.75), 3),
    p90      = round(wq(0.90), 3),
    p95      = round(wq(0.95), 3),
    p99      = round(wq(0.99), 3),
    max      = round(max(tenure_vals), 3),
    sd       = round(sqrt(weighted.mean((tenure_vals - weighted.mean(tenure_vals, tenure_wts))^2,
                                         tenure_wts)), 3),
    skewness = round({
      m <- weighted.mean(tenure_vals, tenure_wts)
      s <- sqrt(weighted.mean((tenure_vals - m)^2, tenure_wts))
      weighted.mean(((tenure_vals - m) / s)^3, tenure_wts)
    }, 3)
  )
}

#' Extract pooled tenure values and weights for employed obs (for plotting)
#'
#' @param df Data frame.
#' @param label Character label.
#' @return data.table with columns: tenure, weight, dataset.
extract_tenure_pool <- function(df, label = "data") {
  dt <- as.data.table(df)
  data.table(
    tenure  = c(dt[y1 == 1, tenure1], dt[y2 == 1, tenure2], dt[y3 == 1, tenure3]),
    weight  = c(dt[y1 == 1, weight],  dt[y2 == 1, weight],  dt[y3 == 1, weight]),
    dataset = label
  )
}


# --- (d) Timegap category distribution ---------------------------------------

#' Compute timegap category proportions for nonemployed observations
#'
#' @param df Data frame with timegap_cat1-3, y1-3, weight.
#' @param label Character label.
#' @return data.table with columns: dataset, cat, n_weighted, pct.
compute_timegap_cat_dist <- function(df, label = "data") {
  dt <- as.data.table(df)
  cats <- c(
    dt[y1 == 0, timegap_cat1],
    dt[y2 == 0, timegap_cat2],
    dt[y3 == 0, timegap_cat3]
  )
  wts <- c(
    dt[y1 == 0, weight],
    dt[y2 == 0, weight],
    dt[y3 == 0, weight]
  )
  tab <- data.table(cat = cats, weight = wts)
  tab <- tab[!is.na(cat), .(n_weighted = sum(weight)), by = cat]
  tab[, pct := round(n_weighted / sum(n_weighted) * 100, 2)]
  tab[, dataset := label]
  setcolorder(tab, c("dataset", "cat", "n_weighted", "pct"))
  tab[order(cat)]
}


# --- (e) Tenure increments ---------------------------------------------------

#' Compute tenure increment statistics for emp->emp pairs
#'
#' Increments: (tenure_t - tenure_{t-1} - 0.25). Under the model, these
#' should be ~N(0, 2*sigma2_g).
#'
#' @param df Data frame.
#' @param label Character label.
#' @return List with $stats (data.table) and $increments (numeric vector).
compute_tenure_increments <- function(df, label = "data") {
  dt <- as.data.table(df)
  inc_12 <- dt[y1 == 1 & y2 == 1, tenure2 - tenure1 - 0.25]
  inc_23 <- dt[y2 == 1 & y3 == 1, tenure3 - tenure2 - 0.25]
  inc_all <- c(inc_12, inc_23)

  stats <- data.table(
    dataset         = label,
    n_pairs         = length(inc_all),
    mean            = round(mean(inc_all, na.rm = TRUE), 4),
    sd              = round(sd(inc_all, na.rm = TRUE), 4),
    implied_sigma_g = round(sd(inc_all, na.rm = TRUE) / sqrt(2), 4),
    median          = round(median(inc_all, na.rm = TRUE), 4),
    p05             = round(quantile(inc_all, 0.05, na.rm = TRUE), 4),
    p95             = round(quantile(inc_all, 0.95, na.rm = TRUE), 4)
  )

  list(stats = stats, increments = inc_all)
}


# --- (f) Timegap category transitions ----------------------------------------

#' Compute observed 7x7 timegap category transition matrix (nonemp->nonemp)
#'
#' @param df Data frame with timegap_cat1-3, y1-3.
#' @param label Character label.
#' @return data.table in wide format (rows = from_cat, cols = to_cat pct).
compute_timegap_transitions <- function(df, label = "data") {
  dt <- as.data.table(df)
  pairs <- rbind(
    dt[y1 == 0 & y2 == 0, .(from = timegap_cat1, to = timegap_cat2)],
    dt[y2 == 0 & y3 == 0, .(from = timegap_cat2, to = timegap_cat3)]
  )
  pairs <- pairs[!is.na(from) & !is.na(to)]

  if (nrow(pairs) == 0) {
    return(data.table(dataset = label, note = "no nonemp-nonemp pairs"))
  }

  counts <- pairs[, .N, by = .(from, to)]
  counts[, row_total := sum(N), by = from]
  counts[, pct := round(N / row_total * 100, 1)]

  wide <- dcast(counts, from ~ to, value.var = "pct", fill = 0)
  wide[, dataset := label]
  setcolorder(wide, c("dataset", "from"))
  wide
}


# --- (g) Three-wave sequence distribution ------------------------------------

#' Compute distribution of 8 binary patterns (y1, y2, y3)
#'
#' @param df Data frame with y1, y2, y3, weight.
#' @param label Character label.
#' @return data.table with columns: dataset, pattern, n_weighted, pct.
compute_sequence_dist <- function(df, label = "data") {
  dt <- as.data.table(df)
  dt[, pattern := paste0(y1, y2, y3)]
  tab <- dt[, .(n_weighted = sum(weight)), by = pattern]
  tab[, pct := round(n_weighted / sum(n_weighted) * 100, 2)]
  tab[, dataset := label]
  setcolorder(tab, c("dataset", "pattern", "n_weighted", "pct"))
  tab[order(pattern)]
}


# --- Divergence metrics -------------------------------------------------------

#' Compute Jensen-Shannon divergence between two discrete distributions
#'
#' @param p,q Numeric vectors of probabilities (same length, sum to 1).
#' @return Scalar JSD in [0, log(2)].
jsd <- function(p, q) {
  # Normalise

  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2
  # KL(p||m) + KL(q||m), with 0*log(0) = 0
  kl <- function(a, b) {
    idx <- a > 0 & b > 0
    sum(a[idx] * log(a[idx] / b[idx]))
  }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}


#' Compute summary divergence table between real and simulated data
#'
#' @param real_df,sim_df Data frames.
#' @param real_label,sim_label Character labels.
#' @return data.table with one row per comparison dimension.
compute_divergence_summary <- function(real_df, sim_df,
                                       real_label = "Real",
                                       sim_label = "Simulated") {
  rows <- list()

  # Employment rate difference (avg across waves)
  er_real <- compute_employment_rates(real_df, real_label)
  er_sim  <- compute_employment_rates(sim_df, sim_label)
  rows[[1]] <- data.table(
    dimension = "Employment rate (mean abs diff)",
    value     = round(mean(abs(er_real$emp_rate - er_sim$emp_rate)), 4)
  )

  # Transition probability difference (P(y2=1|y1=1) and P(y2=1|y1=0))
  tr_real <- compute_transitions(real_df, label = real_label)
  tr_sim  <- compute_transitions(sim_df, label = sim_label)
  # row_pct for (from=1, to=1) = theta1_naive
  t1_real <- tr_real[from_state == 1 & to_state == 1, row_pct]
  t1_sim  <- tr_sim[from_state == 1 & to_state == 1, row_pct]
  t0_real <- tr_real[from_state == 0 & to_state == 1, row_pct]
  t0_sim  <- tr_sim[from_state == 0 & to_state == 1, row_pct]
  rows[[2]] <- data.table(
    dimension = "Theta1 naive (ppt diff)",
    value     = round(abs(t1_real - t1_sim), 2)
  )
  rows[[3]] <- data.table(
    dimension = "Theta0 naive (ppt diff)",
    value     = round(abs(t0_real - t0_sim), 2)
  )

  # Sequence JSD
  seq_real <- compute_sequence_dist(real_df, real_label)
  seq_sim  <- compute_sequence_dist(sim_df, sim_label)
  merged <- merge(seq_real, seq_sim, by = "pattern", all = TRUE,
                  suffixes = c("_real", "_sim"))
  merged[is.na(pct_real), pct_real := 0]
  merged[is.na(pct_sim),  pct_sim  := 0]
  rows[[4]] <- data.table(
    dimension = "Sequence pattern JSD",
    value     = round(jsd(merged$pct_real, merged$pct_sim), 6)
  )

  # Timegap category JSD
  tg_real <- compute_timegap_cat_dist(real_df, real_label)
  tg_sim  <- compute_timegap_cat_dist(sim_df, sim_label)
  merged_tg <- merge(tg_real, tg_sim, by = "cat", all = TRUE,
                     suffixes = c("_real", "_sim"))
  merged_tg[is.na(pct_real), pct_real := 0]
  merged_tg[is.na(pct_sim),  pct_sim  := 0]
  rows[[5]] <- data.table(
    dimension = "Timegap category JSD",
    value     = round(jsd(merged_tg$pct_real, merged_tg$pct_sim), 6)
  )

  # Tenure: abs diff in mean and median
  ten_real <- compute_tenure_stats(real_df, real_label)
  ten_sim  <- compute_tenure_stats(sim_df, sim_label)
  rows[[6]] <- data.table(
    dimension = "Tenure mean (abs diff, yrs)",
    value     = round(abs(ten_real$mean - ten_sim$mean), 3)
  )
  rows[[7]] <- data.table(
    dimension = "Tenure median (abs diff, yrs)",
    value     = round(abs(ten_real$median - ten_sim$median), 3)
  )

  # Tenure increment: abs diff in mean and SD
  ti_real <- compute_tenure_increments(real_df, real_label)$stats
  ti_sim  <- compute_tenure_increments(sim_df, sim_label)$stats
  rows[[8]] <- data.table(
    dimension = "Tenure increment mean (abs diff)",
    value     = round(abs(ti_real$mean - ti_sim$mean), 4)
  )
  rows[[9]] <- data.table(
    dimension = "Tenure increment SD (abs diff)",
    value     = round(abs(ti_real$sd - ti_sim$sd), 4)
  )

  out <- rbindlist(rows)
  out[, comparison := sim_label]
  setcolorder(out, c("comparison", "dimension", "value"))
  out
}
