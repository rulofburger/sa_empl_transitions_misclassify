# Audit the reported current-job start dates underlying the QLFS tenure
# variable.  This is a descriptive diagnostic only: it does not alter the
# estimation data or any fitted model.

args <- commandArgs(trailingOnly = TRUE)
panel_file <- if (length(args)) args[1] else
  Sys.getenv("QLFS_PANEL_FILE", "data/raw/df_qlfs_A.rds")
if (!file.exists(panel_file))
  stop("Run from the project root; panel file was not found: ", panel_file)
panel_name <- tools::file_path_sans_ext(basename(panel_file))
outroot <- "diagnostics/output/job_start_dates"
outdir <- file.path(outroot, panel_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

panel <- readRDS(panel_file)
needed <- unlist(lapply(1:4, function(t)
  paste0(c("period", "employed", "tenure", "weight"), t)))
needed <- c("hhnr", "pnr", "age1", needed)
if (!all(needed %in% names(panel)))
  stop("Panel is missing required columns: ",
       paste(setdiff(needed, names(panel)), collapse = ", "))

# Retain the paper's age range.  Pair-specific calculations below require
# observed employment and tenure only for the waves used in that pair.
panel <- panel[is.finite(panel$age1) & panel$age1 > 17 & panel$age1 < 56, ]
panel$row_id <- seq_len(nrow(panel))

nominal_interview <- function(period) {
  period <- as.integer(period)
  year <- 2008L + (period - 1L) %/% 4L
  quarter <- (period - 1L) %% 4L + 1L
  month <- 3L * quarter
  # A one-based serial month: 12*year + month.
  list(year = year, month = month, serial = 12L * year + month)
}

records <- do.call(rbind, lapply(1:4, function(t) {
  interview <- nominal_interview(panel[[paste0("period", t)]])
  tenure <- as.numeric(panel[[paste0("tenure", t)]])
  employed <- as.numeric(panel[[paste0("employed", t)]])
  weight <- as.numeric(panel[[paste0("weight", t)]])
  valid <- employed == 1 & is.finite(tenure) & tenure >= 0 &
    abs(tenure - round(tenure)) < 1e-8
  start_serial <- rep(NA_integer_, length(tenure))
  start_serial[valid] <- interview$serial[valid] - as.integer(round(tenure[valid]))
  start_month <- ifelse(valid, (start_serial - 1L) %% 12L + 1L, NA_integer_)
  start_year <- ifelse(valid, (start_serial - start_month) %/% 12L, NA_integer_)
  data.frame(
    row_id = panel$row_id, hhnr = panel$hhnr, pnr = panel$pnr,
    panel_wave = t, period = as.integer(panel[[paste0("period", t)]]),
    interview_year = interview$year, interview_month = interview$month,
    interview_serial = interview$serial, employed = employed,
    tenure_months = tenure, weight = weight, valid_start = valid,
    start_serial = start_serial, start_year = start_year,
    start_month = start_month
  )
}))

weighted_share <- function(condition, weight) {
  use <- !is.na(condition) & is.finite(weight) & weight > 0
  if (!any(use)) return(NA_real_)
  sum(weight[use] * condition[use]) / sum(weight[use])
}

valid_records <- records[records$valid_start, ]
month_n <- table(factor(valid_records$start_month, levels = 1:12))
month_weight <- tapply(valid_records$weight, factor(valid_records$start_month,
  levels = 1:12), sum, na.rm = TRUE)
month_table <- data.frame(
  start_month = 1:12,
  month_name = month.abb,
  observations = as.integer(month_n),
  unweighted_share = as.integer(month_n) / sum(month_n),
  weighted_share = as.numeric(month_weight) / sum(month_weight)
)

# Consecutive employed-wave pairs.  A reported start after the preceding
# nominal interview month can represent a genuine job-to-job move.  A start
# in the preceding interview month is separated because actual interview days
# are unavailable and its interpretation is therefore ambiguous.
pairs <- do.call(rbind, lapply(1:3, function(t) {
  a <- records[records$panel_wave == t, ]
  b <- records[records$panel_wave == t + 1L, ]
  stopifnot(identical(a$row_id, b$row_id))
  use <- a$valid_start & b$valid_start
  change <- b$start_serial - a$start_serial
  category <- rep(NA_character_, nrow(a))
  category[use & change == 0] <- "Same reported start month"
  category[use & change != 0 & b$start_serial > a$interview_serial] <-
    "Plausible job-to-job reset"
  category[use & change != 0 & b$start_serial == a$interview_serial] <-
    "Boundary-month reset"
  residual <- use & is.na(category)
  category[residual & abs(change) <= 3] <- "Revision: 1-3 months"
  category[residual & abs(change) >= 4 & abs(change) <= 12] <-
    "Revision: 4-12 months"
  category[residual & abs(change) > 12] <- "Revision: over 12 months"
  data.frame(
    row_id = a$row_id, pair = paste0(t, "-", t + 1L),
    previous_interview_serial = a$interview_serial,
    previous_start_serial = a$start_serial,
    current_start_serial = b$start_serial,
    start_change_months = ifelse(use, change, NA_integer_),
    category = category,
    weight = sqrt(pmax(a$weight, 0) * pmax(b$weight, 0))
  )
}))
pairs_valid <- pairs[!is.na(pairs$category) & is.finite(pairs$weight) &
  pairs$weight > 0, ]

category_order <- c(
  "Same reported start month", "Plausible job-to-job reset",
  "Boundary-month reset", "Revision: 1-3 months",
  "Revision: 4-12 months", "Revision: over 12 months")
pair_category <- do.call(rbind, lapply(split(pairs_valid, pairs_valid$pair),
  function(z) {
    n <- table(factor(z$category, levels = category_order))
    w <- tapply(z$weight, factor(z$category, levels = category_order),
      sum, na.rm = TRUE)
    data.frame(pair = z$pair[1], category = category_order,
      observations = as.integer(n),
      unweighted_share = as.integer(n) / sum(n),
      weighted_share = as.numeric(w) / sum(w))
  }))
pair_category_all <- within(do.call(rbind, list(pair_category,
  transform(data.frame(
    pair = "All",
    category = category_order,
    observations = as.integer(table(factor(pairs_valid$category,
      levels = category_order))),
    unweighted_share = as.integer(table(factor(pairs_valid$category,
      levels = category_order))) / nrow(pairs_valid),
    weighted_share = as.numeric(tapply(pairs_valid$weight,
      factor(pairs_valid$category, levels = category_order), sum,
      na.rm = TRUE)) / sum(pairs_valid$weight)
  )))), row.names <- NULL)

change_freq <- aggregate(weight ~ start_change_months, pairs_valid, sum)
change_n <- aggregate(row_id ~ start_change_months, pairs_valid, length)
names(change_n)[2] <- "observations"
change_freq <- merge(change_n, change_freq, by = "start_change_months")
change_freq$weighted_share <- change_freq$weight / sum(change_freq$weight)
change_freq <- change_freq[order(change_freq$start_change_months), ]

# Does a changed date persist into the following wave?  Persistence is much
# more consistent with a new job or a stable recalled start date than with an
# independent marginal draw at every wave.
triples <- do.call(rbind, lapply(1:2, function(t) {
  p1 <- pairs[pairs$pair == paste0(t, "-", t + 1L), ]
  p2 <- pairs[pairs$pair == paste0(t + 1L, "-", t + 2L), ]
  stopifnot(identical(p1$row_id, p2$row_id))
  use <- !is.na(p1$category) & !is.na(p2$category)
  data.frame(
    triple = paste0(t, "-", t + 1L, "-", t + 2L),
    initial_category = p1$category[use],
    next_report_same_as_middle =
      p2$current_start_serial[use] == p1$current_start_serial[use],
    weight = sqrt(p1$weight[use] * p2$weight[use])
  )
}))
triple_persistence <- do.call(rbind, lapply(
  split(triples, interaction(triples$triple, triples$initial_category,
    drop = TRUE)), function(z) data.frame(
      triple = z$triple[1], initial_category = z$initial_category[1],
      observations = nrow(z),
      next_report_same_unweighted = mean(z$next_report_same_as_middle),
      next_report_same_weighted = weighted_share(
        z$next_report_same_as_middle, z$weight))))
triple_persistence <- triple_persistence[
  order(triple_persistence$triple,
    match(triple_persistence$initial_category, category_order)), ]

summary_table <- data.frame(
  panel = panel_name,
  quantity = c(
    "Valid employed-wave start dates",
    "January share of reported start months",
    "Round-year share (year ending 0 or 5)",
    "Consecutive employed-wave pairs",
    "Same reported start month",
    "Plausible job-to-job reset",
    "Any revised pre-existing start date",
    "Revisions that are exact whole-year shifts",
    "Revised date repeated in the next employed wave",
    "Plausible reset repeated in the next employed wave"),
  estimate = c(
    nrow(valid_records),
    weighted_share(valid_records$start_month == 1, valid_records$weight),
    weighted_share(valid_records$start_year %% 5 == 0, valid_records$weight),
    nrow(pairs_valid),
    weighted_share(pairs_valid$category == "Same reported start month",
      pairs_valid$weight),
    weighted_share(pairs_valid$category == "Plausible job-to-job reset",
      pairs_valid$weight),
    weighted_share(grepl("^Revision", pairs_valid$category),
      pairs_valid$weight),
    weighted_share(abs(pairs_valid$start_change_months) %% 12 == 0,
      pairs_valid$weight * grepl("^Revision", pairs_valid$category)),
    weighted_share(triples$next_report_same_as_middle[
      grepl("^Revision", triples$initial_category)],
      triples$weight[grepl("^Revision", triples$initial_category)]),
    weighted_share(triples$next_report_same_as_middle[
      triples$initial_category == "Plausible job-to-job reset"],
      triples$weight[triples$initial_category ==
        "Plausible job-to-job reset"])),
  unit = c("count", rep("share", 2), "count", rep("share", 6))
)

write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
write.csv(month_table, file.path(outdir, "start_month_heaping.csv"),
  row.names = FALSE)
write.csv(pair_category_all,
  file.path(outdir, "consecutive_employed_pair_categories.csv"),
  row.names = FALSE)
write.csv(change_freq,
  file.path(outdir, "reported_start_month_changes.csv"), row.names = FALSE)
write.csv(triple_persistence,
  file.path(outdir, "changed_date_persistence.csv"), row.names = FALSE)

# Refresh cross-panel comparisons whenever more than one panel audit has been
# run.  This keeps the role of the A/B/C matching restrictions transparent.
summary_files <- file.path(outroot, paste0("df_qlfs_", c("A", "B", "C")),
  "summary.csv")
available_summary <- summary_files[file.exists(summary_files)]
if (length(available_summary)) {
  comparison <- do.call(rbind, lapply(available_summary, read.csv,
    stringsAsFactors = FALSE))
  write.csv(comparison, file.path(outroot, "matching_comparison.csv"),
    row.names = FALSE)
}

png(file.path(outdir, "job_start_date_audit.png"), width = 1800,
  height = 1300, res = 170)
op <- par(mfrow = c(2, 2), mar = c(5, 4.5, 3, 1))
barplot(100 * month_table$weighted_share, names.arg = month.abb,
  col = "#3973AC", border = NA, ylab = "Weighted share (%)",
  main = "Reported current-job start month")
abline(h = 100 / 12, lty = 2, col = "grey40")
plot_changes <- change_freq[abs(change_freq$start_change_months) <= 24, ]
plot(plot_changes$start_change_months, 100 * plot_changes$weighted_share,
  type = "h", lwd = 3, col = "#A23E48", xlab = "Change in reported start month",
  ylab = "Weighted share of pairs (%)",
  main = "Revisions across employed waves")
all_categories <- pair_category_all[pair_category_all$pair == "All", ]
barplot(100 * all_categories$weighted_share,
  names.arg = c("Same", "Reset", "Boundary", "1-3m", "4-12m", ">12m"),
  col = "#4C956C", border = NA, las = 2,
  ylab = "Weighted share (%)", main = "Pair classification")
persistence_plot <- aggregate(
  cbind(success = triples$weight * triples$next_report_same_as_middle,
        total = triples$weight) ~ initial_category, triples, sum)
persistence_plot$share <- persistence_plot$success / persistence_plot$total
persistence_plot <- persistence_plot[match(category_order,
  persistence_plot$initial_category, nomatch = 0), ]
barplot(100 * persistence_plot$share,
  names.arg = c("Same", "Reset", "Boundary", "1-3m", "4-12m", ">12m"),
  col = "#8C6BB1", border = NA, las = 2,
  ylab = "Next report repeats middle date (%)",
  main = "Persistence into a third employed wave")
par(op)
dev.off()

cat("\nJob-start-date audit\n")
print(summary_table, row.names = FALSE, digits = 5)
cat("\nReported start-month distribution\n")
print(month_table, row.names = FALSE, digits = 5)
cat("\nConsecutive employed-wave pair classifications\n")
print(pair_category_all, row.names = FALSE, digits = 5)
cat("\nPersistence of the middle-wave date into the next employed wave\n")
print(triple_persistence, row.names = FALSE, digits = 5)
cat("\nOutputs written to", outdir, "\n")
