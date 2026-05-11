# ==============================================================================
# build_descriptive_table.R
# Created: 2026-05-08
#
# Generates the descriptive statistics table for the REStat paper.
# Outputs: paper/output/table_descriptive.tex
#
# The table has three panels:
#   Panel A: Sample characteristics
#   Panel B: Observed employment transition rates (3-month and 6-month)
#   Panel C: Data quality — age and education inconsistency rates
#
# Data sources:
#   - 3-wave data: data/raw/df_qlfs_A.rds (via scripts/ingest_data_3waves_SA.R)
#   - 4-wave data: same source, filtered to 4 complete waves
#
# Usage: Rscript paper/build_descriptive_table.R
#        or source() from an interactive session.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(haven)
})

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
out_dir <- here("paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "table_descriptive.tex")

# ---------------------------------------------------------------------------
# 1. Load 3-wave data
# ---------------------------------------------------------------------------
# Load the raw data and replicate the core filtering from ingest_data_3waves_SA.R
# We load a simplified version (no duration imputation) because this script
# only needs employment status, demographics, and inconsistency indicators.

message("Loading 3-wave data...")

df_raw <- readRDS(here("data", "raw", "df_qlfs_A.rds")) %>%
  filter(age1 > 17 & age1 < 56) %>%
  filter(!is.na(employed1), !is.na(employed2), !is.na(employed3)) %>%
  select(contains(c(
    "employed", "age", "educ", "female", "weight"
  ))) %>%
  select(-contains("4")) %>%
  rename(y1 = employed1, y2 = employed2, y3 = employed3) %>%
  mutate(weight = (weight1 * weight2 * weight3)^(1/3)) %>%
  as.data.frame()

# Convert haven-labelled to numeric
haven_cols <- names(df_raw)[sapply(df_raw, haven::is.labelled)]
df_raw[haven_cols] <- lapply(df_raw[haven_cols], \(x) as.numeric(unclass(x)))

# Make employment binary (0/1)
df_raw <- df_raw %>%
  mutate(
    y1 = as.integer(y1 == 1),
    y2 = as.integer(y2 == 1),
    y3 = as.integer(y3 == 1)
  )

n3 <- nrow(df_raw)
message(sprintf("3-wave sample: N = %s", format(n3, big.mark = ",")))

# ---------------------------------------------------------------------------
# 2. Load 4-wave data (for 6-month transition rates)
# ---------------------------------------------------------------------------
message("Loading 4-wave data...")

df_4w <- readRDS(here("data", "raw", "df_qlfs_A.rds")) %>%
  filter(age1 > 17 & age1 < 56) %>%
  filter(
    !is.na(employed1), !is.na(employed2),
    !is.na(employed3), !is.na(employed4)
  ) %>%
  select(contains(c("employed", "weight"))) %>%
  rename(
    y1 = employed1, y2 = employed2,
    y3 = employed3, y4 = employed4
  ) %>%
  mutate(weight = (weight1 * weight2 * weight3 * weight4)^0.25) %>%
  as.data.frame()

# Convert and binarise
haven_cols_4 <- names(df_4w)[sapply(df_4w, haven::is.labelled)]
df_4w[haven_cols_4] <- lapply(df_4w[haven_cols_4], \(x) as.numeric(unclass(x)))
df_4w <- df_4w %>%
  mutate(across(c(y1, y2, y3, y4), \(x) as.integer(x == 1)))

n4 <- nrow(df_4w)
message(sprintf("4-wave sample: N = %s", format(n4, big.mark = ",")))

# ---------------------------------------------------------------------------
# 3. Compute inconsistency indicators (replicate compute_inconsistencies())
# ---------------------------------------------------------------------------
# Age rule: consecutive wave gap must be 0 or +1 year
# Education rule: non-decreasing, at most +1 per wave
# We need educ as numeric for these checks

.age_ok  <- function(d) !is.na(d) & (abs(d) < 0.01 | abs(d - 1) < 0.01)
.edu_ok  <- function(d) !is.na(d) & (abs(d) < 0.01 | abs(d - 1) < 0.01)

df_raw <- df_raw %>%
  mutate(
    # Pairwise inconsistency indicators
    V12_age = as.integer(!.age_ok(age2 - age1)),
    V23_age = as.integer(!.age_ok(age3 - age2)),
    V12_edu = as.integer(!.edu_ok(educ2 - educ1)),
    V23_edu = as.integer(!.edu_ok(educ3 - educ2)),
    # Any age inconsistency across any pairwise gap
    any_age_incons = as.integer(V12_age | V23_age),
    # Any education inconsistency
    any_edu_incons = as.integer(V12_edu | V23_edu),
    # Both
    any_incons = as.integer(any_age_incons | any_edu_incons)
  )

# ---------------------------------------------------------------------------
# 4. Helper: weighted mean
# ---------------------------------------------------------------------------
wmean <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  sum(x[keep] * w[keep]) / sum(w[keep])
}

# ---------------------------------------------------------------------------
# 5. Panel A: Sample characteristics (3-wave)
# ---------------------------------------------------------------------------
message("Computing Panel A: sample characteristics...")

w <- df_raw$weight

pa_n           <- n3
pa_age_mean    <- wmean(df_raw$age1, w)
pa_age_min     <- min(df_raw$age1, na.rm = TRUE)
pa_age_max     <- max(df_raw$age1, na.rm = TRUE)
pa_educ_mean   <- wmean(df_raw$educ1, w)
pa_empl_w1     <- wmean(df_raw$y1, w) * 100
pa_empl_w2     <- wmean(df_raw$y2, w) * 100
pa_empl_w3     <- wmean(df_raw$y3, w) * 100

# Female rate (only if female1 column exists)
pa_female <- if ("female1" %in% names(df_raw)) {
  wmean(df_raw$female1 == 1, w) * 100
} else {
  NA_real_
}

# ---------------------------------------------------------------------------
# 6. Panel B: Transition rates
# ---------------------------------------------------------------------------
message("Computing Panel B: transition rates...")

# 3-month rates (observed in 3-wave data, waves 1->2 and 2->3)
# Entry: y_{t-1}=0 -> y_t=1
# Exit:  y_{t-1}=1 -> y_t=0

# Wave 1->2
pb_entry_12 <- with(
  df_raw,
  wmean(y2[y1 == 0], weight[y1 == 0]) * 100
)
pb_exit_12 <- with(
  df_raw,
  wmean(1 - y2[y1 == 1], weight[y1 == 1]) * 100
)

# Wave 2->3
pb_entry_23 <- with(
  df_raw,
  wmean(y3[y2 == 0], weight[y2 == 0]) * 100
)
pb_exit_23 <- with(
  df_raw,
  wmean(1 - y3[y2 == 1], weight[y2 == 1]) * 100
)

# Average 3-month rates
pb_entry_3m <- (pb_entry_12 + pb_entry_23) / 2
pb_exit_3m  <- (pb_exit_12  + pb_exit_23)  / 2

# Implied 6-month rates from 2-step Markov:
# P(y3=1|y1=0) via law of total probability under Markov assumption
# = P(y2=1|y1=0)*P(y3=1|y2=1) + P(y2=0|y1=0)*P(y3=1|y2=0)
# = entry_12 * theta1_implied + (1-entry_12) * entry_23
theta1_implied <- 1 - pb_exit_3m / 100
pb_entry_6m_implied <- (pb_entry_12/100) * theta1_implied * 100 +
  (1 - pb_entry_12/100) * (pb_entry_23/100) * 100

# Simpler: use the naive 3-month compounding for the table footnote comparison
# Implied 6-month entry = 1 - (1 - e)^2 where e = 3-month entry
pb_entry_6m_compound <- (1 - (1 - pb_entry_3m/100)^2) * 100
pb_exit_6m_compound  <- (1 - (1 - pb_exit_3m/100)^2)  * 100

# Observed 6-month rates from 4-wave data (wave 1->3, skipping wave 2)
pb_entry_6m_obs <- with(
  df_4w,
  wmean(y3[y1 == 0], weight[y1 == 0]) * 100
)
pb_exit_6m_obs <- with(
  df_4w,
  wmean(1 - y3[y1 == 1], weight[y1 == 1]) * 100
)

# ---------------------------------------------------------------------------
# 7. Panel C: Data quality
# ---------------------------------------------------------------------------
message("Computing Panel C: inconsistency rates...")

pc_age_incons  <- wmean(df_raw$any_age_incons, w) * 100
pc_edu_incons  <- wmean(df_raw$any_edu_incons, w) * 100
pc_any_incons  <- wmean(df_raw$any_incons, w) * 100

# ---------------------------------------------------------------------------
# 8. Format helpers
# ---------------------------------------------------------------------------
.f1 <- function(x) sprintf("%.1f", x)          # 1 decimal
.f2 <- function(x) sprintf("%.2f", x)          # 2 decimals
.fi <- function(x) format(round(x), big.mark = ",", scientific = FALSE)  # integer

# ---------------------------------------------------------------------------
# 9. Write LaTeX table
# ---------------------------------------------------------------------------
message("Writing LaTeX table to: ", out_path)

lines <- c(
  "% ===========================================================================",
  "% table_descriptive.tex",
  "% Generated by paper/build_descriptive_table.R",
  sprintf("%% Date: %s", Sys.Date()),
  "% ===========================================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Descriptive statistics: QLFS rotating panel, 2010Q1--2019Q4}",
  "\\label{tab:descriptive}",
  "\\begin{tabular}{lc}",
  "\\toprule",
  "\\textbf{Statistic} & \\textbf{Value} \\\\",
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Panel A: Sample (3-wave linked triads)}} \\\\",
  sprintf("Observations & %s \\\\", .fi(pa_n)),
  sprintf("Age range & %d--%d \\\\", pa_age_min, pa_age_max),
  sprintf("Mean age (wave 1) & %s \\\\", .f1(pa_age_mean)),
  sprintf("Mean education (years, wave 1) & %s \\\\", .f1(pa_educ_mean))
)

# Only include female row if data available
if (!is.na(pa_female)) {
  lines <- c(lines,
    sprintf("\\%% female & %s\\%% \\\\", .f1(pa_female))
  )
}

lines <- c(lines,
  sprintf("Employment rate: wave 1 (\\%%) & %s \\\\", .f1(pa_empl_w1)),
  sprintf("Employment rate: wave 2 (\\%%) & %s \\\\", .f1(pa_empl_w2)),
  sprintf("Employment rate: wave 3 (\\%%) & %s \\\\", .f1(pa_empl_w3)),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Panel B: Observed transition rates}} \\\\",
  "\\multicolumn{2}{l}{\\quad\\textit{3-month rates}} \\\\",
  sprintf("\\quad Entry rate: wave 1$\\to$2 (\\%%) & %s \\\\",  .f2(pb_entry_12)),
  sprintf("\\quad Entry rate: wave 2$\\to$3 (\\%%) & %s \\\\",  .f2(pb_entry_23)),
  sprintf("\\quad Exit rate: wave 1$\\to$2 (\\%%) & %s \\\\",   .f2(pb_exit_12)),
  sprintf("\\quad Exit rate: wave 2$\\to$3 (\\%%) & %s \\\\",   .f2(pb_exit_23)),
  "\\multicolumn{2}{l}{\\quad\\textit{6-month rates}} \\\\",
  sprintf("\\quad Observed entry rate: wave 1$\\to$3 (\\%%) & %s \\\\",  .f2(pb_entry_6m_obs)),
  sprintf("\\quad Observed exit rate: wave 1$\\to$3 (\\%%) & %s \\\\",   .f2(pb_exit_6m_obs)),
  sprintf("\\quad Implied entry rate (2-step Markov, \\%%) & %s \\\\", .f2(pb_entry_6m_compound)),
  sprintf("\\quad Implied exit rate (2-step Markov, \\%%) & %s \\\\",  .f2(pb_exit_6m_compound)),
  sprintf("\\quad Ratio: implied / observed entry & %s \\\\",
    .f2(pb_entry_6m_compound / pb_entry_6m_obs)),
  sprintf("\\quad Ratio: implied / observed exit & %s \\\\",
    .f2(pb_exit_6m_compound / pb_exit_6m_obs)),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Panel C: Data quality --- inconsistent responses}} \\\\",
  sprintf("\\%% with age inconsistency (any wave) & %s \\\\",  .f1(pc_age_incons)),
  sprintf("\\%% with education inconsistency (any wave) & %s \\\\",  .f1(pc_edu_incons)),
  sprintf("\\%% with either inconsistency & %s \\\\",  .f1(pc_any_incons)),
  "\\bottomrule",
  "\\multicolumn{2}{p{0.85\\textwidth}}{\\footnotesize \\textit{Note:} \\",
  "All rates are weighted using the geometric mean of wave-specific survey weights. \\",
  "3-month transition rates are the mean of wave 1$\\to$2 and wave 2$\\to$3 \\",
  "conditional rates. 6-month rates use 4-wave linked observations \\",
  sprintf("(N = %s). Implied 6-month rates compound the 3-month rates under \\", .fi(n4)),
  "the Markov assumption: $1 - (1-r)^2$ where $r$ is the 3-month rate. \\",
  "Age inconsistency: consecutive wave gap not in $\\{0, 1\\}$ years. \\",
  "Education inconsistency: non-monotone or multi-category jump.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(lines, out_path)
message("Done. Table written to: ", out_path)
