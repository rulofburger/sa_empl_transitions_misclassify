# ==============================================================================
# build_descriptive_table.R
# Created: 2026-05-08
#
# Generates data-section tables for the REStat paper.
# Outputs:
#   - paper/output/table_descriptive.tex
#   - paper/output/table_tenure_diagnostics.tex
#
# Data source:
#   - data/raw/df_qlfs_A.rds
#
# Usage: Rscript paper/build_descriptive_table.R
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(haven)
})

# -----------------------------------------------------------------------------
# 0. Paths and helpers
# -----------------------------------------------------------------------------
out_dir <- here("paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

descriptive_out_path <- file.path(out_dir, "table_descriptive.tex")
tenure_out_path <- file.path(out_dir, "table_tenure_diagnostics.tex")

weighted_mean <- function(values, weights) {
  keep <- !is.na(values) & !is.na(weights) & weights > 0
  if (!any(keep)) {
    return(NA_real_)
  }

  result <- sum(values[keep] * weights[keep]) / sum(weights[keep])
  return(result)
}

weighted_percent <- function(values, weights) {
  result <- weighted_mean(as.integer(values), weights) * 100
  return(result)
}

format_1 <- function(value) {
  return(sprintf("%.1f", value))
}

format_2 <- function(value) {
  return(sprintf("%.2f", value))
}

format_int <- function(value) {
  result <- format(round(value), big.mark = ",", scientific = FALSE)
  return(result)
}

latex_bs <- intToUtf8(92L)
line_break <- paste0(" ", latex_bs, latex_bs)

latex_cmd <- function(command) {
  return(paste0(latex_bs, command))
}

latex_env <- function(environment, action = "begin") {
  return(paste0(latex_cmd(action), "{", environment, "}"))
}

latex_textbf <- function(text) {
  return(paste0(latex_cmd("textbf"), "{", text, "}"))
}

latex_textit <- function(text) {
  return(paste0(latex_cmd("textit"), "{", text, "}"))
}

latex_multicolumn <- function(n_columns, spec, text) {
  result <- paste0(
    latex_cmd("multicolumn"), "{", n_columns, "}{", spec, "}{", text, "}"
  )
  return(result)
}

latex_row <- function(...) {
  return(paste0(..., line_break))
}

convert_haven_labels <- function(data) {
  labelled_columns <- names(data)[
    vapply(data, haven::is.labelled, logical(1))
  ]

  for (column in labelled_columns) {
    set(data, j = column, value = as.numeric(unclass(data[[column]])))
  }

  return(data)
}

keep_matching_columns <- function(data, patterns) {
  keep_columns <- unique(unlist(lapply(
    patterns,
    function(pattern) grep(pattern, names(data), fixed = TRUE, value = TRUE)
  )))

  result <- data[, keep_columns, with = FALSE]
  return(result)
}

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------
message("Loading raw QLFS panel data...")

raw_data <- as.data.table(readRDS(here("data", "raw", "df_qlfs_A.rds")))
raw_data <- convert_haven_labels(raw_data)

df_three <- copy(raw_data)[
  age1 > 17 & age1 < 56 &
    !is.na(employed1) & !is.na(employed2) & !is.na(employed3)
]
df_three <- keep_matching_columns(
  df_three,
  c("employed", "age", "educ", "female", "weight", "tenure")
)

wave4_columns <- grep("4", names(df_three), fixed = TRUE, value = TRUE)
if (length(wave4_columns) > 0L) {
  df_three[, (wave4_columns) := NULL]
}

setnames(
  df_three,
  old = c("employed1", "employed2", "employed3"),
  new = c("y1", "y2", "y3")
)

df_three[, weight := (weight1 * weight2 * weight3)^(1 / 3)]
df_three[, c("y1", "y2", "y3") := lapply(
  .SD,
  function(status) as.integer(status == 1)
), .SDcols = c("y1", "y2", "y3")]

n_three <- nrow(df_three)
message(sprintf("3-wave sample: N = %s", format_int(n_three)))

df_four <- copy(raw_data)[
  age1 > 17 & age1 < 56 &
    !is.na(employed1) & !is.na(employed2) &
    !is.na(employed3) & !is.na(employed4)
]
df_four <- keep_matching_columns(df_four, c("employed", "weight"))

setnames(
  df_four,
  old = c("employed1", "employed2", "employed3", "employed4"),
  new = c("y1", "y2", "y3", "y4")
)

df_four[, weight := (weight1 * weight2 * weight3 * weight4)^0.25]
df_four[, c("y1", "y2", "y3", "y4") := lapply(
  .SD,
  function(status) as.integer(status == 1)
), .SDcols = c("y1", "y2", "y3", "y4")]

n_four <- nrow(df_four)
message(sprintf("4-wave sample: N = %s", format_int(n_four)))

# -----------------------------------------------------------------------------
# 2. Data quality indicators
# -----------------------------------------------------------------------------
age_ok <- function(delta) {
  result <- !is.na(delta) & (abs(delta) < 0.01 | abs(delta - 1) < 0.01)
  return(result)
}

education_ok <- function(delta) {
  result <- !is.na(delta) & (abs(delta) < 0.01 | abs(delta - 1) < 0.01)
  return(result)
}

df_three[, `:=`(
  v12_age = as.integer(!age_ok(age2 - age1)),
  v23_age = as.integer(!age_ok(age3 - age2)),
  v12_edu = as.integer(!education_ok(educ2 - educ1)),
  v23_edu = as.integer(!education_ok(educ3 - educ2))
)]

df_three[, `:=`(
  any_age_incons = as.integer(v12_age == 1L | v23_age == 1L),
  any_edu_incons = as.integer(v12_edu == 1L | v23_edu == 1L)
)]
df_three[, any_incons := as.integer(
  any_age_incons == 1L | any_edu_incons == 1L
)]

# -----------------------------------------------------------------------------
# 3. Descriptive table moments
# -----------------------------------------------------------------------------
message("Computing descriptive table moments...")

survey_weights <- df_three$weight

sample_n <- n_three
age_mean <- weighted_mean(df_three$age1, survey_weights)
age_min <- min(df_three$age1, na.rm = TRUE)
age_max <- max(df_three$age1, na.rm = TRUE)
education_mean <- weighted_mean(df_three$educ1, survey_weights)
employment_wave1 <- weighted_mean(df_three$y1, survey_weights) * 100
employment_wave2 <- weighted_mean(df_three$y2, survey_weights) * 100
employment_wave3 <- weighted_mean(df_three$y3, survey_weights) * 100

female_share <- if ("female1" %in% names(df_three)) {
  weighted_mean(df_three$female1 == 1, survey_weights) * 100
} else {
  NA_real_
}

entry_12 <- df_three[y1 == 0L, weighted_mean(y2, weight) * 100]
entry_23 <- df_three[y2 == 0L, weighted_mean(y3, weight) * 100]
exit_12 <- df_three[y1 == 1L, weighted_mean(1 - y2, weight) * 100]
exit_23 <- df_three[y2 == 1L, weighted_mean(1 - y3, weight) * 100]

entry_3m <- (entry_12 + entry_23) / 2
exit_3m <- (exit_12 + exit_23) / 2

entry_6m_compound <- (1 - (1 - entry_3m / 100)^2) * 100
exit_6m_compound <- (1 - (1 - exit_3m / 100)^2) * 100

entry_6m_observed <- df_four[y1 == 0L, weighted_mean(y3, weight) * 100]
exit_6m_observed <- df_four[y1 == 1L, weighted_mean(1 - y3, weight) * 100]

age_incons <- weighted_mean(df_three$any_age_incons, survey_weights) * 100
education_incons <- weighted_mean(df_three$any_edu_incons, survey_weights) * 100
any_incons <- weighted_mean(df_three$any_incons, survey_weights) * 100

# -----------------------------------------------------------------------------
# 4. Tenure diagnostics table moments
# -----------------------------------------------------------------------------
message("Computing tenure diagnostics...")

tenure_dt <- df_three[, .(y1, y2, y3, tenure1, tenure2, tenure3, weight)]
tenure_dt[, `:=`(
  tenure1_mo = fifelse(y1 == 1L & tenure1 >= 0, tenure1, NA_real_),
  tenure2_mo = fifelse(y2 == 1L & tenure2 >= 0, tenure2, NA_real_),
  tenure3_mo = fifelse(y3 == 1L & tenure3 >= 0, tenure3, NA_real_)
)]

ee_pairs <- rbindlist(list(
  tenure_dt[
    y1 == 1L & y2 == 1L &
      !is.na(tenure1_mo) & !is.na(tenure2_mo),
    .(dt_mo = tenure2_mo - tenure1_mo, weight)
  ],
  tenure_dt[
    y2 == 1L & y3 == 1L &
      !is.na(tenure2_mo) & !is.na(tenure3_mo),
    .(dt_mo = tenure3_mo - tenure2_mo, weight)
  ]
))

ee_exact <- weighted_percent(abs(ee_pairs$dt_mo - 3) < 0.01, ee_pairs$weight)
ee_non_clock <- 100 - ee_exact
ee_large_deviation <- weighted_percent(
  abs(ee_pairs$dt_mo - 3) > 12,
  ee_pairs$weight
)

s101 <- tenure_dt[
  y1 == 1L & y2 == 0L & y3 == 1L &
    !is.na(tenure1_mo) & !is.na(tenure3_mo),
  .(gap = tenure3_mo - tenure1_mo, tenure3_mo, weight)
]
s101[, classify := fcase(
  abs(gap - 6) <= 2,
  "continuous_employment_clock",
  tenure3_mo <= 3 & gap < 0,
  "fresh_start_clock",
  default = "other"
)]

s101_continuous <- weighted_percent(
  s101$classify == "continuous_employment_clock",
  s101$weight
)
s101_fresh_start <- weighted_percent(
  s101$classify == "fresh_start_clock",
  s101$weight
)

fresh_starts <- rbindlist(list(
  tenure_dt[
    y1 == 0L & y2 == 1L & !is.na(tenure2_mo),
    .(tenure_mo = tenure2_mo, weight)
  ],
  tenure_dt[
    y2 == 0L & y3 == 1L & !is.na(tenure3_mo),
    .(tenure_mo = tenure3_mo, weight)
  ]
))

fresh_start_short <- weighted_percent(
  fresh_starts$tenure_mo <= 3,
  fresh_starts$weight
)
fresh_start_long <- weighted_percent(
  fresh_starts$tenure_mo > 12,
  fresh_starts$weight
)

# -----------------------------------------------------------------------------
# 5. Write descriptive table
# -----------------------------------------------------------------------------
message("Writing LaTeX table to: ", descriptive_out_path)

descriptive_lines <- c(
  "% ===========================================================================",
  "% table_descriptive.tex",
  "% Generated by paper/build_descriptive_table.R",
  sprintf("%% Date: %s", Sys.Date()),
  "% ===========================================================================",
  paste0(latex_env("table"), "[htbp]"),
  latex_cmd("centering"),
  paste0(
    latex_cmd("caption"),
    "{Descriptive statistics: QLFS rotating panel, 2010Q1--2019Q4}"
  ),
  paste0(latex_cmd("label"), "{tab:descriptive}"),
  paste0(latex_env("tabular"), "{lc}"),
  latex_cmd("toprule"),
  latex_row(latex_textbf("Statistic"), " & ", latex_textbf("Value")),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(
    2,
    "l",
    latex_textit("Panel A: Sample (3-wave linked triads)")
  )),
  latex_row("Observations & ", format_int(sample_n)),
  latex_row("Age range & ", as.integer(age_min), "--", as.integer(age_max)),
  latex_row("Mean age (wave 1) & ", format_1(age_mean)),
  latex_row("Mean education (years, wave 1) & ", format_1(education_mean))
)

if (!is.na(female_share)) {
  descriptive_lines <- c(
    descriptive_lines,
    latex_row(latex_cmd("%"), " female & ", format_1(female_share), latex_cmd("%"))
  )
}

descriptive_lines <- c(
  descriptive_lines,
  latex_row("Employment rate: wave 1 (", latex_cmd("%"), ") & ", format_1(employment_wave1)),
  latex_row("Employment rate: wave 2 (", latex_cmd("%"), ") & ", format_1(employment_wave2)),
  latex_row("Employment rate: wave 3 (", latex_cmd("%"), ") & ", format_1(employment_wave3)),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(2, "l", latex_textit("Panel B: Observed transition rates"))),
  latex_row(latex_multicolumn(2, "l", paste0(latex_cmd("quad"), latex_textit("3-month rates")))),
  latex_row(latex_cmd("quad"), " Entry rate: wave 1$", latex_cmd("to"), "$2 (", latex_cmd("%"), ") & ", format_2(entry_12)),
  latex_row(latex_cmd("quad"), " Entry rate: wave 2$", latex_cmd("to"), "$3 (", latex_cmd("%"), ") & ", format_2(entry_23)),
  latex_row(latex_cmd("quad"), " Exit rate: wave 1$", latex_cmd("to"), "$2 (", latex_cmd("%"), ") & ", format_2(exit_12)),
  latex_row(latex_cmd("quad"), " Exit rate: wave 2$", latex_cmd("to"), "$3 (", latex_cmd("%"), ") & ", format_2(exit_23)),
  latex_row(latex_multicolumn(2, "l", paste0(latex_cmd("quad"), latex_textit("6-month rates")))),
  latex_row(latex_cmd("quad"), " Observed entry rate: wave 1$", latex_cmd("to"), "$3 (", latex_cmd("%"), ") & ", format_2(entry_6m_observed)),
  latex_row(latex_cmd("quad"), " Observed exit rate: wave 1$", latex_cmd("to"), "$3 (", latex_cmd("%"), ") & ", format_2(exit_6m_observed)),
  latex_row(latex_cmd("quad"), " Implied entry rate (2-step Markov, ", latex_cmd("%"), ") & ", format_2(entry_6m_compound)),
  latex_row(latex_cmd("quad"), " Implied exit rate (2-step Markov, ", latex_cmd("%"), ") & ", format_2(exit_6m_compound)),
  latex_row(latex_cmd("quad"), " Ratio: implied / observed entry & ", format_2(entry_6m_compound / entry_6m_observed)),
  latex_row(latex_cmd("quad"), " Ratio: implied / observed exit & ", format_2(exit_6m_compound / exit_6m_observed)),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(2, "l", latex_textit("Panel C: Data quality -- inconsistent responses"))),
  latex_row(latex_cmd("%"), " with age inconsistency (any wave) & ", format_1(age_incons)),
  latex_row(latex_cmd("%"), " with education inconsistency (any wave) & ", format_1(education_incons)),
  latex_row(latex_cmd("%"), " with either inconsistency & ", format_1(any_incons)),
  latex_cmd("bottomrule"),
  latex_row(latex_multicolumn(
    2,
    paste0("p{0.85", latex_cmd("textwidth"), "}"),
    paste0(
      latex_cmd("footnotesize"), " ", latex_textit("Note:"),
      " All rates are weighted using the geometric mean of wave-specific ",
      "survey weights. 3-month transition rates are the mean of wave ",
      "1$", latex_cmd("to"), "$2 and wave 2$", latex_cmd("to"), "$3 ",
      "conditional rates. 6-month rates use 4-wave linked observations ",
      "(N = ", format_int(n_four), "). Implied 6-month rates compound ",
      "the 3-month rates under the Markov assumption: $1 - (1-r)^2$. ",
      "Age inconsistency: consecutive wave gap not 0 or 1 years. ",
      "Education inconsistency: non-monotone or multi-category jump."
    )
  )),
  latex_env("tabular", "end"),
  latex_env("table", "end")
)

writeLines(descriptive_lines, descriptive_out_path)
message("Done. Table written to: ", descriptive_out_path)

# -----------------------------------------------------------------------------
# 6. Write tenure diagnostics table
# -----------------------------------------------------------------------------
message("Writing LaTeX table to: ", tenure_out_path)

tenure_lines <- c(
  "% ===========================================================================",
  "% table_tenure_diagnostics.tex",
  "% Generated by paper/build_descriptive_table.R",
  sprintf("%% Date: %s", Sys.Date()),
  "% ===========================================================================",
  paste0(latex_env("table"), "[htbp]"),
  latex_cmd("centering"),
  latex_env("threeparttable"),
  paste0(
    latex_cmd("caption"),
    "{Tenure-based evidence of spurious employment transitions}"
  ),
  paste0(latex_cmd("label"), "{tab:tenure_diagnostics}"),
  paste0(latex_env("tabular"), paste0("{p{0.58", latex_cmd("textwidth"), "}rr}")),
  latex_cmd("toprule"),
  latex_row(
    latex_textbf("Diagnostic moment"), " & ",
    latex_textbf("N"), " & ",
    latex_textbf(paste0("Weighted share (", latex_cmd("%"), ")"))
  ),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(
    3,
    "l",
    latex_textit("Panel A: Employed--employed continuation pairs")
  )),
  latex_row(
    "Exact tenure clock: $g_t-g_{t-1}=3$ months & ",
    format_int(nrow(ee_pairs)), " & ", format_1(ee_exact)
  ),
  latex_row(
    "Non-clock increment: $g_t-g_{t-1}", latex_cmd("neq"),
    " 3$ months & ", format_int(nrow(ee_pairs)), " & ",
    format_1(ee_non_clock)
  ),
  latex_row(
    "Large non-clock deviation: $|g_t-g_{t-1}-3|>12$ months & ",
    format_int(nrow(ee_pairs)), " & ", format_1(ee_large_deviation)
  ),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(
    3,
    "l",
    latex_textit("Panel B: Employed--non-employed--employed sequences")
  )),
  latex_row(
    "Tenure flank follows continuous-employment clock: $g_3",
    latex_cmd("approx"), " g_1+6$ months & ", format_int(nrow(s101)),
    " & ", format_1(s101_continuous)
  ),
  latex_row(
    "Wave-3 tenure consistent with a genuinely new job: $g_3",
    latex_cmd("leq"), " 3$ months & ", format_int(nrow(s101)),
    " & ", format_1(s101_fresh_start)
  ),
  latex_cmd("midrule"),
  latex_row(latex_multicolumn(
    3,
    "l",
    latex_textit("Panel C: Apparent within-panel job starts")
  )),
  latex_row(
    "Short reported tenure at start wave: $g_t", latex_cmd("leq"),
    " 3$ months & ", format_int(nrow(fresh_starts)), " & ",
    format_1(fresh_start_short)
  ),
  latex_row(
    "Impossible long tenure at start wave: $g_t>12$ months & ",
    format_int(nrow(fresh_starts)), " & ", format_1(fresh_start_long)
  ),
  latex_cmd("bottomrule"),
  latex_env("tabular", "end"),
  paste0(latex_env("tablenotes"), "[flushleft]"),
  paste0(
    latex_cmd("item"), " ", latex_textit("Note:"),
    " Rows report weighted shares among observations with the status pattern ",
    "shown in each panel. ", latex_textit("N"),
    " is the unweighted denominator: adjacent pairs in Panels A and C, ",
    "and three-wave triads in Panel B. Tenure $g$ is measured in months. ",
    "For $(1,0,1)$ sequences, the continuous-employment clock is defined ",
    "as $|(g_3-g_1)-6|", latex_cmd("leq"), " 2$ months. Apparent starts ",
    "pool observed $0", latex_cmd("to"), "1$ transitions over waves ",
    "1$", latex_cmd("to"), "$2 and 2$", latex_cmd("to"), "$3. ",
    "All shares use the geometric mean of wave-specific survey weights."
  ),
  latex_env("tablenotes", "end"),
  latex_env("threeparttable", "end"),
  latex_env("table", "end")
)

writeLines(tenure_lines, tenure_out_path)
message("Done. Table written to: ", tenure_out_path)