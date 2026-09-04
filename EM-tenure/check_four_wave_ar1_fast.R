# Lightweight regression checks; no model estimation or bootstrap.
testthat::test_dir("EM-tenure/tests/testthat",
  filter = "four-wave|piecewise-risk-exact|tenure-monthly-eps",
  reporter = "summary", stop_on_failure = TRUE)

# Check the audit's result table against the saved output without rendering PDF.
audit <- readLines("paper/tenure_start_date_model_audit.Rmd", warn = FALSE)
audit_env <- new.env(parent = globalenv())
for (label in c("four-wave-bridge-results", "four-wave-recovery-score",
    "four-wave-recovery-mass", "four-wave-clock-repair-checks",
    "four-wave-clock-recovery-fits", "four-wave-clock-recovery-hazards",
    "four-wave-clock-empirical", "four-wave-clock-empirical-hazards",
    "four-wave-ar2-recovery", "four-wave-ar2-empirical", "four-wave-ar2-implications",
    "four-wave-ar2-hazards")) {
  start <- grep(paste0("^```\\{r ",label,"[,}]"), audit)
  stopifnot(length(start) == 1L)
  end <- start + which(audit[seq.int(start + 1L, length(audit))] == "```")[1L]
  eval(parse(text = audit[seq.int(start + 1L, end - 1L)]), envir = audit_env)
}
