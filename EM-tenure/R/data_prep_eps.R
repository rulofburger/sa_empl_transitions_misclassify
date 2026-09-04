# Tenure-contamination-specific recoding applied after shared ingestion.
prepare_eps_estimation_data <- function(df, allow_zero_tenure = isTRUE(
    getOption("sa_empl_transitions.preserve_zero_tenure", FALSE))) {
  if (requireNamespace("haven", quietly=TRUE)) {
    labelled <- vapply(df, haven::is.labelled, logical(1))
    df[labelled] <- lapply(df[labelled], as.numeric)
  }
  nw_category <- function(age_years) as.integer(cut(
    pmax((age_years-16)*12,0), breaks=c(0,3,6,9,12,36,60,Inf),
    right=FALSE, include.lowest=TRUE))
  for (wave in 1:3) {
    nw <- paste0("neverworked",wave); tc <- paste0("timegap_cat",wave)
    y <- paste0("y",wave); tenure <- paste0("tenure",wave)
    age <- paste0("age",wave)
    if (all(c(nw,tc,age) %in% names(df))) {
      use <- !is.na(df[[nw]]) & as.numeric(df[[nw]]) == 1
      df[[tc]][use] <- nw_category(df[[age]][use])
    }
    # Durations belonging to the unreported state are unavailable, not zero
    # and not observations that can be filled from another wave.  Preserve
    # that distinction for latent histories in which y_t is misclassified.
    if (all(c(y,tenure) %in% names(df))) df[[tenure]][df[[y]] == 0] <- NA_real_
    if (all(c(y,tc) %in% names(df))) df[[tc]][df[[y]] == 1] <- NA_integer_
  }
  if (!"weight" %in% names(df)) df$weight <- df$weight1
  validate_df_eps(df, allow_zero_tenure = allow_zero_tenure)
  df
}

# Attach the nominal QLFS interview month implied by the period index.  Exact
# interview dates are unavailable; the survey quarters are represented by
# March, June, September, and December, matching the start-date audit.
add_nominal_interview_months <- function(df) {
  period_cols <- paste0("period",1:3)
  if (length(setdiff(period_cols,names(df))))
    stop("Nominal interview months require period1-period3")
  for (wave in 1:3) {
    period <- as.integer(df[[period_cols[wave]]])
    if (any(!is.finite(period) | period < 1L))
      stop("period columns must contain positive integers")
    df[[paste0("interview_month",wave)]] <-
      3L*((period-1L)%%4L+1L)
  }
  df
}
