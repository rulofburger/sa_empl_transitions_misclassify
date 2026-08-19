# Tenure-contamination-specific recoding applied after shared ingestion.
prepare_eps_estimation_data <- function(df) {
  if (requireNamespace("haven", quietly=TRUE)) {
    labelled <- vapply(df, haven::is.labelled, logical(1))
    df[labelled] <- lapply(df[labelled], as.numeric)
  }
  nw_category <- function(age_years) as.integer(cut(
    pmax((age_years-16)*12,0), breaks=c(0,3,6,9,12,36,60,Inf),
    right=FALSE, include.lowest=TRUE))
  for (wave in 1:3) {
    nw <- paste0("neverworked",wave); tc <- paste0("timegap_cat",wave)
    age <- paste0("age",wave)
    if (all(c(nw,tc,age) %in% names(df))) {
      use <- !is.na(df[[nw]]) & as.numeric(df[[nw]]) == 1
      df[[tc]][use] <- nw_category(df[[age]][use])
    }
  }
  if (!"weight" %in% names(df)) df$weight <- df$weight1
  validate_df_eps(df)
  df
}
