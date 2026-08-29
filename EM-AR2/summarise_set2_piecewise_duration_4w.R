# Reference duration-hazard profiles for the Table 5 piecewise-duration model.
if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source this script from the project root")
fit_path <- paste0("EM-AR2/output/results/set2_piecewise_duration/",
                   "fit_ar2_set2_piecewise_latest.rds")
fit <- readRDS(fit_path)
eta <- fit$eta
V <- fit$analytical_inference$vcov_eta

duration_labels <- c("0--3 months", "4--6 months", "7--12 months",
  "13--24 months", "25--60 months", "Over 60 months")
entry_terms <- c(NA, paste0("entry_timegap_",
  c("4_6m", "7_12m", "13_24m", "25_60m", "over_60m")))
exit_terms <- c(NA, paste0("persistence_tenure_",
  c("4_6m", "7_12m", "13_24m", "25_60m", "over_60m")))

profile_row <- function(equation, lag2, duration, duration_term) {
  x <- setNames(numeric(length(eta)), names(eta))
  if (equation == "Entry") {
    x["entry_intercept"] <- 1
    x["entry_lag2"] <- lag2
    if (!is.na(duration_term)) x[duration_term] <- 1
    index <- sum(x * eta)
    estimate <- pnorm(index)
    gradient <- dnorm(index) * x
  } else {
    x["persistence_intercept"] <- 1
    x["persistence_lag2"] <- lag2
    if (!is.na(duration_term)) x[duration_term] <- 1
    index <- sum(x * eta)
    estimate <- 1 - pnorm(index)
    gradient <- -dnorm(index) * x
  }
  data.frame(equation = equation, second_lag_employed = lag2,
    duration = duration, estimate = estimate,
    se = sqrt(as.numeric(t(gradient) %*% V %*% gradient)))
}

profile <- do.call(rbind, c(
  lapply(0:1, function(lag) do.call(rbind,
    Map(function(label, term) profile_row("Entry", lag, label, term),
        duration_labels, entry_terms))),
  lapply(0:1, function(lag) do.call(rbind,
    Map(function(label, term) profile_row("Exit", lag, label, term),
        duration_labels, exit_terms)))))
output_path <- paste0("EM-AR2/output/results/set2_piecewise_duration/",
                      "piecewise_duration_hazard_profile.csv")
write.csv(profile, output_path, row.names = FALSE)
print(transform(profile, estimate = 100 * estimate, se = 100 * se),
      row.names = FALSE, digits = 5)
