# Collect completed four-wave runs without evaluating any three-wave model.
source("EM-tenure/R/source_all.R")
outdir <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
df <- readRDS(file.path(outdir, "prepared_cells_latest.rds"))$df4
paths <- list.files(outdir, "^fit_fast_.*_latest[.]rds$", full.names = TRUE)
if (!length(paths)) stop("No completed fast refinement runs to summarize")
fits <- lapply(paths, readRDS)
names(fits) <- sub("_latest[.]rds$", "", sub("^fit_fast_", "", basename(paths)))
# This is the historical pre-repair comparison. Do not mix the repaired
# likelihood into its starting-point comparison or overwrite its interpretation.
legacy <- vapply(fits,function(fit)
  is.null(fit$params$timegap_clock_model) ||
    identical(fit$params$timegap_clock_model,"legacy_pairwise"),logical(1))
fits <- fits[legacy]; paths <- paths[legacy]
if (!length(fits)) stop("No legacy fits; use a separate repaired-model comparison")
preliminary <- readRDS(file.path(outdir, "fit_reporting2_latest.rds"))
all_fits <- c(list(preliminary = preliminary), fits)
summary <- do.call(rbind, lapply(names(all_fits), function(label) {
  fit <- all_fits[[label]]
  rates <- duration_weighted_transition_rates_4w(df, fit)[1L, ]
  data.frame(stage = label, observations = attr(df, "n_original"),
    loglik = fit$loglik, gain_over_preliminary = fit$loglik - preliminary$loglik,
    entry_rate = rates$entry_rate, exit_rate = rates$exit_rate,
    status_misclassification = fit$params$pi,
    initial_employment = fit$params$alpha,
    projected_score = if (is.null(fit$projected_score)) NA_real_ else
      max(abs(fit$projected_score)),
    optimizer_code = fit$convergence,
    independently_converged = isTRUE(fit$converged))
}))
eligible <- which(vapply(fits, function(fit) isTRUE(fit$converged), logical(1)))
if (!length(eligible)) stop("No run passes the independent convergence checks")
best_index <- eligible[which.max(vapply(fits[eligible], function(fit)
  fit$loglik, numeric(1)))]
selected <- fits[[best_index]]
multistart <- data.frame(start = names(fits),
  loglik_gap = selected$loglik - vapply(fits, function(fit) fit$loglik, numeric(1)),
  max_parameter_gap = vapply(fits, function(fit)
    max(abs(fit$par_unconstrained - selected$par_unconstrained)), numeric(1)),
  independently_converged = vapply(fits, function(fit) isTRUE(fit$converged), logical(1)))
selected_row <- summary[summary$stage == names(fits)[best_index], ]
for (rate in c("entry_rate", "exit_rate", "status_misclassification"))
  multistart[[paste0(rate, "_gap")]] <-
    summary[[rate]][match(names(fits), summary$stage)] - selected_row[[rate]]
write.csv(summary, file.path(outdir, "convergence_summary_latest.csv"), row.names=FALSE)
write.csv(multistart, file.path(outdir, "multistart_summary_latest.csv"), row.names=FALSE)
saveRDS(list(best = selected, best_label = names(fits)[best_index],
  summary = summary, multistart = multistart, input_paths = paths,
  multistart_completed = length(fits) >= 3L &&
    all(vapply(fits, function(fit) isTRUE(fit$converged), logical(1))),
  same_person_three_wave_reestimated = FALSE),
  file.path(outdir, "converged_comparison_latest.rds"))
print(summary, row.names = FALSE, digits = 9)
print(multistart, row.names = FALSE, digits = 9)
