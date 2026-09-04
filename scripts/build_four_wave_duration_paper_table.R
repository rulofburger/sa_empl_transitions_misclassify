# Rebuild a paper-facing table from saved fits only; no estimation or inference.
# Run from the repository root.
source("EM-tenure/R/source_all.R")
paths <- c(
  ar2="EM-tenure/output/results/job_change_monthly/four_wave_ar2/empirical/comparison_latest.rds",
  ar1="EM-tenure/output/results/job_change_monthly/four_wave_ar1/continuous_clock_empirical/comparison_latest.rds",
  data="EM-tenure/output/results/job_change_monthly/four_wave_ar1/prepared_cells_latest.rds")
saved <- readRDS(paths["ar2"])
stopifnot(identical(saved$source_md5, four_wave_ar2_source_fingerprint()),
          !isTRUE(saved$standard_errors_estimated))
fits <- list(AR1=readRDS(paths["ar1"])$best, AR2=saved$best)
stopifnot(all(vapply(fits, function(f) isTRUE(f$converged), logical(1))))
df <- readRDS(paths["data"])$df4
data <- prepare_four_wave_kernel_data(df)
h <- latent_histories_eps_4w()

# Condition on latent origin and preceding transition, not reported switching.
# All four reports inform posterior risk-set membership. Pool intervals 2:3.
history_rates <- function(fit) {
  risks <- .ar2_baseline_risks(data, fit$params)
  effects <- .ar2_effects(fit$params)
  rows <- list()
  for (origin in 0:1) for (recent in c(FALSE, TRUE)) {
    name <- if (origin == 0) "entry" else "exit"
    numerator <- denominator <- 0
    for (t in 2:3) {
      selected <- h[,t] == origin & (h[,t] != h[,t-1]) == recent
      w <- df$weight * rowSums(fit$gamma[,selected,drop=FALSE])
      p <- risks[[name]][,t]
      if (recent) p <- .ar2_shift(p, effects[name])
      numerator <- numerator + sum(w*p)
      denominator <- denominator + sum(w)
    }
    rows[[length(rows)+1L]] <- data.frame(rate=name, recent=recent,
      numerator=numerator, denominator=denominator,
      probability=numerator/denominator)
  }
  out <- do.call(rbind, rows)
  stopifnot(all(is.finite(out$probability)),
            all(out$probability > 0 & out$probability < 1))
  # Grouped risk sets must reproduce the existing interval-rate calculation.
  overall <- duration_weighted_transition_rates_ar2(df, fit)
  for (name in c("entry", "exit")) {
    group <- out[out$rate == name,]
    interval <- overall[overall$wave %in% 2:3,]
    stopifnot(abs(sum(group$numerator) -
      sum(interval[[paste0(name,"_numerator")]])) < 1e-6,
      abs(sum(group$denominator) -
      sum(interval[[paste0(name,"_denominator")]])) < 1e-6)
  }
  out
}
groups <- lapply(fits, history_rates)
labels <- c("Entry rate", "Exit rate", "Misclassification rate",
  "Initial employment rate", "Entry after a recent exit",
  "Entry after continued nonemployment", "Exit after a recent entry",
  "Exit after continued employment", "Log likelihood", "N", "Parameters")
values <- lapply(names(fits), function(model) {
  fit <- fits[[model]]
  rates <- duration_weighted_transition_rates_ar2(df, fit)[1,]
  group <- groups[[model]]
  get_rate <- function(name, recent) group$probability[
    group$rate == name & group$recent == recent]
  c(rates$entry_rate, rates$exit_rate, fit$params$pi, fit$params$alpha,
    get_rate("entry", TRUE), get_rate("entry", FALSE),
    get_rate("exit", TRUE), get_rate("exit", FALSE),
    fit$loglik, attr(df,"n_original"),
    length(fit$par_unconstrained))
})
names(values) <- names(fits)
numeric_table <- data.frame(measure=labels, AR1=values$AR1, AR2=values$AR2)
stopifnot(numeric_table$AR1[10] == numeric_table$AR2[10],
          numeric_table$AR2[11] - numeric_table$AR1[11] == 2)
# Check against the consolidated main rates, rather than copying rounded text.
for (model in names(fits)) {
  key <- if (model == "AR1") "AR1" else saved$best_label
  reference <- saved$summary[saved$summary$model == key,]
  stopifnot(max(abs(values[[model]][1:4] - unlist(reference[
    c("entry","exit","misclassification","initial_employment")]))) < 1e-10)
}
display <- numeric_table
for (model in names(fits)) {
  display[[model]] <- c(sprintf("%.3f", 100*values[[model]][1:8]),
    formatC(values[[model]][9],format="f",digits=2,big.mark=","),
    formatC(values[[model]][10:11],format="f",digits=0,big.mark=","))
}
notes <- paste(
  "Notes: All rates are percentages. Panel A; four waves; identical samples.",
  "Entry and exit rates are survey- and posterior-risk-weighted fitted probabilities",
  "pooled across all three quarterly transitions. History-specific rates pool only",
  "transitions 2--3 and 3--4; recent means a latent status change in the preceding",
  "interval, while continued means no such change. All four reports inform the",
  "posterior risk sets; these are not raw observed-switcher rates or real-time forecasts.",
  "Groups have different duration compositions. Both models retain the same",
  "piecewise duration hazards, duration-reporting mechanisms, and AR(1) first-transition",
  "initialization. The proposed new-nonemployment onset restriction and alternative",
  "same-employer return treatment have not been implemented. Likelihood weights are",
  "normalized to N. Standard errors and calibrated significance tests are pending;",
  "no significance stars are reported.")
out <- "paper/generated/four_wave_duration"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(numeric_table,file.path(out,"table_four_wave_duration_values.csv"),row.names=FALSE)
group_table <- do.call(rbind,lapply(names(groups),function(model)
  data.frame(model=model,groups[[model]])))
write.csv(group_table,file.path(out,"history_specific_rates.csv"),row.names=FALSE)
body <- vapply(seq_along(labels),function(i) paste0(
  if (i %in% c(5,9)) "\\addlinespace\n" else "",
  paste(display[i,],collapse=" & ")," \\\\"),character(1))
tex <- c("% Generated by scripts/build_four_wave_duration_paper_table.R; do not edit.",
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{Four-wave duration model: transition rates and recent reversals}",
  "\\label{tab:four-wave-duration-ar2}", "\\begin{tabular}{lrr}",
  "\\toprule", " & AR(1) & AR(2) \\\\", "\\midrule", body,
  "\\bottomrule", "\\end{tabular}", "\\par\\vspace{0.5em}",
  "\\begin{minipage}{0.98\\linewidth}\\footnotesize", notes,
  "\\end{minipage}", "\\end{table}")
writeLines(tex,file.path(out,"table_four_wave_duration.tex"))
markdown <- c("| Measure | AR(1) | AR(2) |", "|:--|--:|--:|",
  vapply(seq_along(labels),function(i)
    paste0("| ",paste(display[i,],collapse=" | ")," |"),character(1)),
  "", notes)
writeLines(markdown,file.path(out,"table_four_wave_duration.md"))
provenance_paths <- c(paths, builder="scripts/build_four_wave_duration_paper_table.R")
write.csv(data.frame(input=names(provenance_paths),path=unname(provenance_paths),
  md5=unname(tools::md5sum(provenance_paths))),file.path(out,"provenance.csv"),row.names=FALSE)
cat(markdown,sep="\n")
cat("\nAll table and grouped-rate checks passed. No models re-estimated.\n")
