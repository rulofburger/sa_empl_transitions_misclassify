#!/usr/bin/env Rscript

# Copy compact, paper-facing results out of ignored working-output directories.
# Run this from anywhere after estimating the models and rebuilding the tables:
#   Rscript --vanilla scripts/snapshot_paper_artifacts.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
  project_root <- normalizePath(file.path(dirname(script_path), ".."),
                                winslash = "/")
} else {
  project_root <- normalizePath(getwd(), winslash = "/")
}

snapshot_root <- file.path(project_root, "paper", "generated")

table_files <- c(
  "EM-baseline/output/tables/table_baseline_implied.tex",
  "EM-baseline/output/tables/table_baseline_params.tex",
  "EM-baseline/output/tables/table_matching_implied.tex",
  "EM-baseline-ext/output/tables/table_inconsistency.tex",
  "EM-baseline-ext/output/tables/table_inconsistency_appendix.tex",
  "EM-baseline-ext/output/tables/table_inconsistency_reliability.tex",
  "EM-baseline-ext/output/tables/table_inconsistency_extent.tex"
)

result_files <- c(
  "EM-baseline/output/results/run_summary.csv",
  "EM-baseline/output/results/matching_sensitivity_summary.csv",
  "EM-baseline-ext/output/results/table6_inconsistency_audit.csv",
  "EM-baseline-ext/output/results/table6_inconsistency_extent.csv",
  "EM-baseline-ext/output/results/table6_inconsistency_extent_audit.csv",
  "EM-baseline-ext/output/results/table6_inconsistency_prevalence.csv",
  "EM-baseline-ext/output/results/table6_inconsistency_severity_prevalence.csv",
  "EM-baseline-ext/output/results/table6_matching_membership_prevalence.csv",
  "EM-baseline-ext/output/results/table6_matching_quality_audit.csv",
  "EM-baseline-ext/output/results/table6_reliability_group_robustness.csv"
)

sources <- file.path(project_root, c(table_files, result_files))
missing_sources <- sources[!file.exists(sources)]
if (length(missing_sources)) {
  stop("Missing paper artifacts:\n", paste(" -", missing_sources,
                                            collapse = "\n"),
       "\nRun the Table 1--3 estimation and table pipeline first.")
}

relative_sources <- substring(
  normalizePath(sources, winslash = "/"),
  nchar(project_root) + 2L
)
destinations <- file.path(snapshot_root, relative_sources)

# Remove stale files from the generated snapshot, while preserving its README.
if (dir.exists(snapshot_root)) {
  existing <- list.files(snapshot_root, recursive = TRUE, full.names = TRUE,
                         all.files = TRUE, no.. = TRUE)
  existing <- existing[file.info(existing)$isdir %in% FALSE]
  keep <- normalizePath(file.path(snapshot_root, "README.md"), winslash = "/",
                        mustWork = FALSE)
  stale <- existing[normalizePath(existing, winslash = "/",
                                  mustWork = FALSE) != keep]
  if (length(stale)) unlink(stale)
}

for (i in seq_along(sources)) {
  dir.create(dirname(destinations[i]), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(sources[i], destinations[i], overwrite = TRUE,
                      copy.mode = FALSE, copy.date = FALSE)
  if (!copied) stop("Could not copy ", relative_sources[i])
  if (grepl("[.]tex$", destinations[i], ignore.case = TRUE)) {
    lines <- readLines(destinations[i], warn = FALSE)
    writeLines(sub("[ \\t]+$", "", lines), destinations[i], useBytes = TRUE)
  }
}

manifest <- data.frame(
  source = relative_sources,
  snapshot = file.path("paper", "generated", relative_sources),
  bytes = unname(file.info(destinations)$size),
  md5 = unname(tools::md5sum(destinations)),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(snapshot_root, "manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE, na = "")

cat(sprintf("Updated %d paper artifacts in %s\n",
            nrow(manifest), snapshot_root))
