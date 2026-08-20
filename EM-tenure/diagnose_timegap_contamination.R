# Decompose timegap inconsistencies under the fitted contamination model.
# This script does not estimate or modify the model.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
fit_file <- "EM-tenure/output/results/timegap_contamination/fits_latest.rds"
if (!file.exists(fit_file)) stop("Run EM-tenure/estimate_timegap_contamination.R first")
saved <- readRDS(fit_file); fit <- saved$fit
if (nrow(fit$gamma)!=nrow(df_fit)) stop("Saved fit does not match current data")

decomp <- timegap_contamination_decomposition(df_fit,fit)
summary <- decomp$summary
summary$share_percent <- 100*summary$share_of_latent_continuations
summary$contamination_percent <- 100*summary$posterior_contamination_share
top_cells <- do.call(rbind,lapply(split(decomp$cells,decomp$cells$transition),
  function(x) head(x,12L)))
top_cells$share_percent <- 100*top_cells$share_within_transition
top_cells$contamination_percent <- 100*top_cells$posterior_contamination_share

cat("\nPosterior decomposition of reported N-to-N timegap transitions\n")
print(summary[,c("transition","mechanism","posterior_weight","share_percent",
  "contamination_percent")],row.names=FALSE,digits=7)
cat("\nLargest category-to-category cells within each interview interval\n")
print(top_cells[,c("transition","prev_category","curr_category","mechanism",
  "posterior_weight","share_percent","contamination_percent")],
  row.names=FALSE,digits=7)

outdir <- "EM-tenure/output/results/timegap_contamination"
write.csv(summary,file.path(outdir,"mechanism_decomposition_latest.csv"),
  row.names=FALSE)
write.csv(decomp$cells,file.path(outdir,"category_cell_decomposition_latest.csv"),
  row.names=FALSE)
