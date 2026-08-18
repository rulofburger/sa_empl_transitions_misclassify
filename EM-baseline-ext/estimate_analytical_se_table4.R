library(here)
library(dplyr)
library(ggplot2)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))

source(here::here("scripts", "ingest_data_3waves_SA.R"))
sector_source <- readRDS(here::here("data", "raw", "QLFSmerged_mapped.rds"))
df_qlfs <- attach_transition_informal_sector(df_qlfs, sector_source)
rm(sector_source)
keep <- complete.cases(df_qlfs[, c("y1", "y2", "y3", "weight",
                                    "age1", "age2", "age3",
                                    "educ1", "educ2", "educ3")]) &
  df_qlfs$weight > 0
df_ext <- as.data.frame(df_qlfs[keep, , drop = FALSE])
df_ext$y1 <- as.integer(df_ext$y1); df_ext$y2 <- as.integer(df_ext$y2)
df_ext$y3 <- as.integer(df_ext$y3); df_ext$weight <- as.numeric(df_ext$weight)
for (nm in c("contracttype1", "contracttype2"))
  df_ext[[nm]] <- ifelse(is.na(df_ext[[nm]]), 0L, as.integer(df_ext[[nm]]))
rm(df_qlfs)

cv1 <- prepare_covariate_matrix(df_ext, 1L)
cv2 <- prepare_covariate_matrix(df_ext, 2L)
cv3 <- prepare_covariate_matrix(df_ext, 3L)
designs <- list(s1 = cv1$X, s2 = cv2$X, s3 = cv3$X_transition)
labels <- c("cov_s1_non_free", "cov_s2_non_free", "cov_s3_non_free",
            "cov_s1_sym_free", "cov_s2_sym_free", "cov_s3_sym_free")
results_dir <- here::here("EM-baseline-ext", "output", "results")

for (label in labels) {
  set_name <- sub("^cov_(s[123])_.*$", "\\1", label)
  model_type <- if (grepl("_sym_", label)) "symmetric" else "none"
  fit <- readRDS(file.path(results_dir, paste0("fit_", label, ".rds")))
  cat(sprintf("[%s] analytical sandwich/delta SE...\n", label))
  t0 <- proc.time()[["elapsed"]]
  out <- analytical_se_covariates(df_ext, designs[[set_name]], fit, model_type)
  out$label <- label; out$model_type <- model_type
  out$elapsed_seconds <- proc.time()[["elapsed"]] - t0
  saveRDS(out, file.path(results_dir, paste0("analytical_se_", label, ".rds")))
  cat(sprintf("  cells=%d K=%d condition=%.3g max|score|=%.3g elapsed=%.1fs\n",
              out$diagnostics$cells, out$diagnostics$K,
              out$diagnostics$condition_number,
              out$diagnostics$max_abs_score, out$elapsed_seconds))
  print(out$summary, row.names = FALSE, digits = 5)
}
