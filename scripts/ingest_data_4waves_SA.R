
# Four-wave QLFS panel used by the AR(2) model. This script is deliberately
# dependency-free so it behaves identically from R, Rscript, and tests.
df_qlfs <- readRDS("data/raw/df_qlfs_A.rds")

keep <- df_qlfs$age1 > 17 & df_qlfs$age1 < 56 &
  stats::complete.cases(df_qlfs[, paste0("employed", 1:4), drop = FALSE])
df_qlfs <- df_qlfs[keep, , drop = FALSE]

name_pattern <- paste(c("employed", "age", "educ", "race", "female", "weight",
                        "tenure", "timegap", "formal", "period", "contracttype",
                        "neverworked"), collapse = "|")
df_qlfs <- df_qlfs[, grepl(name_pattern, names(df_qlfs)), drop = FALSE]
names(df_qlfs)[match(paste0("employed", 1:4), names(df_qlfs))] <- paste0("y", 1:4)
df_qlfs$weight <- with(df_qlfs,
  (weight1 * weight2 * weight3 * weight4) ^ 0.25)
df_qlfs <- as.data.frame(df_qlfs)

