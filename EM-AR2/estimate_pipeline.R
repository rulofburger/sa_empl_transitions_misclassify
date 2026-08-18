# ==============================================================================
# EM-AR2: full-sample estimation with analytical standard errors
# ==============================================================================
# Uses every complete four-wave record in df_qlfs_A.rds. Identical observed
# histories are collapsed for EM estimation, without changing the likelihood.
# Analytical SEs use the observed-data Hessian, individual weighted scores,
# and the delta method. No bootstrap is run by this script.

if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source EM-AR2/estimate_pipeline.R from the project root")

source("EM-AR2/R/source_all.R")
source("scripts/ingest_data_4waves_SA.R")

required <- c("y1", "y2", "y3", "y4", "weight")
missing <- setdiff(required, names(df_qlfs))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

# All complete four-wave observations are retained; ingestion applies only the
# paper's age restriction and complete-case requirement for employment status.
df_4w <- df_qlfs[, required]
n_obs <- nrow(df_4w)
weight_total <- sum(df_4w$weight)
df_4w$weight <- n_obs * df_4w$weight / weight_total
message(sprintf("EM-AR2 full sample: N = %s; periods %d-%d.",
                format(n_obs, big.mark = ","), min(df_qlfs$period1), max(df_qlfs$period1)))

settings <- list(max_iter = 2000L, tol = 1e-10, verbose = 1L,
                 collapse_cells = TRUE)
fit_sym <- do.call(em_fit_ar2, c(list(df = df_4w, estimate_pi = TRUE), settings))
fit_nopi <- do.call(em_fit_ar2, c(list(df = df_4w, estimate_pi = FALSE,
                                       fixed_pi = 0), settings))
fit_asym <- do.call(em_fit_ar2, c(list(df = df_4w, asymmetric = TRUE), settings))

fits <- list(nopi = fit_nopi, sym = fit_sym, asym = fit_asym)
model_types <- c(nopi = "none", sym = "symmetric", asym = "asymmetric")
inference <- Map(function(fit, type) analytical_se_ar2(df_4w, fit, type),
                 fits, model_types)

results_dir <- "EM-AR2/output/results"
tables_dir <- "EM-AR2/output/tables"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

for (nm in names(fits)) {
  fits[[nm]]$analytical_inference <- inference[[nm]]
  saveRDS(fits[[nm]], file.path(results_dir,
    sprintf("em_ar2_%s_%s.rds", nm, run_id)))
}

summary_rows <- do.call(rbind, lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  inf <- inference[[nm]]
  p <- fit$params
  data.frame(
    run_id = run_id, model = nm, converged = fit$converged,
    iterations = fit$iterations, loglik = fit$loglik,
    theta0 = p$theta0, theta01 = p$theta01,
    theta1 = p$theta1, theta10 = p$theta10,
    pi = p$pi %||% NA_real_, pi0 = p$pi0 %||% NA_real_,
    pi1 = p$pi1 %||% NA_real_, n_obs = n_obs,
    max_abs_score = inf$diagnostics$max_abs_average_score,
    information_condition = inf$diagnostics$information_condition,
    stringsAsFactors = FALSE)
}))
write.csv(summary_rows, file.path(results_dir,
  sprintf("run_summary_%s.csv", run_id)), row.names = FALSE)

analytical_rows <- do.call(rbind, lapply(names(inference), function(nm) {
  cbind(model = nm, inference[[nm]]$summary, row.names = NULL)
}))
write.csv(analytical_rows, file.path(results_dir,
  sprintf("analytical_se_%s.csv", run_id)), row.names = FALSE)

# Stable latest files are intentionally overwritten by each complete pipeline
# run; timestamped fit objects retain an auditable history.
write.csv(summary_rows, file.path(results_dir, "run_summary_latest.csv"), row.names = FALSE)
write.csv(analytical_rows, file.path(results_dir, "analytical_se_latest.csv"), row.names = FALSE)

cat("\nAR(2) full-sample estimates with analytical standard errors\n")
print(analytical_rows[analytical_rows$quantity %in%
  c("p_00", "p_10", "p_01", "p_11", "employment_rate", "pi", "pi0", "pi1"), ],
  row.names = FALSE, digits = 5)
cat("\nDiagnostics\n")
print(summary_rows[, c("model", "converged", "iterations", "loglik", "n_obs",
                       "max_abs_score", "information_condition")],
      row.names = FALSE, digits = 5)
