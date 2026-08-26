# ==============================================================================
# EM-baseline: exact eight-cell MLE and analytical inference pipeline
#
# Uses data/raw/df_qlfs_A.rds through the existing ingestion script. Alternative
# B/C matching panels are deliberately reserved for a separate sensitivity run.
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2) # required by diagnostic plots in the shared ingestion script

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

RANDOM_SEED <- 1234L
N_STARTS <- 12L
set.seed(RANDOM_SEED)

results_dir <- here::here("EM-baseline", "output", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

for (nm in c("y1", "y2", "y3")) {
  if (is.factor(df_qlfs[[nm]])) stop("Employment status is factor encoded: ", nm)
  values <- df_qlfs[[nm]][!is.na(df_qlfs[[nm]])]
  if (!all(values %in% c(0, 1))) stop("Non-binary employment status: ", nm)
}

df_baseline <- df_qlfs |>
  filter(!is.na(y1), !is.na(y2), !is.na(y3),
         !is.na(weight), is.finite(weight), weight > 0) |>
  transmute(y1 = as.integer(y1), y2 = as.integer(y2), y3 = as.integer(y3),
            weight = as.numeric(weight)) |>
  as.data.frame()
rm(df_qlfs)

cells <- collapse_baseline_cells(df_baseline)
cat(sprintf("Default panel: df_qlfs_A.rds | N=%s | sum(weights)=%.3f\n",
            format(cells$n, big.mark = ","), cells$weight_sum))

configs <- list(
  list(model_type = "none", stationary = TRUE, label = "none_stat"),
  list(model_type = "none", stationary = FALSE, label = "none_free"),
  list(model_type = "symmetric", stationary = TRUE, label = "sym_stat"),
  list(model_type = "symmetric", stationary = FALSE, label = "sym_free"),
  list(model_type = "asymmetric", stationary = TRUE, label = "asym_stat"),
  list(model_type = "asymmetric", stationary = FALSE, label = "asym_free")
)

.observed_start <- function(model_type, stationary) {
  w <- df_baseline$weight
  y1 <- df_baseline$y1; y2 <- df_baseline$y2; y3 <- df_baseline$y3
  denom0 <- sum(w * (y1 == 0)) + sum(w * (y2 == 0))
  denom1 <- sum(w * (y1 == 1)) + sum(w * (y2 == 1))
  p <- list(
    alpha = sum(w * y1) / sum(w),
    theta0 = (sum(w * (y1 == 0 & y2 == 1)) + sum(w * (y2 == 0 & y3 == 1))) / denom0,
    theta1 = (sum(w * (y1 == 1 & y2 == 1)) + sum(w * (y2 == 1 & y3 == 1))) / denom1
  )
  if (stationary) p$alpha <- stationary_alpha(p$theta0, p$theta1)
  if (model_type == "symmetric") p$pi <- 0.025
  if (model_type == "asymmetric") { p$pi0 <- 0.025; p$pi1 <- 0.025 }
  p
}

.nested_start <- function(cfg, fits) {
  suffix <- if (cfg$stationary) "stat" else "free"
  if (cfg$model_type == "symmetric") {
    p <- fits[[paste0("none_", suffix)]]$params
    p$pi <- 0.01
    return(p)
  }
  if (cfg$model_type == "asymmetric") {
    p0 <- fits[[paste0("sym_", suffix)]]$params
    return(list(alpha = p0$alpha, theta0 = p0$theta0, theta1 = p0$theta1,
                pi0 = p0$pi, pi1 = p0$pi))
  }
  if (!cfg$stationary && !is.null(fits$none_stat)) {
    p <- fits$none_stat$params
    p$alpha <- sum(df_baseline$weight * df_baseline$y1) / sum(df_baseline$weight)
    return(p)
  }
  NULL
}

.make_starts <- function(cfg, fits) {
  seeds <- list(init_params(cfg$model_type, cfg$stationary),
                .observed_start(cfg$model_type, cfg$stationary))
  nested <- .nested_start(cfg, fits)
  if (!is.null(nested)) seeds <- c(list(nested), seeds)
  base <- seeds[[1L]]
  eta <- pack_baseline_params(base, cfg$model_type, cfg$stationary)
  while (length(seeds) < N_STARTS) {
    seeds[[length(seeds) + 1L]] <- unpack_baseline_eta(
      eta + rnorm(length(eta), 0, 0.8), cfg$model_type, cfg$stationary
    )
  }
  seeds
}

fits <- list()
for (cfg in configs) {
  cat("\n--- ", cfg$label, " ---\n", sep = "")
  fit <- fit_baseline_mle(
    df_baseline, cfg$model_type, cfg$stationary,
    starts = .make_starts(cfg, fits), compute_gamma = FALSE, verbose = 1L,
    source_panel = "df_qlfs_A.rds"
  )
  fit$label <- cfg$label
  if (!fit$converged) stop(cfg$label, " failed convergence diagnostics")
  fits[[cfg$label]] <- fit
}

check_baseline_nesting(fits)
cat("\nAll likelihood nesting checks passed.\n")

analytical <- list()
for (cfg in configs) {
  analytical[[cfg$label]] <- analytical_se_baseline(df_baseline, fits[[cfg$label]])
}

# Save only after every fit and inference calculation has passed.
for (cfg in configs) {
  label <- cfg$label
  saveRDS(fits[[label]], file.path(results_dir, paste0("fit_", label, ".rds")))
  saveRDS(analytical[[label]],
          file.path(results_dir, paste0("analytical_se_", label, ".rds")))
}

run_rows <- lapply(configs, function(cfg) {
  fit <- fits[[cfg$label]]; p <- fit$params
  data.frame(
    timestamp = run_ts, source_panel = fit$sample$source_panel,
    sample_signature = fit$sample$signature, N = fit$sample$n,
    weight_sum = fit$sample$weight_sum, model_type = cfg$model_type,
    stationary = cfg$stationary, label = cfg$label,
    converged = fit$converged, optimizer_code = fit$diagnostics$optimizer_code,
    max_abs_score = fit$diagnostics$max_abs_score,
    information_min_eigenvalue = fit$diagnostics$information_min_eigenvalue,
    loglik = fit$loglik, alpha = p$alpha, theta0 = p$theta0, theta1 = p$theta1,
    pi = p$pi %||% NA_real_, pi0 = p$pi0 %||% NA_real_, pi1 = p$pi1 %||% NA_real_,
    stringsAsFactors = FALSE
  )
})
run_summary <- do.call(rbind, run_rows)
write.csv(run_summary, file.path(results_dir, "run_summary.csv"), row.names = FALSE)

cat("\n--- Validated baseline estimates ---\n")
display <- do.call(rbind, lapply(configs, function(cfg) {
  fit <- fits[[cfg$label]]; imp <- implied_baseline(fit$params, cfg$model_type)
  se <- analytical[[cfg$label]]$summary
  lookup <- setNames(se$se, se$quantity)
  data.frame(
    model = cfg$label,
    entry_pct = 100 * imp$entry_rate, entry_se = 100 * lookup[["entry_rate"]],
    exit_pct = 100 * imp$exit_rate, exit_se = 100 * lookup[["exit_rate"]],
    employment_pct = 100 * imp$employment_rate,
    employment_se = 100 * lookup[["employment_rate"]],
    pi_pct = ifelse(is.na(imp$pi), NA, 100 * imp$pi),
    loglik_m = fit$loglik / 1e6,
    max_score = fit$diagnostics$max_abs_score
  )
}))
print(display, row.names = FALSE, digits = 6)
cat("\nSaved validated fits, analytical SEs, and current-run summary.\n")
