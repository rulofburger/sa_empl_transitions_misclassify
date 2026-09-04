# Estimate the preferred Table 7 model with exact, local, and gross tenure
# reports. Timegap contamination retains the joint per-wave marginal model.
# No bootstrap is run.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
baseline_file <- "EM-tenure/output/results/timegap_contamination_robustness/fits_latest.rds"
if (!file.exists(baseline_file)) stop("Run the Table 7 robustness estimator first")
robust_saved <- readRDS(baseline_file)
baseline <- robust_saved$joint$fit
if (nrow(baseline$gamma) != nrow(df_fit)) stop("Saved fit does not match data")

outdir <- "EM-tenure/output/results/tenure_local_gross"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
resume <- identical(tolower(Sys.getenv("TENURE_HYBRID_RESUME", "true")), "true")
prelim_maxit <- as.integer(Sys.getenv("TENURE_HYBRID_PRELIM_MAXIT", "300"))
refine_maxit <- as.integer(Sys.getenv("TENURE_HYBRID_REFINE_MAXIT", "500"))

# Conditional profile at the saved optimum verifies the direction and selects
# empirically useful starting shares without optimizing the nuisance block.
profile_shares <- c(.001, .01, .03, .06, .10, .20, .35, .50, .70, .90)
profile <- do.call(rbind, lapply(profile_shares, function(share) {
  p <- baseline$params
  p$eps_local <- p$eps * share
  p$eps_gross <- p$eps * (1 - share)
  p$tenure_contamination_model <- "local_gross"
  value <- e_step_eps(df_fit, p, check_df = FALSE, suff_stats = FALSE)$loglik
  data.frame(local_share = share, local_rate = p$eps_local,
    gross_rate = p$eps_gross, loglik = value)
}))
profile$loglik_difference <- profile$loglik - baseline$loglik
write.csv(profile, file.path(outdir, "conditional_profile_latest.csv"),
  row.names = FALSE)

# A fixed dispersed screen is used instead of an open-ended derivative-free
# optimizer. At the boundary the 14 nuisance parameters are already at their
# accepted optimum. Alternative nuisance blocks come from the independent and
# local-timegap fits, with their pair probabilities converted to per-wave
# probabilities before imposing the preferred joint timegap model.
candidate_blocks <- list(
  preferred_joint = baseline$params,
  independent_timegap = robust_saved$marginal$params,
  local_timegap = robust_saved$local$fit$params)
for (nm in setdiff(names(candidate_blocks), "preferred_joint")) {
  candidate_blocks[[nm]]$eps_d <- 1 - sqrt(1 - candidate_blocks[[nm]]$eps_d)
  candidate_blocks[[nm]]$timegap_contamination_model <- "joint_marginal"
}
screen_shares <- c(.001, .01, .05, .15, .35, .60)
message("Evaluating a fixed dispersed interior screen")
screen <- do.call(rbind, lapply(names(candidate_blocks), function(block) {
  do.call(rbind, lapply(screen_shares, function(share) {
    p <- candidate_blocks[[block]]
    p$eps_local <- p$eps * share
    p$eps_gross <- p$eps * (1 - share)
    p$tenure_contamination_model <- "local_gross"
    p$timegap_contamination_model <- "joint_marginal"
    value <- e_step_eps(df_fit, p, check_df = FALSE, suff_stats = FALSE)$loglik
    data.frame(block = block, local_share = share,
      local_rate = p$eps_local, gross_rate = p$eps_gross, loglik = value,
      loglik_difference = value - baseline$loglik)
  }))
}))
write.csv(screen, file.path(outdir, "dispersed_screen_latest.csv"),
  row.names = FALSE)
if (max(screen$loglik) > baseline$loglik + 1e-4)
  stop("An interior screen candidate improves on the boundary; full refinement is required")

message("No interior candidate improves on the gross-only boundary")
fit <- baseline
fit$params$eps_local <- 0
fit$params$eps_gross <- fit$params$eps
fit$params$tenure_contamination_model <- "local_gross"
fit$convergence <- 0L
fit$iterations <- 0L
fit$boundary_solution <- TRUE
gradient <- c(local_probability_boundary_score =
  (profile$loglik[1L] - baseline$loglik) / profile$local_rate[1L] /
    sum(df_fit$weight))
fit$gradient <- gradient
fit$gradient_max <- NA_real_
rates <- duration_weighted_transition_rates(df_fit, fit)[1L, ]
p <- fit$params
comparison <- data.frame(
  model = c("Joint per-wave marginal tenure", "Local plus gross tenure"),
  parameters = c(14L, 15L),
  loglik = c(baseline$loglik, fit$loglik),
  AIC = c(-2 * baseline$loglik + 28, -2 * fit$loglik + 30),
  alpha = c(baseline$params$alpha, p$alpha),
  pi = c(baseline$params$pi, p$pi),
  eps_total = c(baseline$params$eps, p$eps),
  eps_local = c(0, p$eps_local),
  eps_gross = c(baseline$params$eps, p$eps_gross),
  eps_d = c(baseline$params$eps_d, p$eps_d),
  effective_pair_contamination = c(1 - (1 - baseline$params$eps_d)^2,
    1 - (1 - p$eps_d)^2),
  weighted_exit = c(
    duration_weighted_transition_rates(df_fit, baseline)[1L, "exit_rate"],
    rates$exit_rate),
  weighted_entry = c(
    duration_weighted_transition_rates(df_fit, baseline)[1L, "entry_rate"],
    rates$entry_rate))
lr_stat <- max(0, 2 * (fit$loglik - baseline$loglik))
lr <- data.frame(comparison = "Local plus gross versus gross-only tenure error",
  LR = lr_stat, df = 1L,
  p_chibar2 = if (lr_stat == 0) 1 else .5 * pchisq(lr_stat, 1,
    lower.tail = FALSE))
multistart <- data.frame(start = seq_len(nrow(screen) + 1L),
  start_local_share = c(screen$local_share, 0),
  convergence = c(rep(NA_integer_, nrow(screen)), 0L),
  loglik = c(screen$loglik, fit$loglik), evaluations = 1L,
  error = NA_character_, phase = c(rep("fixed dispersed screen", nrow(screen)),
    "accepted boundary"))

cat("\nConditional local-share profile at the old optimum\n")
print(profile, row.names = FALSE, digits = 8)
cat("\nHybrid tenure model comparison\n")
print(comparison, row.names = FALSE, digits = 8)
cat("\nDispersed interior screen\n")
print(screen, row.names = FALSE, digits = 8)
cat("\nBoundary likelihood-ratio diagnostic\n")
print(lr, row.names = FALSE, digits = 8)
cat("\nMultistart diagnostics\n")
print(multistart, row.names = FALSE, digits = 8)
cat("\nScaled numerical gradient\n")
print(gradient, digits = 8)

write.csv(comparison, file.path(outdir, "model_comparison_latest.csv"),
  row.names = FALSE)
write.csv(lr, file.path(outdir, "boundary_lr_latest.csv"), row.names = FALSE)
write.csv(multistart, file.path(outdir, "screen_diagnostics_latest.csv"),
  row.names = FALSE)
fit$objective_function <- NULL
saveRDS(list(fit = fit, baseline = baseline, comparison = comparison,
  profile = profile, lr = lr, multistart = multistart),
  file.path(outdir, "fits_latest.rds"))
