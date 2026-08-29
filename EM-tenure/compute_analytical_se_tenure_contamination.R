# Numerical observed-information and delta-method inference for the six
# tenure-contamination specifications reported in Table 7.
#
# This script does not re-estimate point estimates. It reconstructs each saved
# observed-likelihood objective at its optimum, computes the full numerical
# Hessian on the unconstrained optimizer scale, and propagates its inverse to
# the paper-facing probabilities. The calculation follows the scaling used in
# estimate_timegap_contamination.R: the objective is the average negative
# survey-weighted log likelihood, so the inverse Hessian is divided by the
# total survey weight.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
total_weight <- sum(df_fit$weight)

duration_saved <- readRDS(
  "EM-tenure/output/results/duration_hazard/fits_latest.rds")
piecewise_saved <- readRDS(
  "EM-tenure/output/results/piecewise_hazard/fits_latest.rds")
robust_saved <- readRDS(
  "EM-tenure/output/results/timegap_contamination_robustness/fits_latest.rds")

model_specs <- list(
  list(label = "Constant, exact", fit = duration_saved$constant_direct,
       family = "duration", beta_fixed = c(beta_g = 0, beta_d = 0)),
  list(label = "Power law, exact", fit = duration_saved$duration,
       family = "duration", beta_fixed = NULL),
  list(label = "Piecewise, exact", fit = piecewise_saved$piecewise,
       family = "piecewise", timegap_model = NULL),
  list(label = "Independent error", fit = robust_saved$marginal,
       family = "piecewise", timegap_model = "marginal"),
  list(label = "Local error", fit = robust_saved$local$fit,
       family = "piecewise", timegap_model = "local"),
  list(label = "Joint per-wave error", fit = robust_saved$joint$fit,
       family = "piecewise", timegap_model = "joint_marginal")
)

make_transform <- function(spec) {
  if (identical(spec$family, "duration")) {
    beta_fixed <- spec$beta_fixed
    z0 <- .duration_eps_pack(spec$fit$params, beta_fixed = beta_fixed)
    unpack <- function(z) .duration_eps_unpack(z, beta_fixed = beta_fixed)
  } else {
    contaminated <- !is.null(spec$timegap_model)
    z0 <- .piecewise_eps_pack(spec$fit$params,
      timegap_contamination = contaminated)
    unpack <- function(z) {
      p <- .piecewise_eps_unpack(z,
        timegap_contamination = contaminated)
      if (contaminated) {
        p$timegap_contamination_model <- spec$timegap_model
        p$timegap_local_decay <- if (identical(spec$timegap_model, "local")) 1 else 1
      }
      p
    }
  }
  list(z0 = z0, unpack = unpack)
}

evaluate_model <- function(z, unpack, quantities_only = FALSE) {
  p <- unpack(z)
  estep <- e_step_eps(df_fit, p, check_df = FALSE, suff_stats = FALSE)
  if (!quantities_only) return(-estep$loglik / total_weight)
  fit <- list(params = p, gamma = estep$gamma)
  rates <- duration_weighted_transition_rates(df_fit, fit)[1L, ]
  eps_d <- if (is.null(p$eps_d)) 0 else p$eps_d
  pair <- if (identical(p$timegap_contamination_model, "joint_marginal"))
    1 - (1 - eps_d)^2 else eps_d
  c(entry_rate = rates$entry_rate, exit_rate = rates$exit_rate,
    misclassification_rate = p$pi, tenure_contamination_rate = p$eps,
    timegap_contamination_rate = eps_d,
    comparable_pair_contamination_rate = pair,
    initial_employment_rate = p$alpha)
}

numeric_jacobian <- function(fn, x, step = 1e-3) {
  base <- fn(x)
  out <- matrix(NA_real_, nrow = length(base), ncol = length(x),
                dimnames = list(names(base), names(x)))
  for (j in seq_along(x)) {
    xp <- xm <- x
    xp[j] <- xp[j] + step
    xm[j] <- xm[j] - step
    out[, j] <- (fn(xp) - fn(xm)) / (2 * step)
  }
  out
}

infer_one <- function(spec) {
  transform <- make_transform(spec)
  z0 <- transform$z0
  objective <- function(z) evaluate_model(z, transform$unpack, FALSE)
  quantities <- function(z) evaluate_model(z, transform$unpack, TRUE)
  message("Computing full observed-information Hessian: ", spec$label,
          " (", length(z0), " parameters)")
  hessian <- optimHess(z0, objective,
    control = list(ndeps = rep(1e-3, length(z0))))
  hessian <- (hessian + t(hessian)) / 2
  eigenvalues <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- max(abs(eigenvalues)) * sqrt(.Machine$double.eps)
  if (sum(eigenvalues > tolerance) != length(z0) || min(eigenvalues) <= 0)
    stop("Observed information is not positive definite for ", spec$label)
  vcov_z <- solve(hessian) / total_weight
  jacobian <- numeric_jacobian(quantities, z0)
  estimates <- quantities(z0)
  variance <- jacobian %*% vcov_z %*% t(jacobian)
  standard_errors <- sqrt(pmax(diag(variance), 0))
  if (is.null(spec$timegap_model)) {
    standard_errors["timegap_contamination_rate"] <- NA_real_
    standard_errors["comparable_pair_contamination_rate"] <- NA_real_
  }
  list(
    summary = data.frame(model = spec$label, quantity = names(estimates),
      estimate = unname(estimates), se = unname(standard_errors)),
    diagnostics = data.frame(model = spec$label, parameters = length(z0),
      rank = sum(eigenvalues > tolerance),
      minimum_eigenvalue = min(eigenvalues),
      condition_number = max(eigenvalues) / min(eigenvalues)),
    vcov_optimizer = vcov_z,
    jacobian = jacobian
  )
}

fits <- lapply(model_specs, infer_one)
summary_table <- do.call(rbind, lapply(fits, `[[`, "summary"))
diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))

outdir <- "EM-tenure/output/results/analytical_inference"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(outdir, "table7_analytical_se_latest.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(outdir, "table7_hessian_diagnostics_latest.csv"),
          row.names = FALSE)
saveRDS(list(summary = summary_table, diagnostics = diagnostics,
             models = fits),
        file.path(outdir, "table7_analytical_inference_latest.rds"))

cat("\nTable 7 observed-information/delta-method inference\n")
print(transform(summary_table, estimate = 100 * estimate, se = 100 * se),
      row.names = FALSE, digits = 6)
cat("\nHessian diagnostics\n")
print(diagnostics, row.names = FALSE, digits = 6)
