# Nested Table 5 extension: AR(2) + Set 4 transitions + reliability-dependent error.
if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source EM-AR2/estimate_set4_reliability_4w.R from the project root")
source("EM-AR1-4W/R/source_all.R")
source("EM-AR2/R/source_all.R")

set.seed(20260827L)
results_dir <- "EM-AR2/output/results/set4_reliability"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

full_data <- prepare_ar2_set4_reliability_4w()
cat(sprintf("Table 5 extension sample: N=%s; exact cells=%s\n",
  format(full_data$n_original, big.mark=","),
  format(full_data$n_cells, big.mark=",")))

duration_fit <- readRDS("EM-AR2/output/results/ar2_duration_covariates_4w_latest.rds")
table4_fit <- readRDS("EM-baseline-ext/output/results/fit_cov_s4_reliability_free.rds")

start_from_duration <- function(data) {
  out <- setNames(rep(0, length(ar2_reliability_names(data))),
                  ar2_reliability_names(data))
  out[1:3] <- duration_fit$eta[1:3]
  for (block in c("entry", "persistence")) {
    old <- duration_fit$params[[if (block == "entry") "beta0" else "beta1"]]
    target <- paste0(block, "_", names(old))
    common <- intersect(target, names(out))
    out[common] <- old[sub(paste0("^", block, "_"), "", common)]
  }
  out["error_intercept"] <- qlogis(2 * duration_fit$params$pi)
  out
}

expand_start <- function(old_eta, data, reliability_warm = FALSE) {
  out <- setNames(rep(0, length(ar2_reliability_names(data))),
                  ar2_reliability_names(data))
  common <- intersect(names(old_eta), names(out)); out[common] <- old_eta[common]
  if (reliability_warm) {
    common_delta <- intersect(names(table4_fit$params$delta), names(out))
    out[common_delta] <- table4_fit$params$delta[common_delta]
  }
  out
}

stage_names <- c("constant", "inconsistency", "reliability")
requested <- Sys.getenv("AR2_EXTENSION_STAGES", paste(stage_names, collapse=","))
stage_names <- intersect(stage_names, trimws(strsplit(requested, ",", fixed=TRUE)[[1]]))
if (!length(stage_names)) stop("AR2_EXTENSION_STAGES selected no valid stages")
n_starts <- as.integer(Sys.getenv("AR2_EXTENSION_STARTS", "3"))
if (!is.finite(n_starts) || n_starts < 1L || n_starts > 5L)
  stop("AR2_EXTENSION_STARTS must be between 1 and 5")
resume <- identical(Sys.getenv("AR2_EXTENSION_RESUME"), "1")

best <- list(); inference <- list(); previous <- NULL
for (stage in stage_names) {
  data <- subset_ar2_reliability_stage(full_data, stage)
  accepted_checkpoint <- file.path(results_dir,
    sprintf("fit_ar2_set4_%s_latest.rds", stage))
  if (resume && file.exists(accepted_checkpoint)) {
    fit <- readRDS(accepted_checkpoint)
    if (length(fit$eta) != length(ar2_reliability_names(data)))
      stop("Accepted checkpoint has the wrong parameter dimension for ", stage)
    inf <- fit$analytical_inference
    if (is.null(inf)) stop("Accepted checkpoint lacks analytical inference: ", stage)
    best[[stage]] <- fit; inference[[stage]] <- inf; previous <- fit
    cat(sprintf("\n=== Stage: %s | accepted checkpoint ===\n", stage))
    cat(sprintf("resumed accepted %s: ll=%.3f score=%.3e rank=%d/%d\n",
      stage, fit$loglik, fit$max_abs_score,
      inf$diagnostics$information_rank, inf$diagnostics$parameter_count))
    next
  }
  nested <- if (is.null(previous)) start_from_duration(data) else
    expand_start(previous$eta, data, FALSE)
  warm <- if (stage == "constant") nested else expand_start(previous$eta, data, TRUE)
  starts <- list(nested = nested)
  if (n_starts >= 2L) {
    if (stage == "constant") {
      scale <- ifelse(grepl("alpha|intercept|lag2", names(warm)), .08, .04)
      starts$perturbed_1 <- warm + rnorm(length(warm), 0, scale)
    } else starts$reliability_warm <- warm
  }
  if (n_starts >= 3L) for (j in 3:n_starts) {
    scale <- ifelse(grepl("inconsistency|panel_|error_", names(warm)), .12, .05)
    starts[[paste0("perturbed_", length(starts))]] <-
      warm + rnorm(length(warm), 0, scale)
  }
  fits <- vector("list", length(starts)); names(fits) <- names(starts)
  elapsed <- numeric(length(starts))
  cat(sprintf("\n=== Stage: %s | K=%d | starts=%d ===\n",
              stage, length(warm), length(starts)))
  for (i in seq_along(starts)) {
    checkpoint <- file.path(results_dir,
      sprintf("fit_ar2_set4_%s_start_%s.rds", stage, names(starts)[i]))
    if (resume && file.exists(checkpoint)) {
      fits[[i]] <- readRDS(checkpoint)
      cat(sprintf("resumed %-18s ll=%.3f score=%.3e\n", names(starts)[i],
                  fits[[i]]$loglik, fits[[i]]$max_abs_score))
      next
    }
    preliminary_maxit <- as.integer(Sys.getenv("AR2_EXTENSION_PRELIM_MAXIT", "250"))
    tm <- system.time(fits[[i]] <- fit_ar2_set4_reliability(
      data, starts[[i]], maxit=preliminary_maxit, reltol=1e-7))
    elapsed[i] <- unname(tm["elapsed"]); saveRDS(fits[[i]], checkpoint)
    cat(sprintf("%-25s ll=%.3f code=%d eval=%d score=%.3e time=%.1fs\n",
      names(starts)[i], fits[[i]]$loglik, fits[[i]]$optimizer_code,
      fits[[i]]$iterations, fits[[i]]$max_abs_score, elapsed[i]))
  }
  tab <- data.frame(stage=stage, start=names(fits),
    loglik=vapply(fits, `[[`, numeric(1), "loglik"),
    code=vapply(fits, `[[`, integer(1), "optimizer_code"),
    evaluations=vapply(fits, `[[`, numeric(1), "iterations"),
    max_abs_score=vapply(fits, `[[`, numeric(1), "max_abs_score"),
    elapsed_seconds=elapsed)
  write.csv(tab, file.path(results_dir,
    sprintf("multistart_ar2_set4_%s.csv", stage)), row.names=FALSE)
  winner <- which.max(tab$loglik); fit <- fits[[winner]]
  if (!is.null(previous)) {
    nested_ll <- sum(data$weight * .ar2r_components(nested, data)$loglik_i)
    if (abs(nested_ll - previous$loglik) > 1e-4)
      stop("Nested likelihood does not reproduce the preceding stage")
    if (fit$loglik < previous$loglik - 1e-4)
      stop("Likelihood nesting failed at stage ", stage)
  }
  cat("Computing analytical inference for ", stage, "...\n", sep="")
  inf <- analytical_se_ar2_set4_reliability(data, fit)
  correction <- inf$diagnostics$newton_correction
  if (max(abs(correction)) > 1e-3) {
    fractions <- c(1, .5, .25, .125, 0)
    trials <- lapply(fractions, function(frac) fit$eta + frac * correction)
    ll <- vapply(trials, function(z)
      sum(data$weight * .ar2r_components(z, data)$loglik_i), numeric(1))
    j <- which.max(ll)
    if (ll[j] > fit$loglik + 1e-5) {
      cat(sprintf("Newton polish fraction %.3f; gain %.3f\n",
                  fractions[j], ll[j] - fit$loglik))
      saveRDS(list(eta=trials[[j]], fraction=fractions[j], loglik=ll[j],
                   source_fit=fit$eta, inference=inf),
        file.path(results_dir,
          sprintf("newton_start_ar2_set4_%s_latest.rds", stage)))
      fit <- fit_ar2_set4_reliability(data, trials[[j]], 400L, 1e-10)
      inf <- analytical_se_ar2_set4_reliability(data, fit)
    }
  }
  if (fit$max_abs_score > 2e-5)
    stop("Final normalized score remains too large at stage ", stage,
         "; resume from the saved best start with a longer polish")
  if (inf$diagnostics$information_rank < inf$diagnostics$parameter_count)
    stop("Local identification failed at stage ", stage)
  fit$analytical_inference <- inf; fit$multistart <- tab
  saveRDS(fit, file.path(results_dir,
    sprintf("fit_ar2_set4_%s_latest.rds", stage)))
  write.csv(inf$coefficient_summary, file.path(results_dir,
    sprintf("coefficients_ar2_set4_%s.csv", stage)), row.names=FALSE)
  write.csv(inf$probability_summary, file.path(results_dir,
    sprintf("probabilities_ar2_set4_%s.csv", stage)), row.names=FALSE)
  best[[stage]] <- fit; inference[[stage]] <- inf; previous <- fit
  cat(sprintf("accepted %s: ll=%.3f score=%.3e rank=%d/%d cond=%.3g\n",
    stage, fit$loglik, fit$max_abs_score,
    inf$diagnostics$information_rank, inf$diagnostics$parameter_count,
    inf$diagnostics$information_condition))
  print(inf$probability_summary, row.names=FALSE, digits=6)
}

summary <- do.call(rbind, lapply(names(best), function(stage) {
  q <- inference[[stage]]$probability_summary
  data.frame(stage=stage, quantity=q$quantity, estimate=q$estimate, se=q$se,
    loglik=best[[stage]]$loglik, N=best[[stage]]$n_obs,
    K=inference[[stage]]$diagnostics$parameter_count,
    rank=inference[[stage]]$diagnostics$information_rank,
    condition=inference[[stage]]$diagnostics$information_condition)
}))
write.csv(summary, file.path(results_dir, "table5_extension_summary.csv"), row.names=FALSE)
cat("\n=== Table 5 extension summary ===\n")
print(summary, row.names=FALSE, digits=6)
