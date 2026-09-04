# Efficient outer-product-of-gradients (OPG) analytical inference for the
# 30-parameter exact-retention model. This uses record-level likelihood scores
# and therefore needs 61, rather than 1,801, likelihood evaluations. It is a
# diagnostic analytical approximation, not the exact observed Hessian.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full, allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month", 1:3))

fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/fits_latest.rds")
fit <- readRDS(fit_file)$extension
z0 <- fit$par_unconstrained
if (length(z0) != 30L) stop("Expected a 30-parameter fit")
unpacked0 <- .piecewise_calendar_revision_monthly_unpack(z0,
  "joint_marginal")
if (!isTRUE(all.equal(unpacked0$tenure_exact_anchor_retention_prob,
    fit$params$tenure_exact_anchor_retention_prob, tolerance=1e-10)))
  stop("Saved constrained and unconstrained retention estimates disagree")
workers <- as.integer(Sys.getenv("EXACT_ANCHOR_OPG_WORKERS", "1"))
step_scale <- as.numeric(Sys.getenv("EXACT_ANCHOR_OPG_STEP", "0.001"))
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/inference_opg")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cluster <- parallel::makePSOCKcluster(workers)
on.exit(parallel::stopCluster(cluster), add=TRUE)
worker_wd <- getwd()
parallel::clusterCall(cluster, function(path) {
  setwd(path); source("EM-tenure/R/source_all.R"); NULL
}, worker_wd)
df_worker <- df_fit
parallel::clusterExport(cluster, "df_worker", envir=environment())

worker_eval <- function(z) {
  p <- .piecewise_calendar_revision_monthly_unpack(z, "joint_marginal")
  e <- e_step_eps(df_worker, p, check_df=FALSE, suff_stats=FALSE)
  proxy <- list(params=p, gamma=e$gamma,
    job_change_posterior=e$job_change_posterior)
  rates <- duration_weighted_transition_rates(df_worker, proxy)[1, ]
  jp <- e$job_change_posterior
  posterior_q <- sum(df_worker$weight*jp$expected_changes) /
    sum(df_worker$weight*jp$opportunities)
  quantities <- c(entry_rate=rates$entry_rate, exit_rate=rates$exit_rate,
    initial_employment=p$alpha, status_misclassification=p$pi,
    tenure_contamination=p$eps, timegap_contamination=p$eps_d,
    job_change_prob=p$job_change_prob,
    posterior_job_change_prob=posterior_q,
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
    setNames(p$lambda_g, paste0("exit_hazard_", 1:5)),
    setNames(p$lambda_d, paste0("entry_hazard_", 1:5)),
    setNames(p$tenure_start_month_probs, paste0("start_month_", 1:12)))
  list(row_loglik=e$row_loglik, quantities=quantities)
}
parallel::clusterExport(cluster, "worker_eval", envir=environment())

step <- step_scale*pmax(1, abs(z0))
p <- length(z0)
tasks <- vector("list", 2L*p+1L)
tasks[[1L]] <- z0
plus <- minus <- integer(p)
for (j in seq_len(p)) {
  zp <- zm <- z0
  zp[j] <- zp[j]+step[j]; zm[j] <- zm[j]-step[j]
  plus[j] <- 2L*j; minus[j] <- 2L*j+1L
  tasks[[plus[j]]] <- zp; tasks[[minus[j]]] <- zm
}
message("Evaluating ", length(tasks), " OPG score points")
values <- parallel::parLapplyLB(cluster, tasks, worker_eval)

ncell <- nrow(df_fit)
scores <- matrix(NA_real_, ncell, p,
  dimnames=list(NULL, names(z0)))
jacobian <- matrix(NA_real_, length(values[[1L]]$quantities), p,
  dimnames=list(names(values[[1L]]$quantities), names(z0)))
for (j in seq_len(p)) {
  scores[,j] <- (values[[plus[j]]]$row_loglik-
    values[[minus[j]]]$row_loglik)/(2*step[j])
  jacobian[,j] <- (values[[plus[j]]]$quantities-
    values[[minus[j]]]$quantities)/(2*step[j])
}
weighted_score <- scores*sqrt(df_fit$weight)
information <- crossprod(weighted_score)
information <- (information+t(information))/2
gradient <- colSums(scores*df_fit$weight)
eigenvalues <- eigen(information, symmetric=TRUE, only.values=TRUE)$values
tolerance <- max(abs(eigenvalues))*sqrt(.Machine$double.eps)
positive_definite <- min(eigenvalues) > tolerance
diagnostics <- data.frame(parameters=p,
  cells=ncell, rank=sum(eigenvalues>tolerance),
  positive_definite=positive_definite,
  minimum_eigenvalue=min(eigenvalues),
  maximum_eigenvalue=max(eigenvalues),
  condition_number=if (positive_definite)
    max(eigenvalues)/min(eigenvalues) else Inf,
  maximum_absolute_total_score=max(abs(gradient)),
  total_score_l2=sqrt(sum(gradient^2)), step_scale=step_scale)
if (!positive_definite) stop("OPG information is not positive definite")

vcov_z <- solve(information)
vcov_quantities <- jacobian%*%vcov_z%*%t(jacobian)
estimates <- values[[1L]]$quantities
se <- sqrt(pmax(diag(vcov_quantities), 0))
summary <- data.frame(quantity=names(estimates), estimate=unname(estimates),
  se=unname(se), lower=unname(estimates-1.96*se),
  upper=unname(estimates+1.96*se))
summary$near_boundary <- summary$estimate < 3*summary$se |
  summary$estimate > 1-3*summary$se
reporting_names <- c("status_misclassification", "tenure_contamination",
  "timegap_contamination", "job_change_prob",
  "gross_january_heaping_prob", "gross_anchor_revision_prob",
  "clean_anchor_revision_prob", "exact_anchor_retention_prob")
reporting_vcov <- vcov_quantities[reporting_names, reporting_names,
  drop=FALSE]
reporting_cor <- cov2cor(reporting_vcov)

result <- list(summary=summary, diagnostics=diagnostics, scores=scores,
  gradient=gradient, information=information, vcov_optimizer=vcov_z,
  jacobian=jacobian, vcov_quantities=vcov_quantities,
  reporting_vcov=reporting_vcov, reporting_cor=reporting_cor)
saveRDS(result, file.path(outdir, "opg_inference_latest.rds"))
write.csv(summary, file.path(outdir, "analytical_se_latest.csv"),
  row.names=FALSE)
write.csv(diagnostics, file.path(outdir, "information_diagnostics_latest.csv"),
  row.names=FALSE)
write.csv(data.frame(parameter=names(gradient), total_score=gradient),
  file.path(outdir, "total_score_latest.csv"), row.names=FALSE)
write.csv(reporting_cor,
  file.path(outdir, "reporting_probability_correlations_latest.csv"))

cat("\nExact-retention OPG analytical inference\n")
print(transform(summary, estimate=100*estimate, se=100*se,
  lower=100*lower, upper=100*upper), row.names=FALSE, digits=7)
cat("\nOPG and score diagnostics\n")
print(diagnostics, row.names=FALSE, digits=7)
cat("\nReporting-probability correlations\n")
print(round(reporting_cor, 3))
