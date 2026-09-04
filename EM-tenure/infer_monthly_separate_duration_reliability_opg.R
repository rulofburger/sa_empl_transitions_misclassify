# Full-parameter outer-product-of-gradients (OPG) inference for the model with
# separate tenure and timegap reliability dispersions. This is a provisional
# analytical approximation: it adjusts the focal reliability parameters for
# all fitted nuisance parameters but is not the exact observed Hessian.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month",1:3))

root <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability")
saved <- readRDS(file.path(root,"fits_latest.rds"))
fit <- saved$extension
z0 <- fit$par_unconstrained
unpacked0 <- .piecewise_calendar_revision_monthly_unpack(z0,
  "joint_marginal")
if (!isTRUE(all.equal(unpacked0$tenure_reliability_shift,
    fit$params$tenure_reliability_shift,tolerance=1e-10)) ||
    !isTRUE(all.equal(unpacked0$timegap_reliability_shift,
    fit$params$timegap_reliability_shift,tolerance=1e-10)))
  stop("Saved constrained and unconstrained reliability estimates disagree")

workers <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_OPG_WORKERS","8"))
step_scale <- as.numeric(Sys.getenv("SEPARATE_RELIABILITY_OPG_STEP","0.001"))
outdir <- file.path(root,"inference_opg")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

cluster <- parallel::makePSOCKcluster(max(1L,min(workers,8L)))
on.exit(parallel::stopCluster(cluster),add=TRUE)
worker_wd <- getwd()
parallel::clusterCall(cluster,function(path) {
  setwd(path); source("EM-tenure/R/source_all.R"); NULL
},worker_wd)
df_worker <- df_fit
parallel::clusterExport(cluster,"df_worker",envir=environment())

worker_eval <- function(z) {
  tryCatch({
    p <- .piecewise_calendar_revision_monthly_unpack(z,"joint_marginal")
    e <- e_step_eps(df_worker,p,check_df=FALSE,suff_stats=FALSE)
    proxy <- list(params=p,gamma=e$gamma,
      job_change_posterior=e$job_change_posterior)
    rates <- duration_weighted_transition_rates(df_worker,proxy)[1,]
    components <- c(
      reliable_tenure=.duration_reliability_component_params(p,
        "reliable")$eps,
      unreliable_tenure=.duration_reliability_component_params(p,
        "unreliable")$eps,
      reliable_timegap=.duration_reliability_component_params(p,
        "reliable")$eps_d,
      unreliable_timegap=.duration_reliability_component_params(p,
        "unreliable")$eps_d)
    jp <- e$job_change_posterior
    posterior_q <- sum(df_worker$weight*jp$expected_changes)/
      sum(df_worker$weight*jp$opportunities)
    quantities <- c(entry_rate=rates$entry_rate,exit_rate=rates$exit_rate,
      initial_employment=p$alpha,status_misclassification=p$pi,
      tenure_contamination_midpoint=p$eps,
      timegap_contamination_midpoint=p$eps_d,components,
      tenure_reliability_shift=p$tenure_reliability_shift,
      timegap_reliability_shift=p$timegap_reliability_shift,
      shift_difference=p$timegap_reliability_shift-
        p$tenure_reliability_shift,
      job_change_prob=p$job_change_prob,
      posterior_job_change_prob=posterior_q,
      gross_january_heaping_prob=p$tenure_heaping_prob,
      gross_anchor_revision_prob=p$tenure_year_revision_prob,
      clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
      exact_anchor_retention_prob=p$tenure_exact_anchor_retention_prob,
      local_anchor_revision_prob=p$tenure_local_revision_prob,
      setNames(p$lambda_g,paste0("exit_hazard_",1:5)),
      setNames(p$lambda_d,paste0("entry_hazard_",1:5)))
    list(row_loglik=e$row_loglik,quantities=quantities,error=NA_character_)
  },error=function(error) list(row_loglik=NULL,quantities=NULL,
    error=conditionMessage(error)))
}
parallel::clusterExport(cluster,"worker_eval",envir=environment())

step <- step_scale*pmax(1,abs(z0))
p <- length(z0)
tasks <- vector("list",2L*p+1L); tasks[[1L]] <- z0
plus <- minus <- integer(p)
for (j in seq_len(p)) {
  zp <- zm <- z0
  zp[j] <- zp[j]+step[j]; zm[j] <- zm[j]-step[j]
  plus[j] <- 2L*j; minus[j] <- 2L*j+1L
  tasks[[plus[j]]] <- zp; tasks[[minus[j]]] <- zm
}
point_file <- file.path(outdir,"score_points_latest.rds")
if (file.exists(point_file)) {
  cached <- readRDS(point_file)
  cache_matches <- isTRUE(all.equal(cached$z0,z0,tolerance=0)) &&
    isTRUE(all.equal(cached$step,step,tolerance=0))
} else cache_matches <- FALSE
if (cache_matches) {
  message("Loading ",length(tasks)," cached OPG score points")
  values <- cached$values
} else {
  message("Evaluating ",length(tasks)," OPG score points")
  values <- parallel::parLapplyLB(cluster,tasks,worker_eval)
  saveRDS(list(z0=z0,step=step,values=values),point_file)
}
valid <- vapply(values,function(value) is.null(value$error) ||
  is.na(value$error),logical(1))
if (!valid[1L]) stop("The central likelihood evaluation failed: ",
  values[[1L]]$error)

ncell <- nrow(df_fit)
scores <- matrix(NA_real_,ncell,p,dimnames=list(NULL,names(z0)))
jacobian <- matrix(NA_real_,length(values[[1L]]$quantities),p,
  dimnames=list(names(values[[1L]]$quantities),names(z0)))
adaptive_scale <- rep(1,p); names(adaptive_scale) <- names(z0)
derivative_scheme <- rep("central",p); names(derivative_scheme) <- names(z0)
for (j in seq_len(p)) {
  plus_value <- values[[plus[j]]]; minus_value <- values[[minus[j]]]
  plus_valid <- valid[plus[j]]; minus_valid <- valid[minus[j]]
  h <- step[j]
  if (!plus_valid && !minus_valid) {
    for (scale in c(.1,.01,.001)) {
      zp <- zm <- z0; h <- step[j]*scale
      zp[j] <- zp[j]+h; zm[j] <- zm[j]-h
      adapted <- parallel::parLapplyLB(cluster,list(zp,zm),worker_eval)
      plus_value <- adapted[[1L]]; minus_value <- adapted[[2L]]
      plus_valid <- is.null(plus_value$error) || is.na(plus_value$error)
      minus_valid <- is.null(minus_value$error) || is.na(minus_value$error)
      if (plus_valid || minus_valid) {
        adaptive_scale[j] <- scale
        break
      }
    }
  }
  if (plus_valid && minus_valid) {
    if (adaptive_scale[j]<1) derivative_scheme[j] <- "adaptive central"
    scores[,j] <- (plus_value$row_loglik-minus_value$row_loglik)/(2*h)
    jacobian[,j] <- (plus_value$quantities-minus_value$quantities)/(2*h)
  } else if (plus_valid) {
    derivative_scheme[j] <- "forward"
    scores[,j] <- (plus_value$row_loglik-
      values[[1L]]$row_loglik)/h
    jacobian[,j] <- (plus_value$quantities-
      values[[1L]]$quantities)/h
  } else if (minus_valid) {
    derivative_scheme[j] <- "backward"
    scores[,j] <- (values[[1L]]$row_loglik-
      minus_value$row_loglik)/h
    jacobian[,j] <- (values[[1L]]$quantities-
      minus_value$quantities)/h
  } else stop("Both finite-difference sides failed for ",names(z0)[j],
    " after adaptive step reduction")
}
failed_points <- which(!valid)
failed_parameters <- if (length(failed_points)) unique(vapply(failed_points,
  function(index) names(z0)[ceiling((index-1L)/2)],character(1))) else
  character()
weighted_scores <- scores*sqrt(df_fit$weight)
information <- crossprod(weighted_scores)
information <- (information+t(information))/2
gradient <- colSums(scores*df_fit$weight)
eigenvalues <- eigen(information,symmetric=TRUE,only.values=TRUE)$values
tolerance <- max(abs(eigenvalues))*sqrt(.Machine$double.eps)
rank <- sum(eigenvalues>tolerance)
positive_definite <- min(eigenvalues)>tolerance
diagnostics <- data.frame(parameters=p,cells=ncell,rank=rank,
  positive_definite=positive_definite,minimum_eigenvalue=min(eigenvalues),
  maximum_eigenvalue=max(eigenvalues),
  condition_number=if (positive_definite)
    max(eigenvalues)/min(eigenvalues) else Inf,
  maximum_absolute_total_score=max(abs(gradient)),
  total_score_l2=sqrt(sum(gradient^2)),step_scale=step_scale,
  failed_perturbations=length(failed_points),
  originally_failed_parameters=paste(failed_parameters,collapse=";"),
  one_sided_parameters=paste(names(derivative_scheme)[
    derivative_scheme %in% c("forward","backward")],collapse=";"),
  adaptively_reduced_parameters=paste(names(adaptive_scale)[
    adaptive_scale<1],collapse=";"),
  minimum_step_scale=min(adaptive_scale))
if (!positive_definite) stop("OPG information is not positive definite")

vcov_z <- solve(information)
vcov_quantities <- jacobian%*%vcov_z%*%t(jacobian)
estimates <- values[[1L]]$quantities
se <- sqrt(pmax(diag(vcov_quantities),0))
summary <- data.frame(quantity=names(estimates),estimate=unname(estimates),
  se=unname(se),lower=unname(estimates-1.96*se),
  upper=unname(estimates+1.96*se))
focal <- c("tenure_contamination_midpoint",
  "timegap_contamination_midpoint","reliable_tenure",
  "unreliable_tenure","reliable_timegap","unreliable_timegap",
  "tenure_reliability_shift","timegap_reliability_shift",
  "shift_difference")
focal_cor <- cov2cor(vcov_quantities[focal,focal,drop=FALSE])

result <- list(summary=summary,diagnostics=diagnostics,
  gradient=gradient,information=information,vcov_optimizer=vcov_z,
  jacobian=jacobian,vcov_quantities=vcov_quantities,
  focal_correlations=focal_cor,step=step,
  adaptive_scale=adaptive_scale,derivative_scheme=derivative_scheme)
saveRDS(result,file.path(outdir,"opg_inference_latest.rds"))
write.csv(summary,file.path(outdir,"analytical_se_latest.csv"),row.names=FALSE)
write.csv(diagnostics,file.path(outdir,
  "information_diagnostics_latest.csv"),row.names=FALSE)
write.csv(data.frame(parameter=names(gradient),total_score=gradient),
  file.path(outdir,"total_score_latest.csv"),row.names=FALSE)
write.csv(focal_cor,file.path(outdir,
  "focal_parameter_correlations_latest.csv"))
cat("\nSeparate-duration-reliability OPG analytical inference\n")
print(summary,row.names=FALSE,digits=8)
cat("\nOPG and score diagnostics\n")
print(diagnostics,row.names=FALSE,digits=8)
cat("\nFocal-parameter correlations\n"); print(round(focal_cor,3))
