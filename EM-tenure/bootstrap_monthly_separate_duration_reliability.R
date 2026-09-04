# Parametric-bootstrap calibration and full-nuisance recovery for the monthly
# two-clock reliability model. Null samples compare the common- and
# separate-dispersion fits; unrestricted samples re-estimate all 33 parameters.
# Every data set and optimizer stage is cached so larger runs are resumable.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

root <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability")
fit_file <- file.path(root,"fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate and polish the empirical models first")
empirical <- readRDS(fit_file)
truth_null <- empirical$common$params
truth_alt <- empirical$extension$params
truth_null$tenure_measurement_model <- "monthly"
truth_alt$tenure_measurement_model <- "monthly"

n <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_N","10000"))
reps <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_REPS","10"))
workers <- max(1L,min(8L,as.integer(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_WORKERS","8"))))
screen_maxit <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_SCREEN_MAXIT","1"))
refine_maxit <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_REFINE_MAXIT","5"))
refine_chunks <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_REFINE_CHUNKS","3"))
reltol <- as.numeric(Sys.getenv("RELIABILITY_BOOTSTRAP_RELTOL","1e-9"))
pgtol <- as.numeric(Sys.getenv("RELIABILITY_BOOTSTRAP_PGTOL","1e-7"))
resume <- identical(tolower(Sys.getenv("RELIABILITY_BOOTSTRAP_RESUME","true")),
  "true")
run_null <- identical(tolower(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_RUN_NULL","true")),"true")
run_recovery <- identical(tolower(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_RUN_RECOVERY","true")),"true")
recovery_blockwise <- identical(tolower(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_RECOVERY_BLOCKWISE","false")),"true")
block_cycles <- as.integer(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_BLOCK_CYCLES","1"))
block_maxit <- as.integer(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_BLOCK_MAXIT","8"))
block_full_maxit <- as.integer(Sys.getenv(
  "RELIABILITY_BOOTSTRAP_BLOCK_FULL_MAXIT","3"))
if (!run_null && !run_recovery) stop("At least one exercise must be enabled")
base_seed <- as.integer(Sys.getenv("RELIABILITY_BOOTSTRAP_SEED","261000"))
outdir <- file.path(root,"bootstrap_full_nuisance",
  paste0("n",n,"_b",reps))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

.compact_fit <- function(fit) {
  fit$gamma <- NULL
  fit$job_change_posterior <- NULL
  fit$objective_function <- NULL
  fit
}
.as_separate <- function(p,difference=0) {
  common <- p$duration_reliability_shift
  separate <- c(p$tenure_reliability_shift,p$timegap_reliability_shift)
  separate <- separate[is.finite(separate)]
  if (is.null(common) || !is.finite(common) ||
      (common==0 && length(separate)==2L)) common <- mean(separate)
  p$duration_reliability_shift <- 0
  p$tenure_reliability_shift <- max(0,common-difference/2)
  p$timegap_reliability_shift <- max(0,common+difference/2)
  p
}
.as_common <- function(p) {
  shifts <- c(p$tenure_reliability_shift,p$timegap_reliability_shift)
  shifts <- shifts[is.finite(shifts)]
  common <- p$duration_reliability_shift
  if (is.null(common) || !is.finite(common) || common==0)
    common <- if (length(shifts)) mean(shifts) else 1
  p$tenure_reliability_shift <- NULL
  p$timegap_reliability_shift <- NULL
  p$duration_reliability_shift <- common
  p
}
.perturb_start <- function(p,model) {
  p$eps <- min(.45,max(.01,p$eps*.75+.03))
  p$eps_d <- min(.45,max(.005,p$eps_d*1.35+.01))
  p$pi <- min(.20,max(.005,p$pi*1.25))
  p$job_change_prob <- min(.10,max(.001,p$job_change_prob*1.5))
  if (model=="common") {
    p$duration_reliability_shift <- min(5,
      max(.2,p$duration_reliability_shift+.7))
  } else {
    average <- mean(c(p$tenure_reliability_shift,
      p$timegap_reliability_shift))
    p$tenure_reliability_shift <- max(.1,average-.5)
    p$timegap_reliability_shift <- min(5,average+.5)
  }
  p
}
.fit_stage <- function(df_fit,start,model,path,label,maxit,free_names=NULL,
    gradient_scheme="forward") {
  if (resume && file.exists(path)) return(readRDS(path))
  message(label)
  p <- if (model=="common") .as_common(start) else .as_separate(start,
    if (!is.null(start$timegap_reliability_shift) &&
        !is.null(start$tenure_reliability_shift))
      start$timegap_reliability_shift-start$tenure_reliability_shift else 0)
  args <- list(df=df_fit,start=p,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    q_start=p$job_change_prob,maxit=maxit,reltol=reltol,pgtol=pgtol,
    workers=workers,verbose=0L,gradient_step=1e-4,
    gradient_scheme=gradient_scheme,free_names=free_names)
  if (model=="common") {
    args$reliability_shift_start <- p$duration_reliability_shift
  } else {
    args$tenure_reliability_shift_start <- p$tenure_reliability_shift
    args$timegap_reliability_shift_start <- p$timegap_reliability_shift
  }
  fit <- tryCatch(do.call(fit_eps_piecewise_calendar_revision_monthly,args),
    error=function(error) structure(list(error=conditionMessage(error),
      label=label,model=model),class="failed_full_recovery_fit"))
  if (!inherits(fit,"failed_full_recovery_fit")) {
    fit <- .compact_fit(fit); fit$label <- label; fit$model <- model
    fit$optimizer_control <- list(maxit=maxit,reltol=reltol,pgtol=pgtol,
      workers=workers,gradient_step=1e-4,
      gradient_scheme=gradient_scheme,free_names=free_names)
  }
  saveRDS(fit,path)
  fit
}
.fit_multistart <- function(df_fit,starts,model,prefix) {
  screens <- lapply(seq_along(starts),function(j) .fit_stage(df_fit,
    starts[[j]],model,file.path(outdir,paste0(prefix,"_screen",j,".rds")),
    paste(prefix,"screen",j),screen_maxit))
  ok <- !vapply(screens,inherits,logical(1),"failed_full_recovery_fit")
  if (!any(ok)) stop("All screens failed for ",prefix)
  candidates <- screens[ok]
  current <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),
    "loglik"))]]
  stages <- screens
  if (!identical(current$convergence,0L)) {
    for (chunk in seq_len(refine_chunks)) {
      refined <- .fit_stage(df_fit,current$params,model,
        file.path(outdir,paste0(prefix,"_refine",chunk,".rds")),
        paste(prefix,"refine",chunk),refine_maxit)
      stages[[length(stages)+1L]] <- refined
      if (inherits(refined,"failed_full_recovery_fit")) next
      gain <- refined$loglik-current$loglik
      if (refined$loglik>=current$loglik-1e-6) current <- refined
      if (identical(current$convergence,0L) || gain<1e-4) break
    }
  }
  list(best=current,stages=stages,overall_convergence=current$convergence)
}
.fit_recovery_blocks <- function(df_fit,result,prefix) {
  calendar_names <- c("tenure_heaping","tenure_year_revision",
    "tenure_clean_anchor_revision",paste0("start_month_logit",1:11),
    "tenure_exact_anchor_retention","tenure_local_revision")
  reliability_names <- c("pi","eps","eps_d","job_change",
    "tenure_reliability_shift","timegap_reliability_shift")
  transition_names <- c("alpha",paste0("log_hg",1:5),paste0("log_hd",1:5))
  current <- result$best
  stages <- result$stages
  for (cycle in seq_len(block_cycles)) {
    specifications <- list(calendar=list(names=calendar_names,scheme="forward"),
      reliability=list(names=reliability_names,scheme="central"),
      transition=list(names=transition_names,scheme="central"))
    for (block in names(specifications)) {
      specification <- specifications[[block]]
      fit <- .fit_stage(df_fit,current$params,"separate",
        file.path(outdir,paste0(prefix,"_block",cycle,"_",block,".rds")),
        paste(prefix,"block",cycle,block),block_maxit,
        free_names=specification$names,
        gradient_scheme=specification$scheme)
      stages[[length(stages)+1L]] <- fit
      if (!inherits(fit,"failed_full_recovery_fit") &&
          fit$loglik>=current$loglik-1e-6) current <- fit
    }
  }
  final <- .fit_stage(df_fit,current$params,"separate",
    file.path(outdir,paste0(prefix,"_block",block_cycles,"_full.rds")),
    paste(prefix,"block",block_cycles,"full"),block_full_maxit)
  stages[[length(stages)+1L]] <- final
  final_accepted <- !inherits(final,"failed_full_recovery_fit") &&
    final$loglik>=current$loglik-1e-6
  final_gain <- if (inherits(final,"failed_full_recovery_fit")) NA_real_ else
    final$loglik-current$loglik
  if (final_accepted) current <- final
  overall_convergence <- if (!final_accepted) 1L else final$convergence
  list(best=current,stages=stages,
    overall_convergence=overall_convergence,last_gain=final_gain)
}
.last_gain <- function(result) {
  ok <- !vapply(result$stages,inherits,logical(1),
    "failed_full_recovery_fit")
  values <- vapply(result$stages[ok],`[[`,numeric(1),"loglik")
  if (length(values)<2L) return(NA_real_)
  tail(values,1L)-tail(values,2L)[1L]
}
.simulate_cached <- function(kind,replication,truth,seed) {
  path <- file.path(outdir,sprintf("%s_rep%02d_data.rds",kind,replication))
  if (resume && file.exists(path)) return(readRDS(path))
  d <- simulate_eps_piecewise_job_change(n,truth,seed=seed)
  saveRDS(d,path); d
}
.stage_rows <- function(result,kind,replication,model) {
  do.call(rbind,lapply(seq_along(result$stages),function(j) {
    fit <- result$stages[[j]]
    if (inherits(fit,"failed_full_recovery_fit")) return(data.frame(
      kind=kind,replication=replication,model=model,stage=j,
      label=fit$label,loglik=NA_real_,convergence=NA_integer_,
      iterations=NA_integer_,maxit=NA_integer_,reltol=NA_real_,
      pgtol=NA_real_,error=fit$error))
    control <- fit$optimizer_control
    data.frame(kind=kind,replication=replication,model=model,stage=j,
      label=fit$label,loglik=fit$loglik,convergence=fit$convergence,
      iterations=fit$iterations,
      maxit=if (is.null(control)) NA_integer_ else control$maxit,
      reltol=if (is.null(control)) NA_real_ else control$reltol,
      pgtol=if (is.null(control)) NA_real_ else control$pgtol,
      error=NA_character_)
  }))
}
.natural_parameters <- function(p,df_fit) {
  evaluation <- e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)
  rates <- duration_weighted_transition_rates(df_fit,
    list(params=p,gamma=evaluation$gamma))[1,]
  shifts <- .duration_reliability_shifts(p)
  values <- c(entry_rate=rates$entry_rate,exit_rate=rates$exit_rate,
    initial_employment=p$alpha,status_error=p$pi,
    tenure_midpoint=p$eps,timegap_midpoint=p$eps_d,
    job_change=p$job_change_prob,
    tenure_heaping=p$tenure_heaping_prob,
    year_revision=p$tenure_year_revision_prob,
    clean_anchor_revision=p$tenure_clean_anchor_revision_prob,
    exact_anchor_retention=p$tenure_exact_anchor_retention_prob,
    local_revision=p$tenure_local_revision_prob,
    tenure_shift=unname(shifts["tenure"]),
    timegap_shift=unname(shifts["timegap"]),
    shift_difference=unname(shifts["timegap"]-shifts["tenure"]),
    setNames(p$lambda_g,paste0("exit_hazard_",seq_along(p$lambda_g))),
    setNames(p$lambda_d,paste0("entry_hazard_",seq_along(p$lambda_d))),
    setNames(p$tenure_start_month_probs,paste0("start_month_",1:12)))
  values
}

null_results <- vector("list",reps)
alt_results <- vector("list",reps)
stage_rows <- list()
for (replication in seq_len(reps)) {
  if (run_null) {
  message("Null bootstrap replication ",replication," of ",reps)
  d0 <- .simulate_cached("null",replication,truth_null,
    base_seed+replication)
  df0 <- collapse_eps_cells(d0,allow_zero_tenure=TRUE,
    extra_cols=paste0("interview_month",1:3))
  common_start <- .as_common(truth_null)
  common_result <- .fit_multistart(df0,list(common_start,
    .perturb_start(common_start,"common")),"common",
    sprintf("null_rep%02d_common",replication))
  separate_start <- .as_separate(common_result$best$params,0)
  separate_result <- .fit_multistart(df0,list(separate_start,
    .perturb_start(separate_start,"separate")),"separate",
    sprintf("null_rep%02d_separate",replication))
  LR <- max(0,2*(separate_result$best$loglik-common_result$best$loglik))
  null_results[[replication]] <- list(common=common_result$best,
    separate=separate_result$best,LR=LR,
    common_last_gain=.last_gain(common_result),
    separate_last_gain=.last_gain(separate_result))
  stage_rows[[length(stage_rows)+1L]] <- .stage_rows(common_result,"null",
    replication,"common")
  stage_rows[[length(stage_rows)+1L]] <- .stage_rows(separate_result,"null",
    replication,"separate")
  }

  if (run_recovery) {
  message("Unrestricted recovery replication ",replication," of ",reps)
  d1 <- .simulate_cached("alternative",replication,truth_alt,
    base_seed+10000L+replication)
  df1 <- collapse_eps_cells(d1,allow_zero_tenure=TRUE,
    extra_cols=paste0("interview_month",1:3))
  equal_start <- .as_separate(.as_common(truth_alt),0)
  alt_result <- .fit_multistart(df1,list(truth_alt,equal_start),"separate",
    sprintf("alternative_rep%02d_separate",replication))
  if (recovery_blockwise) alt_result <- .fit_recovery_blocks(df1,alt_result,
    sprintf("alternative_rep%02d_separate",replication))
  fitted <- .natural_parameters(alt_result$best$params,df1)
  generating <- .natural_parameters(truth_alt,df1)
  alt_results[[replication]] <- list(best=alt_result$best,
    fitted=fitted,generating=generating,
    last_gain=if (is.null(alt_result$last_gain)) .last_gain(alt_result) else
      alt_result$last_gain,
    convergence=alt_result$overall_convergence)
  stage_rows[[length(stage_rows)+1L]] <- .stage_rows(alt_result,
    "alternative",replication,"separate")
  }
}

observed_LR <- 2*(empirical$extension$loglik-empirical$common$loglik)
null_index <- which(!vapply(null_results,is.null,logical(1)))
null_draws <- if (length(null_index)) do.call(rbind,lapply(null_index,function(i) {
  result <- null_results[[i]]
  data.frame(replication=i,common_loglik=result$common$loglik,
    separate_loglik=result$separate$loglik,LR=result$LR,
    common_convergence=result$common$convergence,
    separate_convergence=result$separate$convergence,
    common_last_gain=result$common_last_gain,
    separate_last_gain=result$separate_last_gain,
    valid=identical(result$common$convergence,0L) &&
      identical(result$separate$convergence,0L))
})) else data.frame(replication=integer(),common_loglik=double(),
  separate_loglik=double(),LR=double(),common_convergence=integer(),
  separate_convergence=integer(),common_last_gain=double(),
  separate_last_gain=double(),valid=logical())
valid_null <- null_draws[null_draws$valid,,drop=FALSE]
valid_reps <- nrow(valid_null)
null_summary <- data.frame(reps_requested=reps,valid_reps=valid_reps,n=n,
  observed_LR=observed_LR,
  exceedances=if (valid_reps) sum(valid_null$LR>=observed_LR) else NA_integer_,
  bootstrap_p=if (valid_reps)
    (1+sum(valid_null$LR>=observed_LR))/(valid_reps+1) else NA_real_,
  LR_mean=if (valid_reps) mean(valid_null$LR) else NA_real_,
  LR_median=if (valid_reps) median(valid_null$LR) else NA_real_,
  LR_95=if (valid_reps) unname(quantile(valid_null$LR,.95)) else NA_real_,
  LR_max=if (valid_reps) max(valid_null$LR) else NA_real_)

alt_index <- which(!vapply(alt_results,is.null,logical(1)))
recovery_draws <- if (length(alt_index)) do.call(rbind,lapply(alt_index,function(i) {
  data.frame(replication=i,parameter=names(alt_results[[i]]$fitted),
    truth=unname(alt_results[[i]]$generating),
    estimate=unname(alt_results[[i]]$fitted),
    fit_convergence=alt_results[[i]]$convergence,
    last_loglik_gain=alt_results[[i]]$last_gain,
    gain_per_observation=alt_results[[i]]$last_gain/n,
    valid=identical(alt_results[[i]]$convergence,0L))
})) |> mutate(error=estimate-truth) else data.frame(replication=integer(),
  parameter=character(),truth=double(),estimate=double(),
  fit_convergence=integer(),last_loglik_gain=double(),
  gain_per_observation=double(),valid=logical(),error=double())
recovery_summary <- recovery_draws |>
  group_by(parameter) |>
  summarise(reps=n(),valid_reps=sum(valid),truth=first(truth),
    mean_estimate=if (any(valid)) mean(estimate[valid]) else NA_real_,
    bias=if (any(valid)) mean(error[valid]) else NA_real_,
    rmse=if (any(valid)) sqrt(mean(error[valid]^2)) else NA_real_,
    minimum=if (any(valid)) min(estimate[valid]) else NA_real_,
    maximum=if (any(valid)) max(estimate[valid]) else NA_real_,
    .groups="drop")
optimizer_stages <- if (length(stage_rows)) do.call(rbind,stage_rows) else
  data.frame(kind=character(),replication=integer(),model=character(),
    stage=integer(),label=character(),loglik=double(),convergence=integer(),
    iterations=integer(),maxit=integer(),reltol=double(),pgtol=double(),
    error=character())

write.csv(null_draws,file.path(outdir,"null_lr_draws_latest.csv"),
  row.names=FALSE)
write.csv(null_summary,file.path(outdir,"null_lr_summary_latest.csv"),
  row.names=FALSE)
write.csv(recovery_draws,file.path(outdir,"full_recovery_draws_latest.csv"),
  row.names=FALSE)
write.csv(recovery_summary,file.path(outdir,
  "full_recovery_summary_latest.csv"),row.names=FALSE)
write.csv(optimizer_stages,file.path(outdir,"optimizer_stages_latest.csv"),
  row.names=FALSE)
saveRDS(list(requested_configuration=list(n=n,reps=reps,seed=base_seed,
  workers=workers,screen_maxit=screen_maxit,refine_maxit=refine_maxit,
  refine_chunks=refine_chunks,reltol=reltol,pgtol=pgtol,resume=resume,
  run_null=run_null,run_recovery=run_recovery,
  recovery_blockwise=recovery_blockwise,block_cycles=block_cycles,
  block_maxit=block_maxit,block_full_maxit=block_full_maxit),
  empirical_LR=observed_LR,null_draws=null_draws,null_summary=null_summary,
  recovery_draws=recovery_draws,recovery_summary=recovery_summary,
  optimizer_stages=optimizer_stages,null_results=null_results,
  alt_results=alt_results),file.path(outdir,"bootstrap_latest.rds"))

cat("\nNull LR bootstrap pilot\n")
print(null_draws,row.names=FALSE,digits=8)
cat("\nNull LR summary\n")
print(null_summary,row.names=FALSE,digits=8)
cat("\nFull-nuisance recovery summary\n")
print(recovery_summary,row.names=FALSE,digits=8)
cat("\nOptimizer stages\n")
print(optimizer_stages,row.names=FALSE,digits=8)
if (run_null && valid_reps<reps) warning(reps-valid_reps,
  " null bootstrap fit(s) did not converge and were excluded")
invalid_recovery <- length(unique(recovery_draws$replication[!recovery_draws$valid]))
if (run_recovery && invalid_recovery) warning(invalid_recovery,
  " recovery fit(s) did not converge; aggregate recovery statistics are NA")
