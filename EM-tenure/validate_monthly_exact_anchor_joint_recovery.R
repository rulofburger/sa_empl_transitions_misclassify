# Full joint simulation recovery for the 30-parameter exact-anchor model.
# Every transition, reporting, duration-hazard, and baseline-month parameter is
# re-estimated. Two starts per replication diagnose mode dependence; completed
# fits are cached so the grid is resumable.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")
fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate exact-anchor retention first")
truth <- readRDS(fit_file)$extension$params
sample_sizes <- as.integer(strsplit(Sys.getenv(
  "EXACT_ANCHOR_JOINT_RECOVERY_N","5000"),",",fixed=TRUE)[[1]])
reps <- as.integer(Sys.getenv("EXACT_ANCHOR_JOINT_RECOVERY_REPS","2"))
workers <- as.integer(Sys.getenv("EXACT_ANCHOR_JOINT_RECOVERY_WORKERS","4"))
maxit <- as.integer(Sys.getenv("EXACT_ANCHOR_JOINT_RECOVERY_MAXIT","40"))
central_maxit <- as.integer(Sys.getenv(
  "EXACT_ANCHOR_JOINT_RECOVERY_CENTRAL_MAXIT","8"))
resume <- identical(tolower(Sys.getenv(
  "EXACT_ANCHOR_JOINT_RECOVERY_RESUME","true")),"true")
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/validation_joint")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

fit_one <- function(df_fit,start,label,path,scheme="forward",iterations=maxit) {
  if (resume && file.exists(path)) return(readRDS(path))
  message("Joint recovery fit: ",label)
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start,
    heaping_start=max(start$tenure_heaping_prob,1e-6),
    revision_start=start$tenure_year_revision_prob,
    clean_anchor_revision_start=start$tenure_clean_anchor_revision_prob,
    start_month_probs_start=start$tenure_start_month_probs,
    exact_anchor_retention_start=start$tenure_exact_anchor_retention_prob,
    q_start=start$job_change_prob,maxit=iterations,reltol=1e-9,pgtol=1e-7,
    workers=max(1L,min(workers,8L)),verbose=1L,gradient_step=1e-3,
    gradient_scheme=scheme)
  fit$objective_function <- NULL
  fit$recovery_label <- label
  saveRDS(fit,path)
  fit
}

parameter_vector <- function(p) c(alpha=p$alpha,pi=p$pi,eps=p$eps,
  eps_d=p$eps_d,job_change_prob=p$job_change_prob,
  january_heaping=p$tenure_heaping_prob,
  gross_anchor_revision=p$tenure_year_revision_prob,
  clean_anchor_revision=p$tenure_clean_anchor_revision_prob,
  exact_anchor_retention=p$tenure_exact_anchor_retention_prob,
  setNames(p$lambda_g,paste0("exit_hazard_",1:5)),
  setNames(p$lambda_d,paste0("entry_hazard_",1:5)),
  setNames(p$tenure_start_month_probs,paste0("start_month_",1:12)))
truth_vector <- parameter_vector(truth)

fit_records <- list(); parameter_records <- list(); summary_records <- list()
counter <- 0L
for (n in sample_sizes) for (replication in seq_len(reps)) {
  seed <- 260918L+1000L*n+replication
  d <- simulate_eps_piecewise_job_change(n,truth,seed=seed)
  df_fit <- collapse_eps_cells(d,allow_zero_tenure=TRUE,
    extra_cols=paste0("interview_month",1:3))
  truth_eval <- e_step_eps(df_fit,truth,check_df=FALSE,suff_stats=FALSE)
  truth_fit <- list(params=truth,gamma=truth_eval$gamma,
    job_change_posterior=truth_eval$job_change_posterior,
    loglik=truth_eval$loglik)
  truth_rates <- duration_weighted_transition_rates(df_fit,truth_fit)[1,]
  repdir <- file.path(outdir,paste0("n",n,"_rep",replication))
  dir.create(repdir,recursive=TRUE,showWarnings=FALSE)

  near <- truth
  dispersed <- truth
  dispersed$tenure_exact_anchor_retention_prob <-
    if (replication%%2L) .04 else .40
  dispersed$eps <- if (replication%%2L) .12 else .38
  dispersed$tenure_clean_anchor_revision_prob <-
    if (replication%%2L) .03 else .30
  dispersed$tenure_year_revision_prob <-
    if (replication%%2L) .02 else .20
  dispersed$tenure_start_month_probs <- .5*truth$tenure_start_month_probs+
    .5/12
  dispersed$tenure_start_month_probs <-
    dispersed$tenure_start_month_probs/
      sum(dispersed$tenure_start_month_probs)

  near_fit <- fit_one(df_fit,near,paste0("n",n," rep",replication,
    " near"),file.path(repdir,"fit_near_latest.rds"))
  dispersed_fit <- fit_one(df_fit,dispersed,paste0("n",n," rep",replication,
    " dispersed"),file.path(repdir,"fit_dispersed_latest.rds"))
  candidates <- list(near=near_fit,dispersed=dispersed_fit)
  selected_name <- names(candidates)[which.max(vapply(candidates,`[[`,
    numeric(1),"loglik"))]
  selected <- candidates[[selected_name]]
  central <- fit_one(df_fit,selected$params,paste0("n",n," rep",replication,
    " central"),file.path(repdir,"fit_central_latest.rds"),"central",
    central_maxit)
  candidates$central <- central
  best_name <- names(candidates)[which.max(vapply(candidates,`[[`,numeric(1),
    "loglik"))]
  best <- candidates[[best_name]]
  saveRDS(list(truth=truth,data_seed=seed,near=near_fit,
    dispersed=dispersed_fit,central=central,best_name=best_name,best=best),
    file.path(repdir,"recovery_latest.rds"))

  for (label in names(candidates)) {
    counter <- counter+1L
    f <- candidates[[label]]
    v <- parameter_vector(f$params)
    parameter_records[[counter]] <- data.frame(n=n,replication=replication,
      start=label,selected=label==best_name,parameter=names(v),
      truth=unname(truth_vector[names(v)]),estimate=unname(v),
      bias=unname(v-truth_vector[names(v)]))
    rates <- duration_weighted_transition_rates(df_fit,f)[1,]
    summary_records[[counter]] <- data.frame(n=n,replication=replication,
      start=label,selected=label==best_name,loglik=f$loglik,
      convergence=f$convergence,iterations=f$iterations,
      entry_rate=rates$entry_rate,exit_rate=rates$exit_rate,
      true_entry_rate=truth_rates$entry_rate,
      true_exit_rate=truth_rates$exit_rate,
      entry_rate_bias=rates$entry_rate-truth_rates$entry_rate,
      exit_rate_bias=rates$exit_rate-truth_rates$exit_rate,
      exact_anchor_retention=f$params$tenure_exact_anchor_retention_prob,
      tenure_contamination=f$params$eps,
      status_misclassification=f$params$pi,
      month_rmse=sqrt(mean((f$params$tenure_start_month_probs-
        truth$tenure_start_month_probs)^2)),
      exit_hazard_rmse=sqrt(mean((f$params$lambda_g-truth$lambda_g)^2)),
      entry_hazard_rmse=sqrt(mean((f$params$lambda_d-truth$lambda_d)^2)))
  }
}

parameters <- do.call(rbind,parameter_records)
fits <- do.call(rbind,summary_records)
selected_parameters <- parameters[parameters$selected,]
selected_fits <- fits[fits$selected,]
recovery_summary <- do.call(rbind,lapply(split(selected_parameters,
  selected_parameters$parameter),function(z) data.frame(
    parameter=z$parameter[1],truth=z$truth[1],mean_estimate=mean(z$estimate),
    mean_bias=mean(z$bias),rmse=sqrt(mean(z$bias^2)),
    minimum=min(z$estimate),maximum=max(z$estimate))))

write.csv(fits,file.path(outdir,"fit_summary_latest.csv"),row.names=FALSE)
write.csv(parameters,file.path(outdir,"parameter_draws_latest.csv"),
  row.names=FALSE)
write.csv(recovery_summary,file.path(outdir,"recovery_summary_latest.csv"),
  row.names=FALSE)
saveRDS(list(truth=truth,fits=fits,parameters=parameters,
  recovery_summary=recovery_summary),
  file.path(outdir,"validation_joint_latest.rds"))
cat("\nFull joint recovery fits\n")
print(fits,row.names=FALSE,digits=8)
cat("\nSelected-fit parameter recovery\n")
print(recovery_summary,row.names=FALSE,digits=8)
