# Simulation-recovery checks for the model with separate tenure and timegap
# reliability dispersions. The primary check estimates the four parameters
# that distinguish the two reliability components while holding the remaining
# nuisance parameters at their generating values. Multiple starts and cached
# replication-level results make the check resumable.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")

fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate separate duration reliability first")
truth <- readRDS(fit_file)$extension$params
truth$tenure_measurement_model <- "monthly"

n <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_RECOVERY_N","20000"))
reps <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_RECOVERY_REPS","3"))
maxit <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_RECOVERY_MAXIT","40"))
workers <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_RECOVERY_WORKERS","8"))
resume <- identical(tolower(Sys.getenv(
  "SEPARATE_RELIABILITY_RECOVERY_RESUME","true")),"true")
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

free_names <- c("eps","eps_d","tenure_reliability_shift",
  "timegap_reliability_shift")
quantity_names <- c("eps","eps_d","tenure_reliability_shift",
  "timegap_reliability_shift")

fit_once <- function(replication) {
  path <- file.path(outdir,sprintf("replication_%02d_latest.rds",replication))
  if (resume && file.exists(path)) {
    cached <- readRDS(path)
    if (identical(cached$n,n) && identical(cached$maxit,maxit)) return(cached)
  }
  d <- simulate_eps_piecewise_job_change(n,truth,
    seed=260905L+replication)
  class <- d$duration_reliability_class
  df_fit <- collapse_eps_cells(d,allow_zero_tenure=TRUE,
    extra_cols=paste0("interview_month",1:3))
  starts <- list(
    truth=c(eps=unname(truth$eps),eps_d=unname(truth$eps_d),
      tenure=unname(truth$tenure_reliability_shift),
      timegap=unname(truth$timegap_reliability_shift)),
    compressed=c(eps=.10,eps_d=.10,tenure=1,timegap=1),
    dispersed=c(eps=.20,eps_d=.08,tenure=3,timegap=3))
  fits <- lapply(names(starts),function(label) {
    message("Recovery replication ",replication,", start ",label)
    s <- starts[[label]]
    start <- truth; start$eps <- unname(s["eps"])
    start$eps_d <- unname(s["eps_d"])
    start$tenure_reliability_shift <- unname(s["tenure"])
    start$timegap_reliability_shift <- unname(s["timegap"])
    fit <- tryCatch(fit_eps_piecewise_calendar_revision_monthly(df_fit,start,
        heaping_start=start$tenure_heaping_prob,
        revision_start=start$tenure_year_revision_prob,
        clean_anchor_revision_start=start$tenure_clean_anchor_revision_prob,
        start_month_probs_start=start$tenure_start_month_probs,
        exact_anchor_retention_start=
          start$tenure_exact_anchor_retention_prob,
        local_revision_start=start$tenure_local_revision_prob,
        tenure_reliability_shift_start=
          start$tenure_reliability_shift,
        timegap_reliability_shift_start=
          start$timegap_reliability_shift,
        q_start=start$job_change_prob,maxit=maxit,reltol=1e-10,
        pgtol=1e-7,workers=max(1L,min(workers,8L)),verbose=0L,
        gradient_step=5e-4,free_names=free_names,
        gradient_scheme="central"),error=function(error)
          structure(list(start_label=label,error=conditionMessage(error)),
            class="failed_recovery_fit"))
    if (inherits(fit,"failed_recovery_fit")) return(fit)
    fit$objective_function <- NULL
    fit$start_label <- label
    fit
  })
  successful <- !vapply(fits,inherits,logical(1),"failed_recovery_fit")
  if (!any(successful)) stop("All recovery starts failed in replication ",
    replication,": ",paste(vapply(fits,`[[`,character(1),"error"),
      collapse="; "))
  successful_fits <- fits[successful]
  best <- successful_fits[[which.max(vapply(successful_fits,`[[`,numeric(1),
    "loglik"))]]
  # The simulator preserves the generating class label. Compare it with the
  # fitted posterior before collapsing, because class membership is individual.
  full_eval <- e_step_eps(d,best$params,check_df=FALSE,suff_stats=FALSE)
  posterior <- full_eval$duration_reliability_posterior
  unreliable <- as.integer(class=="unreliable")
  ranks <- rank(posterior,ties.method="average")
  n1 <- sum(unreliable); n0 <- length(unreliable)-n1
  auc <- (sum(ranks[unreliable==1])-n1*(n1+1)/2)/(n1*n0)
  result <- list(replication=replication,n=n,maxit=maxit,best=best,fits=fits,
    class_diagnostic=data.frame(
      posterior_mean_reliable=mean(posterior[unreliable==0]),
      posterior_mean_unreliable=mean(posterior[unreliable==1]),AUC=auc))
  saveRDS(result,path)
  result
}

results <- lapply(seq_len(reps),fit_once)
draws <- do.call(rbind,lapply(results,function(result) {
  p <- result$best$params
  data.frame(replication=result$replication,loglik=result$best$loglik,
    convergence=result$best$convergence,iterations=result$best$iterations,
    eps=p$eps,eps_d=p$eps_d,
    tenure_reliability_shift=p$tenure_reliability_shift,
    timegap_reliability_shift=p$timegap_reliability_shift,
    result$class_diagnostic)
}))
truth_values <- unlist(truth[quantity_names])
summary <- do.call(rbind,lapply(quantity_names,function(quantity) {
  estimate <- draws[[quantity]]
  data.frame(parameter=quantity,truth=truth_values[[quantity]],
    mean_estimate=mean(estimate),bias=mean(estimate)-truth_values[[quantity]],
    rmse=sqrt(mean((estimate-truth_values[[quantity]])^2)),
    minimum=min(estimate),maximum=max(estimate))
}))
start_summary <- do.call(rbind,lapply(results,function(result)
  do.call(rbind,lapply(result$fits,function(fit) {
    if (inherits(fit,"failed_recovery_fit")) return(data.frame(
      replication=result$replication,start=fit$start_label,
      loglik=NA_real_,convergence=NA_integer_,iterations=NA_integer_,
      error=fit$error))
    data.frame(replication=result$replication,start=fit$start_label,
      loglik=fit$loglik,convergence=fit$convergence,
      iterations=fit$iterations,error=NA_character_)
  }))))

write.csv(draws,file.path(outdir,"recovery_draws_latest.csv"),row.names=FALSE)
write.csv(summary,file.path(outdir,"recovery_summary_latest.csv"),
  row.names=FALSE)
write.csv(start_summary,file.path(outdir,"multistart_summary_latest.csv"),
  row.names=FALSE)
saveRDS(list(truth=truth,draws=draws,summary=summary,
  start_summary=start_summary,results=results),
  file.path(outdir,"validation_latest.rds"))
cat("\nSeparate-duration-reliability conditional recovery draws\n")
print(draws,row.names=FALSE,digits=8)
cat("\nRecovery summary\n"); print(summary,row.names=FALSE,digits=8)
cat("\nMultiple-start comparison\n")
print(start_summary,row.names=FALSE,digits=8)
