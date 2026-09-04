# Add a saturated, shared distribution for the implied calendar month in which
# a job began.  Eleven month logits are estimated relative to December.  The
# clean-anchor calendar-revision model is nested when all 12 probabilities are
# equal.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
interview_cols <- paste0("interview_month",1:3)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=interview_cols)
screen_n <- min(as.integer(Sys.getenv("START_MONTH_SCREEN_N","50000")),
  nrow(df_full))
set.seed(260914L)
screen_index <- if (screen_n<nrow(df_full))
  sample.int(nrow(df_full),screen_n) else seq_len(nrow(df_full))
df_screen <- collapse_eps_cells(df_full[screen_index,],
  allow_zero_tenure=TRUE,extra_cols=interview_cols)

base_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_clean_anchor_revision/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the clean-anchor model first")
base <- readRDS(base_file)$extension
base$params$tenure_start_month_probs <- rep(1/12,12L)
base_eval <- e_step_eps(df_fit,base$params,check_df=FALSE,suff_stats=FALSE)
base$loglik <- base_eval$loglik
base$gamma <- base_eval$gamma
base$job_change_posterior <- base_eval$job_change_posterior

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_start_month_baseline")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("START_MONTH_BASELINE_WORKERS","1"))
maxit <- as.integer(Sys.getenv("START_MONTH_BASELINE_MAXIT","120"))
resume <- identical(tolower(Sys.getenv("START_MONTH_BASELINE_RESUME","true")),
  "true")

nested_loglik <- e_step_eps(df_fit,base$params,check_df=FALSE,
  suff_stats=FALSE)$loglik
nesting <- data.frame(base_loglik=base$loglik,
  uniform_month_loglik=nested_loglik,difference=nested_loglik-base$loglik)
if (abs(nesting$difference)>1e-6)
  stop("uniform start months do not reproduce the clean-anchor model")

# Use the empirical distribution only to disperse starts.  Shrink it toward
# uniform because its January excess also contains gross-report heaping.
month_count <- numeric(12L)
for (wave in 1:3) {
  keep <- df_full[[paste0("y",wave)]]==1L
  tenure_month <- round(12*df_full[[paste0("tenure",wave)]][keep])
  start_month <- ((df_full[[paste0("interview_month",wave)]][keep]-1L-
    tenure_month%%12L)%%12L)+1L
  month_count <- month_count+vapply(1:12,function(m)
    sum(df_full$weight[keep][start_month==m]),numeric(1))
}
empirical <- month_count/sum(month_count)
smoothed <- .5*empirical+.5/12
january_start <- c(.20,rep(.80/11,11))

starts <- list(uniform=rep(1/12,12L),smoothed_empirical=smoothed,
  january_heavy=january_start)

# First isolate the 11 new month effects while holding the already-estimated
# parameters fixed. This gives a stable, inexpensive mode screen before the
# full 29-parameter refinement.
month_pack <- function(probs) log(probs[1:11]/probs[12])
month_unpack <- function(z) {
  value <- c(unname(z),0); value <- value-max(value)
  exp(value)/sum(exp(value))
}
total_weight <- sum(df_fit$weight)
screen_weight <- sum(df_screen$weight)
conditional_workers <- max(1L,min(workers,4L))
joint_workers <- max(1L,min(workers,8L))
cluster <- if (conditional_workers>1L)
  parallel::makePSOCKcluster(conditional_workers) else NULL
if (!is.null(cluster)) {
  worker_wd <- getwd()
  parallel::clusterCall(cluster,function(path) {
    setwd(path); source("EM-tenure/R/source_all.R"); NULL
  },worker_wd)
  df_worker <- df_screen; base_worker <- base$params
  total_weight_worker <- screen_weight
  parallel::clusterExport(cluster,c("df_worker","base_worker",
    "total_weight_worker"),envir=environment())
}
conditional_objective <- function(z) {
  p <- base$params; p$tenure_start_month_probs <- month_unpack(z)
  -e_step_eps(df_screen,p,check_df=FALSE,suff_stats=FALSE)$loglik/screen_weight
}
conditional_worker <- function(z) {
  value <- c(unname(z),0); value <- value-max(value)
  p <- base_worker; p$tenure_start_month_probs <- exp(value)/sum(exp(value))
  -e_step_eps(df_worker,p,check_df=FALSE,suff_stats=FALSE)$loglik/
    total_weight_worker
}
conditional_gradient <- function(z) {
  step <- 1e-3*pmax(1,abs(z))
  points <- lapply(seq_along(z),function(j) {
    value <- z; value[j] <- value[j]+step[j]; value
  })
  plus <- if (is.null(cluster)) vapply(points,conditional_objective,numeric(1))
    else unlist(parallel::parLapplyLB(cluster,points,conditional_worker),
      use.names=FALSE)
  out <- (plus-conditional_objective(z))/step
  out[!is.finite(out)] <- 0
  out
}
conditional_fits <- lapply(names(starts),function(label) {
  message("Screening conditional month effects from ",label," start")
  z0 <- month_pack(starts[[label]])
  opt <- optim(z0,conditional_objective,gr=conditional_gradient,
    method="L-BFGS-B",lower=rep(-8,11),upper=rep(8,11),
    control=list(maxit=maxit,factr=1e7,pgtol=1e-7,trace=1L))
  list(label=label,opt=opt,probs=month_unpack(opt$par),
    screen_loglik=-opt$value*screen_weight)
})
names(conditional_fits) <- names(starts)
if (!is.null(cluster)) parallel::stopCluster(cluster)
for (label in names(conditional_fits)) {
  p <- base$params
  p$tenure_start_month_probs <- conditional_fits[[label]]$probs
  conditional_fits[[label]]$loglik <- e_step_eps(df_fit,p,check_df=FALSE,
    suff_stats=FALSE)$loglik
}
best_conditional <- conditional_fits[[which.max(vapply(conditional_fits,
  `[[`,numeric(1),"loglik"))]]
conditional_params <- base$params
conditional_params$tenure_start_month_probs <- best_conditional$probs

fit_joint <- function(label,data,start_params,gradient_scheme,iterations) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Joint refinement: ",label)
  fit <- fit_eps_piecewise_calendar_revision_monthly(data,start_params,
    heaping_start=start_params$tenure_heaping_prob,
    revision_start=start_params$tenure_year_revision_prob,
    clean_anchor_revision_start=
      start_params$tenure_clean_anchor_revision_prob,
    start_month_probs_start=start_params$tenure_start_month_probs,
    q_start=start_params$job_change_prob,maxit=iterations,reltol=1e-9,
    pgtol=1e-7,workers=joint_workers,verbose=1L,
    gradient_step=1e-3,gradient_scheme=gradient_scheme)
  fit$objective_function <- NULL
  saveRDS(fit,path)
  fit
}
fit_screen_forward <- fit_joint("screen_forward",df_screen,
  conditional_params,"forward",min(maxit,35L))
fit_screen_central <- fit_joint("screen_central",df_screen,
  fit_screen_forward$params,"central",min(maxit,10L))
fit_forward <- fit_joint("full_forward",df_fit,fit_screen_central$params,
  "forward",min(maxit,8L))
fit_central <- fit_joint("full_central",df_fit,fit_forward$params,"central",
  min(maxit,3L))
fits <- list(conditional=structure(list(params=conditional_params,
    loglik=best_conditional$loglik,convergence=best_conditional$opt$convergence,
    iterations=unname(best_conditional$opt$counts["function"])),
    class="conditional_month_fit"),screen_forward=fit_screen_forward,
    screen_central=fit_screen_central,forward=fit_forward,central=fit_central)
joint_fits <- fits[c("forward","central")]
best <- joint_fits[[which.max(vapply(joint_fits,`[[`,numeric(1),"loglik"))]]

summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior
  posterior_q <- sum(df_fit$weight*jp$expected_changes)/
    sum(df_fit$weight*jp$opportunities)
  p <- fit$params
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=posterior_q,
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
comparison <- rbind(
  summarise_fit("Uniform baseline start months",base,18L),
  summarise_fit("Flexible baseline start months",best,29L))
start_comparison <- rbind(do.call(rbind,lapply(names(conditional_fits),
  function(label) {
    f <- conditional_fits[[label]]
    data.frame(stage="Conditional month effects",start=label,
      screen_loglik=f$screen_loglik,loglik=f$loglik,
      convergence=f$opt$convergence,
      iterations=unname(f$opt$counts["function"]),january=f$probs[1],
      gross_january_heaping_prob=base$params$tenure_heaping_prob,
      gross_anchor_revision_prob=base$params$tenure_year_revision_prob,
      clean_anchor_revision_prob=
        base$params$tenure_clean_anchor_revision_prob,
      tenure_contamination=base$params$eps,
      status_misclassification=base$params$pi)
  })),do.call(rbind,lapply(c("forward","central"),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(stage="Full-sample joint refinement",start=label,
    screen_loglik=NA_real_,loglik=f$loglik,
    convergence=f$convergence,iterations=f$iterations,
    january=p$tenure_start_month_probs[1],
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi)
})))
month_distribution <- data.frame(month=month.abb,
  observed_report_share=empirical,estimated_baseline=
    best$params$tenure_start_month_probs,uniform=1/12)
LR <- max(0,2*(best$loglik-base$loglik))
lr <- data.frame(comparison="Equal baseline start-month probabilities",
  LR=LR,df=11L,p_value=pchisq(LR,11,lower.tail=FALSE))

cat("\nExact nesting check\n"); print(nesting,row.names=FALSE,digits=12)
cat("\nFlexible start-month comparison\n")
print(comparison,row.names=FALSE,digits=8)
cat("\nStarting-value comparison\n")
print(start_comparison,row.names=FALSE,digits=8)
cat("\nBaseline start-month distribution\n")
print(month_distribution,row.names=FALSE,digits=8)
cat("\nLikelihood-ratio test\n"); print(lr,row.names=FALSE,digits=8)

write.csv(nesting,file.path(outdir,"nesting_latest.csv"),row.names=FALSE)
write.csv(comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(start_comparison,file.path(outdir,"starting_values_latest.csv"),
  row.names=FALSE)
write.csv(month_distribution,file.path(outdir,
  "start_month_distribution_latest.csv"),row.names=FALSE)
write.csv(lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
saveRDS(list(base=base,seasonal=best,fits=fits,nesting=nesting,
  comparison=comparison,start_comparison=start_comparison,
  month_distribution=month_distribution,lr=lr),
  file.path(outdir,"fits_latest.rds"))
