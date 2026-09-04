# Joint simulation-recovery check for the 11-parameter baseline start-month
# distribution and the existing transition/reporting parameters.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")
fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_start_month_baseline/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate the start-month model first")
truth <- readRDS(fit_file)$seasonal$params
n <- as.integer(Sys.getenv("START_MONTH_RECOVERY_N","10000"))
workers <- as.integer(Sys.getenv("START_MONTH_RECOVERY_WORKERS","4"))
maxit <- as.integer(Sys.getenv("START_MONTH_RECOVERY_MAXIT","50"))
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_start_month_baseline/validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

d <- simulate_eps_piecewise_job_change(n,truth,seed=260916L)
df_fit <- collapse_eps_cells(d,allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month",1:3))
starts <- list(uniform=rep(1/12,12L),reversed=
  rev(truth$tenure_start_month_probs))
pack_month <- function(probs) log(probs[1:11]/probs[12])
unpack_month <- function(z) {
  value <- c(z,0); value <- value-max(value)
  exp(value)/sum(exp(value))
}
month_objective <- function(z) {
  p <- truth; p$tenure_start_month_probs <- unpack_month(z)
  -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik/nrow(d)
}
month_fits <- lapply(names(starts),function(label) {
  opt <- optim(pack_month(starts[[label]]),month_objective,method="BFGS",
    control=list(maxit=min(maxit,40L),reltol=1e-8,trace=1L,REPORT=10L))
  list(label=label,opt=opt,probs=unpack_month(opt$par),
    loglik=-opt$value*nrow(d))
})
names(month_fits) <- names(starts)
best_month <- month_fits[[which.max(vapply(month_fits,`[[`,numeric(1),
  "loglik"))]]
saveRDS(month_fits,file.path(outdir,"conditional_month_fits_latest.rds"))
joint_start <- truth
joint_start$tenure_start_month_probs <- best_month$probs
forward <- fit_eps_piecewise_calendar_revision_monthly(df_fit,joint_start,
  heaping_start=max(joint_start$tenure_heaping_prob,.002),
  revision_start=joint_start$tenure_year_revision_prob,
  clean_anchor_revision_start=joint_start$tenure_clean_anchor_revision_prob,
  start_month_probs_start=joint_start$tenure_start_month_probs,
  q_start=joint_start$job_change_prob,maxit=maxit,reltol=1e-9,pgtol=1e-7,
  workers=workers,verbose=1L,gradient_step=1e-3,
  gradient_scheme="forward")
forward$objective_function <- NULL
saveRDS(forward,file.path(outdir,"fit_forward_latest.rds"))
best <- forward
central <- fit_eps_piecewise_calendar_revision_monthly(df_fit,best$params,
  heaping_start=max(best$params$tenure_heaping_prob,1e-5),
  revision_start=best$params$tenure_year_revision_prob,
  clean_anchor_revision_start=best$params$tenure_clean_anchor_revision_prob,
  start_month_probs_start=best$params$tenure_start_month_probs,
  q_start=best$params$job_change_prob,maxit=10L,reltol=1e-9,pgtol=1e-7,
  workers=workers,verbose=1L,gradient_step=1e-3,
  gradient_scheme="central")
central$objective_function <- NULL
if (central$loglik>=best$loglik) best <- central

fit_summary <- do.call(rbind,lapply(names(month_fits),function(label) {
  f <- month_fits[[label]]
  data.frame(stage="Conditional months",start=label,loglik=f$loglik,
    convergence=f$opt$convergence,
    iterations=unname(f$opt$counts["function"]),january=f$probs[1],
    month_rmse=sqrt(mean((f$probs-
      truth$tenure_start_month_probs)^2)),
    gross_january_heaping_prob=truth$tenure_heaping_prob,
    gross_anchor_revision_prob=truth$tenure_year_revision_prob,
    clean_anchor_revision_prob=truth$tenure_clean_anchor_revision_prob,
    tenure_contamination=truth$eps,status_misclassification=truth$pi)
}))
fit_summary <- rbind(fit_summary,data.frame(stage="Joint",
  start="forward",loglik=forward$loglik,convergence=forward$convergence,
  iterations=forward$iterations,
  january=forward$params$tenure_start_month_probs[1],
  month_rmse=sqrt(mean((forward$params$tenure_start_month_probs-
    truth$tenure_start_month_probs)^2)),
  gross_january_heaping_prob=forward$params$tenure_heaping_prob,
  gross_anchor_revision_prob=forward$params$tenure_year_revision_prob,
  clean_anchor_revision_prob=forward$params$tenure_clean_anchor_revision_prob,
  tenure_contamination=forward$params$eps,
  status_misclassification=forward$params$pi),data.frame(stage="Joint",
  start="central_refinement",
  loglik=central$loglik,convergence=central$convergence,
  iterations=central$iterations,
  january=central$params$tenure_start_month_probs[1],
  month_rmse=sqrt(mean((central$params$tenure_start_month_probs-
    truth$tenure_start_month_probs)^2)),
  gross_january_heaping_prob=central$params$tenure_heaping_prob,
  gross_anchor_revision_prob=central$params$tenure_year_revision_prob,
  clean_anchor_revision_prob=
    central$params$tenure_clean_anchor_revision_prob,
  tenure_contamination=central$params$eps,
  status_misclassification=central$params$pi))
month_recovery <- data.frame(month=month.abb,
  truth=truth$tenure_start_month_probs,
  estimate=best$params$tenure_start_month_probs)
month_recovery$bias <- month_recovery$estimate-month_recovery$truth
quantities <- c("pi","eps","eps_d","job_change_prob",
  "tenure_heaping_prob","tenure_year_revision_prob",
  "tenure_clean_anchor_revision_prob")
parameter_recovery <- data.frame(parameter=quantities,
  truth=unlist(truth[quantities]),estimate=unlist(best$params[quantities]))
parameter_recovery$bias <- parameter_recovery$estimate-
  parameter_recovery$truth

write.csv(fit_summary,file.path(outdir,"fits_latest.csv"),row.names=FALSE)
write.csv(month_recovery,file.path(outdir,"month_recovery_latest.csv"),
  row.names=FALSE)
write.csv(parameter_recovery,file.path(outdir,
  "parameter_recovery_latest.csv"),row.names=FALSE)
saveRDS(list(truth=truth,month_fits=month_fits,forward=forward,
  central=central,best=best,
  fit_summary=fit_summary,month_recovery=month_recovery,
  parameter_recovery=parameter_recovery),
  file.path(outdir,"validation_latest.rds"))
cat("\nStart-month recovery fits\n")
print(fit_summary,row.names=FALSE,digits=8)
cat("\nMonth distribution recovery\n")
print(month_recovery,row.names=FALSE,digits=8)
cat("\nOther parameter recovery\n")
print(parameter_recovery,row.names=FALSE,digits=8)
