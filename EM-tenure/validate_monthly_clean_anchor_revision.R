# Recovery diagnostic for the clean-anchor whole-year revision probability.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")
base_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_revision/fits_latest.rds")
if (!file.exists(base_file)) stop("Estimate the calendar-revision model first")
truth <- readRDS(base_file)$revision$params
truth$tenure_clean_anchor_revision_prob <- as.numeric(Sys.getenv(
  "CLEAN_ANCHOR_RECOVERY_TRUTH",".35"))
n <- as.integer(Sys.getenv("CLEAN_ANCHOR_RECOVERY_N","8000"))
workers <- as.integer(Sys.getenv("CLEAN_ANCHOR_RECOVERY_WORKERS","1"))
maxit <- as.integer(Sys.getenv("CLEAN_ANCHOR_RECOVERY_MAXIT","70"))
seed <- 260912L
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_clean_anchor_revision/validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

d <- simulate_eps_piecewise_job_change(n,truth,seed=seed)
df_fit <- collapse_eps_cells(d,allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month",1:3))

grid <- seq(.02,.70,length.out=29L)
conditional_loglik <- vapply(grid,function(eta) {
  p <- truth; p$tenure_clean_anchor_revision_prob <- eta
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
},numeric(1))
conditional <- data.frame(clean_anchor_revision_prob=grid,
  loglik=conditional_loglik,
  relative_loglik=conditional_loglik-max(conditional_loglik))
conditional_optimum <- optimize(function(eta) {
  p <- truth; p$tenure_clean_anchor_revision_prob <- eta
  -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
},c(.001,.90),tol=1e-5)

fit_one <- function(label,eta_start) {
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,truth,
    heaping_start=truth$tenure_heaping_prob,
    revision_start=truth$tenure_year_revision_prob,
    clean_anchor_revision_start=eta_start,q_start=truth$job_change_prob,
    maxit=maxit,reltol=1e-9,pgtol=1e-7,workers=workers,verbose=1L)
  fit$objective_function <- NULL
  saveRDS(fit,file.path(outdir,paste0("fit_",label,"_latest.rds")))
  fit
}
fits <- list(low=fit_one("low",.08),high=fit_one("high",.65))
fit_summary <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,clean_anchor_revision_prob=
      p$tenure_clean_anchor_revision_prob,gross_anchor_revision_prob=
      p$tenure_year_revision_prob,january_heaping_prob=
      p$tenure_heaping_prob,tenure_contamination=p$eps,
      job_change_prob=p$job_change_prob,status_misclassification=p$pi)
}))
best <- fits[[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]]
quantities <- c("tenure_clean_anchor_revision_prob",
  "tenure_year_revision_prob","tenure_heaping_prob","eps",
  "job_change_prob","pi")
recovery <- data.frame(parameter=quantities,
  truth=unlist(truth[quantities]),estimate=unlist(best$params[quantities]))
recovery$bias <- recovery$estimate-recovery$truth

cat("\nConditional clean-anchor revision profile\n")
print(conditional,row.names=FALSE,digits=7)
cat(sprintf("\nConditional optimum %.5f (truth %.5f)\n",
  conditional_optimum$minimum,truth$tenure_clean_anchor_revision_prob))
cat("\nJoint recovery starts\n"); print(fit_summary,row.names=FALSE,digits=8)
cat("\nPreferred recovery fit\n"); print(recovery,row.names=FALSE,digits=8)

write.csv(conditional,file.path(outdir,"conditional_profile_latest.csv"),
  row.names=FALSE)
write.csv(fit_summary,file.path(outdir,"fits_latest.csv"),row.names=FALSE)
write.csv(recovery,file.path(outdir,"recovery_latest.csv"),row.names=FALSE)
saveRDS(list(truth=truth,data=d,conditional=conditional,
  conditional_optimum=conditional_optimum,fits=fits,recovery=recovery),
  file.path(outdir,"validation_latest.rds"))
