# Simulation-recovery check for the exact reported-start-date retention
# probability. Other parameters are held at their generating values so this
# isolates whether the new channel is identified by the intended continuation
# pattern rather than by numerical movement elsewhere in the model.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")
fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate exact-anchor retention first")
truth <- readRDS(fit_file)$extension$params
n <- as.integer(Sys.getenv("EXACT_ANCHOR_RECOVERY_N","5000"))
reps <- as.integer(Sys.getenv("EXACT_ANCHOR_RECOVERY_REPS","20"))
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

recover_once <- function(replication) {
  d <- simulate_eps_piecewise_job_change(n,truth,
    seed=260917L+replication)
  df_fit <- collapse_eps_cells(d,allow_zero_tenure=TRUE,
    extra_cols=paste0("interview_month",1:3))
  objective <- function(rho) {
    p <- truth; p$tenure_exact_anchor_retention_prob <- rho
    -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
  }
  opt <- optimize(objective,c(1e-7,.95),tol=1e-6)
  data.frame(replication=replication,n=n,
    truth=truth$tenure_exact_anchor_retention_prob,
    estimate=opt$minimum,bias=opt$minimum-
      truth$tenure_exact_anchor_retention_prob,
    loglik=-opt$objective,loglik_zero=-objective(0),
    LR_zero=max(0,2*(-opt$objective+objective(0))))
}
recovery <- do.call(rbind,lapply(seq_len(reps),recover_once))
summary <- data.frame(reps=reps,n=n,truth=unique(recovery$truth),
  mean_estimate=mean(recovery$estimate),sd_estimate=sd(recovery$estimate),
  rmse=sqrt(mean((recovery$estimate-recovery$truth)^2)),
  minimum=min(recovery$estimate),maximum=max(recovery$estimate),
  mean_LR_zero=mean(recovery$LR_zero))
write.csv(recovery,file.path(outdir,"recovery_draws_latest.csv"),
  row.names=FALSE)
write.csv(summary,file.path(outdir,"recovery_summary_latest.csv"),
  row.names=FALSE)
saveRDS(list(truth=truth,recovery=recovery,summary=summary),
  file.path(outdir,"validation_latest.rds"))
cat("\nExact-anchor retention recovery\n")
print(summary,row.names=FALSE,digits=8)
cat("\nReplication estimates\n")
print(recovery,row.names=FALSE,digits=8)
