# Recovery diagnostic for the persistent reported-start-date extension. Data
# are simulated from a positive persistence parameter and re-estimated using
# the same full observed likelihood. Outputs are resumable and separate from
# empirical fits.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
source("EM-tenure/R/source_all.R")

base_file <- "EM-tenure/output/results/job_change_monthly/fits_latest.rds"
if (!file.exists(base_file)) stop("Estimate the corrected Panel A model first")
truth <- readRDS(base_file)$reset$params
truth$tenure_measurement_model <- "monthly"
truth$tenure_report_persistence <- as.numeric(Sys.getenv(
  "PERSISTENT_VALIDATION_RHO","0.60"))

n <- as.integer(Sys.getenv("PERSISTENT_VALIDATION_N","12000"))
workers <- as.integer(Sys.getenv("PERSISTENT_VALIDATION_WORKERS","1"))
maxit <- as.integer(Sys.getenv("PERSISTENT_VALIDATION_MAXIT","60"))
resume <- identical(tolower(Sys.getenv("PERSISTENT_VALIDATION_RESUME","true")),
  "true")
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "persistent_reporting/validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

sim_file <- file.path(outdir,"simulated_data_latest.rds")
if (resume && file.exists(sim_file)) sim <- readRDS(sim_file) else {
  sim <- simulate_eps_piecewise_job_change(n,truth,seed=260831L)
  saveRDS(sim,sim_file)
}
validate_df_eps(sim,allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(sim,allow_zero_tenure=TRUE)

evaluate_rho <- function(rho) {
  p <- truth; p$tenure_report_persistence <- rho
  e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
}
rho_grid <- seq(0,.95,by=.05)
conditional <- data.frame(rho=rho_grid,
  loglik=vapply(rho_grid,evaluate_rho,numeric(1)))
conditional$relative_loglik <- conditional$loglik-max(conditional$loglik)
conditional_optimum <- optimize(function(rho) -evaluate_rho(rho),
  interval=c(1e-6,.999),tol=1e-5)

fit_one <- function(label,rho_start) {
  path <- file.path(outdir,paste0("fit_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  start <- truth
  start$eps <- min(.45,max(.05,truth$eps + if (rho_start<.5) -.04 else .04))
  start$job_change_prob <- max(.003,truth$job_change_prob*
    if (rho_start<.5) .7 else 1.4)
  fit <- fit_eps_piecewise_persistent_monthly(df_fit,start,
    persistence_start=rho_start,q_start=start$job_change_prob,maxit=maxit,
    reltol=1e-9,pgtol=1e-7,workers=workers,verbose=1L)
  fit$objective_function <- NULL
  saveRDS(fit,path)
  fit
}
fits <- list(low=fit_one("low",.20),high=fit_one("high",.85))
fit_summary <- do.call(rbind,lapply(names(fits),function(label) {
  f <- fits[[label]]; p <- f$params
  data.frame(start=label,loglik=f$loglik,convergence=f$convergence,
    iterations=f$iterations,tenure_report_persistence=
      p$tenure_report_persistence,job_change_prob=p$job_change_prob,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d)
}))
best <- fits[[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]]
quantities <- c("tenure_report_persistence","job_change_prob","pi","eps",
  "eps_d")
recovery <- data.frame(parameter=quantities,
  truth=unlist(truth[quantities]),estimate=unlist(best$params[quantities]))
recovery$bias <- recovery$estimate-recovery$truth

cat("\nConditional persistence profile at generating nuisance parameters\n")
print(conditional,row.names=FALSE,digits=8)
cat(sprintf("\nConditional optimum: rho = %.6f (truth %.6f)\n",
  conditional_optimum$minimum,truth$tenure_report_persistence))
cat("\nJoint full-likelihood starts\n")
print(fit_summary,row.names=FALSE,digits=8)
cat("\nBest-fit recovery\n"); print(recovery,row.names=FALSE,digits=8)

write.csv(conditional,file.path(outdir,"conditional_profile_latest.csv"),
  row.names=FALSE)
write.csv(fit_summary,file.path(outdir,"fits_latest.csv"),row.names=FALSE)
write.csv(recovery,file.path(outdir,"recovery_latest.csv"),row.names=FALSE)
saveRDS(list(truth=truth,conditional=conditional,
  conditional_optimum=conditional_optimum,fit_summary=fit_summary,
  recovery=recovery,fits=fits),file.path(outdir,"validation_latest.rds"))
