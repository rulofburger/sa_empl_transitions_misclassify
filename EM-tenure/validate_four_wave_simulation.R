# Gate simulation recovery on agreement between generator and likelihood.
# This diagnoses the sequential simulator; it does not refit empirical data.
source("EM-tenure/R/source_all.R")
base <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
clock_model <- match.arg(Sys.getenv("FOUR_WAVE_TIMEGAP_CLOCK","continuous_joint"),
  c("continuous_joint","legacy_pairwise"))
outdir <- file.path(base,if (clock_model=="continuous_joint")
  "recovery_continuous_clock" else "recovery")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
truth <- readRDS(file.path(base,"converged_comparison_latest.rds"))$best$params
n <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_SCORE_N","20000"))
reps <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_REPS","3"))
stopifnot(n>=1000L,reps>=1L)
z <- .pack_four_wave_preferred(truth)
records <- vector("list",reps)
for (r in seq_len(reps)) {
  seed <- 260903L+r
  message("Sequential-generator score check: replication ",r,"; n=",n)
  d <- simulate_eps_piecewise_job_change(n,truth,seed=seed,waves=4,exact_risk=TRUE)
  data <- prepare_four_wave_kernel_data(d)
  eval_z <- function(v) evaluate_four_wave_monthly_fast(data,
    .piecewise_calendar_revision_monthly_unpack(v,"joint_marginal"),
    posterior=FALSE,threads=8L,timegap_clock=clock_model)$row_loglik
  scores <- matrix(NA_real_,n,length(z),dimnames=list(NULL,names(z)))
  for (j in seq_along(z)) {
    step <- 1e-5*max(1,abs(z[j]))
    plus <- minus <- z
    plus[j] <- plus[j]+step
    minus[j] <- minus[j]-step
    scores[,j] <- (eval_z(plus)-eval_z(minus))/(2*step)
  }
  stopifnot(all(is.finite(scores)))
  mean_score <- colMeans(scores)
  mcse <- apply(scores,2,sd)/sqrt(n)
  records[[r]] <- data.frame(replication=r,n=n,seed=seed,parameter=names(z),
    mean_score=mean_score,mcse=mcse,z_score=mean_score/mcse)
  saveRDS(list(data=d,truth=truth,seed=seed,score=scores,
    summary=records[[r]],source_md5=four_wave_fast_source_fingerprint()),
    file.path(outdir,sprintf("sequential_score_n%d_rep%d.rds",n,r)))
  print(records[[r]][order(-abs(records[[r]]$z_score)),][1:10,],row.names=FALSE)
}
draws <- do.call(rbind,records)
pooled <- do.call(rbind,lapply(split(draws,draws$parameter),function(x)
  data.frame(parameter=x$parameter[1],n_total=sum(x$n),
    mean_score=mean(x$mean_score),mcse=sqrt(sum(x$mcse^2))/nrow(x))))
pooled$z_score <- pooled$mean_score/pooled$mcse
# Deliberately conservative screening threshold across 33 coordinates.
passed <- all(abs(pooled$z_score)<5)
write.csv(draws,file.path(outdir,"sequential_score_draws_latest.csv"),row.names=FALSE)
write.csv(pooled,file.path(outdir,"sequential_score_summary_latest.csv"),row.names=FALSE)
saveRDS(list(passed=passed,draws=draws,summary=pooled,truth=truth,timegap_clock=clock_model,
  source_md5=four_wave_fast_source_fingerprint(),
  interpretation="Generator-likelihood agreement gate, not parameter recovery"),
  file.path(outdir,"sequential_score_gate_latest.rds"))
print(pooled[order(-abs(pooled$z_score)),],row.names=FALSE)
message("Generator-likelihood score gate passed: ",passed)
if (!passed) warning("Do not treat fits to this generator as validated recovery")
