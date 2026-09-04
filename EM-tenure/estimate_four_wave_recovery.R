# Small full-nuisance recovery exercise for the repaired four-wave likelihood.
# All 33 parameters are free. No empirical panel is re-estimated.
source("EM-tenure/R/source_all.R")
base <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
outdir <- file.path(base,"recovery_continuous_clock")
gate <- readRDS(file.path(outdir,"recovery_status_latest.rds"))
validation <- readRDS(file.path(base,"compiled_equivalence_latest.rds"))
fingerprint <- four_wave_fast_source_fingerprint()
stopifnot(isTRUE(gate$ready_for_recovery_fits),
  identical(gate$source_md5,fingerprint),
  identical(validation$source_md5,fingerprint),
  identical(validation$timegap_clock,"continuous_joint"))
empirical <- readRDS(file.path(base,"converged_comparison_latest.rds"))$best
truth <- empirical$params
truth$timegap_clock_model <- "continuous_joint"
ztruth <- .pack_four_wave_preferred(truth)
bounds <- .four_wave_parameter_bounds(ztruth)
parscale <- empirical$controls$parscale
if (is.null(names(parscale))) {
  stopifnot(length(parscale)==length(ztruth),
    identical(names(empirical$par_unconstrained),names(ztruth)))
  names(parscale) <- names(ztruth)
} else parscale <- parscale[names(ztruth)]
stopifnot(identical(names(parscale),names(ztruth)),all(is.finite(parscale)))
n <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_N","5000"))
reps <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_FIT_REPS","2"))
maxit <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_MAXIT","180"))
workers <- as.integer(Sys.getenv("FOUR_WAVE_RECOVERY_WORKERS","4"))
stopifnot(n>=1000,reps>=1,workers>=1,workers<=8)
score_tolerance <- 1e-5
cluster <- parallel::makePSOCKcluster(workers)
tryCatch(for (node in seq_along(cluster)) parallel::clusterCall(cluster[node],
  function(path) {setwd(path);source("EM-tenure/R/source_all.R");
    load_four_wave_monthly_kernel();NULL}, getwd()),
  error=function(e) {parallel::stopCluster(cluster);stop(e)})
worker_objective <- function(z) {
  p <- .piecewise_calendar_revision_monthly_unpack(z,"joint_marginal")
  -evaluate_four_wave_monthly_fast(data_worker,p,posterior=FALSE,threads=1L,
    timegap_clock="continuous_joint")$loglik/weight_worker
}
parallel::clusterExport(cluster,"worker_objective")
results <- list(); details <- list()
tryCatch({
  for (replication in seq_len(reps)) {
    seed <- 261003L+replication
    d <- simulate_eps_piecewise_job_change(n,truth,seed=seed,waves=4,exact_risk=TRUE)
    df <- collapse_eps_cells_4w(d,extra_cols=paste0("interview_month",1:4))
    data_worker <- prepare_four_wave_kernel_data(df)
    weight_worker <- sum(df$weight)
    parallel::clusterExport(cluster,c("data_worker","weight_worker"))
    truth_eval <- evaluate_four_wave_monthly_fast(data_worker,truth,threads=4L)
    target <- duration_weighted_transition_rates_4w(df,list(params=truth,
      gamma=truth_eval$gamma,integration_method="exact_piecewise"))[1,]
    starts <- list(near=ztruth+.12*sin(seq_along(ztruth)))
    if (replication==1L) starts$alternative <- ztruth-.12*cos(seq_along(ztruth))
    for (label in names(starts)) {
      id <- sprintf("n%d_rep%d_%s",n,replication,label)
      path <- file.path(outdir,paste0("fit_",id,".rds"))
      checkpoint <- file.path(outdir,paste0("checkpoint_",id,".rds"))
      cached <- if (file.exists(path)) readRDS(path) else NULL
      if (!is.null(cached) && identical(cached$source_md5,fingerprint) &&
          isTRUE(cached$converged)) {
        fit <- cached
      } else {
        z0 <- pmin(bounds$upper,pmax(bounds$lower,starts[[label]]))
        if (file.exists(checkpoint)) {
          old <- readRDS(checkpoint)
          if (identical(old$source_md5,fingerprint)) z0 <- old$par
        }
        best <- Inf; gradients <- 0L; history <- data.frame()
        objective <- function(z) {
          value <- worker_objective(z)
          if (!is.finite(value)) stop("Nonfinite recovery objective")
          if (value<best) {
            best <<- value
            history <<- rbind(history,data.frame(gradient=gradients,
              loglik=-value*weight_worker))
            saveRDS(list(par=z,loglik=-value*weight_worker,history=history,
              source_md5=fingerprint),checkpoint)
          }
          value
        }
        gradient <- function(z) {
          step <- 1e-5*pmax(1,abs(z))
          plus <- pmin(step,bounds$upper-z)
          minus <- pmin(step,z-bounds$lower)
          points <- vector("list",2*length(z))
          for (j in seq_along(z)) {
            points[[2*j-1]] <- points[[2*j]] <- z
            points[[2*j-1]][j] <- z[j]+plus[j]
            points[[2*j]][j] <- z[j]-minus[j]
          }
          values <- unlist(parallel::parLapplyLB(cluster,points,worker_objective,
            chunk.size=1L),use.names=FALSE)
          g <- (values[seq(1,length(values),2)]-values[seq(2,length(values),2)])/
            (plus+minus)
          if (any(!is.finite(g))) stop("Nonfinite recovery score")
          gradients <<- gradients+1L
          if (gradients%%5==0) message(id,": gradient ",gradients,
            ", max score ",signif(max(abs(g)),4),", best LL ",signif(-best*weight_worker,10))
          g
        }
        message("Full 33-parameter recovery fit: ",id)
        opt <- optim(z0,objective,gradient,method="L-BFGS-B",
          lower=bounds$lower,upper=bounds$upper,control=list(maxit=maxit,
            factr=1e5,pgtol=score_tolerance*min(parscale),lmm=33,parscale=parscale))
        score <- gradient(opt$par)
        projected <- score
        projected[opt$par<=bounds$lower+1e-8 & score>0] <- 0
        projected[opt$par>=bounds$upper-1e-8 & score<0] <- 0
        p <- .piecewise_calendar_revision_monthly_unpack(opt$par,"joint_marginal")
        p$timegap_clock_model <- "continuous_joint"
        evaluated <- evaluate_four_wave_monthly_fast(data_worker,p,threads=4L)
        fit <- c(list(params=p,par_unconstrained=opt$par,
          integration_method="exact_piecewise",timegap_clock="continuous_joint",
          convergence=opt$convergence,message=opt$message,
          converged=opt$convergence==0L && max(abs(projected))<=score_tolerance,
          projected_score=projected,central_score=score,history=history,
          n=n,replication=replication,seed=seed,start=label,
          truth=truth,target_rates=target,source_md5=fingerprint,
          controls=list(maxit=maxit,score_tolerance=score_tolerance,
            parscale=parscale,workers=workers),session_info=sessionInfo()),evaluated)
        saveRDS(fit,path)
      }
      rates <- duration_weighted_transition_rates_4w(df,fit)[1,]
      summary <- data.frame(replication=replication,start=label,n=n,
        converged=fit$converged,score=max(abs(fit$projected_score)),loglik=fit$loglik,
        entry_truth=target$entry_rate,entry_estimate=rates$entry_rate,
        exit_truth=target$exit_rate,exit_estimate=rates$exit_rate,
        pi_truth=truth$pi,pi_estimate=fit$params$pi,
        alpha_truth=truth$alpha,alpha_estimate=fit$params$alpha,
        exit_hazard_rmse=sqrt(mean((fit$params$lambda_g-truth$lambda_g)^2)),
        entry_hazard_rmse=sqrt(mean((fit$params$lambda_d-truth$lambda_d)^2)))
      results[[id]] <- summary
      coefficient <- function(p) c(pi=p$pi,alpha=p$alpha,
        setNames(p$lambda_g,paste0("exit_hazard_",1:5)),
        setNames(p$lambda_d,paste0("entry_hazard_",1:5)))
      details[[id]] <- data.frame(replication=replication,start=label,
        parameter=names(coefficient(truth)),truth=unname(coefficient(truth)),
        estimate=unname(coefficient(fit$params)))
      write.csv(do.call(rbind,results),file.path(outdir,"fit_summary_latest.csv"),row.names=FALSE)
      write.csv(do.call(rbind,details),file.path(outdir,"coefficient_recovery_latest.csv"),row.names=FALSE)
      print(summary,row.names=FALSE,digits=8)
    }
  }
},finally=parallel::stopCluster(cluster))
summary <- do.call(rbind,results)
saveRDS(list(summary=summary,coefficients=do.call(rbind,details),truth=truth,
  all_converged=all(summary$converged),n=n,replications=reps,
  empirical_fit_changed=FALSE,source_md5=fingerprint),
  file.path(outdir,"recovery_fits_latest.rds"))
print(summary,row.names=FALSE,digits=8)
