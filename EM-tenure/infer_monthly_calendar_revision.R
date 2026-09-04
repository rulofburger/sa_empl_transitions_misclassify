# Full observed-information and delta-method inference for the preferred
# January-heaping plus whole-year revision models. The gross-anchor model has
# 579 finite-difference points for 17 parameters; the clean-anchor extension
# has 649 points for 18 parameters. Both are evaluated in resumable chunks.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
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
total_weight <- sum(df_fit$weight)

inference_model <- tolower(Sys.getenv("REVISION_INFERENCE_MODEL","gross_anchor"))
if (!inference_model %in% c("gross_anchor","clean_anchor"))
  stop("REVISION_INFERENCE_MODEL must be gross_anchor or clean_anchor")
model_dir <- if (inference_model=="clean_anchor")
  "calendar_clean_anchor_revision" else "calendar_revision"
fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  model_dir,"/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate the requested revision model first")
fit <- if (inference_model=="clean_anchor") readRDS(fit_file)$extension else
  readRDS(fit_file)$revision
z0 <- fit$par_unconstrained
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  model_dir,"/inference")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("REVISION_INFERENCE_WORKERS","1"))
chunk_multiplier <- as.integer(Sys.getenv("REVISION_INFERENCE_CHUNK","2"))
resume <- identical(tolower(Sys.getenv("REVISION_INFERENCE_RESUME","true")),
  "true")

cluster <- parallel::makePSOCKcluster(workers)
on.exit(parallel::stopCluster(cluster),add=TRUE)
worker_wd <- getwd()
parallel::clusterCall(cluster,function(path) {
  setwd(path); source("EM-tenure/R/source_all.R"); NULL
},worker_wd)
df_worker <- df_fit
total_weight_worker <- total_weight
parallel::clusterExport(cluster,c("df_worker","total_weight_worker"),
  envir=environment())

unpack_worker <- function(z) .piecewise_calendar_revision_monthly_unpack(z,
  "joint_marginal")
worker_eval <- function(task) {
  p <- unpack_worker(task$z)
  e <- e_step_eps(df_worker,p,check_df=FALSE,suff_stats=FALSE)
  ans <- list(objective=-e$loglik/total_weight_worker,quantities=NULL)
  if (task$quantities) {
    proxy <- list(params=p,gamma=e$gamma)
    rates <- duration_weighted_transition_rates(df_worker,proxy)[1,]
    jp <- e$job_change_posterior
    posterior_q <- sum(df_worker$weight*jp$expected_changes)/
      sum(df_worker$weight*jp$opportunities)
    ans$quantities <- c(entry_rate=rates$entry_rate,
      exit_rate=rates$exit_rate,status_misclassification=p$pi,
      tenure_contamination=p$eps,timegap_contamination=p$eps_d,
      initial_employment=p$alpha,job_change_prob=p$job_change_prob,
      posterior_job_change_prob=posterior_q,
      gross_january_heaping_prob=p$tenure_heaping_prob,
      whole_year_revision_prob=p$tenure_year_revision_prob,
      clean_anchor_revision_prob=if (is.null(
        p$tenure_clean_anchor_revision_prob)) 0 else
        p$tenure_clean_anchor_revision_prob)
  }
  ans
}
parallel::clusterExport(cluster,c("unpack_worker","worker_eval"),
  envir=environment())

step <- 1e-3*pmax(1,abs(z0))
p <- length(z0)
points <- list(list(z=z0,quantities=TRUE))
add_point <- function(z,quantities=FALSE) {
  points[[length(points)+1L]] <<- list(z=z,quantities=quantities)
  length(points)
}
plus <- minus <- integer(p)
for (j in seq_len(p)) {
  zp <- zm <- z0
  zp[j] <- zp[j]+step[j]; zm[j] <- zm[j]-step[j]
  plus[j] <- add_point(zp,TRUE); minus[j] <- add_point(zm,TRUE)
}
pp <- pm <- mp <- mm <- matrix(NA_integer_,p,p)
for (j in seq_len(p-1L)) for (k in (j+1L):p) {
  zpp <- zpm <- zmp <- zmm <- z0
  zpp[j] <- zpp[j]+step[j]; zpp[k] <- zpp[k]+step[k]
  zpm[j] <- zpm[j]+step[j]; zpm[k] <- zpm[k]-step[k]
  zmp[j] <- zmp[j]-step[j]; zmp[k] <- zmp[k]+step[k]
  zmm[j] <- zmm[j]-step[j]; zmm[k] <- zmm[k]-step[k]
  pp[j,k] <- add_point(zpp); pm[j,k] <- add_point(zpm)
  mp[j,k] <- add_point(zmp); mm[j,k] <- add_point(zmm)
}

partial_file <- file.path(outdir,"observed_information_points_partial.rds")
values <- vector("list",length(points))
if (resume && file.exists(partial_file)) {
  partial <- readRDS(partial_file)
  if (identical(names(partial$z0),names(z0)) &&
      isTRUE(all.equal(unname(partial$z0),unname(z0),tolerance=1e-10)) &&
      isTRUE(all.equal(partial$step,step,tolerance=1e-12)) &&
      length(partial$values)==length(values)) values <- partial$values
}
chunk_size <- max(workers,workers*chunk_multiplier)
repeat {
  missing_index <- which(vapply(values,is.null,logical(1)))
  if (!length(missing_index)) break
  index <- head(missing_index,chunk_size)
  message(sprintf("Evaluating information points %d-%d of %d (%d saved)",
    min(index),max(index),length(points),length(points)-length(missing_index)))
  fresh <- parallel::parLapplyLB(cluster,points[index],worker_eval)
  values[index] <- fresh
  saveRDS(list(z0=z0,step=step,values=values),partial_file)
}

f0 <- values[[1L]]$objective
hessian <- matrix(0,p,p,dimnames=list(names(z0),names(z0)))
for (j in seq_len(p)) hessian[j,j] <-
  (values[[plus[j]]]$objective-2*f0+values[[minus[j]]]$objective)/step[j]^2
for (j in seq_len(p-1L)) for (k in (j+1L):p) {
  hessian[j,k] <- hessian[k,j] <-
    (values[[pp[j,k]]]$objective-values[[pm[j,k]]]$objective-
     values[[mp[j,k]]]$objective+values[[mm[j,k]]]$objective)/
    (4*step[j]*step[k])
}
hessian <- (hessian+t(hessian))/2
eigenvalues <- eigen(hessian,symmetric=TRUE,only.values=TRUE)$values
tolerance <- max(abs(eigenvalues))*sqrt(.Machine$double.eps)
diagnostics <- data.frame(parameters=p,rank=sum(eigenvalues>tolerance),
  minimum_eigenvalue=min(eigenvalues),maximum_eigenvalue=max(eigenvalues),
  condition_number=if (min(eigenvalues)>0)
    max(eigenvalues)/min(eigenvalues) else Inf,step_scale=1e-3)
saveRDS(list(z0=z0,step=step,hessian=hessian,
  eigenvalues=eigenvalues,diagnostics=diagnostics),
  file.path(outdir,"observed_information_raw_latest.rds"))
write.csv(diagnostics,file.path(outdir,"hessian_diagnostics_latest.csv"),
  row.names=FALSE)
if (min(eigenvalues)<=tolerance) stop(
  "Revision observed information is not positive definite; min eigenvalue=",
  signif(min(eigenvalues),6))
vcov_z <- solve(hessian)/total_weight

estimates <- values[[1L]]$quantities
jacobian <- matrix(NA_real_,length(estimates),p,
  dimnames=list(names(estimates),names(z0)))
for (j in seq_len(p)) jacobian[,j] <-
  (values[[plus[j]]]$quantities-values[[minus[j]]]$quantities)/(2*step[j])
vcov_quantities <- jacobian%*%vcov_z%*%t(jacobian)
se <- sqrt(pmax(diag(vcov_quantities),0))
summary <- data.frame(quantity=names(estimates),estimate=unname(estimates),
  se=unname(se),lower=unname(estimates-1.96*se),
  upper=unname(estimates+1.96*se))
inference <- list(summary=summary,diagnostics=diagnostics,hessian=hessian,
  vcov_optimizer=vcov_z,jacobian=jacobian,
  vcov_quantities=vcov_quantities)
saveRDS(inference,file.path(outdir,"analytical_inference_latest.rds"))
write.csv(summary,file.path(outdir,"analytical_se_latest.csv"),row.names=FALSE)

cat("\nCalendar-revision observed-information and delta-method inference\n")
print(transform(summary,estimate=100*estimate,se=100*se,
  lower=100*lower,upper=100*upper),row.names=FALSE,digits=7)
cat("\nHessian diagnostics\n")
print(diagnostics,row.names=FALSE,digits=7)
