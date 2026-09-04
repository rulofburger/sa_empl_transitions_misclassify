# Observed-information, delta-method, and nuisance-adjusted profile inference
# for the corrected Panel A discrete-month job-change model. Calculations are
# resumable; no bootstrap is run and no Table 7 output is overwritten.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(df_qlfs, allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full, allow_zero_tenure=TRUE)
total_weight <- sum(df_fit$weight)

saved <- readRDS("EM-tenure/output/results/job_change_monthly/fits_latest.rds")
fit <- saved$reset
z0 <- fit$par_unconstrained
q_hat <- fit$params$job_change_prob
outdir <- "EM-tenure/output/results/job_change_monthly/inference"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("MONTHLY_INFERENCE_WORKERS", "1"))
profile_maxit <- as.integer(Sys.getenv("MONTHLY_PROFILE_MAXIT", "40"))
resume <- identical(tolower(Sys.getenv("MONTHLY_INFERENCE_RESUME", "true")),
  "true")

inference_file <- file.path(outdir, "analytical_inference_latest.rds")
cluster <- parallel::makePSOCKcluster(workers)
on.exit(parallel::stopCluster(cluster), add=TRUE)
worker_wd <- getwd()
parallel::clusterCall(cluster, function(path) {
  setwd(path)
  source("EM-tenure/R/source_all.R")
  NULL
}, worker_wd)
df_worker <- df_fit
total_weight_worker <- total_weight
parallel::clusterExport(cluster, c("df_worker", "total_weight_worker"),
  envir=environment())

if (resume && file.exists(inference_file)) {
  message("Loading saved observed-information calculation")
  inference <- readRDS(inference_file)
} else {
  unpack_worker <- function(z) {
    p <- .piecewise_job_change_unpack(z,
      timegap_contamination_model="joint_marginal")
    p$tenure_measurement_model <- "monthly"
    p
  }
  worker_eval <- function(task) {
    p <- unpack_worker(task$z)
    e <- e_step_eps(df_worker, p, check_df=FALSE, suff_stats=FALSE)
    ans <- list(objective=-e$loglik/total_weight_worker, quantities=NULL)
    if (task$quantities) {
      proxy <- list(params=p, gamma=e$gamma)
      rates <- duration_weighted_transition_rates(df_worker, proxy)[1, ]
      jp <- e$job_change_posterior
      posterior_q <- sum(df_worker$weight*jp$expected_changes) /
        sum(df_worker$weight*jp$opportunities)
      ans$quantities <- c(entry_rate=rates$entry_rate,
        exit_rate=rates$exit_rate, status_misclassification=p$pi,
        tenure_contamination=p$eps, timegap_contamination=p$eps_d,
        initial_employment=p$alpha, job_change_prob=p$job_change_prob,
        posterior_job_change_prob=posterior_q)
    }
    ans
  }
  parallel::clusterExport(cluster, c("unpack_worker", "worker_eval"),
    envir=environment())

  step <- 1e-3*pmax(1, abs(z0))
  p <- length(z0)
  points <- list(list(z=z0, quantities=TRUE))
  add_point <- function(z, quantities=FALSE) {
    points[[length(points)+1L]] <<- list(z=z, quantities=quantities)
    length(points)
  }
  plus <- minus <- integer(p)
  for (j in seq_len(p)) {
    zp <- zm <- z0
    zp[j] <- zp[j]+step[j]; zm[j] <- zm[j]-step[j]
    plus[j] <- add_point(zp, TRUE)
    minus[j] <- add_point(zm, TRUE)
  }
  pp <- pm <- mp <- mm <- matrix(NA_integer_, p, p)
  for (j in seq_len(p-1L)) for (k in (j+1L):p) {
    zpp <- zpm <- zmp <- zmm <- z0
    zpp[j] <- zpp[j]+step[j]; zpp[k] <- zpp[k]+step[k]
    zpm[j] <- zpm[j]+step[j]; zpm[k] <- zpm[k]-step[k]
    zmp[j] <- zmp[j]-step[j]; zmp[k] <- zmp[k]+step[k]
    zmm[j] <- zmm[j]-step[j]; zmm[k] <- zmm[k]-step[k]
    pp[j,k] <- add_point(zpp); pm[j,k] <- add_point(zpm)
    mp[j,k] <- add_point(zmp); mm[j,k] <- add_point(zmm)
  }
  message(sprintf("Evaluating %d observed-information points", length(points)))
  values <- parallel::parLapply(cluster, points, worker_eval)
  f0 <- values[[1L]]$objective
  hessian <- matrix(0, p, p, dimnames=list(names(z0), names(z0)))
  for (j in seq_len(p)) hessian[j,j] <-
    (values[[plus[j]]]$objective-2*f0+values[[minus[j]]]$objective)/step[j]^2
  for (j in seq_len(p-1L)) for (k in (j+1L):p) {
    hessian[j,k] <- hessian[k,j] <-
      (values[[pp[j,k]]]$objective-values[[pm[j,k]]]$objective-
       values[[mp[j,k]]]$objective+values[[mm[j,k]]]$objective) /
      (4*step[j]*step[k])
  }
  eigenvalues <- eigen(hessian, symmetric=TRUE, only.values=TRUE)$values
  tolerance <- max(abs(eigenvalues))*sqrt(.Machine$double.eps)
  if (min(eigenvalues) <= tolerance) stop(
    "Monthly observed information is not positive definite; min eigenvalue=",
    signif(min(eigenvalues), 6))
  vcov_z <- solve(hessian)/total_weight

  estimates <- values[[1L]]$quantities
  jacobian <- matrix(NA_real_, length(estimates), p,
    dimnames=list(names(estimates), names(z0)))
  for (j in seq_len(p)) jacobian[,j] <-
    (values[[plus[j]]]$quantities-values[[minus[j]]]$quantities)/(2*step[j])
  vcov_quantities <- jacobian %*% vcov_z %*% t(jacobian)
  se <- sqrt(pmax(diag(vcov_quantities), 0))
  summary <- data.frame(quantity=names(estimates), estimate=unname(estimates),
    se=unname(se), lower=unname(estimates-1.96*se),
    upper=unname(estimates+1.96*se))
  diagnostics <- data.frame(parameters=p, rank=sum(eigenvalues>tolerance),
    minimum_eigenvalue=min(eigenvalues),
    maximum_eigenvalue=max(eigenvalues),
    condition_number=max(eigenvalues)/min(eigenvalues), step_scale=1e-3)
  inference <- list(summary=summary, diagnostics=diagnostics,
    hessian=hessian, vcov_optimizer=vcov_z, jacobian=jacobian,
    vcov_quantities=vcov_quantities)
  saveRDS(inference, inference_file)
  write.csv(summary, file.path(outdir, "analytical_se_latest.csv"),
    row.names=FALSE)
  write.csv(diagnostics, file.path(outdir, "hessian_diagnostics_latest.csv"),
    row.names=FALSE)
}

# Nuisance-adjusted profile points are chosen from the analytical q standard
# error, then each point reoptimizes all 14 nuisance parameters.
q_se <- inference$summary$se[inference$summary$quantity=="job_change_prob"]
if (!is.finite(q_se) || q_se<=0) stop("Invalid analytical reset SE")
q_grid <- sort(unique(pmax(1e-7, c(q_hat,
  q_hat+c(-2.1,-1.5,1.5,2.75)*q_se))))
profile_paths <- file.path(outdir,
  paste0("profile_q_", formatC(q_grid, digits=8, format="f"), ".rds"))
profile_fits <- vector("list", length(q_grid))
missing_profile <- logical(length(q_grid))
for (j in seq_along(q_grid)) {
  if (abs(q_grid[j]-q_hat)<1e-12) {
    profile_fits[[j]] <- fit
  } else if (resume && file.exists(profile_paths[j])) {
    cached <- readRDS(profile_paths[j])
    if (identical(cached$convergence, 0L)) profile_fits[[j]] <- cached else
      missing_profile[j] <- TRUE
  } else missing_profile[j] <- TRUE
}
if (any(missing_profile)) {
  start_worker <- fit$params
  maxit_worker <- profile_maxit
  parallel::clusterExport(cluster, c("start_worker", "maxit_worker"),
    envir=environment())
  message("Reoptimizing nuisance parameters at ", sum(missing_profile),
    " fixed reset probabilities")
  new_fits <- parallel::parLapply(cluster, q_grid[missing_profile], function(q) {
    ans <- fit_eps_piecewise_job_change(df_worker, start_worker, q_fixed=q,
      maxit=maxit_worker, reltol=1e-10, pgtol=1e-7, workers=1L,
      tenure_measurement_model="monthly")
    ans$objective_function <- NULL
    ans
  })
  which_missing <- which(missing_profile)
  for (k in seq_along(which_missing)) {
    j <- which_missing[k]
    profile_fits[[j]] <- new_fits[[k]]
    saveRDS(profile_fits[[j]], profile_paths[j])
  }
}
profile <- data.frame(q=q_grid,
  loglik=vapply(profile_fits, `[[`, numeric(1), "loglik"),
  convergence=vapply(profile_fits, `[[`, integer(1), "convergence"),
  iterations=vapply(profile_fits, function(x) as.integer(x$iterations),
    integer(1)))
profile$LR <- 2*(max(profile$loglik)-profile$loglik)
if (any(profile$convergence != 0L)) stop(
  "At least one nuisance-adjusted profile point did not converge")

# A local quadratic through all reoptimized points provides a smooth profile
# interval; the grid must bracket both 95% roots.
x <- profile$q-q_hat
quad <- lm(loglik~x+I(x^2), data=profile)
target <- max(profile$loglik)-qchisq(.95,1)/2
co <- coef(quad)
profile_roots <- polyroot(c(co[1]-target,co[2],co[3]))
roots <- sort(Re(profile_roots[abs(Im(profile_roots))<1e-8]))+q_hat
if (length(roots)!=2L || roots[1]<min(q_grid) || roots[2]>max(q_grid))
  stop("Profile grid does not bracket the quadratic 95% interval")
profile_interval <- data.frame(estimate=q_hat, analytical_se=q_se,
  lower=roots[1], upper=roots[2], quadratic_R2=summary(quad)$r.squared)

write.csv(profile, file.path(outdir, "profile_latest.csv"), row.names=FALSE)
write.csv(profile_interval, file.path(outdir, "profile_interval_latest.csv"),
  row.names=FALSE)
saveRDS(list(inference=inference, profile=profile,
  profile_interval=profile_interval),
  file.path(outdir, "inference_and_profile_latest.rds"))

cat("\nDiscrete-month observed-information and delta-method inference\n")
print(transform(inference$summary, estimate=100*estimate, se=100*se,
  lower=100*lower, upper=100*upper), row.names=FALSE, digits=7)
cat("\nHessian diagnostics\n")
print(inference$diagnostics, row.names=FALSE, digits=7)
cat("\nNuisance-adjusted reset profile\n")
print(profile, row.names=FALSE, digits=9)
print(transform(profile_interval, estimate=100*estimate,
  analytical_se=100*analytical_se, lower=100*lower, upper=100*upper),
  row.names=FALSE, digits=7)
