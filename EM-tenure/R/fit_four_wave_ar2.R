# Full observed-likelihood optimization; neither the AR1 nuisance vector nor
# its posterior is held fixed. Forward differences during search; an independent
# central score is required for acceptance. All files are restartable by label.
fit_four_wave_ar2 <- function(df,start,parscale,prefix,workers=8L,maxit=150L,
    input_md5,score_tolerance=1e-5) {
  fingerprint <- four_wave_ar2_source_fingerprint()
  output <- paste0(prefix,"_fit.rds")
  checkpoint <- paste0(prefix,"_checkpoint.rds")
  if (file.exists(output)) {
    old <- readRDS(output)
    if (identical(old$source_md5,fingerprint) &&
        identical(old$input_md5,input_md5) && isTRUE(old$converged)) return(old)
  }
  z0 <- .pack_four_wave_ar2(start)
  initial_start <- z0
  bounds <- .four_wave_ar2_bounds(z0)
  if (file.exists(checkpoint)) {
    old <- readRDS(checkpoint)
    if (!identical(old$source_md5,fingerprint) || !identical(old$input_md5,input_md5))
      stop("Stale AR2 checkpoint: use a fresh output label")
    z0 <- old$par
  }
  z0 <- pmin(bounds$upper,pmax(bounds$lower,z0))
  stopifnot(length(z0)==35L,length(parscale)==35L,
    all(is.finite(parscale)),all(parscale>0),workers>=1,workers<=8)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  data_worker <- prepare_four_wave_kernel_data(df)
  weight_worker <- sum(df$weight)
  objective_worker <- function(z) {
    -evaluate_four_wave_ar2(data_worker,.unpack_four_wave_ar2(z),
      posterior=FALSE,threads=1L)$loglik/weight_worker
  }
  # Do not capture this fit's cluster/socket environment in serialized closures.
  environment(objective_worker) <- .GlobalEnv
  cluster <- parallel::makePSOCKcluster(workers)
  tryCatch(for (j in seq_along(cluster)) parallel::clusterCall(cluster[j],function(path) {
    setwd(path);source("EM-tenure/R/source_all.R");load_four_wave_monthly_kernel();NULL
  },getwd()),error=function(e) {parallel::stopCluster(cluster);stop(e)})
  parallel::clusterExport(cluster,c("data_worker","weight_worker","objective_worker"),
    envir=environment())
  best <- Inf; evaluations <- gradients <- 0L; history <- data.frame()
  objective <- function(z) {
    value <- -evaluate_four_wave_ar2(data_worker,.unpack_four_wave_ar2(z),
      posterior=FALSE,threads=workers)$loglik/weight_worker
    if (!is.finite(value)) stop("Nonfinite AR2 objective")
    evaluations <<- evaluations+1L
    if (value<best) {
      best <<- value
      history <<- rbind(history,data.frame(evaluation=evaluations,gradient=gradients,
        loglik=-value*weight_worker,timestamp=as.character(Sys.time())))
      saveRDS(list(par=z,loglik=-value*weight_worker,history=history,
        source_md5=fingerprint,input_md5=input_md5),checkpoint)
      message(basename(prefix),": LL ",sprintf("%.6f",-value*weight_worker))
    }
    value
  }
  gradient <- function(z,central=FALSE) {
    step <- 1e-5*pmax(1,abs(z))
    plus <- pmin(step,bounds$upper-z); minus <- pmin(step,z-bounds$lower)
    if (central) {
      points <- vector("list",2*length(z))
      for (j in seq_along(z)) {
        points[[2*j-1]] <- points[[2*j]] <- z
        points[[2*j-1]][j] <- z[j]+plus[j];points[[2*j]][j] <- z[j]-minus[j]
      }
    } else {
      step <- ifelse(plus>=step,step,-minus)
      points <- lapply(seq_along(z),function(j) {v<-z;v[j]<-v[j]+step[j];v})
    }
    values <- unlist(parallel::parLapplyLB(cluster,points,objective_worker,chunk.size=1L),
      use.names=FALSE)
    g <- if (central) (values[seq(1,length(values),2)]-values[seq(2,length(values),2)])/
      (plus+minus) else (values-objective(z))/step
    if (any(!is.finite(g))) stop("Nonfinite AR2 score")
    names(g) <- names(z)
    gradients <<- gradients+1L
    message(basename(prefix),": gradient ",gradients,if (central) " central" else " forward",
      "; max score ",signif(max(abs(g)),5)," (",names(which.max(abs(g))),")")
    g
  }
  tryCatch({
    opt <- optim(z0,objective,gradient,method="L-BFGS-B",lower=bounds$lower,
      upper=bounds$upper,control=list(maxit=maxit,lmm=35,factr=1e5,
        pgtol=score_tolerance*min(parscale),parscale=parscale))
    score <- gradient(opt$par,central=TRUE)
    projected <- score
    projected[opt$par<=bounds$lower+1e-8 & score>0] <- 0
    projected[opt$par>=bounds$upper-1e-8 & score<0] <- 0
    p <- .unpack_four_wave_ar2(opt$par)
    evaluated <- evaluate_four_wave_ar2(data_worker,p,threads=workers)
    fit <- c(list(params=p,par_unconstrained=opt$par,
      integration_method="exact_piecewise",convergence=opt$convergence,message=opt$message,
      converged=opt$convergence==0 && max(abs(projected))<=score_tolerance,
      independent_score_pass=max(abs(projected))<=score_tolerance,
      central_score=score,projected_score=projected,gradient_max=max(abs(projected)),
      initial_start=initial_start,history=history,counts=opt$counts,
      source_md5=fingerprint,input_md5=input_md5,n=attr(df,"n_original") %||% nrow(df),
      controls=list(workers=workers,maxit=maxit,parscale=parscale,step=1e-5,
        score_tolerance=score_tolerance,search_gradient="forward",check_gradient="central"),
      session_info=sessionInfo()),evaluated)
    saveRDS(fit,output)
    print(c(loglik=fit$loglik,score=fit$gradient_max,converged=fit$converged,
      ar2_entry=p$ar2_entry,ar2_exit=p$ar2_exit))
    fit
  },finally=parallel::stopCluster(cluster))
}
