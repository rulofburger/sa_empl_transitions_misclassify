if (!file.exists("EM-AR1-4W/R/source_all.R")) stop("Run from project root")
source("EM-AR1-4W/R/source_all.R")

path <- "EM-AR1-4W/output/results/fmm_type_error_table3_column1_4w_latest.rds"
fit <- readRDS(path)
if (is.null(fit$analytical_inference$bread)) stop("Finalize the fit first")
cat("Preparing full Table 3-predictor design...\n")
data <- prepare_fmm_covariates_inconsistency_4w(error_design="table3_column1")
objective <- .direct_fmm_covinc_objective(data,fit$params)
z <- .pack_fmm_covinc_4w(fit$params)
ee <- eigen(fit$analytical_inference$bread,symmetric=TRUE)
direction <- ee$vectors[,which.min(ee$values)]
steps <- c(0, .05, -.05, .10, -.10, .20, -.20, .40, -.40, .80, -.80, 1.20, -1.20)
profile <- data.frame(step=steps,value=vapply(steps,function(s)objective$fn(z+s*direction),numeric(1)))
profile$loglik <- -profile$value*sum(data$weight)
print(profile,row.names=FALSE,digits=10)
best <- profile$step[which.min(profile$value)]
if (best == 0) stop("No improving saddle-profile step found")
cat(sprintf("Restarting from eigenvector step %.3f\n",best))
maxit <- as.integer(Sys.getenv("FMM_POLISH_MAXIT","200"))
opt <- optim(z+best*direction,objective$fn,objective$gr,method="BFGS",
  control=list(maxit=maxit,reltol=1e-12,trace=1,REPORT=20))
detail <- objective$details(opt$par)
p <- .normalize_fmm_covinc_labels(detail$params,data)
e <- e_step_fmm_covinc_4w(data,p)
objective2 <- .direct_fmm_covinc_objective(data,p)
z2 <- .pack_fmm_covinc_4w(p); score <- objective2$gr(z2)
fit$params <- p; fit$e <- e; fit$loglik <- e$loglik
fit$optimizer_code <- opt$convergence
fit$converged <- opt$convergence == 0L && max(abs(score)) < 1e-6
fit$fixed_point_error <- max(abs(score))
fit$analytical_inference <- NULL; fit$derived_inference <- NULL
fit$candidates <- rbind(fit$candidates,
  data.frame(loglik=fit$loglik,converged=fit$converged))
saveRDS(fit,path)
cat(sprintf("After saddle escape: ll=%.3f max|score|=%.3e code=%d\n",
  fit$loglik,fit$fixed_point_error,opt$convergence))
