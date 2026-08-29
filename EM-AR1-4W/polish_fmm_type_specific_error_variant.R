if (!file.exists("EM-AR1-4W/R/source_all.R")) stop("Run from project root")
source("EM-AR1-4W/R/source_all.R")

args <- commandArgs(trailingOnly=TRUE)
variant <- if (length(args)) args[[1]] else "intercept_only"
lookup <- list(
  intercept_only=list(stem="type_error_intercept_only",design="intercept_only"),
  table3_column1=list(stem="type_error_table3_column1",design="table3_column1"))
if (!variant %in% names(lookup)) stop("Unknown variant: ",variant)
cfg <- lookup[[variant]]
path <- file.path("EM-AR1-4W/output/results",
                  paste0("fmm_",cfg$stem,"_4w_latest.rds"))
if (!file.exists(path)) stop("Missing initial fit: ",path)

cat("Preparing data for ",variant," polish...\n",sep="")
data <- prepare_fmm_covariates_inconsistency_4w(error_design=cfg$design)
fit <- readRDS(path); objective <- .direct_fmm_covinc_objective(data,fit$params)
z0 <- .pack_fmm_covinc_4w(fit$params)
g0 <- objective$gr(z0); f0 <- objective$fn(z0)
cat(sprintf("Before: ll=%.3f max|score|=%.3e\n",-f0*sum(data$weight),max(abs(g0))))

maxit <- as.integer(Sys.getenv("FMM_POLISH_MAXIT","120"))
opt <- optim(z0,objective$fn,objective$gr,method="BFGS",
             control=list(maxit=maxit,reltol=1e-12,trace=1,REPORT=20))
detail <- objective$details(opt$par)
p <- .normalize_fmm_covinc_labels(detail$params,data)
e <- e_step_fmm_covinc_4w(data,p)
# Rebuild after label normalization before evaluating and storing diagnostics.
objective2 <- .direct_fmm_covinc_objective(data,p)
z <- .pack_fmm_covinc_4w(p); score <- objective2$gr(z)
fit$params <- p; fit$e <- e; fit$loglik <- e$loglik
fit$iterations <- (fit$iterations %||% 0L) + unname(opt$counts[["function"]])
fit$optimizer_code <- opt$convergence
fit$converged <- opt$convergence == 0L && max(abs(score)) < 1e-6
fit$fixed_point_error <- max(abs(score))
fit$candidates <- rbind(fit$candidates,
  data.frame(loglik=fit$loglik,converged=fit$converged))
fit$analytical_inference <- NULL; fit$derived_inference <- NULL
saveRDS(fit,path)
cat(sprintf("After: ll=%.3f max|score|=%.3e code=%d converged=%s\n",
  fit$loglik,fit$fixed_point_error,opt$convergence,fit$converged))
