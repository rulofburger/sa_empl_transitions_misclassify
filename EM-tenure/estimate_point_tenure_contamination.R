# Point estimates only for the unchanged epsilon tenure-contamination model.
# This runner deliberately does not source or execute bootstrap code.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from the project root")
ingest_packages <- c("dplyr", "ggplot2")
missing_packages <- ingest_packages[
  !vapply(ingest_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Packages required by scripts/ingest_data_3waves_SA.R are missing: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")
df_point <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
message(sprintf("Epsilon point sample: N=%s; exact cells=%s; weights normalized to N",
  format(attr(df_point,"n_original"),big.mark=","),
  format(attr(df_point,"n_cells"),big.mark=",")))

specs <- list(eps_free=list(stationary=FALSE,linked=FALSE),
  eps_stationary=list(stationary=TRUE,linked=FALSE),
  eps_linked=list(stationary=FALSE,linked=TRUE),
  eps_stationary_linked=list(stationary=TRUE,linked=TRUE))
results <- lapply(names(specs),function(nm){
  s <- specs[[nm]]; message("Estimating ",nm)
  ans <- em_multistart_eps(df_point,stationary=s$stationary,linked=s$linked,
    K=3L,max_iter=1000L,tol=1e-10,param_tol=1e-7,verbose=0L)
  message(sprintf("  best start=%d; LL=%.6f; fixed residual=%.3g",
    ans$best_start,ans$best$loglik,ans$best$diagnostics$fixedpoint_residual))
  ans
}); names(results) <- names(specs)

rows <- do.call(rbind,lapply(names(results),function(nm){
  fit<-results[[nm]]$best;p<-fit$params
  data.frame(model=nm,status=fit$status,iterations=fit$iterations,loglik=fit$loglik,
    alpha=p$alpha,entry=p$theta0,exit=1-p$theta1,pi=p$pi,eps=p$eps,
    lambda_g=p$lambda_g,lambda_d=p$lambda_d,mean_E_months=12/p$lambda_g,
    mean_N_months=12/p$lambda_d,fixedpoint_residual=fit$diagnostics$fixedpoint_residual)
}))
cat("\nRepaired epsilon point estimates (no bootstrap)\n")
print(rows,row.names=FALSE,digits=7)
cat("\nMulti-start diagnostics\n")
for(nm in names(results)){cat("\n",nm,"\n",sep="");print(results[[nm]]$summary,row.names=FALSE,digits=8)}
outdir <- "EM-tenure/output/results/point";dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(rows,file.path(outdir,"point_estimates_latest.csv"),row.names=FALSE)
saveRDS(results,file.path(outdir,"point_fits_latest.rds"))
