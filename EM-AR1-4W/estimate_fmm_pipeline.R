if(!file.exists("EM-AR1-4W/R/source_all.R"))stop("Run from project root")
source("EM-AR1-4W/R/source_all.R");source("scripts/ingest_data_4waves_SA.R")
df<-df_qlfs[,c("y1","y2","y3","y4","weight")];df$weight<-nrow(df)*df$weight/sum(df$weight)
configs<-list(none_stat=c("none",TRUE),sym_stat=c("symmetric",TRUE),
              none_free=c("none",FALSE),sym_free=c("symmetric",FALSE))
fits<-lapply(configs,function(x)fit_fmm_ar1_4w(df,x[[1]],as.logical(x[[2]]),random_starts=60L))
inference<-lapply(fits,function(x)analytical_se_fmm_ar1_4w(df,x))
outdir<-"EM-AR1-4W/output/results";dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
run_id<-format(Sys.time(),"%Y%m%d_%H%M%S")
for(nm in names(fits)){fits[[nm]]$analytical_inference<-inference[[nm]]
 saveRDS(fits[[nm]],file.path(outdir,paste0("fmm_ar1_4w_",nm,"_",run_id,".rds")))}
rows<-do.call(rbind,lapply(names(fits),function(nm)cbind(model=nm,inference[[nm]]$summary)))
diag<-do.call(rbind,lapply(names(fits),function(nm)data.frame(model=nm,loglik=fits[[nm]]$loglik,
 converged=fits[[nm]]$converged,identified=fits[[nm]]$identified,n_obs=fits[[nm]]$sample$n,
 max_abs_score=fits[[nm]]$diagnostics$max_abs_score,
 information_condition=fits[[nm]]$diagnostics$information_condition)))
write.csv(rows,file.path(outdir,"fmm_analytical_se_latest.csv"),row.names=FALSE)
write.csv(diag,file.path(outdir,"fmm_run_summary_latest.csv"),row.names=FALSE)
cat("\nTwo-type four-wave AR(1) FMM estimates (percent)\n")
print(transform(rows,estimate=100*estimate,se=100*se),row.names=FALSE,digits=5)
cat("\nDiagnostics\n");print(diag,row.names=FALSE,digits=5)
