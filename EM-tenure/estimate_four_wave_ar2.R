source("EM-tenure/R/source_all.R")
base<-"EM-tenure/output/results/job_change_monthly/four_wave_ar2"
gate<-readRDS(file.path(base,"validation/gates_latest.rds"))
recovery<-readRDS(file.path(base,"recovery/recovery_latest.rds"))
fingerprint<-four_wave_ar2_source_fingerprint()
stopifnot(isTRUE(gate$passed),isTRUE(recovery$all_converged),
  identical(gate$source_md5,fingerprint),identical(recovery$source_md5,fingerprint),
  min(recovery$information_eigenvalues)>1e-6)
ar1base<-"EM-tenure/output/results/job_change_monthly/four_wave_ar1"
old<-readRDS(file.path(ar1base,"continuous_clock_empirical/comparison_latest.rds"))$best
data_path<-file.path(ar1base,"prepared_cells_latest.rds")
df<-readRDS(data_path)$df4
data<-prepare_four_wave_kernel_data(df)
p<-old$params
scale<-c(unname(old$controls$parscale),5,5)
outdir<-file.path(base,"empirical");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# A cheap conditional warm start only. The following fits free all 35 parameters.
conditional<-optim(c(0,0),function(b) {
  p$ar2_entry<-b[1];p$ar2_exit<-b[2]
  -evaluate_four_wave_ar2(data,p,posterior=FALSE,base_fit=old)$loglik/sum(df$weight)
},method="L-BFGS-B",lower=c(-8,-8),upper=c(8,8),control=list(factr=1e5))
p$ar2_entry<-conditional$par[1];p$ar2_exit<-conditional$par[2]
saveRDS(conditional,file.path(outdir,"conditional_warm_start.rds"))
print(c(conditional_entry=p$ar2_entry,conditional_exit=p$ar2_exit))
fits<-list()
fits$primary<-fit_four_wave_ar2(df,p,scale,file.path(outdir,"primary"),
  input_md5=tools::md5sum(data_path))
z<-.pack_four_wave_ar2(fits$primary$params)
bounds<-.four_wave_ar2_bounds(z)
for (sgn in c(-1,1)) {
  start<-z+sgn*.08*sin(seq_along(z))
  start[c("ar2_entry","ar2_exit")]<-start[c("ar2_entry","ar2_exit")]+sgn*.4
  start<-pmin(bounds$upper,pmax(bounds$lower,start))
  label<-if (sgn<0) "minus" else "plus"
  fits[[label]]<-fit_four_wave_ar2(df,.unpack_four_wave_ar2(start),scale,
    file.path(outdir,label),input_md5=tools::md5sum(data_path))
}
stopifnot(all(vapply(fits,function(f) isTRUE(f$converged),logical(1))))
saveRDS(list(paths=file.path(outdir,paste0(names(fits),"_fit.rds")),
  all_converged=TRUE,source_md5=fingerprint),file.path(outdir,"completed_starts.rds"))
