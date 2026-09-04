source("EM-tenure/R/source_all.R")
base<-"EM-tenure/output/results/job_change_monthly/four_wave_ar2"
gate<-readRDS(file.path(base,"validation/gates_latest.rds"))
stopifnot(isTRUE(gate$passed),identical(gate$source_md5,four_wave_ar2_source_fingerprint()))
old<-readRDS("EM-tenure/output/results/job_change_monthly/four_wave_ar1/continuous_clock_empirical/comparison_latest.rds")$best
scale<-c(unname(old$controls$parscale),5,5)
outdir<-file.path(base,"recovery");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
truth<-gate$truth;z<-.pack_four_wave_ar2(truth)
rows<-list()
for (rep in 1:2) {
  d<-simulate_four_wave_ar2(5000,truth,seed=280904L+rep,exact_risk=TRUE)
  df<-collapse_eps_cells_4w(d,extra_cols=paste0("interview_month",1:4))
  data_path<-file.path(outdir,paste0("sample",rep,".rds"))
  if (!file.exists(data_path)) saveRDS(df,data_path)
  stopifnot(identical(df,readRDS(data_path)))
  target_eval<-evaluate_four_wave_ar2(prepare_four_wave_kernel_data(df),truth)
  target<-duration_weighted_transition_rates_ar2(df,list(params=truth,gamma=target_eval$gamma))[1,]
  starts<-list(near=z+.08*sin(seq_along(z)))
  starts$near[c("ar2_entry","ar2_exit")]<-0
  if (rep==1) {
    starts$alternative<-z-.08*cos(seq_along(z))
    starts$alternative[c("ar2_entry","ar2_exit")]<-c(1.8,1.3)
  }
  for (label in names(starts)) {
    f<-fit_four_wave_ar2(df,.unpack_four_wave_ar2(starts[[label]]),scale,
      file.path(outdir,paste0("rep",rep,"_",label)),workers=4,maxit=180,
      input_md5=tools::md5sum(data_path))
    rate<-duration_weighted_transition_rates_ar2(df,f)[1,]
    rows[[paste(rep,label)]]<-data.frame(replication=rep,start=label,converged=f$converged,
      loglik=f$loglik,score=f$gradient_max,
      ar2_entry_truth=truth$ar2_entry,ar2_entry=f$params$ar2_entry,
      ar2_exit_truth=truth$ar2_exit,ar2_exit=f$params$ar2_exit,
      entry_truth=target$entry_rate,entry=rate$entry_rate,
      exit_truth=target$exit_rate,exit=rate$exit_rate,
      pi_truth=truth$pi,pi=f$params$pi,alpha_truth=truth$alpha,alpha=f$params$alpha)
    write.csv(do.call(rbind,rows),file.path(outdir,"summary.csv"),row.names=FALSE)
  }
}
# Residualized score information: do the two AR2 scores contain independent
# variation after projecting out ALL 33 nuisance score directions?
scores<-do.call(rbind,lapply(1:3,function(r)
  readRDS(file.path(base,paste0("validation/score_rep",r,".rds")))$scores))
nuisance<-scores[,1:33,drop=FALSE]
sdx<-sqrt(colMeans(nuisance^2));nuisance<-sweep(nuisance,2,pmax(sdx,1e-12),"/")
q<-qr(nuisance,tol=1e-8)
residual<-qr.resid(q,scores[,34:35,drop=FALSE])
info<-crossprod(residual)
eig<-eigen(info,symmetric=TRUE)$values
summary<-do.call(rbind,rows)
saveRDS(list(summary=summary,all_converged=all(summary$converged),truth=truth,
  efficient_score_information=info,information_eigenvalues=eig,nuisance_rank=q$rank,
  information_condition=max(eig)/min(eig),source_md5=four_wave_ar2_source_fingerprint()),
  file.path(outdir,"recovery_latest.rds"))
print(summary,row.names=FALSE,digits=7)
print(info);print(eig)
stopifnot(all(summary$converged),min(eig)>1e-6)
