source("EM-tenure/R/source_all.R")
base <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
outdir <- "EM-tenure/output/results/job_change_monthly/four_wave_ar2/validation"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
old <- readRDS(file.path(base,"continuous_clock_empirical/comparison_latest.rds"))
stopifnot(isTRUE(old$best$converged),identical(old$source_md5,four_wave_fast_source_fingerprint()))
p <- old$best$params
df <- readRDS(file.path(base,"prepared_cells_latest.rds"))$df4
data <- prepare_four_wave_kernel_data(df)
zero <- evaluate_four_wave_ar2(data,p)
nesting <- c(loglik=zero$loglik-old$best$loglik,
  posterior=max(abs(zero$gamma-old$best$gamma)))
stopifnot(abs(nesting[1])<1e-6,nesting[2]<1e-10)
message("Full-sample zero-effect nesting passed")
p$ar2_entry <- log(3);p$ar2_exit <- log(2)
fast <- evaluate_four_wave_ar2(data,p)
ref <- e_step_four_wave_ar2_reference(df,p)
equivalence <- c(loglik=fast$loglik-ref$loglik,posterior=max(abs(fast$gamma-ref$gamma)))
stopifnot(abs(equivalence[1])<1e-6,equivalence[2]<1e-10)
message("Full-sample nonzero-effect reference check passed")
rm(data,df,fast,ref,zero);gc()

# Exhaustive joint probabilities, including wrong status and duration reports.
mass_p <- p;mass_p$lambda_g<-rep(20,5);mass_p$eps<-.2;mass_p$eps_d<-.3;mass_p$pi<-.15
mass_p$job_change_prob<-0;mass_p$tenure_start_month_probs<-rep(1/12,12)
for (nm in c("tenure_reliability_shift","timegap_reliability_shift","tenure_heaping_prob",
    "tenure_year_revision_prob","tenure_clean_anchor_revision_prob",
    "tenure_exact_anchor_retention_prob","tenure_local_revision_prob")) mass_p[[nm]]<-0
max_month<-24L
grid<-as.matrix(expand.grid(rep(list(seq_len(7L+max_month+1L)),4)))
df<-data.frame(weight=rep(1,nrow(grid)))
for (t in 1:4) {
  df[[paste0("y",t)]]<-as.integer(grid[,t]>7)
  df[[paste0("tenure",t)]]<-ifelse(grid[,t]>7,(grid[,t]-8)/12,NA_real_)
  df[[paste0("timegap_cat",t)]]<-ifelse(grid[,t]<=7,grid[,t],NA_integer_)
  df[[paste0("interview_month",t)]]<-3L*t
}
data<-prepare_four_wave_kernel_data(df);rm(grid,df);gc()
cases<-expand.grid(varying_hazards=c(FALSE,TRUE),opposite_effects=c(FALSE,TRUE))
cases$mass<-NA_real_
for (i in seq_len(nrow(cases))) {
  mass_p$lambda_d<-if (cases$varying_hazards[i]) c(4,1,.5,.3,.2) else rep(.4,5)
  mass_p$ar2_entry<-.8;mass_p$ar2_exit<-if (cases$opposite_effects[i]) -.6 else .6
  cases$mass[i]<-sum(exp(evaluate_four_wave_ar2(data,mass_p,posterior=FALSE)$row_loglik))
  print(cases[i,],digits=12)
}
tail_bound<-64*exp(-20*(max_month+1-9)/12)
stopifnot(all(abs(cases$mass-1)<tail_bound+1e-10))
rm(data);gc()

# Scores at known nonzero generating effects; all 35 coordinates are checked.
n<-10000L;reps<-3L;records<-list()
z<-.pack_four_wave_ar2(p)
for (r in 1:reps) {
  message("AR2 generator score check ",r,"/",reps)
  d<-simulate_four_wave_ar2(n,p,seed=270904L+r,exact_risk=TRUE)
  data<-prepare_four_wave_kernel_data(d)
  scores<-matrix(NA_real_,n,35,dimnames=list(NULL,names(z)))
  for (j in seq_along(z)) {
    step<-1e-5*max(1,abs(z[j]));plus<-minus<-z
    plus[j]<-plus[j]+step;minus[j]<-minus[j]-step
    scores[,j]<-(evaluate_four_wave_ar2(data,.unpack_four_wave_ar2(plus),posterior=FALSE)$row_loglik-
      evaluate_four_wave_ar2(data,.unpack_four_wave_ar2(minus),posterior=FALSE)$row_loglik)/(2*step)
  }
  stopifnot(all(is.finite(scores)))
  records[[r]]<-data.frame(replication=r,parameter=names(z),n=n,
    mean_score=colMeans(scores),mcse=apply(scores,2,sd)/sqrt(n))
  saveRDS(list(data=d,truth=p,scores=scores,seed=270904L+r,
    source_md5=four_wave_ar2_source_fingerprint()),file.path(outdir,paste0("score_rep",r,".rds")))
}
draws<-do.call(rbind,records)
pooled<-do.call(rbind,lapply(split(draws,draws$parameter),function(x)
  data.frame(parameter=x$parameter[1],mean_score=mean(x$mean_score),
    mcse=sqrt(sum(x$mcse^2))/nrow(x))))
pooled$z_score<-pooled$mean_score/pooled$mcse
score_pass<-all(abs(pooled$z_score)<5)
write.csv(pooled,file.path(outdir,"score_summary.csv"),row.names=FALSE)
write.csv(cases,file.path(outdir,"probability_mass.csv"),row.names=FALSE)
saveRDS(list(passed=score_pass,nesting=nesting,equivalence=equivalence,mass=cases,
  tail_bound=tail_bound,scores=pooled,truth=p,n_simulated=n*reps,
  source_md5=four_wave_ar2_source_fingerprint()),file.path(outdir,"gates_latest.rds"))
print(pooled[order(-abs(pooled$z_score)),],row.names=FALSE)
stopifnot(score_pass)
