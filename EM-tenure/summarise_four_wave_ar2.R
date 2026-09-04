source("EM-tenure/R/source_all.R")
base<-"EM-tenure/output/results/job_change_monthly/four_wave_ar2"
outdir<-file.path(base,"empirical")
completed<-readRDS(file.path(outdir,"completed_starts.rds"))
fingerprint<-four_wave_ar2_source_fingerprint()
stopifnot(isTRUE(completed$all_converged),identical(completed$source_md5,fingerprint))
fits<-lapply(completed$paths,readRDS)
names(fits)<-sub("_fit.rds","",basename(completed$paths),fixed=TRUE)
stopifnot(all(vapply(fits,function(f) isTRUE(f$converged) && f$gradient_max<=1e-5 &&
  identical(f$source_md5,fingerprint),logical(1))))
best_label<-names(fits)[which.max(vapply(fits,`[[`,numeric(1),"loglik"))]
best<-fits[[best_label]]
bounds<-.four_wave_ar2_bounds(best$par_unconstrained)
boundary_parameters<-names(best$par_unconstrained)[
  best$par_unconstrained<=bounds$lower+1e-6 | best$par_unconstrained>=bounds$upper-1e-6]
ar2_boundary<-any(c("ar2_entry","ar2_exit") %in% boundary_parameters)
if (ar2_boundary) warning("An AR2 effect is at its numerical bound; do not interpret it as well determined")
ar1base<-"EM-tenure/output/results/job_change_monthly/four_wave_ar1"
old<-readRDS(file.path(ar1base,"continuous_clock_empirical/comparison_latest.rds"))$best
df<-readRDS(file.path(ar1base,"prepared_cells_latest.rds"))$df4
data<-prepare_four_wave_kernel_data(df)
ref<-e_step_four_wave_ar2_reference(df,best$params)
reference_checks<-c(loglik=best$loglik-ref$loglik,
  posterior=max(abs(best$gamma-ref$gamma)))
stopifnot(abs(reference_checks[1])<1e-6,reference_checks[2]<1e-10,
  best$loglik>=old$loglik-1e-6)
# These aggregates must come from the AR2 reference, not the AR1 reweighting input.
best$duration_reliability_posterior<-ref$duration_reliability_posterior
best$duration_reliability_component_probabilities<-ref$duration_reliability_component_probabilities
best$job_change_posterior<-ref$job_change_posterior
rm(ref);gc()
h<-latent_histories_eps_4w()
post_employment<-best$gamma %*% h
observed<-as.matrix(df[paste0("y",1:4)])
post_errors<-ifelse(observed==1,1-post_employment,post_employment)
posterior_status_error<-sum(df$weight*rowSums(post_errors))/(4*sum(df$weight))
risks<-.ar2_baseline_risks(data,best$params)
exact_ar2_score<-c(ar2_entry=0,ar2_exit=0)
for (t in 2:3) for (j in 1:16) if (h[j,t]!=h[j,t-1]) {
  origin<-if (h[j,t]) "exit" else "entry"
  name<-paste0("ar2_",origin)
  risk<-.ar2_shift(risks[[origin]][,t],best$params[[name]])
  exact_ar2_score[name]<-exact_ar2_score[name]+sum(df$weight*best$gamma[,j]*
    (as.numeric(h[j,t+1]!=h[j,t])-risk))/sum(df$weight)
}
score_identity<-exact_ar2_score+best$central_score[names(exact_ar2_score)]
stopifnot(max(abs(score_identity))<1e-7,
  abs(posterior_status_error-best$params$pi)<5e-6)
all_fits<-c(list(AR1=old),fits)
summary<-do.call(rbind,lapply(names(all_fits),function(label) {
  f<-all_fits[[label]];r<-duration_weighted_transition_rates_ar2(df,f)[1,]
  data.frame(model=label,n=attr(df,"n_original"),entry=r$entry_rate,exit=r$exit_rate,
    misclassification=f$params$pi,initial_employment=f$params$alpha,
    ar2_entry=f$params$ar2_entry %||% 0,ar2_exit=f$params$ar2_exit %||% 0,
    loglik=f$loglik,score=f$gradient_max,converged=f$converged)
}))
intervals<-do.call(rbind,lapply(c("AR1",best_label),function(label)
  data.frame(model=label,duration_weighted_transition_rates_ar2(df,all_fits[[label]]))))
coef<-data.frame(parameter=names(best$par_unconstrained),
  AR1=.pack_four_wave_ar2(old$params),AR2=best$par_unconstrained,row.names=NULL)
hazards<-data.frame(hazard=rep(c("Entry","Exit"),each=5),
  duration=rep(c("0--3 months","3--12 months","1--3 years","3--5 years","5+ years"),2),
  AR1=c(old$params$lambda_d,old$params$lambda_g),AR2=c(best$params$lambda_d,best$params$lambda_g))

# Local identification diagnostic, not a survey-design covariance estimate.
z<-best$par_unconstrained
scores<-matrix(NA_real_,nrow(df),35,dimnames=list(NULL,names(z)))
for (j in seq_along(z)) {
  step<-1e-5*max(1,abs(z[j]));plus<-minus<-z
  plus[j]<-plus[j]+step;minus[j]<-minus[j]-step
  scores[,j]<-(evaluate_four_wave_ar2(data,.unpack_four_wave_ar2(plus),posterior=FALSE)$row_loglik-
    evaluate_four_wave_ar2(data,.unpack_four_wave_ar2(minus),posterior=FALSE)$row_loglik)/(2*step)
  if (j%%5==0) message("Identification scores: ",j,"/35")
}
stopifnot(all(is.finite(scores)))
weighted<-scores*sqrt(df$weight)
scale<-sqrt(colSums(weighted[,1:33,drop=FALSE]^2))
x<-sweep(weighted[,1:33,drop=FALSE],2,pmax(scale,1e-12),"/")
identification<-do.call(rbind,lapply(c(1e-6,1e-8,1e-10),function(tol) {
  q<-qr(x,tol=tol)
  r<-qr.resid(q,weighted[,34:35,drop=FALSE]);info<-crossprod(r)
  eig<-eigen(info,symmetric=TRUE)$values
  data.frame(qr_tolerance=tol,nuisance_rank=q$rank,min_eigenvalue=min(eig),
    max_eigenvalue=max(eig),condition=max(eig)/min(eig),
    entry_information=info[1,1],exit_information=info[2,2],cross_information=info[1,2])
}))
saveRDS(list(scores=scores,identification=identification,source_md5=fingerprint),
  file.path(outdir,"identification_scores.rds"))

# Descriptive posterior decomposition of apparent transitions with a follow-up
# wave. All four reports inform the posterior; this is not a two-report forecast.
decompose<-function(f,label) {
  h<-latent_histories_eps_4w();rows<-list()
  for (origin in 0:1) {
    denom<-raw_reverse<-confirmed_reverse<-confirmed_persist<-0
    for (t in 1:2) {
      observed<-df[[paste0("y",t)]]==origin & df[[paste0("y",t+1)]]==1-origin
      w<-df$weight*observed;denom<-denom+sum(w)
      raw_reverse<-raw_reverse+sum(w*(df[[paste0("y",t+2)]]==origin))
      same_direction<-h[,t]==origin & h[,t+1]==1-origin
      rev<-same_direction & h[,t+2]==origin
      persist<-same_direction & h[,t+2]==1-origin
      confirmed_reverse<-confirmed_reverse+sum(w*rowSums(f$gamma[,rev,drop=FALSE]))
      confirmed_persist<-confirmed_persist+sum(w*rowSums(f$gamma[,persist,drop=FALSE]))
    }
    rows[[origin+1]]<-data.frame(model=label,transition=if (origin) "Exit" else "Entry",
      apparent_weight=denom,raw_reversal=raw_reverse/denom,
      not_confirmed=1-(confirmed_reverse+confirmed_persist)/denom,
      genuine_reversed=confirmed_reverse/denom,genuine_persistent=confirmed_persist/denom)
  }
  do.call(rbind,rows)
}
implications<-rbind(decompose(old,"AR1"),decompose(best,"AR2"))
stopifnot(max(abs(rowSums(implications[c("not_confirmed","genuine_reversed","genuine_persistent")])-1))<1e-10)
write.csv(summary,file.path(outdir,"fit_summary.csv"),row.names=FALSE)
write.csv(intervals,file.path(outdir,"interval_rates.csv"),row.names=FALSE)
write.csv(hazards,file.path(outdir,"hazards.csv"),row.names=FALSE)
write.csv(coef,file.path(outdir,"transformed_coefficients.csv"),row.names=FALSE)
write.csv(identification,file.path(outdir,"identification.csv"),row.names=FALSE)
write.csv(implications,file.path(outdir,"transition_implications.csv"),row.names=FALSE)
saveRDS(list(best=best,best_label=best_label,summary=summary,intervals=intervals,
  hazards=hazards,coefficients=coef,identification=identification,implications=implications,
  reference_checks=reference_checks,likelihood_gain=best$loglik-old$loglik,
  exact_ar2_score=exact_ar2_score,score_identity_difference=score_identity,
  posterior_status_error=posterior_status_error,
  boundary_parameters=boundary_parameters,ar2_boundary=ar2_boundary,
  lr_statistic=2*(best$loglik-old$loglik),initial_pair="alpha_plus_AR1_first_transition",
  standard_errors_estimated=FALSE,source_md5=fingerprint,
  all_R_source_md5=tools::md5sum(list.files("EM-tenure/R",pattern="[.]R$",full.names=TRUE)),
  runner_md5=tools::md5sum(c("EM-tenure/validate_four_wave_ar2.R",
    "EM-tenure/recover_four_wave_ar2.R","EM-tenure/estimate_four_wave_ar2.R",
    "EM-tenure/summarise_four_wave_ar2.R")),
  input_paths=completed$paths,input_md5=tools::md5sum(completed$paths),session_info=sessionInfo()),
  file.path(outdir,"comparison_latest.rds"))
print(summary,row.names=FALSE,digits=8)
print(identification,row.names=FALSE)
print(implications,row.names=FALSE,digits=7)
print(reference_checks)
