# Estimate the tail-safe piecewise-hazard model with contaminated follow-up
# timegap categories. Point estimates and numerical observed-information
# diagnostics only; no bootstrap.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
baseline_file <- "EM-tenure/output/results/piecewise_hazard/fits_latest.rds"
if (!file.exists(baseline_file)) stop(
  "Run EM-tenure/estimate_piecewise_hazard_tenure_contamination.R first")
baseline_saved <- readRDS(baseline_file)
baseline <- baseline_saved$piecewise
if (nrow(baseline$gamma) != nrow(df_fit)) stop(
  "Saved piecewise fit does not match the current exact-cell data")

starts <- lapply(c(.02,.10,.30),function(e) {
  p <- baseline$params; p$eps_d <- e; p
})
starts[[3]]$lambda_d <- starts[[3]]$lambda_d*c(1.10,.90,1.08,.92,1.05)
starts[[3]]$lambda_g <- starts[[3]]$lambda_g*c(.95,1.06,.94,1.05,1.02)

message(sprintf("Timegap-contamination sample: N=%s; exact cells=%s",
  format(attr(df_fit,"n_original"),big.mark=","),
  format(attr(df_fit,"n_cells"),big.mark=",")))
message("Estimating the 14-parameter model from three starts")
preliminary <- fit_eps_piecewise_multistart(df_fit,starts,
  maxit=450L,reltol=1e-9,pgtol=1e-6,verbose=0L,
  timegap_contamination=TRUE)
message("Refining the best solution with BFGS")
refined <- fit_eps_piecewise_multistart(df_fit,
  list(preliminary$best$params),maxit=700L,reltol=1e-11,pgtol=1e-7,
  method="BFGS",verbose=0L,timegap_contamination=TRUE)
fit <- refined$best
multistart <- rbind(transform(preliminary$summary,phase="preliminary"),
  transform(refined$summary,phase="refinement"))

fits <- list("Piecewise, timegap exact"=baseline,
             "Piecewise, timegap contaminated"=fit)
model_comparison <- do.call(rbind,lapply(seq_along(fits),function(k) {
  f <- fits[[k]]; rates <- duration_weighted_transition_rates(df_fit,f)[1,]
  p <- f$params
  data.frame(model=names(fits)[k],parameters=12L+k,loglik=f$loglik,
    AIC=-2*f$loglik+2*(12L+k),alpha=p$alpha,pi=p$pi,eps_g=p$eps,
    eps_d=if(is.null(p$eps_d)) 0 else p$eps_d,
    weighted_exit=rates$exit_rate,weighted_entry=rates$entry_rate,
    mean_employment_months=12*duration_mean_years(p$lambda_g,p$beta_g),
    mean_nonemployment_months=12*duration_mean_years(p$lambda_d,p$beta_d))
}))
lr <- data.frame(comparison="Timegap contamination versus boundary eps_d=0",
  LR=2*(fit$loglik-baseline$loglik),df=1L,
  p_chibar2=.5*pchisq(2*(fit$loglik-baseline$loglik),1,lower.tail=FALSE))
hazards <- piecewise_hazard_table(fit)
representative <- c(.125,.625,2,4,7.5)
hazards$quarterly_exit <- .duration_transition_probability(
  representative,fit$params$lambda_g,0)
hazards$quarterly_entry <- .duration_transition_probability(
  representative,fit$params$lambda_d,0)
contamination <- timegap_contamination_diagnostics(df_fit,fit)

# Conditional likelihood slice: useful for detecting a boundary solution or a
# locally flat contamination parameter without repeatedly optimizing 13
# nuisance parameters.
slice_eps <- sort(unique(c(.001,.005,.01,.025,.05,.10,.20,.35,
                           fit$params$eps_d)))
slice <- do.call(rbind,lapply(slice_eps,function(e) {
  p <- fit$params; p$eps_d <- e
  z <- .piecewise_eps_pack(p,timegap_contamination=TRUE)
  data.frame(eps_d=e,loglik=-fit$objective_function(z)*sum(df_fit$weight))
}))
slice$loglik_difference <- slice$loglik-max(slice$loglik)

# Numerical observed information on the optimizer scale. Because the objective
# is the average negative log likelihood, divide its inverse by total weight.
message("Computing numerical observed-information diagnostics")
hessian <- tryCatch(optimHess(fit$par_unconstrained,
  fit$objective_function,control=list(ndeps=rep(1e-3,
    length(fit$par_unconstrained)))),error=function(e) e)
if (inherits(hessian,"error")) {
  hessian_diagnostics <- data.frame(rank=NA_integer_,parameters=14L,
    minimum_eigenvalue=NA_real_,condition_number=NA_real_,eps_d_se=NA_real_,
    error=conditionMessage(hessian))
} else {
  hs <- (hessian+t(hessian))/2
  ev <- eigen(hs,symmetric=TRUE,only.values=TRUE)$values
  tol <- max(abs(ev))*sqrt(.Machine$double.eps)
  vcov_z <- tryCatch(solve(hs)/sum(df_fit$weight),error=function(e) NULL)
  d_eps <- with(fit$params,eps_d*(1-eps_d/.95))
  eps_se <- if(is.null(vcov_z)) NA_real_ else
    d_eps*sqrt(max(vcov_z["eps_d","eps_d"],0))
  hessian_diagnostics <- data.frame(rank=sum(ev>tol),parameters=length(ev),
    minimum_eigenvalue=min(ev),condition_number=max(ev)/min(ev),
    eps_d_se=eps_se,error=NA_character_)
}

cat("\nModel comparison\n")
print(model_comparison,row.names=FALSE,digits=8)
cat("\nBoundary likelihood-ratio diagnostic\n")
print(lr,row.names=FALSE,digits=8)
cat("\nPiecewise hazards and implied quarterly risks\n")
print(hazards,row.names=FALSE,digits=7)
cat("\nPosterior timegap-contamination diagnostics\n")
print(contamination,row.names=FALSE,digits=7)
cat("\nMultistart diagnostics\n")
print(multistart,row.names=FALSE,digits=8)
cat("\nConditional likelihood slice\n")
print(slice,row.names=FALSE,digits=8)
cat("\nObserved-information diagnostics\n")
print(hessian_diagnostics,row.names=FALSE,digits=8)
cat("\nRefined scaled gradient\n")
print(fit$gradient,digits=8)

outdir <- "EM-tenure/output/results/timegap_contamination"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(model_comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(hazards,file.path(outdir,"hazards_latest.csv"),row.names=FALSE)
write.csv(contamination,file.path(outdir,"contamination_diagnostics_latest.csv"),
  row.names=FALSE)
write.csv(multistart,file.path(outdir,"multistart_latest.csv"),row.names=FALSE)
write.csv(slice,file.path(outdir,"likelihood_slice_latest.csv"),row.names=FALSE)
write.csv(hessian_diagnostics,file.path(outdir,"hessian_latest.csv"),row.names=FALSE)
fit$objective_function <- NULL
saveRDS(list(fit=fit,baseline=baseline,comparison=model_comparison,lr=lr,
  hazards=hazards,contamination=contamination,multistart=multistart,
  slice=slice,hessian=hessian_diagnostics),file.path(outdir,"fits_latest.rds"))
