# Robustness of the piecewise model to the timegap contamination distribution.
# Estimates (i) local category contamination and (ii) joint per-wave marginal
# contamination. No bootstrap or full Hessian.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
exact_file <- "EM-tenure/output/results/piecewise_hazard/fits_latest.rds"
marginal_file <- "EM-tenure/output/results/timegap_contamination/fits_latest.rds"
if (!file.exists(exact_file) || !file.exists(marginal_file)) stop(
  "Run the piecewise and timegap-contamination estimators first")
exact <- readRDS(exact_file)$piecewise
marginal <- readRDS(marginal_file)$fit
if (nrow(marginal$gamma)!=nrow(df_fit)) stop("Saved fit does not match data")

outdir <- "EM-tenure/output/results/timegap_contamination_robustness"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

estimate_variant <- function(label,model,starts,local_decay=1) {
  cache <- file.path(outdir,paste0(label,"_fit.rds"))
  if (file.exists(cache) &&
      identical(tolower(Sys.getenv("TIMEGAP_ROBUST_RESUME")),"true")) {
    message("Loading saved ",label," fit")
    return(readRDS(cache))
  }
  message("Estimating ",label," model from two starts")
  preliminary <- fit_eps_piecewise_multistart(df_fit,starts,maxit=350L,
    reltol=1e-9,pgtol=1e-6,method="L-BFGS-B",verbose=0L,
    timegap_contamination=TRUE,timegap_contamination_model=model,
    timegap_local_decay=local_decay)
  message("Refining ",label," optimum")
  refined <- fit_eps_piecewise_multistart(df_fit,list(preliminary$best$params),
    maxit=600L,reltol=1e-11,pgtol=1e-7,method="BFGS",verbose=0L,
    timegap_contamination=TRUE,timegap_contamination_model=model,
    timegap_local_decay=local_decay)
  report <- rbind(transform(preliminary$summary,phase="preliminary"),
    transform(refined$summary,phase="refinement"))
  ans <- list(fit=refined$best,multistart=report)
  ans$fit$objective_function <- NULL
  saveRDS(ans,cache)
  ans
}

local_starts <- list(marginal$params,marginal$params)
local_starts[[1]]$eps_d <- .20
local_starts[[2]]$eps_d <- .40
local_starts[[2]]$lambda_d <- local_starts[[2]]$lambda_d*
  c(1.05,.95,1.04,.96,1.02)
local_result <- estimate_variant("local","local",local_starts,local_decay=1)

joint_start <- marginal$params
joint_start$eps_d <- 1-sqrt(1-marginal$params$eps_d)
joint_starts <- list(joint_start,joint_start)
joint_starts[[1]]$eps_d <- .08
joint_starts[[2]]$eps_d <- .25
joint_starts[[2]]$lambda_d <- joint_starts[[2]]$lambda_d*
  c(.96,1.05,.95,1.04,.98)
joint_result <- estimate_variant("joint_marginal","joint_marginal",joint_starts)

fits <- list("Timegap exact"=exact,"Independent marginal follow-up"=marginal,
  "Local category follow-up"=local_result$fit,
  "Joint per-wave marginal"=joint_result$fit)
model_comparison <- do.call(rbind,lapply(seq_along(fits),function(k) {
  f <- fits[[k]]; p <- f$params
  npar <- if(k==1)13L else 14L
  rates <- duration_weighted_transition_rates(df_fit,f)[1,]
  eps_d <- if(is.null(p$eps_d)) 0 else p$eps_d
  joint <- identical(p$timegap_contamination_model,"joint_marginal")
  data.frame(model=names(fits)[k],parameters=npar,
    loglik=f$loglik,AIC=-2*f$loglik+2*npar,
    alpha=p$alpha,pi=p$pi,eps_g=p$eps,eps_d=eps_d,
    effective_pair_contamination=if(joint)1-(1-eps_d)^2 else eps_d,
    weighted_exit=rates$exit_rate,weighted_entry=rates$entry_rate,
    mean_employment_months=12*duration_mean_years(p$lambda_g,p$beta_g),
    mean_nonemployment_months=12*duration_mean_years(p$lambda_d,p$beta_d))
}))
hazards <- do.call(rbind,lapply(names(fits)[-1],function(nm) {
  f <- fits[[nm]]; tab <- piecewise_hazard_table(f)
  tab$model <- nm
  tab$quarterly_exit <- .duration_transition_probability(
    c(.125,.625,2,4,7.5),f$params$lambda_g,0)
  tab$quarterly_entry <- .duration_transition_probability(
    c(.125,.625,2,4,7.5),f$params$lambda_d,0)
  tab
}))
multistart <- rbind(transform(local_result$multistart,model="local"),
  transform(joint_result$multistart,model="joint_marginal"))

cat("\nTimegap-contamination robustness comparison\n")
print(model_comparison,row.names=FALSE,digits=8)
cat("\nDuration-specific hazards\n")
print(hazards[,c("model","interval","exit_hazard","entry_hazard",
  "quarterly_exit","quarterly_entry")],row.names=FALSE,digits=7)
cat("\nMultistart and refinement diagnostics\n")
print(multistart,row.names=FALSE,digits=8)
cat("\nLocal-model scaled gradient\n")
print(local_result$fit$gradient,digits=8)
cat("\nJoint-model scaled gradient\n")
print(joint_result$fit$gradient,digits=8)

write.csv(model_comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(hazards,file.path(outdir,"hazards_latest.csv"),row.names=FALSE)
write.csv(multistart,file.path(outdir,"multistart_latest.csv"),row.names=FALSE)
saveRDS(list(local=local_result,joint=joint_result,exact=exact,
  marginal=marginal,comparison=model_comparison,hazards=hazards,
  multistart=multistart),file.path(outdir,"fits_latest.rds"))
