# Estimate the tail-safe five-interval piecewise duration-hazard epsilon model.
# Reuses the corrected constant and power-law fits; no bootstrap.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
previous_file <- "EM-tenure/output/results/duration_hazard/fits_latest.rds"
if (!file.exists(previous_file)) stop(
  "Run EM-tenure/estimate_duration_hazard_tenure_contamination.R first")
previous <- readRDS(previous_file)
constant_fit <- previous$constant_direct
power_fit <- previous$duration
if (nrow(power_fit$gamma) != nrow(df_fit)) stop(
  "Saved duration fit does not match the current exact-cell data")

start_power <- piecewise_start_from_power(power_fit$params)
start_constant <- constant_fit$params
start_constant$lambda_g <- start_constant$lambda_g*c(1.10,1,.90,.85,.80)
start_constant$lambda_d <- start_constant$lambda_d*c(1.20,1,.85,.80,.75)
start_constant$beta_g <- start_constant$beta_d <- 0
start_perturbed <- start_power
start_perturbed$lambda_g <- start_power$lambda_g*c(1.20,.85,1.12,.88,1.08)
start_perturbed$lambda_d <- start_power$lambda_d*c(.85,1.18,.88,1.15,.92)

message(sprintf("Piecewise sample: N=%s; exact cells=%s",
  format(attr(df_fit,"n_original"),big.mark=","),
  format(attr(df_fit,"n_cells"),big.mark=",")))
outdir <- "EM-tenure/output/results/piecewise_hazard"
resume_file <- file.path(outdir,"fits_latest.rds")
refine_only <- identical(tolower(Sys.getenv("PIECEWISE_REFINE_ONLY")),"true") &&
  file.exists(resume_file)
if (refine_only) {
  message("Refining saved piecewise optimum with tight gradient controls")
  saved_piecewise <- readRDS(resume_file)
  piecewise_multi <- fit_eps_piecewise_multistart(df_fit,
    list(saved_piecewise$piecewise$params),maxit=500L,reltol=1e-11,
    pgtol=1e-7,method="BFGS",verbose=0L)
  multistart_report <- rbind(
    transform(saved_piecewise$multistart,phase="preliminary"),
    transform(piecewise_multi$summary,phase="refinement"))
} else {
  message("Estimating five-interval piecewise hazards from three starts")
  preliminary <- fit_eps_piecewise_multistart(df_fit,
    list(start_power,start_constant,start_perturbed),maxit=350L,
    reltol=1e-8,pgtol=1e-5,verbose=0L)
  message("Refining best piecewise optimum with tight gradient controls")
  piecewise_multi <- fit_eps_piecewise_multistart(df_fit,
    list(preliminary$best$params),maxit=500L,reltol=1e-11,
    pgtol=1e-7,method="BFGS",verbose=0L)
  multistart_report <- rbind(
    transform(preliminary$summary,phase="preliminary"),
    transform(piecewise_multi$summary,phase="refinement"))
}
piecewise_fit <- piecewise_multi$best

fits <- list("Constant hazard"=constant_fit,"Power-law hazard"=power_fit,
             "Piecewise hazard"=piecewise_fit)
npar <- c(5L,7L,13L)
model_comparison <- do.call(rbind,lapply(seq_along(fits),function(k) {
  f <- fits[[k]]; p <- f$params
  rates <- duration_weighted_transition_rates(df_fit,f)[1,]
  data.frame(model=names(fits)[k],parameters=npar[k],loglik=f$loglik,
    AIC=-2*f$loglik+2*npar[k],alpha=p$alpha,pi=p$pi,eps=p$eps,
    weighted_exit=rates$exit_rate,weighted_entry=rates$entry_rate,
    mean_employment_months=duration_mean_years(p$lambda_g,p$beta_g)*12,
    mean_nonemployment_months=duration_mean_years(p$lambda_d,p$beta_d)*12,
    median_employment_months=.duration_inverse_cumhaz(log(2),p$lambda_g,p$beta_g)*12,
    median_nonemployment_months=.duration_inverse_cumhaz(log(2),p$lambda_d,p$beta_d)*12)
}))

hazards <- piecewise_hazard_table(piecewise_fit)
representative <- c(.125,.625,2,4,7.5)
hazards$representative_years <- representative
hazards$quarterly_exit <- .duration_transition_probability(
  representative,piecewise_fit$params$lambda_g,0)
hazards$quarterly_entry <- .duration_transition_probability(
  representative,piecewise_fit$params$lambda_d,0)
missing_shares <- duration_missing_risk_shares(df_fit,piecewise_fit)
lr_constant <- data.frame(comparison="Piecewise versus constant",
  LR=2*(piecewise_fit$loglik-constant_fit$loglik),df=8L,
  p_value=pchisq(2*(piecewise_fit$loglik-constant_fit$loglik),8,lower.tail=FALSE))

cat("\nModel comparison\n")
print(model_comparison,row.names=FALSE,digits=8)
cat("\nPiecewise hazards and implied quarterly risks\n")
print(hazards,row.names=FALSE,digits=7)
cat("\nPosterior risk-set shares with unavailable clocks\n")
print(missing_shares[,c("wave","exit_clock_missing","entry_clock_missing")],
  row.names=FALSE,digits=7)
cat("\nMultistart diagnostics\n")
print(multistart_report,row.names=FALSE,digits=8)
cat("\nNested comparison\n")
print(lr_constant,row.names=FALSE,digits=8)
cat("\nRefined scaled gradient by transformed parameter\n")
print(piecewise_fit$gradient,digits=8)

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(model_comparison,file.path(outdir,"model_comparison_latest.csv"),row.names=FALSE)
write.csv(hazards,file.path(outdir,"hazards_latest.csv"),row.names=FALSE)
write.csv(missing_shares,file.path(outdir,"missing_clock_shares_latest.csv"),row.names=FALSE)
write.csv(multistart_report,file.path(outdir,"multistart_latest.csv"),row.names=FALSE)
write.csv(lr_constant,file.path(outdir,"nested_comparison_latest.csv"),row.names=FALSE)
piecewise_fit$objective_function <- NULL
saveRDS(list(piecewise=piecewise_fit,multistart=multistart_report,
  comparison=model_comparison),file.path(outdir,"fits_latest.rds"))
