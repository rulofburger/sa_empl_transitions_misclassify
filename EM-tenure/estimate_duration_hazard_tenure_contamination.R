# Compare the unchanged constant-hazard linked epsilon model with a nested
# log-duration-hazard extension. Point estimates only; no bootstrap.
if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from the project root")
required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_3waves_SA.R")

df_fit <- collapse_eps_cells(prepare_eps_estimation_data(df_qlfs))
message(sprintf("Duration-hazard sample: N=%s; exact cells=%s",
  format(attr(df_fit,"n_original"),big.mark=","),
  format(attr(df_fit,"n_cells"),big.mark=",")))

message("Re-estimating unchanged nonstationary CTMC-linked model")
constant_em <- em_multistart_eps(df_fit, stationary=FALSE, linked=TRUE,
  K=3L, max_iter=1000L, tol=1e-10, param_tol=1e-7, verbose=0L)$best

message("Validating direct duration likelihood at beta_g=beta_d=0")
constant_direct <- fit_eps_duration_hazard(
  df_fit, constant_em$params, beta_fixed=0, maxit=250L, reltol=1e-11)

message("Estimating log-duration hazards from three starts")
duration_multi <- fit_eps_duration_multistart(
  df_fit, constant_direct$params, maxit=500L, reltol=1e-10)
duration_fit <- duration_multi$best
p0 <- constant_direct$params; p1 <- duration_fit$params

summary_table <- data.frame(
  model=c("Constant hazard","Log-duration hazard"),
  loglik=c(constant_direct$loglik,duration_fit$loglik),
  alpha=c(p0$alpha,p1$alpha), pi=c(p0$pi,p1$pi), eps=c(p0$eps,p1$eps),
  lambda_g=c(p0$lambda_g,p1$lambda_g), beta_g=c(0,p1$beta_g),
  lambda_d=c(p0$lambda_d,p1$lambda_d), beta_d=c(0,p1$beta_d),
  mean_employment_months=c(duration_mean_years(p0$lambda_g,0),
                            duration_mean_years(p1$lambda_g,p1$beta_g))*12,
  mean_nonemployment_months=c(duration_mean_years(p0$lambda_d,0),
                               duration_mean_years(p1$lambda_d,p1$beta_d))*12)
profile <- duration_transition_profile(p1)
weighted_rates <- rbind(
  cbind(model="Constant hazard",
        duration_weighted_transition_rates(df_fit,constant_direct)),
  cbind(model="Log-duration hazard",
        duration_weighted_transition_rates(df_fit,duration_fit)))
missing_clock_shares <- duration_missing_risk_shares(df_fit,duration_fit)
lr <- 2*(duration_fit$loglik-constant_direct$loglik)
lr_table <- data.frame(LR=lr,df=2,p_value=pchisq(lr,df=2,lower.tail=FALSE),
  constant_EM_LL=constant_em$loglik,
  constant_direct_LL=constant_direct$loglik,
  nested_LL_difference=constant_direct$loglik-constant_em$loglik,
  optimizer_code=duration_fit$convergence,
  scaled_gradient_max=duration_fit$gradient_max)

cat("\nConstant versus log-duration hazards\n")
print(summary_table,row.names=FALSE,digits=8)
cat("\nImplied quarterly transition probabilities\n")
print(profile,row.names=FALSE,digits=7)
cat("\nPosterior-risk and survey-weighted transition probabilities\n")
print(weighted_rates[,c("model","wave","exit_rate","entry_rate")],
      row.names=FALSE,digits=7)
cat("\nPosterior risk-set shares with unavailable origin clocks\n")
print(missing_clock_shares[,c("wave","exit_clock_missing","entry_clock_missing")],
      row.names=FALSE,digits=7)
cat("\nDuration-dependent multi-start diagnostics\n")
print(duration_multi$summary,row.names=FALSE,digits=8)
cat("\nNested-model and optimizer diagnostics\n")
print(lr_table,row.names=FALSE,digits=8)

outdir <- "EM-tenure/output/results/duration_hazard"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(summary_table,file.path(outdir,"model_comparison_latest.csv"),row.names=FALSE)
write.csv(profile,file.path(outdir,"transition_profile_latest.csv"),row.names=FALSE)
write.csv(weighted_rates,file.path(outdir,"weighted_rates_latest.csv"),row.names=FALSE)
write.csv(missing_clock_shares,file.path(outdir,"missing_clock_shares_latest.csv"),row.names=FALSE)
write.csv(lr_table,file.path(outdir,"diagnostics_latest.csv"),row.names=FALSE)
constant_direct$objective_function <- NULL
duration_fit$objective_function <- NULL
saveRDS(list(constant_em=constant_em,constant_direct=constant_direct,
             duration=duration_fit,multistart=duration_multi$summary),
        file.path(outdir,"fits_latest.rds"))
