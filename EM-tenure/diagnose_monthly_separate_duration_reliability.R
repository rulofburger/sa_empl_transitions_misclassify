# Predictive and likelihood-gain diagnostics for separate tenure and timegap
# reliability dispersions. No parameters are estimated here.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=paste0("interview_month",1:3))

fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate separate duration reliability first")
saved <- readRDS(fit_file)
models <- list(`Common dispersion`=saved$common$params,
  `Separate dispersions`=saved$extension$params)
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/diagnostics")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

observed_patterns <- .duration_reporting_patterns(df_full,models[[2]])
observed <- observed_patterns$metrics
reps <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_PREDICTIVE_REPS","20"))
sim_n <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_PREDICTIVE_N","50000"))
workers <- as.integer(Sys.getenv("SEPARATE_RELIABILITY_PREDICTIVE_WORKERS","4"))
tasks <- do.call(c,lapply(names(models),function(model)
  lapply(seq_len(reps),function(replication)
    list(model=model,replication=replication,
      seed=260904L+replication+1000L*match(model,names(models))))))
task_eval <- function(task) {
  p <- models[[task$model]]
  d <- simulate_eps_piecewise_job_change(sim_n,p,seed=task$seed)
  value <- .duration_reporting_patterns(d,p)$metrics
  data.frame(model=task$model,replication=task$replication,
    metric=names(value),estimate=unname(value))
}
if (workers>1L) {
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster),add=TRUE)
  worker_wd <- getwd()
  parallel::clusterCall(cluster,function(path) {
    setwd(path); source("EM-tenure/R/source_all.R"); NULL
  },worker_wd)
  parallel::clusterExport(cluster,c("sim_n","models"),envir=environment())
  draws <- do.call(rbind,parallel::parLapplyLB(cluster,tasks,task_eval))
} else draws <- do.call(rbind,lapply(tasks,task_eval))
predictive <- draws |>
  group_by(model,metric) |>
  summarise(predicted_mean=mean(estimate),predicted_sd=sd(estimate),
    predicted_lower=quantile(estimate,.025),
    predicted_upper=quantile(estimate,.975),.groups="drop") |>
  mutate(observed=unname(observed[metric]),
    standardized_difference=(observed-predicted_mean)/predicted_sd,
    within_95=observed>=predicted_lower & observed<=predicted_upper)

common_eval <- e_step_eps(df_fit,models[[1]],check_df=FALSE,suff_stats=FALSE)
separate_eval <- e_step_eps(df_fit,models[[2]],check_df=FALSE,suff_stats=FALSE)
fit_patterns <- .duration_reporting_patterns(df_fit,models[[2]])
gain <- df_fit$weight*(separate_eval$row_loglik-common_eval$row_loglik)
posterior <- separate_eval$duration_reliability_posterior
decomposition <- data.frame(pattern=fit_patterns$pattern,
  weight=df_fit$weight,gain=gain,posterior=posterior) |>
  group_by(pattern) |>
  summarise(population_weight=sum(weight),loglik_gain=sum(gain),
    mean_gain=loglik_gain/population_weight,
    mean_posterior_unreliable=weighted.mean(posterior,weight),.groups="drop") |>
  mutate(population_share=population_weight/sum(population_weight),
    gain_share=loglik_gain/sum(loglik_gain))

write.csv(draws,file.path(outdir,"predictive_draws_latest.csv"),row.names=FALSE)
write.csv(predictive,file.path(outdir,"posterior_predictive_latest.csv"),
  row.names=FALSE)
write.csv(decomposition,file.path(outdir,
  "likelihood_gain_by_duration_pattern_latest.csv"),row.names=FALSE)
cat("\nSeparate-reliability predictive checks\n")
print(predictive,n=Inf,width=Inf)
cat("\nLikelihood gain and posterior class by observed duration pattern\n")
print(decomposition,n=Inf,width=Inf)
