# Predictive and likelihood-gain diagnostics for the shared person-level
# duration-reporting reliability mixture. No parameters are estimated here.

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
  "shared_duration_reliability/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate shared duration reliability first")
saved <- readRDS(fit_file)
models <- list(`Independent reliability`=saved$base$params,
  `Shared reliability`=saved$extension$params)
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "shared_duration_reliability/diagnostics")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

.weighted_share <- function(value,weight) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  if (!any(keep)) return(NA_real_)
  weighted.mean(as.numeric(value[keep]),weight[keep])
}
.duration_patterns <- function(d,params) {
  n <- nrow(d)
  tenure_deviation <- timegap_inconsistent <- opportunities <- integer(n)
  exact <- local <- yearly <- other <- logical()
  calendar_weight <- numeric()
  timegap_bad <- logical(); timegap_weight <- numeric()
  w <- if ("weight" %in% names(d)) d$weight else rep(1,n)
  for (j in 1:2) {
    employed <- d[[paste0("y",j)]]==1L &
      d[[paste0("y",j+1L)]]==1L
    delta <- round(12*d[[paste0("tenure",j+1L)]])-
      round(12*d[[paste0("tenure",j)]])-3L
    tenure_deviation <- tenure_deviation+as.integer(employed & delta!=0L)
    opportunities <- opportunities+as.integer(employed)
    exact <- c(exact,delta[employed]==0L)
    local <- c(local,delta[employed]!=0L & abs(delta[employed])<=6L)
    yearly <- c(yearly,delta[employed]!=0L & delta[employed]%%12L==0L)
    other <- c(other,delta[employed]!=0L & abs(delta[employed])>6L &
      delta[employed]%%12L!=0L)
    calendar_weight <- c(calendar_weight,w[employed])

    nonemployed <- d[[paste0("y",j)]]==0L &
      d[[paste0("y",j+1L)]]==0L
    clean_lp <- log_emission_transition_d(
      d[[paste0("timegap_cat",j+1L)]][nonemployed],
      d[[paste0("timegap_cat",j)]][nonemployed],params$lambda_d,
      if (is.null(params$beta_d)) 0 else params$beta_d)
    bad <- !is.finite(clean_lp)
    timegap_inconsistent[nonemployed] <-
      timegap_inconsistent[nonemployed]+as.integer(bad)
    opportunities <- opportunities+as.integer(nonemployed)
    timegap_bad <- c(timegap_bad,bad)
    timegap_weight <- c(timegap_weight,w[nonemployed])
  }
  irregular <- tenure_deviation+timegap_inconsistent
  two <- opportunities==2L
  metrics <- c(exact_continuation=.weighted_share(exact,calendar_weight),
    local_revision_1_to_6_months=.weighted_share(local,calendar_weight),
    whole_year_revision=.weighted_share(yearly,calendar_weight),
    other_revision=.weighted_share(other,calendar_weight),
    timegap_inconsistent=.weighted_share(timegap_bad,timegap_weight),
    two_opportunities_no_irregularity=.weighted_share(two & irregular==0L,
      w)/.weighted_share(two,w),
    two_opportunities_one_irregularity=.weighted_share(two & irregular==1L,
      w)/.weighted_share(two,w),
    two_opportunities_two_irregularities=.weighted_share(two & irregular==2L,
      w)/.weighted_share(two,w))
  pattern <- ifelse(opportunities==0L,"No repeated-state duration",
    ifelse(tenure_deviation>0L,"Tenure deviation",
      ifelse(timegap_inconsistent>0L,"Timegap inconsistency",
        "Repeated state; no irregularity")))
  list(metrics=metrics,pattern=pattern,opportunities=opportunities,
    irregular=irregular,tenure_deviation=tenure_deviation,
    timegap_inconsistent=timegap_inconsistent)
}

observed_patterns <- .duration_patterns(df_full,models[[2]])
observed <- observed_patterns$metrics
reps <- as.integer(Sys.getenv("SHARED_RELIABILITY_PREDICTIVE_REPS","20"))
sim_n <- as.integer(Sys.getenv("SHARED_RELIABILITY_PREDICTIVE_N","50000"))
workers <- as.integer(Sys.getenv("SHARED_RELIABILITY_PREDICTIVE_WORKERS","4"))
tasks <- do.call(c,lapply(names(models),function(model)
  lapply(seq_len(reps),function(replication)
    list(model=model,replication=replication,
      seed=260903L+replication+1000L*match(model,names(models))))))
task_eval <- function(task) {
  p <- models[[task$model]]
  d <- simulate_eps_piecewise_job_change(sim_n,p,seed=task$seed)
  value <- .duration_patterns(d,p)$metrics
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
  parallel::clusterExport(cluster,c("sim_n","models",".duration_patterns",
    ".weighted_share"),envir=environment())
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

old_eval <- e_step_eps(df_fit,models[[1]],check_df=FALSE,suff_stats=FALSE)
new_eval <- e_step_eps(df_fit,models[[2]],check_df=FALSE,suff_stats=FALSE)
fit_patterns <- .duration_patterns(df_fit,models[[2]])
gain <- df_fit$weight*(new_eval$row_loglik-old_eval$row_loglik)
posterior <- new_eval$duration_reliability_posterior
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
cat("\nShared-reliability predictive checks\n")
print(predictive,n=Inf,width=Inf)
cat("\nLikelihood gain and posterior class by observed duration pattern\n")
print(decomposition,n=Inf,width=Inf)
cat("\nIdentification note\n")
cat(paste("With three waves, one record cannot contain both a repeated",
  "employment duration and a repeated nonemployment duration. The common",
  "cross-clock shift is therefore a parsimonious restriction supported",
  "indirectly by within-clock serial clustering, not a directly observed",
  "within-person tenure-timegap correlation.\n"))
