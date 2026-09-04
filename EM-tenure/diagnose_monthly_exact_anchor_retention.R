# Predictive checks and likelihood-gain decomposition for the exact reported
# start-date retention extension. No parameters are fitted here.

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
  "calendar_exact_anchor_retention/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate exact-anchor retention first")
saved <- readRDS(fit_file)
models <- list(`Flexible baseline`=saved$base$params,
  `Exact-anchor retention`=saved$extension$params)
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention/diagnostics")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

.weighted_share <- function(value,weight) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  if (!any(keep)) return(NA_real_)
  weighted.mean(as.numeric(value[keep]),weight[keep])
}
.calendar_metrics <- function(d) {
  w <- if ("weight" %in% names(d)) d$weight else rep(1,nrow(d))
  exact <- yearly <- other <- weights <- list()
  for (j in 1:2) {
    keep <- d[[paste0("y",j)]]==1L & d[[paste0("y",j+1L)]]==1L
    delta <- round(12*d[[paste0("tenure",j+1L)]][keep])-
      round(12*d[[paste0("tenure",j)]][keep])-3L
    exact[[j]] <- delta==0L
    yearly[[j]] <- delta!=0L & delta%%12L==0L
    other[[j]] <- !exact[[j]] & !yearly[[j]]
    weights[[j]] <- w[keep]
  }
  pooled_weight <- unlist(weights)
  c(exact_continuation=.weighted_share(unlist(exact),pooled_weight),
    whole_year_revision=.weighted_share(unlist(yearly),pooled_weight),
    other_revision=.weighted_share(unlist(other),pooled_weight))
}

observed <- .calendar_metrics(df_full)
reps <- as.integer(Sys.getenv("EXACT_ANCHOR_PREDICTIVE_REPS","30"))
sim_n <- as.integer(Sys.getenv("EXACT_ANCHOR_PREDICTIVE_N","50000"))
workers <- as.integer(Sys.getenv("EXACT_ANCHOR_PREDICTIVE_WORKERS","1"))
tasks <- do.call(c,lapply(names(models),function(model)
  lapply(seq_len(reps),function(replication)
    list(model=model,replication=replication,
      seed=260916L+replication+1000L*match(model,names(models))))))
task_eval <- function(task) {
  d <- simulate_eps_piecewise_job_change(sim_n,models[[task$model]],
    seed=task$seed)
  value <- .calendar_metrics(d)
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
  parallel::clusterExport(cluster,c("sim_n","models",".calendar_metrics",
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
gain <- df_fit$weight*(new_eval$row_loglik-old_eval$row_loglik)
reported <- as.matrix(df_fit[paste0("y",1:3)])==1L
tenure <- round(12*as.matrix(df_fit[paste0("tenure",1:3)]))
exact_count <- yearly_count <- other_count <- integer(nrow(df_fit))
for (j in 1:2) {
  keep <- reported[,j] & reported[,j+1L]
  delta <- tenure[,j+1L]-tenure[,j]-3L
  exact_count <- exact_count+as.integer(keep & delta==0L)
  yearly_count <- yearly_count+as.integer(keep & delta!=0L & delta%%12L==0L)
  other_count <- other_count+as.integer(keep & delta!=0L & delta%%12L!=0L)
}
pattern <- ifelse(exact_count+yearly_count+other_count==0L,
  "No consecutive employed reports",ifelse(exact_count>0L &
  yearly_count+other_count==0L,"Exact continuation only",
  ifelse(exact_count==0L & yearly_count>0L & other_count==0L,
    "Whole-year revision only",ifelse(exact_count==0L &
    yearly_count==0L & other_count>0L,"Other revision only","Mixed"))))
decomposition <- data.frame(pattern=pattern,weight=df_fit$weight,gain=gain) |>
  group_by(pattern) |>
  summarise(population_weight=sum(weight),loglik_gain=sum(gain),
    mean_gain=loglik_gain/population_weight,.groups="drop") |>
  mutate(population_share=population_weight/sum(population_weight),
    gain_share=loglik_gain/sum(loglik_gain))

write.csv(draws,file.path(outdir,"predictive_draws_latest.csv"),row.names=FALSE)
write.csv(predictive,file.path(outdir,"posterior_predictive_latest.csv"),
  row.names=FALSE)
write.csv(decomposition,file.path(outdir,
  "likelihood_gain_by_continuation_pattern_latest.csv"),row.names=FALSE)
cat("\nExact-anchor predictive checks\n")
print(predictive,n=Inf,width=Inf)
cat("\nLikelihood-gain decomposition\n")
print(decomposition,n=Inf,width=Inf)
