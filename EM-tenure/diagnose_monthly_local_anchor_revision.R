# Predictive checks, likelihood-gain decomposition, and conditional curvature
# for the short-range reported-start-date revision extension. No parameters are
# fitted here.

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
  "calendar_local_anchor_revision/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate local anchor revisions first")
saved <- readRDS(fit_file)
models <- list(`Exact-anchor model`=saved$base$params,
  `Local anchor revisions`=saved$extension$params)
outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_local_anchor_revision/diagnostics")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

.weighted_share <- function(value,weight) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  if (!any(keep)) return(NA_real_)
  weighted.mean(as.numeric(value[keep]),weight[keep])
}
.calendar_metrics <- function(d) {
  w <- if ("weight" %in% names(d)) d$weight else rep(1,nrow(d))
  exact <- local <- yearly <- other <- weights <- list()
  for (j in 1:2) {
    keep <- d[[paste0("y",j)]]==1L & d[[paste0("y",j+1L)]]==1L
    delta <- round(12*d[[paste0("tenure",j+1L)]][keep])-
      round(12*d[[paste0("tenure",j)]][keep])-3L
    exact[[j]] <- delta==0L
    local[[j]] <- delta!=0L & abs(delta)<=6L
    yearly[[j]] <- delta!=0L & delta%%12L==0L
    other[[j]] <- !exact[[j]] & !local[[j]] & !yearly[[j]]
    weights[[j]] <- w[keep]
  }
  pooled_weight <- unlist(weights)
  c(exact_continuation=.weighted_share(unlist(exact),pooled_weight),
    local_revision_1_to_6_months=
      .weighted_share(unlist(local),pooled_weight),
    whole_year_revision=.weighted_share(unlist(yearly),pooled_weight),
    other_revision=.weighted_share(unlist(other),pooled_weight))
}

observed <- .calendar_metrics(df_full)
reps <- as.integer(Sys.getenv("LOCAL_ANCHOR_PREDICTIVE_REPS","20"))
sim_n <- as.integer(Sys.getenv("LOCAL_ANCHOR_PREDICTIVE_N","50000"))
workers <- as.integer(Sys.getenv("LOCAL_ANCHOR_PREDICTIVE_WORKERS","4"))
tasks <- do.call(c,lapply(names(models),function(model)
  lapply(seq_len(reps),function(replication)
    list(model=model,replication=replication,
      seed=260902L+replication+1000L*match(model,names(models))))))
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
exact_count <- local_count <- yearly_count <- other_count <-
  integer(nrow(df_fit))
for (j in 1:2) {
  keep <- reported[,j] & reported[,j+1L]
  delta <- tenure[,j+1L]-tenure[,j]-3L
  exact_count <- exact_count+as.integer(keep & delta==0L)
  local_count <- local_count+as.integer(keep & delta!=0L & abs(delta)<=6L)
  yearly_count <- yearly_count+as.integer(keep & delta!=0L &
    delta%%12L==0L)
  other_count <- other_count+as.integer(keep & delta!=0L & abs(delta)>6L &
    delta%%12L!=0L)
}
total_count <- exact_count+local_count+yearly_count+other_count
nonzero_types <- (exact_count>0)+(local_count>0)+(yearly_count>0)+
  (other_count>0)
pattern <- ifelse(total_count==0L,"No consecutive employed reports",
  ifelse(nonzero_types>1L,"Mixed",
    ifelse(exact_count>0L,"Exact continuation only",
      ifelse(local_count>0L,"Local revision only",
        ifelse(yearly_count>0L,"Whole-year revision only",
          "Other revision only")))))
decomposition <- data.frame(pattern=pattern,weight=df_fit$weight,gain=gain) |>
  group_by(pattern) |>
  summarise(population_weight=sum(weight),loglik_gain=sum(gain),
    mean_gain=loglik_gain/population_weight,.groups="drop") |>
  mutate(population_share=population_weight/sum(population_weight),
    gain_share=loglik_gain/sum(loglik_gain))

# Conditional observed-curvature approximation for kappa. This is deliberately
# labelled conditional because the remaining parameters are held at their
# retained blockwise estimates.
p <- saved$extension$params
z0 <- qlogis(p$tenure_local_revision_prob)
h <- as.numeric(Sys.getenv("LOCAL_ANCHOR_CURVATURE_STEP","0.005"))
loglik_at_z <- function(z) {
  value <- p
  value$tenure_local_revision_prob <- plogis(z)
  e_step_eps(df_fit,value,check_df=FALSE,suff_stats=FALSE)$loglik
}
ll <- c(minus=loglik_at_z(z0-h),centre=loglik_at_z(z0),
  plus=loglik_at_z(z0+h))
information_z <- -(ll["plus"]-2*ll["centre"]+ll["minus"])/h^2
se_z <- if (is.finite(information_z) && information_z>0)
  sqrt(1/information_z) else NA_real_
curvature <- data.frame(estimate=p$tenure_local_revision_prob,
  logit_estimate=z0,step=h,information_logit=information_z,se_logit=se_z,
  se_probability=p$tenure_local_revision_prob*
    (1-p$tenure_local_revision_prob)*se_z,
  lower_probability=plogis(z0-qnorm(.975)*se_z),
  upper_probability=plogis(z0+qnorm(.975)*se_z),
  inference="Conditional curvature; other parameters fixed")

write.csv(draws,file.path(outdir,"predictive_draws_latest.csv"),row.names=FALSE)
write.csv(predictive,file.path(outdir,"posterior_predictive_latest.csv"),
  row.names=FALSE)
write.csv(decomposition,file.path(outdir,
  "likelihood_gain_by_continuation_pattern_latest.csv"),row.names=FALSE)
write.csv(curvature,file.path(outdir,"conditional_curvature_latest.csv"),
  row.names=FALSE)
cat("\nLocal-anchor predictive checks\n")
print(predictive,n=Inf,width=Inf)
cat("\nLikelihood-gain decomposition\n")
print(decomposition,n=Inf,width=Inf)
cat("\nConditional curvature for local revision probability\n")
print(curvature,row.names=FALSE,digits=10)
