# Observable predictive checks and likelihood-gain decomposition for the
# preferred whole-year revision model. No parameter re-estimation is done.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
required <- c("dplyr","ggplot2")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")

options(sa_empl_transitions.qlfs_3wave_panel="df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure=TRUE)
source("scripts/ingest_data_3waves_SA.R")
df_full <- prepare_eps_estimation_data(add_nominal_interview_months(df_qlfs),
  allow_zero_tenure=TRUE)
interview_cols <- paste0("interview_month",1:3)
df_fit <- collapse_eps_cells(df_full,allow_zero_tenure=TRUE,
  extra_cols=interview_cols)

fit_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_revision/fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate the calendar-revision model first")
fits <- readRDS(fit_file)
base_params <- fits$base$params
revision_params <- fits$revision$params
base_params$tenure_year_revision_prob <- 0

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_revision/diagnostics")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

.weighted_share <- function(value,weight) {
  keep <- is.finite(value) & is.finite(weight) & weight>0
  if (!any(keep)) return(NA_real_)
  weighted.mean(as.numeric(value[keep]),weight[keep])
}

.calendar_metrics <- function(d) {
  w <- if ("weight" %in% names(d)) d$weight else rep(1,nrow(d))
  ans <- numeric(0)
  pair_values <- list()
  for (pair in list(c(1L,2L),c(2L,3L),c(1L,3L))) {
    a <- pair[1L]; b <- pair[2L]; gap <- b-a
    keep <- d[[paste0("y",a)]]==1L & d[[paste0("y",b)]]==1L
    prior <- round(12*d[[paste0("tenure",a)]][keep])
    current <- round(12*d[[paste0("tenure",b)]][keep])
    delta <- current-prior-gap*.TENURE_MONTHS_PER_WAVE
    exact <- delta==0L
    whole_year <- delta!=0L & delta%%12L==0L
    other <- !exact & !whole_year
    label <- paste0("waves_",a,b,"_")
    ans[paste0(label,"exact_continuation")] <-
      .weighted_share(exact,w[keep])
    ans[paste0(label,"whole_year_revision")] <-
      .weighted_share(whole_year,w[keep])
    ans[paste0(label,"other_revision")] <-
      .weighted_share(other,w[keep])
    if (gap==1L) pair_values[[length(pair_values)+1L]] <-
      list(exact=exact,whole_year=whole_year,other=other,weight=w[keep])
  }
  pooled_weight <- unlist(lapply(pair_values,`[[`,"weight"))
  for (metric in c("exact","whole_year","other"))
    ans[paste0("consecutive_pairs_",metric)] <- .weighted_share(
      unlist(lapply(pair_values,`[[`,metric)),pooled_weight)
  january_values <- january_weights <- list()
  for (wave in 1:3) {
    keep <- d[[paste0("y",wave)]]==1L
    month <- round(12*d[[paste0("tenure",wave)]][keep])
    january <- month%%12L==d[[paste0("interview_month",wave)]][keep]-1L
    ans[paste0("wave_",wave,"_january_start")] <-
      .weighted_share(january,w[keep])
    january_values[[wave]] <- january; january_weights[[wave]] <- w[keep]
  }
  ans["all_employed_reports_january_start"] <- .weighted_share(
    unlist(january_values),unlist(january_weights))
  ans
}

observed <- .calendar_metrics(df_full)
reps <- as.integer(Sys.getenv("REVISION_PREDICTIVE_REPS","40"))
sim_n <- as.integer(Sys.getenv("REVISION_PREDICTIVE_N","50000"))
workers <- as.integer(Sys.getenv("REVISION_PREDICTIVE_WORKERS","1"))
resume <- identical(tolower(Sys.getenv("REVISION_PREDICTIVE_RESUME","true")),
  "true")
sim_file <- file.path(outdir,"predictive_draws_latest.rds")

models <- list(`January heaping`=base_params,
  `January + whole-year revisions`=revision_params)
tasks <- do.call(c,lapply(names(models),function(model)
  lapply(seq_len(reps),function(replication)
    list(model=model,replication=replication,
      seed=260910L+replication+1000L*match(model,names(models))))))
draws <- NULL
if (resume && file.exists(sim_file)) {
  saved <- readRDS(sim_file)
  if (identical(saved$reps,reps) && identical(saved$sim_n,sim_n))
    draws <- saved$draws
}
if (is.null(draws)) {
  task_eval <- function(task) {
    d <- simulate_eps_piecewise_job_change(sim_n,models[[task$model]],
      seed=task$seed)
    metric_value <- .calendar_metrics(d)
    data.frame(model=task$model,replication=task$replication,
      metric=names(metric_value),estimate=unname(metric_value))
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
  saveRDS(list(reps=reps,sim_n=sim_n,draws=draws),sim_file)
}

predictive <- draws |>
  group_by(model,metric) |>
  summarise(predicted_mean=mean(estimate),predicted_sd=sd(estimate),
    predicted_lower=quantile(estimate,.025),
    predicted_upper=quantile(estimate,.975),.groups="drop") |>
  mutate(observed=unname(observed[metric]),
    standardized_difference=(observed-predicted_mean)/predicted_sd,
    within_95=observed>=predicted_lower & observed<=predicted_upper)
write.csv(predictive,file.path(outdir,"posterior_predictive_latest.csv"),
  row.names=FALSE)

# Decompose the exact observed-data likelihood improvement by mutually
# exclusive observable calendar patterns on the collapsed likelihood cells.
base_eval <- e_step_eps(df_fit,base_params,check_df=FALSE,suff_stats=FALSE)
revision_eval <- e_step_eps(df_fit,revision_params,check_df=FALSE,
  suff_stats=FALSE)
gain <- df_fit$weight*(revision_eval$row_loglik-base_eval$row_loglik)
if (abs(sum(gain)-(revision_eval$loglik-base_eval$loglik))>1e-5)
  stop("Row-level likelihood decomposition does not reproduce total gain")

reported_employed <- as.matrix(df_fit[paste0("y",1:3)])==1L
tenure_month <- round(12*as.matrix(df_fit[paste0("tenure",1:3)]))
pair_exact <- pair_year <- pair_other <- matrix(FALSE,nrow(df_fit),3L)
pair_names <- c("Waves 1-2","Waves 2-3","Waves 1-3")
pair_index <- list(c(1L,2L),c(2L,3L),c(1L,3L))
for (j in seq_along(pair_index)) {
  a <- pair_index[[j]][1L]; b <- pair_index[[j]][2L]; gap <- b-a
  eligible <- reported_employed[,a] & reported_employed[,b]
  delta <- tenure_month[,b]-tenure_month[,a]-
    gap*.TENURE_MONTHS_PER_WAVE
  pair_exact[,j] <- eligible & delta==0L
  pair_year[,j] <- eligible & delta!=0L & delta%%12L==0L
  pair_other[,j] <- eligible & delta!=0L & delta%%12L!=0L
}
pattern <- ifelse(rowSums(pair_year)>0,"At least one whole-year revision",
  ifelse(rowSums(pair_exact)>0,"At least one exact continuation",
    ifelse(rowSums(pair_other)>0,"Other multiple employed reports",
      "Fewer than two employed reports")))
decomposition <- data.frame(pattern=pattern,weight=df_fit$weight,gain=gain) |>
  group_by(pattern) |>
  summarise(population_weight=sum(weight),loglik_gain=sum(gain),
    mean_gain=loglik_gain/population_weight,.groups="drop") |>
  mutate(population_share=population_weight/sum(population_weight),
    gain_share=loglik_gain/sum(loglik_gain))
write.csv(decomposition,file.path(outdir,
  "likelihood_gain_by_pattern_latest.csv"),row.names=FALSE)

pair_decomposition <- do.call(rbind,lapply(seq_along(pair_index),function(j) {
  category <- ifelse(pair_year[,j],"Whole-year revision",
    ifelse(pair_exact[,j],"Exact continuation",
      ifelse(pair_other[,j],"Other revision","Not jointly employed")))
  data.frame(pair=pair_names[j],category=category,weight=df_fit$weight,
    gain=gain) |>
    group_by(pair,category) |>
    summarise(population_weight=sum(weight),loglik_gain=sum(gain),
      .groups="drop")
})) |>
  group_by(pair) |>
  mutate(population_share=population_weight/sum(population_weight),
    gain_share=loglik_gain/sum(loglik_gain)) |>
  ungroup()
write.csv(pair_decomposition,file.path(outdir,
  "likelihood_gain_by_pair_latest.csv"),row.names=FALSE)

cat("\nObservable posterior-predictive checks\n")
print(predictive |> filter(grepl("consecutive_pairs|all_employed",metric)),
  n=Inf,width=Inf)
cat("\nLikelihood-gain decomposition\n")
print(decomposition,n=Inf,width=Inf)
