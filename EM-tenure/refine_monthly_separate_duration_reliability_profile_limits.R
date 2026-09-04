# Add one tightly optimized point near each provisional profile-likelihood
# confidence limit, then recompute the interpolated crossings. The outer
# points have already passed two-start screening in the main profile script.

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

root <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability")
saved <- readRDS(file.path(root,"fits_latest.rds"))
unrestricted <- saved$extension
outdir <- file.path(root,"profile_difference")
profile_file <- file.path(outdir,"profile_latest.rds")
if (!file.exists(profile_file)) stop("Run the main difference profile first")
profile <- readRDS(profile_file)
workers <- max(1L,min(8L,as.integer(Sys.getenv(
  "RELIABILITY_PROFILE_WORKERS","8"))))
maxit <- as.integer(Sys.getenv("RELIABILITY_PROFILE_LIMIT_MAXIT","2"))
resume <- identical(tolower(Sys.getenv("RELIABILITY_PROFILE_RESUME","true")),
  "true")
threshold <- qchisq(.95,1)
d_hat <- with(unrestricted$params,timegap_reliability_shift-
  tenure_reliability_shift)

outer <- profile$profile |>
  filter(abs(difference-d_hat)>1e-10,LR>=threshold) |>
  mutate(side=ifelse(difference<d_hat,"lower","upper")) |>
  group_by(side) |>
  slice_min(abs(difference-d_hat),n=1,with_ties=FALSE) |>
  ungroup()
if (!all(c("lower","upper") %in% outer$side))
  stop("The existing profile does not bracket both limits")
targets <- outer |>
  transmute(side,difference=d_hat+(difference-d_hat)*sqrt(threshold/LR))

.safe_name <- function(x) {
  prefix <- if (x<0) "m" else "p"
  paste0(prefix,gsub("\\.","_",formatC(abs(x),digits=9,format="f")))
}
.compact_fit <- function(fit) {
  fit$gamma <- NULL; fit$job_change_posterior <- NULL
  fit$objective_function <- NULL; fit
}
.fit_target <- function(side,difference) {
  path <- file.path(outdir,paste0("profile_",.safe_name(difference),"_",
    side,"_limit.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Refining %s profile limit at d=%+.9f",side,difference))
  p <- unrestricted$params
  fit <- fit_eps_piecewise_calendar_revision_monthly_fixed_difference(
    df_fit,p,difference=difference,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    q_start=p$job_change_prob,maxit=maxit,reltol=1e-12,pgtol=1e-10,
    workers=workers,verbose=1L,gradient_step=5e-5,
    gradient_scheme="forward")
  fit <- .compact_fit(fit)
  fit$profile_start <- paste0(side,"_limit")
  saveRDS(fit,path)
  fit
}

limit_fits <- lapply(seq_len(nrow(targets)),function(i)
  .fit_target(targets$side[i],targets$difference[i]))
limit_rows <- do.call(rbind,lapply(limit_fits,function(fit) data.frame(
  difference=fit$profile_difference,loglik=fit$loglik,
  convergence=fit$convergence,iterations=fit$iterations,
  tenure_shift=fit$params$tenure_reliability_shift,
  timegap_shift=fit$params$timegap_reliability_shift,
  relative_loglik=fit$loglik-unrestricted$loglik,
  LR=max(0,2*(unrestricted$loglik-fit$loglik)))))

combined <- bind_rows(profile$profile,limit_rows) |>
  arrange(difference) |>
  distinct(difference,.keep_all=TRUE)
.profile_crossing <- function(x,y,target=threshold) {
  hit <- which(diff(sign(y-target))!=0)
  if (!length(hit)) return(NA_real_)
  j <- hit[1L]
  x[j]+(target-y[j])*(x[j+1L]-x[j])/(y[j+1L]-y[j])
}
left <- combined[combined$difference<=d_hat,] |> arrange(desc(difference))
right <- combined[combined$difference>=d_hat,] |> arrange(difference)
interval <- data.frame(estimate=d_hat,
  profile_lower=.profile_crossing(left$difference,left$LR),
  profile_upper=.profile_crossing(right$difference,right$LR),
  threshold=threshold,
  equality_LR=2*(unrestricted$loglik-saved$common$loglik),
  equality_p_value=pchisq(2*(unrestricted$loglik-saved$common$loglik),1,
    lower.tail=FALSE))
profile$profile <- combined
profile$interval <- interval
profile$limit_points <- limit_rows
for (i in seq_along(limit_fits))
  profile$fits[[.safe_name(limit_fits[[i]]$profile_difference)]] <-
    limit_fits[[i]]
saveRDS(profile,profile_file)
write.csv(combined,file.path(outdir,"profile_likelihood_latest.csv"),
  row.names=FALSE)
write.csv(interval,file.path(outdir,"profile_interval_latest.csv"),
  row.names=FALSE)
write.csv(limit_rows,file.path(outdir,"profile_limit_checks_latest.csv"),
  row.names=FALSE)

cat("\nTargeted profile-limit checks\n")
print(limit_rows,row.names=FALSE,digits=12)
cat("\nRefined nuisance-adjusted profile interval\n")
print(interval,row.names=FALSE,digits=12)
if (any(limit_rows$convergence!=0L))
  warning("At least one tightened limit-point optimization did not converge")
