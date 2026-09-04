# Nuisance-adjusted profile likelihood for the difference between the
# timegap and tenure reliability dispersions in the preferred monthly model.
# Each fixed-difference point is screened from two nuisance starts and then
# continued from the better solution. Completed fits are cached.

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
fit_file <- file.path(root,"fits_latest.rds")
if (!file.exists(fit_file))
  stop("Estimate and polish the separate-reliability model first")
saved <- readRDS(fit_file)
unrestricted <- saved$extension
common <- saved$common
if (is.null(unrestricted$params$tenure_reliability_shift) ||
    is.null(unrestricted$params$timegap_reliability_shift))
  stop("The saved preferred fit does not contain separate reliability shifts")

outdir <- file.path(root,"profile_difference")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- max(1L,min(8L,as.integer(Sys.getenv(
  "RELIABILITY_PROFILE_WORKERS","8"))))
screen_maxit <- as.integer(Sys.getenv("RELIABILITY_PROFILE_SCREEN_MAXIT","1"))
refine_maxit <- as.integer(Sys.getenv("RELIABILITY_PROFILE_REFINE_MAXIT","5"))
refine_chunks <- as.integer(Sys.getenv("RELIABILITY_PROFILE_REFINE_CHUNKS","3"))
resume <- identical(tolower(Sys.getenv("RELIABILITY_PROFILE_RESUME","true")),
  "true")

d_hat <- with(unrestricted$params,timegap_reliability_shift-
  tenure_reliability_shift)
offsets <- c(-.04,-.025,-.015,0,.015,.025,.04)
if (nzchar(Sys.getenv("RELIABILITY_PROFILE_OFFSETS")))
  offsets <- as.numeric(strsplit(Sys.getenv("RELIABILITY_PROFILE_OFFSETS"),
    ",",fixed=TRUE)[[1L]])
grid <- sort(unique(c(0,d_hat+offsets)))
write.csv(data.frame(difference=grid,offset_from_mle=grid-d_hat),
  file.path(outdir,"profile_grid_latest.csv"),row.names=FALSE)

.safe_name <- function(x) {
  prefix <- if (x<0) "m" else "p"
  paste0(prefix,gsub("\\.","_",formatC(abs(x),digits=9,format="f")))
}
.compact_fit <- function(fit) {
  fit$gamma <- NULL
  fit$job_change_posterior <- NULL
  fit$objective_function <- NULL
  fit
}
.fit_point <- function(difference,start,label,maxit) {
  path <- file.path(outdir,paste0("profile_",.safe_name(difference),"_",
    label,".rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Profile d=%+.9f, start=%s, maxit=%d",
    difference,label,maxit))
  p <- start$params
  fit <- fit_eps_piecewise_calendar_revision_monthly_fixed_difference(
    df_fit,p,difference=difference,
    heaping_start=max(p$tenure_heaping_prob,1e-8),
    revision_start=p$tenure_year_revision_prob,
    clean_anchor_revision_start=p$tenure_clean_anchor_revision_prob,
    start_month_probs_start=p$tenure_start_month_probs,
    exact_anchor_retention_start=p$tenure_exact_anchor_retention_prob,
    local_revision_start=p$tenure_local_revision_prob,
    q_start=p$job_change_prob,maxit=maxit,reltol=1e-11,pgtol=1e-8,
    workers=workers,verbose=1L,gradient_step=5e-5,
    gradient_scheme="forward")
  fit <- .compact_fit(fit)
  fit$profile_start <- label
  saveRDS(fit,path)
  fit
}

# Exact central and equality-restriction fits anchor the profile. The equality
# point is algebraically identical to the common-dispersion model.
exact <- list()
exact[[.safe_name(d_hat)]] <- .compact_fit(unrestricted)
common_separate <- common
common_separate$params$tenure_reliability_shift <-
  common$params$duration_reliability_shift
common_separate$params$timegap_reliability_shift <-
  common$params$duration_reliability_shift
common_separate$params$duration_reliability_shift <- 0
exact[[.safe_name(0)]] <- .compact_fit(common_separate)

all_starts <- list()
selected <- exact
.run_direction <- function(values,direction_label) {
  previous <- unrestricted
  output <- list()
  for (difference in values) {
    key <- .safe_name(difference)
    if (key %in% names(exact)) {
      previous <- exact[[key]]
      next
    }
    # The warm continuation and common-dispersion nuisance estimates are
    # genuinely distinct starts except at the equality point itself.
    warm <- .fit_point(difference,previous,
      paste0(direction_label,"_warm_screen"),screen_maxit)
    alternative <- .fit_point(difference,common,
      paste0(direction_label,"_common_screen"),screen_maxit)
    candidates <- list(warm=warm,common=alternative)
    current <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),
      "loglik"))]]
    all_starts[[paste0(key,"_",direction_label,"_warm")]] <<- warm
    all_starts[[paste0(key,"_",direction_label,"_common")]] <<- alternative
    if (!identical(current$convergence,0L)) {
      for (chunk in seq_len(refine_chunks)) {
        refined <- .fit_point(difference,current,paste0(direction_label,
          "_refine_",chunk),refine_maxit)
        all_starts[[paste0(key,"_",direction_label,"_refine_",chunk)]] <<-
          refined
        gain <- refined$loglik-current$loglik
        if (refined$loglik>=current$loglik-1e-6) current <- refined
        if (identical(refined$convergence,0L) || gain<1e-4) break
      }
    }
    output[[key]] <- current
    previous <- current
  }
  output
}

lower_values <- sort(grid[grid<d_hat & abs(grid)>1e-12],decreasing=TRUE)
upper_values <- sort(grid[grid>d_hat])
selected <- c(selected,.run_direction(lower_values,"lower"),
  .run_direction(upper_values,"upper"))
# Duplicate names arise only if an exact anchor was also encountered while
# walking a direction; retain the best likelihood for each fixed difference.
selected <- split(selected,names(selected)) |>
  lapply(function(x) x[[which.max(vapply(x,`[[`,numeric(1),"loglik"))]])

selected_rows <- do.call(rbind,lapply(names(selected),function(key) {
  fit <- selected[[key]]
  d <- if (!is.null(fit$profile_difference)) fit$profile_difference else
    with(fit$params,timegap_reliability_shift-tenure_reliability_shift)
  data.frame(difference=d,loglik=fit$loglik,
    convergence=fit$convergence,iterations=fit$iterations,
    tenure_shift=fit$params$tenure_reliability_shift,
    timegap_shift=fit$params$timegap_reliability_shift)
})) |>
  arrange(difference) |>
  mutate(relative_loglik=loglik-unrestricted$loglik,
    LR=pmax(0,-2*relative_loglik))

start_rows <- if (length(all_starts)) do.call(rbind,
  lapply(names(all_starts),function(label) {
    fit <- all_starts[[label]]
    data.frame(label=label,difference=fit$profile_difference,
      loglik=fit$loglik,convergence=fit$convergence,
      iterations=fit$iterations)
  })) else data.frame()

.profile_crossing <- function(x,y,target=qchisq(.95,1)) {
  if (length(x)<2L) return(NA_real_)
  hit <- which(diff(sign(y-target))!=0)
  if (!length(hit)) return(NA_real_)
  j <- hit[1L]
  x[j]+(target-y[j])*(x[j+1L]-x[j])/(y[j+1L]-y[j])
}
left <- selected_rows[selected_rows$difference<=d_hat,] |>
  arrange(desc(difference))
right <- selected_rows[selected_rows$difference>=d_hat,] |>
  arrange(difference)
lower <- .profile_crossing(left$difference,left$LR)
upper <- .profile_crossing(right$difference,right$LR)
interval <- data.frame(estimate=d_hat,profile_lower=lower,
  profile_upper=upper,threshold=qchisq(.95,1),
  equality_LR=2*(unrestricted$loglik-common$loglik),
  equality_p_value=pchisq(2*(unrestricted$loglik-common$loglik),1,
    lower.tail=FALSE))

write.csv(start_rows,file.path(outdir,"all_starts_latest.csv"),
  row.names=FALSE)
write.csv(selected_rows,file.path(outdir,"profile_likelihood_latest.csv"),
  row.names=FALSE)
write.csv(interval,file.path(outdir,"profile_interval_latest.csv"),
  row.names=FALSE)
saveRDS(list(profile=selected_rows,interval=interval,
  all_starts=start_rows,fits=selected),file.path(outdir,
  "profile_latest.rds"))

cat("\nSelected nuisance-adjusted dispersion-difference profile\n")
print(selected_rows,row.names=FALSE,digits=12)
cat("\nProfile-likelihood interval and equality test\n")
print(interval,row.names=FALSE,digits=12)
if (any(selected_rows$loglik>unrestricted$loglik+.05))
  warning("A fixed-difference point materially exceeds the saved unrestricted fit")
if (any(!is.finite(c(lower,upper))))
  stop("Profile grid did not bracket both 95% interval crossings; extend offsets")
