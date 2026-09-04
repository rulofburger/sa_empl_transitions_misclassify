# Nuisance-adjusted multistart profile likelihood for the whole-year revision
# probability. The grid is centered and scaled using analytical inference;
# each noncentral point is optimized from two distinct nuisance starts.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
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

root <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_revision")
fit_file <- file.path(root,"fits_latest.rds")
inference_file <- file.path(root,"inference","analytical_inference_latest.rds")
if (!file.exists(fit_file) || !file.exists(inference_file))
  stop("Estimate and run analytical inference for calendar revisions first")
fits <- readRDS(fit_file)
preferred <- fits$revision
alternative <- fits$fits$main
inference <- readRDS(inference_file)
row <- inference$summary[inference$summary$quantity==
  "whole_year_revision_prob",]
if (nrow(row)!=1L || !is.finite(row$se) || row$se<=0)
  stop("Valid analytical revision-probability standard error not found")

outdir <- file.path(root,"inference","profile")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
workers <- as.integer(Sys.getenv("REVISION_PROFILE_WORKERS","1"))
maxit <- as.integer(Sys.getenv("REVISION_PROFILE_MAXIT","40"))
resume <- identical(tolower(Sys.getenv("REVISION_PROFILE_RESUME","true")),
  "true")
multipliers <- c(-2,-1,0,1,2)
grid <- pmin(.95,pmax(1e-5,row$estimate+multipliers*row$se))
grid <- sort(unique(grid))
grid_definition <- data.frame(revision_prob=grid,
  analytical_se_units=(grid-row$estimate)/row$se)
write.csv(grid_definition,file.path(outdir,"profile_grid_latest.csv"),
  row.names=FALSE)

.safe_name <- function(x) gsub("\\.","_",formatC(x,digits=8,format="f"))
.fit_point <- function(omega,start,label) {
  path <- file.path(outdir,paste0("profile_",.safe_name(omega),"_",label,
    ".rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message(sprintf("Profiling omega=%.8f from %s",omega,label))
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start$params,
    heaping_start=start$params$tenure_heaping_prob,
    revision_start=omega,q_start=start$params$job_change_prob,
    revision_fixed=omega,maxit=maxit,reltol=1e-10,pgtol=1e-7,
    workers=workers,verbose=1L)
  fit$objective_function <- NULL
  fit$profile_revision_prob <- omega
  fit$profile_start <- label
  saveRDS(fit,path)
  fit
}

centre <- which.min(abs(grid-row$estimate))
results <- list()
# The exact unrestricted solution supplies the central point.
results[[paste0(centre,"_preferred")]] <- preferred
results[[paste0(centre,"_alternative")]] <- preferred

# Move outward from the MLE so the first start at every point is the nearest
# already optimized profile fit. The second always returns to the distinct
# central local-mode nuisance estimates as a safeguard.
for (direction in c(-1L,1L)) {
  previous <- preferred
  index <- centre+direction
  while (index>=1L && index<=length(grid)) {
    omega <- grid[index]
    warm <- .fit_point(omega,previous,paste0("warm_",ifelse(direction<0,
      "lower","upper")))
    alt <- .fit_point(omega,alternative,"alternative")
    results[[paste0(index,"_warm")]] <- warm
    results[[paste0(index,"_alternative")]] <- alt
    candidates <- list(warm=warm,alternative=alt)
    previous <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),
      "loglik"))]]
    index <- index+direction
  }
}

rows <- do.call(rbind,lapply(names(results),function(key) {
  split <- strsplit(key,"_",fixed=TRUE)[[1L]]
  index <- as.integer(split[1L]); fit <- results[[key]]
  data.frame(revision_prob=grid[index],start=paste(split[-1L],collapse="_"),
    loglik=fit$loglik,convergence=fit$convergence,
    iterations=fit$iterations)
}))
best <- rows |>
  group_by(revision_prob) |>
  slice_max(loglik,n=1,with_ties=FALSE) |>
  ungroup() |>
  arrange(revision_prob) |>
  mutate(relative_loglik=loglik-max(loglik),LR=-2*relative_loglik)

.profile_crossing <- function(x,y,target=3.841459) {
  hit <- which(diff(sign(y-target))!=0)
  if (!length(hit)) return(NA_real_)
  j <- hit[1L]
  x[j]+(target-y[j])*(x[j+1L]-x[j])/(y[j+1L]-y[j])
}
left <- best[best$revision_prob<=preferred$params$tenure_year_revision_prob,]
right <- best[best$revision_prob>=preferred$params$tenure_year_revision_prob,]
lower <- .profile_crossing(rev(left$revision_prob),rev(left$LR))
upper <- .profile_crossing(right$revision_prob,right$LR)
interval <- data.frame(estimate=preferred$params$tenure_year_revision_prob,
  analytical_lower=row$lower,analytical_upper=row$upper,
  profile_lower=lower,profile_upper=upper,threshold=3.841459)
write.csv(rows,file.path(outdir,"all_starts_latest.csv"),row.names=FALSE)
write.csv(best,file.path(outdir,"profile_likelihood_latest.csv"),
  row.names=FALSE)
write.csv(interval,file.path(outdir,"profile_interval_latest.csv"),
  row.names=FALSE)
saveRDS(list(grid=grid_definition,all_starts=rows,profile=best,
  interval=interval,results=results),file.path(outdir,"profile_latest.rds"))

cat("\nAll calendar-revision profile starts\n")
print(rows,row.names=FALSE,digits=9)
cat("\nSelected nuisance-adjusted profile\n")
print(best,row.names=FALSE,digits=9)
cat("\nAnalytical and profile intervals\n")
print(interval,row.names=FALSE,digits=9)
if (any(best$convergence!=0L)) stop("A selected profile point did not converge")
if (any(!is.finite(c(lower,upper))))
  stop("Profile grid did not bracket both 95% interval crossings")
