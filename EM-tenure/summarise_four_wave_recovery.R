# Consolidate failed or successful validation gates, without running a fit.
source("EM-tenure/R/source_all.R")
clock_model <- match.arg(Sys.getenv("FOUR_WAVE_TIMEGAP_CLOCK","continuous_joint"),
  c("continuous_joint","legacy_pairwise"))
base <- file.path("EM-tenure/output/results/job_change_monthly/four_wave_ar1",
  if (clock_model=="continuous_joint") "recovery_continuous_clock" else "recovery")
score <- readRDS(file.path(base,"sequential_score_gate_latest.rds"))
mass <- readRDS(file.path(base,"probability_mass_m24.rds"))
reference <- readRDS(file.path(base,"probability_mass_m12.rds"))
adequate_tail <- mass$conservative_tail_bound < 1e-8
normalization_passed <- adequate_tail && all(abs(mass$cases$mass_minus_one)<1e-8)
reference_agrees <- is.finite(reference$reference_mass) &&
  abs(reference$reference_mass-tail(reference$cases$probability_mass,1))<1e-10
fit_path <- file.path(base,"recovery_fits_latest.rds")
fits <- if (file.exists(fit_path)) readRDS(fit_path) else NULL
if (!is.null(fits) && !identical(fits$source_md5,four_wave_fast_source_fingerprint()))
  fits <- NULL
status <- list(timegap_clock=clock_model,
  score_passed=score$passed,normalization_passed=normalization_passed,
  reference_agrees=reference_agrees,ready_for_recovery_fits=
    isTRUE(score$passed)&&normalization_passed&&reference_agrees,
  recovery_fits_run=if (is.null(fits)) 0L else nrow(fits$summary),
  recovery_fits_all_converged=if (is.null(fits)) NA else fits$all_converged,
  empirical_fit_changed=FALSE,
  n_simulated=sum(unique(score$draws[c("replication","n")])$n),
  score=score$summary,mass=mass$cases,tail_bound=mass$conservative_tail_bound,
  prefix=mass$prefix,clock=mass$clock,
  source_md5=four_wave_fast_source_fingerprint())
saveRDS(status,file.path(base,"recovery_status_latest.rds"))
print(status[c("score_passed","normalization_passed","reference_agrees",
  "ready_for_recovery_fits","recovery_fits_run","empirical_fit_changed")])
print(mass$cases,row.names=FALSE,digits=12)
message("Conservative omitted-tail bound: ",format(mass$conservative_tail_bound))
