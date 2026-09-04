# Exhaustive finite-support checks in a nested, rapidly decaying tenure case.
# No empirical model is changed or re-estimated.
source("EM-tenure/R/source_all.R")
base <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
clock_model <- match.arg(Sys.getenv("FOUR_WAVE_TIMEGAP_CLOCK","continuous_joint"),
  c("continuous_joint","legacy_pairwise"))
outdir <- file.path(base,if (clock_model=="continuous_joint")
  "recovery_continuous_clock" else "recovery")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
p <- readRDS(file.path(base,"converged_comparison_latest.rds"))$best$params
p$alpha <- .5
p$lambda_g <- rep(20,5)
p$eps <- .2
p$job_change_prob <- 0
p$tenure_start_month_probs <- rep(1/12,12)
for (name in c("tenure_reliability_shift","timegap_reliability_shift",
    "tenure_heaping_prob","tenure_year_revision_prob",
    "tenure_clean_anchor_revision_prob","tenure_exact_anchor_retention_prob",
    "tenure_local_revision_prob")) p[[name]] <- 0
max_month <- as.integer(Sys.getenv("FOUR_WAVE_MASS_MAX_MONTH","12"))
# Categories 1:7 denote reported nonemployment; the remaining codes encode
# reported employment with tenure 0:max_month months.
support <- seq_len(7L+max_month+1L)
grid <- as.matrix(expand.grid(rep(list(support),4)))
d <- data.frame(weight=rep(1,nrow(grid)))
for (j in 1:4) {
  d[[paste0("y",j)]] <- as.integer(grid[,j]>7L)
  d[[paste0("timegap_cat",j)]] <- ifelse(grid[,j]<=7L,grid[,j],NA_integer_)
  d[[paste0("tenure",j)]] <- ifelse(grid[,j]>7L,(grid[,j]-8L)/12,NA_real_)
  d[[paste0("interview_month",j)]] <- 3L*j
}
data <- prepare_four_wave_kernel_data(d)
cases <- expand.grid(varying_hazard=c(FALSE,TRUE),status_error=c(0,.15),
  timegap_error=c(1e-8,.3))
cases$probability_mass <- NA_real_
for (i in seq_len(nrow(cases))) {
  p$lambda_d <- if (cases$varying_hazard[i]) c(4,1,.5,.3,.2) else rep(.4,5)
  p$pi <- cases$status_error[i]
  p$eps_d <- cases$timegap_error[i]
  cases$probability_mass[i] <- sum(exp(evaluate_four_wave_monthly_fast(data,p,
    posterior=FALSE,threads=8L,timegap_clock=clock_model)$row_loglik))
  print(cases[i,],row.names=FALSE,digits=12)
}
cases$mass_minus_one <- cases$probability_mass-1
reference_mass <- NA_real_
if (identical(Sys.getenv("FOUR_WAVE_MASS_REFERENCE","false"),"true")) {
  reference <- e_step_eps_4w(d,p,exact_risk=TRUE,suff_stats=FALSE,timegap_clock=clock_model)
  reference_mass <- sum(exp(reference$row_loglik))
  print(c(reference_mass=reference_mass,
    compiled_mass=tail(cases$probability_mass,1)),digits=14)
}
write.csv(cases,file.path(outdir,paste0("probability_mass_m",max_month,".csv")),
  row.names=FALSE)
# Separate deterministic illustration of the future-report-dependent timegap
# rule: a fourth false employment report switches a 0000 spell from the joint
# kernel to pairwise fallback even for its first three observed categories.
cats <- as.matrix(expand.grid(rep(list(1:7),3)))
hazards <- c(4,1,.5,.3,.2)
eps <- .3
joint <- exp(log_emission_timegap_spell_joint(cats,hazards,eps_d=eps))
pairwise <- exp(log_emission_interval_d(cats[,1],hazards)+
  log_emission_transition_d_contaminated(cats[,2],cats[,1],hazards,
    eps_d=1-(1-eps)^2)+
  log_emission_transition_d_contaminated(cats[,3],cats[,2],hazards,
    eps_d=1-(1-eps)^2))
prefix <- data.frame(joint_mass=sum(joint),pairwise_mass=sum(pairwise),
  total_variation=sum(abs(joint-pairwise))/2,
  max_cell_difference=max(abs(joint-pairwise)))
print(prefix,digits=12)
# A clean continuous clock retained across all waves is not a first-order
# Markov process after categorization. Compare its exact intersection masses
# with the likelihood's product of pairwise category transitions.
cats4 <- as.matrix(expand.grid(rep(list(1:7),4)))
exact_clock <- apply(cats4,1,function(path) {
  intervals <- vapply(path,.timegap_interval,numeric(2))
  lower <- max(intervals[1,]-(0:3)*.QUARTER_YEARS)
  upper <- min(intervals[2,]-(0:3)*.QUARTER_YEARS)
  if (lower>=upper) 0 else exp(.log_duration_interval_prob(lower,upper,hazards))
})
markov_clock <- exp(log_emission_timegap_spell_joint(cats4,hazards,eps_d=1e-12))
clock <- data.frame(exact_clock_mass=sum(exact_clock),markov_mass=sum(markov_clock),
  total_variation=sum(abs(exact_clock-markov_clock))/2)
print(clock,digits=12)
# Every tenure component in this submodel is either a marginal exponential
# draw or such a draw shifted by at most nine months (resets are shorter).
# Summing a union bound over four reports and sixteen histories is conservative.
tail_bound <- 64*exp(-20*(max_month+1-9)/12)
saveRDS(list(cases=cases,max_month=max_month,n_cells=nrow(grid),prefix=prefix,
  timegap_clock=clock_model,
  clock=clock,reference_mass=reference_mass,conservative_tail_bound=tail_bound,
  note="Use increasing tenure cutoffs and the tail bound to distinguish missing tail mass",
  source_md5=four_wave_fast_source_fingerprint()),
  file.path(outdir,paste0("probability_mass_m",max_month,".rds")))
