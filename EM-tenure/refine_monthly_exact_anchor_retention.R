# Continue the full-sample exact-anchor model from its current polished fit,
# and screen dispersed retention starts without overwriting a superior saved
# solution. Every completed stage is cached and resumable.

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

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "calendar_exact_anchor_retention")
fit_file <- file.path(outdir,"fits_latest.rds")
if (!file.exists(fit_file)) stop("Estimate exact-anchor retention first")
saved <- readRDS(fit_file)
incumbent <- saved$extension
workers <- as.integer(Sys.getenv("EXACT_ANCHOR_REFINE_WORKERS","8"))
main_maxit <- as.integer(Sys.getenv("EXACT_ANCHOR_REFINE_MAXIT","8"))
screen_maxit <- as.integer(Sys.getenv("EXACT_ANCHOR_SCREEN_MAXIT","3"))
resume <- identical(tolower(Sys.getenv("EXACT_ANCHOR_REFINE_RESUME","true")),
  "true")

fit_stage <- function(label,start,scheme,iterations) {
  path <- file.path(outdir,paste0("refine_",label,"_latest.rds"))
  if (resume && file.exists(path)) return(readRDS(path))
  message("Exact-anchor refinement: ",label)
  fit <- fit_eps_piecewise_calendar_revision_monthly(df_fit,start$params,
    heaping_start=max(start$params$tenure_heaping_prob,1e-8),
    revision_start=start$params$tenure_year_revision_prob,
    clean_anchor_revision_start=
      start$params$tenure_clean_anchor_revision_prob,
    start_month_probs_start=start$params$tenure_start_month_probs,
    exact_anchor_retention_start=
      start$params$tenure_exact_anchor_retention_prob,
    q_start=start$params$job_change_prob,maxit=iterations,reltol=1e-10,
    pgtol=1e-8,workers=max(1L,min(workers,8L)),verbose=1L,
    gradient_step=5e-4,gradient_scheme=scheme)
  fit$objective_function <- NULL
  fit$stage <- label
  saveRDS(fit,path)
  fit
}

# Continue the incumbent with the more accurate central gradient.
continued <- fit_stage("incumbent_central",incumbent,"central",main_maxit)

# Dispersed rho starts retain the incumbent nuisance values. Short forward
# screens are enough to reveal a competing basin; only a superior basin is
# promoted to the longer central refinement.
make_start <- function(rho,label) {
  value <- incumbent
  value$params$tenure_exact_anchor_retention_prob <- rho
  eval <- e_step_eps(df_fit,value$params,check_df=FALSE,suff_stats=FALSE)
  value$loglik <- eval$loglik
  value$gamma <- eval$gamma
  value$job_change_posterior <- eval$job_change_posterior
  value$stage <- label
  value
}
low <- fit_stage("low_rho_forward",make_start(.03,"low_rho"),"forward",
  screen_maxit)
high <- fit_stage("high_rho_forward",make_start(.40,"high_rho"),"forward",
  screen_maxit)
screened <- list(incumbent=incumbent,continued=continued,low=low,high=high)
best_screen_name <- names(screened)[which.max(vapply(screened,`[[`,numeric(1),
  "loglik"))]
best_screen <- screened[[best_screen_name]]
promoted_path <- file.path(outdir,"refine_promoted_central_latest.rds")
if (resume && file.exists(promoted_path)) {
  promoted <- readRDS(promoted_path)
} else if (best_screen_name %in% c("low","high") &&
    best_screen$loglik>max(incumbent$loglik,continued$loglik)+1e-6) {
  promoted <- fit_stage("promoted_central",best_screen,"central",main_maxit)
} else promoted <- continued
best <- list(incumbent=incumbent,continued=continued,promoted=promoted,
  low=low,high=high)[[which.max(vapply(list(incumbent=incumbent,
    continued=continued,promoted=promoted,low=low,high=high),`[[`,
    numeric(1),"loglik"))]]

# Finish with an exact conditional maximization of rho. A synchronized cached
# polish is reusable because all nuisance parameters and the empirical data are
# unchanged on a resumed run.
polished_path <- file.path(outdir,"refine_rho_polished_latest.rds")
polished <- if (resume && file.exists(polished_path))
  readRDS(polished_path) else NULL
polished_valid <- !is.null(polished) &&
  "tenure_exact_anchor_retention" %in% names(polished$par_unconstrained) &&
  isTRUE(all.equal(plogis(unname(polished$par_unconstrained[
    "tenure_exact_anchor_retention"])),
    polished$params$tenure_exact_anchor_retention_prob,tolerance=1e-10)) &&
  isTRUE(all.equal(polished$params[names(best$params)],best$params,
    tolerance=2e-4))
if (!polished_valid) {
  rho_objective <- function(rho) {
    p <- best$params; p$tenure_exact_anchor_retention_prob <- rho
    -e_step_eps(df_fit,p,check_df=FALSE,suff_stats=FALSE)$loglik
  }
  rho_opt <- optimize(rho_objective,c(1e-7,.95),tol=1e-7)
  polished_params <- best$params
  polished_params$tenure_exact_anchor_retention_prob <- rho_opt$minimum
  polished_eval <- e_step_eps(df_fit,polished_params,check_df=FALSE,
    suff_stats=FALSE)
  polished <- best
  polished$params <- polished_params
  polished$par_unconstrained["tenure_exact_anchor_retention"] <-
    qlogis(rho_opt$minimum)
  polished$loglik <- polished_eval$loglik
  polished$gamma <- polished_eval$gamma
  polished$job_change_posterior <- polished_eval$job_change_posterior
  polished$rho_polished <- TRUE
  polished$stage <- "refined_rho_polished"
  saveRDS(polished,polished_path)
}
if (polished$loglik>=best$loglik) best <- polished

all_fits <- list(incumbent=incumbent,continued=continued,low=low,high=high,
  promoted=promoted,polished=polished)
refinement <- do.call(rbind,lapply(names(all_fits),function(label) {
  f <- all_fits[[label]]; p <- f$params
  data.frame(stage=label,loglik=f$loglik,
    gain_over_incumbent=f$loglik-incumbent$loglik,
    convergence=f$convergence,iterations=f$iterations,
    exact_anchor_retention=p$tenure_exact_anchor_retention_prob,
    tenure_contamination=p$eps,status_misclassification=p$pi,
    gross_anchor_revision=p$tenure_year_revision_prob,
    clean_anchor_revision=p$tenure_clean_anchor_revision_prob,
    january_heaping=p$tenure_heaping_prob)
}))

if (best$loglik>=saved$extension$loglik-1e-7) {
  saved$fits$long_refinement <- best
  saved$extension <- best
}
saved$refinement <- refinement
final <- saved$extension
summarise_fit <- function(model,fit,k) {
  rates <- duration_weighted_transition_rates(df_fit,fit)[1,]
  jp <- fit$job_change_posterior; p <- fit$params
  data.frame(model=model,loglik=fit$loglik,parameters=k,
    AIC=-2*fit$loglik+2*k,entry_rate=rates$entry_rate,
    exit_rate=rates$exit_rate,initial_employment=p$alpha,
    status_misclassification=p$pi,tenure_contamination=p$eps,
    timegap_contamination=p$eps_d,job_change_prob=p$job_change_prob,
    posterior_job_change_rate=sum(df_fit$weight*jp$expected_changes)/
      sum(df_fit$weight*jp$opportunities),
    gross_january_heaping_prob=p$tenure_heaping_prob,
    gross_anchor_revision_prob=p$tenure_year_revision_prob,
    clean_anchor_revision_prob=p$tenure_clean_anchor_revision_prob,
    exact_anchor_retention_prob=
      if (is.null(p$tenure_exact_anchor_retention_prob)) 0 else
        p$tenure_exact_anchor_retention_prob,
    convergence=fit$convergence,iterations=fit$iterations)
}
saved$comparison <- rbind(
  summarise_fit("Flexible baseline; no exact retention",saved$base,29L),
  summarise_fit("Conditional exact retention only",saved$conditional,30L),
  summarise_fit("Joint exact-anchor retention",final,30L))
LR <- max(0,2*(final$loglik-saved$base$loglik))
saved$lr <- data.frame(comparison="No exact-anchor retention",LR=LR,df=1L,
  p_value_chisq1=pchisq(LR,1,lower.tail=FALSE),
  p_value_boundary_mixture=.5*pchisq(LR,1,lower.tail=FALSE))
saveRDS(saved,fit_file)
write.csv(saved$comparison,file.path(outdir,"model_comparison_latest.csv"),
  row.names=FALSE)
write.csv(saved$lr,file.path(outdir,"lr_latest.csv"),row.names=FALSE)
write.csv(refinement,file.path(outdir,"refinement_summary_latest.csv"),
  row.names=FALSE)
cat("\nExact-anchor longer refinement\n")
print(refinement,row.names=FALSE,digits=10)
cat("\nRetained best stage\n")
print(refinement[which.max(refinement$loglik),],row.names=FALSE,digits=10)
