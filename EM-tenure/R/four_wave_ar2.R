# Two AR(2) return effects; reporting and baseline duration laws are unchanged.
# Source after source_all.R. The validated AR(1) implementation is not edited.
.ar2_effects <- function(p) {
  b <- c(entry = p$ar2_entry %||% 0, exit = p$ar2_exit %||% 0)
  if (length(b) != 2L || any(!is.finite(b))) stop("Invalid AR2 effects")
  b
}
.ar2_shift <- function(risk, effect) {
  if (length(effect) == 1L && effect == 0) return(risk)
  plogis(qlogis(risk) + effect)
}
.pack_four_wave_ar2 <- function(p) c(.pack_four_wave_preferred(p),
  ar2_entry = unname(.ar2_effects(p)[1]), ar2_exit = unname(.ar2_effects(p)[2]))
.unpack_four_wave_ar2 <- function(z) {
  p <- .piecewise_calendar_revision_monthly_unpack(
    z[!names(z) %in% c("ar2_entry", "ar2_exit")], "joint_marginal")
  p$ar2_entry <- unname(z["ar2_entry"])
  p$ar2_exit <- unname(z["ar2_exit"])
  p$timegap_clock_model <- "continuous_joint"
  p$employment_model <- "AR2_return_logit"
  p$initial_pair_model <- "alpha_plus_AR1_first_transition"
  p
}
.four_wave_ar2_bounds <- function(z) {
  old <- .four_wave_parameter_bounds(z[!names(z) %in% c("ar2_entry", "ar2_exit")])
  list(lower=c(old$lower, ar2_entry=-8, ar2_exit=-8),
    upper=c(old$upper, ar2_entry=8, ar2_exit=8))
}

# Original risks condition on reports; missing wrong-state clocks are integrated
# using the existing baseline distribution. Shift the resulting transition ODDS,
# not the hazards defining that duration distribution.
.ar2_baseline_risks <- function(data, p) {
  entry_lut <- .piecewise_category_transition_risk_exact(1:7, p$lambda_d)
  entry_missing <- .piecewise_mean_transition_risk_exact(p$lambda_d)
  exit_missing <- .piecewise_mean_transition_risk_exact(p$lambda_g)
  entry <- exit <- matrix(NA_real_, nrow(data$y), 3L)
  for (t in 1:3) {
    exit[,t] <- .duration_transition_probability(pmax(0,data$month[,t])/12,p$lambda_g)
    exit[data$y[,t]==0,t] <- exit_missing
    entry[,t] <- entry_missing
    observed <- data$y[,t]==0
    entry[observed,t] <- entry_lut[data$category[observed,t]]
  }
  list(entry=entry,exit=exit)
}

.ar2_logratio <- function(data, p, risks = .ar2_baseline_risks(data,p)) {
  h <- latent_histories_eps_4w()
  b <- .ar2_effects(p)
  out <- matrix(0,nrow(data$y),16L)
  if (all(b==0)) return(out)
  for (t in 2:3) for (j in 1:16) {
    if (h[j,t]==h[j,t-1]) next
    origin <- if (h[j,t]) "exit" else "entry"
    r <- risks[[origin]][,t]
    r2 <- .ar2_shift(r,b[origin])
    out[,j] <- out[,j] + if (h[j,t+1]!=h[j,t])
      log(r2)-log(r) else log1p(-r2)-log1p(-r)
  }
  out
}

# Exact identity: L2 = L1 * sum_h posterior1(h) * prior2(h)/prior1(h).
# Nuisance parameters are re-estimated in L1 at EVERY objective evaluation.
# This is not a two-step/fixed-posterior estimator. AR1 class/job aggregates
# cannot be reused after reweighting; obtain those from the R reference if needed.
evaluate_four_wave_ar2 <- function(data, p, posterior=TRUE, threads=8L,
    base_fit=NULL) {
  if (!identical(p$timegap_clock_model %||% "continuous_joint","continuous_joint"))
    stop("AR2 requires the repaired continuous clock")
  b <- .ar2_effects(p)
  if (is.null(base_fit)) base_fit <- evaluate_four_wave_monthly_fast(data,p,
    posterior=TRUE,threads=threads,timegap_clock="continuous_joint")
  if (any(!is.finite(base_fit$row_loglik)) || nrow(base_fit$gamma)!=nrow(data$y))
    stop("Invalid AR1 likelihood supplied to AR2 evaluator")
  ratio <- .ar2_logratio(data,p)
  logterms <- log(base_fit$gamma)+ratio
  adjustment <- .row_logsumexp(logterms)
  rowll <- base_fit$row_loglik+adjustment
  if (any(!is.finite(rowll))) stop("Nonfinite AR2 likelihood")
  list(loglik=sum(data$weight*rowll),row_loglik=rowll,
    gamma=if (posterior) exp(logterms-adjustment) else NULL)
}

# Independent direct-prior R reference. Rebind only the history-prior function
# in a private environment; retain the original reporting and mixture recursion.
e_step_four_wave_ar2_reference <- function(df,params) {
  ref <- new.env(parent=environment(e_step_eps_4w))
  ref$e_step_eps_4w <- e_step_eps_4w
  environment(ref$e_step_eps_4w) <- ref
  ref$.log_duration_history_prior_eps_4w <- function(hmat,alpha,s,g,d,
      lambda_g,beta_g,lambda_d,beta_d,exact_risk=FALSE) {
    stopifnot(exact_risk)
    n <- length(s[[1]])
    out <- matrix(0,n,nrow(hmat))
    exit <- entry <- vector("list",3)
    for (t in 1:3) {
      exit[[t]] <- .duration_transition_probability(g[[t]],lambda_g,beta_g)
      miss <- s[[t]]==0 | !is.finite(exit[[t]])
      exit[[t]][miss] <- .piecewise_mean_transition_risk_exact(lambda_g)
      entry[[t]] <- .piecewise_category_transition_risk_exact(d[[t]],lambda_d)
      miss <- s[[t]]==1 | !is.finite(entry[[t]])
      entry[[t]][miss] <- .piecewise_mean_transition_risk_exact(lambda_d)
    }
    for (j in seq_len(nrow(hmat))) {
      h <- hmat[j,]
      value <- rep(if (h[1]) log(alpha) else log1p(-alpha),n)
      for (t in 1:3) {
        risk <- if (h[t]) exit[[t]] else entry[[t]]
        if (t>=2 && h[t]!=h[t-1]) {
          effect <- if (h[t]) params$ar2_exit %||% 0 else params$ar2_entry %||% 0
          risk <- plogis(qlogis(risk)+effect)
        }
        value <- value + if (h[t+1]!=h[t]) log(risk) else log1p(-risk)
      }
      out[,j] <- value
    }
    out
  }
  ref$e_step_eps_4w(df,params,exact_risk=TRUE,timegap_clock="continuous_joint")
}

duration_weighted_transition_rates_ar2 <- function(df,fit) {
  data <- prepare_four_wave_kernel_data(df)
  risks <- .ar2_baseline_risks(data,fit$params)
  b <- .ar2_effects(fit$params)
  h <- latent_histories_eps_4w()
  rows <- lapply(1:3,function(t) {
    sums <- c(entry_numerator=0,entry_denominator=0,exit_numerator=0,exit_denominator=0)
    for (j in 1:16) {
      origin <- if (h[j,t]) "exit" else "entry"
      risk <- risks[[origin]][,t]
      if (t>=2 && h[j,t]!=h[j,t-1]) risk <- .ar2_shift(risk,b[origin])
      w <- df$weight*fit$gamma[,j]
      sums[paste0(origin,"_numerator")] <- sums[paste0(origin,"_numerator")]+sum(w*risk)
      sums[paste0(origin,"_denominator")] <- sums[paste0(origin,"_denominator")]+sum(w)
    }
    data.frame(wave=t,as.list(sums))
  })
  out <- do.call(rbind,rows)
  out <- rbind(data.frame(wave=0,as.list(colSums(out[,-1]))),out)
  out$entry_rate <- out$entry_numerator/out$entry_denominator
  out$exit_rate <- out$exit_numerator/out$exit_denominator
  out
}

four_wave_ar2_source_fingerprint <- function() c(four_wave_fast_source_fingerprint(),
  tools::md5sum(c("EM-tenure/R/four_wave_ar2.R","EM-tenure/R/simulate_four_wave_ar2.R",
    "EM-tenure/R/fit_four_wave_ar2.R")))
