# Diagnostics and exact-cell acceleration for the epsilon contamination model.

collapse_eps_cells <- function(df, normalize_weights = TRUE,
    allow_zero_tenure = FALSE, extra_cols = character()) {
  validate_df_eps(df, allow_zero_tenure = allow_zero_tenure)
  if (length(setdiff(extra_cols, names(df))))
    stop("Unknown extra collapse columns: ",
      paste(setdiff(extra_cols, names(df)), collapse = ", "))
  cols <- c(paste0("y", 1:3), paste0("tenure", 1:3),
    paste0("timegap_cat", 1:3), extra_cols)
  key <- do.call(paste, c(df[cols], sep = "\r"))
  first <- !duplicated(key)
  group <- match(key, key[first])
  out <- df[first, c(cols, "weight"), drop = FALSE]
  out$weight <- as.vector(rowsum(df$weight, group, reorder = FALSE))
  if (normalize_weights) out$weight <- nrow(df) * out$weight / sum(out$weight)
  attr(out, "n_original") <- nrow(df)
  attr(out, "n_cells") <- nrow(out)
  out
}

.perturb_eps_start <- function(params, linked, stationary, sd, seed) {
  set.seed(seed)
  out <- params
  for (nm in c("alpha", "theta0", "theta1", "pi", "eps"))
    out[[nm]] <- plogis(qlogis(.bound01(out[[nm]])) + rnorm(1, 0, sd))
  if (!linked) {
    out$lambda_g <- out$lambda_g * exp(rnorm(1, 0, sd))
    out$lambda_d <- out$lambda_d * exp(rnorm(1, 0, sd))
  } else {
    out$lambda_g <- ctmc_lambda_from_persistence(out$theta1)
    out$lambda_d <- ctmc_lambda_from_transition(out$theta0)
  }
  if (stationary) out$alpha <- out$theta0 / (out$theta0 + 1 - out$theta1)
  out
}

em_multistart_eps <- function(df, stationary = FALSE, linked = FALSE,
                              K = 3L, perturb_sd = 0.45, seed = 2718L,
                              max_iter = 1000L, tol = 1e-10,
                              param_tol = 1e-7, verbose = 0L) {
  if (K < 1L) stop("K must be positive")
  base <- init_params_eps(df, linked = linked)
  starts <- c(list(base), lapply(seq_len(K - 1L), function(k)
    .perturb_eps_start(base, linked, stationary, perturb_sd, seed + k)))
  fits <- lapply(seq_along(starts), function(k) tryCatch(
    em_fit_tenure_eps(df, params0=starts[[k]], stationary=stationary,
      linked=linked, max_iter=max_iter, tol=tol, param_tol=param_tol,
      verbose=verbose),
    error=function(e) structure(list(error=conditionMessage(e)), class="eps_fit_error")))
  summary <- do.call(rbind, lapply(seq_along(fits), function(k) {
    fit <- fits[[k]]
    if (inherits(fit, "eps_fit_error")) return(data.frame(start=k,status="error",
      converged=FALSE,iterations=NA,loglik=NA,fixedpoint_residual=NA,error=fit$error))
    data.frame(start=k,status=fit$status,converged=fit$converged,
      iterations=fit$iterations,loglik=fit$loglik,
      fixedpoint_residual=fit$diagnostics$fixedpoint_residual,error=NA_character_)
  }))
  eligible <- which(summary$converged & is.finite(summary$loglik))
  if (!length(eligible)) stop("No epsilon-model start genuinely converged: ",
                              paste(summary$status, collapse=", "))
  best_index <- eligible[which.max(summary$loglik[eligible])]
  list(best=fits[[best_index]], summary=summary, best_start=best_index)
}
# Weighted duration-reporting patterns used to assess serial contamination
# models. "Irregularity" means a tenure report that does not advance by three
# months within employment, or a timegap transition with zero clean-process
# probability within nonemployment.
.duration_reporting_patterns <- function(df,params) {
  n <- nrow(df)
  tenure_deviation <- timegap_inconsistent <- opportunities <- integer(n)
  exact <- local <- yearly <- other <- logical()
  calendar_weight <- numeric()
  timegap_bad <- logical(); timegap_weight <- numeric()
  weight <- if ("weight" %in% names(df)) df$weight else rep(1,n)
  weighted_share <- function(value,w) {
    keep <- is.finite(value) & is.finite(w) & w>0
    if (!any(keep)) return(NA_real_)
    weighted.mean(as.numeric(value[keep]),w[keep])
  }
  for (j in 1:2) {
    employed <- df[[paste0("y",j)]]==1L &
      df[[paste0("y",j+1L)]]==1L
    delta <- round(12*df[[paste0("tenure",j+1L)]])-
      round(12*df[[paste0("tenure",j)]])-3L
    tenure_deviation <- tenure_deviation+as.integer(employed & delta!=0L)
    opportunities <- opportunities+as.integer(employed)
    exact <- c(exact,delta[employed]==0L)
    local <- c(local,delta[employed]!=0L & abs(delta[employed])<=6L)
    yearly <- c(yearly,delta[employed]!=0L & delta[employed]%%12L==0L)
    other <- c(other,delta[employed]!=0L & abs(delta[employed])>6L &
      delta[employed]%%12L!=0L)
    calendar_weight <- c(calendar_weight,weight[employed])

    nonemployed <- df[[paste0("y",j)]]==0L &
      df[[paste0("y",j+1L)]]==0L
    clean_lp <- log_emission_transition_d(
      df[[paste0("timegap_cat",j+1L)]][nonemployed],
      df[[paste0("timegap_cat",j)]][nonemployed],params$lambda_d,
      if (is.null(params$beta_d)) 0 else params$beta_d)
    bad <- !is.finite(clean_lp)
    timegap_inconsistent[nonemployed] <-
      timegap_inconsistent[nonemployed]+as.integer(bad)
    opportunities <- opportunities+as.integer(nonemployed)
    timegap_bad <- c(timegap_bad,bad)
    timegap_weight <- c(timegap_weight,weight[nonemployed])
  }
  irregular <- tenure_deviation+timegap_inconsistent
  two <- opportunities==2L
  metrics <- c(exact_continuation=weighted_share(exact,calendar_weight),
    local_revision_1_to_6_months=weighted_share(local,calendar_weight),
    whole_year_revision=weighted_share(yearly,calendar_weight),
    other_revision=weighted_share(other,calendar_weight),
    timegap_inconsistent=weighted_share(timegap_bad,timegap_weight),
    two_opportunities_no_irregularity=weighted_share(two & irregular==0L,
      weight)/weighted_share(two,weight),
    two_opportunities_one_irregularity=weighted_share(two & irregular==1L,
      weight)/weighted_share(two,weight),
    two_opportunities_two_irregularities=weighted_share(two & irregular==2L,
      weight)/weighted_share(two,weight))
  pattern <- ifelse(opportunities==0L,"No repeated-state duration",
    ifelse(tenure_deviation>0L,"Tenure deviation",
      ifelse(timegap_inconsistent>0L,"Timegap inconsistency",
        "Repeated state; no irregularity")))
  list(metrics=metrics,pattern=pattern,opportunities=opportunities,
    irregular=irregular,tenure_deviation=tenure_deviation,
    timegap_inconsistent=timegap_inconsistent)
}
