# Diagnostics and exact-cell acceleration for the epsilon contamination model.

collapse_eps_cells <- function(df, normalize_weights = TRUE) {
  validate_df_eps(df)
  cols <- c(paste0("y", 1:3), paste0("tenure", 1:3), paste0("timegap_cat", 1:3))
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
