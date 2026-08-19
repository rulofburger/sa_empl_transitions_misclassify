# Duration-dependent linked epsilon model.
#
# The continuous-time hazards are
#   h_g(x) = lambda_g * (1 + x)^beta_g
#   h_d(x) = lambda_d * (1 + x)^beta_d.
# beta = 0 nests the constant-hazard CTMC-linked model. The same hazard enters
# the spell-duration emissions and the history transition probabilities.

.duration_eps_fixed_beta <- function(beta_fixed) {
  if (is.null(beta_fixed)) return(c(beta_g = NA_real_, beta_d = NA_real_))
  if (length(beta_fixed) == 1L) return(c(beta_g = beta_fixed, beta_d = beta_fixed))
  if (is.null(names(beta_fixed))) names(beta_fixed) <- c("beta_g", "beta_d")
  beta_fixed[c("beta_g", "beta_d")]
}

.duration_eps_unpack <- function(z, beta_fixed = NULL,
                                 pi_cap = 0.49, eps_floor = 1e-4,
                                 eps_cap = 0.95) {
  fixed <- .duration_eps_fixed_beta(beta_fixed)
  beta_g <- if (is.na(fixed["beta_g"])) unname(z["beta_g"]) else fixed["beta_g"]
  beta_d <- if (is.na(fixed["beta_d"])) unname(z["beta_d"]) else fixed["beta_d"]
  lambda_g <- exp(unname(z["log_lambda_g"]))
  lambda_d <- exp(unname(z["log_lambda_d"]))
  p_exit0 <- .duration_transition_probability(0, lambda_g, beta_g)
  p_entry0 <- .duration_transition_probability(0, lambda_d, beta_d)
  list(
    alpha = plogis(unname(z["alpha"])),
    theta0 = p_entry0,
    theta1 = 1 - p_exit0,
    pi = pi_cap * plogis(unname(z["pi"])),
    eps = eps_floor + (eps_cap - eps_floor) * plogis(unname(z["eps"])),
    lambda_g = lambda_g,
    beta_g = beta_g,
    lambda_d = lambda_d,
    beta_d = beta_d
  )
}

.duration_eps_pack <- function(params, beta_fixed = NULL,
                               pi_cap = 0.49, eps_floor = 1e-4,
                               eps_cap = 0.95) {
  z <- c(
    alpha = qlogis(.bound01(params$alpha)),
    pi = qlogis(.bound01(params$pi / pi_cap)),
    eps = qlogis(.bound01((params$eps - eps_floor) / (eps_cap - eps_floor))),
    log_lambda_g = log(params$lambda_g),
    log_lambda_d = log(params$lambda_d)
  )
  fixed <- .duration_eps_fixed_beta(beta_fixed)
  if (is.na(fixed["beta_g"]))
    z <- c(z, beta_g = if (is.null(params$beta_g)) 0 else params$beta_g)
  if (is.na(fixed["beta_d"]))
    z <- c(z, beta_d = if (is.null(params$beta_d)) 0 else params$beta_d)
  z
}

fit_eps_duration_hazard <- function(df, start, beta_fixed = NULL,
                                    maxit = 400L, reltol = 1e-10,
                                    verbose = 0L) {
  validate_df_eps(df)
  z0 <- .duration_eps_pack(start, beta_fixed)
  lower <- rep(-Inf, length(z0)); upper <- rep(Inf, length(z0))
  names(lower) <- names(upper) <- names(z0)
  lower[c("log_lambda_g", "log_lambda_d")] <- log(1e-5)
  upper[c("log_lambda_g", "log_lambda_d")] <- log(50)
  beta_names <- intersect(c("beta_g", "beta_d"), names(z0))
  lower[beta_names] <- -0.95
  upper[beta_names] <- 3
  total_weight <- sum(df$weight)
  objective <- function(z) {
    params <- .duration_eps_unpack(z, beta_fixed)
    value <- tryCatch(e_step_eps(df, params, check_df = FALSE,
                                 suff_stats = FALSE)$loglik,
                      error = function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value / total_weight
  }
  opt <- optim(z0, objective, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = maxit, factr = reltol / .Machine$double.eps,
                              trace = as.integer(verbose)))
  params <- .duration_eps_unpack(opt$par, beta_fixed)
  estep <- e_step_eps(df, params, check_df = FALSE, suff_stats = FALSE)
  list(params = params, loglik = estep$loglik, gamma = estep$gamma,
       convergence = opt$convergence, message = opt$message,
       iterations = unname(opt$counts["function"]), objective = opt$value,
       gradient_max = NA_real_, par_unconstrained = opt$par,
       objective_function = objective)
}

.duration_eps_numeric_gradient <- function(fn, x, step = 1e-5) {
  vapply(seq_along(x), function(j) {
    xp <- xm <- x; xp[j] <- xp[j] + step; xm[j] <- xm[j] - step
    (fn(xp) - fn(xm)) / (2 * step)
  }, numeric(1))
}

fit_eps_duration_multistart <- function(df, start, beta_starts = list(
                                          c(-0.4, -0.8), c(0, 0),
                                          c(0.25, -0.50)), ...) {
  fits <- lapply(beta_starts, function(b) {
    s <- start; s$beta_g <- b[1]; s$beta_d <- b[2]
    tryCatch(fit_eps_duration_hazard(df, s, ...),
             error = function(e) structure(list(error = conditionMessage(e)),
                                          class = "duration_eps_error"))
  })
  tab <- do.call(rbind, lapply(seq_along(fits), function(k) {
    x <- fits[[k]]
    if (inherits(x, "duration_eps_error")) return(data.frame(
      start = k, convergence = NA_integer_, loglik = NA_real_,
      gradient_max = NA_real_, beta_g = NA_real_, beta_d = NA_real_,
      error = x$error))
    data.frame(start = k, convergence = x$convergence, loglik = x$loglik,
      gradient_max = x$gradient_max, beta_g = x$params$beta_g,
      beta_d = x$params$beta_d, error = NA_character_)
  }))
  eligible <- which(tab$convergence == 0 & is.finite(tab$loglik))
  if (!length(eligible)) stop("No duration-dependent fit converged")
  best <- eligible[which.max(tab$loglik[eligible])]
  grad <- .duration_eps_numeric_gradient(fits[[best]]$objective_function,
                                         fits[[best]]$par_unconstrained)
  fits[[best]]$gradient_max <- max(abs(grad))
  tab$gradient_max[best] <- fits[[best]]$gradient_max
  list(best = fits[[best]], best_start = best, fits = fits, summary = tab)
}

duration_transition_profile <- function(params,
                                        durations = c(0, .25, .5, 1, 2, 5, 10)) {
  data.frame(
    duration_years = durations,
    exit_probability = .duration_transition_probability(
      durations, params$lambda_g, params$beta_g),
    entry_probability = .duration_transition_probability(
      durations, params$lambda_d, params$beta_d)
  )
}

duration_mean_years <- function(lambda, beta) {
  integrate(function(x) exp(-.duration_cumhaz(x, lambda, beta)),
            lower = 0, upper = Inf, rel.tol = 1e-8,
            subdivisions = 1000L)$value
}

duration_weighted_transition_rates <- function(df, fit) {
  h <- latent_histories()
  p <- fit$params
  gamma <- fit$gamma
  rows <- lapply(1:2, function(t) {
    pr_emp <- as.vector(gamma %*% h[, t])
    pr_nonemp <- 1 - pr_emp
    g <- df[[paste0("tenure", t)]]
    cat_d <- df[[paste0("timegap_cat", t)]]
    p_exit <- .duration_transition_probability(g, p$lambda_g, p$beta_g)
    p_entry <- .duration_category_transition_probability(
      cat_d, p$lambda_d, p$beta_d)
    data.frame(wave = t,
      exit_rate = weighted.mean(p_exit, df$weight * pr_emp),
      entry_rate = weighted.mean(p_entry, df$weight * pr_nonemp),
      exit_numerator = sum(df$weight * pr_emp * p_exit),
      exit_denominator = sum(df$weight * pr_emp),
      entry_numerator = sum(df$weight * pr_nonemp * p_entry),
      entry_denominator = sum(df$weight * pr_nonemp))
  })
  out <- do.call(rbind, rows)
  pooled <- data.frame(wave = 0,
    exit_rate = sum(out$exit_numerator) / sum(out$exit_denominator),
    entry_rate = sum(out$entry_numerator) / sum(out$entry_denominator),
    exit_numerator = sum(out$exit_numerator),
    exit_denominator = sum(out$exit_denominator),
    entry_numerator = sum(out$entry_numerator),
    entry_denominator = sum(out$entry_denominator))
  rbind(pooled, out)
}
