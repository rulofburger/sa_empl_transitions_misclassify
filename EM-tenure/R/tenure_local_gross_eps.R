# Local-plus-gross tenure contamination for the piecewise duration model.
#
# Each observed tenure report belongs to one of three branches:
#   C: exact clock consistency;
#   L: a bounded discrete month error around the latent clock; or
#   G: an independent draw from the model-implied marginal tenure density.
#
# The local support deliberately excludes zero, so the exact branch remains a
# separate probability mass.  A fixed kernel avoids the variance-collapse
# pathology of an unrestricted continuous measurement-error density.

.TENURE_LOCAL_MONTHS <- c(-6:-1, 1:6)
.TENURE_LOCAL_DECAY_MONTHS <- 3

.tenure_local_logweights <- function(
    support_months = .TENURE_LOCAL_MONTHS,
    decay_months = .TENURE_LOCAL_DECAY_MONTHS) {
  if (!length(support_months) || any(support_months == 0L) ||
      any(!is.finite(support_months)))
    stop("Local tenure-error support must be finite and exclude zero")
  if (!is.finite(decay_months) || decay_months <= 0)
    stop("Local tenure-error decay must be positive")
  z <- -abs(support_months) / decay_months
  z - (max(z) + log(sum(exp(z - max(z)))))
}

.tenure_local_match_logprob <- function(
    error_years, support_months = .TENURE_LOCAL_MONTHS,
    logweights = .tenure_local_logweights(support_months), tol = 1e-7) {
  error_months <- 12 * error_years
  rounded <- round(error_months)
  valid <- is.finite(error_months) &
    abs(error_months - rounded) < tol & rounded %in% support_months
  out <- rep(-Inf, length(error_years))
  if (any(valid))
    out[valid] <- logweights[match(rounded[valid], support_months)]
  out
}

.row_logsumexp <- function(x) {
  if (!nrow(x)) return(numeric(0))
  mx <- do.call(pmax, as.data.frame(x))
  out <- mx
  finite <- is.finite(mx)
  out[!finite] <- -Inf
  if (any(finite))
    out[finite] <- mx[finite] + log(rowSums(exp(x[finite, , drop = FALSE] -
      mx[finite])))
  out
}

.log_prob_count <- function(probability, count) {
  if (!count) return(0)
  if (!is.finite(probability) || probability <= 0) return(-Inf)
  count * log(probability)
}

.hybrid_tenure_group_loglik <- function(
    g, offsets, lambda_g, beta_g, p_clean, p_local, p_gross,
    support_months = .TENURE_LOCAL_MONTHS,
    decay_months = .TENURE_LOCAL_DECAY_MONTHS, tol = 1e-7) {
  K <- ncol(g); n <- nrow(g)
  if (!K) return(numeric(n))
  support_years <- support_months / 12
  log_q <- .tenure_local_logweights(support_months, decay_months)
  types <- as.matrix(expand.grid(rep(list(0:2), K)))
  # 0 = clean, 1 = local, 2 = gross
  terms <- matrix(-Inf, n, nrow(types))
  fg <- lapply(seq_len(K), function(j)
    .log_duration_density(g[, j], lambda_g, beta_g))
  delta <- .QUARTER_YEARS

  for (r in seq_len(nrow(types))) {
    type <- types[r, ]
    clean <- which(type == 0L)
    local <- which(type == 1L)
    gross <- which(type == 2L)
    log_mix <- .log_prob_count(p_clean, length(clean)) +
      .log_prob_count(p_local, length(local)) +
      .log_prob_count(p_gross, length(gross))
    if (!is.finite(log_mix)) next
    gross_term <- if (length(gross))
      Reduce(`+`, fg[gross]) else rep(0, n)

    if (!length(clean) && !length(local)) {
      terms[, r] <- log_mix + gross_term
      next
    }

    if (length(clean)) {
      anchor <- clean[1L]
      T0 <- g[, anchor] - offsets[anchor] * delta
      lp <- log_mix + gross_term +
        .log_duration_density(T0, lambda_g, beta_g)
      valid <- is.finite(lp) & T0 > 0
      if (length(clean) > 1L) for (j in clean[-1L])
        valid <- valid & abs(g[, j] - (T0 + offsets[j] * delta)) < tol
      if (length(local)) for (j in local)
        lp <- lp + .tenure_local_match_logprob(
          g[, j] - (T0 + offsets[j] * delta), support_months, log_q, tol)
      lp[!valid] <- -Inf
      terms[, r] <- lp
      next
    }

    # No clean report: integrate the latent spell start by summing over the
    # finite error support for the first local report.
    anchor <- local[1L]
    candidates <- matrix(-Inf, n, length(support_years))
    for (q in seq_along(support_years)) {
      T0 <- g[, anchor] - offsets[anchor] * delta - support_years[q]
      lp <- log_mix + gross_term + log_q[q] +
        .log_duration_density(T0, lambda_g, beta_g)
      if (length(local) > 1L) for (j in local[-1L])
        lp <- lp + .tenure_local_match_logprob(
          g[, j] - (T0 + offsets[j] * delta), support_months, log_q, tol)
      lp[!is.finite(T0) | T0 <= 0] <- -Inf
      candidates[, q] <- lp
    }
    terms[, r] <- .row_logsumexp(candidates)
  }
  .row_logsumexp(terms)
}

log_emission_spell_g_local_gross <- function(
    g_mat, s_mat, t_offsets, lambda_g, eps_local, eps_gross,
    beta_g = 0, support_months = .TENURE_LOCAL_MONTHS,
    decay_months = .TENURE_LOCAL_DECAY_MONTHS, tol = 1e-7) {
  if (!is.matrix(g_mat) || !is.matrix(s_mat) ||
      !identical(dim(g_mat), dim(s_mat)) || !is.logical(s_mat))
    stop("g_mat and logical s_mat must be conformable matrices")
  if (length(t_offsets) != ncol(g_mat) || any(t_offsets < 0L))
    stop("Invalid tenure-clock offsets")
  if (any(!is.finite(lambda_g)) || any(lambda_g <= 0))
    stop("Tenure hazards must be positive")
  if (!is.finite(eps_local) || !is.finite(eps_gross) ||
      eps_local < 0 || eps_gross < 0 || eps_local + eps_gross >= 1)
    stop("Local and gross tenure-error probabilities must be nonnegative and sum below one")
  if (any(!is.finite(g_mat[s_mat])))
    stop("Observed tenure values must be finite")

  N <- nrow(g_mat); K <- as.integer(rowSums(s_mat))
  out <- numeric(N)
  code <- as.vector(s_mat %*% (2^(seq_len(ncol(s_mat)) - 1L)))
  for (key in sort(unique(code[code > 0]))) {
    rows <- which(code == key)
    cols <- which(as.logical(intToBits(key)[seq_len(ncol(s_mat))]))
    out[rows] <- .hybrid_tenure_group_loglik(
      g_mat[rows, cols, drop = FALSE], t_offsets[cols], lambda_g, beta_g,
      1 - eps_local - eps_gross, eps_local, eps_gross,
      support_months, decay_months, tol)
  }
  # Direct-likelihood estimation does not use EM sufficient statistics, but
  # return the common structure expected by e_step_eps().
  list(loglik = out, tau_sum = numeric(N), K = K,
       eps_informative = rep(FALSE, N), lambda_count = numeric(N),
       lambda_xsum = numeric(N))
}

.piecewise_hybrid_unpack <- function(z, timegap_contamination = TRUE) {
  base_names <- setdiff(names(z), "tenure_local_share")
  p <- .piecewise_eps_unpack(z[base_names],
    timegap_contamination = timegap_contamination)
  share <- plogis(unname(z["tenure_local_share"]))
  p$eps_local <- p$eps * share
  p$eps_gross <- p$eps * (1 - share)
  p$tenure_contamination_model <- "local_gross"
  p$tenure_local_support_months <- .TENURE_LOCAL_MONTHS
  p$tenure_local_decay_months <- .TENURE_LOCAL_DECAY_MONTHS
  p
}

.piecewise_hybrid_pack <- function(params, timegap_contamination = TRUE) {
  total <- params$eps
  local <- if (is.null(params$eps_local)) 0.25 * total else params$eps_local
  share <- min(max(local / total, 1e-6), 1 - 1e-6)
  c(.piecewise_eps_pack(params,
      timegap_contamination = timegap_contamination),
    tenure_local_share = qlogis(share))
}

fit_eps_piecewise_hybrid_tenure <- function(
    df, start, maxit = 500L, reltol = 1e-9, pgtol = 1e-7,
    method = "L-BFGS-B", verbose = 0L,
    timegap_contamination_model = "joint_marginal") {
  validate_df_eps(df)
  z0 <- .piecewise_hybrid_pack(start, timegap_contamination = TRUE)
  lower <- rep(-Inf, length(z0)); upper <- rep(Inf, length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- c(paste0("log_hg", 1:5), paste0("log_hd", 1:5))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  lower["tenure_local_share"] <- -12
  upper["tenure_local_share"] <- 12
  total_weight <- sum(df$weight)
  objective <- function(z) {
    p <- .piecewise_hybrid_unpack(z, timegap_contamination = TRUE)
    p$timegap_contamination_model <- timegap_contamination_model
    value <- tryCatch(e_step_eps(df, p, check_df = FALSE,
      suff_stats = FALSE)$loglik, error = function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value / total_weight
  }
  control <- list(maxit = maxit, reltol = reltol,
    ndeps = rep(1e-3, length(z0)), trace = as.integer(verbose))
  if (identical(method, "L-BFGS-B")) {
    control$factr <- reltol / .Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0, objective, method = method, lower = lower, upper = upper,
      control = control)
  } else {
    opt <- optim(z0, objective, method = method, control = control)
  }
  params <- .piecewise_hybrid_unpack(opt$par, timegap_contamination = TRUE)
  params$timegap_contamination_model <- timegap_contamination_model
  estep <- e_step_eps(df, params, check_df = FALSE, suff_stats = FALSE)
  list(params = params, loglik = estep$loglik, gamma = estep$gamma,
    convergence = opt$convergence, message = opt$message,
    iterations = unname(opt$counts["function"]), objective = opt$value,
    par_unconstrained = opt$par, objective_function = objective)
}
