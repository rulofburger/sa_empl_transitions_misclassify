# Current-job changes within a latent employment spell.
#
# A maximal latent employment spell may contain a job-to-job transition at
# either quarterly boundary. A reset partitions the employment spell into
# current-job spells without changing latent employment status. The initial
# job spell retains the existing marginal duration density. A job spell that
# begins at an E->E reset has clean initial tenure U in (0, Delta], with a
# uniform conditional event-time density; subsequent clean reports advance by
# Delta. Gross tenure contamination remains the existing marginal draw.

.job_segment_emission <- function(
    g_mat, s_mat, t_offsets, lambda_g, eps, beta_g = 0,
    start_model = c("marginal", "within_interval"), tol = 1e-7) {
  start_model <- match.arg(start_model)
  if (!is.matrix(g_mat) || !is.matrix(s_mat) ||
      !identical(dim(g_mat), dim(s_mat)) || !is.logical(s_mat))
    stop("g_mat and logical s_mat must be conformable matrices")
  if (length(t_offsets) != ncol(g_mat))
    stop("t_offsets must have length ncol(g_mat)")
  if (!is.finite(eps) || eps <= 0 || eps >= 1)
    stop("eps must be in (0,1)")

  n <- nrow(g_mat)
  K <- as.integer(rowSums(s_mat))
  ans <- list(loglik = numeric(n), tau_sum = numeric(n), K = K,
    eps_informative = logical(n), lambda_count = numeric(n),
    lambda_xsum = numeric(n))
  if (!n || !ncol(g_mat)) return(ans)

  code <- as.vector(s_mat %*% (2^(seq_len(ncol(s_mat)) - 1L)))
  for (key in sort(unique(code[code > 0]))) {
    rows <- which(code == key)
    cols <- which(as.logical(intToBits(key)[seq_len(ncol(s_mat))]))
    k <- length(cols)
    types <- as.matrix(expand.grid(rep(list(0:1), k)))
    # 0 = clean, 1 = gross marginal contamination
    lp <- matrix(-Inf, length(rows), nrow(types))
    tau <- matrix(0, length(rows), nrow(types))
    lc <- matrix(0, length(rows), nrow(types))
    lx <- matrix(0, length(rows), nrow(types))
    fg <- lapply(cols, function(j)
      .log_duration_density(g_mat[rows, j], lambda_g, beta_g))

    for (r in seq_len(nrow(types))) {
      contaminated <- which(types[r, ] == 1L)
      clean <- which(types[r, ] == 0L)
      n_x <- length(contaminated)
      log_mix <- n_x * log(eps) + (k - n_x) * log1p(-eps)
      gross_term <- if (n_x) Reduce(`+`, fg[contaminated]) else
        rep(0, length(rows))
      gross_xsum <- if (n_x) rowSums(g_mat[rows,
        cols[contaminated], drop = FALSE]) else rep(0, length(rows))

      if (!length(clean)) {
        lp[, r] <- log_mix + gross_term
        tau[, r] <- n_x
        lc[, r] <- n_x
        lx[, r] <- gross_xsum
        next
      }

      anchor <- clean[1L]
      U <- g_mat[rows, cols[anchor]] - t_offsets[cols[anchor]] *
        .QUARTER_YEARS
      valid <- is.finite(U) & U > 0
      if (length(clean) > 1L) for (j in clean[-1L])
        valid <- valid & abs(g_mat[rows, cols[j]] -
          (U + t_offsets[cols[j]] * .QUARTER_YEARS)) < tol

      if (identical(start_model, "within_interval")) {
        valid <- valid & U <= .QUARTER_YEARS + tol
        base <- rep(-log(.QUARTER_YEARS), length(rows))
        base_count <- 0
        base_xsum <- 0
      } else {
        base <- .log_duration_density(U, lambda_g, beta_g)
        base_count <- 1
        base_xsum <- pmax(U, 0)
      }
      value <- log_mix + gross_term + base
      value[!valid] <- -Inf
      lp[, r] <- value
      tau[, r] <- n_x
      lc[, r] <- n_x + base_count
      lx[, r] <- gross_xsum + base_xsum
    }

    ll <- .row_logsumexp(lp)
    post <- exp(lp - ll)
    ans$loglik[rows] <- ll
    ans$tau_sum[rows] <- rowSums(post * tau)
    ans$lambda_count[rows] <- rowSums(post * lc)
    ans$lambda_xsum[rows] <- rowSums(post * lx)
    ans$eps_informative[rows] <- if (identical(start_model, "marginal"))
      k >= 2L || any(t_offsets[cols] > 0L) else TRUE
  }
  ans
}

.job_change_partitions <- function(length_spell) {
  if (length_spell <= 1L) return(list(integer(0)))
  grid <- expand.grid(rep(list(0:1), length_spell - 1L))
  lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, ]))
}

log_emission_spell_g_job_change <- function(
    g_mat, s_mat, t_offsets, lambda_g, eps, job_change_prob,
    beta_g = 0, tol = 1e-7) {
  if (!is.finite(job_change_prob) || job_change_prob < 0 ||
      job_change_prob >= 1)
    stop("job_change_prob must be in [0,1)")
  L <- ncol(g_mat)
  if (L <= 1L || job_change_prob == 0) {
    out <- log_emission_spell_g(g_mat, s_mat, t_offsets, lambda_g, eps,
      tol = tol, beta_g = beta_g)
    out$job_changes <- numeric(nrow(g_mat))
    out$job_change_opportunities <- rep(max(L - 1L, 0L), nrow(g_mat))
    return(out)
  }

  patterns <- .job_change_partitions(L)
  n <- nrow(g_mat)
  lp <- matrix(-Inf, n, length(patterns))
  tau <- lc <- lx <- matrix(0, n, length(patterns))
  changes <- numeric(length(patterns))

  for (r in seq_along(patterns)) {
    reset <- patterns[[r]]
    changes[r] <- sum(reset)
    lp[, r] <- changes[r] * log(job_change_prob) +
      (length(reset) - changes[r]) * log1p(-job_change_prob)
    starts <- c(1L, which(reset == 1L) + 1L)
    ends <- c(which(reset == 1L), L)
    for (b in seq_along(starts)) {
      cols <- starts[b]:ends[b]
      segment <- .job_segment_emission(
        g_mat[, cols, drop = FALSE], s_mat[, cols, drop = FALSE],
        as.integer(seq_along(cols) - 1L), lambda_g, eps, beta_g,
        start_model = if (b == 1L) "marginal" else "within_interval",
        tol = tol)
      lp[, r] <- lp[, r] + segment$loglik
      tau[, r] <- tau[, r] + segment$tau_sum
      lc[, r] <- lc[, r] + segment$lambda_count
      lx[, r] <- lx[, r] + segment$lambda_xsum
    }
  }

  ll <- .row_logsumexp(lp)
  post <- exp(lp - ll)
  list(loglik = ll, tau_sum = rowSums(post * tau),
    K = as.integer(rowSums(s_mat)),
    eps_informative = rowSums(s_mat) >= 2L,
    lambda_count = rowSums(post * lc),
    lambda_xsum = rowSums(post * lx),
    job_changes = as.vector(post %*% changes),
    job_change_opportunities = rep(L - 1L, n))
}

# Simulate from the intended factorisation of the piecewise job-change model.
# The transition risk at each boundary is evaluated from the duration report
# available to the likelihood (or its model-implied marginal risk when that
# state-specific duration is unobserved). This makes the simulator suitable
# for end-to-end recovery diagnostics of the implemented observed likelihood.
# A recovery failure is informative: in particular, it can reveal a mismatch
# between the intended discrete reporting process and the continuous-density
# approximation used by the estimator.
.duration_reliability_shifts <- function(params) {
  common <- if (is.null(params$duration_reliability_shift)) 0 else
    params$duration_reliability_shift
  tenure <- if (is.null(params$tenure_reliability_shift)) common else
    params$tenure_reliability_shift
  timegap <- if (is.null(params$timegap_reliability_shift)) common else
    params$timegap_reliability_shift
  shifts <- c(tenure=unname(tenure),timegap=unname(timegap))
  if (any(!is.finite(shifts)) || any(shifts<0))
    stop("duration reliability shifts must be finite and nonnegative")
  shifts
}

.duration_reliability_component_params <- function(params,
    class=c("reliable","unreliable")) {
  class <- match.arg(class)
  shifts <- .duration_reliability_shifts(params)
  direction <- if (identical(class,"reliable")) -1 else 1
  shift_probability <- function(probability,shift) {
    if (!is.finite(probability) || probability<0 || probability>1)
      stop("duration contamination probabilities must be in [0,1]")
    if (probability==0 || probability==1) return(probability)
    plogis(qlogis(probability)+direction*shift)
  }
  out <- params
  out$eps <- shift_probability(params$eps,unname(shifts["tenure"]))
  out$eps_d <- shift_probability(params$eps_d,unname(shifts["timegap"]))
  out$duration_reliability_shift <- 0
  out$tenure_reliability_shift <- 0
  out$timegap_reliability_shift <- 0
  out$duration_reliability_class <- class
  out
}

simulate_eps_piecewise_job_change <- function(n, params, seed = NULL,
    waves = 3L, exact_risk = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  if (length(n) != 1L || !is.finite(n) || n < 1) stop("n must be positive")
  n <- as.integer(n)
  if (length(waves) != 1L || !waves %in% c(3L, 4L))
    stop("waves must be 3 or 4")
  if (exact_risk && (length(params$lambda_g) != 5L ||
      length(params$lambda_d) != 5L)) stop("Exact risks require piecewise hazards")
  required <- c("alpha", "pi", "eps", "eps_d", "job_change_prob",
    "lambda_g", "lambda_d")
  absent <- required[!vapply(required, function(z) !is.null(params[[z]]),
    logical(1))]
  if (length(absent)) stop("Missing simulation parameters: ",
    paste(absent, collapse=", "))
  reliability_shifts <- .duration_reliability_shifts(params)
  if (any(reliability_shifts>0)) {
    unreliable <- runif(n)<.5
    component_seeds <- sample.int(.Machine$integer.max,2L)
    pieces <- list()
    if (any(!unreliable)) {
      pieces[["reliable"]] <- simulate_eps_piecewise_job_change(
        sum(!unreliable),.duration_reliability_component_params(params,
          "reliable"),seed=component_seeds[1L],waves=waves,exact_risk=exact_risk)
      pieces[["reliable"]]$duration_reliability_class <- "reliable"
    }
    if (any(unreliable)) {
      pieces[["unreliable"]] <- simulate_eps_piecewise_job_change(
        sum(unreliable),.duration_reliability_component_params(params,
          "unreliable"),seed=component_seeds[2L],waves=waves,exact_risk=exact_risk)
      pieces[["unreliable"]]$duration_reliability_class <- "unreliable"
    }
    result <- do.call(rbind,pieces)
    result <- result[sample.int(nrow(result)),,drop=FALSE]
    rownames(result) <- NULL
    return(result)
  }
  draw_duration <- function(lambda, size) {
    y <- rexp(size)
    .duration_inverse_cumhaz(y, lambda, 0)
  }
  monthly_tenure <- identical(params$tenure_measurement_model, "monthly")
  tenure_persistence <- if (is.null(params$tenure_report_persistence)) 0 else
    params$tenure_report_persistence
  tenure_heaping <- if (is.null(params$tenure_heaping_prob)) 0 else
    params$tenure_heaping_prob
  tenure_year_revision <- if (is.null(params$tenure_year_revision_prob)) 0 else
    params$tenure_year_revision_prob
  tenure_clean_anchor_revision <-
    if (is.null(params$tenure_clean_anchor_revision_prob)) 0 else
      params$tenure_clean_anchor_revision_prob
  tenure_exact_anchor_retention <-
    if (is.null(params$tenure_exact_anchor_retention_prob)) 0 else
      params$tenure_exact_anchor_retention_prob
  tenure_local_revision <-
    if (is.null(params$tenure_local_revision_prob)) 0 else
      params$tenure_local_revision_prob
  tenure_start_month_probs <- .validate_start_month_probs(
    params$tenure_start_month_probs)
  if (!is.finite(tenure_persistence) || tenure_persistence < 0 ||
      tenure_persistence >= 1)
    stop("tenure_report_persistence must be in [0,1)")
  if (tenure_persistence > 0 && !monthly_tenure)
    stop("tenure-report persistence requires monthly tenure")
  if (!is.finite(tenure_heaping) || tenure_heaping < 0 ||
      tenure_heaping >= 1)
    stop("tenure_heaping_prob must be in [0,1)")
  if (tenure_heaping > 0 && !monthly_tenure)
    stop("tenure heaping requires monthly tenure")
  if (tenure_heaping > 0 && tenure_persistence > 0)
    stop("persistence and calendar heaping are not jointly implemented")
  if (!is.finite(tenure_year_revision) || tenure_year_revision < 0 ||
      tenure_year_revision >= 1)
    stop("tenure_year_revision_prob must be in [0,1)")
  if (tenure_year_revision > 0 && !monthly_tenure)
    stop("whole-year revisions require monthly tenure")
  if (tenure_year_revision > 0 && tenure_persistence > 0)
    stop("persistence and whole-year revisions are not jointly implemented")
  if (!is.finite(tenure_clean_anchor_revision) ||
      tenure_clean_anchor_revision < 0 || tenure_clean_anchor_revision >= 1)
    stop("tenure_clean_anchor_revision_prob must be in [0,1)")
  if (tenure_clean_anchor_revision > 0 && !monthly_tenure)
    stop("clean-anchor revisions require monthly tenure")
  if (tenure_clean_anchor_revision > 0 && tenure_persistence > 0)
    stop("persistence and clean-anchor revisions are not jointly implemented")
  if (!is.finite(tenure_exact_anchor_retention) ||
      tenure_exact_anchor_retention < 0 ||
      tenure_exact_anchor_retention >= 1)
    stop("tenure_exact_anchor_retention_prob must be in [0,1)")
  if (tenure_exact_anchor_retention > 0 && !monthly_tenure)
    stop("exact-anchor retention requires monthly tenure")
  if (tenure_exact_anchor_retention > 0 && tenure_persistence > 0)
    stop("persistence and exact-anchor retention are not jointly implemented")
  if (!is.finite(tenure_local_revision) || tenure_local_revision < 0 ||
      tenure_local_revision >= 1)
    stop("tenure_local_revision_prob must be in [0,1)")
  if (tenure_local_revision > 0 && !monthly_tenure)
    stop("local tenure revisions require monthly tenure")
  if (tenure_local_revision > 0 && tenure_persistence > 0)
    stop("persistence and local revisions are not jointly implemented")
  report_tenure <- function(x) if (monthly_tenure) floor(12 * x) / 12 else x
  draw_tenure <- function(lambda, size)
    report_tenure(draw_duration(lambda, size))
  draw_calendar_tenure <- function(lambda,interview_month) {
    if (!monthly_tenure ||
        max(abs(tenure_start_month_probs-1/12))<=1e-12)
      return(draw_tenure(lambda,length(interview_month)))
    out <- rep(NA_integer_,length(interview_month))
    pending <- seq_along(out)
    ceiling_prob <- max(tenure_start_month_probs)
    while (length(pending)) {
      candidate <- floor(12*draw_duration(lambda,length(pending)))
      start_idx <- ((as.integer(interview_month[pending])-1L-
        candidate%%12L)%%12L)+1L
      accept <- runif(length(pending)) <
        tenure_start_month_probs[start_idx]/ceiling_prob
      out[pending[accept]] <- candidate[accept]
      pending <- pending[!accept]
    }
    out/12
  }
  draw_january_tenure <- function(lambda, interview_month) {
    out <- rep(NA_integer_, length(interview_month))
    pending <- seq_along(out)
    while (length(pending)) {
      candidate <- floor(12*draw_duration(lambda,length(pending)))
      accept <- candidate %% 12L == as.integer(interview_month[pending])-1L
      out[pending[accept]] <- candidate[accept]
      pending <- pending[!accept]
    }
    out/12
  }
  draw_gross_tenure <- function(lambda, interview_month) {
    value <- draw_calendar_tenure(lambda,interview_month)
    heaped <- runif(length(interview_month)) < tenure_heaping
    if (any(heaped)) value[heaped] <- draw_january_tenure(lambda,
      interview_month[heaped])
    value
  }
  draw_year_revision_tenure <- function(
      lambda, interview_month, continuation) {
    out <- rep(NA_integer_, length(interview_month))
    continuation_month <- as.integer(round(12 * continuation))
    target_residue <- continuation_month %% 12L
    pending <- seq_along(out)
    while (length(pending)) {
      candidate <- floor(12 * draw_duration(lambda, length(pending)))
      accept <- candidate %% 12L == target_residue[pending] &
        candidate != continuation_month[pending]
      out[pending[accept]] <- candidate[accept]
      pending <- pending[!accept]
    }
    out / 12
  }
  draw_local_revision_tenure <- function(continuation) {
    continuation_month <- as.integer(round(12*continuation))
    vapply(continuation_month,function(anchor) {
      feasible <- anchor+.TENURE_LOCAL_MONTHS>=0L
      support <- .TENURE_LOCAL_MONTHS[feasible]
      weights <- exp(.tenure_local_logweights()[feasible])
      anchor+sample(support,1L,prob=weights)
    },integer(1))/12
  }
  draw_cat <- function(lambda, size)
    .continuous_to_cat(draw_duration(lambda, size))

  h <- y <- matrix(NA_integer_, n, waves)
  interview_month <- matrix(NA_integer_,n,waves)
  interview_month[,1L] <- sample(c(3L,6L,9L,12L),n,replace=TRUE)
  for (tt in 2:waves) interview_month[,tt] <-
    (interview_month[,tt-1L]+2L) %% 12L+1L
  tenure <- matrix(NA_real_, n, waves)
  timegap <- matrix(NA_integer_, n, waves)
  reset <- matrix(FALSE, n, waves-1L)
  clean_g <- clean_d <- rep(NA_real_, n)
  gross_active <- rep(FALSE, n)
  gross_value <- rep(NA_real_, n)
  gross_wave <- rep(NA_integer_, n)
  anchor_was_gross <- rep(FALSE,n)
  h[, 1L] <- as.integer(runif(n) < params$alpha)
  emp <- h[, 1L] == 1L
  clean_g[emp] <- draw_calendar_tenure(params$lambda_g,
    interview_month[emp,1L])
  clean_d[!emp] <- draw_duration(params$lambda_d, sum(!emp))
  p_exit_missing <- if (exact_risk)
    .piecewise_mean_transition_risk_exact(params$lambda_g) else
    .duration_marginal_transition_probability(params$lambda_g, 0)
  p_entry_missing <- if (exact_risk)
    .piecewise_mean_transition_risk_exact(params$lambda_d) else
    .duration_marginal_transition_probability(params$lambda_d, 0)

  for (tt in seq_len(waves)) {
    latent_emp <- h[, tt] == 1L
    y[, tt] <- ifelse(runif(n) < params$pi, 1L - h[, tt], h[, tt])
    reported_emp <- y[, tt] == 1L

    both_emp <- latent_emp & reported_emp
    if (any(both_emp)) {
      idx <- which(both_emp)
      gross <- runif(length(idx)) < params$eps
      value <- clean_g[idx]
      if (any(gross)) {
        gi <- idx[gross]
        retain_prob <- numeric(length(gi))
        if (tenure_exact_anchor_retention > 0) {
          can_retain <- gross_active[gi] & is.finite(gross_value[gi]) &
            !is.na(gross_wave[gi])
          retain_prob[can_retain] <- tenure_exact_anchor_retention
        } else {
          can_retain <- gross_active[gi] & anchor_was_gross[gi] &
            !is.na(gross_wave[gi])
          retain_prob[can_retain] <- tenure_persistence^(
            tt-gross_wave[gi[can_retain]])
        }
        retain <- runif(length(gi)) < retain_prob
        if (any(retain)) value[which(gross)[retain]] <- gross_value[gi[retain]]
        local <- rep(FALSE,length(gi))
        can_local <- !retain & gross_active[gi] & is.finite(gross_value[gi])
        local[can_local] <- runif(sum(can_local)) < tenure_local_revision
        if (any(local)) value[which(gross)[local]] <-
          draw_local_revision_tenure(gross_value[gi[local]])
        revise <- rep(FALSE, length(gi))
        can_revise_gross <- !retain & !local & gross_active[gi] &
          anchor_was_gross[gi] & is.finite(gross_value[gi])
        revise[can_revise_gross] <- runif(sum(can_revise_gross)) <
          tenure_year_revision
        can_revise_clean <- !retain & !local & !revise & gross_active[gi] &
          !anchor_was_gross[gi] & is.finite(gross_value[gi])
        revise[can_revise_clean] <- runif(sum(can_revise_clean)) <
          tenure_clean_anchor_revision
        if (any(revise)) value[which(gross)[revise]] <-
          draw_year_revision_tenure(params$lambda_g,
            interview_month[gi[revise],tt],gross_value[gi[revise]])
        redraw <- !retain & !local & !revise
        if (any(redraw)) value[which(gross)[redraw]] <-
          draw_gross_tenure(params$lambda_g,interview_month[gi[redraw],tt])
        gross_active[gi] <- TRUE
        gross_value[gi] <- value[which(gross)]
        gross_wave[gi] <- tt
        anchor_was_gross[gi] <- TRUE
      }
      if (any(!gross)) {
        clean_idx <- idx[!gross]
        gross_active[clean_idx] <- TRUE
        gross_value[clean_idx] <- clean_g[clean_idx]
        gross_wave[clean_idx] <- tt
        anchor_was_gross[clean_idx] <- FALSE
      }
      tenure[idx, tt] <- value
    }
    false_emp <- !latent_emp & reported_emp
    if (any(false_emp)) {
      tenure[false_emp, tt] <- draw_gross_tenure(params$lambda_g,
        interview_month[false_emp,tt])
      gross_active[false_emp] <- FALSE
      gross_value[false_emp] <- NA_real_
      gross_wave[false_emp] <- NA_integer_
      anchor_was_gross[false_emp] <- FALSE
    }

    both_non <- !latent_emp & !reported_emp
    if (any(both_non)) {
      gross <- runif(sum(both_non)) < params$eps_d
      value <- .continuous_to_cat(clean_d[both_non])
      value[gross] <- draw_cat(params$lambda_d, sum(gross))
      timegap[both_non, tt] <- value
    }
    false_non <- latent_emp & !reported_emp
    if (any(false_non)) timegap[false_non, tt] <-
      draw_cat(params$lambda_d, sum(false_non))

    if (tt == waves) next
    risk <- numeric(n)
    use_g <- latent_emp & reported_emp
    risk[use_g] <- .duration_transition_probability(tenure[use_g, tt],
      params$lambda_g, 0)
    risk[latent_emp & !reported_emp] <- p_exit_missing
    use_d <- !latent_emp & !reported_emp
    risk[use_d] <- if (exact_risk)
      .piecewise_category_transition_risk_exact(timegap[use_d, tt], params$lambda_d)
      else .duration_category_transition_probability(
        timegap[use_d, tt], params$lambda_d, 0)
    risk[!latent_emp & reported_emp] <- p_entry_missing
    change <- runif(n) < risk
    h[, tt + 1L] <- ifelse(change, 1L - h[, tt], h[, tt])
    next_emp <- h[, tt + 1L] == 1L

    ee <- latent_emp & next_emp
    reset[ee, tt] <- runif(sum(ee)) < params$job_change_prob
    same_job <- ee & !reset[, tt]
    new_job <- ee & reset[, tt]
    clean_g[same_job] <- clean_g[same_job] + .QUARTER_YEARS
    gross_value[same_job & gross_active] <-
      gross_value[same_job & gross_active] + .QUARTER_YEARS
    draw_reset_month <- function(idx) {
      if (!length(idx)) return(numeric())
      vapply(idx,function(ii) {
        support <- .RESET_TENURE_MONTHS
        start_idx <- ((interview_month[ii,tt+1L]-1L-support)%%12L)+1L
        sample(support,1L,prob=tenure_start_month_probs[start_idx])
      },numeric(1))/12
    }
    clean_g[new_job] <- if (monthly_tenure)
      draw_reset_month(which(new_job)) else
      runif(sum(new_job), 1e-8, .QUARTER_YEARS)
    enter_emp <- !latent_emp & next_emp
    clean_g[enter_emp] <- if (monthly_tenure)
      draw_reset_month(which(enter_emp)) else
      draw_tenure(params$lambda_g, sum(enter_emp))
    clean_g[!next_emp] <- NA_real_
    break_anchor <- !next_emp | new_job | enter_emp
    gross_active[break_anchor] <- FALSE
    gross_value[break_anchor] <- NA_real_
    gross_wave[break_anchor] <- NA_integer_
    anchor_was_gross[break_anchor] <- FALSE

    nn <- !latent_emp & !next_emp
    clean_d[nn] <- clean_d[nn] + .QUARTER_YEARS
    enter_non <- latent_emp & !next_emp
    clean_d[enter_non] <- draw_duration(params$lambda_d, sum(enter_non))
    clean_d[next_emp] <- NA_real_
  }

  result <- data.frame(y1=y[, 1], y2=y[, 2], y3=y[, 3],
    tenure1=tenure[, 1], tenure2=tenure[, 2], tenure3=tenure[, 3],
    timegap_cat1=timegap[, 1], timegap_cat2=timegap[, 2],
    timegap_cat3=timegap[, 3],
    interview_month1=interview_month[,1],
    interview_month2=interview_month[,2],
    interview_month3=interview_month[,3],weight=1,
    h1=h[, 1], h2=h[, 2], h3=h[, 3],
    reset12=reset[, 1], reset23=reset[, 2])
  if (waves == 4L) {
    result$y4 <- y[,4L]
    result$tenure4 <- tenure[,4L]
    result$timegap_cat4 <- timegap[,4L]
    result$interview_month4 <- interview_month[,4L]
    result$h4 <- h[,4L]
    result$reset34 <- reset[,3L]
  }
  result
}

.piecewise_job_change_unpack <- function(
    z, q_cap = .50, timegap_contamination_model = "joint_marginal") {
  base_names <- setdiff(names(z), "job_change")
  p <- .piecewise_eps_unpack(z[base_names], timegap_contamination = TRUE)
  p$job_change_prob <- q_cap * plogis(unname(z["job_change"]))
  p$job_change_model <- "within_interval_uniform"
  p$timegap_contamination_model <- timegap_contamination_model
  p
}

.piecewise_job_change_pack <- function(params, q_start = NULL, q_cap = .50) {
  q <- if (is.null(q_start)) {
    if (is.null(params$job_change_prob)) .02 else params$job_change_prob
  } else q_start
  q <- min(max(q, 1e-8), q_cap - 1e-8)
  c(.piecewise_eps_pack(params, timegap_contamination = TRUE),
    job_change = qlogis(q / q_cap))
}

fit_eps_piecewise_job_change <- function(
    df, start, q_start = NULL, maxit = 400L, reltol = 1e-9,
    pgtol = 1e-7, method = "L-BFGS-B", verbose = 0L,
    timegap_contamination_model = "joint_marginal", workers = 1L,
    gradient_step = 1e-4, q_fixed = NULL,
    tenure_measurement_model = c("continuous", "monthly")) {
  tenure_measurement_model <- match.arg(tenure_measurement_model)
  validate_df_eps(df, allow_zero_tenure =
    identical(tenure_measurement_model, "monthly"))
  if (!is.null(q_fixed) && (!is.finite(q_fixed) || q_fixed < 0 ||
      q_fixed >= .50)) stop("q_fixed must be in [0, .50)")
  z0 <- if (is.null(q_fixed))
    .piecewise_job_change_pack(start, q_start) else
    .piecewise_eps_pack(start, timegap_contamination = TRUE)
  unpack <- function(z) {
    if (is.null(q_fixed)) {
      p <- .piecewise_job_change_unpack(z,
        timegap_contamination_model = timegap_contamination_model)
    } else {
      p <- .piecewise_eps_unpack(z, timegap_contamination = TRUE)
      p$job_change_prob <- q_fixed
      p$job_change_model <- "within_interval_uniform"
      p$timegap_contamination_model <- timegap_contamination_model
    }
    p$tenure_measurement_model <- tenure_measurement_model
    p
  }
  lower <- rep(-Inf, length(z0)); upper <- rep(Inf, length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- c(paste0("log_hg", 1:5), paste0("log_hd", 1:5))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  if ("job_change" %in% names(z0)) {
    lower["job_change"] <- -16; upper["job_change"] <- 8
  }
  total_weight <- sum(df$weight)
  objective <- function(z) {
    p <- unpack(z)
    value <- tryCatch(e_step_eps(df, p, check_df = FALSE,
      suff_stats = FALSE)$loglik, error = function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value / total_weight
  }
  cluster <- NULL
  gradient <- NULL
  workers <- as.integer(workers)
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster, function(path) {
      setwd(path)
      source("EM-tenure/R/source_all.R")
      NULL
    }, worker_wd)
    df_worker <- df
    total_weight_worker <- total_weight
    timegap_worker <- timegap_contamination_model
    tenure_measurement_worker <- tenure_measurement_model
    q_fixed_worker <- q_fixed
    parallel::clusterExport(cluster,
      c("df_worker", "total_weight_worker", "timegap_worker",
        "q_fixed_worker", "tenure_measurement_worker"),
      envir = environment())
    worker_objective <- function(z) {
      if (is.null(q_fixed_worker)) {
        p <- .piecewise_job_change_unpack(z,
          timegap_contamination_model = timegap_worker)
      } else {
        p <- .piecewise_eps_unpack(z, timegap_contamination = TRUE)
        p$job_change_prob <- q_fixed_worker
        p$job_change_model <- "within_interval_uniform"
        p$timegap_contamination_model <- timegap_worker
      }
      p$tenure_measurement_model <- tenure_measurement_worker
      value <- tryCatch(e_step_eps(df_worker, p, check_df = FALSE,
        suff_stats = FALSE)$loglik, error = function(e) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value / total_weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step * pmax(1, abs(z))
      points <- vector("list", 2L * length(z))
      for (j in seq_along(z)) {
        points[[2L*j - 1L]] <- points[[2L*j]] <- z
        points[[2L*j - 1L]][j] <- z[j] + step[j]
        points[[2L*j]][j] <- z[j] - step[j]
      }
      values <- unlist(parallel::parLapply(cluster, points,
        worker_objective), use.names = FALSE)
      plus <- values[seq(1L, length(values), by=2L)]
      minus <- values[seq(2L, length(values), by=2L)]
      centre <- objective(z)
      # Central differences are preferred. At a numerical boundary one of the
      # perturbed likelihoods can be invalid even though the current point is
      # valid; use the derivative from the valid side rather than handing
      # L-BFGS-B an infinite gradient.
      valid_plus <- is.finite(plus) & plus < 1e50
      valid_minus <- is.finite(minus) & minus < 1e50
      out <- numeric(length(z))
      both <- valid_plus & valid_minus
      out[both] <- (plus[both] - minus[both]) / (2 * step[both])
      plus_only <- valid_plus & !valid_minus
      out[plus_only] <- (plus[plus_only] - centre) / step[plus_only]
      minus_only <- !valid_plus & valid_minus
      out[minus_only] <- (centre - minus[minus_only]) / step[minus_only]
      out[!is.finite(out)] <- 0
      out
    }
  }
  control <- list(maxit = maxit, reltol = reltol,
    ndeps = rep(1e-3, length(z0)), trace = as.integer(verbose))
  if (identical(method, "L-BFGS-B")) {
    control$factr <- reltol / .Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0, objective, gr = gradient, method = method,
      lower = lower, upper = upper, control = control)
  } else {
    opt <- optim(z0, objective, gr = gradient, method = method,
      control = control)
  }
  params <- unpack(opt$par)
  estep <- e_step_eps(df, params, check_df = FALSE, suff_stats = FALSE)
  list(params = params, loglik = estep$loglik, gamma = estep$gamma,
    job_change_posterior = estep$job_change_posterior,
    convergence = opt$convergence, message = opt$message,
    iterations = unname(opt$counts["function"]), objective = opt$value,
    par_unconstrained = opt$par, objective_function = objective)
}
