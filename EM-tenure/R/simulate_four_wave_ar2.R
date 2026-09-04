# Sequential AR2 simulator derived from simulate_eps_piecewise_job_change.
# Reporting/duration draws are identical; only the transition draw is extended.
# A seeded zero-effect regression test guards synchronization with AR1.
simulate_four_wave_ar2 <- function(n, params, seed = NULL,
    waves = 4L, exact_risk = FALSE) {
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
      pieces[["reliable"]] <- simulate_four_wave_ar2(
        sum(!unreliable),.duration_reliability_component_params(params,
          "reliable"),seed=component_seeds[1L],waves=waves,exact_risk=exact_risk)
      pieces[["reliable"]]$duration_reliability_class <- "reliable"
    }
    if (any(unreliable)) {
      pieces[["unreliable"]] <- simulate_four_wave_ar2(
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
    # AR2 changes only the two later transitions, conditional on latent history.
    if (tt >= 2L) {
      recent <- h[,tt] != h[,tt-1L]
      effect <- ifelse(latent_emp, params$ar2_exit %||% 0, params$ar2_entry %||% 0)
      risk[recent] <- plogis(qlogis(risk[recent])+effect[recent])
    }
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
