# Discrete-month tenure emissions.
#
# QLFS tenure is constructed from the reported month in which the current job
# started. Its likelihood contribution is therefore a probability mass, not a
# continuous density. Clean clocks advance by three months per panel wave.

.TENURE_MONTHS_PER_WAVE <- 3L
.RESET_TENURE_MONTHS <- 0:2

.tenure_month_index <- function(g, tol = 1e-6) {
  raw <- 12 * g
  month <- round(raw)
  valid <- is.finite(raw) & month >= 0 & abs(raw - month) < tol
  list(month=as.integer(month), valid=valid)
}

.log_duration_month_mass <- function(month, lambda, beta = 0) {
  out <- rep(-Inf, length(month))
  valid <- is.finite(month) & month >= 0 &
    abs(month - round(month)) < 1e-8
  if (!any(valid)) return(out)
  a <- round(month[valid]) / 12
  b <- (round(month[valid]) + 1) / 12
  Ha <- .duration_cumhaz(a, lambda, beta)
  span <- .duration_cumhaz(b, lambda, beta) - Ha
  out[valid] <- -Ha + log(-expm1(-span))
  out
}

# A January-heaped report is a draw from the same marginal duration law,
# conditional on the implied reported job-start month being January.  QLFS
# interview months are known, so a tenure of m months is January-heaped when
# m mod 12 equals interview_month - 1.  For the piecewise-constant model the
# normalizing probability has an exact geometric tail beyond the final knot.
.log_january_duration_normalizers <- function(lambda, beta = 0) {
  if (abs(beta) > 1e-12)
    stop("January heaping currently requires piecewise-constant hazards")
  if (length(lambda) != 5L || any(!is.finite(lambda)) || any(lambda <= 0))
    stop("January heaping requires five positive duration hazards")
  out <- numeric(12L)
  for (residue in 0:11) {
    early <- seq.int(residue, 59L, by = 12L)
    early_mass <- if (length(early))
      sum(exp(.log_duration_month_mass(early, lambda, beta))) else 0
    first_tail <- residue + 12L * ceiling((60L - residue) / 12L)
    first_mass <- exp(.log_duration_month_mass(first_tail, lambda, beta))
    tail_mass <- first_mass / (-expm1(-lambda[5L]))
    out[residue + 1L] <- log(early_mass + tail_mass)
  }
  out
}

.log_january_duration_month_mass <- function(
    month, interview_month, lambda, beta = 0,
    log_normalizers = .log_january_duration_normalizers(lambda, beta)) {
  if (length(interview_month) != length(month))
    stop("interview_month and month must have equal length")
  valid_interview <- is.finite(interview_month) &
    interview_month %in% 1:12
  if (!all(valid_interview))
    stop("interview months must be integers from 1 to 12")
  base <- .log_duration_month_mass(month, lambda, beta)
  residue <- as.integer(interview_month) - 1L
  january <- is.finite(month) & month >= 0 &
    (as.integer(round(month)) %% 12L == residue)
  out <- rep(-Inf, length(month))
  out[january] <- base[january] - log_normalizers[residue[january] + 1L]
  out
}

.validate_start_month_probs <- function(probs) {
  if (is.null(probs)) return(rep(1/12,12L))
  if (length(probs)!=12L || any(!is.finite(probs)) || any(probs<=0))
    stop("start-month probabilities must be 12 finite positive values")
  probs/sum(probs)
}

# Combine the duration law with a flexible calendar distribution for the
# implied job-start month.  The multiplier is normalized separately for each
# interview month, so uniform month probabilities reproduce the unmodified
# duration mass exactly.  The same mass is used for clean and gross marginal
# reports; gross-specific January heaping remains a separate reporting branch.
.log_calendar_duration_month_mass <- function(
    month, interview_month, lambda, beta=0,
    start_month_probs=rep(1/12,12L),
    log_residue_masses=.log_january_duration_normalizers(lambda,beta)) {
  if (length(interview_month)!=length(month))
    stop("interview_month and month must have equal length")
  if (any(!interview_month %in% 1:12))
    stop("interview months must be integers from 1 to 12")
  probs <- .validate_start_month_probs(start_month_probs)
  base <- .log_duration_month_mass(month,lambda,beta)
  residue_mass <- exp(log_residue_masses)
  normalizer <- vapply(1:12,function(i) {
    start_index <- ((i-1L-(0:11))%%12L)+1L
    sum(residue_mass*12*probs[start_index])
  },numeric(1))
  valid <- is.finite(month) & month>=0 &
    abs(month-round(month))<1e-8
  out <- rep(-Inf,length(month))
  if (any(valid)) {
    residue <- as.integer(round(month[valid]))%%12L
    start_index <- ((as.integer(interview_month[valid])-1L-residue)%%12L)+1L
    out[valid] <- base[valid]+log(12*probs[start_index])-
      log(normalizer[as.integer(interview_month[valid])])
  }
  out
}

.log_within_interval_start_month_mass <- function(
    month, interview_month, start_month_probs=rep(1/12,12L)) {
  if (length(interview_month)!=length(month))
    stop("interview_month and month must have equal length")
  probs <- .validate_start_month_probs(start_month_probs)
  out <- rep(-Inf,length(month))
  valid <- is.finite(month) & month %in% .RESET_TENURE_MONTHS &
    interview_month %in% 1:12
  if (!any(valid)) return(out)
  normalizers <- vapply(1:12,function(i) {
    idx <- ((i-1L-.RESET_TENURE_MONTHS)%%12L)+1L
    sum(probs[idx])
  },numeric(1))
  possible <- normalizers[as.integer(interview_month[valid])]
  idx <- ((as.integer(interview_month[valid])-1L-
    as.integer(month[valid]))%%12L)+1L
  out[valid] <- log(probs[idx])-log(possible)
  out
}

# A whole-year revision preserves the previously reported job-start month but
# changes its year.  This is a proper conditional mass: the marginal duration
# law is restricted to the required start-month residue and the exact
# continuation of the previous reported anchor is removed from its support.
.log_whole_year_revision_month_mass <- function(
    month, previous_month, interview_month, previous_interview_month,
    gap, lambda, beta = 0,
    log_normalizers = .log_january_duration_normalizers(lambda, beta)) {
  lens <- c(length(month), length(previous_month), length(interview_month),
    length(previous_interview_month), length(gap))
  if (length(unique(lens)) != 1L)
    stop("whole-year revision inputs must have equal length")
  if (any(!interview_month %in% 1:12) ||
      any(!previous_interview_month %in% 1:12))
    stop("interview months must be integers from 1 to 12")
  continuation <- previous_month + gap * .TENURE_MONTHS_PER_WAVE
  previous_start <- (as.integer(previous_interview_month) - 1L -
    as.integer(previous_month)) %% 12L
  required_residue <- (as.integer(interview_month) - 1L -
    previous_start) %% 12L
  base <- .log_duration_month_mass(month, lambda, beta)
  continuation_mass <- .log_duration_month_mass(continuation, lambda, beta)
  log_z <- log_normalizers[required_residue + 1L]
  continuation_in_support <- continuation %% 12L == required_residue
  removed_share <- ifelse(continuation_in_support,
    pmin(exp(continuation_mass - log_z),1-1e-15),0)
  log_z_without_continuation <- log_z + log1p(-removed_share)
  valid <- is.finite(month) & month >= 0 &
    as.integer(round(month)) %% 12L == required_residue &
    as.integer(round(month)) != as.integer(round(continuation))
  out <- rep(-Inf, length(month))
  out[valid] <- base[valid] - log_z_without_continuation[valid]
  out
}

.log_probability_mixture <- function(log_a, log_b, probability) {
  if (probability == 0) return(log_a)
  if (probability == 1) return(log_b)
  top <- pmax(log1p(-probability) + log_a, log(probability) + log_b)
  out <- top + log(exp(log1p(-probability) + log_a - top) +
    exp(log(probability) + log_b - top))
  out[!is.finite(top)] <- -Inf
  out
}

# A normalized local revision changes the implied reported start date by one
# to six months around the preceding anchor. Zero is excluded because exact
# retention is a separate channel. The support is truncated and renormalized
# when a negative current tenure would otherwise result.
.log_local_revision_month_mass <- function(current_month, previous_month, gap,
    support_months=.TENURE_LOCAL_MONTHS,
    decay_months=.TENURE_LOCAL_DECAY_MONTHS) {
  if (length(gap)==1L) gap <- rep(gap,length(current_month))
  if (length(current_month)!=length(previous_month) ||
      length(gap)!=length(current_month))
    stop("Local-revision inputs must have equal lengths")
  logweights <- .tenure_local_logweights(support_months,decay_months)
  expected <- previous_month+gap*.TENURE_MONTHS_PER_WAVE
  delta <- current_month-expected
  out <- rep(-Inf,length(current_month))
  matched <- match(delta,support_months)
  valid <- is.finite(current_month) & is.finite(expected) &
    !is.na(matched) & current_month>=0L
  if (!any(valid)) return(out)

  # Only expected tenures below the largest backward revision require a
  # truncated normalizer. Compute those few denominators once rather than
  # looping over every collapsed cell in every likelihood evaluation.
  logsum <- function(x) {
    maximum <- max(x)
    maximum+log(sum(exp(x-maximum)))
  }
  normalizer <- rep(logsum(logweights),length(out))
  cutoff <- max(0L,ceiling(-min(support_months))-1L)
  for (value in 0L:cutoff) {
    rows <- is.finite(expected) & expected==value
    if (any(rows))
      normalizer[rows] <- logsum(logweights[value+support_months>=0L])
  }
  out[valid] <- logweights[matched[valid]]-normalizer[valid]
  out
}

# Gross tenure reports normally draw a fresh reported start-date anchor from
# the marginal duration distribution. The legacy persistence branch can link
# uninterrupted gross runs. The exact-retention branch instead lets a gross
# current report retain either a clean or gross preceding observed anchor, so
# tenure advances by exactly three months per intervening wave. Both branches
# reduce exactly to independent gross draws at zero.
.gross_month_emission_persistent <- function(
    month_mat, cols, contaminated, t_offsets, month_mass,
    persistence = 0, heaped_month_mass = NULL, heaping = 0,
    year_revision_month_mass = NULL, year_revision = 0,
    clean_anchor_revision = 0, exact_anchor_retention = 0,
    local_revision = 0) {
  n <- nrow(month_mat)
  if (!length(contaminated)) return(list(loglik=rep(0, n),
    lambda_count=numeric(n), lambda_xsum=numeric(n)))
  if (!is.finite(persistence) || persistence < 0 || persistence >= 1)
    stop("tenure-report persistence must be in [0,1)")
  if (!is.finite(heaping) || heaping < 0 || heaping >= 1)
    stop("tenure-report heaping must be in [0,1)")
  if (!is.finite(year_revision) || year_revision < 0 || year_revision >= 1)
    stop("whole-year revision probability must be in [0,1)")
  if (!is.finite(clean_anchor_revision) || clean_anchor_revision < 0 ||
      clean_anchor_revision >= 1)
    stop("clean-anchor revision probability must be in [0,1)")
  if (!is.finite(exact_anchor_retention) || exact_anchor_retention < 0 ||
      exact_anchor_retention >= 1)
    stop("exact-anchor retention probability must be in [0,1)")
  if (!is.finite(local_revision) || local_revision < 0 ||
      local_revision >= 1)
    stop("local revision probability must be in [0,1)")
  if (persistence > 0 && (heaping > 0 || year_revision > 0 ||
      clean_anchor_revision > 0 || exact_anchor_retention > 0 ||
      local_revision > 0))
    stop("persistence is not jointly implemented with calendar revisions")
  if (heaping > 0 && (is.null(heaped_month_mass) ||
      length(heaped_month_mass) != length(month_mass)))
    stop("heaped month masses are required for positive heaping")

  gross_mass <- if (heaping == 0) month_mass else
    Map(function(a, b) .log_probability_mixture(a, b, heaping),
      month_mass, heaped_month_mass)

  # Separate gross and clean preceding anchors because whole-year revision has
  # type-specific probabilities. Exact retention applies to their union.
  candidate_gross <- contaminated[contaminated > 1L &
    (contaminated - 1L) %in% contaminated]
  candidate_clean <- contaminated[contaminated > 1L &
    !((contaminated - 1L) %in% contaminated)]
  candidate <- union(candidate_gross, candidate_clean)
  if ((year_revision > 0 || clean_anchor_revision > 0) &&
      (is.null(year_revision_month_mass) ||
      length(year_revision_month_mass) != length(month_mass)))
    stop("whole-year revision masses are required for positive revision probability")
  if (!length(candidate) || (persistence == 0 && year_revision == 0 &&
      clean_anchor_revision == 0 && exact_anchor_retention == 0 &&
      local_revision == 0)) {
    return(list(
      loglik=Reduce(`+`, gross_mass[contaminated]),
      lambda_count=rep(length(contaminated), n),
      lambda_xsum=rowSums((month_mat[, cols[contaminated], drop=FALSE] + .5)/12)))
  }

  # Conditional revision choices factor by current report once the observed
  # preceding report is fixed. Evaluate the product of two-component mixtures
  # directly instead of enumerating 2^J link patterns.
  if (year_revision > 0 || clean_anchor_revision > 0 ||
      exact_anchor_retention > 0 || local_revision > 0) {
    ll <- rep(0,n)
    for (pos in contaminated) {
      revision_probability <- if (pos %in% candidate_gross) year_revision else
        if (pos %in% candidate_clean) clean_anchor_revision else 0
      contribution <- if (revision_probability > 0)
        .log_probability_mixture(gross_mass[[pos]],
          year_revision_month_mass[[pos]],revision_probability) else
        gross_mass[[pos]]
      if (exact_anchor_retention > 0 && pos %in% candidate) {
        current_col <- cols[pos]
        previous_col <- cols[pos-1L]
        gap <- t_offsets[current_col]-t_offsets[previous_col]
        if (local_revision > 0) {
          local_mass <- .log_local_revision_month_mass(
            month_mat[,current_col],month_mat[,previous_col],gap)
          contribution <- .log_probability_mixture(contribution,local_mass,
            local_revision)
        }
        exact_mass <- ifelse(month_mat[,current_col] ==
          month_mat[,previous_col]+gap*.TENURE_MONTHS_PER_WAVE,0,-Inf)
        contribution <- .log_probability_mixture(contribution,exact_mass,
          exact_anchor_retention)
      } else if (local_revision > 0 && pos %in% candidate) {
        current_col <- cols[pos]
        previous_col <- cols[pos-1L]
        gap <- t_offsets[current_col]-t_offsets[previous_col]
        local_mass <- .log_local_revision_month_mass(
          month_mat[,current_col],month_mat[,previous_col],gap)
        contribution <- .log_probability_mixture(contribution,local_mass,
          local_revision)
      }
      ll <- ll+contribution
    }
    return(list(loglik=ll,
      lambda_count=rep(length(contaminated),n),
      lambda_xsum=rowSums((month_mat[,cols[contaminated],drop=FALSE]+.5)/12)))
  }

  choices <- as.matrix(expand.grid(rep(list(0:1), length(candidate))))
  lp <- matrix(-Inf, n, nrow(choices))
  lc <- lx <- matrix(0, n, nrow(choices))
  for (r in seq_len(nrow(choices))) {
    linked <- candidate[choices[r, ] == 1L]
    redraw <- setdiff(contaminated, linked)
    value <- if (length(redraw)) Reduce(`+`, gross_mass[redraw]) else
      rep(0, n)
    valid <- rep(TRUE, n)
    for (pos in candidate) {
      current_col <- cols[pos]
      previous_col <- cols[pos - 1L]
      gap <- t_offsets[current_col] - t_offsets[previous_col]
      retain <- if (year_revision > 0) year_revision else persistence^gap
      if (pos %in% linked) {
        value <- value + log(retain)
        if (year_revision > 0) {
          value <- value + year_revision_month_mass[[pos]]
        } else {
          valid <- valid & month_mat[, current_col] ==
            month_mat[, previous_col] + gap*.TENURE_MONTHS_PER_WAVE
        }
      } else value <- value + log1p(-retain)
    }
    value[!valid] <- -Inf
    lp[, r] <- value
    lc[, r] <- length(redraw)
    lx[, r] <- if (length(redraw))
      rowSums((month_mat[, cols[redraw], drop=FALSE] + .5)/12) else 0
  }
  ll <- .row_logsumexp(lp)
  post <- exp(lp - ll)
  post[!is.finite(post)] <- 0
  list(loglik=ll, lambda_count=rowSums(post*lc),
    lambda_xsum=rowSums(post*lx))
}

.job_segment_emission_monthly <- function(
    g_mat, s_mat, t_offsets, lambda_g, eps, beta_g = 0,
    start_model = c("marginal", "within_interval"),
    tenure_report_persistence = 0, tenure_heaping_prob = 0,
    tenure_year_revision_prob = 0,
    tenure_clean_anchor_revision_prob = 0,
    tenure_exact_anchor_retention_prob = 0,
    tenure_local_revision_prob = 0,
    tenure_start_month_probs = NULL,
    interview_month_mat = NULL) {
  start_model <- match.arg(start_model)
  if (!is.matrix(g_mat) || !is.matrix(s_mat) ||
      !identical(dim(g_mat), dim(s_mat)) || !is.logical(s_mat))
    stop("g_mat and logical s_mat must be conformable matrices")
  if (length(t_offsets) != ncol(g_mat))
    stop("t_offsets must have length ncol(g_mat)")
  if (!is.finite(eps) || eps <= 0 || eps >= 1)
    stop("eps must be in (0,1)")
  if (!is.finite(tenure_report_persistence) ||
      tenure_report_persistence < 0 || tenure_report_persistence >= 1)
    stop("tenure_report_persistence must be in [0,1)")
  if (!is.finite(tenure_heaping_prob) || tenure_heaping_prob < 0 ||
      tenure_heaping_prob >= 1)
    stop("tenure_heaping_prob must be in [0,1)")
  if (!is.finite(tenure_year_revision_prob) ||
      tenure_year_revision_prob < 0 || tenure_year_revision_prob >= 1)
    stop("tenure_year_revision_prob must be in [0,1)")
  if (!is.finite(tenure_clean_anchor_revision_prob) ||
      tenure_clean_anchor_revision_prob < 0 ||
      tenure_clean_anchor_revision_prob >= 1)
    stop("tenure_clean_anchor_revision_prob must be in [0,1)")
  if (!is.finite(tenure_exact_anchor_retention_prob) ||
      tenure_exact_anchor_retention_prob < 0 ||
      tenure_exact_anchor_retention_prob >= 1)
    stop("tenure_exact_anchor_retention_prob must be in [0,1)")
  if (!is.finite(tenure_local_revision_prob) ||
      tenure_local_revision_prob < 0 || tenure_local_revision_prob >= 1)
    stop("tenure_local_revision_prob must be in [0,1)")
  start_month_probs <- .validate_start_month_probs(tenure_start_month_probs)
  seasonal_start_month <- max(abs(start_month_probs-1/12))>1e-12
  if (tenure_report_persistence > 0 &&
      (tenure_heaping_prob > 0 || tenure_year_revision_prob > 0 ||
      tenure_clean_anchor_revision_prob > 0 ||
      tenure_exact_anchor_retention_prob > 0 ||
      tenure_local_revision_prob > 0))
    stop("persistence is not jointly implemented with calendar revisions")
  if ((tenure_heaping_prob > 0 || tenure_year_revision_prob > 0 ||
      tenure_clean_anchor_revision_prob > 0 ||
      tenure_exact_anchor_retention_prob > 0 ||
      tenure_local_revision_prob > 0 || seasonal_start_month) &&
      (is.null(interview_month_mat) ||
      !is.matrix(interview_month_mat) ||
      !identical(dim(interview_month_mat), dim(g_mat))))
    stop("calendar revisions require interview months conformable with g_mat")

  indexed <- .tenure_month_index(g_mat)
  month_mat <- matrix(indexed$month, nrow(g_mat), ncol(g_mat))
  valid_mat <- matrix(indexed$valid, nrow(g_mat), ncol(g_mat))
  if (any(s_mat & !valid_mat))
    stop("Observed tenure is not on the non-negative monthly grid")

  n <- nrow(g_mat)
  K <- as.integer(rowSums(s_mat))
  ans <- list(loglik=numeric(n), tau_sum=numeric(n), K=K,
    eps_informative=logical(n), lambda_count=numeric(n),
    lambda_xsum=numeric(n))
  if (!n || !ncol(g_mat)) return(ans)

  code <- as.vector(s_mat %*% (2^(seq_len(ncol(s_mat)) - 1L)))
  for (key in sort(unique(code[code > 0]))) {
    rows <- which(code == key)
    cols <- which(as.logical(intToBits(key)[seq_len(ncol(s_mat))]))
    k <- length(cols)
    types <- as.matrix(expand.grid(rep(list(0:1), k)))
    lp <- matrix(-Inf, length(rows), nrow(types))
    tau <- lc <- lx <- matrix(0, length(rows), nrow(types))
    fm <- lapply(cols,function(j) if (seasonal_start_month)
      .log_calendar_duration_month_mass(month_mat[rows,j],
        interview_month_mat[rows,j],lambda_g,beta_g,start_month_probs) else
      .log_duration_month_mass(month_mat[rows,j],lambda_g,beta_g))
    hm <- NULL
    log_norm <- NULL
    if (tenure_heaping_prob > 0 || tenure_year_revision_prob > 0 ||
        tenure_clean_anchor_revision_prob > 0)
      log_norm <- .log_january_duration_normalizers(lambda_g, beta_g)
    if (tenure_heaping_prob > 0) {
      hm <- lapply(cols, function(j)
        .log_january_duration_month_mass(month_mat[rows, j],
          interview_month_mat[rows, j], lambda_g, beta_g, log_norm))
    }
    ym <- vector("list", length(cols))
    if ((tenure_year_revision_prob > 0 ||
        tenure_clean_anchor_revision_prob > 0) && length(cols) > 1L) {
      for (pos in 2:length(cols)) {
        current_col <- cols[pos]
        previous_col <- cols[pos - 1L]
        gap <- t_offsets[current_col] - t_offsets[previous_col]
        ym[[pos]] <- .log_whole_year_revision_month_mass(
          month_mat[rows, current_col], month_mat[rows, previous_col],
          interview_month_mat[rows, current_col],
          interview_month_mat[rows, previous_col], rep(gap, length(rows)),
          lambda_g, beta_g, log_norm)
      }
    }

    for (r in seq_len(nrow(types))) {
      contaminated <- which(types[r, ] == 1L)
      clean <- which(types[r, ] == 0L)
      n_x <- length(contaminated)
      log_mix <- n_x * log(eps) + (k - n_x) * log1p(-eps)
      gross <- .gross_month_emission_persistent(
        month_mat[rows, , drop=FALSE], cols, contaminated, t_offsets, fm,
        tenure_report_persistence, hm, tenure_heaping_prob, ym,
        tenure_year_revision_prob,tenure_clean_anchor_revision_prob,
        tenure_exact_anchor_retention_prob,tenure_local_revision_prob)

      if (!length(clean)) {
        lp[, r] <- log_mix + gross$loglik
        tau[, r] <- n_x
        lc[, r] <- gross$lambda_count
        lx[, r] <- gross$lambda_xsum
        next
      }

      anchor <- clean[1L]
      M0 <- month_mat[rows, cols[anchor]] -
        t_offsets[cols[anchor]] * .TENURE_MONTHS_PER_WAVE
      valid <- M0 >= 0
      if (length(clean) > 1L) for (j in clean[-1L])
        valid <- valid & month_mat[rows, cols[j]] ==
          M0 + t_offsets[cols[j]] * .TENURE_MONTHS_PER_WAVE

      if (identical(start_model, "within_interval")) {
        valid <- valid & M0 %in% .RESET_TENURE_MONTHS
        base <- if (seasonal_start_month)
          .log_within_interval_start_month_mass(M0,
            interview_month_mat[rows,cols[anchor]],start_month_probs) else
          rep(-log(length(.RESET_TENURE_MONTHS)),length(rows))
        base_count <- 0
        base_xsum <- 0
      } else {
        reference_interview <- if (seasonal_start_month)
          ((as.integer(interview_month_mat[rows,cols[anchor]])-1L-
            t_offsets[cols[anchor]]*.TENURE_MONTHS_PER_WAVE)%%12L)+1L else
          NULL
        base <- if (seasonal_start_month)
          .log_calendar_duration_month_mass(M0,reference_interview,
            lambda_g,beta_g,start_month_probs) else
          .log_duration_month_mass(M0,lambda_g,beta_g)
        base_count <- 1
        base_xsum <- (M0 + .5) / 12
      }
      value <- log_mix + gross$loglik + base
      value[!valid] <- -Inf
      lp[, r] <- value
      tau[, r] <- n_x
      lc[, r] <- gross$lambda_count + base_count
      lx[, r] <- gross$lambda_xsum + base_xsum
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

log_emission_spell_g_monthly <- function(
    g_mat, s_mat, t_offsets, lambda_g, eps, job_change_prob = 0,
    beta_g = 0, initial_model = c("marginal", "within_interval"),
    tenure_report_persistence = 0, tenure_heaping_prob = 0,
    tenure_year_revision_prob = 0,
    tenure_clean_anchor_revision_prob = 0,
    tenure_exact_anchor_retention_prob = 0,
    tenure_local_revision_prob = 0,
    tenure_start_month_probs = NULL,
    interview_month_mat = NULL, segment_cache = NULL,
    segment_ids = seq_len(ncol(g_mat)), compress_segments = TRUE) {
  initial_model <- match.arg(initial_model)
  if (is.null(segment_cache)) segment_cache <- new.env(parent = emptyenv())
  segment_value <- function(cols, start_model) {
    key <- paste(c(start_model, segment_ids[cols]), collapse = "_")
    if (exists(key, segment_cache, inherits = FALSE))
      return(get(key, segment_cache, inherits = FALSE))
    gg <- g_mat[, cols, drop = FALSE]
    ss <- s_mat[, cols, drop = FALSE]
    ii <- if (is.null(interview_month_mat)) NULL else
      interview_month_mat[, cols, drop = FALSE]
    rows <- seq_len(nrow(gg))
    index <- rows
    if (compress_segments && length(rows)) {
      cell_key <- do.call(paste, c(as.data.frame(cbind(gg, ss, ii)),
        sep = "\r"))
      first <- !duplicated(cell_key)
      rows <- which(first)
      index <- match(cell_key, cell_key[first])
    }
    out <- .job_segment_emission_monthly(
      gg[rows, , drop = FALSE], ss[rows, , drop = FALSE],
      as.integer(seq_along(cols) - 1L), lambda_g, eps, beta_g,
      start_model = start_model,
      tenure_report_persistence = tenure_report_persistence,
      tenure_heaping_prob = tenure_heaping_prob,
      tenure_year_revision_prob = tenure_year_revision_prob,
      tenure_clean_anchor_revision_prob = tenure_clean_anchor_revision_prob,
      tenure_exact_anchor_retention_prob = tenure_exact_anchor_retention_prob,
      tenure_local_revision_prob = tenure_local_revision_prob,
      tenure_start_month_probs = tenure_start_month_probs,
      interview_month_mat = if (is.null(ii)) NULL else ii[rows, , drop = FALSE])
    out <- lapply(out, function(value) value[index])
    assign(key, out, segment_cache)
    out
  }
  if (!is.finite(job_change_prob) || job_change_prob < 0 ||
      job_change_prob >= 1)
    stop("job_change_prob must be in [0,1)")
  L <- ncol(g_mat)
  if (L <= 1L || job_change_prob == 0) {
    out <- .job_segment_emission_monthly(g_mat, s_mat, t_offsets,
      lambda_g, eps, beta_g, start_model=initial_model,
      tenure_report_persistence=tenure_report_persistence,
      tenure_heaping_prob=tenure_heaping_prob,
      tenure_year_revision_prob=tenure_year_revision_prob,
      tenure_clean_anchor_revision_prob=tenure_clean_anchor_revision_prob,
      tenure_exact_anchor_retention_prob=
        tenure_exact_anchor_retention_prob,
      tenure_local_revision_prob=tenure_local_revision_prob,
      tenure_start_month_probs=tenure_start_month_probs,
      interview_month_mat=interview_month_mat)
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
      segment <- segment_value(cols,
        if (b == 1L) initial_model else "within_interval")
      lp[, r] <- lp[, r] + segment$loglik
      tau[, r] <- tau[, r] + segment$tau_sum
      lc[, r] <- lc[, r] + segment$lambda_count
      lx[, r] <- lx[, r] + segment$lambda_xsum
    }
  }
  ll <- .row_logsumexp(lp)
  post <- exp(lp - ll)
  list(loglik=ll, tau_sum=rowSums(post*tau),
    K=as.integer(rowSums(s_mat)), eps_informative=rowSums(s_mat)>=2L,
    lambda_count=rowSums(post*lc), lambda_xsum=rowSums(post*lx),
    job_changes=as.vector(post %*% changes),
    job_change_opportunities=rep(L-1L, n))
}

.piecewise_persistent_monthly_pack <- function(
    params, persistence_start=.5, q_start=NULL) {
  rho <- min(max(persistence_start, 1e-8), 1-1e-8)
  c(.piecewise_job_change_pack(params, q_start=q_start),
    tenure_persistence=qlogis(rho))
}

.piecewise_persistent_monthly_unpack <- function(
    z, timegap_contamination_model="joint_marginal") {
  base <- z[setdiff(names(z), "tenure_persistence")]
  p <- .piecewise_job_change_unpack(base,
    timegap_contamination_model=timegap_contamination_model)
  p$tenure_measurement_model <- "monthly"
  p$tenure_report_persistence <-
    plogis(unname(z["tenure_persistence"]))
  p
}

# Jointly estimate the corrected monthly model and the probability that a
# consecutive gross tenure report preserves the preceding reported start-date
# anchor. This is a direct observed-likelihood fit with parallel central
# differences, matching the numerical approach used for the job-reset model.
fit_eps_piecewise_persistent_monthly <- function(
    df, start, persistence_start=.5, q_start=NULL, maxit=400L,
    reltol=1e-9, pgtol=1e-7, method="L-BFGS-B", verbose=0L,
    timegap_contamination_model="joint_marginal", workers=1L,
    gradient_step=1e-4) {
  validate_df_eps(df, allow_zero_tenure=TRUE)
  z0 <- .piecewise_persistent_monthly_pack(start, persistence_start, q_start)
  unpack <- function(z) .piecewise_persistent_monthly_unpack(z,
    timegap_contamination_model)
  lower <- rep(-Inf, length(z0)); upper <- rep(Inf, length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- c(paste0("log_hg",1:5), paste0("log_hd",1:5))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  lower["job_change"] <- -16; upper["job_change"] <- 8
  lower["tenure_persistence"] <- -12
  upper["tenure_persistence"] <- 12
  total_weight <- sum(df$weight)
  objective <- function(z) {
    value <- tryCatch(e_step_eps(df, unpack(z), check_df=FALSE,
      suff_stats=FALSE)$loglik, error=function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value/total_weight
  }

  cluster <- NULL; gradient <- NULL
  workers <- as.integer(workers)
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add=TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster, function(path) {
      setwd(path); source("EM-tenure/R/source_all.R"); NULL
    }, worker_wd)
    df_worker <- df
    total_weight_worker <- total_weight
    timegap_worker <- timegap_contamination_model
    parallel::clusterExport(cluster,
      c("df_worker", "total_weight_worker", "timegap_worker"),
      envir=environment())
    worker_objective <- function(z) {
      p <- .piecewise_persistent_monthly_unpack(z, timegap_worker)
      value <- tryCatch(e_step_eps(df_worker, p, check_df=FALSE,
        suff_stats=FALSE)$loglik, error=function(e) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value/total_weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step*pmax(1,abs(z))
      points <- vector("list",2L*length(z))
      for (j in seq_along(z)) {
        points[[2L*j-1L]] <- points[[2L*j]] <- z
        points[[2L*j-1L]][j] <- z[j]+step[j]
        points[[2L*j]][j] <- z[j]-step[j]
      }
      values <- unlist(parallel::parLapply(cluster,points,
        worker_objective),use.names=FALSE)
      plus <- values[seq(1L,length(values),by=2L)]
      minus <- values[seq(2L,length(values),by=2L)]
      centre <- objective(z)
      valid_plus <- is.finite(plus) & plus<1e50
      valid_minus <- is.finite(minus) & minus<1e50
      out <- numeric(length(z))
      both <- valid_plus & valid_minus
      out[both] <- (plus[both]-minus[both])/(2*step[both])
      plus_only <- valid_plus & !valid_minus
      out[plus_only] <- (plus[plus_only]-centre)/step[plus_only]
      minus_only <- !valid_plus & valid_minus
      out[minus_only] <- (centre-minus[minus_only])/step[minus_only]
      out[!is.finite(out)] <- 0
      out
    }
  }
  control <- list(maxit=maxit, reltol=reltol,
    ndeps=rep(1e-3,length(z0)), trace=as.integer(verbose))
  if (identical(method,"L-BFGS-B")) {
    control$factr <- reltol/.Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0,objective,gr=gradient,method=method,
      lower=lower,upper=upper,control=control)
  } else opt <- optim(z0,objective,gr=gradient,method=method,control=control)
  params <- unpack(opt$par)
  estep <- e_step_eps(df,params,check_df=FALSE,suff_stats=FALSE)
  list(params=params,loglik=estep$loglik,gamma=estep$gamma,
    job_change_posterior=estep$job_change_posterior,
    convergence=opt$convergence,message=opt$message,
    iterations=unname(opt$counts["function"]),objective=opt$value,
    par_unconstrained=opt$par,objective_function=objective)
}

.piecewise_heaping_monthly_pack <- function(
    params, heaping_start=.2, q_start=NULL) {
  heap <- min(max(heaping_start, 1e-8), 1-1e-8)
  c(.piecewise_job_change_pack(params, q_start=q_start),
    tenure_heaping=qlogis(heap))
}

.piecewise_heaping_monthly_unpack <- function(
    z, timegap_contamination_model="joint_marginal") {
  base <- z[setdiff(names(z), "tenure_heaping")]
  p <- .piecewise_job_change_unpack(base,
    timegap_contamination_model=timegap_contamination_model)
  p$tenure_measurement_model <- "monthly"
  p$tenure_report_persistence <- 0
  p$tenure_heaping_prob <- plogis(unname(z["tenure_heaping"]))
  p
}

# Jointly estimate the corrected monthly model and the conditional probability
# that a gross tenure report is drawn from the marginal duration distribution
# restricted to an implied January job-start month.  The independent-gross
# model is nested exactly at h=0.
fit_eps_piecewise_heaping_monthly <- function(
    df, start, heaping_start=.2, q_start=NULL, maxit=400L,
    reltol=1e-9, pgtol=1e-7, method="L-BFGS-B", verbose=0L,
    timegap_contamination_model="joint_marginal", workers=1L,
    gradient_step=1e-4, heaping_fixed=NULL) {
  validate_df_eps(df, allow_zero_tenure=TRUE)
  interview_cols <- paste0("interview_month",1:3)
  if (length(setdiff(interview_cols,names(df))) ||
      any(!as.matrix(df[interview_cols]) %in% 1:12))
    stop("heaping fit requires interview_month1-3 in 1:12")
  if (!is.null(heaping_fixed) && (!is.finite(heaping_fixed) ||
      heaping_fixed<0 || heaping_fixed>=1))
    stop("heaping_fixed must be in [0,1)")
  z0 <- if (is.null(heaping_fixed))
    .piecewise_heaping_monthly_pack(start,heaping_start,q_start) else
    .piecewise_job_change_pack(start,q_start=q_start)
  unpack <- function(z) {
    if (is.null(heaping_fixed)) return(.piecewise_heaping_monthly_unpack(z,
      timegap_contamination_model))
    p <- .piecewise_job_change_unpack(z,
      timegap_contamination_model=timegap_contamination_model)
    p$tenure_measurement_model <- "monthly"
    p$tenure_report_persistence <- 0
    p$tenure_heaping_prob <- heaping_fixed
    p
  }
  lower <- rep(-Inf,length(z0)); upper <- rep(Inf,length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- c(paste0("log_hg",1:5),paste0("log_hd",1:5))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  lower["job_change"] <- -16; upper["job_change"] <- 8
  if ("tenure_heaping" %in% names(z0)) {
    lower["tenure_heaping"] <- -12; upper["tenure_heaping"] <- 12
  }
  total_weight <- sum(df$weight)
  objective <- function(z) {
    value <- tryCatch(e_step_eps(df,unpack(z),check_df=FALSE,
      suff_stats=FALSE)$loglik,error=function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value/total_weight
  }

  cluster <- NULL; gradient <- NULL
  workers <- as.integer(workers)
  if (workers>1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster,function(path) {
      setwd(path); source("EM-tenure/R/source_all.R"); NULL
    },worker_wd)
    df_worker <- df
    total_weight_worker <- total_weight
    timegap_worker <- timegap_contamination_model
    heaping_fixed_worker <- heaping_fixed
    parallel::clusterExport(cluster,
      c("df_worker","total_weight_worker","timegap_worker",
        "heaping_fixed_worker"),
      envir=environment())
    worker_objective <- function(z) {
      if (is.null(heaping_fixed_worker)) {
        p <- .piecewise_heaping_monthly_unpack(z,timegap_worker)
      } else {
        p <- .piecewise_job_change_unpack(z,
          timegap_contamination_model=timegap_worker)
        p$tenure_measurement_model <- "monthly"
        p$tenure_report_persistence <- 0
        p$tenure_heaping_prob <- heaping_fixed_worker
      }
      value <- tryCatch(e_step_eps(df_worker,p,check_df=FALSE,
        suff_stats=FALSE)$loglik,error=function(e) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value/total_weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step*pmax(1,abs(z))
      points <- vector("list",2L*length(z))
      for (j in seq_along(z)) {
        points[[2L*j-1L]] <- points[[2L*j]] <- z
        points[[2L*j-1L]][j] <- z[j]+step[j]
        points[[2L*j]][j] <- z[j]-step[j]
      }
      values <- unlist(parallel::parLapply(cluster,points,
        worker_objective),use.names=FALSE)
      plus <- values[seq(1L,length(values),by=2L)]
      minus <- values[seq(2L,length(values),by=2L)]
      centre <- objective(z)
      valid_plus <- is.finite(plus) & plus<1e50
      valid_minus <- is.finite(minus) & minus<1e50
      out <- numeric(length(z))
      both <- valid_plus & valid_minus
      out[both] <- (plus[both]-minus[both])/(2*step[both])
      plus_only <- valid_plus & !valid_minus
      out[plus_only] <- (plus[plus_only]-centre)/step[plus_only]
      minus_only <- !valid_plus & valid_minus
      out[minus_only] <- (centre-minus[minus_only])/step[minus_only]
      out[!is.finite(out)] <- 0
      out
    }
  }
  control <- list(maxit=maxit,reltol=reltol,
    ndeps=rep(1e-3,length(z0)),trace=as.integer(verbose))
  if (identical(method,"L-BFGS-B")) {
    control$factr <- reltol/.Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0,objective,gr=gradient,method=method,
      lower=lower,upper=upper,control=control)
  } else opt <- optim(z0,objective,gr=gradient,method=method,control=control)
  params <- unpack(opt$par)
  estep <- e_step_eps(df,params,check_df=FALSE,suff_stats=FALSE)
  list(params=params,loglik=estep$loglik,gamma=estep$gamma,
    job_change_posterior=estep$job_change_posterior,
    convergence=opt$convergence,message=opt$message,
    iterations=unname(opt$counts["function"]),objective=opt$value,
    par_unconstrained=opt$par,objective_function=objective)
}

.piecewise_calendar_revision_monthly_pack <- function(
    params, heaping_start=.2, revision_start=.1, q_start=NULL,
    clean_anchor_revision_start=NULL,start_month_probs_start=NULL,
    exact_anchor_retention_start=NULL,local_revision_start=NULL,
    reliability_shift_start=NULL,tenure_reliability_shift_start=NULL,
    timegap_reliability_shift_start=NULL) {
  heap <- min(max(heaping_start,1e-8),1-1e-8)
  revision <- min(max(revision_start,1e-8),1-1e-8)
  out <- c(.piecewise_job_change_pack(params,q_start=q_start),
    tenure_heaping=qlogis(heap),tenure_year_revision=qlogis(revision))
  if (!is.null(clean_anchor_revision_start)) {
    clean_revision <- min(max(clean_anchor_revision_start,1e-8),1-1e-8)
    out <- c(out,tenure_clean_anchor_revision=qlogis(clean_revision))
  }
  if (!is.null(start_month_probs_start)) {
    probs <- .validate_start_month_probs(start_month_probs_start)
    out <- c(out,setNames(log(probs[1:11]/probs[12]),
      paste0("start_month_logit",1:11)))
  }
  if (!is.null(exact_anchor_retention_start)) {
    retention <- min(max(exact_anchor_retention_start,1e-8),1-1e-8)
    out <- c(out,tenure_exact_anchor_retention=qlogis(retention))
  }
  if (!is.null(local_revision_start)) {
    local <- min(max(local_revision_start,1e-8),1-1e-8)
    out <- c(out,tenure_local_revision=qlogis(local))
  }
  if (!is.null(reliability_shift_start)) {
    if (!is.finite(reliability_shift_start) || reliability_shift_start<0)
      stop("reliability_shift_start must be finite and nonnegative")
    out <- c(out,duration_reliability_shift=reliability_shift_start)
  }
  separate_starts <- c(tenure=tenure_reliability_shift_start,
    timegap=timegap_reliability_shift_start)
  if (length(separate_starts)) {
    if (!is.null(reliability_shift_start))
      stop("common and separate reliability shifts cannot both be packed")
    if (length(separate_starts)!=2L || any(!is.finite(separate_starts)) ||
        any(separate_starts<0))
      stop("both separate reliability shifts must be finite and nonnegative")
    out <- c(out,tenure_reliability_shift=unname(separate_starts["tenure"]),
      timegap_reliability_shift=unname(separate_starts["timegap"]))
  }
  out
}

.piecewise_calendar_revision_monthly_unpack <- function(
    z, timegap_contamination_model="joint_marginal") {
  calendar_names <- c("tenure_heaping","tenure_year_revision",
    "tenure_clean_anchor_revision","tenure_exact_anchor_retention",
    "tenure_local_revision","duration_reliability_shift",
    "tenure_reliability_shift","timegap_reliability_shift",
    paste0("start_month_logit",1:11))
  base <- z[setdiff(names(z),calendar_names)]
  p <- .piecewise_job_change_unpack(base,
    timegap_contamination_model=timegap_contamination_model)
  p$tenure_measurement_model <- "monthly"
  p$tenure_report_persistence <- 0
  p$tenure_heaping_prob <- plogis(unname(z["tenure_heaping"]))
  p$tenure_year_revision_prob <-
    plogis(unname(z["tenure_year_revision"]))
  p$tenure_clean_anchor_revision_prob <-
    if ("tenure_clean_anchor_revision" %in% names(z))
      plogis(unname(z["tenure_clean_anchor_revision"])) else 0
  p$tenure_exact_anchor_retention_prob <-
    if ("tenure_exact_anchor_retention" %in% names(z))
      plogis(unname(z["tenure_exact_anchor_retention"])) else 0
  p$tenure_local_revision_prob <-
    if ("tenure_local_revision" %in% names(z))
      plogis(unname(z["tenure_local_revision"])) else 0
  p$duration_reliability_shift <-
    if ("duration_reliability_shift" %in% names(z))
      unname(z["duration_reliability_shift"]) else 0
  p$tenure_reliability_shift <-
    if ("tenure_reliability_shift" %in% names(z))
      unname(z["tenure_reliability_shift"]) else NULL
  p$timegap_reliability_shift <-
    if ("timegap_reliability_shift" %in% names(z))
      unname(z["timegap_reliability_shift"]) else NULL
  month_names <- paste0("start_month_logit",1:11)
  if (all(month_names %in% names(z))) {
    logits <- c(unname(z[month_names]),0)
    logits <- logits-max(logits)
    p$tenure_start_month_probs <- exp(logits)/sum(exp(logits))
  } else p$tenure_start_month_probs <- rep(1/12,12L)
  p
}

# Extend January heaping with a correlated date-revision branch. Conditional
# on two consecutive gross reports for the same latent job, omega is the
# probability that the new report preserves the previous reported start month
# but selects a different year. The January-heaping model is nested at omega=0.
fit_eps_piecewise_calendar_revision_monthly <- function(
    df, start, heaping_start=.2, revision_start=.1, q_start=NULL,
    maxit=400L, reltol=1e-9, pgtol=1e-7, method="L-BFGS-B", verbose=0L,
    timegap_contamination_model="joint_marginal", workers=1L,
    gradient_step=1e-4, revision_fixed=NULL,
    clean_anchor_revision_start=NULL,start_month_probs_start=NULL,
    exact_anchor_retention_start=NULL,
    local_revision_start=NULL,
    reliability_shift_start=NULL,
    tenure_reliability_shift_start=NULL,
    timegap_reliability_shift_start=NULL,
    free_names=NULL,
    gradient_scheme=c("central","forward")) {
  gradient_scheme <- match.arg(gradient_scheme)
  validate_df_eps(df,allow_zero_tenure=TRUE)
  interview_cols <- paste0("interview_month",1:3)
  if (length(setdiff(interview_cols,names(df))) ||
      any(!as.matrix(df[interview_cols]) %in% 1:12))
    stop("calendar-revision fit requires interview_month1-3 in 1:12")
  if (!is.null(revision_fixed) && (!is.finite(revision_fixed) ||
      revision_fixed<0 || revision_fixed>=1))
    stop("revision_fixed must be in [0,1)")
  if (!is.null(revision_fixed) && !is.null(clean_anchor_revision_start))
    stop("fixed gross-anchor revisions and free clean-anchor revisions are not jointly implemented")
  if (!is.null(revision_fixed) && !is.null(exact_anchor_retention_start))
    stop("fixed gross-anchor revisions and free exact retention are not jointly implemented")
  if (!is.null(revision_fixed) && !is.null(local_revision_start))
    stop("fixed gross-anchor revisions and free local revisions are not jointly implemented")
  if (!is.null(revision_fixed) && (!is.null(reliability_shift_start) ||
      !is.null(tenure_reliability_shift_start) ||
      !is.null(timegap_reliability_shift_start)))
    stop("fixed gross-anchor revisions and reliability mixtures are not jointly implemented")
  z_full0 <- if (is.null(revision_fixed))
    .piecewise_calendar_revision_monthly_pack(start,heaping_start,
      revision_start,q_start,clean_anchor_revision_start,
      start_month_probs_start,exact_anchor_retention_start,
      local_revision_start,reliability_shift_start,
      tenure_reliability_shift_start,timegap_reliability_shift_start) else
    .piecewise_heaping_monthly_pack(start,heaping_start,q_start)
  unpack_full <- function(z) {
    if (is.null(revision_fixed))
      return(.piecewise_calendar_revision_monthly_unpack(z,
        timegap_contamination_model))
    p <- .piecewise_heaping_monthly_unpack(z,timegap_contamination_model)
    p$tenure_year_revision_prob <- revision_fixed
    p
  }
  lower_full <- rep(-Inf,length(z_full0)); upper_full <- rep(Inf,length(z_full0))
  names(lower_full) <- names(upper_full) <- names(z_full0)
  hz <- c(paste0("log_hg",1:5),paste0("log_hd",1:5))
  lower_full[hz] <- log(1e-4); upper_full[hz] <- log(20)
  lower_full["job_change"] <- -16; upper_full["job_change"] <- 8
  lower_full["tenure_heaping"] <- -12
  upper_full["tenure_heaping"] <- 12
  if ("tenure_year_revision" %in% names(z_full0)) {
    lower_full["tenure_year_revision"] <- -12
    upper_full["tenure_year_revision"] <- 12
  }
  if ("tenure_clean_anchor_revision" %in% names(z_full0)) {
    lower_full["tenure_clean_anchor_revision"] <- -12
    upper_full["tenure_clean_anchor_revision"] <- 12
  }
  if ("tenure_exact_anchor_retention" %in% names(z_full0)) {
    lower_full["tenure_exact_anchor_retention"] <- -12
    upper_full["tenure_exact_anchor_retention"] <- 12
  }
  if ("tenure_local_revision" %in% names(z_full0)) {
    lower_full["tenure_local_revision"] <- -12
    upper_full["tenure_local_revision"] <- 12
  }
  if ("duration_reliability_shift" %in% names(z_full0)) {
    lower_full["duration_reliability_shift"] <- 0
    upper_full["duration_reliability_shift"] <- 6
  }
  for (shift_name in intersect(c("tenure_reliability_shift",
      "timegap_reliability_shift"),names(z_full0))) {
    lower_full[shift_name] <- 0
    upper_full[shift_name] <- 6
  }
  month_names <- intersect(paste0("start_month_logit",1:11),names(z_full0))
  lower_full[month_names] <- -8
  upper_full[month_names] <- 8
  if (is.null(free_names)) free_names <- names(z_full0)
  if (!length(free_names) || anyDuplicated(free_names) ||
      length(setdiff(free_names,names(z_full0))))
    stop("free_names must be unique names in the packed parameter vector")
  z0 <- z_full0[free_names]
  lower <- lower_full[free_names]
  upper <- upper_full[free_names]
  expand_z <- function(z) {
    full <- z_full0
    full[names(z)] <- z
    full
  }
  unpack <- function(z) unpack_full(expand_z(z))
  total_weight <- sum(df$weight)
  objective <- function(z) {
    value <- tryCatch(e_step_eps(df,unpack(z),check_df=FALSE,
      suff_stats=FALSE)$loglik,error=function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value/total_weight
  }

  cluster <- NULL; gradient <- NULL
  workers <- as.integer(workers)
  if (workers>1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster,function(path) {
      setwd(path); source("EM-tenure/R/source_all.R"); NULL
    },worker_wd)
    df_worker <- df
    total_weight_worker <- total_weight
    timegap_worker <- timegap_contamination_model
    revision_fixed_worker <- revision_fixed
    z_full0_worker <- z_full0
    free_names_worker <- free_names
    parallel::clusterExport(cluster,
      c("df_worker","total_weight_worker","timegap_worker",
        "revision_fixed_worker","z_full0_worker","free_names_worker"),
      envir=environment())
    worker_objective <- function(z) {
      full_z <- z_full0_worker
      full_z[free_names_worker] <- z
      if (is.null(revision_fixed_worker)) {
        p <- .piecewise_calendar_revision_monthly_unpack(full_z,timegap_worker)
      } else {
        p <- .piecewise_heaping_monthly_unpack(full_z,timegap_worker)
        p$tenure_year_revision_prob <- revision_fixed_worker
      }
      value <- tryCatch(e_step_eps(df_worker,p,check_df=FALSE,
        suff_stats=FALSE)$loglik,error=function(e) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value/total_weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step*pmax(1,abs(z))
      if (identical(gradient_scheme,"forward")) {
        points <- lapply(seq_along(z),function(j) {
          point <- z; point[j] <- point[j]+step[j]; point
        })
        plus <- unlist(parallel::parLapplyLB(cluster,points,
          worker_objective),use.names=FALSE)
        centre <- objective(z)
        valid_plus <- is.finite(plus) & plus<1e50
        out <- numeric(length(z))
        out[valid_plus] <- (plus[valid_plus]-centre)/step[valid_plus]
        invalid <- which(!valid_plus)
        if (length(invalid)) {
          minus_points <- lapply(invalid,function(j) {
            point <- z; point[j] <- point[j]-step[j]; point
          })
          minus <- unlist(parallel::parLapplyLB(cluster,minus_points,
            worker_objective),use.names=FALSE)
          valid_minus <- is.finite(minus) & minus<1e50
          out[invalid[valid_minus]] <-
            (centre-minus[valid_minus])/step[invalid[valid_minus]]
        }
        out[!is.finite(out)] <- 0
        out <- pmin(pmax(out,-1e6),1e6)
        return(out)
      }
      points <- vector("list",2L*length(z))
      for (j in seq_along(z)) {
        points[[2L*j-1L]] <- points[[2L*j]] <- z
        points[[2L*j-1L]][j] <- z[j]+step[j]
        points[[2L*j]][j] <- z[j]-step[j]
      }
      values <- unlist(parallel::parLapply(cluster,points,
        worker_objective),use.names=FALSE)
      plus <- values[seq(1L,length(values),by=2L)]
      minus <- values[seq(2L,length(values),by=2L)]
      centre <- objective(z)
      valid_plus <- is.finite(plus) & plus<1e50
      valid_minus <- is.finite(minus) & minus<1e50
      out <- numeric(length(z))
      both <- valid_plus & valid_minus
      out[both] <- (plus[both]-minus[both])/(2*step[both])
      plus_only <- valid_plus & !valid_minus
      out[plus_only] <- (plus[plus_only]-centre)/step[plus_only]
      minus_only <- !valid_plus & valid_minus
      out[minus_only] <- (centre-minus[minus_only])/step[minus_only]
      out[!is.finite(out)] <- 0
      out
    }
  }
  control <- list(maxit=maxit,reltol=reltol,
    ndeps=rep(1e-3,length(z0)),trace=as.integer(verbose))
  if (identical(method,"L-BFGS-B")) {
    control$factr <- reltol/.Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0,objective,gr=gradient,method=method,
      lower=lower,upper=upper,control=control)
  } else opt <- optim(z0,objective,gr=gradient,method=method,control=control)
  full_opt <- expand_z(opt$par)
  params <- unpack_full(full_opt)
  estep <- e_step_eps(df,params,check_df=FALSE,suff_stats=FALSE)
  list(params=params,loglik=estep$loglik,gamma=estep$gamma,
    job_change_posterior=estep$job_change_posterior,
    convergence=opt$convergence,message=opt$message,
    iterations=unname(opt$counts["function"]),objective=opt$value,
    par_unconstrained=full_opt,objective_function=objective,
    free_names=free_names)
}

# Reparameterize the two nonnegative duration-reliability shifts by their
# average and a fixed timegap-minus-tenure difference. This is used for a
# nuisance-adjusted likelihood profile of the difference without changing the
# underlying likelihood or the interpretation of either shift.
.piecewise_calendar_revision_monthly_difference_reduce <- function(z_full) {
  shift_names <- c("tenure_reliability_shift",
    "timegap_reliability_shift")
  if (length(setdiff(shift_names,names(z_full))))
    stop("both separate reliability shifts are required")
  level <- mean(unname(z_full[shift_names]))
  reduced <- z_full[setdiff(names(z_full),shift_names)]
  c(reduced,reliability_shift_level=level)
}

.piecewise_calendar_revision_monthly_difference_expand <- function(
    z_reduced, difference) {
  if (!is.finite(difference) || abs(difference)>6)
    stop("difference must be finite and lie in [-6,6]")
  if (!"reliability_shift_level" %in% names(z_reduced))
    stop("reliability_shift_level is required")
  level <- unname(z_reduced["reliability_shift_level"])
  z_full <- z_reduced[names(z_reduced)!="reliability_shift_level"]
  c(z_full,tenure_reliability_shift=level-difference/2,
    timegap_reliability_shift=level+difference/2)
}

fit_eps_piecewise_calendar_revision_monthly_fixed_difference <- function(
    df, start, difference, heaping_start=.2, revision_start=.1,
    q_start=NULL, maxit=40L, reltol=1e-9, pgtol=1e-7,
    method="L-BFGS-B", verbose=0L,
    timegap_contamination_model="joint_marginal", workers=1L,
    gradient_step=1e-4, clean_anchor_revision_start=NULL,
    start_month_probs_start=NULL, exact_anchor_retention_start=NULL,
    local_revision_start=NULL, gradient_scheme=c("central","forward")) {
  gradient_scheme <- match.arg(gradient_scheme)
  validate_df_eps(df,allow_zero_tenure=TRUE)
  interview_cols <- paste0("interview_month",1:3)
  if (length(setdiff(interview_cols,names(df))) ||
      any(!as.matrix(df[interview_cols]) %in% 1:12))
    stop("fixed-difference fit requires interview_month1-3 in 1:12")
  if (!is.finite(difference) || abs(difference)>6)
    stop("difference must be finite and lie in [-6,6]")

  tenure_start <- start$tenure_reliability_shift
  timegap_start <- start$timegap_reliability_shift
  if (is.null(tenure_start) || is.null(timegap_start)) {
    common <- start$duration_reliability_shift
    if (is.null(common)) common <- 0
    tenure_start <- timegap_start <- common
  }
  z_full0 <- .piecewise_calendar_revision_monthly_pack(start,
    heaping_start,revision_start,q_start,clean_anchor_revision_start,
    start_month_probs_start,exact_anchor_retention_start,
    local_revision_start,NULL,tenure_start,timegap_start)
  z0 <- .piecewise_calendar_revision_monthly_difference_reduce(z_full0)
  level_lower <- abs(difference)/2
  level_upper <- 6-abs(difference)/2
  z0["reliability_shift_level"] <- min(level_upper,max(level_lower,
    mean(c(tenure_start,timegap_start))))

  lower <- rep(-Inf,length(z0)); upper <- rep(Inf,length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- intersect(c(paste0("log_hg",1:5),paste0("log_hd",1:5)),names(z0))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  if ("job_change" %in% names(z0)) {
    lower["job_change"] <- -16; upper["job_change"] <- 8
  }
  probability_logits <- intersect(c("tenure_heaping",
    "tenure_year_revision","tenure_clean_anchor_revision",
    "tenure_exact_anchor_retention","tenure_local_revision"),names(z0))
  lower[probability_logits] <- -12; upper[probability_logits] <- 12
  month_names <- intersect(paste0("start_month_logit",1:11),names(z0))
  lower[month_names] <- -8; upper[month_names] <- 8
  lower["reliability_shift_level"] <- level_lower
  upper["reliability_shift_level"] <- level_upper

  expand_z <- function(z) {
    full <- .piecewise_calendar_revision_monthly_difference_expand(z,
      difference)
    # Restore the original packed order so saved parameter vectors remain
    # directly comparable with the unrestricted estimates.
    full[names(z_full0)]
  }
  unpack <- function(z) .piecewise_calendar_revision_monthly_unpack(
    expand_z(z),timegap_contamination_model)
  total_weight <- sum(df$weight)
  objective <- function(z) {
    value <- tryCatch(e_step_eps(df,unpack(z),check_df=FALSE,
      suff_stats=FALSE)$loglik,error=function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value/total_weight
  }

  cluster <- NULL; gradient <- NULL
  workers <- as.integer(workers)
  if (workers>1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
    worker_wd <- getwd()
    parallel::clusterCall(cluster,function(path) {
      setwd(path); source("EM-tenure/R/source_all.R"); NULL
    },worker_wd)
    df_worker <- df
    total_weight_worker <- total_weight
    timegap_worker <- timegap_contamination_model
    difference_worker <- difference
    full_names_worker <- names(z_full0)
    parallel::clusterExport(cluster,c("df_worker","total_weight_worker",
      "timegap_worker","difference_worker","full_names_worker"),
      envir=environment())
    worker_objective <- function(z) {
      full_z <- .piecewise_calendar_revision_monthly_difference_expand(z,
        difference_worker)[full_names_worker]
      p <- .piecewise_calendar_revision_monthly_unpack(full_z,timegap_worker)
      value <- tryCatch(e_step_eps(df_worker,p,check_df=FALSE,
        suff_stats=FALSE)$loglik,error=function(e) NA_real_)
      if (!is.finite(value)) return(1e100)
      -value/total_weight_worker
    }
    gradient <- function(z) {
      step <- gradient_step*pmax(1,abs(z))
      centre <- objective(z)
      if (identical(gradient_scheme,"forward")) {
        points <- lapply(seq_along(z),function(j) {
          point <- z; point[j] <- min(upper[j],point[j]+step[j]); point
        })
        plus <- unlist(parallel::parLapplyLB(cluster,points,
          worker_objective),use.names=FALSE)
        actual_plus <- vapply(seq_along(z),function(j)
          points[[j]][j]-z[j],numeric(1))
        valid_plus <- is.finite(plus) & plus<1e50 & actual_plus>0
        out <- numeric(length(z))
        out[valid_plus] <- (plus[valid_plus]-centre)/actual_plus[valid_plus]
        invalid <- which(!valid_plus)
        if (length(invalid)) {
          minus_points <- lapply(invalid,function(j) {
            point <- z; point[j] <- max(lower[j],point[j]-step[j]); point
          })
          minus <- unlist(parallel::parLapplyLB(cluster,minus_points,
            worker_objective),use.names=FALSE)
          actual_minus <- vapply(seq_along(invalid),function(k) {
            j <- invalid[k]; z[j]-minus_points[[k]][j]
          },numeric(1))
          valid_minus <- is.finite(minus) & minus<1e50 & actual_minus>0
          out[invalid[valid_minus]] <- (centre-minus[valid_minus])/
            actual_minus[valid_minus]
        }
        out[!is.finite(out)] <- 0
        return(pmin(pmax(out,-1e6),1e6))
      }
      points <- vector("list",2L*length(z))
      for (j in seq_along(z)) {
        points[[2L*j-1L]] <- points[[2L*j]] <- z
        points[[2L*j-1L]][j] <- min(upper[j],z[j]+step[j])
        points[[2L*j]][j] <- max(lower[j],z[j]-step[j])
      }
      values <- unlist(parallel::parLapplyLB(cluster,points,
        worker_objective),use.names=FALSE)
      plus <- values[seq(1L,length(values),by=2L)]
      minus <- values[seq(2L,length(values),by=2L)]
      plus_step <- vapply(seq_along(z),function(j)
        points[[2L*j-1L]][j]-z[j],numeric(1))
      minus_step <- vapply(seq_along(z),function(j)
        z[j]-points[[2L*j]][j],numeric(1))
      valid_plus <- is.finite(plus) & plus<1e50 & plus_step>0
      valid_minus <- is.finite(minus) & minus<1e50 & minus_step>0
      out <- numeric(length(z)); both <- valid_plus & valid_minus
      out[both] <- (plus[both]-minus[both])/
        (plus_step[both]+minus_step[both])
      plus_only <- valid_plus & !valid_minus
      out[plus_only] <- (plus[plus_only]-centre)/plus_step[plus_only]
      minus_only <- !valid_plus & valid_minus
      out[minus_only] <- (centre-minus[minus_only])/minus_step[minus_only]
      out[!is.finite(out)] <- 0
      out
    }
  }
  control <- list(maxit=maxit,reltol=reltol,
    ndeps=rep(1e-3,length(z0)),trace=as.integer(verbose))
  if (identical(method,"L-BFGS-B")) {
    control$factr <- reltol/.Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0,objective,gr=gradient,method=method,
      lower=lower,upper=upper,control=control)
  } else opt <- optim(z0,objective,gr=gradient,method=method,
    control=control)
  full_opt <- expand_z(opt$par)
  params <- .piecewise_calendar_revision_monthly_unpack(full_opt,
    timegap_contamination_model)
  estep <- e_step_eps(df,params,check_df=FALSE,suff_stats=FALSE)
  list(params=params,loglik=estep$loglik,gamma=estep$gamma,
    job_change_posterior=estep$job_change_posterior,
    convergence=opt$convergence,message=opt$message,
    iterations=unname(opt$counts["function"]),objective=opt$value,
    par_unconstrained=full_opt,profile_difference=difference,
    profile_reduced_par=opt$par,objective_function=objective)
}
