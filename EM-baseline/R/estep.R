# ==============================================================================
# EM-baseline: E-step — responsibilities and sufficient statistics
# Created: 2026-05-05
#
# Computes posterior probabilities gamma_{ih} = P(H_i = h | y_i; Theta)
# over the 8 latent employment histories, and accumulates the sufficient
# statistics needed by the M-step.
#
# TeX reference: EM baseline.tex, Eqs (12)-(16) for E-step and Eq (23)
# for the asymmetric extension.
# ==============================================================================

#' Compute log misclassification probabilities for a wave
#'
#' Returns the log probability log q(s_t | h_t) for each (observation,
#' history) pair at a single wave.
#'
#' @param s_t Integer vector of length N: observed employment state at wave t.
#' @param h_t Integer vector of length H: latent employment state at wave t
#'   (one element per history).
#' @param model_type Character: \code{"symmetric"}, \code{"asymmetric"}, or
#'   \code{"none"}.
#' @param pi Symmetric misclassification probability (used when
#'   \code{model_type = "symmetric"}).
#' @param pi0 False-positive rate P(s=1|h=0) (used when
#'   \code{model_type = "asymmetric"}).
#' @param pi1 False-negative rate P(s=0|h=1) (used when
#'   \code{model_type = "asymmetric"}).
#' @return N x H matrix of log q(s_t | h_t) values.
#' @keywords internal
.log_misclass_wave <- function(s_t, h_t, model_type, pi = 0, pi0 = 0, pi1 = 0) {
  .validate_model_type(model_type)  # secondary guard (nocov); primary validation is in e_step() # nolint
  N <- length(s_t)
  H <- length(h_t)

  # Create N x H matrix: h_mat[i,j] = h_t[j] (each row is h_t).
  # Vectorised ==: match[i,j] = (s_t[i] == h_t[j]) — no intermediate s_mat needed.
  h_mat <- matrix(h_t, nrow = N, ncol = H, byrow = TRUE)  # each row = h_t
  match <- (s_t == h_mat)

  if (model_type == "none") {
    # q(s|h) = 1{s=h}; log = 0 if match, -Inf if mismatch
    log_q <- ifelse(match, 0, -Inf)

  } else if (model_type == "symmetric") {
    # q(s|h) = (1-pi) if s=h, pi if s!=h
    pi <- .bound01(pi, eps = 1e-12)
    log_q <- ifelse(match, log(1 - pi), log(pi))

  } else if (model_type == "asymmetric") {
    # q(s=0|h=0) = 1-pi0, q(s=1|h=0) = pi0
    # q(s=0|h=1) = pi1,   q(s=1|h=1) = 1-pi1
    pi0 <- .bound01(pi0, eps = 1e-12)
    pi1 <- .bound01(pi1, eps = 1e-12)

    # h=0: correct (s=0) -> log(1-pi0); incorrect (s=1) -> log(pi0)
    # h=1: correct (s=1) -> log(1-pi1); incorrect (s=0) -> log(pi1)
    h_is_0 <- (h_mat == 0)
    log_q <- ifelse(h_is_0,
                    ifelse(match, log(1 - pi0), log(pi0)),
                    ifelse(match, log(1 - pi1), log(pi1)))
  } else {
    # Unreachable: .validate_model_type() above already stopped for unknown values.
    stop(sprintf("Unknown model_type: '%s'. Use 'symmetric', 'asymmetric', or 'none'.", model_type))  # nocov
  }

  log_q
}

#' E-step: compute responsibilities and sufficient statistics
#'
#' Given current parameter estimates, computes the posterior probability
#' gamma_{ih} that individual i has latent history h, conditional on the
#' observed data. Accumulates the sufficient statistics needed by
#' \code{\link{m_step}}.
#'
#' TeX ref: Eqs (12)-(16) for symmetric; Eq (23) for asymmetric extension.
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer, 0/1) and \code{weight} (positive numeric). N rows.
#' @param params Named list of current parameter estimates. Required fields:
#'   \code{alpha}, \code{theta0}, \code{theta1}. Optional fields:
#'   \code{pi} (symmetric), \code{pi0} and \code{pi1} (asymmetric).
#' @param model_type Character: \code{"symmetric"} (default),
#'   \code{"asymmetric"}, or \code{"none"}.
#' @param validate Logical (default \code{TRUE}). If \code{TRUE}, validates
#'   column types, NA values, binary range, and \code{params} values on each
#'   call. Set to \code{FALSE} only when data and params have been fully
#'   validated upstream (e.g., inside \code{\link{em_fit_baseline}} after the
#'   pre-loop checks).
#' @return List with:
#'   \describe{
#'     \item{gamma}{N x 8 matrix of posterior responsibilities. Rows sum to 1.}
#'     \item{loglik}{Scalar: observed-data weighted log-likelihood.}
#'     \item{suff}{Named list of sufficient statistics (see Details).}
#'   }
#'
#' @details
#' Sufficient statistics returned in \code{suff}:
#' \describe{
#'   \item{C1, C0}{Weighted mass at h1=1 vs h1=0. TeX Eq (13).}
#'   \item{D1, D0}{Person-periods from state 1 vs 0. TeX Eq (14).}
#'   \item{T11, T01}{Stay-employed and find-job transition counts. TeX Eq (15).}
#'   \item{M}{Total misclassified wave-observations (symmetric). TeX Eq (16).}
#'   \item{M0}{False positives (h=0, s=1) (asymmetric). TeX Eq (23).}
#'   \item{M1}{False negatives (h=1, s=0) (asymmetric). TeX Eq (23).}
#'   \item{H0}{Total latent-nonemployed wave-observations (denominator for pi0).}
#'   \item{H1}{Total latent-employed wave-observations (denominator for pi1).}
#'   \item{total_weight}{Sum of survey weights W = sum(w_i).}
#' }
#' @examples
#' df <- data.frame(
#'   y1 = c(1L, 0L, 1L), y2 = c(1L, 0L, 0L), y3 = c(1L, 1L, 0L),
#'   weight = c(1, 1, 1)
#' )
#' params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
#' result <- e_step(df, params, model_type = "symmetric")
#' result$loglik
#' # str(result$suff)  # inspect sufficient statistics
#' @export
e_step <- function(df, params, model_type = "symmetric", validate = TRUE) {
  # --- Validate inputs -------------------------------------------------------
  .validate_panel_df(df)
  .validate_model_type(model_type)

  # Per-value checks gated by validate. Set validate = FALSE in em_fit_baseline()
  # to avoid repeating these on every EM iteration (data is immutable between calls).
  if (validate) {
    # Factor columns: as.integer(factor) returns 1/2 codes, not 0/1
    for (.col in c("y1", "y2", "y3")) {
      if (is.factor(df[[.col]]))
        stop(sprintf(
          "Column '%s' is a factor; coerce to integer 0/1 before calling e_step()",
          .col
        ))
    }
    # NA checks (before binary check so NA rows get the right error message)
    if (anyNA(df$y1) || anyNA(df$y2) || anyNA(df$y3))
      stop("e_step: y1/y2/y3 must not contain NA values")
    if (anyNA(df$weight))
      stop("e_step: weight must not contain NA values")
    if (any(df$weight <= 0))
      stop("e_step: all weights must be strictly positive")
    # Binary check on raw pre-coercion values — catches non-integer numerics
    # (e.g. 0.5) that as.integer() would silently truncate to 0L or 1L.
    for (.col in c("y1", "y2", "y3")) {
      vals <- df[[.col]]
      if (!all(vals %in% c(0, 1)))
        stop(sprintf("e_step: column '%s' must be binary (0/1 only)", .col))
    }
    # Params validation: catch invalid values before .bound01() silently clamps them
    for (pname in c("alpha", "theta0", "theta1")) {
      pval <- params[[pname]]
      if (is.null(pval) || !is.numeric(pval) || length(pval) != 1L ||
          !is.finite(pval) || pval <= 0 || pval >= 1)
        stop(sprintf(
          "e_step: params$%s must be a single finite numeric in (0, 1); got: %s",
          pname,
          if (is.numeric(pval) && length(pval) == 1L && !is.na(pval))
            format(pval) else class(pval)
        ))
    }
  }

  N <- nrow(df)
  w <- df$weight

  # --- Extract params --------------------------------------------------------
  alpha  <- .bound01(params$alpha)
  theta1 <- .bound01(params$theta1)
  theta0 <- .bound01(params$theta0)

  pi  <- params$pi  %||% 0
  pi0 <- params$pi0 %||% 0
  pi1 <- params$pi1 %||% 0

  # --- Build latent histories and prior -------------------------------------
  hmat  <- latent_histories()    # 8 x 3
  prior <- prior_over_histories(hmat, theta1, theta0, alpha)  # length-8

  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]

  s1 <- as.integer(df$y1)
  s2 <- as.integer(df$y2)
  s3 <- as.integer(df$y3)

  # --- Compute N x 8 log-joint matrix ----------------------------------------
  # log_joint[i,j] = log P(h=j) + sum_t log q(s_{it} | h_{jt})
  log_prior <- log(prior)  # length-8; log(0) = -Inf for impossible histories

  log_emit_1 <- .log_misclass_wave(s1, h1, model_type, pi, pi0, pi1)  # N x 8
  log_emit_2 <- .log_misclass_wave(s2, h2, model_type, pi, pi0, pi1)  # N x 8
  log_emit_3 <- .log_misclass_wave(s3, h3, model_type, pi, pi0, pi1)  # N x 8

  # Each log_prior element is broadcast across N rows
  log_joint <- sweep(log_emit_1 + log_emit_2 + log_emit_3,
                     MARGIN = 2, STATS = log_prior, FUN = "+")
  # log_joint: N x 8

  # --- Normalise rows to get responsibilities --------------------------------
  # H=8 is a contract for the 3-wave model; assert to catch future extensions early.
  stopifnot(ncol(log_joint) == 8L)
  # Vectorised row logsumexp: avoids N R-level function calls per E-step.
  # Direct pmax() over H=8 named columns — no as.data.frame() copy of N×8 matrix.
  # Subtracting log_max before exp() keeps values in (-Inf, 0] — numerically safe.
  log_max     <- pmax(log_joint[, 1], log_joint[, 2], log_joint[, 3], log_joint[, 4],
                      log_joint[, 5], log_joint[, 6], log_joint[, 7], log_joint[, 8])
  log_row_sum <- log_max + log(rowSums(exp(log_joint - log_max)))
  # Guard: all-impossible observations (log_joint all -Inf) produce NaN silently.
  # Triggers when model_type='none' and data contain any observation-history mismatch.
  if (anyNA(log_row_sum))
    stop(paste0(
      "e_step: at least one observation has zero probability under all 8 latent ",
      "histories. Check for non-binary values in y1/y2/y3, or use model_type != 'none'."
    ))
  log_gamma   <- log_joint - log_row_sum           # N x 8 (log-scale)
  gamma       <- exp(log_gamma)                    # N x 8

  # --- Observed-data log-likelihood (survey-weighted) ------------------------
  loglik <- sum(w * log_row_sum)

  # --- Accumulate sufficient statistics -------------------------------------
  # Weighted responsibility: w_i * gamma_{ih} for each (i,h)
  # wg[i,j] = w[i] * gamma[i,j]
  wg <- gamma * w  # N x 8 (broadcasting w as column)

  # C: initial state mass  (TeX Eq 13)
  C1 <- sum(wg[, h1 == 1, drop = FALSE])
  C0 <- sum(wg[, h1 == 0, drop = FALSE])

  # Transition counts, summed over transitions t=2,3  (TeX Eqs 14-15)
  # D1: person-periods where h_{t-1} = 1  (from state: employed)
  # D0: person-periods where h_{t-1} = 0  (from state: nonemployed)
  # T11: h_{t-1}=1 -> h_t=1
  # T01: h_{t-1}=0 -> h_t=1
  # For 3-wave data there are 2 transitions: (h1->h2) and (h2->h3)
  # Indicator vectors for the 8 histories: whether each from/to combination occurs
  from1_to2_is_11 <- as.integer(h1 == 1 & h2 == 1)  # T11 contribution from t=2
  from1_to2_is_01 <- as.integer(h1 == 0 & h2 == 1)  # T01 contribution from t=2
  from1_to2_from1 <- as.integer(h1 == 1)             # D1  contribution from t=2
  from1_to2_from0 <- as.integer(h1 == 0)             # D0  contribution from t=2

  from2_to3_is_11 <- as.integer(h2 == 1 & h3 == 1)
  from2_to3_is_01 <- as.integer(h2 == 0 & h3 == 1)
  from2_to3_from1 <- as.integer(h2 == 1)
  from2_to3_from0 <- as.integer(h2 == 0)

  # colSums(wg) = sum_i w_i * gamma_{ih}: length-8 vector
  wg_sum <- colSums(wg)

  T11 <- sum(wg_sum * (from1_to2_is_11 + from2_to3_is_11))
  T01 <- sum(wg_sum * (from1_to2_is_01 + from2_to3_is_01))
  D1  <- sum(wg_sum * (from1_to2_from1 + from2_to3_from1))
  D0  <- sum(wg_sum * (from1_to2_from0 + from2_to3_from0))

  # Misclassification counts  (TeX Eqs 16, 23)
  # Gate M computation to symmetric only: m_step() reads suff$M only in the
  # symmetric branch (model_type='asymmetric' uses M0/M1 instead; 'none' uses neither).
  # This avoids 6 outer() N×8 allocations on every asymmetric E-step call.
  if (model_type == "symmetric") {
    mm1 <- outer(s1, h1, "!=")  # N x 8: s1 != h1
    mm2 <- outer(s2, h2, "!=")  # N x 8: s2 != h2
    mm3 <- outer(s3, h3, "!=")  # N x 8: s3 != h3
    # Separate sums avoid an additional N×8 intermediate matrix
    M <- sum(wg * mm1) + sum(wg * mm2) + sum(wg * mm3)
  } else {
    M <- 0  # unused for 'asymmetric'; zero for 'none'
  }

  if (model_type == "asymmetric") {
    # Asymmetric: separate counts by true state (TeX Eq 23).
    # Index-based approach: h_*_zero_idx / h_*_one_idx are length-4 fixed integer
    # vectors (same for every call). Avoids allocating 9 N×8 logical matrices.
    h1_zero_idx <- which(h1 == 0)   # column indices of histories where h1=0 (length 4)
    h2_zero_idx <- which(h2 == 0)
    h3_zero_idx <- which(h3 == 0)
    h1_one_idx  <- which(h1 == 1)   # column indices where h1=1 (length 4)
    h2_one_idx  <- which(h2 == 1)
    h3_one_idx  <- which(h3 == 1)

    # M0: false positives — h=0 and s=1
    M0 <- sum(wg[s1 == 1, h1_zero_idx, drop = FALSE]) +
          sum(wg[s2 == 1, h2_zero_idx, drop = FALSE]) +
          sum(wg[s3 == 1, h3_zero_idx, drop = FALSE])

    # M1: false negatives — h=1 and s=0
    M1 <- sum(wg[s1 == 0, h1_one_idx, drop = FALSE]) +
          sum(wg[s2 == 0, h2_one_idx, drop = FALSE]) +
          sum(wg[s3 == 0, h3_one_idx, drop = FALSE])

    # H0, H1: total latent wave-observations by true state (denominators for pi0, pi1)
    H0 <- sum(wg_sum * (as.integer(h1 == 0) + as.integer(h2 == 0) + as.integer(h3 == 0)))
    H1 <- sum(wg_sum * (as.integer(h1 == 1) + as.integer(h2 == 1) + as.integer(h3 == 1)))
  } else {
    M0 <- 0; M1 <- 0; H0 <- 0; H1 <- 0
  }

  suff <- list(
    C1           = C1,
    C0           = C0,
    D1           = D1,
    D0           = D0,
    T11          = T11,
    T01          = T01,
    M            = M,
    M0           = M0,
    M1           = M1,
    H0           = H0,
    H1           = H1,
    total_weight = sum(w)
  )

  list(gamma = gamma, loglik = loglik, suff = suff)
}
