# ==============================================================================
# EM-baseline-ext: E-step for Extension I (observable heterogeneity)
# Created: 2026-05-06
#
# Computes responsibilities gamma_{ih} and sufficient statistics for the
# covariate extension of the baseline EM model (TeX Section 5, Eq 11).
#
# Key difference from baseline e_step(): each individual has their own
# Markov prior computed from individual-specific transition probabilities
# theta1_i = Phi(X_i' beta1), theta0_i = Phi(X_i' beta0).
#
# The misclassification block is unchanged: scalar pi applied uniformly.
# ==============================================================================

#' E-step for the covariate (observable heterogeneity) extension
#'
#' Computes posterior responsibilities \eqn{\gamma_{ih}} and sufficient
#' statistics for the GEM M-step, using individual-specific Markov priors
#' derived from probit-linked covariates.
#'
#' TeX ref: EM baseline.tex Section 5, Eqs (10)--(11).
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric). N rows.
#' @param X N × p numeric design matrix from
#'   \code{\link{prepare_covariate_matrix}}.
#' @param params Named list with:
#'   \code{beta0} (length-p numeric), \code{beta1} (length-p numeric),
#'   and optionally \code{pi} (scalar, used when \code{model_type = "symmetric"}).
#' @param model_type Character: \code{"symmetric"} (default) or \code{"none"}.
#' @param validate Logical (default \code{TRUE}). Set \code{FALSE} inside
#'   the EM loop after initial validation.
#' @return List with:
#'   \describe{
#'     \item{gamma}{N × 8 matrix of posterior responsibilities (rows sum to 1).}
#'     \item{loglik}{Scalar weighted observed-data log-likelihood.}
#'     \item{suff}{Named list of sufficient statistics for \code{\link{m_step_covariates}}.}
#'   }
#' @details
#' Sufficient statistics in \code{suff}:
#' \describe{
#'   \item{XtWy_0, XtWX_0}{Cross-products for GEM probit update of \eqn{\beta_0}
#'     (job-finding from nonemployment). Accumulated from person-periods where
#'     \eqn{h_{t-1}=0}.}
#'   \item{XtWy_1, XtWX_1}{Cross-products for GEM probit update of \eqn{\beta_1}
#'     (employment persistence). Accumulated from person-periods where
#'     \eqn{h_{t-1}=1}.}
#'   \item{eff_w_0, eff_w_1}{Total effective weight at each individual for
#'     transitions from state 0 vs state 1 (used for IRLS normalisation).}
#'   \item{M}{Total responsibility-weighted misclassifications (for \eqn{\pi}).}
#'   \item{C1, C0}{Initial-state responsibility mass (for free \eqn{\alpha}).}
#'   \item{total_weight}{Sum of survey weights W.}
#' }
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   X      <- matrix(1, nrow = 50, ncol = 1)
#'   params <- init_params_covariates(p = 1L)
#'   e_step_covariates(df, X, params)
#' }
#' @export
e_step_covariates <- function(df, X, params, model_type = "symmetric",
                               validate = TRUE, mm_precomp = NULL,
                               stationary = TRUE) {
  if (!model_type %in% c("symmetric", "none"))
    stop("e_step_covariates: model_type must be 'symmetric' or 'none'")

  N <- nrow(df)
  X_transition <- .as_transition_design(X, N)
  if (stationary && !.transition_design_is_time_invariant(X_transition))
    stop("e_step_covariates: stationarity is undefined with time-varying transition covariates; use stationary=FALSE")

  if (validate) {
    .validate_panel_df(df)
    if (is.null(params$beta0) || is.null(params$beta1))
      stop("e_step_covariates: params must contain 'beta0' and 'beta1'")
    if (length(params$beta0) != ncol(X_transition$X12) ||
        length(params$beta1) != ncol(X_transition$X12))
      stop("e_step_covariates: length(beta0) and length(beta1) must equal ncol(X)")
    for (.col in c("y1", "y2", "y3")) {
      if (!all(df[[.col]] %in% c(0, 1)))
        stop(sprintf("e_step_covariates: column '%s' must be binary (0/1 only)", .col))
    }
  }

  w  <- df$weight
  p  <- ncol(X_transition$X12)
  s1 <- as.integer(df$y1)
  s2 <- as.integer(df$y2)
  s3 <- as.integer(df$y3)

  beta0 <- params$beta0
  beta1 <- params$beta1
  pi    <- params$pi %||% 0

  # ---- Individual-specific transition probabilities ------------------------
  theta0_12 <- pmin(pmax(pnorm(as.vector(X_transition$X12 %*% beta0)), 1e-6), 1 - 1e-6)
  theta0_23 <- pmin(pmax(pnorm(as.vector(X_transition$X23 %*% beta0)), 1e-6), 1 - 1e-6)
  theta1_12 <- pmin(pmax(pnorm(as.vector(X_transition$X12 %*% beta1)), 1e-6), 1 - 1e-6)
  theta1_23 <- pmin(pmax(pnorm(as.vector(X_transition$X23 %*% beta1)), 1e-6), 1 - 1e-6)

  # Initial employment probability. Under stationarity it is implied by the
  # individual transition probabilities; otherwise alpha is a free scalar.
  alpha_i <- if (stationary) {
    theta0_12 / (theta0_12 + 1 - theta1_12)
  } else {
    if (is.null(params$alpha) || !is.finite(params$alpha))
      stop("e_step_covariates: free-alpha model requires finite params$alpha")
    rep(params$alpha, N)
  }
  alpha_i <- pmin(pmax(alpha_i, 1e-6), 1 - 1e-6)

  # ---- Latent histories (8 x 3) -------------------------------------------
  hmat <- .hmat_cache()
  h1   <- hmat[, 1L]
  h2   <- hmat[, 2L]
  h3   <- hmat[, 3L]

  # ---- Individual-specific N x 8 prior matrix ------------------------------
  # For individual i and history h:
  # log P(h | theta0_i, theta1_i) = log P(h1) + log P(h2|h1) + log P(h3|h2)
  # Vectorised: use outer products.

  # Log initial probability: P(h1=1)=alpha_i, P(h1=0)=1-alpha_i
  # log_p_h1[i, j] = log P(h1[j] | alpha_i[i])
  log_p_h1 <- matrix(NA_real_, nrow = N, ncol = 8L)
  log_p_h1[, h1 == 1L] <- log(alpha_i)
  log_p_h1[, h1 == 0L] <- log(1 - alpha_i)

  # Log transition 1->2 for each individual
  log_p_12 <- .log_markov_trans_indiv(h1, h2, theta1_12, theta0_12)
  log_p_23 <- .log_markov_trans_indiv(h2, h3, theta1_23, theta0_23)

  log_prior <- log_p_h1 + log_p_12 + log_p_23  # N x 8

  # ---- Log-emission probabilities ------------------------------------------
  log_emit_1 <- .log_misclass_wave_ext(s1, h1, model_type, pi)  # N x 8
  log_emit_2 <- .log_misclass_wave_ext(s2, h2, model_type, pi)
  log_emit_3 <- .log_misclass_wave_ext(s3, h3, model_type, pi)

  # ---- Log-joint and responsibilities (N x 8) ------------------------------
  log_joint <- log_prior + log_emit_1 + log_emit_2 + log_emit_3

  # Row logsumexp (H=8 columns) — lapply avoids as.data.frame overhead
  log_max     <- Reduce(pmax, lapply(seq_len(8L), function(j) log_joint[, j]))
  log_row_sum <- log_max + log(rowSums(exp(log_joint - log_max)))
  if (anyNA(log_row_sum))
    stop(paste0(
      "e_step_covariates: at least one observation has zero probability ",
      "under all 8 histories. Check model_type or covariate values."
    ))

  log_gamma <- log_joint - log_row_sum   # N x 8
  gamma     <- exp(log_gamma)

  loglik <- sum(w * log_row_sum)

  # ---- Sufficient statistics for GEM M-step --------------------------------
  wg     <- gamma * w  # N x 8
  wg_sum <- colSums(wg)  # length 8

  # For GEM probit update, accumulate IRLS cross-products.
  # For each transition (h_{t-1}, h_t), the effective weight and outcome
  # at individual i for beta1 (from h=1) or beta0 (from h=0).
  #
  # Effective weight for beta_k at (i, t):
  #   eff_w_k[i] = sum_{h: h_{t-1}=k} w_i * gamma_{ih}   (summed over t=2,3)
  # Effective outcome (fractional):
  #   eff_y_k[i] = sum_{h: h_{t-1}=k, h_t=1} w_i * gamma_{ih} / eff_w_k[i]

  # Indicator vectors over 8 histories (used at both transitions)
  from1_t2  <- as.integer(h1 == 1L)   # h1=1: from employment at t=2
  from0_t2  <- as.integer(h1 == 0L)
  to1_t2_1  <- as.integer(h1 == 1L & h2 == 1L)   # from 1, go to 1 at t=2
  to1_t2_0  <- as.integer(h1 == 0L & h2 == 1L)

  from1_t3  <- as.integer(h2 == 1L)
  from0_t3  <- as.integer(h2 == 0L)
  to1_t3_1  <- as.integer(h2 == 1L & h3 == 1L)
  to1_t3_0  <- as.integer(h2 == 0L & h3 == 1L)

  # Keep interval-specific sufficient statistics because origin-wave contract
  # type and sector can differ between transitions.
  eff_w_1_12  <- as.vector(wg %*% from1_t2)
  eff_w_0_12  <- as.vector(wg %*% from0_t2)
  eff_wy_1_12 <- as.vector(wg %*% to1_t2_1)
  eff_wy_0_12 <- as.vector(wg %*% to1_t2_0)
  eff_w_1_23  <- as.vector(wg %*% from1_t3)
  eff_w_0_23  <- as.vector(wg %*% from0_t3)
  eff_wy_1_23 <- as.vector(wg %*% to1_t3_1)
  eff_wy_0_23 <- as.vector(wg %*% to1_t3_0)

  eff_w_1 <- eff_w_1_12 + eff_w_1_23
  eff_w_0 <- eff_w_0_12 + eff_w_0_23
  eff_wy_1 <- eff_wy_1_12 + eff_wy_1_23
  eff_wy_0 <- eff_wy_0_12 + eff_wy_0_23

  # Misclassification count M (for pi update)
  # mm_precomp = outer(s1,h1,"!=") + outer(s2,h2,"!=") + outer(s3,h3,"!=")
  # can be precomputed once in the driver and passed in.
  if (model_type == "symmetric") {
    if (!is.null(mm_precomp)) {
      M <- sum(wg * mm_precomp)
    } else {
      M <- sum(wg * outer(s1, h1, "!=")) +
           sum(wg * outer(s2, h2, "!=")) +
           sum(wg * outer(s3, h3, "!="))
    }
  } else {
    M <- 0
  }

  # Initial-state posterior mass.  Retain individual contributions because the
  # stationary initial probability alpha_i depends jointly on beta0 and beta1.
  init_w1 <- as.vector(wg %*% as.integer(h1 == 1L))
  init_w0 <- as.vector(wg %*% as.integer(h1 == 0L))
  C1 <- sum(init_w1)
  C0 <- sum(init_w0)

  suff <- list(
    # GEM probit cross-products — computed lazily in m_step_covariates
    eff_w_1  = eff_w_1,
    eff_wy_1 = eff_wy_1,
    eff_w_0  = eff_w_0,
    eff_wy_0 = eff_wy_0,
    transition = list(
      X12 = list(eff_w_1 = eff_w_1_12, eff_wy_1 = eff_wy_1_12,
                 eff_w_0 = eff_w_0_12, eff_wy_0 = eff_wy_0_12),
      X23 = list(eff_w_1 = eff_w_1_23, eff_wy_1 = eff_wy_1_23,
                 eff_w_0 = eff_w_0_23, eff_wy_0 = eff_wy_0_23)
    ),
    M            = M,
    C1           = C1,
    C0           = C0,
    init_w1      = init_w1,
    init_w0      = init_w0,
    total_weight = sum(w)
  )

  list(gamma = gamma, loglik = loglik, suff = suff)
}


# NOTE: .log_markov_trans_indiv() and .log_misclass_wave_ext() are defined in
# helpers_ext.R (sourced first) and shared across all extension E-steps.
