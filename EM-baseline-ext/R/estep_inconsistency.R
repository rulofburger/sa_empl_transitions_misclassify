# ==============================================================================
# EM-baseline-ext: E-step for Extension IV (inconsistency-augmented model)
# Created: 2026-05-06
#
# Models wave-specific misclassification probability as
#   pi_it = 0.5 * sigma(delta_0 + delta_1 * Y_age_it + delta_2 * Y_edu_it)
# where sigma() is the logistic function, ensuring pi_it in (0, 0.5).
#
# TeX ref: EM baseline.tex Section 8, Eqs (26)--(30).
# ==============================================================================

#' Compute posterior mismatch probabilities from gamma and observed states
#'
#' Internal helper. For each individual and wave, sums gamma over histories
#' where the latent state disagrees with the observed state.
#'
#' @param gamma N x 8 responsibility matrix.
#' @param s1,s2,s3 Integer 0/1 vectors of length N.
#' @param hmat 8 x 3 matrix of latent histories (from \code{latent_histories()}).
#' @return N x 3 matrix of posterior mismatch probabilities.
.compute_p_mat <- function(gamma, s1, s2, s3, hmat) {
  h1 <- hmat[, 1L]
  h2 <- hmat[, 2L]
  h3 <- hmat[, 3L]
  # For each wave t: P(h_t != s_it | y_i) = sum_{h: h_t != s_it} gamma_ih
  .pmis <- function(st, ht) {
    # indices of histories where h_t == 0 and h_t == 1
    idx0 <- which(ht == 0L)
    idx1 <- which(ht == 1L)
    # p_mismatch_i = sum_{h: h_t != s_it} gamma_ih
    # = I(s_it == 1) * rowSums(gamma[, idx0]) + I(s_it == 0) * rowSums(gamma[, idx1])
    as.numeric(
      (st == 1L) * rowSums(gamma[, idx0, drop = FALSE]) +
      (st == 0L) * rowSums(gamma[, idx1, drop = FALSE])
    )
  }
  cbind(
    p1 = .pmis(s1, h1),
    p2 = .pmis(s2, h2),
    p3 = .pmis(s3, h3)
  )
}


#' Compute individual- and wave-specific log-emissions under inconsistency model
#'
#' Internal helper.  Returns an N x 8 matrix for one wave.
#'
#' @param st Integer 0/1 vector (length N): observed employment at wave t.
#' @param ht Integer 0/1 vector (length 8): latent state for all histories at t.
#' @param pi_t Numeric vector (length N): individual misclassification probs at t.
#' @return N x 8 matrix of log P(s_it | h_t, pi_it).
.log_emit_indiv <- function(st, ht, pi_t) {
  # Clamp pi_t for numerical safety
  pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8)
  # match_mask[i, h]: 1 if s_it == h_t, 0 otherwise
  match_mask <- outer(st, ht, "==")
  # log P = match_mask * log(1-pi_t) + (1-match_mask) * log(pi_t)
  match_mask * log(1 - pi_t) + (1 - match_mask) * log(pi_t)
}


#' E-step for the inconsistency-augmented model
#'
#' Uses individual- and wave-specific misclassification probabilities
#' \eqn{\pi_{it} = \tfrac{1}{2}\,\sigma(\delta_0 + \delta_1 Y_{it}^{\text{age}} +
#' \delta_2 Y_{it}^{\text{edu}})}.
#'
#' TeX ref: EM baseline.tex Section 8, Eq (29).
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric).
#' @param incons_mat N × 6 numeric matrix with columns
#'   \code{Y_age_1}, \code{Y_age_2}, \code{Y_age_3},
#'   \code{Y_edu_1}, \code{Y_edu_2}, \code{Y_edu_3}.
#'   Typically produced by \code{\link{compute_inconsistencies}}.
#' @param params Named list with:
#'   \code{theta0}, \code{theta1}, \code{alpha} (Markov parameters);
#'   \code{delta} (length-3 vector: intercept, age coefficient, edu coefficient).
#' @param validate Logical (default \code{TRUE}).
#' @return List with:
#'   \describe{
#'     \item{gamma}{N × 8 responsibility matrix.}
#'     \item{loglik}{Scalar weighted observed-data log-likelihood.}
#'     \item{suff}{Named list of sufficient statistics including \code{p_mat}.}
#'   }
#' @examples
#' \dontrun{
#'   df      <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                         y3 = rbinom(50,1,.6), weight = rep(1,50),
#'                         age1=25:74, age2=26:75, age3=27:76,
#'                         educ1=rep(3L,50), educ2=rep(3L,50), educ3=rep(3L,50))
#'   inc_mat <- as.matrix(compute_inconsistencies(df)[,
#'     c("Y_age_1","Y_age_2","Y_age_3","Y_edu_1","Y_edu_2","Y_edu_3")])
#'   e_step_inconsistency(df, inc_mat, init_params_inconsistency())
#' }
#' @export
e_step_inconsistency <- function(df, incons_mat, params, validate = TRUE) {
  if (validate) {
    .validate_panel_df(df)
    required_p <- c("theta0", "theta1", "alpha", "delta")
    missing_p  <- setdiff(required_p, names(params))
    if (length(missing_p) > 0L)
      stop(sprintf("e_step_inconsistency: params missing: %s",
                   paste(missing_p, collapse = ", ")))
    if (length(params$delta) != 3L)
      stop("e_step_inconsistency: params$delta must have length 3")
    if (!is.matrix(incons_mat) || ncol(incons_mat) != 6L)
      stop("e_step_inconsistency: incons_mat must be an N x 6 matrix")
    if (nrow(incons_mat) != nrow(df))
      stop("e_step_inconsistency: incons_mat must have the same number of rows as df")
    for (.col in c("y1", "y2", "y3")) {
      if (!all(df[[.col]] %in% c(0, 1)))
        stop(sprintf("e_step_inconsistency: column '%s' must be binary", .col))
    }
    # P1: validate incons_mat column order (wrong order silently corrupts eta_mat)
    if (!is.null(colnames(incons_mat))) {
      expected_cols <- c("Y_age_1", "Y_age_2", "Y_age_3",
                         "Y_edu_1", "Y_edu_2", "Y_edu_3")
      if (!identical(colnames(incons_mat), expected_cols))
        stop(sprintf(
          "e_step_inconsistency: incons_mat columns must be ordered: %s",
          paste(expected_cols, collapse = ", ")
        ))
    }
    if (!all(incons_mat %in% c(0L, 1L, 0, 1)))
      stop("e_step_inconsistency: incons_mat must be binary (0/1 only)")
  }

  w  <- df$weight
  s1 <- as.integer(df$y1)
  s2 <- as.integer(df$y2)
  s3 <- as.integer(df$y3)

  delta <- params$delta
  Y_age <- incons_mat[, 1:3, drop = FALSE]  # columns 1-3: Y_age_1, Y_age_2, Y_age_3
  Y_edu <- incons_mat[, 4:6, drop = FALSE]  # columns 4-6: Y_edu_1, Y_edu_2, Y_edu_3

  # ---- Individual- and wave-specific pi_it ---------------------------------
  # eta_it = delta_0 + delta_1 * Y_age_it + delta_2 * Y_edu_it  (N x 3)
  eta_mat <- delta[1L] + delta[2L] * Y_age + delta[3L] * Y_edu
  sig_mat <- plogis(eta_mat)                  # N x 3: sigma(eta_it)
  pi_mat  <- 0.5 * sig_mat                    # N x 3: pi_it in (0, 0.5)

  # ---- Latent history prior ------------------------------------------------
  hmat    <- .hmat_cache()
  prior_h <- prior_over_histories(hmat,
                                  theta1 = .bound01(params$theta1),
                                  theta0 = .bound01(params$theta0),
                                  alpha  = .bound01(params$alpha))
  log_prior <- log(prior_h)  # length 8

  # ---- Log-emissions (N x 8 per wave, with individual pi_it) --------------
  log_emit_1 <- .log_emit_indiv(s1, hmat[, 1L], pi_mat[, 1L])
  log_emit_2 <- .log_emit_indiv(s2, hmat[, 2L], pi_mat[, 2L])
  log_emit_3 <- .log_emit_indiv(s3, hmat[, 3L], pi_mat[, 3L])
  log_emit   <- log_emit_1 + log_emit_2 + log_emit_3  # N x 8

  # ---- Log joint + normalisation -------------------------------------------
  log_joint <- sweep(log_emit, 2L, log_prior, FUN = "+")  # N x 8
  # Row logsumexp (H=8 columns) — lapply avoids as.data.frame overhead
  log_max   <- Reduce(pmax, lapply(seq_len(8L), function(j) log_joint[, j]))
  log_rs    <- log_max + log(rowSums(exp(log_joint - log_max)))

  if (anyNA(log_rs))
    stop(paste0(
      "e_step_inconsistency: at least one observation has zero probability ",
      "under all 8 histories. Check model parameters."
    ))

  gamma  <- exp(log_joint - log_rs)
  loglik <- sum(w * log_rs)

  # ---- Sufficient statistics -----------------------------------------------
  wg <- gamma * w  # N x 8
  h1 <- hmat[, 1L]
  h2 <- hmat[, 2L]
  h3 <- hmat[, 3L]

  wg_col <- colSums(wg)  # length 8: avoids 4 N×8 outer-product allocations
  suff <- list(
    T11 = sum(wg_col * (as.integer(h1 == 1L & h2 == 1L) + as.integer(h2 == 1L & h3 == 1L))),
    T01 = sum(wg_col * (as.integer(h1 == 0L & h2 == 1L) + as.integer(h2 == 0L & h3 == 1L))),
    D1  = sum(wg_col * (as.integer(h1 == 1L) + as.integer(h2 == 1L))),
    D0  = sum(wg_col * (as.integer(h1 == 0L) + as.integer(h2 == 0L))),
    C1  = sum(wg[, h1 == 1L, drop = FALSE]),
    C0  = sum(wg[, h1 == 0L, drop = FALSE]),
    # Posterior mismatch probabilities (N x 3) — used for NR delta update
    p_mat        = .compute_p_mat(gamma, s1, s2, s3, hmat),
    # Diagnostic: pi_mat and sig_mat (not consumed by m_step; useful for review)
    pi_mat       = pi_mat,
    sig_mat      = sig_mat,
    # Survey weights — used by m_step_inconsistency for NR step
    weights      = w,
    total_weight = sum(w)
  )

  list(gamma = gamma, loglik = loglik, suff = suff)
}
