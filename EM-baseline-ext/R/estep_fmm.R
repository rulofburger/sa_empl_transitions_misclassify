# ==============================================================================
# EM-baseline-ext: E-step for Extension III (2-type FMM)
# Created: 2026-05-06
#
# Extends the baseline E-step to a 2-type finite mixture model. The latent
# space has 16 states: 8 histories x 2 types. Responsibilities are stored
# as two N x 8 matrices (gamma_A, gamma_B) rather than a single N x 16 matrix
# for cleaner downstream code.
#
# TeX ref: EM baseline.tex Section 7, Eqs (19)--(22).
# ==============================================================================

#' E-step for the 2-type FMM extension
#'
#' Computes joint history-type responsibilities
#' \eqn{\gamma_{ihk} = P(H_i = h, K_i = k \mid y_i; \Theta)}
#' for both types \eqn{k \in \{A, B\}}.
#'
#' TeX ref: EM baseline.tex Section 7, Eq (21).
#'
#' @param df Data frame with columns \code{y1}, \code{y2}, \code{y3}
#'   (integer 0/1) and \code{weight} (positive numeric). N rows.
#' @param params Named list with:
#'   \code{theta0_A}, \code{theta1_A}, \code{alpha_A} (type A Markov params);
#'   \code{theta0_B}, \code{theta1_B}, \code{alpha_B} (type B Markov params);
#'   \code{phi} (mixing weight for type A);
#'   and optionally \code{pi} (shared symmetric misclassification).
#' @param model_type Character: \code{"symmetric"} (default) or \code{"none"}.
#' @param validate Logical (default \code{TRUE}).
#' @return List with:
#'   \describe{
#'     \item{gamma_A}{N × 8 matrix: \eqn{P(H=h, K=A \mid y_i)}.}
#'     \item{gamma_B}{N × 8 matrix: \eqn{P(H=h, K=B \mid y_i)}.}
#'     \item{loglik}{Scalar weighted observed-data log-likelihood.}
#'     \item{suff}{Named list of sufficient statistics.}
#'   }
#' @examples
#' \dontrun{
#'   df     <- data.frame(y1 = rbinom(50,1,.6), y2 = rbinom(50,1,.6),
#'                        y3 = rbinom(50,1,.6), weight = rep(1,50))
#'   params <- init_params_fmm()
#'   e_step_fmm(df, params)
#' }
#' @export
e_step_fmm <- function(df, params, model_type = "symmetric", validate = TRUE,
                        mm_precomp = NULL) {
  if (!model_type %in% c("symmetric", "none"))
    stop("e_step_fmm: model_type must be 'symmetric' or 'none'")

  if (validate) {
    .validate_panel_df(df)
    required <- c("theta0_A", "theta1_A", "alpha_A",
                  "theta0_B", "theta1_B", "alpha_B", "phi")
    missing_p <- setdiff(required, names(params))
    if (length(missing_p) > 0L)
      stop(sprintf("e_step_fmm: params missing fields: %s",
                   paste(missing_p, collapse = ", ")))
    for (.col in c("y1", "y2", "y3")) {
      if (!all(df[[.col]] %in% c(0, 1)))
        stop(sprintf("e_step_fmm: column '%s' must be binary (0/1 only)", .col))
    }
  }

  N  <- nrow(df)
  w  <- df$weight
  s1 <- as.integer(df$y1)
  s2 <- as.integer(df$y2)
  s3 <- as.integer(df$y3)

  phi    <- .bound01(params$phi)
  pi     <- params$pi %||% 0

  hmat <- .hmat_cache()
  h1   <- hmat[, 1L]
  h2   <- hmat[, 2L]
  h3   <- hmat[, 3L]

  # ---- Type-specific priors (length-8 vectors) -----------------------------
  prior_A <- prior_over_histories(hmat,
                                  theta1 = .bound01(params$theta1_A),
                                  theta0 = .bound01(params$theta0_A),
                                  alpha  = .bound01(params$alpha_A))
  prior_B <- prior_over_histories(hmat,
                                  theta1 = .bound01(params$theta1_B),
                                  theta0 = .bound01(params$theta0_B),
                                  alpha  = .bound01(params$alpha_B))

  # Scale by mixture weight: phi * prior_A, (1-phi) * prior_B
  log_joint_A_prior <- log(phi)       + log(prior_A)  # length 8
  log_joint_B_prior <- log(1 - phi)   + log(prior_B)

  # ---- Log-emissions (shared across types) ---------------------------------
  log_emit_1 <- .log_misclass_wave_ext(s1, h1, model_type, pi)  # N x 8
  log_emit_2 <- .log_misclass_wave_ext(s2, h2, model_type, pi)
  log_emit_3 <- .log_misclass_wave_ext(s3, h3, model_type, pi)
  log_emit   <- log_emit_1 + log_emit_2 + log_emit_3              # N x 8

  # ---- Log-joint for each type (N x 8) ------------------------------------
  # log_joint_k[i, h] = log phi_k + log prior_k(h) + sum_t log q(s_{it}|h_t)
  log_jA <- sweep(log_emit, 2L, log_joint_A_prior, FUN = "+")  # N x 8
  log_jB <- sweep(log_emit, 2L, log_joint_B_prior, FUN = "+")  # N x 8

  # ---- Normalise over all 16 states ----------------------------------------
  # lapply avoids as.data.frame overhead for row-max
  log_max16   <- pmax(Reduce(pmax, lapply(seq_len(8L), function(j) log_jA[, j])),
                      Reduce(pmax, lapply(seq_len(8L), function(j) log_jB[, j])))
  log_row_sum <- log_max16 + log(
    rowSums(exp(log_jA - log_max16)) + rowSums(exp(log_jB - log_max16))
  )

  if (anyNA(log_row_sum))
    stop(paste0(
      "e_step_fmm: at least one observation has zero probability under all ",
      "16 (history x type) states. Check model parameters or model_type."
    ))

  # Responsibilities for each type (N x 8)
  gamma_A <- exp(log_jA - log_row_sum)
  gamma_B <- exp(log_jB - log_row_sum)

  loglik <- sum(w * log_row_sum)

  # ---- Sufficient statistics -----------------------------------------------
  wgA <- gamma_A * w   # N x 8: w_i * gamma_{ihA}
  wgB <- gamma_B * w

  # Transition indicator vectors (shared across types)
  from1_t2  <- as.integer(h1 == 1L)
  from0_t2  <- as.integer(h1 == 0L)
  to1_t2_1  <- as.integer(h1 == 1L & h2 == 1L)
  to1_t2_0  <- as.integer(h1 == 0L & h2 == 1L)
  from1_t3  <- as.integer(h2 == 1L)
  from0_t3  <- as.integer(h2 == 0L)
  to1_t3_1  <- as.integer(h2 == 1L & h3 == 1L)
  to1_t3_0  <- as.integer(h2 == 0L & h3 == 1L)

  .trans_suff <- function(wg) {
    wg_sum <- colSums(wg)
    list(
      T11 = sum(wg_sum * (to1_t2_1  + to1_t3_1)),
      T01 = sum(wg_sum * (to1_t2_0  + to1_t3_0)),
      D1  = sum(wg_sum * (from1_t2  + from1_t3)),
      D0  = sum(wg_sum * (from0_t2  + from0_t3)),
      C1  = sum(wg[, h1 == 1L, drop = FALSE]),
      C0  = sum(wg[, h1 == 0L, drop = FALSE])
    )
  }

  sA <- .trans_suff(wgA)
  sB <- .trans_suff(wgB)

  # Type-marginalised responsibility sums (for phi update)
  gamma_iA_sum <- rowSums(gamma_A)  # N-vector: P(K=A | y_i)
  W            <- sum(w)

  # Misclassification: M pooled across types
  # mm_precomp = outer(s1,h1,"!=") + outer(s2,h2,"!=") + outer(s3,h3,"!=")
  # can be precomputed once in the driver and passed in.
  if (model_type == "symmetric") {
    wg_all <- wgA + wgB  # N x 8: total responsibility weight
    if (!is.null(mm_precomp)) {
      M <- sum(wg_all * mm_precomp)
    } else {
      M <- sum(wg_all * outer(s1, h1, "!=")) +
           sum(wg_all * outer(s2, h2, "!=")) +
           sum(wg_all * outer(s3, h3, "!="))
    }
  } else {
    M <- 0
  }

  suff <- list(
    # Type A
    T11_A = sA$T11, T01_A = sA$T01, D1_A = sA$D1, D0_A = sA$D0,
    C1_A  = sA$C1,  C0_A  = sA$C0,
    # Type B
    T11_B = sB$T11, T01_B = sB$T01, D1_B = sB$D1, D0_B = sB$D0,
    C1_B  = sB$C1,  C0_B  = sB$C0,
    # Mixing weight
    phi_num      = sum(w * gamma_iA_sum),  # numerator of phi update
    total_weight = W,
    # Misclassification
    M = M
  )

  list(gamma_A = gamma_A, gamma_B = gamma_B, loglik = loglik, suff = suff)
}
