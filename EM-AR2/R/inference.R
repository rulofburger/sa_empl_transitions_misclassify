# ==============================================================================
# EM-AR2: Post-estimation inference — cell probabilities, transitions, GoF
# Created: 2026-05-05
# ==============================================================================
# After EM convergence, this module computes:
#   1. Model-implied 16 cell probabilities P(y1,y2,y3,y4 | params)
#   2. Implied true transition probabilities (comparable to MLE pipeline output)
#   3. Goodness-of-fit: empirical vs model cell probabilities
#
# The approach mirrors the `df_template` + `group_by` + `summarise` pattern
# from scripts/estimate_models_4waves_SA.R, using the 256 latent
# (y*1,y*2,y*3,y*4) x (y1,y2,y3,y4) combinations.
# ==============================================================================

#' Compute model-implied cell probabilities for all 16 observed sequences
#'
#' Marginalises over all 16 latent histories to compute
#' P(y1,y2,y3,y4 | params) = sum_{h} P(h|params) * P(y|h,pi).
#'
#' @param params Named list with theta0, theta01, theta1, theta10, pi
#'   (and optionally pi0, pi1 if asymmetric).
#' @return Data frame with columns y1, y2, y3, y4, model_prob.
#'   Rows are all 16 observed binary sequences. model_prob sums to 1.
#' @export
model_cell_probs_ar2 <- function(params) {
  hmat  <- latent_histories_ar2()
  prior <- prior_over_histories_ar2(hmat,
    theta0  = params$theta0,
    theta01 = params$theta01,
    theta1  = params$theta1,
    theta10 = params$theta10,
    alpha   = params$alpha
  )

  asymmetric <- isTRUE(params$asymmetric)
  pi  <- params$pi  %||% 0
  pi0 <- params$pi0 %||% pi
  pi1 <- params$pi1 %||% pi

  # All 16 observed sequences
  obs_seqs <- as.matrix(expand.grid(
    y1 = 0:1, y2 = 0:1, y3 = 0:1, y4 = 0:1
  ))

  n_obs <- nrow(obs_seqs)  # 16

  model_prob <- numeric(n_obs)

  for (i in seq_len(n_obs)) {
    y <- obs_seqs[i, ]

    # For each latent history, compute emission probability
    p_y_given_h <- numeric(nrow(hmat))
    for (j in seq_len(nrow(hmat))) {
      h <- hmat[j, ]
      log_p <- 0
      for (t in 1:4) {
        if (!asymmetric) {
          log_p <- log_p +
            if (y[t] == h[t]) log(1 - pi) else log(pi)
        } else {
          if (h[t] == 1) {
            log_p <- log_p +
              if (y[t] == 1) log(1 - pi0) else log(pi0)
          } else {
            log_p <- log_p +
              if (y[t] == 0) log(1 - pi1) else log(pi1)
          }
        }
      }
      p_y_given_h[j] <- exp(log_p)
    }

    # Marginalise over histories
    model_prob[i] <- sum(prior * p_y_given_h)
  }

  data.frame(obs_seqs, model_prob = model_prob)
}

#' Compute implied true (latent) transition probabilities
#'
#' Computes four summary transition quantities for the AR(2) model,
#' matching the MLE pipeline output in estimate_models_4waves_SA.R:
#'   1. P(y*4=1 | y*1=0) — probability employed at wave 4 given non-employed at wave 1
#'   2. P(ever employed | y*1=0) — probability of employment at waves 2, 3, or 4
#'   3. P(y*4=0 | y*1=1) — probability non-employed at wave 4 given employed at wave 1
#'   4. P(ever non-employed | y*1=1) — probability of non-employment at waves 2, 3, or 4
#'
#' @param params Named list with theta0, theta01, theta1, theta10.
#' @return Named numeric vector with four elements:
#'   p_emp_w4_from0, p_ever_emp_from0, p_nonemp_w4_from1, p_ever_nonemp_from1.
#' @examples
#' \dontrun{
#' params <- list(theta0=0.10, theta01=0.15, theta1=0.08, theta10=0.12, pi=0)
#' implied_transitions_ar2(params)
#' }
#' @export
implied_transitions_ar2 <- function(params) {
  hmat  <- latent_histories_ar2()
  # Use params$alpha if provided (e.g. during EM diagnostics or testing with a
  # custom initial distribution). For standard post-estimation inference, pass
  # params without alpha (or alpha=NULL) to use the stationary distribution,
  # which reflects model-implied long-run transitions.
  prior <- prior_over_histories_ar2(hmat,
    theta0  = params$theta0,
    theta01 = params$theta01,
    theta1  = params$theta1,
    theta10 = params$theta10,
    alpha   = params$alpha  # NULL → stationary; non-NULL → use as-is
  )

  h1 <- hmat[, 1]
  h2 <- hmat[, 2]
  h3 <- hmat[, 3]
  h4 <- hmat[, 4]

  denom0 <- sum(prior[h1 == 0])
  denom1 <- sum(prior[h1 == 1])

  if (denom0 < 1e-12) stop("implied_transitions_ar2: no prior mass on h1=0. Model may be degenerate.")
  if (denom1 < 1e-12) stop("implied_transitions_ar2: no prior mass on h1=1. Model may be degenerate.")

  # P(y*4=1 | y*1=0)
  # = sum_{h: h1=0, h4=1} P(h) / sum_{h: h1=0} P(h)
  p_emp_w4_from0 <- sum(prior[h1 == 0 & h4 == 1]) / denom0

  # P(ever employed | y*1=0)
  # = P(h2=1 OR h3=1 OR h4=1 | h1=0)
  ever_emp <- (h2 == 1) | (h3 == 1) | (h4 == 1)
  p_ever_emp_from0 <- sum(prior[h1 == 0 & ever_emp]) / denom0

  # P(y*4=0 | y*1=1)
  p_nonemp_w4_from1 <- sum(prior[h1 == 1 & h4 == 0]) / denom1

  # P(ever non-employed | y*1=1)
  ever_nonemp <- (h2 == 0) | (h3 == 0) | (h4 == 0)
  p_ever_nonemp_from1 <- sum(prior[h1 == 1 & ever_nonemp]) / denom1

  c(
    p_emp_w4_from0      = p_emp_w4_from0,
    p_ever_emp_from0    = p_ever_emp_from0,
    p_nonemp_w4_from1   = p_nonemp_w4_from1,
    p_ever_nonemp_from1 = p_ever_nonemp_from1
  )
}

#' Compute goodness-of-fit by comparing empirical and model cell probabilities
#'
#' Computes the 16 empirical cell proportions (weighted) and the 16 model-
#' implied cell probabilities, returning a data frame with both and the
#' residuals. Also returns the sum of squared residuals (SSR) as an overall
#' fit statistic.
#'
#' @param df Data frame with y1, y2, y3, y4, weight columns.
#' @param params Named list of model parameters.
#' @return Named list:
#'   - table: data.frame with y1,y2,y3,y4, empirical, model, residual
#'   - ssr: scalar sum of squared residuals
#' @examples
#' \dontrun{
#' df <- data.frame(y1=rbinom(200,1,.3), y2=rbinom(200,1,.3),
#'                  y3=rbinom(200,1,.3), y4=rbinom(200,1,.3), weight=rep(1,200))
#' params <- list(theta0=0.10, theta01=0.15, theta1=0.08, theta10=0.12, pi=0.05)
#' gof <- goodness_of_fit_ar2(df, params)
#' gof$ssr
#' }
#' @export
goodness_of_fit_ar2 <- function(df, params) {
  # Empirical cell proportions — single pass using integer cell index
  total_weight <- sum(df$weight)

  # Encode each row as an integer cell ID in 1..16
  cell_id  <- df$y1 + 2L * df$y2 + 4L * df$y3 + 8L * df$y4 + 1L
  wt_sums  <- tapply(df$weight, cell_id, sum)
  # tapply returns a named vector; ensure all 16 cells are represented
  all_ids  <- 1:16
  empirical_vec <- as.numeric(wt_sums[match(all_ids, as.integer(names(wt_sums)))])
  empirical_vec[is.na(empirical_vec)] <- 0
  empirical_vec <- empirical_vec / total_weight

  # Re-order to match expand.grid(y1=0:1, y2=0:1, y3=0:1, y4=0:1)
  cell_keys <- as.matrix(expand.grid(y1=0:1, y2=0:1, y3=0:1, y4=0:1))
  key_ids   <- cell_keys[, "y1"] + 2L * cell_keys[, "y2"] +
               4L * cell_keys[, "y3"] + 8L * cell_keys[, "y4"] + 1L
  empirical <- empirical_vec[key_ids]

  # Model cell probabilities
  model_df <- model_cell_probs_ar2(params)

  # Join by y1,y2,y3,y4 ordering (expand.grid has same ordering)
  model_prob <- model_df$model_prob

  residual <- model_prob - empirical

  gof_table <- data.frame(
    cell_keys,
    empirical  = round(empirical,  6),
    model_prob = round(model_prob, 6),
    residual   = round(residual,   6)
  )

  ssr <- sum(residual^2)

  list(table = gof_table, ssr = ssr)
}
