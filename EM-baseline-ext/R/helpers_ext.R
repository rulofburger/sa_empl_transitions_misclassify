# ==============================================================================
# EM-baseline-ext: Shared helpers used across extension E-steps
# Created: 2026-05-06
#
# Defines internal helpers that are called by more than one E-step file:
#   .hmat_cache()              — memoised latent_histories() (constant 8×3)
#   .log_markov_trans_indiv()  — individual-specific N×8 log-transition matrix
#   .log_misclass_wave_ext()   — N×8 log-emission for symmetric/none models
#
# NOTE: This file MUST be sourced before estep_covariates.R, estep_fmm.R, and
# estep_inconsistency.R. source_all.R and helper-source.R enforce this order.
# ==============================================================================

# ---- Memoised latent history matrix -----------------------------------------

# latent_histories() returns the same constant 8×3 matrix on every call.
# Cache it in a private environment to avoid reconstructing it 1000×/EM run.
.hmat_env <- new.env(parent = emptyenv())
.hmat_env$.cache <- NULL

#' Return cached 8×3 latent history matrix
#'
#' Calls \code{\link{latent_histories}} on the first invocation and caches the
#' result. Subsequent calls return the cached matrix directly.
#'
#' @return 8×3 integer matrix; rows are histories, columns are waves 1–3.
#' @keywords internal
.hmat_cache <- function() {
  if (is.null(.hmat_env$.cache))
    .hmat_env$.cache <- latent_histories()
  .hmat_env$.cache
}


# ---- Log Markov transition probabilities for individual-specific rates -------

#' Log Markov transition probabilities for individual-specific rates
#'
#' Vectorised computation: precomputes four log-probability N-vectors and fills
#' matrix columns by logical index, avoiding the serial for-loop over H=8
#' histories.
#'
#' @param h_from Integer vector of length H (8): latent state at wave t-1.
#' @param h_to   Integer vector of length H (8): latent state at wave t.
#' @param theta1_i N-vector: employment persistence for each individual.
#' @param theta0_i N-vector: job-finding rate for each individual.
#' @return N × H matrix of log transition probabilities.
#' @keywords internal
.log_markov_trans_indiv <- function(h_from, h_to, theta1_i, theta0_i) {
  log_t1  <- log(theta1_i)       # stay employed
  log_1t1 <- log(1 - theta1_i)  # separate
  log_t0  <- log(theta0_i)       # find job
  log_1t0 <- log(1 - theta0_i)  # stay nonemployed

  N   <- length(theta1_i)
  H   <- length(h_from)
  out <- matrix(0, nrow = N, ncol = H)

  # Fill by logical index — vectorised, no per-column loop
  if (any(h_from == 1L & h_to == 1L)) out[, h_from == 1L & h_to == 1L] <- log_t1
  if (any(h_from == 1L & h_to == 0L)) out[, h_from == 1L & h_to == 0L] <- log_1t1
  if (any(h_from == 0L & h_to == 1L)) out[, h_from == 0L & h_to == 1L] <- log_t0
  if (any(h_from == 0L & h_to == 0L)) out[, h_from == 0L & h_to == 0L] <- log_1t0
  out
}


# ---- Log misclassification emission for extensions --------------------------

#' Log misclassification emission for extension E-steps (symmetric or none)
#'
#' Variant of the baseline \code{.log_misclass_wave} restricted to
#' \code{"symmetric"} and \code{"none"} model types for use in extension
#' E-steps.  Uses logical-index assignment instead of \code{ifelse()} to avoid
#' evaluating both branches and allocating redundant intermediate vectors.
#'
#' @param s_t Integer N-vector: observed state at wave t.
#' @param h_t Integer H-vector: latent state for all histories at t.
#' @param model_type Character: \code{"symmetric"} or \code{"none"}.
#' @param pi Scalar misclassification probability (used when \code{model_type = "symmetric"}).
#' @return N × H matrix of log emission probabilities.
#' @keywords internal
.log_misclass_wave_ext <- function(s_t, h_t, model_type, pi = 0) {
  N          <- length(s_t)
  H          <- length(h_t)
  match_mask <- outer(s_t, h_t, "==")  # N × H logical

  if (model_type == "none") {
    out             <- matrix(-Inf, N, H)
    out[match_mask] <- 0
    return(out)
  }

  # Symmetric misclassification: P(s=h) = 1-pi, P(s!=h) = pi
  pi              <- max(min(pi, 1 - 1e-12), 1e-12)
  out             <- matrix(log(pi), N, H)         # non-match
  out[match_mask] <- log(1 - pi)                   # match
  out
}
