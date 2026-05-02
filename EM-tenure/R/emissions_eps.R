# ==============================================================================
# EM-tenure: Emission log-densities for the eps (Spec I) model
# ==============================================================================
# Created: 2026-04-30
# Updated: 2026-05-01 — replaced all-or-nothing (AON) approximation with full
#   per-wave Bernoulli contamination enumeration (2^K patterns for K <= 3).
# Updated: 2026-05-01 — fixed K=1 handler for t_1 > a (shifted Exp mixture).
#
# Implements the spell-pair joint tenure emission of Spec I:
#   - Spell length T_g ~ Exp(lambda_g) (no sigma).
#   - Each observed tenure report is independently clock-consistent
#     (probability 1 - eps) or contaminated (probability eps, iid Exp draw).
#   - For a maximal latent E-spell with K observed tenures, the joint emission
#     is the sum over all 2^K contamination patterns m in {0,1}^K:
#
#       p(m) = eps^|X| * (1-eps)^|C|
#       T_g pinned by first clean wave j*: T_g = g_{j*} - d_{j*}*Delta
#       All other clean waves must satisfy g_j = T_g + d_j*Delta (within tol).
#       Log density = p(m) * lambda_g * exp(-lambda_g * T_g)   [if any clean]
#                               * prod_{j in X} lambda_g * exp(-lambda_g * g_j)
#       If all contaminated (|C|=0): density = prod_j lambda_g*exp(-lambda_g*g_j)
#       If T_g <= 0 (impossible spell start), pattern density = 0.
#
#   - For K = 0 spells: contributes nothing to the tenure emission.
#   - For K = 1 spells with offset = 0 (t_1 = a): plain Exp; eps drops out.
#   - For K = 1 spells with offset > 0 (t_1 > a): 2-pattern mixture;
#     clean branch is shifted Exp, contaminated branch is unshifted Exp.
#     eps is identifiable; tau_sum = P(contaminated | g).
#   - For K = 2: full 4-pattern enumeration (CC, CX, XC, XX).
#   - For K = 3: full 8-pattern enumeration (CCC, ..., XXX).
#
# Design rationale: The original AON approximation ("all clean OR all contam")
# is inconsistent with the per-wave Bernoulli M-step used in mstep_eps.R.
# Full enumeration makes the emission and M-step consistent (same Q function),
# restoring EM monotonicity. See brainstorm 2026-05-01-eps-model-convergence-fix.md.
#
# Return value change from original: 'tau' is now 'tau_sum' with different
# semantics — it is E[# contaminated waves | g], not P(at least one contam).
#
# Companion: documents/EM tenure epsilon.tex (Section 3.4-3.5).
# ==============================================================================

#' Log Exp(lambda_g) emission density (vectorised)
#'
#' Returns log f_Exp(g; lambda_g) = log(lambda_g) - lambda_g * g for g > 0.
#' Used for K = 1 spell emissions and for miscl-employed singletons.
#'
#' @param g Numeric vector of observed tenures (years), > 0.
#' @param lambda_g Exponential rate (scalar, > 0).
#' @return Numeric vector of log-densities. -Inf for g <= 0 or non-finite g.
#' @details Used for K = 1 spell emissions and miscl-employed singleton waves.
#'   Returns -Inf for any non-positive or non-finite g (avoids NaN in downstream
#'   log-sum-exp operations).
#' @keywords internal
.log_exp_g <- function(g, lambda_g) {
  out <- log(lambda_g) - lambda_g * g
  out[!is.finite(g) | g <= 0] <- -Inf
  out
}


#' Spell-pair joint tenure emission for a single E-spell (vectorised)
#'
#' For a maximal latent E-spell with relative wave offsets \code{t_offsets}
#' and per-individual tenure observations \code{g_mat} and observation
#' indicators \code{s_mat}, computes the joint log-emission and the
#' auxiliary quantities needed for the M-step.
#'
#' Algorithm (per row \eqn{i}):
#' \enumerate{
#'   \item \eqn{K_i = \sum_t s_{it}} (number of observed tenures in spell).
#'   \item If \eqn{K_i = 0}: log-emission = 0; no contribution to suff stats.
#'   \item If \eqn{K_i = 1} and offset = 0 (\eqn{t_1 = a}): both patterns
#'     give the same density so \eqn{\varepsilon} drops out.
#'     log-emission = \eqn{\log\lambda_g - \lambda_g g}. tau_sum = 0.
#'   \item If \eqn{K_i = 1} and offset > 0 (\eqn{t_1 > a}): clean branch
#'     gives shifted Exp, contaminated gives unshifted. 2-pattern mixture;
#'     tau_sum = P(contaminated | g) > 0.
#'   \item If \eqn{K_i \ge 2}: full \eqn{2^{K_i}} contamination-pattern
#'     enumeration (4 patterns for K=2; 8 for K=3). For each pattern
#'     \eqn{m = (m_1, \ldots, m_K) \in \{0,1\}^K}:
#'     \itemize{
#'       \item eps weight: \eqn{\varepsilon^{|X|} (1-\varepsilon)^{|C|}}.
#'       \item \eqn{T_g} pinned by first clean wave \eqn{j^*}:
#'         \eqn{T_g = g_{j^*} - d_{j^*}\Delta}. If \eqn{T_g \le 0},
#'         pattern density = 0.
#'       \item All subsequent clean waves must satisfy
#'         \eqn{g_j = T_g + d_j\Delta} within \code{tol}. If not, density = 0.
#'       \item Log-density = \eqn{\log\lambda_g - \lambda_g T_g}
#'         (if any clean) + \eqn{\sum_{j \in X}(\log\lambda_g - \lambda_g g_j)}.
#'       \item If all contaminated: \eqn{T_g} is integrated out;
#'         density = \eqn{\prod_j \lambda_g e^{-\lambda_g g_j}}.
#'     }
#'     The log-emission is log-sum-exp over all \eqn{2^K} patterns.
#' }
#'
#' Auxiliary outputs (for the M-step):
#' \itemize{
#'   \item \code{tau_sum}: \eqn{E[\text{# contaminated waves} | g]}, in
#'         \eqn{[0, K]}. Always 0 for \eqn{K \le 1}.
#'   \item \code{lambda_count}: \eqn{E[\text{# Exp evaluations} | g]}.
#'         For K=1: 1. For K>=2: 1 (T_g) + E[# contaminated waves].
#'   \item \code{lambda_xsum}: \eqn{E[\text{sum of Exp arguments} | g]}.
#'         Includes T_g (if any clean wave) plus all contaminated g_j.
#' }
#'
#' @param g_mat Numeric \code{N x L} matrix of tenures at each spell wave.
#'   Values at unobserved waves (where \code{s_mat} is FALSE) are not used.
#' @param s_mat Logical \code{N x L} matrix; TRUE where tenure is observed.
#' @param t_offsets Integer vector of length \code{L}: relative wave offsets
#'   within the spell (e.g., \code{c(0, 1, 2)} for a spell spanning waves
#'   \code{a, a+1, a+2}; origin is spell-start wave).
#' @param lambda_g Exponential rate (> 0).
#' @param eps Per-wave contamination probability in (0, 1).
#' @param tol Tolerance for clock-consistency check (default 1e-6 years).
#' @return Named list with N-vectors: \code{loglik}, \code{tau_sum}
#'   (0 for K=0 or K=1/offset=0; positive for K=1/offset>0 and K>=2),
#'   \code{eps_informative} (logical; TRUE when spell contributes to eps
#'   sufficient statistics: K>=2 or K=1 with offset>0),
#'   \code{K}, \code{lambda_count}, \code{lambda_xsum}.
#' @references TeX: \emph{EM tenure epsilon.tex}, Section 3.5.
#' @examples
#' \dontrun{
#' g_mat <- matrix(c(2.0, 2.25), nrow = 1)
#' s_mat <- matrix(c(TRUE, TRUE), nrow = 1)
#' log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = 0.15, eps = 0.20)
#' }
#' @export
log_emission_spell_g <- function(g_mat, s_mat, t_offsets,
                                 lambda_g, eps, tol = 1e-6) {
  if (!is.matrix(g_mat) || !is.matrix(s_mat)) {
    stop("log_emission_spell_g: g_mat and s_mat must be matrices.")
  }
  if (!identical(dim(g_mat), dim(s_mat))) {
    stop("log_emission_spell_g: g_mat and s_mat must have identical dims.")
  }
  if (!is.logical(s_mat)) {
    stop("log_emission_spell_g: s_mat must be a logical matrix (not numeric/integer).")
  }
  if (length(t_offsets) != ncol(s_mat)) {
    stop("log_emission_spell_g: t_offsets must have length ncol(s_mat).")
  }
  if (!is.integer(t_offsets) || any(is.na(t_offsets)) || any(t_offsets < 0L)) {
    stop("log_emission_spell_g: t_offsets must be a non-negative integer vector (use 0L, 1L, ...).")
  }
  if (!is.finite(tol) || tol <= 0) {
    stop(sprintf("log_emission_spell_g: tol must be > 0; got %.4g", tol))
  }
  if (!is.finite(lambda_g) || lambda_g <= 0) {
    stop(sprintf("log_emission_spell_g: lambda_g must be > 0; got %.4g", lambda_g))
  }
  if (!is.finite(eps) || eps <= 0 || eps >= 1) {
    stop(sprintf("log_emission_spell_g: eps must be in (0, 1); got %.4g", eps))
  }
  if (any(!is.finite(g_mat[s_mat]))) {
    stop("log_emission_spell_g: NA/Inf tenure at an observed (s=TRUE) position.")
  }

  N <- nrow(s_mat)
  L <- ncol(s_mat)
  K <- as.integer(rowSums(s_mat))

  if (any(K > 3L)) {
    stop(sprintf(
      "log_emission_spell_g: K > 3 not supported (%d rows affected; max K = %d).",
      sum(K > 3L), max(K)
    ))
  }

  loglik          <- numeric(N)
  tau_sum         <- numeric(N)    # E[# contaminated waves | g]; 0 for K=0 or K=1/offset=0
  eps_informative <- logical(N)    # TRUE iff spell contributes to eps sufficient stats
  lambda_count    <- numeric(N)
  lambda_xsum     <- numeric(N)

  if (N == 0L) {
    return(list(loglik = loglik, tau_sum = tau_sum, K = K,
                eps_informative = eps_informative,
                lambda_count = lambda_count,
                lambda_xsum = lambda_xsum))
  }

  Delta   <- .QUARTER_YEARS
  if (!is.finite(Delta) || Delta <= 0) {
    stop(sprintf("log_emission_spell_g: .QUARTER_YEARS is not a positive finite value; got %.4g", Delta))
  }
  log_lam <- log(lambda_g)
  log_1me <- log1p(-eps)   # log(1 - eps)
  log_e   <- log(eps)

  # --- K = 1 -----------------------------------------------------------------
  # When t_1 = a (offset == 0): both contamination patterns yield
  #   lambda_g * exp(-lambda_g * g_obs), so eps drops out. Plain Exp eval.
  # When t_1 > a (offset > 0): the clean branch gives a *shifted* Exp density
  #   lambda_g * exp(-lambda_g * (g - d*Delta)), while the contaminated branch
  #   gives the unshifted lambda_g * exp(-lambda_g * g). They differ, so eps
  #   is identifiable and a 2-pattern mixture is required.
  # See EM tenure epsilon.tex, Section 3.4 (Singletons paragraph).
  mask1 <- K == 1L
  if (any(mask1)) {
    s_sub1  <- s_mat[mask1, , drop = FALSE]
    g_sub1  <- g_mat[mask1, , drop = FALSE]
    n1      <- sum(mask1)
    col_obs <- max.col(s_sub1, ties.method = "first")
    g_obs   <- g_sub1[cbind(seq_len(n1), col_obs)]
    d_obs   <- t_offsets[col_obs]   # relative wave offset from spell start

    # Split: at-start (d=0) vs shifted (d>0)
    at_start <- d_obs == 0L

    # Precompute indices to avoid double-subscript intermediate copies
    idx1   <- which(mask1)
    idx_at <- idx1[at_start]
    idx_sh <- idx1[!at_start]

    # --- K=1, t_1 = a (d=0): eps drops out, plain Exp ---
    if (any(at_start)) {
      loglik[idx_at]       <- log_lam - lambda_g * g_obs[at_start]
      lambda_count[idx_at] <- 1L
      lambda_xsum[idx_at]  <- g_obs[at_start]
      # tau_sum and eps_informative remain FALSE/0: no eps info from offset-0 singletons
    }

    # --- K=1, t_1 > a (d>0): 2-pattern mixture (clean vs contaminated) ---
    if (any(!at_start)) {
      g1s  <- g_obs[!at_start]
      d1s  <- d_obs[!at_start]
      T_g  <- g1s - d1s * Delta        # T_g implied by clean branch
      v    <- T_g > 0                  # T_g must be positive; invalid → clean branch = -Inf

      # log-densities for the two patterns
      lp_clean  <- rep(-Inf, sum(!at_start))
      if (any(v)) lp_clean[v]  <- log_1me + log_lam - lambda_g * T_g[v]
      lp_contam <- log_e + log_lam - lambda_g * g1s   # always finite (g > 0)

      # log-sum-exp (lp_contam is always finite, so lp_mx is always finite)
      lp_mx   <- pmax(lp_clean, lp_contam)
      ll_k1   <- lp_mx + log(exp(lp_clean - lp_mx) + exp(lp_contam - lp_mx))

      nu       <- exp(lp_contam - ll_k1)  # P(contaminated | g): weight on contam branch
      T_g_safe <- pmax(T_g, 0)            # guard 0 * (-Inf) NaN when T_g <= 0 (nu=1 there)

      loglik[idx_sh]          <- ll_k1
      tau_sum[idx_sh]         <- nu
      eps_informative[idx_sh] <- TRUE     # K=1/offset>0: eps is identifiable
      lambda_count[idx_sh]    <- 1L
      lambda_xsum[idx_sh]     <- nu * g1s + (1 - nu) * T_g_safe
    }
  }

  # --- K = 2 -----------------------------------------------------------------
  # Full 4-pattern enumeration: CC, CX, XC, XX
  mask2 <- K == 2L
  if (any(mask2)) {
    g_sub2 <- g_mat[mask2, , drop = FALSE]
    s_sub2 <- s_mat[mask2, , drop = FALSE]
    n2     <- nrow(g_sub2)
    row2   <- seq_len(n2)

    # Identify the two observed columns per row.
    j1 <- max.col(s_sub2, ties.method = "first")   # first observed column
    j2 <- max.col(s_sub2, ties.method = "last")    # second observed column

    d1 <- t_offsets[j1]  # relative offset at first observed wave
    d2 <- t_offsets[j2]  # relative offset at second observed wave

    g1 <- g_sub2[cbind(row2, j1)]
    g2 <- g_sub2[cbind(row2, j2)]

    # T_g implied by each wave being clean
    T_from_j1 <- g1 - d1 * Delta   # valid only when > 0
    T_from_j2 <- g2 - d2 * Delta

    v1 <- T_from_j1 > 0  # CC and CX patterns require T_g from j1 > 0
    v2 <- T_from_j2 > 0  # XC pattern requires T_g from j2 > 0

    # Full 4-pattern enumeration: CC, CX, XC, XX. Each tenure observation is
    # independently clean (C, prob 1-eps) or contaminated (X, iid Exp, prob eps).
    # Pre-fill with -Inf; assign only valid rows to avoid ifelse() evaluating
    # both branches (saves ~2 n2-length temporary vectors per pattern).

    # CC: both clean; consistency g2 = g1 + (d2-d1)*Delta
    consistent_k2 <- abs(g2 - g1 - (d2 - d1) * Delta) < tol
    lp_cc <- rep(-Inf, n2)
    m_cc  <- v1 & consistent_k2
    if (any(m_cc)) lp_cc[m_cc] <- 2 * log_1me + log_lam - lambda_g * T_from_j1[m_cc]
    # CX: j1 clean, j2 contaminated
    lp_cx <- rep(-Inf, n2)
    if (any(v1)) lp_cx[v1] <- log_1me + log_e + 2 * log_lam - lambda_g * (T_from_j1[v1] + g2[v1])
    # XC: j1 contaminated, j2 clean
    lp_xc <- rep(-Inf, n2)
    if (any(v2)) lp_xc[v2] <- log_e + log_1me + 2 * log_lam - lambda_g * (g1[v2] + T_from_j2[v2])
    # XX: both contaminated (T_g integrated out) — always finite
    lp_xx <- 2 * log_e + 2 * log_lam - lambda_g * (g1 + g2)

    # Log-sum-exp over 4 patterns (XX is always finite, so mx2 is always finite)
    mx2 <- pmax(lp_cc, lp_cx, lp_xc, lp_xx)
    ll2 <- mx2 + log(exp(lp_cc - mx2) + exp(lp_cx - mx2) +
                     exp(lp_xc - mx2) + exp(lp_xx - mx2))

    w_cc <- exp(lp_cc - ll2)
    w_cx <- exp(lp_cx - ll2)
    w_xc <- exp(lp_xc - ll2)
    w_xx <- exp(lp_xx - ll2)

    # tau_sum = E[# contaminated waves]: CC=0, CX=1, XC=1, XX=2
    tau2 <- w_cx + w_xc + 2 * w_xx

    # lambda_count = E[# Exp evaluations]: CC=1, CX=2, XC=2, XX=2
    lc2 <- w_cc + 2 * (w_cx + w_xc + w_xx)

    # lambda_xsum = E[sum of Exp args].
    # Guard against 0 * (-Inf) NaN when T <= 0 (weight is 0 there anyway).
    T1s <- pmax(T_from_j1, 0)   # safe T_g from j1 (0 when invalid, w=0)
    T2s <- pmax(T_from_j2, 0)   # safe T_g from j2 (0 when invalid, w=0)
    lx2 <- T1s         * w_cc +
           (T1s + g2)  * w_cx +
           (g1  + T2s) * w_xc +
           (g1  + g2)  * w_xx

    loglik[mask2]           <- ll2
    tau_sum[mask2]          <- tau2
    eps_informative[mask2]  <- TRUE   # K=2: always eps-informative
    lambda_count[mask2]     <- lc2
    lambda_xsum[mask2]      <- lx2
  }

  # --- K = 3 -----------------------------------------------------------------
  # Full 8-pattern enumeration: CCC, CCX, CXC, CXX, XCC, XCX, XXC, XXX
  # In a 3-wave model K=3 implies L=3 and all spell columns observed.
  mask3 <- K == 3L
  if (any(mask3)) {
    g_sub3 <- g_mat[mask3, , drop = FALSE]
    n3     <- nrow(g_sub3)

    # All 3 columns are observed; extract each wave's tenure.
    g1 <- g_sub3[, 1L]
    g2 <- g_sub3[, 2L]
    g3 <- g_sub3[, 3L]

    d1 <- t_offsets[1L]; d2 <- t_offsets[2L]; d3 <- t_offsets[3L]

    # T_g implied by each wave being the first (and determining) clean wave
    T1 <- g1 - d1 * Delta   # from wave 1  (always > 0 since g1 > 0 and d1 >= 0)
    T2 <- g2 - d2 * Delta   # from wave 2
    T3 <- g3 - d3 * Delta   # from wave 3

    v2 <- T2 > 0   # patterns requiring T_g from wave 2
    v3 <- T3 > 0   # patterns requiring T_g from wave 3
    # v1 (T1 > 0) is practically always true since g1 > 0; handle via pmax below

    # Consistency indicators
    D21 <- (d2 - d1) * Delta   # expected g increment wave1 -> wave2
    D31 <- (d3 - d1) * Delta   # expected g increment wave1 -> wave3
    D32 <- (d3 - d2) * Delta   # expected g increment wave2 -> wave3

    cons12 <- abs(g2 - g1 - D21) < tol   # wave2 consistent from T_g at wave1
    cons13 <- abs(g3 - g1 - D31) < tol   # wave3 consistent from T_g at wave1
    cons23 <- abs(g3 - g2 - D32) < tol   # wave3 consistent from T_g at wave2

    # Pre-compute shared log-weight offsets
    l1 <- log_1me; l2 <- 2 * log_1me; l3 <- 3 * log_1me
    e1 <- log_e;   e2 <- 2 * log_e;   e3 <- 3 * log_e
    ll3 <- 3 * log_lam   # 3 * log(lambda_g)
    ll2 <- 2 * log_lam
    ll1 <- log_lam
    sum123 <- g1 + g2 + g3

    # Full 8-pattern enumeration packed into an n3×8 matrix to reduce
    # allocations from 16×n3 (8 pattern vecs + 8 weight vecs) to 3×n3
    # (lp_mat, w_mat, one rowSums buffer). [P2-11]
    lp_mat <- matrix(-Inf, n3, 8L)
    # CCC=col1, CCX=col2, CXC=col3, CXX=col4, XCC=col5, XCX=col6, XXC=col7, XXX=col8
    m_CCC <- cons12 & cons13
    if (any(m_CCC)) lp_mat[m_CCC, 1L] <- l3 + ll1 - lambda_g * T1[m_CCC]
    if (any(cons12)) lp_mat[cons12, 2L] <- l2 + e1 + ll2 - lambda_g * (T1[cons12] + g3[cons12])
    if (any(cons13)) lp_mat[cons13, 3L] <- l2 + e1 + ll2 - lambda_g * (T1[cons13] + g2[cons13])
    m_CXX <- T1 > 0
    if (any(m_CXX)) lp_mat[m_CXX, 4L] <- l1 + e2 + ll3 - lambda_g * (T1[m_CXX] + g2[m_CXX] + g3[m_CXX])
    m_XCC <- v2 & cons23
    if (any(m_XCC)) lp_mat[m_XCC, 5L] <- e1 + l2 + ll2 - lambda_g * (g1[m_XCC] + T2[m_XCC])
    if (any(v2)) lp_mat[v2, 6L] <- e2 + l1 + ll3 - lambda_g * (g1[v2] + T2[v2] + g3[v2])
    if (any(v3)) lp_mat[v3, 7L] <- e2 + l1 + ll3 - lambda_g * (g1[v3] + g2[v3] + T3[v3])
    lp_mat[, 8L] <- e3 + ll3 - lambda_g * sum123  # XXX: always finite

    # Log-sum-exp + normalise in one pass (XXX col finite ⇒ mx3 always finite)
    mx3   <- pmax(lp_mat[,1], lp_mat[,2], lp_mat[,3], lp_mat[,4],
                  lp_mat[,5], lp_mat[,6], lp_mat[,7], lp_mat[,8])
    w_mat <- exp(lp_mat - mx3)
    rs3   <- rowSums(w_mat)
    ll3_out <- mx3 + log(rs3)
    w_mat   <- w_mat / rs3   # normalise: w_mat[i,] now sums to 1

    # tau_sum = E[# contaminated waves]: col1=0, cols2-4,5,6,7 depend on # X
    # CCC=0, CCX/CXC/XCC=1, CXX/XCX/XXC=2, XXX=3
    tau3 <- (w_mat[,2] + w_mat[,3] + w_mat[,5]) +
             2 * (w_mat[,4] + w_mat[,6] + w_mat[,7]) +
             3 * w_mat[,8]

    # lambda_count = E[# Exp evaluations]: CCC=1, one-X patterns=2, rest=3
    lc3 <- w_mat[,1] + 2*(w_mat[,2]+w_mat[,3]+w_mat[,5]) +
            3*(w_mat[,4]+w_mat[,6]+w_mat[,7]+w_mat[,8])

    # lambda_xsum = E[sum of Exp args]; safe T values (0 when physically impossible)
    T1s <- pmax(T1, 0); T2s <- pmax(T2, 0); T3s <- pmax(T3, 0)
    lx3 <- T1s              * w_mat[,1] +
           (T1s + g3)       * w_mat[,2] +
           (T1s + g2)       * w_mat[,3] +
           (T1s + g2 + g3)  * w_mat[,4] +
           (g1  + T2s)      * w_mat[,5] +
           (g1  + T2s + g3) * w_mat[,6] +
           (g1  + g2  + T3s)* w_mat[,7] +
           sum123            * w_mat[,8]

    loglik[mask3]           <- ll3_out
    tau_sum[mask3]          <- tau3
    eps_informative[mask3]  <- TRUE   # K=3: always eps-informative
    lambda_count[mask3]     <- lc3
    lambda_xsum[mask3]      <- lx3
  }

  list(loglik = loglik, tau_sum = tau_sum, K = K,
       eps_informative = eps_informative,
       lambda_count = lambda_count, lambda_xsum = lambda_xsum)
}


#' Enumerate maximal E-spells in a 3-wave latent history
#'
#' Returns a list of integer vectors, one per maximal contiguous run of
#' \code{h_t = 1} in the history.
#'
#' Examples:
#' \itemize{
#'   \item \code{(0,0,0)} -> \code{list()}
#'   \item \code{(1,0,1)} -> \code{list(1L, 3L)}
#'   \item \code{(1,1,1)} -> \code{list(1:3)}
#'   \item \code{(0,1,1)} -> \code{list(2:3)}
#' }
#'
#' @param h Integer vector of length 3 with values in \{0, 1\}.
#' @return List of integer vectors of wave indices (each of length >= 1).
#' @keywords internal
.maximal_e_spells <- function(h) {
  if (length(h) != 3L) stop(".maximal_e_spells: h must have length 3.")
  out <- list()
  i <- 1L
  while (i <= 3L) {
    if (h[i] == 1L) {
      j <- i
      while (j + 1L <= 3L && h[j + 1L] == 1L) j <- j + 1L
      out[[length(out) + 1L]] <- i:j
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}
