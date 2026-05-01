---
date: 2026-04-30
title: "Spec I (epsilon model): design rationale, identification, and implementation plan"
status: implemented
chosen-approach: "Point-mass + Exp contamination with spell-pair joint emission"
tags: [em-algorithm, identification, tenure, contamination, story-A, spec-i]
predecessor-brainstorms:
  - 2026-04-28-mismatch-model-sigma-collapse-fix.md
  - 2026-04-29-duration-contamination-rho-model.md
---

# Spec I (epsilon model): design rationale and implementation plan

## TL;DR

The base EM-tenure spec collapses sigma_g^2 -> 0 and inflates pi to 26%
because its Gaussian measurement model is structurally wrong for tenure data.
The rho extension fixes pi (~4%) but routes ~33% of observations through rho,
which absorbs all tenure inconsistencies and so prevents tenure data from
identifying pi beyond what state-Markov violations alone provide. Spec I
replaces the Gaussian-EMG with a structurally correct point-mass +
Exponential-contamination mixture and adds a spell-pair joint emission so
that tenure flanks discriminate between (1,0,1) hypotheses where the same
spell continues vs. fresh job entry. This is the structure required to make
Story A in the paper -- "tenure data corroborates and refines the simple
model's pi estimate" -- empirically active in the model.

## Headline empirical evidence (from `diagnose_tenure_consistency_SA.R`)

| Diagnostic | N | Share |
|---|---|---|
| E->E continuation increments exactly +3mo | 343,850 | **63.2%** |
| (1,0,1) panels with tenure flanks corroborating miscl (g3 ~ g1+6mo) | 5,997 | **27.4%** |
| (1,0,1) panels with tenure flanks consistent with fresh start (g3 <= 3mo) | 5,997 | 18.3% |
| Within-panel "fresh starts" reporting tenure > 12 months (impossible) | 35,389 | **36.7%** |
| (0,1,0) panels with identical timegap_cat at flanks (corroborates miscl) | 6,445 | **63.7%** |

The 36.7% and 63.7% numbers are the descriptive workhorses for the paper:
- A "fresh start" with 12+ months of tenure is structurally impossible under
  genuine job entry. 36.7% of observed (0,1) transitions thus imply a lower
  bound on misclassification at the previous wave.
- Identical timegap categories at the flanks of a (0,1,0) panel are exactly
  what a misclassification-at-wave-2 generates: continuous nonemployment with
  the same elapsed nonemployment duration on both sides.

These are both descriptive (no model required) and complementary to the
structural Spec I estimate of pi.

## Story selection (from conversation)

The user wants the paper to lead with **Story A** (tenure data corroborates
and amplifies pi via flank consistency) but retain Story B (tenure data is a
robustness check that produces a similar pi to the simple model) for the
robustness section. Story A is the headline; the rho model is preserved in
the appendix as the Story B analogue.

## Why the base / rho models cannot tell Story A

In the base spec, sigma_g^2 -> 0 is the correct MLE because 63.2% of E->E
increments are exactly +3 months. With sigma_g^2 -> 0 the Gaussian increment
density becomes a Dirac spike at +3mo, which assigns zero likelihood to the
36.8% of non-standard increments. The model then routes those into
misclassified histories, inflating pi to 26%.

In the rho spec, the contamination probability rho absorbs the non-standard
tenure increments before they affect pi. By construction (under the spec's
emissions), the (1,1,1) hypothesis at s=(1,0,1) treats the wave-3 tenure as
an *independent* EMG draw -- discarding the structural information that
g_1 and g_3 belong to the same spell. So tenure flanks do not discriminate
hypotheses, and pi is identified from state-Markov violations alone.

## Spec I design

### Tenure measurement model

Replace the Gaussian increment + EMG level with a point-mass + Exp
contamination mixture:

- Spell length T_g ~ Exp(lambda_g) (no sigma).
- Each tenure report is either *clock-consistent* (probability 1-eps,
  reports the deterministic clock value exactly) or *contaminated*
  (probability eps, reports an independent Exp(lambda_g) draw).

### Spell-pair joint emission (the key innovation)

For each maximal latent E-spell [a, b] under history h with K observed
tenures (s_t = 1 at K of the b-a+1 spell waves):

```
ell_spell(g) = (1-eps)^K * Exp(g_a) * 1{g_t = g_a + (t-a)*Delta for all observed t}
             + (1 - (1-eps)^K) * prod_t Exp(g_t)
```

The first branch fires only when all K reports are exactly clock-consistent;
the second branch is approximate (treats all K reports as independent Exp).

For the critical case h=(1,1,1) under s=(1,0,1): K=2 (waves 1 and 3 observed,
same spell). Spell-pair emission has a point mass at {g_3 = g_1 + 0.5}, which
fires for the 27.4% of observations that exhibit this pattern, pushing
gamma toward (1,1,1) and increasing pi.

For h=(1,0,1) under s=(1,0,1): K=1 for each of the two singleton E-spells.
No point mass; standard Exp evaluation per wave. Loses the corroboration
factor.

### Parameters

`(alpha, theta_0, theta_1, pi, eps, lambda_g, lambda_d)` -- 7 parameters,
same as the base model.

### Identification (block separation)

1. **Markov prior** (alpha, theta_0, theta_1): from latent transition counts.
2. **State misclassification** (pi): from state-Markov violations *weighted
   by tenure-flank consistency through gamma*. This is the corroboration
   channel.
3. **Tenure contamination** (eps): from the empirical fraction of
   spell-pairs satisfying the exact clock relation. Method-of-moments
   estimates: ~0.37 from adjacent E->E pairs, ~0.48 from (1,0,1) bridges.
4. **Spell-length rates** (lambda_g, lambda_d): from inverse means of
   tenure / interval-censored timegap.

### Predicted empirical results

- **eps**: ~0.37 to 0.50 depending on which observations dominate the M-step.
- **pi**: 3% to 5% -- modestly above the simple-model 2.98%, reflecting
  tenure-flank corroboration.
- **theta_1**: ~0.97 (matching simple model and rho model).
- **theta_0**: ~0.03.
- **lambda_g**: ~0.17 (mean spell ~6 years).
- **lambda_d**: ~0.13 (mean nonemp spell ~7-8 years).
- **No sigma_g**: removed entirely.

If realised, these support the headline:

> Two independent identification channels -- state-Markov violations alone
> (simple model, pi=2.98%) and state-Markov violations + tenure flank
> consistency (Spec I, pi=4%) -- both produce a structural misclassification
> rate around 3-5%, supporting the paper's central claim that observed
> employment churn is partly a measurement artefact.

## Answers to user's six points

1. **Duration floor**: keep at one week (`0.25/12 = 0.0208 years`). The
   Spec I model treats the floor as the imputed value for unobserved-tenure
   waves; the spell-pair structure is unaffected because the floor only
   applies to nonemployment waves where tenure does not enter the emission
   anyway. (If desired, change to `1/12` (one month) for cleanliness; this
   changes nothing structural.)

2. **Tenure granularity**: confirmed in the data. 100% integer-month;
   only 6.2% multiples of 12 (year-rounding is a minor not dominant feature).
   The "63.2% exactly +3mo and ~7% small errors" pattern is captured by the
   point-mass + Exp mixture exactly.

3. **(1,0,1) diagnostic script**: delivered as
   `scripts/diagnose_tenure_consistency_SA.R`. Outputs five tables and
   three figures characterising the corroboration evidence. **Headline
   finding: 27.4% of (1,0,1) panels exhibit tenure flanks within 2 months
   of `g_1 + 6mo` -- the strongest possible evidence of miscl-at-wave-2 in
   the data.** Also: 36.7% of fresh starts have tenure > 12mo (impossible
   under genuine entry); 63.7% of (0,1,0) panels have identical timegap
   categories at flanks.

4. **Critical distinction (pi vs rho roles)**: in the unified mismatch
   spec the user proposes,
   - pi handles `s=(1,0,1)` with `g=(10mo, [n/a], 16mo)` -- "(1,1,1) miscl
     at wave 2 with consistent flanks". Spec I delivers this through the
     spell-pair point mass.
   - rho handles `s=(1,1,1)` with `g=(10, 72, 16)` -- "wave-2 tenure is a
     random draw, person was employed throughout". Spec I delivers this
     through the eps contamination probability on the wave-2 increment.
   Spec I unifies both cases through eps with appropriate weights, without
   needing rho. (A pure rho model misses the corroboration channel; a Spec II
   with both eps and rho is documented but not implemented in this round --
   the marginal value of rho on top of eps is small for QLFS.)

5. **Drop sigma in EMG -> use Exp**: yes. In the conversation we
   established that
   - The 63.2% point-mass at +3mo *is* the data-generating clock; not noise.
   - The 36.8% non-standard residual *is* contamination, not Gaussian noise.
   - Removing sigma simplifies identification (lambda_g is identified
     immediately from the inverse-mean of tenure under matched-employed) and
     removes the sigma -> 0 collapse pathology.

6. **Will pi be informed by tenure under Spec I?**: yes, but **only because
   of the spell-pair joint emission**. Without spell-pair emission, pi is
   identified from state-Markov violations alone (same as the simple model).
   The spell-pair emission causes responsibilities gamma_(1,1,1) to be
   pushed up (or down) for individual (1,0,1) observations based on whether
   their tenure flanks satisfy `g_3 = g_1 + 0.5`. This re-routing flows
   through the M-step formula `pi_hat = M / (3W)` to refine pi.

   **Critical implementation detail**: under the *base / rho* spec the
   wave-3 emission for h=(1,1,1) at s=(1,0,1) is `f_EMG(g_3)` with no
   reference to g_1 -- so tenure flanks do not enter pi identification.
   Spec I's spell-pair emission is the structural fix that makes this work.

## Implementation plan

The Spec I implementation lives alongside the base and rho code, not
replacing it. New files / functions:

```
EM-tenure/R/
  emissions_eps.R      (DONE)  log_emission_spell_g(); spell-pair point mass
  estep_eps.R          (DONE)  computes gamma + tau_spell
  mstep_eps.R          (DONE)  closed-form eps update + Exp-MLE for lambda
  em_driver_eps.R      (DONE)  em_fit_tenure_eps(); top-level driver
  init_params_eps.R    (DONE)  init_params_eps() -- delegates to init_params()
  source_all.R         (DONE)  new sources added

scripts/
  diagnose_tenure_consistency_SA.R  (DONE)

documents/
  EM tenure epsilon.tex (DONE)

EM-tenure/feedback/
  2026-04-30-epsilon-spec-design.md (THIS FILE)

EM-tenure/estimate_pipeline.R (DONE)  eps fit calls appended below rho block
```

120 eps-model tests pass; 31 rho-model tests pass (no regression). Two
thorough review rounds completed (2026-04-30 to 2026-05-01).

### Naming conventions

- Driver: `em_fit_tenure_eps(df, params0, stationary, linked, verbose)`.
- Init: `init_params_eps(df, linked = FALSE)`.
- Saved fits: `fit_eps_<run_id>.rds`,
  `fit_eps_stationary_<run_id>.rds`,
  `fit_eps_linked_<run_id>.rds`,
  `fit_eps_stationary_linked_<run_id>.rds`.
- run_summary.csv: 4 new rows per run, model in
  `{eps_free, eps_stationary, eps_linked, eps_stationary_linked}`.

### High-level pseudocode for the spell-pair emission

```r
log_emission_spell_g <- function(g, s_obs_mask, K, lambda_g, eps) {
  # g: numeric vector of tenures at all spell waves (length n_spell_waves)
  # s_obs_mask: logical vector, TRUE where tenure is observed (s = 1)
  # K = sum(s_obs_mask)
  if (K == 0L) return(0)            # no tenure observed in this spell
  if (K == 1L) {                     # singleton: just Exp evaluation
    g_obs <- g[s_obs_mask]
    return(log(lambda_g) - lambda_g * g_obs)
  }
  # K >= 2: mixture
  obs_idx <- which(s_obs_mask)
  g_obs   <- g[obs_idx]
  # Clock-consistency: all observed tenures equal g_first + (t - t_first) * Delta
  Delta    <- 0.25
  expected <- g_obs[1L] + (obs_idx - obs_idx[1L]) * Delta
  consistent <- all(abs(g_obs - expected) < 1e-8)

  log_p_clean   <- K * log1p(-eps) +
                   log(lambda_g) - lambda_g * g_obs[1L]    # only g_first contributes
  log_p_contam  <- log1p(-(1 - eps)^K) +
                   sum(log(lambda_g) - lambda_g * g_obs)

  if (consistent) {
    matrixStats::logSumExp(c(log_p_clean, log_p_contam))
  } else {
    log_p_contam   # clean branch has zero density off the clock
  }
}
```

### M-step for eps (closed-form)

```r
update_eps <- function(gamma, tau_spell_list) {
  # tau_spell_list[[i]][[h_idx]][[spell_idx]] = posterior P(contaminated | obs)
  # Sum responsibilities * tau across all (i, h, spell) where K >= 2.
  num <- 0; den <- 0
  for (i in seq_len(N)) {
    for (h_idx in seq_along(H)) {
      g <- gamma[i, h_idx]; if (g < 1e-12) next
      for (sp in tau_spell_list[[i]][[h_idx]]) {
        K <- sp$K
        if (K < 2L) next
        num <- num + w[i] * g * sp$tau          # contaminated weight
        den <- den + w[i] * g * 1               # total weight
      }
    }
  }
  num / den   # fraction of spell-pair observations that are contaminated
}
```

The closed form for eps via the K-th-root inversion in the .tex (eq 11)
is the formally correct version when spell sizes K vary across observations.
The simpler formula above is exact when K is constant within a homogeneous
sub-sample.

### Update to estimate_pipeline.R

The user has explicitly requested that the eps fit calls go *underneath*
the existing rho block, preserving the existing pipeline. Sketch
(append at end of `estimate_pipeline.R`):

```r
# ##############################################################################
# EPSILON-AUGMENTED ESTIMATION (Spec I: spell-pair point-mass + Exp contam.)
# Created: 2026-04-30
# ##############################################################################
custom_init_eps <- init_params_eps(df_qlfs, linked = FALSE)

message("\n=== Estimating eps model (free, non-stationary) ===")
fit_eps <- em_fit_tenure_eps(
  df = df_qlfs, params0 = ws$eps %||% custom_init_eps,
  stationary = FALSE, linked = FALSE, verbose = 2L
)

message("\n=== Estimating eps model (free, stationary) ===")
fit_eps_stationary <- em_fit_tenure_eps(
  df = df_qlfs, params0 = ws$eps_stationary %||% custom_init_eps,
  stationary = TRUE, linked = FALSE, verbose = 2L
)

message("\n=== Estimating eps model (linked, non-stationary) ===")
fit_eps_linked <- em_fit_tenure_eps(
  df = df_qlfs, params0 = ws$eps_linked %||% custom_init_eps,
  stationary = FALSE, linked = TRUE, verbose = 2L
)

message("\n=== Estimating eps model (linked, stationary) ===")
fit_eps_stationary_linked <- em_fit_tenure_eps(
  df = df_qlfs, params0 = ws$eps_stationary_linked %||% custom_init_eps,
  stationary = TRUE, linked = TRUE, verbose = 2L
)

saveRDS(fit_eps,                  file.path(results_dir, sprintf("fit_eps_%s.rds", run_id)))
saveRDS(fit_eps_stationary,       file.path(results_dir, sprintf("fit_eps_stationary_%s.rds", run_id)))
saveRDS(fit_eps_linked,           file.path(results_dir, sprintf("fit_eps_linked_%s.rds", run_id)))
saveRDS(fit_eps_stationary_linked,file.path(results_dir, sprintf("fit_eps_stationary_linked_%s.rds", run_id)))

# Append summary rows analogous to rho block, with eps in place of rho.
```

The `.make_summary_row_eps()` function is **not** separate in the final implementation.
The pipeline uses a unified `.make_summary_row()` that handles all three model families
(base, rho, eps) via `%||% NA_real_` for model-specific fields. The LR test
`fit_eps$loglik vs fit_em$loglik` follows the rho block pattern.

## Risks and open questions

1. **Spell-pair joint emission with K=3** (h=(1,1,1) under s=(1,1,1)):
   the formula treats the "at least one contaminated" branch as
   approximately independent. This is exact only when *all* K are
   contaminated. For K=3 with one contaminated the exact joint is more
   complex (involves a conditional Exp + delta over the clean pair). The
   approximation should be adequate for QLFS but a referee may push back.
   *Mitigation*: write up the exact joint as a remark and note the
   approximation is asymptotically negligible.

2. **eps_d for timegap**: not introduced. If a referee asks why timegap
   does not get an analogous error parameter, the answer is that timegap is
   already discretised into 7 categories which provides intrinsic robustness.
   *Mitigation*: add eps_d as a robustness check in the appendix if requested.

3. **Effect of duration floor (1 week vs 1 month)**: the floor only enters
   for unobserved-tenure waves under nonemployment hypotheses, where tenure
   is not used in the emission. Should not affect Spec I results.

4. **Within-panel start emission**: Spec I uses Exp(lambda_g) for fresh
   starts (no point mass). This is justified by the empirical finding that
   only 39% of "fresh starts" have tenure <= 3mo; the rest are likely
   miscl-at-t-1 recovery. The Exp marginal is the right population
   distribution. A refinement would model fresh starts as Uniform[0, 0.25]
   for the genuine-entry component, but this adds complexity.

## Suggested next steps

1. ~~(User confirms) the Spec I plan and the implementation skeleton.~~ ✓ done
2. ~~Implement the new R files following the skeleton in this document.~~ ✓ done (2026-04-30 to 2026-05-01)
3. Run on QLFS data; expect convergence in ~10-30 iterations from cold start.
4. Compare fit_eps vs fit_em vs fit_rho via run_summary.csv.
5. Generate a paper table comparing simple, EM-tenure, rho, and eps models
   on (theta_0, theta_1, pi, eps/rho, log-likelihood).
6. Write the spec into the paper main text (section: "Auxiliary identification
   from tenure consistency").

## Implementation design choices (2026-04-30, post-coding)

Three small deviations from the .tex spec were adopted during
implementation. They are noted here so the documentation, paper, and code
remain in sync.

### 1. eps M-step uses per-wave Bernoulli, not the K-th-root form

The .tex (eq. eps_mstep) writes
`eps_hat = 1 - (sum w*gamma*(1-tau) / sum w*gamma*K)^(1/Kbar)`,
which is the spell-level "at least one contaminated" Bernoulli MLE under
the assumption that K is constant. The implementation instead uses the
per-wave Bernoulli mixture-weight estimator

```
eps_hat = sum_{K>=2} w_i * gamma_ih * tau * K
        / sum_{K>=2} w_i * gamma_ih * K
        = Eps_num / Eps_den
```

which is the standard EM update under the same all-or-nothing
approximation already baked into the contamination branch of the spell
emission. The two estimators coincide when K is constant; they differ at
order K^{-1} when K varies (1, 2, 3 across spells). The per-wave form is
unambiguous and avoids a small upward bias in the K-th-root form when
spells of different lengths are pooled. The choice has no effect on the
EM fixed point at the true parameter and is a numerics-only refinement.

The .tex has been updated to present the per-wave form as primary and the
K-th-root form as a remark.

### 2. eps numerical bounds

`eps_floor = 1e-4` and `eps_cap = 0.95` (vs. the `pi_cap = 0.49` used for
pi and rho). Rationale:

- The floor avoids log(0) in the clean-branch log-density when eps is
  near zero on a particular E-step iteration (e.g., during initialisation
  or after a warm start from rho).
- The cap is permissive (0.95, not 0.49) because eps is not restricted
  to be small under the model. Empirical method-of-moments estimates are
  already in the 0.37-0.50 range; capping at 0.49 would be active and
  would distort the fixed point.

### 3. Timegap (lambda_d) emissions reused verbatim from base / rho

Spec I changes the tenure block only. The full timegap emission
machinery is imported from the base implementation:

- `log_emission_interval_d()`: marginal interval probability for fresh
  starts and miscl-as-nonemployed.
- `log_emission_transition_d()`: conditional transition probability for
  nonemployment continuations with both flanks observed.
- `log_emission_start_d_cat()`: start emission (cat == 1) -- not used
  in the current eps model because we do not introduce eps_d.

This preserves apples-to-apples comparability of lambda_d across the
base, rho, and eps specifications and concentrates Spec I's structural
change in the tenure block, exactly as the .tex states.

### 4. K>3 guard in log_emission_spell_g()

During the second review round (2026-05-01) an explicit guard was added:

```r
if (any(K > 3L)) stop("K > 3 not supported: only 3-wave panels are expected")
```

With 3 observed waves, K ∈ {0,1,2,3} is guaranteed by the data structure,
so this guard is defensive only. It makes the supported domain of the
function explicit and prevents silent wrong results if the function is
reused in a future 4-wave extension without updating the enumeration logic.
