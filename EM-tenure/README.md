# EM-tenure

EM algorithm for estimating employment transition probabilities with tenure and nonemployment durations, allowing for symmetric misclassification. This folder implements three model specifications; the **contamination model (epsilon)** is the primary specification of interest.

## Overview for reviewers

If you are reviewing this code alongside the paper, start here:

1. **Theory**: `documents/EM tenure epsilon.tex` — the self-contained companion paper describing the contamination model.
2. **Estimation**: `estimate_pipeline.R` — end-to-end script that ingests data and runs all model variants (base EM-tenure, epsilon). The epsilon section starts at the comment block `# EPSILON-AUGMENTED ESTIMATION`.
3. **Epsilon model code**: The five `*_eps.R` files in `R/` (see file structure below).
4. **Tests**: `tests/testthat/test-eps-model.R`.

## Why the contamination (epsilon) model?

The base EM-tenure model (`documents/EM tenure.tex`) uses Gaussian measurement error (σ²_g) for tenure and an EMG spell-level distribution. Three empirical facts in the QLFS data rule this out:

1. **Tenure is integer-month, not continuous Gaussian.** Reported tenures are discrete month values.
2. **63.2% of E→E continuation increments are exactly 0.25 years.** The Gaussian MLE for σ²_g collapses to zero. The remaining 36.8% include wildly off values (negatives, year-shifted) — a contamination distribution, not Gaussian noise.
3. **36.7% of within-panel new job entrants report tenure > 12 months.** A genuine job entry cannot have 12+ months of tenure — these must be misclassifications or tenure reporting errors.

The base EM-tenure model cannot accommodate these facts: σ²_g collapses to zero, and all inconsistency is routed through π (inflating it to over 25%). The contamination model replaces the Gaussian + EMG measurement model with a **point-mass + Exponential contamination mixture**, eliminating σ²_g entirely, and introduces a **spell-pair joint emission** that uses tenure flanks to identify misclassification.

## Contamination model (epsilon) — specification

### Parameters

| Parameter | Description |
|-----------|-------------|
| α | Initial employment probability |
| θ₁ | P(employed \| employed last period) |
| θ₀ | P(employed \| nonemployed last period) |
| π | Symmetric misclassification probability |
| ε | Per-wave tenure contamination probability |
| λ_g | Exponential rate for employment spells |
| λ_d | Exponential rate for nonemployment spells |

Note: σ²_g is **eliminated** — there is no Gaussian component. The parameter ε replaces it.

### Tenure measurement model

Each observed tenure report is independently either:
- **Clock-consistent** (probability 1 − ε): the observation equals the deterministic clock value exactly (point mass).
- **Contaminated** (probability ε): an independent draw from Exp(λ_g), unrelated to the clock.

For a maximal latent employment spell with K observed tenure waves, the joint tenure emission sums over all 2^K contamination patterns. This is the **spell-pair joint emission** — the key innovation over the base EM-tenure model. When a spell spans an intermediate misclassified wave (e.g., h = (1,1,1) but s = (1,0,1)), the tenure at waves 1 and 3 are "flanks" of the same spell. If both are clock-consistent, they must satisfy g₃ = g₁ + 2Δ. This structural constraint lets tenure directly inform misclassification identification.

### Nonemployment duration

Nonemployment duration D is observed as an integer category code 1–7 (interval-censored). The model treats the true D as Exp(λ_d). This component is shared with the base EM-tenure model.

| Code | Interval (months) | Interval (years) |
|:----:|:-----------------:|:----------------:|
| 1 | [0, 3) | [0, 0.25) |
| 2 | [3, 6) | [0.25, 0.50) |
| 3 | [6, 9) | [0.50, 0.75) |
| 4 | [9, 12) | [0.75, 1.00) |
| 5 | [12, 36) | [1.00, 3.00) |
| 6 | [36, 60) | [3.00, 5.00) |
| 7 | [60, ∞) | [5.00, ∞) |

### Model variants

Each specification is estimated in four variants via two orthogonal switches:

| Variant | Stationarity | CTMC link | Free parameters |
|---------|:---:|:---:|:---:|
| Free, non-stationary | — | — | 7 (α, θ₀, θ₁, π, ε, λ_g, λ_d) |
| Free, stationary | α = θ₀/(θ₀ + 1 − θ₁) | — | 6 |
| Linked, non-stationary | — | λ_g = −log(θ₁)/Δ, λ_d = −log(1−θ₀)/Δ | 5 |
| Linked, stationary | both | both | 4 (θ₀, θ₁, π, ε) |

The **linked + stationary** variant is the most parsimonious. `estimate_pipeline.R` runs all four and performs nested LR tests.

## File structure

```
EM-tenure/
├── R/
│   ├── utils.R              # Constants, .bound01, .logsumexp, timegap helpers
│   ├── transforms.R         # Logit, CTMC link functions
│   ├── latent_histories.R   # 8 latent histories, clocks, prior
│   │
│   │  ── Base EM-tenure model ──
│   ├── emissions.R          # EMG, Normal increment, Normal start, discrete interval
│   ├── estep.R              # E-step (base): responsibilities + sufficient statistics
│   ├── mstep.R              # M-step (base): Brent theta update + closed-form others
│   ├── em_driver.R          # em_fit_tenure() driver + init_params()
│   │
│   │  ── Contamination (epsilon) model ──
│   ├── emissions_eps.R      # Spell-pair point-mass + Exp emission (2^K enumeration)
│   ├── estep_eps.R          # E-step (eps): responsibilities + tau / spell statistics
│   ├── mstep_eps.R          # M-step (eps): closed-form eps, lambda_g; Brent for linked
│   ├── init_params_eps.R    # init_params_eps(): starting values (no sigma2_g)
│   ├── em_driver_eps.R      # em_fit_tenure_eps() driver
│   │
│   │  ── Shared utilities ──
│   ├── simulate.R           # simulate_panel() synthetic data generator
│   ├── diagnostics.R        # Post-estimation: multistart, profile LL, Hessian
│   ├── compare_distributions.R  # Real vs simulated comparison helpers
│   └── source_all.R         # Source all files in dependency order
├── estimate_pipeline.R      # End-to-end estimation (base + eps, all 4 variants each)
├── simulation_diagnostic.R  # Distributional comparison: real vs simulated data
├── tests/
│   └── testthat/
│       ├── helper-source.R
│       ├── test-utils.R
│       ├── test-transforms.R
│       ├── test-latent_histories.R
│       ├── test-emissions.R
│       ├── test-emissions-discrete.R
│       ├── test-estep-mstep.R
│       ├── test-eps-model.R       # ← Tests for the contamination model
│       ├── test-free-spec.R
│       └── test-integration.R
├── CONCERNS.md
├── SIMULATION_DIAGNOSTIC_REPORT.md
└── README.md
```

### Where to look for the epsilon model

| What | File(s) |
|------|---------|
| Theory / derivation | `documents/EM tenure epsilon.tex` |
| Spell-pair emission densities | `R/emissions_eps.R` |
| E-step (responsibilities, tau, sufficient stats) | `R/estep_eps.R` |
| M-step (ε, λ_g closed-form; Brent for linked θ) | `R/mstep_eps.R` |
| Starting values | `R/init_params_eps.R` |
| EM loop / convergence / monotonicity guard | `R/em_driver_eps.R` |
| Estimation pipeline (data → fit → save) | `estimate_pipeline.R` (search for `EPSILON`) |
| Tests | `tests/testthat/test-eps-model.R` |

## Usage

### Epsilon model quick start (synthetic data)

```r
source("EM-tenure/R/source_all.R")

# Initialise parameters
p0 <- init_params_eps(df, linked = FALSE, eps_init = 0.20)

# Fit the contamination model (free, non-stationary)
fit <- em_fit_tenure_eps(df, params0 = p0, stationary = FALSE,
                         linked = FALSE, verbose = 2L)
fit$params
# Returns: alpha, theta0, theta1, pi, eps, lambda_g, lambda_d
```

### With real data (full pipeline)

```r
source("EM-tenure/estimate_pipeline.R")
```

For epsilon point estimates only (three starts per specification, exact-cell
collapsing, weights normalized to the original sample size, and no bootstrap),
run `EM-tenure/estimate_point_tenure_contamination.R`.

This runs all model variants (base EM-tenure × 4 variants, epsilon × 4 variants), saves each fit as a timestamped `.rds` in `output/results/`, appends a summary row to `output/results/run_summary.csv`, and performs LR tests.

### Running tests

From the project root:

```r
testthat::test_dir("EM-tenure/tests/testthat")
```

## Key epsilon model functions

| Function | File | Description |
|----------|------|-------------|
| `em_fit_tenure_eps()` | `em_driver_eps.R` | Top-level EM driver (mirrors `em_fit_tenure()`) |
| `init_params_eps()` | `init_params_eps.R` | Build starting parameter vector (no σ²_g; adds ε) |
| `e_step_eps()` | `estep_eps.R` | E-step: responsibilities γ, spell-pair tau, sufficient stats |
| `m_step_eps()` | `mstep_eps.R` | M-step: closed-form ε, λ_g; Brent solver for linked θ |
| `log_spell_emission_eps()` | `emissions_eps.R` | Spell-pair joint tenure emission (2^K enumeration) |

The point-only runner `estimate_duration_hazard_tenure_contamination.R`
compares the constant-hazard linked epsilon model with the nested hazards
`h(x) = lambda * (1 + x)^beta`. It uses three starts, integrates transition
risk over interval-censored timegap categories, prints duration-specific and
posterior-risk-weighted rates, and does not run the bootstrap.

The follow-on runner `estimate_piecewise_hazard_tenure_contamination.R`
estimates positive piecewise-constant entry and exit hazards over 0--3 months,
3--12 months, 1--3 years, 3--5 years, and 5+ years. The final positive hazard
is maintained over the open tail, guaranteeing a proper distribution and a
finite mean duration. It reuses the corrected constant and power-law fits,
runs multiple starts plus a tighter refinement, and does not bootstrap.

`estimate_timegap_contamination.R` adds a parsimonious error mechanism for
follow-up timegap categories to that piecewise model. During a latent
nonemployment continuation, the new report follows the clock-consistent
category transition with probability `1 - eps_d`; with probability `eps_d` it
is an independent draw from the model-implied duration-category distribution.
The extension exactly nests the timegap-error-free model at `eps_d = 0`.
Marginal and isolated timegap reports do not identify `eps_d`. The runner uses
three starts, BFGS refinement, a conditional likelihood slice, and a numerical
observed-information diagnostic; it does not bootstrap.

`diagnose_timegap_contamination.R` decomposes posterior-weighted latent
nonemployment continuations without changing or re-estimating the model. It
separates clock-feasible reports, impossible transitions ending in the
0--3-month category that could reflect an unobserved employment spell and
clock reset, and impossible transitions that cannot be generated by such a
within-quarter interruption.

`estimate_timegap_contamination_robustness.R` compares two alternative error
distributions. The local model assigns exponentially less probability to
categories farther from the clean clock-reachable set. The joint model assigns
a clean/contaminated indicator to each report, permits clean reports to remain
linked around a contaminated middle wave, and reports both the per-wave error
probability and its implied two-report probability. Both models retain the same
piecewise entry and exit hazards. Generated fits are resumable and no bootstrap
is run.

`estimate_local_gross_tenure_contamination.R` extends the preferred joint
timegap specification by splitting tenure contamination into a fixed local
month-error kernel (symmetric errors of one to six months, with exponentially
declining weights) and the existing gross marginal draw. The local kernel is
discrete and fixed, avoiding variance collapse. A boundary profile and fixed
dispersed interior screen are reported because the local-error probability is
estimated at zero. Analytical errors for common quantities therefore inherit
the gross-only fit, while the boundary probability itself has no regular Wald
standard error.

`estimate_job_change_tenure_contamination.R` implements the first start-date
improvement without overwriting Table 7. For each latent employed-to-employed
boundary, it mixes a continuation of the current-job tenure clock with a new
job whose start lies uniformly within the quarterly interval. The single reset
probability is jointly estimated with the existing piecewise transition
hazards and status, tenure, and timegap contamination rates. The likelihood
exactly equals the established likelihood when the reset probability is zero;
the runner reports a conditional profile, the nested model comparison, and
posterior reset counts. Generated fits are stored under
`output/results/job_change/` and no bootstrap is run.

`estimate_job_change_matching_robustness.R` keeps Panel A as the main dataset
and re-estimates that identical 15-parameter extension on stricter matching
Panels B and C. It starts from the converged Panel A optimum and writes the
panel fits and comparison under `output/results/job_change/matching_robustness/`
without replacing any main result.

`validate_job_change_model.R` runs dispersed starts, tighter refinement, a
nuisance-adjusted reset-probability profile, multistart recovery simulations,
and a large-sample likelihood-score diagnostic. The current continuous-density
version is numerically stable on the real data but fails full-model recovery:
it understates simulated resets and overstates tenure contamination. Because
tenure is reported through a job-start month, the next required repair is a
discrete-month probability-mass likelihood. The job-change estimates should
remain diagnostic until that repair and its recovery checks pass.

`R/tenure_monthly_eps.R` implements that repair. Marginal tenure reports use
monthly interval probabilities, clean clocks advance by three months, and
latent employment entries and current-job resets have support at 0, 1, or 2
months. `validate_monthly_tenure_model.R` verifies normalization and performs a
full-likelihood recovery exercise. `estimate_monthly_job_change_tenure_contamination.R`
retains genuine zero-month reports, compares the no-reset and reset models on
Panel A, and screens three dispersed reset starts. Its outputs are stored under
`output/results/job_change_monthly/`; they do not overwrite Table 7.
`validate_monthly_tenure_recovery_grid.R` repeats the full-likelihood recovery
exercise at configurable sample sizes, caches every fit, screens a dispersed
start at each size, and reports bias, RMSE, hazard recovery, and scores at the
generating parameters. In the completed 10,000/40,000-observation grid, all six
fits converge and bias contracts sharply at the larger sample.
`infer_monthly_job_change_tenure_contamination.R` computes the full numerical
observed-information matrix in parallel, applies the delta method to the
paper-facing rates, and constructs a nuisance-adjusted profile interval for
the reset probability. `refine_monthly_profile_point.R` provides a resumable
single-point refinement when a profile tail needs a closer converged bracket.
`estimate_monthly_job_change_matching_robustness.R` re-estimates the identical
corrected model on Panels B and C, validates their monthly tenure grids, and
screens low, main, and high reset-probability starts. Its cached fits and the
A/B/C comparison are stored under
`output/results/job_change_monthly/matching_robustness/`.
`estimate_monthly_persistent_tenure_reporting.R` adds a one-parameter process
in which consecutive gross reports within a job can retain the preceding
reported start-date anchor. It checks exact nesting at zero, screens three
starts, and retains the nested boundary fit if it dominates all positive
approximations. `validate_monthly_persistent_tenure_reporting.R` checks that
positive persistence is recoverable using a conditional profile and dispersed
full-likelihood starts. Outputs are stored under
`output/results/job_change_monthly/persistent_reporting/`.
`estimate_monthly_calendar_heaping.R` adds a one-parameter reporting component
under which a gross tenure report is drawn from the marginal duration law
conditional on implying a January job-start month. Nominal interview months
are reconstructed from `period1`--`period3` and retained during cell collapse.
The extension normalizes exactly, nests the independent-gross model at zero,
and screens three starts without overwriting the existing monthly fit.
`validate_monthly_calendar_heaping.R` checks a conditional likelihood profile
and joint recovery from dispersed starts. Outputs are stored under
`output/results/job_change_monthly/calendar_heaping/`.
`infer_monthly_calendar_heaping.R` computes the full observed-information
matrix in resumable parallel chunks and propagates it to the reported rates by
the delta method. `estimate_monthly_calendar_revision.R` then adds a nested
whole-year date-revision branch: consecutive gross reports can preserve the
previously reported start month while selecting a different start year.
`validate_monthly_calendar_revision.R` checks normalization, conditional
identification, and joint recovery. Its outputs are stored under
`output/results/job_change_monthly/calendar_revision/`.
`infer_monthly_calendar_revision.R` computes the full 17-parameter observed
information matrix and delta-method uncertainty for this extension.
`profile_monthly_calendar_revision.R` constructs a nuisance-adjusted,
multistart likelihood profile for the whole-year revision probability.
`diagnose_monthly_calendar_revision.R` performs fitted-parameter parametric
predictive checks of January starts and within-person date changes, and
decomposes the likelihood gain into mutually exclusive observed patterns.
`estimate_monthly_clean_anchor_revision.R` adds a separate probability that a
gross current report revises the preceding reported start date when that
preceding report is clean. The earlier gross-to-gross revision model is nested
exactly at zero, and three dispersed starts are cached and compared.
`validate_monthly_clean_anchor_revision.R` checks the new channel with a
conditional likelihood profile and two joint recovery starts.
`diagnose_monthly_clean_anchor_revision.R` compares fitted calendar moments
with the nested model and decomposes the exact likelihood gain by observed
report pattern. `infer_monthly_calendar_revision.R` accepts
`REVISION_INFERENCE_MODEL=clean_anchor` to compute the full 18-parameter
observed information and delta-method uncertainty for this extension.
`estimate_monthly_start_month_baseline.R` adds a saturated 12-month baseline
distribution shared by clean and gross marginal tenure reports and by latent
entries and job resets. Equal month probabilities exactly nest the clean-anchor
model. Estimation uses three conditional month-effect starts on a fixed mode-
search sample, then forward- and central-difference joint refinements on the
full Panel A likelihood. `refine_monthly_start_month_baseline.R` continues a
full-sample central refinement without repeating the mode search.
`diagnose_monthly_start_month_baseline.R` compares fitted calendar moments and
decomposes the likelihood gain by observed start month.
`validate_monthly_start_month_baseline.R` runs dispersed-start and joint
simulation-recovery checks; its forward gradient uses a backward difference
when a positive perturbation enters an invalid penalty region.
`estimate_monthly_exact_anchor_retention.R` adds a shared probability that a
gross current tenure report retains the preceding reported job-start date
exactly, whether that preceding report was clean or gross. Whole-year revision
and marginal redraw probabilities are conditional on not retaining the anchor,
so the flexible baseline model is nested exactly at zero. Estimation combines
a wide conditional profile, bounded forward- and central-difference joint
refinement, and an exact one-dimensional polish of the new parameter.
`diagnose_monthly_exact_anchor_retention.R` compares exact continuations,
whole-year revisions, and other revisions in the observed and simulated data
and decomposes the likelihood gain by continuation pattern.
`validate_monthly_exact_anchor_retention.R` isolates the new channel in a
repeated simulation-recovery exercise with all nuisance parameters held at
their generating values.
`refine_monthly_exact_anchor_retention.R` performs a resumable longer
30-parameter refinement, screens low- and high-retention starts, promotes the
best joint fit, and applies a final scalar polish.
`validate_monthly_exact_anchor_joint_recovery.R` simulates the complete model
and re-estimates all 30 parameters from near and dispersed starts, caching each
replication and reporting recovery of headline rates, duration hazards, the
start-month baseline, and reporting-process probabilities separately.
`infer_monthly_exact_anchor_retention_opg.R` obtains a full-rank analytical OPG
approximation from record-level finite-difference scores and reports
delta-method errors, boundary flags, and reporting-probability correlations.
`infer_monthly_exact_anchor_retention.R` provides the slower, resumable exact
central observed-Hessian calculation for confirmatory inference.
`estimate_monthly_local_anchor_revision.R` adds a normalized non-zero local
revision kernel over plus or minus one to six months, conditional on not
retaining the preceding reported anchor exactly. The fixed kernel has a
three-month discrete-Laplace scale and is renormalized when backward revisions
would imply negative tenure. The estimator verifies exact nesting at zero,
profiles the new probability, screens dispersed starts, jointly refines the
competing reporting parameters, and checks the transition-hazard block. The
calendar fitter's optional `free_names` argument makes these blockwise checks
resumable while leaving full-parameter fitting as the default.
`diagnose_monthly_local_anchor_revision.R` compares exact, one-to-six-month,
whole-year, and residual revisions in observed and simulated data, decomposes
the likelihood gain by continuation pattern, and computes a conditional
curvature diagnostic. Outputs are stored under
`output/results/job_change_monthly/calendar_local_anchor_revision/`.
`estimate_monthly_shared_duration_reliability.R` adds an equal-weight
person-level reliability mixture that shifts tenure and timegap contamination
logits together, verifies exact nesting at zero dispersion, and screens
dispersed starts. `diagnose_monthly_shared_duration_reliability.R` evaluates
duration-pattern predictions, likelihood contributions, and posterior class
probabilities. These provisional outputs are stored under
`output/results/job_change_monthly/shared_duration_reliability/`.
`estimate_monthly_separate_duration_reliability.R` relaxes the common-shift
restriction, verifies exact nesting when the two clock-specific dispersions
are equal, and screens common, tenure-heavy, and timegap-heavy starts.
`diagnose_monthly_separate_duration_reliability.R` compares the nested and
unrestricted predictive distributions and decomposes their likelihood
difference. `refine_monthly_separate_duration_reliability.R`,
`refine_monthly_separate_reliability_blocks.R`, and
`refine_monthly_separate_reliability_joint.R` provide cached focal, blockwise,
and full joint continuations. `refine_monthly_common_duration_reliability.R`
and `polish_monthly_common_duration_reliability.R` put the equality-restricted
fit on the same optimized nuisance base. The retained unrestricted and
restricted fits both converge, and their comparison no longer conflates the
restriction with unrelated optimizer progress.
`validate_monthly_separate_duration_reliability.R` performs multiple-start
conditional recovery and posterior class-separation checks.
`infer_monthly_separate_duration_reliability_opg.R` computes cached
full-nuisance OPG/delta-method inference with adaptive finite differences and
information diagnostics. `profile_monthly_separate_duration_reliability.R`
profiles the timegap-minus-tenure dispersion difference while reoptimizing the
other 32 parameters and screens each outer bracket from two nuisance starts.
`refine_monthly_separate_duration_reliability_profile_limits.R` adds tightly
optimized points near the two likelihood-ratio confidence limits. Outputs are
stored under
`output/results/job_change_monthly/separate_duration_reliability/`.
`bootstrap_monthly_separate_duration_reliability.R` provides a resumable
parametric-bootstrap calibration of the equality-restriction likelihood-ratio
test and a separate full-nuisance recovery exercise. It caches simulated data
and every optimizer stage, reports convergence explicitly, and excludes capped
fits from the bootstrap p-value. Its sample size, replication count, workers,
and continuation budget are controlled by the `RELIABILITY_BOOTSTRAP_*`
environment variables. `RELIABILITY_BOOTSTRAP_RUN_NULL` and
`RELIABILITY_BOOTSTRAP_RUN_RECOVERY` can disable either exercise so that
large recovery and calibration batches can be run separately.
`RELIABILITY_BOOTSTRAP_RELTOL` and `RELIABILITY_BOOTSTRAP_PGTOL` expose the
optimizer tolerances (defaulting to the fitter's standard `1e-9` and `1e-7`).
`RELIABILITY_BOOTSTRAP_RECOVERY_BLOCKWISE=true` continues each unrestricted
recovery fit through calendar-reporting, reliability, and transition blocks
before a final joint optimization. `RELIABILITY_BOOTSTRAP_BLOCK_CYCLES`,
`RELIABILITY_BOOTSTRAP_BLOCK_MAXIT`, and
`RELIABILITY_BOOTSTRAP_BLOCK_FULL_MAXIT` control that continuation.
Draw-level files include the last likelihood gain and gain per observation so
iteration-capped fits can be diagnosed without treating them as converged.
New cache files also record their stage-specific optimizer controls; this
allows later continuation blocks to use different budgets or tolerances
without obscuring how an earlier stage was produced.

`R/four_wave_eps.R` generalizes the corrected discrete-month observed
likelihood from eight three-wave employment histories to all sixteen
four-wave histories. It retains the piecewise duration hazards, current-job
resets, calendar reporting channels, timegap contamination, and the fixed
equal-weight reliability mixture. The established three-wave E-step is left
unchanged. This module supplies a four-wave AR(1) baseline that can be checked
before adding second-order transition dependence.

Run `estimate_monthly_four_wave_ar1.R` from the project root. It saves data,
bounded refinement blocks, sample flow, and a provisional comparison under
`output/results/job_change_monthly/four_wave_ar1/`. Set
`FOUR_WAVE_AR1_REPORT_ONLY=true` to collect existing checkpoints without
estimating another block. The default uses two workers and two cycles of
transition, reporting, and calendar blocks; `FOUR_WAVE_AR1_CYCLES` controls
the number of cycles. Existing checkpoints are reused. A bounded step is
explicitly **not convergence**. `FOUR_WAVE_AR1_FULL_MAXIT` is zero by default;
positive values request an additional, potentially lengthy L-BFGS-B polish.
The same-person three-wave row currently evaluates the inherited parameters,
not a re-estimated three-wave comparator. No AR(2) estimates or analytical
standard errors are produced by this provisional runner.

The monthly tenure emission shares job-segment calculations within each
likelihood evaluation and evaluates identical segment records only once.
This is an exact computational refactor, tested against the uncompressed
calculation and a saved full-sample four-wave posterior; it does not change
the duration or reporting assumptions.

For joint four-wave refinement, run `validate_four_wave_fast.R` followed by
`refine_four_wave_ar1_fast.R`. The optional Rcpp/OpenMP evaluator in
`R/four_wave_fast.R` and `src/four_wave_monthly.cpp` uses probability tables
from the reference R implementation. Validation covers perturbed parameters,
zero status error, disabled job/report-revision channels, serial/parallel
agreement, and the full-sample saved posterior. The refinement runner requires
validation fingerprints matching the current sources. Rcpp and an R-compatible
C++ compiler are required for the first build; compiled files stay in the
output cache, not in the source directory.
Inputs are the cached four-wave cells and preliminary fit produced by
`estimate_monthly_four_wave_ar1.R`, plus the existing preferred three-wave
information matrix for numerical scaling. These are prerequisites, not
additional models re-estimated by the continuation.

The joint runner saves improving parameter checkpoints and reports a separate
central-difference, bound-adjusted score check at termination. Optimizer code
zero alone is not sufficient to flag convergence. `FOUR_WAVE_FAST_LABEL`,
`FOUR_WAVE_FAST_START`, `FOUR_WAVE_FAST_MAXIT`, `FOUR_WAVE_FAST_WORKERS`, and
`FOUR_WAVE_FAST_PERTURB` control reproducible continuations and alternative
starts. Checkpoint resumption defaults to true; set `FOUR_WAVE_FAST_RESUME=false`
to restart a label. It neither estimates nor evaluates a three-wave comparator.

The four-wave continuation uses exact piecewise-exponential integrals for
category-averaged and missing-duration transition risks. Fifty independent
integration/limit checks cover this calculation. The reference four-wave
E-step supports it with `exact_risk=TRUE`; its default legacy quadrature and
the three-wave estimator are unchanged. This eliminates small quadrature
noise from numerical scores without changing the mathematical risk model.
`FOUR_WAVE_FAST_GRADIENT` defaults to `central` (`forward` is also supported).
The saved three-wave information diagonal is used solely for parameter
scaling, with a cap on very weak directions. `summarise_four_wave_ar1_fast.R`
collects completed runs, checks across starts, and consistently calculated
four-wave rates without refitting any other model.
`check_four_wave_ar1_fast.R` runs the focused regression tests and checks the
audit result-table chunk without estimating a model or rendering a PDF.

**Repaired four-wave timegap clock:** four-wave evaluators now default to
`continuous_joint`. `R/timegap_clock_joint.R` integrates one spell-start clock
over the intersection of every clean report's shifted category interval.
It sums over contamination flags and integrates unobserved reports, without
the old joint/pairwise switch. The initial timegap distribution, reporting
error probabilities, transition-risk formulas, tenure model, and three-wave
estimator are unchanged. The existing sequential simulator with
`waves=4, exact_risk=TRUE` implements this continuous-clock law.

`validate_four_wave_simulation.R`, `audit_four_wave_probability_mass.R`, and
`summarise_four_wave_recovery.R` now write to `recovery_continuous_clock`.
The repaired probability-mass and score gates pass. Set
`FOUR_WAVE_MASS_MAX_MONTH=24` for the tight-tail check; use 12 and
`FOUR_WAVE_MASS_REFERENCE=true` for R/C++ comparison. Set
`FOUR_WAVE_TIMEGAP_CLOCK=legacy_pairwise` to reproduce the old diagnostics
in `recovery`; their failed results are retained. Direct evaluator calls
can use `timegap_clock="legacy_pairwise"` to reproduce historical likelihoods.

`estimate_four_wave_recovery.R` requires successful current-source gates and
compiled validation, then estimates all 33 parameters on two 5,000-person
simulations, with an additional start for the first sample. It saves fits,
checkpoints, central-score checks, rates, and hazard coefficients separately.
`FOUR_WAVE_RECOVERY_N`, `FOUR_WAVE_RECOVERY_FIT_REPS`,
`FOUR_WAVE_RECOVERY_MAXIT`, and `FOUR_WAVE_RECOVERY_WORKERS` control this small
exercise. This is not a bootstrap standard-error calculation.

Historical empirical fit objects are retained unchanged. The empirical continuation's
new default label is `clock_joint1`, and it rejects stale-source checkpoints.
The historical summary collector excludes fits explicitly tagged with the
repaired clock, to avoid combining different likelihoods in a start comparison.
For a repaired Panel A continuation, set `FOUR_WAVE_FAST_START` to
`EM-tenure/output/results/job_change_monthly/four_wave_ar1/fit_fast_start_plus_latest.rds`
and run `refine_four_wave_ar1_fast.R`. The start is the best converged historical
fit, not an assertion that its old likelihood is the repaired objective.
After convergence, run `summarise_four_wave_clock_empirical.R`. It requires
current-source recovery gates and an independently converged repaired fit,
checks the full-sample R likelihood and posteriors at the new parameters,
and writes main-rate comparisons, interval rates, annual duration hazards,
transformed coefficients, and a consolidated RDS to
`four_wave_ar1/continuous_clock_empirical/`. Main rates use the same posterior
origin-risk and survey weights as before. Cross-formula likelihood differences
are not LR statistics. This stage does not estimate standard errors or AR(2).

### Four-wave AR(2) duration extension

The separate `R/four_wave_ar2.R` module adds two log-odds return effects to
the repaired four-wave AR(1) law. Effects apply only at the last two transitions
when the preceding latent state differs from the current state. The first
transition retains its AR(1) initialization, with free initial employment.
The duration/reporting laws and equal reliability-class shares are unchanged.
Effects are odds multipliers of the existing integrated quarterly risks, not
hazard multipliers. Both zero exactly nest AR(1).

The evaluator exactly reweights the 16 AR(1) history posteriors at each new
nuisance vector. It is not a fixed-posterior or two-step fit. Final class/job
posterior aggregates are obtained from a direct-prior R reference. The
separate sequential simulator retains the original reporting code; a seeded
zero-effect test checks synchronization with the unchanged AR(1) simulator.

Run from the repository root, in order:

1. `EM-tenure/validate_four_wave_ar2.R` (nesting, reference, normalization, scores);
2. `EM-tenure/recover_four_wave_ar2.R` (all 35 parameters, two simulated samples);
3. `EM-tenure/estimate_four_wave_ar2.R` (Panel A, three starts, independent scores);
4. `EM-tenure/summarise_four_wave_ar2.R` (rates, effects, identification and reversal decomposition).

Results, input/source fingerprints, and resumable checkpoints are stored
separately in `output/results/job_change_monthly/four_wave_ar2/`. A conditional
two-coefficient optimization supplies a warm start only; the empirical fits
free all 35 parameters. Search gradients are forward differences, but acceptance
requires an independent central projected score at most `1e-5` and optimizer
code zero. AR(2) coefficients are bounded at +/-8; boundary solutions should
be flagged rather than interpreted as well-determined effects. No bootstrap
or analytical standard errors are calculated by these runners.

### Key arguments to `em_fit_tenure_eps()`

| Argument | Default | Description |
|----------|---------|-------------|
| `df` | — | Data frame (see data requirements below) |
| `params0` | `NULL` | Starting parameters; if NULL, calls `init_params_eps()` |
| `stationary` | `FALSE` | Impose α = θ₀/(θ₀ + 1 − θ₁) |
| `linked` | `FALSE` | Impose CTMC link: λ_g = −log(θ₁)/Δ, λ_d = −log(1−θ₀)/Δ |
| `max_iter` | `500` | Maximum EM iterations |
| `tol` | `1e-8` | Convergence tolerance (relative ΔLL) |
| `verbose` | `0` | 0 = silent, 1 = final summary, 2 = per-iteration |

## Data requirements

The input data frame must have columns:

| Column | Description | Required |
|--------|-------------|:--------:|
| y1, y2, y3 | Observed employment (0/1) at waves 1-3 | Always |
| tenure1, tenure2, tenure3 | Observed tenure in **years** | Always |
| timegap_cat1, timegap_cat2, timegap_cat3 | Nonemployment category codes (1–7) | Always |
| weight | Survey weight | Always |

### Missing wrong-state durations

The shared ingest script historically fills durations belonging to the
unreported state for compatibility with older estimators. The epsilon model no
longer treats these constructed values as data. `prepare_eps_estimation_data()`
sets tenure to missing at reported-nonemployment waves and timegap to missing
at reported-employment waves. Emissions integrate these unavailable clocks out.
In the duration-dependent transition prior, risk is averaged over the
model-implied duration distribution whenever the latent state differs from the
reported state. Actually observed tenure and timegap reports retain the exact
clock-consistency restrictions.

## TeX references

| Specification | Document |
|---------------|----------|
| Base EM-tenure model | `documents/EM tenure.tex` |
| **Contamination (epsilon) model** | `documents/EM tenure epsilon.tex` |

**Epsilon model TeX cross-references**:
- Contamination model definition: Section 3.2 (Eq. for ε branch)
- Spell-pair joint emission: Section 3.4 (Eq. spell_emission)
- Per-pattern density (clean/contaminated decomposition): Section 3.4 (Eq. per_pattern_density)
- E-step sufficient statistics: Section 4
- M-step updates (ε, λ_g): Section 4.3
