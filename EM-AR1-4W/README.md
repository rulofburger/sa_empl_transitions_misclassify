# Four-wave AR(1) benchmark

Single latent-type first-order Markov employment model estimated on all complete
four-wave QLFS histories. The module supports no measurement error and symmetric
misclassification, with either a stationary or freely estimated initial
employment probability. Point estimates use exact 16-cell observed-data MLE;
standard errors use an individual survey-weighted sandwich covariance and the
delta method. The probability kernel and parameter transforms are intended to
be reused by the subsequent two-type finite-mixture model.

Run from the project root with `source("EM-AR1-4W/estimate_pipeline.R")`.

The two-type finite-mixture extension is in `R/fmm_ar1_4w.R`. It uses
multi-start EM over latent types and employment histories, followed by an exact
observed-likelihood polish, label normalization, a Jacobian-rank identification
test, and analytical sandwich/delta-method inference. Run it with
`source("EM-AR1-4W/estimate_fmm_pipeline.R")`.

The parsimonious controlled-heterogeneity specification is in
`R/fmm_covariates_inconsistency_4w.R`. It retains two type-specific transition
intercepts and free initial conditions, while constraining demographic and job
covariate slopes to be common across types. Origin-wave log tenure enters
persistence only; log time since work and never-worked enter entry only.
Permanent contract and informal sector vary by transition and enter
persistence only. Symmetric
misclassification varies with wave-attributed age and education
inconsistencies. The estimator uses an EM warm start, exact observed-likelihood
analytical scores, multiple starts, and survey-weighted sandwich/delta-method
inference. Run it with
`source("EM-AR1-4W/estimate_fmm_covariates_inconsistency.R")`.

For the robustness model with type-specific baseline misclassification but
common inconsistency slopes, run
`source("EM-AR1-4W/estimate_fmm_type_specific_misclassification.R")`.

For the two Table 6 variants with type-specific error intercepts, run
`source("EM-AR1-4W/estimate_fmm_type_specific_error_variants.R")`. The first
has no observed predictors in the error equation. The second uses the Table 3
column (1) age, education, race, and gender inconsistency indicators and the
exactly-two, exactly-three, and exactly-four indicators, with common slopes and
type-specific intercepts.

The runner writes one saved fit per variant. Use
`polish_fmm_type_specific_error_variant.R <variant>` for an additional exact-
likelihood BFGS refinement and
`finalize_fmm_type_specific_error_variants.R <saved-fit-stem>` to rebuild the
sandwich/delta-method inference and audit the score, rank, and curvature. If the
full Table 3-predictor fit has negative local curvature,
`escape_fmm_table3_saddle.R` profiles the weakest information-matrix direction
and restarts from an improving point before it is finalized again. Table 6 uses
only fits with a maximum normalized score below `1e-5`, full information rank,
and positive minimum curvature.

The uncontrolled reliability-mixture specification is in
`R/fmm_common_transitions_type_error_4w.R`. It holds entry, exit, and initial
employment probabilities common across the two latent groups and allows only
the symmetric misclassification probability to differ. It uses the same
covariate-complete sample as Table 6 for comparison, without putting those
covariates in the likelihood. Run it with
`source("EM-AR1-4W/estimate_fmm_common_transitions_type_error.R")`.
