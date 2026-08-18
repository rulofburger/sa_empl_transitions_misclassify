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
covariate slopes to be common across types. Permanent contract and informal
sector vary by transition and enter persistence only. Symmetric
misclassification varies with wave-attributed age and education
inconsistencies. The estimator uses an EM warm start, exact observed-likelihood
analytical scores, multiple starts, and survey-weighted sandwich/delta-method
inference. Run it with
`source("EM-AR1-4W/estimate_fmm_covariates_inconsistency.R")`.

For the robustness model with type-specific baseline misclassification but
common inconsistency slopes, run
`source("EM-AR1-4W/estimate_fmm_type_specific_misclassification.R")`.
