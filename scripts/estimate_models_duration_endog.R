

df_estimate <- df_qlfs_tenure_timegap

df_covariate_combos <- df_estimate %>%
  select(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3) %>%
  unique() %>%
  na.omit %>%
  data.table()

df_template_employment <- data.table::CJ(
  y1_star = c(0, 1),
  y2_star = c(0, 1),
  y3_star = c(0, 1)
)

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_duration_endog <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_duration_endog[, k := NULL]

# param_init <- data.frame(intercept_0 = -2.524934, timegap = -0.1300759, intercept_1 = 1.02172, tenure = -0.6729012, pi = 0.02924524)
param_init <- data.frame(theta_0 = 0.02881745, theta_1 = 0.02918875, lambda_g = 6.439774, lambda_h = 2.298665, sigma_g = 0.2231302, sigma_h = 0.2231302, pi = 0.02990933)
param_init_transformed <- param_init
param_init_transformed$theta_0 <- logit_transform(param_init$theta_0)
param_init_transformed$theta_1 <- logit_transform(param_init$theta_1)
param_init_transformed$pi <- log(param_init$pi/(0.5 - param_init$pi))
param_init_transformed$lambda_g <- log(param_init$lambda_g)
param_init_transformed$lambda_h <- log(param_init$lambda_h)
param_init_transformed$sigma_g <- log(param_init$sigma_g)
param_init_transformed$sigma_h <- log(param_init$sigma_h)

calc_mle_3waves_ar1_duration_endog(param_init_transformed)
# calc_lli_3waves_ar1_duration_endog(param_init_transformed)
# calc_lli_derivatives_3waves_ar1_duration_endog(param_init_transformed)
# calc_mle_derivatives_3waves_ar1_duration_endog(param_init_transformed)
# 
# param_init_transformed2 <- data.frame(theta_0 = -1.085932, theta_1 = -2.433053, pi = 28.34738, lambda_g = 1.804956, lambda_h = 0.8556096 , sigma_g = -415.6423, sigma_h = -260.3311)
# 
# calc_mle_3waves_ar1_duration_endog(param_init_transformed2)
# llia <- calc_lli_3waves_ar1_duration_endog(param_init_transformed)
# llib <- calc_lli_3waves_ar1_duration_endog(param_init_transformed2)

calc_mle_derivatives_3waves_ar1_duration_endog(param_init_transformed)

pi_fixed <- list(pi = param_init_transformed$pi)
param_init_transformed_pi_fixed <- param_init_transformed
param_init_transformed_pi_fixed$pi <- NULL

model_mle_3w_ar1_duration_endog <- maxLik::maxLik(
  logLik = function(param) calc_mle_3waves_ar1_duration_endog_pi_fixed(param, pi_fixed = pi_fixed),
  grad = function(param) calc_mle_derivatives_3waves_ar1_duration_endog_pi_fixed(param, pi_fixed = pi_fixed),
  start = param_init_transformed_pi_fixed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_duration_endog$estimate
model_mle_3w_ar1_duration_endog$maximum

only_pi <- model_mle_3w_ar1_duration_endog$estimate

param_init_transformed_only_pi <- pi_fixed
pi_fixed <- NULL

model_mle_3w_ar1_duration_endog <- maxLik::maxLik(
  logLik = function(param) calc_mle_3waves_ar1_duration_endog_only_pi(param, only_pi = only_pi),
  grad = function(param) calc_mle_derivatives_3waves_ar1_duration_endog_only_pi(param, only_pi = only_pi),
  start = param_init_transformed_only_pi,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_duration_endog <- maxLik::maxLik(
  logLik = calc_mle_3waves_ar1_duration_endog,
  grad = calc_mle_derivatives_3waves_ar1_duration_endog,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

param_init_transformed_new <- model_mle_3w_ar1_duration_endog$estimate
param_init_transformed_new$pi <- pi_fixed$pi
calc_mle_derivatives_3waves_ar1_duration_endog(param_init_transformed_new)

model_mle_3w_ar1_duration_endog <- maxLik::maxLik(
  calc_mle_3waves_ar1_duration_endog,
  grad = calc_mle_derivatives_3waves_ar1_duration_endog,
  start = model_mle_3w_ar1_duration_endog$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_duration_endog$estimate
model_mle_3w_ar1_duration_endog$maximum
