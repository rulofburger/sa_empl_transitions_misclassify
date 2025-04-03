# METADATA =====================================================================
# DESCRIPTION: Estimate 3 wave models for SA
# CREATED: 2025-02-13 (rulofburger)

# SUMMARY:

# INITIALISE ====

# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)
library(haven)
library(data.table)
library(fastverse)

# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

#> Load estimation functions defined in other scripts ----
source("scripts/define_estimation_functions_3waves_mle_ar1.R")
source("scripts/define_estimation_functions_3waves_mle_ar1_covariates.R")
source("scripts/define_estimation_functions_3waves_mle_ar1_misclassification_symptoms.R")
source("scripts/define_estimation_functions_3waves_mle_ar1_fmm_hessian.R")

#> Define functions to create output tables ----

# These are used because I couldnt find a better way to output maxLik
# coefficients and the other useful statistics into a Latex table. 

# Latex output table for simple models without covariates
create_stargazer_table <- function(model_object, df, formula) {
  # Ensure required libraries are loaded
  if (!requireNamespace("maxLik", quietly = TRUE) || 
      !requireNamespace("stargazer", quietly = TRUE)) {
    stop("Required packages 'maxLik' and 'stargazer' are not installed.")
  }
  
  # Extract coefficients and calculate standard errors
  model_se <- sqrt(((exp(model_object$estimate) / ((1 + exp(model_object$estimate))^2))^2) * 
                     diag(vcov(model_object)))
  
  # Create a formula dynamically
  formula <- as.formula(formula)
  
  # Fit an auxiliary linear model
  lm_model <- lm(data = df, formula)
  
  # Transform coefficients to percentages and standard errors accordingly
  coefficients_transformed <- unlist(logit_inverse(model_object$estimate))*100
  std_errors_transformed <- unlist(model_se*100)
  
  # Update the linear model with transformed coefficients and residuals
  lm_model$coefficients <- coefficients_transformed
  obs <- nrow(df)
  lm_model$residuals <- rnorm(obs)
  lm_model$se_list <- std_errors_transformed
  lm_model$ll <- model_object$maximum
  return(lm_model)
  
}

# Latex output table for simple models with covariates that determine transition rates
create_stargazer_table_covariates <- function(model_object, df, df_transition_probs, mean_transition_rates, formula1, formula2) {
  
  mat_X1 <- model.matrix(as.formula(formula1), df)
  mat_X2 <- model.matrix(as.formula(formula2), df)
  
  k_covariates <- dim(mat_X1)[2]
  
  model_vcov <- vcov(model_object)
  model_se_pi <- sqrt((((logit_inverse(model_object$estimate)$pi)*(1-logit_inverse(model_object$estimate)$pi))^2)*model_vcov["pi", "pi"])
  
  cov_matrix <- model_vcov[-dim(model_vcov)[1],-dim(model_vcov)[1]]
  
  cov_matrix1 <- cov_matrix[seq(1:k_covariates), seq(1:k_covariates)]
  cov_matrix2 <- cov_matrix[seq(k_covariates + 1, 2*k_covariates), seq(k_covariates + 1, 2*k_covariates)]
  
  grad0_1 <- df_transition_probs$theta0_1 * (1 - df_transition_probs$theta0_1) * mat_X1
  grad1_1 <- df_transition_probs$theta1_1 * (1 - df_transition_probs$theta1_1) * mat_X1
  grad0_2 <- df_transition_probs$theta0_2 * (1 - df_transition_probs$theta0_2) * mat_X2
  grad1_2 <- df_transition_probs$theta1_2 * (1 - df_transition_probs$theta1_2) * mat_X2
  
  var_pred_probs0_1 <- rowSums(grad0_1 %*% cov_matrix1 * grad0_1)
  var_pred_probs1_1 <- rowSums(grad1_1 %*% cov_matrix1 * grad1_1)
  var_pred_probs0_2 <- rowSums(grad0_2 %*% cov_matrix2 * grad0_2)
  var_pred_probs1_2 <- rowSums(grad1_2 %*% cov_matrix2 * grad1_2)
  
  # Standard errors of predicted probabilities
  se_pred_probs0_1 <- sqrt(var_pred_probs0_1)
  se_pred_probs1_1 <- sqrt(var_pred_probs1_1)
  se_pred_probs0_2 <- sqrt(var_pred_probs0_2)
  se_pred_probs1_2 <- sqrt(var_pred_probs1_2)
  
  N <- length(df_transition_probs$theta0_1)
  
  # SE of the mean predicted probability
  se_mean_pred0_1 <- sqrt(mean(var_pred_probs0_1) / N)
  se_mean_pred1_1 <- sqrt(mean(var_pred_probs1_1) / N)
  se_mean_pred0_2 <- sqrt(mean(var_pred_probs0_2) / N)
  se_mean_pred1_2 <- sqrt(mean(var_pred_probs1_2) / N)
  
  formula <- "y3 ~ y1 + y2 + educ1 + educ2"
  
  # Create a formula dynamically
  formula <- as.formula(formula)
  
  # Fit an auxiliary linear model
  lm_model <- lm(data = df, formula)
  
  # Transform coefficients to percentages and standard errors accordingly
  coefficients_transformed <- unlist(list(
    theta0_1 = mean_transition_rates$theta0_1,
    theta0_2 = mean_transition_rates$theta0_2,
    theta1_1 = mean_transition_rates$theta1_1,
    theta1_2 = mean_transition_rates$theta1_2,
    pi = logit_inverse(model_object$estimate$pi)
  ))*100
  
  
  std_errors_transformed <- unlist(list(
    theta0_1 = se_mean_pred0_1,
    theta0_2 = se_mean_pred0_2,
    theta1_1 = se_mean_pred1_1,
    theta1_2 = se_mean_pred1_2,
    pi = model_se_pi
  ))*100
  
  
  lm_model$coefficients <- coefficients_transformed
  obs <- nrow(df)
  lm_model$residuals <- rnorm(obs)
  lm_model$se_list <- std_errors_transformed
  lm_model$ll <- model_object$maximum
  return(lm_model)
}

# Latex output table for simple models with covariates that determine misclassification prob
create_stargazer_table_symptoms <- function(model_object, df, formula) {
  
  formula <- as.formula(formula)
  covariate_matrix  <- model.matrix(as.formula(formula), df) %>% unique()
  
  model_vcov <- vcov(model_object)
  model_se_theta_0 <- sqrt((((logit_inverse(model_object$estimate)$theta_0)*(1-logit_inverse(model_object$estimate)$theta_0))^2)*model_vcov["theta_0", "theta_0"])
  model_se_theta_1 <- sqrt((((logit_inverse(model_object$estimate)$theta_1)*(1-logit_inverse(model_object$estimate)$theta_1))^2)*model_vcov["theta_1", "theta_1"])
  
  vcov_matrix <- model_vcov[-c(1:2),-c(1:2)]
  coef_estimates <- model_object$estimate[-c(1:2)]
  k_covariates <- dim(vcov_matrix)[1]
  
  linear_predictor <- covariate_matrix %*% t(as.matrix(coef_estimates))
  
  variances_linear_predictor <- apply(
    covariate_matrix, 
    1, 
    function(row) t(row) %*% vcov_matrix %*% row
  )
  se_linear_predictors <- sqrt(variances_linear_predictor)
  predicted_probabilities <- 1 / (1 + exp(-linear_predictor))
  
  gradient <- predicted_probabilities * (1 - predicted_probabilities)
  variances_probabilities <- (gradient^2) * variances_linear_predictor
  se_probabilities <- sqrt(variances_probabilities)
  
  formula <- update(formula, y3 ~  . + y2 + y1)
  
  
  # Fit an auxiliary linear model
  lm_model <- lm(data = df, formula)
  
  names(lm_model$coefficients)[-c(1:k_covariates)] <- c("theta_0", "theta_1")
  
  lm_model$coefficients["theta_0"] <- logit_inverse(model_object$estimate$theta_0)
  lm_model$coefficients["theta_1"] <- logit_inverse(model_object$estimate$theta_1)
  
  for (j in seq(1, k_covariates)) {
    lm_model$coefficients[j] <- predicted_probabilities[j]
  }
  
  std_errors_transformed <- lm_model$coefficient
  std_errors_transformed["theta_0"] <- model_se_theta_0
  std_errors_transformed["theta_1"] <- model_se_theta_1
  for (j in seq(1, k_covariates)) {
    std_errors_transformed[j] <- se_probabilities[j]
  }  
  
  lm_model$coefficients <- lm_model$coefficients*100
  std_errors_transformed <- std_errors_transformed*100
  
  # Update the linear model with transformed coefficients and residuals
  obs <- nrow(df)
  lm_model$residuals <- rnorm(obs)
  lm_model$se_list <- std_errors_transformed
  lm_model$ll <- model_object$maximum
  return(lm_model)
  
}

# INGEST DATA ====

# Run script that loads 3 wave SA data as df_qlfs
source("scripts/ingest_data_3waves_SA.R")

# Limit survey rounds and calculate weights to be consistent within panel and to sum to 1
df_qlfs <- df_qlfs %>% 
  filter(period1 >= 30 & period1 <= 32) %>% 
  mutate(weight_total = sum(weight))  %>% 
  mutate(weight = dim(df_qlfs)[1]*weight/weight_total) 

# Create data set with covariates for specific models estimated below
df_qlfs_age_educ <- df_qlfs %>% 
  filter(!is.na(age1) & !is.na(age2) & !is.na(age3)) %>% 
  filter(!is.na(educ1) & !is.na(educ2) & !is.na(educ3))

df_qlfs_tenure_timegap <- df_qlfs %>% 
  filter(!is.na(timegap1) & !is.na(timegap2) & !is.na(timegap3)) %>% 
  filter(!is.na(tenure1) & !is.na(tenure2) & !is.na(tenure3))

df_qlfs_age_educ_tenure_timegap <- df_qlfs_age_educ %>% 
  filter(!is.na(timegap1) & !is.na(timegap2) & !is.na(timegap3)) %>% 
  filter(!is.na(tenure1) & !is.na(tenure2) & !is.na(tenure3))

df_qlfs_age_educ_female_race <- df_qlfs_age_educ %>% 
  filter(!is.na(race1) & !is.na(race2) & !is.na(race3)) %>% 
  filter(!is.na(female1) & !is.na(female2) & !is.na(female3))

df_qlfs_age_educ_female_race_contract <- df_qlfs_age_educ_female_race %>%
  mutate(
    contracttype1_missing = if_else(is.na(contracttype1), 1, 0),
    contracttype1 = if_else(is.na(contracttype1), 0, contracttype1),
    contracttype2_missing = if_else(is.na(contracttype2), 1, 0),
    contracttype2 = if_else(contracttype2_missing == 1, 0, contracttype2)
    )

# ESTIMATE MODELS ====

#> 3 waves, AR(1), no covariates ==== 

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 
df_estimate <- df_qlfs

#>> No ME ----

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085)
param_init_transformed <- logit_transform(param_init)
model_mle_3w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_pi0,
  grad = calc_mle_derivatives_3waves_ar1_pi0,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_pi0,
  grad = calc_mle_derivatives_3waves_ar1_pi0,
  start = param_init_transformed,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
print(logit_inverse(model_mle_3w_ar1_pi0$estimate))

#>> ME ----

param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, pi = 0.03)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1 <- maxLik::maxLik(
  calc_mle_3waves_ar1,
  grad = calc_mle_derivatives_3waves_ar1,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1 <- maxLik::maxLik(
  calc_mle_3waves_ar1,
  grad = calc_mle_derivatives_3waves_ar1,
  start = model_mle_3w_ar1$estimate,
  method = "NM",
  reltol = 0,
  gradtol = 0
)

print(logit_inverse(model_mle_3w_ar1$estimate))

#>> Asymmetric ME ----

param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, pi_0 = 0.03, pi_1 = 0.03)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1_asymmetric <- maxLik::maxLik(
  calc_mle_3waves_ar1_asymmetric,
  grad = calc_mle_derivatives_3waves_ar1_asymmetric,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_asymmetric <- maxLik::maxLik(
  calc_mle_3waves_ar1_asymmetric,
  grad = calc_mle_derivatives_3waves_ar1_asymmetric,
  start = model_mle_3w_ar1_asymmetric$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
print(logit_inverse(model_mle_3w_ar1_asymmetric$estimate))

#>> Creat Latex table ====

# AR(1) model, 3 waves, no ME
lm_model_3w_ar1_pi0 <- create_stargazer_table(model_object = model_mle_3w_ar1_pi0, df = df_estimate, formula = "y3 ~ y1 + y2") 

# AR(1) model, 3 waves, ME
lm_model_3w_ar1 <- create_stargazer_table(model_object = model_mle_3w_ar1, df = df_estimate, formula = "y3 ~ y1 + y2 + y3") 

# AR(1) model, 3 waves, Assymetric ME
lm_model_3w_ar1_asymmetric <- create_stargazer_table(model_object = model_mle_3w_ar1_asymmetric, df = df_estimate, formula = "y3 ~ y1 + y2 + y3 + age1") 


# Generate the stargazer table
table_simple_implied <- stargazer::stargazer(
  lm_model_3w_ar1_pi0, lm_model_3w_ar1, lm_model_3w_ar1_asymmetric, 
  type = "latex", 
  no.space = TRUE,
  label = "table_simple_implied",
  digits = 2,
  order = c("theta_1", "theta_10", "theta_0", "theta_01", "pi", "pi_0", "pi_1"),
  se = list(lm_model_3w_ar1_pi0$se_list, lm_model_3w_ar1$se_list, lm_model_3w_ar1_asymmetric$se_list), 
  keep.stat = c("n"), 
  covariate.labels = c("Entry rate", "Exit rate", "Misclassification rate", "Misclass. rate: non-employed", "Misclass. rate: employed"),
  add.lines = list(
    c("Misclassification", "No", "Symmetric", "Asymmetric"),
    c("LL", round(lm_model_3w_ar1_pi0$ll, 1), round(lm_model_3w_ar1$ll, 1), round(lm_model_3w_ar1_asymmetric$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)

output_file <- "./output/tables/SA/table_simple_implied.tex"
cat(table_simple_implied, file = output_file, sep = "\n")

#> Covariates that affect transitions ====

#>> Covariates (age and educ) & ME ----

df_estimate <- df_qlfs_age_educ

df_covariate_combos <- df_qlfs_age_educ %>% 
  select(age1, age2, educ1, educ2) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
  ) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates_age_educ <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates_age_educ[, k := NULL]

param_init <- data.frame(intercept_0 = -9.480411, age_0 = 0.3578881, age2_0 = -0.004708052, educ_0 = 0.007591136, intercept_1 = 0.3971249, age_1 = -0.06975788, age2_1 = 0.0002472373, educ_1 = -0.1404834, pi = 0.02375949)
param_init_transformed <- param_init
param_init_transformed$pi <- logit_transform(param_init$pi)

model_mle_3w_ar1_covariates_age_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ$estimate
model_mle_3w_ar1_covariates_age_educ$maximum

model_mle_3w_ar1_covariates_age_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ,
  start = model_mle_3w_ar1_covariates_age_educ$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ$estimate
model_mle_3w_ar1_covariates_age_educ$maximum

param_temp <- model_mle_3w_ar1_covariates_age_educ$estimate

df_transition_probs <- df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age1 + param_temp$age2_0*age1^2 + param_temp$educ_0*educ1),
    theta0_2 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age2 + param_temp$age2_0*age2^2 + param_temp$educ_0*educ2),
    theta1_1 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age1 + param_temp$age2_1*age1^2 + param_temp$educ_1*educ1),
    theta1_2 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age2 + param_temp$age2_1*age2^2 + param_temp$educ_1*educ2),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  )

mean_transition_rates <- df_transition_probs %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )

logit_inverse(model_mle_3w_ar1_covariates_age_educ$estimate$pi)

lm_model_mle_3w_ar1_covariates_age_educ <- create_stargazer_table_covariates(
  model_object = model_mle_3w_ar1_covariates_age_educ, 
  df = df_estimate, 
  df_transition_probs, 
  mean_transition_rates,
  formula1 = "~ age1 +  I(age1^2) + educ1", 
  formula2 = "~ age2 +  I(age2^2) + educ2"
)


df_educ_age_grid <- data.table::CJ(
  age = seq(18,55),
  educ = seq(0, 17)
) %>%
  mutate(
    theta_0 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age + param_temp$age2_0*age^2 + param_temp$educ_0*educ),
    theta_1 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age + param_temp$age2_1*age^2 + param_temp$educ_1*educ),
    mu = theta_0/(theta_1 + theta_0)
  )

df_educ_age_grid %>%
  mutate(educ = as.factor(educ)) %>%
  ggplot(aes(x = age, y = theta_0, color = educ, group = educ)) +
  geom_line()

df_educ_age_grid %>%
  mutate(educ = as.factor(educ)) %>%
  ggplot(aes(x = age, y = theta_1, color = educ, group = educ)) +
  geom_line()

df_educ_age_grid %>%
  mutate(educ = as.factor(educ)) %>%
  ggplot(aes(x = age, y = mu, color = educ, group = educ)) +
  geom_line()

#>> No ME ----

param_init <- data.frame(intercept_0 = -7.517197, age_0 = 0.263745, age2_0 = -0.003260282, educ_0 = 0.03267523, intercept_1 = 2.445478, age_1 = -0.1652242, age2_1 = 0.001691698, educ_1 = -0.1187962)
param_init_transformed <- param_init

model_mle_3w_ar1_covariates_age_educ_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_pi0,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_pi0$maximum

model_mle_3w_ar1_covariates_age_educ_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_pi0,
  start = model_mle_3w_ar1_covariates_age_educ_pi0$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_pi0$maximum


param_temp <- model_mle_3w_ar1_covariates_age_educ_pi0$estimate

df_transition_probs <- df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age1 + param_temp$age2_0*age1^2 + param_temp$educ_0*educ1),
    theta0_2 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age2 + param_temp$age2_0*age2^2 + param_temp$educ_0*educ2),
    theta1_1 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age1 + param_temp$age2_1*age1^2 + param_temp$educ_1*educ1),
    theta1_2 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age2 + param_temp$age2_1*age2^2 + param_temp$educ_1*educ2),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  ) 

mean_transition_rates <- df_transition_probs %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )

#>> Covariates (age, educ, female & race) & ME ----

df_estimate <- df_qlfs_age_educ_female_race

df_covariate_combos <- df_qlfs_age_educ_female_race %>% 
  select(age1, age2, educ1, educ2, female1, female2, race1, race2) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates_age_educ_female_race <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates_age_educ_female_race[, k := NULL]

param_init <- data.frame(intercept_0 = -9.200884, age_0 = 0.353941, age2_0 = -0.004574609, educ_0 = 0.02526308, female_0 = -0.75322, race2_0 = 0.1216257, race3_0 = -1.188118, race4_0 = -0.7020223, intercept_1 = 0.6947394, age_1 = -0.09688445, age2_1 = 0.000684121, educ_1 = -0.1141063, female_1 = 0.02715901, race2_1 = -0.2864189, race3_1 = -1.388325, race4_1 = -1.507239, pi = 0.0222393)
param_init_transformed <- param_init
param_init_transformed$pi <- logit_transform(param_init$pi)

model_mle_3w_ar1_covariates_age_educ_female_race <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race$estimate
model_mle_3w_ar1_covariates_age_educ_female_race$maximum

model_mle_3w_ar1_covariates_age_educ_female_race <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race,
  start = model_mle_3w_ar1_covariates_age_educ_female_race$estimate,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race$estimate
model_mle_3w_ar1_covariates_age_educ_female_race$maximum

param <- model_mle_3w_ar1_covariates_age_educ_female_race$estimate

df_transition_probs <- df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$female_0*female1 + param$race2_0*(race1 == 2) + param$race3_0*(race1 == 3) + param$race4_0*(race1 == 4)),
    theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$female_0*female2 + param$race2_0*(race2 == 2) + param$race3_0*(race2 == 3) + param$race4_0*(race2 == 4)),
    theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$female_1*female1 + param$race2_1*(race1 == 2) + param$race3_1*(race1 == 3) + param$race4_1*(race1 == 4)),
    theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$female_1*female2 + param$race2_1*(race2 == 2) + param$race3_1*(race2 == 3) + param$race4_1*(race2 == 4)),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  ) 

mean_transition_rates <- df_transition_probs %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )

logit_inverse(model_mle_3w_ar1_covariates_age_educ_female_race$estimate$pi)

lm_model_mle_3w_ar1_covariates_age_educ_female_race <- create_stargazer_table_covariates(
  model_object = model_mle_3w_ar1_covariates_age_educ_female_race, 
  df = df_estimate, 
  df_transition_probs, 
  mean_transition_rates,
  formula1 = "~ age1 +  I(age1^2) + educ1 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)", 
  formula2 = "~ age2 +  I(age2^2) + educ2 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)"
)


#>> No ME ----

param_init <- data.frame(intercept_0 = -7.46008, age_0 = 0.2713965, age2_0 = -0.0033241, educ_0 = 0.0419672, female_0 = -0.6312398, race2_0 = 0.2227918, race3_0 = -0.5986817, race4_0 = -0.07955744, intercept_1 = 2.522993, age_1 = -0.1810162, age2_1 = 0.001934264, educ_1 = -0.1005616, female_1 = 0.1527486, race2_1 = -0.1987458, race3_1 = -0.789203, race4_1 = -0.8971382)

param_init_transformed <- param_init

model_mle_3w_ar1_covariates_age_educ_female_race_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_pi0,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_pi0$maximum

model_mle_3w_ar1_covariates_age_educ_female_race_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_pi0,
  start = model_mle_3w_ar1_covariates_age_educ_female_race_pi0$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_pi0$maximum


param_temp <- model_mle_3w_ar1_covariates_age_educ_female_race_pi0$estimate

df_transition_probs <- df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age1 + param_temp$age2_0*age1^2 + param_temp$educ_0*educ1),
    theta0_2 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age2 + param_temp$age2_0*age2^2 + param_temp$educ_0*educ2),
    theta1_1 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age1 + param_temp$age2_1*age1^2 + param_temp$educ_1*educ1),
    theta1_2 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age2 + param_temp$age2_1*age2^2 + param_temp$educ_1*educ2),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  ) 

mean_transition_rates <- df_transition_probs %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )



# lm_model_mle_3w_ar1_covariates_age_educ_female_race_contract <- create_stargazer_table_covariates(
#   model_object = model_mle_3w_ar1_covariates_age_educ_female_race_contract, 
#   df = df_estimate, 
#   df_transition_probs, 
#   mean_transition_rates,
#   formula1 = "~ age1 +  I(age1^2) + educ1 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)", 
#   formula2 = "~ age2 +  I(age2^2) + educ2 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)"
# )

#>> Covariates (age, educ, female, race & contract) & ME ----

df_estimate <- df_qlfs_age_educ_female_race_contract

df_covariate_combos <- df_qlfs_age_educ_female_race_contract %>% 
  select(age1, age2, educ1, educ2, female1, female2, race1, race2, contracttype1, contracttype2, contracttype1_missing, contracttype2_missing) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates_age_educ_female_race_contract <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates_age_educ_female_race_contract[, k := NULL]

param_init <- data.frame(intercept_0 = -11.09761, age_0 = 0.4816279, age2_0 = -0.006445747, educ_0 = 0.009268962 , female_0 = -0.4586408, race2_0 = 0.006957785, race3_0 = -1.202247, race4_0 = -0.5641202, intercept_1 = 7.782345, age_1 = 0.2402219, age2_1 = -0.003885501, educ_1 = -0.08984728, female_1 = 0.05677661, race2_1 = -0.1199074 , race3_1 = -1.575087, race4_1 = -1.800444 , contract = 0.807705, contract_missing = 4.848624, pi = 0.01719847)
param_init_transformed <- param_init
param_init_transformed$pi <- logit_transform(param_init$pi)

model_mle_3w_ar1_covariates_age_educ_female_race_contract <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_contract,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_contract$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_contract$maximum

model_mle_3w_ar1_covariates_age_educ_female_race_contract <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_contract,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract,
  start = model_mle_3w_ar1_covariates_age_educ_female_race_contract$estimate,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_contract$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_contract$maximum

param <- model_mle_3w_ar1_covariates_age_educ_female_race_contract$estimate

df_transition_probs <- df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$female_0*female1 + param$race2_0*(race1 == 2) + param$race3_0*(race1 == 3) + param$race4_0*(race1 == 4)),
    theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$female_0*female2 + param$race2_0*(race2 == 2) + param$race3_0*(race2 == 3) + param$race4_0*(race2 == 4)),
    theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$female_1*female1 + param$race2_1*(race1 == 2) + param$race3_1*(race1 == 3) + param$race4_1*(race1 == 4) + param$contract*contracttype1 + param$contract*param$contract_missing*contracttype1_missing),
    theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$female_1*female2 + param$race2_1*(race2 == 2) + param$race3_1*(race2 == 3) + param$race4_1*(race2 == 4) + param$contract*contracttype2 + param$contract*param$contract_missing*contracttype2_missing),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  ) 

mean_transition_rates <- df_transition_probs %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )

logit_inverse(model_mle_3w_ar1_covariates_age_educ_female_race_contract$estimate$pi)


# lm_model_mle_3w_ar1_covariates_age_educ_female_race_contract <- create_stargazer_table_covariates(
#   model_object = model_mle_3w_ar1_covariates_age_educ_female_race_contract, 
#   df = df_estimate, 
#   df_transition_probs, 
#   mean_transition_rates,
#   formula1 = "~ age1 +  I(age1^2) + educ1 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)", 
#   formula2 = "~ age2 +  I(age2^2) + educ2 + female1 + I(race1 == 2) + I(race1 == 3) + I(race1 == 4)"
# )

#>> No ME ----

param_init <- data.frame(intercept_0 = -9.480411, age_0 = 0.3694859, age2_0 = -0.004521592, educ_0 = 0.08960009, female_0 = -1.162182 , race2_0 = 0.2146799, race3_0 = -0.0138424, race4_0 = 0.2491305,
                         intercept_1 = 0.3905569, age_1 = -0.1650995, age2_1 = 0.001647857, educ_1 = -0.1005749, female_1 = 0.9313917 , race2_1 = -0.1968393 , race3_1 = -0.06141175, race4_1 = -0.3832552, 
                         contract = 0.4761247, contract_missing = 5.966967)
param_init_transformed <- param_init

model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_contract_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract_pi0,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$maximum

model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates_age_educ_female_race_contract_pi0,
  grad = calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract_pi0,
  start = model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$estimate
model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$maximum


param_temp <- model_mle_3w_ar1_covariates_age_educ_female_race_contract_pi0$estimate

df_estimate %>% 
  mutate(
    theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$female_0*female1 + param$race2_0*(race1 == 2) + param$race3_0*(race1 == 3) + param$race4_0*(race1 == 4)),
    theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$female_0*female2 + param$race2_0*(race2 == 2) + param$race3_0*(race2 == 3) + param$race4_0*(race2 == 4)),
    theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$female_1*female1 + param$race2_1*(race1 == 2) + param$race3_1*(race1 == 3) + param$race4_1*(race1 == 4) + param$contract*contracttype1 + param$contract*param$contract_missing*contracttype1_missing),
    theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$female_1*female2 + param$race2_1*(race2 == 2) + param$race3_1*(race2 == 3) + param$race4_1*(race2 == 4) + param$contract*contracttype2 + param$contract*param$contract_missing*contracttype2_missing),
    mu_1 = theta0_1/(theta1_1 + theta0_1)
  ) %>% 
  summarise(
    theta0_1 = weighted.mean(theta0_1, weight),
    theta0_2 = weighted.mean(theta0_2, weight),
    theta1_1 = weighted.mean(theta1_1, weight),
    theta1_2 = weighted.mean(theta1_2, weight)
  )
#>> Creat Latex table ====

# Generate the stargazer table
table_covariates_implied <- stargazer::stargazer(
  lm_model_mle_3w_ar1_covariates_age_educ, lm_model_mle_3w_ar1_covariates_age_educ_female_race,  
  type = "latex", 
  no.space = TRUE,
  label = "table_covariates_implied",
  digits = 2,
  order = c("theta0_1", "theta0_2", "theta1_1", "theta1_2", "pi"),
  se = list(lm_model_mle_3w_ar1_covariates_age_educ$se_list, lm_model_mle_3w_ar1_covariates_age_educ_female_race$se_list), 
  keep.stat = c("n"), 
  covariate.labels = c("Entry rate: period 2", "Entry rate: period 3", "Exit rate: period 2", "Exit rate: period 3" ,"Misclassification rate"),
  add.lines = list(
    c("Covariates", "Educ + Age", "+ Race + Gender"),
    c("LL", round(lm_model_mle_3w_ar1_covariates_age_educ$ll, 1), round(lm_model_mle_3w_ar1_covariates_age_educ_female_race$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)

output_file <- "./output/tables/SA/table_covariates_implied.tex"
cat(table_covariates_implied, file = output_file, sep = "\n")

model_covariates_age_educ_unrestricted <- lm(data = df_estimate, y3 ~ y1 + y2 + y3 + educ1 + educ2 + educ3 + age1 + age2 + age3)
coefficients_covariates_age_educ_unrestricted <- unlist(model_mle_3w_ar1_covariates_age_educ$estimate)
std_errors_covariates_age_educ_unrestricted <- unlist(sqrt(diag(vcov(model_mle_3w_ar1_covariates_age_educ))))
obs_covariates_age_educ_unrestricted <- 3*nrow(df_estimate)
ll_covariates_age_educ_unrestricted <- model_mle_3w_ar1_covariates_age_educ$maximum
model_covariates_age_educ_unrestricted$coefficients <- coefficients_covariates_age_educ_unrestricted
model_covariates_age_educ_unrestricted$residuals <- rnorm(obs_covariates_age_educ_unrestricted)
model_covariates_age_educ_unrestricted$ll <- ll_covariates_age_educ_unrestricted
stargazer::stargazer(model_covariates_age_educ_unrestricted, type = "text", se = list(std_errors_covariates_age_educ_unrestricted), keep.stat = c("n"))

model_covariates_age_educ_female_race_unrestricted <- lm(data = df_estimate, y3 ~ y1 + y2 + educ1 + educ2 + educ3 + age1 + age2 + age3 + I(race1 == 2)+ I(race1 == 3)+ I(race1 == 4) + I(race2 == 2)+ I(race2 == 3)+ I(race2 == 4) + female1 + female2)
coefficients_covariates_age_educ_female_race_unrestricted <- unlist(model_mle_3w_ar1_covariates_age_educ_female_race$estimate)
std_errors_covariates_age_educ_female_race_unrestricted <- unlist(sqrt(diag(vcov(model_mle_3w_ar1_covariates_age_educ_female_race))))
obs_covariates_age_educ_female_race_unrestricted <- 3*nrow(df_estimate)
model_covariates_age_educ_female_race_unrestricted$coefficients <- coefficients_covariates_age_educ_female_race_unrestricted
model_covariates_age_educ_female_race_unrestricted$residuals <- rnorm(obs_covariates_age_educ_female_race_unrestricted)
model_covariates_age_educ_female_race_unrestricted$ll <- model_mle_3w_ar1_covariates_age_educ_female_race$maximum
stargazer::stargazer(model_covariates_age_educ_female_race_unrestricted, type = "text", se = list(std_errors_covariates_age_educ_female_race_unrestricted), keep.stat = c("n"))

table_covariates_unrestricted <- stargazer::stargazer(
  model_covariates_age_educ_unrestricted, model_covariates_age_educ_female_race_unrestricted,  
  type = "latex", 
  no.space = TRUE,
  label = "table_covariates_unrestricted",
  digits = 2,
  order = c("intercept_0", "age_0", "age2_0", "educ_0", "female_0", "race2_0", "race3_0", "race4_0", "intercept_1", "age_1", "age2_1", "educ_1", "female_1", "race2_1", "race3_1", "race4_1", "pi"),
  se = list(model_covariates_age_educ_unrestricted$se_list, model_covariates_age_educ_female_race_unrestricted$se_list), 
  keep.stat = c("n"), 
  covariate.labels = c("Constant (entry)", "Age (entry)", "Age squared (entry)", "Educ (entry)", "Race = coloured (entry)", "Race = Indian/Asian (entry)", "Race = white (entry)" , "Female (entry)" , "Constant (exit)", "Age (exit)", "Age squared (exit)", "Educ (exit)", "Race = coloured (exit)", "Race = Indian/Asian (exit)", "Race = white (exit)" , "Female (exit)", "Misclassification rate"),
  add.lines = list(
    c("Covariates", "Educ + Age", "+ Race + Gender"),
    c("LL", round(model_covariates_age_educ_unrestricted$ll, 1), round(model_covariates_age_educ_female_race_unrestricted$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)

output_file <- "./output/tables/SA/table_covariates_unrestricted.tex"
cat(table_covariates_unrestricted, file = output_file, sep = "\n")


#>> Duration dependence (timegap and tenure) ----

df_estimate <- df_qlfs_tenure_timegap

df_covariate_combos <- df_estimate %>%
  select(y1, y2, y3, timegap1, timegap2, tenure1, tenure2) %>%
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
df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates[, k := NULL]

# param_init <- data.frame(intercept_0 = -2.524934, timegap = -0.1300759, intercept_1 = 1.02172, tenure = -0.6729012, pi = 0.02924524)
param_init <- data.frame(intercept_0 = -3.47, timegap = 0, intercept_1 = -3.47, tenure = 0, pi = 0.02924524)
param_init_transformed <- param_init
param_init_transformed$pi <- logit_transform(param_init$pi)


model_mle_3w_ar1_covariates2 <- maxLik::maxLik(
  calc_lli_3waves_ar1_covariates_duration,
  grad = calc_mle_derivatives_3waves_ar1_covariates_duration,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_covariates2$estimate
model_mle_3w_ar1_covariates2$maximum

model_mle_3w_ar1_covariates2 <- maxLik::maxLik(
  calc_mle_3waves_ar1_covariates2,
  grad = calc_mle_derivatives_3waves_ar1_covariates2,
  start = model_mle_3w_ar1_covariates2$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
# 
# model_mle_3w_ar1_covariates2$estimate
# model_mle_3w_ar1_covariates2$maximum
# 
# param_temp <- model_mle_3w_ar1_covariates2$estimate
# 
# df_estimate %>% 
#   mutate(
#     theta0_1 = logit_inverse(param$intercept_0 + param$timegap*log(timegap1 + 1.5)),
#     theta0_2 = logit_inverse(param$intercept_0 + param$timegap*log(timegap2 + 1.5)),
#     theta1_1 = logit_inverse(param$intercept_1 + param$tenure*log(tenure1 + 1.5)),
#     theta1_2 = logit_inverse(param$intercept_1 + param$tenure*log(tenure2 + 1.5)),
#     mu_1 = theta0_1/(theta1_1 + theta0_1)
#   ) %>% 
#   summarise(
#     theta0_1 = weighted.mean(theta0_1, weight),
#     theta0_2 = weighted.mean(theta0_2, weight),
#     theta1_1 = weighted.mean(theta1_1, weight),
#     theta1_2 = weighted.mean(theta1_2, weight),
#     mu_1 = weighted.mean(mu_1, weight)
#   )
# 
# 
# 
# 
# 
# #>> Covariates (age, educ, timegap and tenure) ----
# 
# df_covariate_combos <- df_estimate2 %>% 
#   select(y1, y2, y3, age1, age2, educ1, educ2, timegap1, timegap2, tenure1, tenure2) %>% 
#   unique() %>% 
#   na.omit %>% 
#   data.table()
# 
# df_template_employment <- data.table::CJ(
#   y1_star = c(0, 1), 
#   y2_star = c(0, 1), 
#   y3_star = c(0, 1)
# ) 
# 
# # Add a temporary key column without overwriting any existing data
# df_covariate_combos[, k := 1]
# df_template_employment[, k := 1]
# 
# # Perform the join while ensuring no columns are omitted
# df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]
# 
# # Remove the temporary key column `k` from the result
# df_template_covariates[, k := NULL]
# 
# param_init <- data.frame(intercept_0 = -2.524934, age_0 = 0.06983356, age2_0 = -0.001228011, educ_0 = -0.01023191, timegap = -0.1300759, intercept_1 = 1.02172, age_1 = -0.0857186, age2_1 =  0.001056658, educ_1 = -0.07880339, tenure_1 = -0.6729012, pi = 0.02924524)
# param_init_transformed <- param_init
# param_init_transformed$pi <- logit_transform(param_init$pi)
# 
# 
# model_mle_3w_ar1_covariates2 <- maxLik::maxLik(
#   calc_mle_3waves_ar1_covariates2,
#   grad = calc_mle_derivatives_3waves_ar1_covariates2,
#   start = param_init_transformed,
#   method = "BFGS",
#   reltol = 0,
#   gradtol = 0
# )
# 
# model_mle_3w_ar1_covariates2$estimate
# model_mle_3w_ar1_covariates2$maximum
# 
# # model_mle_3w_ar1_covariates1 <- maxLik::maxBHHH(
# #   calc_mle_3waves_ar1_covariates1,
# #   grad = calc_lli_derivatives_3waves_ar1_covariates1,
# #   start = param_init_transformed
# # )
# 
# model_mle_3w_ar1_covariates2 <- maxLik::maxLik(
#   calc_mle_3waves_ar1_covariates2,
#   grad = calc_mle_derivatives_3waves_ar1_covariates2,
#   start = model_mle_3w_ar1_covariates2$estimate,
#   method = "NR",
#   reltol = 0,
#   gradtol = 0
# )
# 
# model_mle_3w_ar1_covariates2$estimate
# model_mle_3w_ar1_covariates2$maximum
# 
# param_temp <- model_mle_3w_ar1_covariates2$estimate
# 
# df_estimate2 %>% 
#   mutate(
#     theta0_1 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age1 + param_temp$age2_0*age1^2 + param_temp$educ_0*educ1),
#     theta0_2 = logit_inverse(param_temp$intercept_0 + param_temp$age_0*age2 + param_temp$age2_0*age2^2 + param_temp$educ_0*educ2),
#     theta1_1 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age1 + param_temp$age2_1*age1^2 + param_temp$educ_1*educ1),
#     theta1_2 = logit_inverse(param_temp$intercept_1 + param_temp$age_1*age2 + param_temp$age2_1*age2^2 + param_temp$educ_1*educ2),
#     mu_1 = theta0_1/(theta1_1 + theta0_1)
#   ) %>% 
#   summarise(
#     theta0_1 = weighted.mean(theta0_1, weight),
#     theta0_2 = weighted.mean(theta0_2, weight),
#     theta1_1 = weighted.mean(theta1_1, weight),
#     theta1_2 = weighted.mean(theta1_2, weight)
#   )
# 
# 
# model_mle_3w_ar1_covariates2 <- maxLik::maxLik(
#   calc_mle_3waves_ar1_covariates2,
#   grad = calc_mle_derivatives_3waves_ar1_covariates2,
#   start = param_init_transformed,
#   method = "NR"
# )

#> Symptoms of mislcassification ----

df_estimate <- df_qlfs_age_educ %>% 
  mutate(age_inconsistent = if_else(
    (abs(age1 - age2) > 1) | (abs(age2 - age3) > 1), 1L, 0L
  )) %>% 
  mutate(educ_inconsistent = if_else(
    (educ2 - educ1 > 1) | (educ2 - educ1 < 0) | (educ3 - educ2 < 0) | (educ3 - educ2 < 0), 1L, 0L
  )) 


#> Symptoms: age inconsistent ----

df_covariate_combos <- df_estimate %>% 
  select(age_inconsistent) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates[, k := NULL]

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, intercept_pi = -7.907681, age_pi = 0.304568)
param_init_transformed <- param_init
param_init_transformed$theta_0 <- logit_transform(param_init$theta_0)
param_init_transformed$theta_1 <- logit_transform(param_init$theta_1)

model_mle_3w_ar1_misclassifiction_symptoms_age <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_age,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_age,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_age$estimate
model_mle_3w_ar1_misclassifiction_symptoms_age$maximum

model_mle_3w_ar1_misclassifiction_symptoms_age <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_age,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_age,
  start = model_mle_3w_ar1_misclassifiction_symptoms_age$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_age$estimate
model_mle_3w_ar1_misclassifiction_symptoms_age$maximum


#> Symptoms: educ inconsistent ----

df_covariate_combos <- df_estimate %>% 
  select(educ_inconsistent) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates[, k := NULL]

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, intercept_pi = -7.907681, educ_pi = 0.304568)
param_init_transformed <- param_init
param_init_transformed$theta_0 <- logit_transform(param_init$theta_0)
param_init_transformed$theta_1 <- logit_transform(param_init$theta_1)

model_mle_3w_ar1_misclassifiction_symptoms_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_educ,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_educ,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_educ$estimate
model_mle_3w_ar1_misclassifiction_symptoms_educ$maximum

model_mle_3w_ar1_misclassifiction_symptoms_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_educ,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_educ,
  start = model_mle_3w_ar1_misclassifiction_symptoms_educ$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_educ$estimate
model_mle_3w_ar1_misclassifiction_symptoms_educ$maximum



#> Symptoms: age & educ inconsistent ----

df_covariate_combos <- df_estimate %>% 
  select(age_inconsistent, educ_inconsistent) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates[, k := NULL]

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, intercept_pi = -7.907681, age_pi = 0.304568, educ_pi = 0.2)
param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, intercept_pi = -3.476099, age_pi = 0, educ_pi = 0)
param_init_transformed <- param_init
param_init_transformed$theta_0 <- logit_transform(param_init$theta_0)
param_init_transformed$theta_1 <- logit_transform(param_init$theta_1)

model_mle_3w_ar1_misclassifiction_symptoms <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms$estimate
model_mle_3w_ar1_misclassifiction_symptoms$maximum

model_mle_3w_ar1_misclassifiction_symptoms <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms,
  start = model_mle_3w_ar1_misclassifiction_symptoms$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms$estimate
model_mle_3w_ar1_misclassifiction_symptoms$maximum


# df_estimate %>% 
#   mutate(
#     pi = logit_inverse(param_temp$intercept_pi + param_temp$age_pi*age_inconsistent + param_temp$educ_pi*educ_inconsistent)
#   ) %>% 
#   summarise(
#     pi = weighted.mean(pi, weight)
#   )
# 
# df_estimate %>% 
#   mutate(
#     pi = logit_inverse(param_temp$intercept_pi + param_temp$age_pi*age_inconsistent + param_temp$educ_pi*educ_inconsistent)
#   ) %>% 
#   group_by(educ_inconsistent, age_inconsistent) %>% 
#   summarise(
#     pi = mean(pi)
#   )




#> Symptoms: age & educ + interaction inconsistent ----

df_covariate_combos <- df_estimate %>% 
  select(age_inconsistent, educ_inconsistent) %>% 
  unique() %>% 
  na.omit %>% 
  data.table()

df_template_employment <- data.table::CJ(
  y1 = c(0, 1), 
  y1_star = c(0, 1), 
  y2 = c(0, 1), 
  y2_star = c(0, 1), 
  y3 = c(0, 1), 
  y3_star = c(0, 1)
) 

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_covariates <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_covariates[, k := NULL]

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, intercept_pi = -7.907681, age_pi = 0.304568, educ_pi = 0.2, age_educ_pi = 0)
# param_init <- data.frame(theta_0 = 0.03, theta_1 = 0.03, intercept_pi = -3.476099, age_pi = 0, educ_pi = 0)
param_init_transformed <- param_init
param_init_transformed$theta_0 <- logit_transform(param_init$theta_0)
param_init_transformed$theta_1 <- logit_transform(param_init$theta_1)

model_mle_3w_ar1_misclassifiction_symptoms_age_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_age_educ,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_age_educ,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_age_educ$estimate
model_mle_3w_ar1_misclassifiction_symptoms_age_educ$maximum

model_mle_3w_ar1_misclassifiction_symptoms_age_educ <- maxLik::maxLik(
  calc_mle_3waves_ar1_misclassifiction_symptoms_age_educ,
  grad = calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms_age_educ,
  start = model_mle_3w_ar1_misclassifiction_symptoms_age_educ$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)

model_mle_3w_ar1_misclassifiction_symptoms_age_educ$estimate
model_mle_3w_ar1_misclassifiction_symptoms_age_educ$maximum

#>> Creat Latex table ====

lm_model_symptoms_age <- create_stargazer_table_symptoms(model_object = model_mle_3w_ar1_misclassifiction_symptoms_age, df = df_estimate, formula = "~ age_inconsistent") 
lm_model_symptoms_educ <- create_stargazer_table_symptoms(model_object = model_mle_3w_ar1_misclassifiction_symptoms_educ, df = df_estimate, formula = "~ educ_inconsistent") 
lm_model_symptoms <- create_stargazer_table_symptoms(model_object = model_mle_3w_ar1_misclassifiction_symptoms, df = df_estimate, formula = "~ age_inconsistent + educ_inconsistent") 
lm_model_symptoms_age_educ <- create_stargazer_table_symptoms(model_object = model_mle_3w_ar1_misclassifiction_symptoms_age_educ, df = df_estimate, formula = "~ age_inconsistent + educ_inconsistent + age_inconsistent*educ_inconsistent") 


# Generate the stargazer table
table_symptoms_implied <- stargazer::stargazer(
  lm_model_symptoms_age, lm_model_symptoms_educ, lm_model_symptoms, 
  type = "latex",
  no.space = TRUE,
  label = "table_symptoms_implied",
  digits = 2,
  order = c("theta_0", "theta_1", "Constant", "age_inconsistent", "educ_inconsistent"),
  se = list(lm_model_symptoms_age$se_list, lm_model_symptoms_educ$se_list, lm_model_symptoms$se_list),
  keep.stat = c("n"),
  covariate.labels = c("Entry rate", "Exit rate", "Miscl. rate: no incon.", "Miscl. rate: age incon.", "Miscl. rate: educ incon."),
  add.lines = list(
    c("Inconsistencies", "Age", "Educ", "Age + Educ"),
    c("LL", round(lm_model_symptoms_age$ll, 1), round(lm_model_symptoms_educ$ll, 1), round(lm_model_symptoms$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)


output_file <- "./output/tables/SA/table_symptoms_implied.tex"
cat(table_symptoms_implied, file = output_file, sep = "\n")

#>> FMM (two groups) ----


# param_init <- data.frame(theta0_1 = 0.1, theta1_1 = 0.1, p_1 = 0.5, theta0_2 = 0.5, theta1_2 = 0.5, pi = 0.01)
param_init <- data.frame(theta0_1 = 0.2946616, theta1_1 = 0.02416365, p_1 = 0.277803, theta0_2 = 0.02341918, theta1_2 = 0.04854722, pi = 0.02714389)


df_estimate <- df_qlfs
param_init <- data.frame(theta0_1 = 0.085, theta1_1 = 0.085, p_1 = 0.1, theta0_2 = 0.03, theta1_2 = 0.03, pi = 0.03)
param_init_transformed <- logit_transform(param_init)
model_mle_3w_ar1_fmm2 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_lli_derivatives_3waves_ar1_fmm2,
  start = param_init_transformed,
  hess = calc_hessian_3waves_ar1_fmm2,
  method = "BFGS"
)

model_mle_3w_ar1_fmm2$estimate
logit_inverse(model_mle_3w_ar1_fmm2$estimate)
model_mle_3w_ar1_fmm2$maximum
calc_mle_derivatives_3waves_ar1_fmm2(model_mle_3w_ar1_fmm2$estimate)
diag(vcov(model_mle_3w_ar1_fmm2))

model_mle_3w_ar1_fmm2 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_lli_derivatives_3waves_ar1_fmm2,
  start = param_init_transformed,
  hess = calc_hessian_3waves_ar1_fmm2,
  method = "NR"
)

model_mle_3w_ar1_fmm2 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_lli_derivatives_3waves_ar1_fmm2,
  start = model_mle_3w_ar1_fmm2$estimate,
  hess = calc_hessian_3waves_ar1_fmm2,
  method = "BFGS"
)

model_mle_3w_ar1_fmm2 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_lli_derivatives_3waves_ar1_fmm2,
  start = model_mle_3w_ar1_fmm2$estimate,
  hess = calc_hessian_3waves_ar1_fmm2,
  method = "CG"
)

model_mle_3w_ar1_fmm2$estimate
logit_inverse(model_mle_3w_ar1_fmm2$estimate)
model_mle_3w_ar1_fmm2$maximum
calc_mle_derivatives_3waves_ar1_fmm2(model_mle_3w_ar1_fmm2$estimate)
diag(vcov(model_mle_3w_ar1_fmm2))
sqrt(((exp(model_mle_3w_ar1_fmm2$estimate) / ((1 + exp(model_mle_3w_ar1_fmm2$estimate))^2))^2) * 
       diag(vcov(model_mle_3w_ar1_fmm2)))

model_mle_3w_ar1_fmm2 <- maxLik::maxNR(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_mle_derivatives_3waves_ar1_fmm2,
  hess = calc_hessian_3waves_ar1_fmm2,
  start = model_mle_3w_ar1_fmm2$estimate,
)

model_mle_3w_ar1_fmm2$estimate
logit_inverse(model_mle_3w_ar1_fmm2$estimate)
model_mle_3w_ar1_fmm2$maximum
hessian_matrix <- model_mle_3w_ar1_fmm2$hessian
var_cov_matrix <- solve(-hessian_matrix) # Invert Hessian for variance-covariance
print(var_cov_matrix)
sqrt(((exp(model_mle_3w_ar1_fmm2$estimate) / ((1 + exp(model_mle_3w_ar1_fmm2$estimate))^2))^2) * 
       diag(var_cov_matrix))


vcov(model_mle_3w_ar1_fmm2)

model_mle_3w_ar1_fmm2 <- maxLik::maxBFGSR(
  calc_mle_3waves_ar1_fmm2,
  grad = calc_mle_derivatives_3waves_ar1_fmm2,
  start = model_mle_3w_ar1_fmm2$estimate,
  finalHessian = TRUE
)

model_mle_3w_ar1_fmm2$estimate
logit_inverse(model_mle_3w_ar1_fmm2$estimate)
model_mle_3w_ar1_fmm2$maximum
calc_mle_derivatives_3waves_ar1_fmm2(model_mle_3w_ar1_fmm2$estimate)
# Extract Hessian and compute variance-covariance matrix
hessian_matrix <- model_mle_3w_ar1_fmm2$hessian
var_cov_matrix <- solve(hessian_matrix) # Invert Hessian for variance-covariance
print(var_cov_matrix)

# Optionally, compute standard errors (square roots of diagonal elements)
std_errors <- sqrt(diag(var_cov_matrix))

model_se <- sqrt(((exp(model_mle_3w_ar1_fmm2$estimate) / ((1 + exp(model_mle_3w_ar1_fmm2$estimate))^2))^2) * 
                   diag(var_cov_matrix))

print(std_errors)
diag(vcov(model_mle_3w_ar1_fmm2))

calc_mle_3waves_ar1_fmm2_help <- function(param_transformed) {
  param$theta0_1 <- param_transformed[1]
  param$theta1_1 <- param_transformed[2]
  param$p_1 <- param_transformed[3]
  param$theta0_2 <- param_transformed[4]
  param$theta1_2 <- param_transformed[5]
  param$pi <- param_transformed[6]
  calc_mle_3waves_ar1_fmm2(param)
}
param_transformed1 <- as.numeric(model_mle_3w_ar1_fmm2$estimate)
library(numDeriv)
  # example parameters
analytic_grad <- calc_mle_derivatives_3waves_ar1_fmm2(model_mle_3w_ar1_fmm2$estimate)
numeric_grad <- grad(calc_mle_3waves_ar1_fmm2_help, param_transformed1)
print(cbind(analytic_grad, numeric_grad))

# model_mle_3w_ar1_fmm2$estimate
# logit_inverse(model_mle_3w_ar1_fmm2$estimate)
# model_mle_3w_ar1_fmm2$maximum
# model_mle_3w_ar1_pi0_se <- sqrt(((exp(model_mle_3w_ar1_pi0$estimate)/((1 + exp(model_mle_3w_ar1_pi0$estimate))^2))^2)*diag(vcov(model_mle_3w_ar1_pi0)))

# print(logit_inverse(model_mle_3w_ar1_pi0$estimate))
# print(model_mle_3w_ar1_pi0_se)
# print(logit_inverse(model_mle_3w_ar1_pi0$estimate)/model_mle_3w_ar1_pi0_se)


#>> FMM (three groups) ----


df_estimate <- df_qlfs

param_init <- data.frame(theta0_1 = 0.1952656, theta1_1 = 0.1040511, p_1 = 0.1058684, theta0_2 = 0.04340665, theta1_2 = 0.03654352, p_2 = 0.1249773, theta0_3 = 0.01580732, theta1_3 = 0.01799972, pi = 0.02826307)
param_init_transformed <- logit_transform(param_init)
model_mle_3w_ar1_fmm3 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm3,
  grad = calc_mle_derivatives_3waves_ar1_fmm3,
  start = param_init_transformed,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_fmm3$estimate
logit_inverse(model_mle_3w_ar1_fmm3$estimate)
model_mle_3w_ar1_fmm3$maximum

model_mle_3w_ar1_fmm3 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm3,
  grad = calc_mle_derivatives_3waves_ar1_fmm3,
  start = model_mle_3w_ar1_fmm3$estimate,
  method = "BFGS",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_fmm3$estimate
logit_inverse(model_mle_3w_ar1_fmm3$estimate)
model_mle_3w_ar1_fmm3$maximum

model_mle_3w_ar1_fmm3 <- maxLik::maxLik(
  calc_mle_3waves_ar1_fmm3,
  grad = calc_mle_derivatives_3waves_ar1_fmm3,
  start = model_mle_3w_ar1_fmm3$estimate,
  method = "NR",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_fmm3$estimate
logit_inverse(model_mle_3w_ar1_fmm3$estimate)
model_mle_3w_ar1_fmm3$maximum

model_mle_3w_ar1_fmm3$estimate
logit_inverse(model_mle_3w_ar1_fmm3$estimate)
model_mle_3w_ar1_fmm3$maximum
calc_mle_derivatives_3waves_ar1_fmm3(model_mle_3w_ar1_fmm3$estimate)
diag(vcov(model_mle_3w_ar1_fmm3))

calc_mle_3waves_ar1_fmm3_help <- function(param_transformed) {
  param$theta0_1 <- param_transformed[1]
  param$theta1_1 <- param_transformed[2]
  param$p_1 <- param_transformed[3]
  param$theta0_2 <- param_transformed[4]
  param$theta1_2 <- param_transformed[5]
  param$p_2 <- param_transformed[6]
  param$theta0_3 <- param_transformed[7]
  param$theta1_3 <- param_transformed[8]
  param$pi <- param_transformed[9]
  calc_mle_3waves_ar1_fmm3(param)
}
param_transformed1 <- as.numeric(param_init_transformed)
library(numDeriv)
# example parameters
analytic_grad <- calc_mle_derivatives_3waves_ar1_fmm3(param_init_transformed)
numeric_grad <- grad(calc_mle_3waves_ar1_fmm3_help, param_transformed1)
print(cbind(analytic_grad, numeric_grad))


#>> Creat Latex table ====

model_mle_3w_ar1_fmm1 <- model_mle_3w_ar1
names(model_mle_3w_ar1_fmm1$estimate) <- c("theta0_1", "theta1_1", "pi")

# FMM: 1 group
lm_model_3w_ar1_fmm1 <- create_stargazer_table(model_object = model_mle_3w_ar1_fmm1, df = df_qlfs, formula = "y3 ~ y1 + y2 + age1") 

# FMM: 2 groups
lm_model_3w_ar1_fmm2 <- create_stargazer_table(model_object = model_mle_3w_ar1_fmm2, df = df_qlfs, formula = "y3 ~ y1 + y2 + age1 + age2 + age3 + educ1") 

# FMM: 3 groups
lm_model_3w_ar1_fmm3 <- create_stargazer_table(model_object = model_mle_3w_ar1_fmm3, df = df_qlfs, formula = "y3 ~ y1 + y2 + age1 + age2 + age3 + educ1 + educ2 + educ3 + female1") 

# Generate the stargazer table
table_fmm_implied <- stargazer::stargazer(
  lm_model_3w_ar1_fmm1, lm_model_3w_ar1_fmm2, lm_model_3w_ar1_fmm3, 
  type = "text", 
  no.space = TRUE,
  label = "table_fmm_implied",
  digits = 2,
  order = c("theta0_1", "theta1_1", "theta0_2", "theta1_2", "theta0_3", "theta1_3","pi", "p_1", "p_2"),
  se = list(lm_model_3w_ar1_fmm1$se_list, lm_model_3w_ar1_fmm2$se_list, lm_model_3w_ar1_fmm3$se_list), 
  keep.stat = c("n"), 
  covariate.labels = c("Entry rate 1", "Exit rate 1", "Entry rate 2", "Exit rate 2", "Entry rate 3", "Exit rate 3", "Misclassification rate", "Prob 1", "Prob 2"),
  add.lines = list(
    c("Misclassification", "No", "Symmetric", "Asymmetric"),
    c("LL", round(lm_model_3w_ar1_fmm1$ll, 1), round(lm_model_3w_ar1_fmm2$ll, 1), round(lm_model_3w_ar1_fmm3$ll, 1))
  ),
  dep.var.labels.include = FALSE,
  dep.var.caption = "",
  title = ""
)

output_file <- "./output/tables/SA/table_simple_implied.tex"
cat(table_simple_implied, file = output_file, sep = "\n")

