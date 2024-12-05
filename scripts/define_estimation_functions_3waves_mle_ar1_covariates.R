# DEFINE FUNCTIONS ====

#> General functions that restrict parameter values to unit interval

logit_transform <- function(param0) {
  log(param0/(1 - param0))
}

logit_inverse <- function(param_input0) {
  1/(1 + exp(-param_input0))
}

# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH MISCLASSIFICATION ERROR AND COVARIATE DEPENDENT TRANSITION RATES (educ + age) ====

calc_lli_3waves_ar1_covariates1 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
      ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, age1, age2, educ1, educ2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'age1', 'age2', 'educ1', 'educ2')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)

}

calc_lli_derivatives_3waves_ar1_covariates1 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    mutate(
      d1_star_theta0_1 = case_when(
        y1_star == 0 ~ -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = case_when(
        y1_star == 0 ~ theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta0_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, age1, age2, educ1, educ2) %>% 
    summarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1),
      joint_d_theta0_2 = sum(joint_d_theta0_2),
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_theta1_2 = sum(joint_d_theta1_2), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      theta0_1 = mean(theta0_1),
      theta0_2 = mean(theta0_2),
      theta1_1 = mean(theta1_1),
      theta1_2 = mean(theta1_2),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta0_1 = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2 = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1 = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2 = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    ) %>%
    mutate(
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0 = joint_d_theta0_1*age1 + joint_d_theta0_2*age2,
      joint_d_age2_0 = joint_d_theta0_1*(age1^2) + joint_d_theta0_2*(age2^2),
      joint_d_educ_0 = joint_d_theta0_1*educ1 + joint_d_theta0_2*educ2,
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1 = joint_d_theta1_1*age1 + joint_d_theta1_2*age2,
      joint_d_age2_1 = joint_d_theta1_1*(age1^2) + joint_d_theta1_2*(age2^2),
      joint_d_educ_1 = joint_d_theta1_1*educ1 + joint_d_theta1_2*educ2,
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'age1', 'age2', 'educ1', 'educ2')) %>% 
    mutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_age_0 = weight*joint_d_age_0/joint_p,
      lgi_age2_0 = weight*joint_d_age2_0/joint_p,
      lgi_educ_0 = weight*joint_d_educ_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_age_1 = weight*joint_d_age_1/joint_p,
      lgi_age2_1 = weight*joint_d_age2_1/joint_p,
      lgi_educ_1 = weight*joint_d_educ_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
      ) %>%
    select(lgi_intercept_0, lgi_age_0, lgi_age2_0, lgi_educ_0, lgi_intercept_1, lgi_age_1, lgi_age2_1, lgi_educ_1, lgi_pi)
  
  if(pi0) df_gi <- df_gi %>% select(-lgi_pi)
  return(df_gi)
  
}
  
calc_mle_3waves_ar1_covariates1 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates1(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates1 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates1(param_transformed))
}

calc_mle_3waves_ar1_covariates1_pi0 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates1(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates1_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates1(param_transformed, pi0 = TRUE))
}


# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH MISCLASSIFICATION ERROR AND COVARIATE DEPENDENT TRANSITION RATES (timegap + tenure) ====


calc_lli_3waves_ar1_covariates2 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, timegap1, timegap2, tenure1, tenure2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'tenure1', 'tenure2')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}

calc_lli_derivatives_3waves_ar1_covariates2 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    mutate(
      d1_star_theta0_1 = case_when(
        y1_star == 0 ~ -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = case_when(
        y1_star == 0 ~ theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta0_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, timegap1, timegap2, tenure1, tenure2) %>% 
    summarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1),
      joint_d_theta0_2 = sum(joint_d_theta0_2),
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_theta1_2 = sum(joint_d_theta1_2), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      theta0_1 = mean(theta0_1),
      theta0_2 = mean(theta0_2),
      theta1_1 = mean(theta1_1),
      theta1_2 = mean(theta1_2),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta0_1 = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2 = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1 = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2 = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    ) %>%
    mutate(
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_timegap_0 = joint_d_theta0_1*log(timegap1 + 1.5) + joint_d_theta0_2*log(timegap2 + 1.5),
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_tenure_1 = joint_d_theta1_1*log(tenure1 + 1.5) + joint_d_theta1_2*log(tenure2 + 1.5)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'tenure1', 'tenure2')) %>% 
    mutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_timegap_0 = weight*joint_d_timegap_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_tenure_1 = weight*joint_d_tenure_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
    ) %>%
    select(lgi_intercept_0, lgi_timegap_0, lgi_intercept_1, lgi_tenure_1, lgi_pi)
}

calc_mle_3waves_ar1_covariates2 <- function(param_transformed) {
  
  sum(calc_lli_3waves_ar1_covariates2(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates2 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates2(param_transformed))
}

calc_mle_3waves_ar1_covariates2_pi0 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates2(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates2_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates2(param_transformed, pi0 = TRUE))
}


# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH MISCLASSIFICATION ERROR AND COVARIATE DEPENDENT TRANSITION RATES (educ + age + timegap + tenure) ====


calc_lli_3waves_ar1_covariates3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, age1, age2, educ1, educ2, timegap1, timegap2, tenure1, tenure2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'age1', 'age2', 'educ1', 'educ2', 'timegap1', 'timegap2', 'tenure1', 'tenure2')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}

calc_lli_derivatives_3waves_ar1_covariates3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - theta0_1,
        y1_star == 0 & y2_star == 1 ~ theta0_1,
        y1_star == 1 & y2_star == 0 ~ theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - theta0_2,
        y2_star == 0 & y3_star == 1 ~ theta0_2,
        y2_star == 1 & y3_star == 0 ~ theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    mutate(
      d1_star_theta0_1 = case_when(
        y1_star == 0 ~ -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = case_when(
        y1_star == 0 ~ theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1 ~ -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta0_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_2 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, age1, age2, educ1, educ2, timegap1, timegap2, tenure1, tenure2) %>% 
    summarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1),
      joint_d_theta0_2 = sum(joint_d_theta0_2),
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_theta1_2 = sum(joint_d_theta1_2), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      theta0_1 = mean(theta0_1),
      theta0_2 = mean(theta0_2),
      theta1_1 = mean(theta1_1),
      theta1_2 = mean(theta1_2),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta0_1 = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2 = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1 = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2 = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    ) %>%
    mutate(
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0 = joint_d_theta0_1*age1 + joint_d_theta0_2*age2,
      joint_d_age2_0 = joint_d_theta0_1*(age1^2) + joint_d_theta0_2*(age2^2),
      joint_d_educ_0 = joint_d_theta0_1*educ1 + joint_d_theta0_2*educ2,
      joint_d_timegap_0 = joint_d_theta0_1*log(timegap1 + 1.5) + joint_d_theta0_2*log(timegap2 + 1.5),
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1 = joint_d_theta1_1*age1 + joint_d_theta1_2*age2,
      joint_d_age2_1 = joint_d_theta1_1*(age1^2) + joint_d_theta1_2*(age2^2),
      joint_d_educ_1 = joint_d_theta1_1*educ1 + joint_d_theta1_2*educ2,
      joint_d_tenure_1 = joint_d_theta1_1*log(tenure1 + 1.5) + joint_d_theta1_2*log(tenure2 + 1.5)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'age1', 'age2', 'educ1', 'educ2', 'timegap1', 'timegap2', 'tenure1', 'tenure2')) %>% 
    mutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_age_0 = weight*joint_d_age_0/joint_p,
      lgi_age2_0 = weight*joint_d_age2_0/joint_p,
      lgi_educ_0 = weight*joint_d_educ_0/joint_p,
      lgi_timegap_0 = weight*joint_d_timegap_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_age_1 = weight*joint_d_age_1/joint_p,
      lgi_age2_1 = weight*joint_d_age2_1/joint_p,
      lgi_educ_1 = weight*joint_d_educ_1/joint_p,
      lgi_tenure_1 = weight*joint_d_tenure_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
    ) %>%
    select(lgi_intercept_0, lgi_age_0, lgi_age2_0, lgi_educ_0, lgi_timegap_0, lgi_intercept_1, lgi_age_1, lgi_age2_1, lgi_educ_1, lgi_tenure_1, lgi_pi)
}

calc_mle_3waves_ar1_covariates3 <- function(param_transformed) {
  sum(calc_lli_3waves_ar3_covariates3(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates3 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates3(param_transformed))
}
