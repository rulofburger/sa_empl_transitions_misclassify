# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH COVARIATE DEPENDENT MISCLASSIFICATION ERROR ====

calc_lli_3waves_ar1_misclassifiction_symptoms <- function(param_transformed) {
  
  param <- param_transformed
  param$theta_0 <- logit_inverse(param_transformed$theta_0)
  param$theta_1 <- logit_inverse(param_transformed$theta_1)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      pi = logit_inverse(param$intercept_pi + param$age_pi*age_inconsistent + param$educ_pi*educ_inconsistent),
      mu = param$theta_0/(param$theta_0 + param$theta_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y1_star == 0 & y2_star == 1 ~ param$theta_0,
        y1_star == 1 & y2_star == 0 ~ param$theta_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
        y2_star == 0 & y3_star == 1 ~ param$theta_0,
        y2_star == 1 & y3_star == 0 ~ param$theta_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
      ),
      p1 = if_else(y1 == y1_star, 1 - pi, pi),
      p2 = if_else(y2 == y2_star, 1 - pi, pi),
      p3 = if_else(y3 == y3_star, 1 - pi, pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, age_inconsistent, educ_inconsistent) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'age_inconsistent', 'educ_inconsistent')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}

calc_lli_derivatives_3waves_ar1_misclassifiction_symptoms <- function(param_transformed) {
  
  param <- param_transformed
  param$theta_0 <- logit_inverse(param_transformed$theta_0)
  param$theta_1 <- logit_inverse(param_transformed$theta_1)
  
  df_probs_temp <- df_template_covariates %>% 
    mutate(
      pi = logit_inverse(param$intercept_pi + param$age_pi*age_inconsistent + param$educ_pi*educ_inconsistent),
      mu = param$theta_0/(param$theta_0 + param$theta_1)
    ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y1_star == 0 & y2_star == 1 ~ param$theta_0,
        y1_star == 1 & y2_star == 0 ~ param$theta_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
        y2_star == 0 & y3_star == 1 ~ param$theta_0,
        y2_star == 1 & y3_star == 0 ~ param$theta_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
      ),
      p1 = if_else(y1 == y1_star, 1 - pi, pi),
      p2 = if_else(y2 == y2_star, 1 - pi, pi),
      p3 = if_else(y3 == y3_star, 1 - pi, pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    mutate(
      d1_star_theta_0 = case_when(
        y1_star == 0 ~ -param$theta_1/((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ param$theta_1/((param$theta_1 + param$theta_0)^2)
      ),
      d1_star_theta_1 = case_when(
        y1_star == 0 ~ param$theta_0/((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ -param$theta_0/((param$theta_1 + param$theta_0)^2)
      ),
      d2_star_theta_0 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta_0 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta_1 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_d_theta_0 = 
        d1_star_theta_0*joint_p/p1_star + 
        d2_star_theta_0*joint_p/p2_star + 
        d3_star_theta_0*joint_p/p3_star,
      joint_d_theta_1 = 
        d1_star_theta_1*joint_p/p1_star + 
        d2_star_theta_1*joint_p/p2_star + 
        d3_star_theta_1*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, age_inconsistent, educ_inconsistent) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0), 
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      pi = mean(pi),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0), 
      joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1),
      joint_d_pi = joint_d_pi*pi*(1 - pi)
    ) %>% 
    mutate(
      joint_d_intercept_pi = joint_d_pi,
      joint_d_age_pi = joint_d_pi*age_inconsistent,
      joint_d_educ_pi = joint_d_pi*educ_inconsistent
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'age_inconsistent', 'educ_inconsistent')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_intercept_pi = weight*joint_d_intercept_pi/joint_p,
      lgi_age_pi = weight*joint_d_age_pi/joint_p,
      lgi_educ_pi = weight*joint_d_educ_pi/joint_p,
    ) %>%
    select(lgi_theta_0, lgi_theta_1, lgi_intercept_pi, lgi_age_pi, lgi_educ_pi)
}

calc_mle_3waves_ar1_misclassifiction_symptoms <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_misclassifiction_symptoms(param_transformed))
}

calc_mle_derivatives_3waves_ar1_misclassifiction_symptoms <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_misclassifiction_symptoms(param_transformed))
}


