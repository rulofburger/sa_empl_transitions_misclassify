
# takes in probability (e.g. job entry rate) and outputs logit transformed
logit_transform <- function(param0) {
  param_input <- log(param0/(1 - param0))
  return(param_input)
}

# produces probability from logit parameters
logit_inverse <- function(param_input0) {
  param <- 1/(1 + exp(-param_input0))
  return(param)
}

logistic_derivative <- function(param) {
  logit_transform(param)*(1 - logit_transform(param))
}

calc_df_probs_cov <- function(df_cov, param_transformed) {
  
  df_probs <- df_cov |>
    mutate(
      p1_star = if_else(y1_star == 1, 
                        logit_inverse(param_transformed$mu), 
                        1 - logit_inverse(param_transformed$mu)),
      
      p2_star = if_else(y1_star == 1, 
                        if_else(y2_star == 1, 
                                logit_inverse(param_transformed$beta01 + 
                                                param_transformed$beta1*age1 + 
                                                param_transformed$beta2*age1*age1 + 
                                                param_transformed$beta3A*educ1), 
                                1 - logit_inverse(param_transformed$beta01 + 
                                                    param_transformed$beta1*age1 + 
                                                    param_transformed$beta2*age1*age1 + 
                                                    param_transformed$beta3A*educ1)), 
                        if_else(y2_star == 1, 
                                logit_inverse(param_transformed$beta02 + 
                                                param_transformed$beta1*age1 + 
                                                param_transformed$beta2*age1*age1 + 
                                                param_transformed$beta3B*educ1), 
                                1 - logit_inverse(param_transformed$beta02 + 
                                                    param_transformed$beta1*age1 + 
                                                    param_transformed$beta2*age1*age1 + 
                                                    param_transformed$beta3B*educ1))),
      
      p3_star = if_else(y2_star == 1, 
                        if_else(y3_star == 1, 
                                logit_inverse(param_transformed$beta01 + 
                                                param_transformed$beta1*age2 + 
                                                param_transformed$beta2*age2*age2 + 
                                                param_transformed$beta3A*educ2),
                                1 - logit_inverse(param_transformed$beta01 + 
                                                    param_transformed$beta1*age2 + 
                                                    param_transformed$beta2*age2*age2 + 
                                                    param_transformed$beta3A*educ2)), 
                        if_else(y3_star == 1, 
                                logit_inverse(param_transformed$beta02 + 
                                                param_transformed$beta1*age2 + 
                                                param_transformed$beta2*age2*age2 + 
                                                param_transformed$beta3B*educ2),
                                1 - logit_inverse(param_transformed$beta02 + 
                                                    param_transformed$beta1*age2 + 
                                                    param_transformed$beta2*age2*age2 + 
                                                    param_transformed$beta3B*educ2))),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3)
  
  df_probs
}



# input logit transformed parameters (i.e. not probabilities)
calc_lli_3waves_ar1_cov <- function(df_template, data, param_transformed) {
  
  # Note: the theta parameters are expanded into linear functions
  #         of the beta parameters to include the covariates
  #         so that:
  #                  1 - theta1 = beta01 + beta1*Age + beta2*Age^2 + beta3A*Educ
  #                      theta2 = beta02 + beta1*Age + beta2*Age^2 + beta3B*Educ
  #         where:
  #                  logit_inverse(theta1) is job exit rate
  #                  logit_inverse(theta2) is job entry rate
  
  
  # Make probabilities
  #param <- logit_inverse(param_transformed)
  #mu    <- param$theta_2/(param$theta_1 + param$theta_2)
  
  df_cov <- joyn::joyn(x          = df_template, 
                       y          = data, 
                       by         = c("y1", "y2", "y3"), 
                       keep       = "left", 
                       match_type = "1:m", 
                       reportvar  = FALSE)
  
  df_probs <- calc_df_probs_cov(df_cov            = df_cov, 
                                param_transformed = param_transformed)
  
  
  df_probs <- df_probs_temp |>
    group_by(y1, y2, y3) |>
    summarise(joint_p = sum(joint_p), 
              .groups = "drop")
  
  df_lli <- df_estimate |>
    joyn::left_join(df_probs, 
                    by = c('y1', 'y2', 'y3')) |>
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
  df_lli
  
}


calc_lli_derivatives_3waves_ar1_cov <- function(param_transformed, df_x) {
  
  
  df_probs <- calc_df_probs_cov(df_cov            = df_cov, 
                                param_transformed = param_transformed)
  
  df_probs <- df_probs |> 
    mutate(
      d1_star_theta_1 = case_when(
        y1_star == 1 ~ logit_inverse(param_transformed$beta02 + 
                                       param_transformed$beta1*age1 + 
                                       param_transformed$beta2*age1*age1 + 
                                       param_transformed$beta3B*educ1) /    # job entry
          ((1 - logit_inverse(param_transformed$beta01 + 
                                param_transformed$beta1*age1 + 
                                param_transformed$beta2*age1*age1 + 
                                param_transformed$beta3A*educ1) +      # job exit
              logit_inverse(param_transformed$beta02 + 
                              param_transformed$beta1*age1 + 
                              param_transformed$beta2*age1*age1 + 
                              param_transformed$beta3B*educ1))^2), # job entry
        
        y1_star == 0 ~ - logit_inverse(param_transformed$beta02 + 
                                         param_transformed$beta1*age1 + 
                                         param_transformed$beta2*age1*age1 + 
                                         param_transformed$beta3B*educ1) /
          ((1 - logit_inverse(param_transformed$beta01 + 
                                param_transformed$beta1*age1 + 
                                param_transformed$beta2*age1*age1 + 
                                param_transformed$beta3A*educ1)
            + logit_inverse(param_transformed$beta02 + 
                              param_transformed$beta1*age1 + 
                              param_transformed$beta2*age1*age1 + 
                              param_transformed$beta3B*educ1))^2)
      ),
      d1_star_theta_2 = case_when(
        y1_star == 1 ~ - logit_inverse(param_transformed$beta02 + 
                                         param_transformed$beta1*age1 + 
                                         param_transformed$beta2*age1*age1 + 
                                         param_transformed$beta3B*educ1) / 
          ((logit_inverse(param_transformed$beta01 + 
                            param_transformed$beta1*age1 + 
                            param_transformed$beta2*age1*age1 + 
                            param_transformed$beta3A*educ1) +
              logit_inverse(param_transformed$beta02 + 
                              param_transformed$beta1*age1 + 
                              param_transformed$beta2*age1*age1 + 
                              param_transformed$beta3B*educ1))^2) + 
          1/(logit_inverse(param_transformed$beta01 + 
                             param_transformed$beta1*age1 + 
                             param_transformed$beta2*age1*age1 + 
                             param_transformed$beta3A*educ1) + 
               logit_inverse(param_transformed$beta02 + 
                               param_transformed$beta1*age1 + 
                               param_transformed$beta2*age1*age1 + 
                               param_transformed$beta3B*educ1)),
        y1_star == 0 ~ param$theta_2/
          ((logit_inverse(param_transformed$beta01 + 
                            param_transformed$beta1*age1 + 
                            param_transformed$beta2*age1*age1 + 
                            param_transformed$beta3A*educ1) + 
              logit_inverse(param_transformed$beta02 + 
                              param_transformed$beta1*age1 + 
                              param_transformed$beta2*age1*age1 + 
                              param_transformed$beta3B*educ1))^2) - 
          1/(logit_inverse(param_transformed$beta01 + 
                             param_transformed$beta1*age1 + 
                             param_transformed$beta2*age1*age1 + 
                             param_transformed$beta3A*educ1) + 
               logit_inverse(param_transformed$beta02 + 
                               param_transformed$beta1*age1 + 
                               param_transformed$beta2*age1*age1 + 
                               param_transformed$beta3B*educ1))
      ),
      d2_star_theta_1 = case_when(
        y2_star == 0 & y1_star == 0 ~ 0,
        y2_star == 0 & y1_star == 1 ~ -1,
        y2_star == 1 & y1_star == 0 ~ 0,
        y2_star == 1 & y1_star == 1 ~ 1
      ),
      d2_star_theta_2 = case_when(
        y2_star == 0 & y1_star == 0 ~ -1,
        y2_star == 0 & y1_star == 1 ~ 0,
        y2_star == 1 & y1_star == 0 ~ 1,
        y2_star == 1 & y1_star == 1 ~ 0
      ),
      d3_star_theta_1 = case_when(
        y3_star == 0 & y2_star == 0 ~ 0,
        y3_star == 0 & y2_star == 1 ~ -1,
        y3_star == 1 & y2_star == 0 ~ 0,
        y3_star == 1 & y2_star == 1 ~ 1
      ),
      d3_star_theta_2 = case_when(
        y3_star == 0 & y2_star == 0 ~ -1,
        y3_star == 0 & y2_star == 1 ~ 0,
        y3_star == 1 & y2_star == 0 ~ 1,
        y3_star == 1 & y2_star == 1 ~ 0
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      
      
      joint_d_theta_1 = 
        d1_star_theta_1*p1*p2_star*p2*p3_star*p3 + 
        p1_star*p1*d2_star_theta_1*p2*p3_star*p3 + 
        p1_star*p1*p2_star*p2*d3_star_theta_1*p3,
      joint_d_theta_2 = 
        d1_star_theta_2*p1*p2_star*p2*p3_star*p3 + 
        p1_star*p1*d2_star_theta_2*p2*p3_star*p3 + 
        p1_star*p1*p2_star*p2*d3_star_theta_2*p3,
      joint_d_pi = 
        p1_star*d1_pi*p2_star*p2*p3_star*p3 + 
        p1_star*p1*p2_star*d2_pi*p3_star*p3 + 
        p1_star*p1*p2_star*p2*p3_star*d3_pi,
      
      joint_d_beta_01 = 
        joint_d_theta_1*1 + 
        joint_d_theta_2*0, 
      
      joint_d_beta_02 = 
        joint_d_theta_1*0 +
        joint_d_theta_2*1 ,
      
      joint_d_beta_1  = 
        joint_d_theta_1*(age1 + age2 + age3) +
        joint_d_theta_2*(age1 + age2 + age3) ,
      
      joint_d_beta_2  = 
        joint_d_theta_1*(age1*age1 + age2*age2 + age3*age3) +
        joint_d_theta_2*(age1*age1 + age2*age2 + age3*age3) ,
      
      joint_d_beta_3A = 
        joint_d_theta_1*(educ1 + educ2 + educ3) +
        joint_d_theta_2*0 ,
      
      joint_d_beta_3B = 
        joint_d_theta_1*0 +
        joint_d_theta_2*(educ1 + educ2 + educ3)) 
  
  df_grad <- df_probs_temp |>
    group_by(y1, y2, y3) |>
    summarise(
      joint_d_beta_01 = sum(joint_d_beta_01),
      joint_d_beta_02 = sum(joint_d_beta_02), 
      joint_d_beta_1  = sum(joint_d_beta_1), 
      joint_d_beta_2  = sum(joint_d_beta_2), 
      joint_d_beta_3A = sum(joint_d_beta_3A), 
      joint_d_beta_3B = sum(joint_d_beta_3B), 
      joint_d_pi      = sum(joint_d_pi),
      joint_p         = sum(joint_p),
      .groups         = "drop")
  
  df_gi <- df_estimate |>
    joyn::left_join(df_grad, 
                    by        = c('y1', 'y2', 'y3'), 
                    reportvar = FALSE) |>
    mutate(
      lgi_beta_01 = weight*joint_d_beta_01/joint_p,
      lgi_beta_02 = weight*lgi_beta_02/joint_p,
      lgi_beta_1  = weight*lgi_beta_1/joint_p,
      lgi_beta_2  = weight*lgi_beta_2/joint_p,
      lgi_beta_3A = weight*lgi_beta_3A/joint_p,
      lgi_beta_3B = weight*lgi_beta_3B/joint_p,
      lgi_pi      = weight*lgi_pi/joint_p) |>
    select(contains("lgi"))
  
}


calc_mle_3waves_ar1_cov <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_cov(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_cov <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_cov(param_transformed))
  return(lg)
}


