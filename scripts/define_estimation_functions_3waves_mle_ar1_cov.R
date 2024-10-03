
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


