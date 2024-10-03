
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

