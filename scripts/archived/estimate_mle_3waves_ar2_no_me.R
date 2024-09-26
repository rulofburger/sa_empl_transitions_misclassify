logit_transform <- function(param0) {
  param_input <- log(param0/(1 - param0))
  return(param_input)
}

logit_inverse <- function(param_input0) {
  param <- 1/(1 + exp(-param_input0))
  return(param)
}


calc_lli <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  param$pi <- 0
  Theta = param$theta_0*(param$theta_10 - 1) + param$theta_1*(param$theta_01 - 1)
  
  df_probs_temp <- df_template %>% 
    mutate(
      p12_star = case_when(
        y2_star == 0 & y1_star == 0 ~ param$theta_1*(param$theta_0 + param$theta_01 - 1)/Theta,
        y2_star == 0 & y1_star == 1 ~ -param$theta_1*param$theta_0/Theta,
        y2_star == 1 & y1_star == 0 ~ -param$theta_1*param$theta_0/Theta,
        y2_star == 1 & y1_star == 1 ~ param$theta_0*(param$theta_1 + param$theta_10 - 1)/Theta
      ),
      p3_star = case_when(
        y3_star == 1 & y2_star == 0 & y1_star == 0 ~ param$theta_0,
        y3_star == 1 & y2_star == 0 & y1_star == 1 ~ param$theta_0 + param$theta_01,
        y3_star == 1 & y2_star == 1 & y1_star == 0 ~ 1 - param$theta_1 - param$theta_10,
        y3_star == 1 & y2_star == 1 & y1_star == 1 ~ 1 - param$theta_1,
        y3_star == 0 & y2_star == 0 & y1_star == 0 ~ 1 - param$theta_0,
        y3_star == 0 & y2_star == 0 & y1_star == 1 ~ 1 - param$theta_0 - param$theta_01,
        y3_star == 0 & y2_star == 1 & y1_star == 0 ~ param$theta_1 + param$theta_10,
        y3_star == 0 & y2_star == 1 & y1_star == 1 ~ param$theta_1,
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p12_star*p1*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>% 
    pull(lli)
  
  # print(df_lli[1])
  # 
  # return(df_lli)
  
}

calc_ll <- function(param_transformed) {
  ll <- sum(calc_lli(param_transformed))
  return(ll)
}

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 

param_init <- data.frame(theta_0 = 0.05, theta_01 = 0.05, theta_10 = 0.05, theta_1 = 0.05)

param_init_transformed <- logit_transform(param_init)

calc_ll(param_init_transformed)




estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = param_init_transformed,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum