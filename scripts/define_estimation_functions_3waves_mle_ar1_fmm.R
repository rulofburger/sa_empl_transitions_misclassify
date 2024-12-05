# DEFINE FUNCTIONS ====

#> General functions that restrict parameter values to unit interval

logit_transform <- function(param0) {
  param_input <- log(param0/(1 - param0))
  return(param_input)
}

logit_inverse <- function(param_input0) {
  param <- 1/(1 + exp(-param_input0))
  return(param)
}

# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH 2 FMM groups ====

calc_lli_3waves_ar1_fmm2 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
  mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
  if(pi0) param$pi <- 0
  
  df_probs_temp1 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_1,
        y1_star == 0 & y2_star == 1 ~ param$theta0_1,
        y1_star == 1 & y2_star == 0 ~ param$theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_1,
        y2_star == 0 & y3_star == 1 ~ param$theta0_1,
        y2_star == 1 & y3_star == 0 ~ param$theta1_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs_temp2 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_2, 1 - mu_2),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_2,
        y1_star == 0 & y2_star == 1 ~ param$theta0_2,
        y1_star == 1 & y2_star == 0 ~ param$theta1_2,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_2
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_2,
        y2_star == 0 & y3_star == 1 ~ param$theta0_2,
        y2_star == 1 & y3_star == 0 ~ param$theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs1 <- df_probs_temp1 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p1 = sum(joint_p), .groups = "drop")
  
  df_probs2 <- df_probs_temp2 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p2 = sum(joint_p), .groups = "drop")
  
  df_probs <- df_probs1 %>% 
    left_join(df_probs2, by = c('y1', 'y2', 'y3')) %>% 
    mutate(joint_p = param$p_1*joint_p1 + (1 - param$p_1)*joint_p2)
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}


calc_lli_derivatives_3waves_ar1_fmm2 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
  mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
  if(pi0) param$pi <- 0
  
  df_probs_temp1 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_1,
        y1_star == 0 & y2_star == 1 ~ param$theta0_1,
        y1_star == 1 & y2_star == 0 ~ param$theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_1,
        y2_star == 0 & y3_star == 1 ~ param$theta0_1,
        y2_star == 1 & y3_star == 0 ~ param$theta1_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_1 = case_when(
        y1_star == 0 ~ -param$theta1_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1 ~ param$theta1_1/((param$theta1_1 + param$theta0_1)^2)
      ),
      d1_star_theta1_1 = case_when(
        y1_star == 0 ~ param$theta0_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1 ~ -param$theta0_1/((param$theta1_1 + param$theta0_1)^2)
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
      d3_star_theta0_1 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_1 = case_when(
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
        d2_star_theta0_1*joint_p/p2_star + 
        d3_star_theta0_1*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star + 
        d3_star_theta1_1*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_probs_temp2 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_2, 1 - mu_2),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_2,
        y1_star == 0 & y2_star == 1 ~ param$theta0_2,
        y1_star == 1 & y2_star == 0 ~ param$theta1_2,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_2
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_2,
        y2_star == 0 & y3_star == 1 ~ param$theta0_2,
        y2_star == 1 & y3_star == 0 ~ param$theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_2 = case_when(
        y1_star == 0 ~ -param$theta1_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1 ~ param$theta1_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d1_star_theta1_2 = case_when(
        y1_star == 0 ~ param$theta0_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1 ~ -param$theta0_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d2_star_theta0_2 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_2 = case_when(
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
      joint_d_theta0_2 = 
        d1_star_theta0_2*joint_p/p1_star + 
        d2_star_theta0_2*joint_p/p2_star + 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_2 = 
        d1_star_theta1_2*joint_p/p1_star + 
        d2_star_theta1_2*joint_p/p2_star + 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    )  
  
  df_grad1 <- df_probs_temp1 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1), 
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_pi1 = sum(joint_d_pi),
      joint_p1 = sum(joint_p),
      .groups = "drop") 
  
  df_grad2 <- df_probs_temp2 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta0_2 = sum(joint_d_theta0_2), 
      joint_d_theta1_2 = sum(joint_d_theta1_2),
      joint_d_pi2 = sum(joint_d_pi),
      joint_p2 = sum(joint_p),
      .groups = "drop") 
  
  df_grad <- df_grad1 %>% 
    left_join(df_grad2, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      joint_p = param$p_1*joint_p1 + (1 - param$p_1)*joint_p2,
      joint_d_pi = param$p_1*joint_d_pi1 + (1 - param$p_1)*joint_d_pi2
      ) %>% 
    mutate(
      joint_d_theta0_1 = param$p_1*joint_d_theta0_2*param$theta0_1*(1 - param$theta0_1), 
      joint_d_theta1_1 = param$p_1*joint_d_theta1_2*param$theta1_1*(1 - param$theta1_1),
      joint_d_theta0_2 = (1 - param$p_1)*joint_d_theta0_2*param$theta0_2*(1 - param$theta0_2), 
      joint_d_theta1_2 = (1 - param$p_1)*joint_d_theta1_2*param$theta1_2*(1 - param$theta1_2),
      joint_d_p_1 = joint_p1*param$p_1*(1 - param$p_1) - joint_p2*param$p_1*(1 - param$p_1),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta0_1 = weight*joint_d_theta0_1/joint_p,
      lgi_theta1_1 = weight*joint_d_theta1_1/joint_p,
      lgi_theta0_2 = weight*joint_d_theta0_2/joint_p,
      lgi_theta1_2 = weight*joint_d_theta1_2/joint_p,
      lgi_p_1 = weight*joint_d_p_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
    ) %>%
    select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_pi)
  
  if(pi0) df_gi <- df_gi %>% select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2)
  return(df_gi)
  
}

calc_mle_3waves_ar1_fmm2 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_fmm2(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_fmm2 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_fmm2(param_transformed))
  return(lg)
}

# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH 3 FMM groups ====

calc_lli_3waves_ar1_fmm3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
  mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
  mu_3 <- param$theta0_3/(param$theta0_3 + param$theta1_3)
  if(pi0) param$pi <- 0
  
  df_probs_temp1 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_1,
        y1_star == 0 & y2_star == 1 ~ param$theta0_1,
        y1_star == 1 & y2_star == 0 ~ param$theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_1,
        y2_star == 0 & y3_star == 1 ~ param$theta0_1,
        y2_star == 1 & y3_star == 0 ~ param$theta1_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs_temp2 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_2, 1 - mu_2),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_2,
        y1_star == 0 & y2_star == 1 ~ param$theta0_2,
        y1_star == 1 & y2_star == 0 ~ param$theta1_2,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_2
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_2,
        y2_star == 0 & y3_star == 1 ~ param$theta0_2,
        y2_star == 1 & y3_star == 0 ~ param$theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs_temp3 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_3, 1 - mu_3),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_3,
        y1_star == 0 & y2_star == 1 ~ param$theta0_3,
        y1_star == 1 & y2_star == 0 ~ param$theta1_3,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_3
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_3,
        y2_star == 0 & y3_star == 1 ~ param$theta0_3,
        y2_star == 1 & y3_star == 0 ~ param$theta1_3,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_3
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  
  df_probs1 <- df_probs_temp1 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p1 = sum(joint_p), .groups = "drop")
  
  df_probs2 <- df_probs_temp2 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p2 = sum(joint_p), .groups = "drop")
  
  df_probs3 <- df_probs_temp3 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p3 = sum(joint_p), .groups = "drop")
  
  
  df_probs <- df_probs1 %>% 
    left_join(df_probs2, by = c('y1', 'y2', 'y3')) %>% 
    left_join(df_probs3, by = c('y1', 'y2', 'y3')) %>% 
    mutate(joint_p = param$p_1*joint_p1 + param$p_2*joint_p2 + (1 - param$p_1 - param$p_2)*joint_p3)
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}


calc_lli_derivatives_3waves_ar1_fmm3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
  mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
  mu_3 <- param$theta0_3/(param$theta0_3 + param$theta1_3)
  if(pi0) param$pi <- 0
  
  df_probs_temp1 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_1,
        y1_star == 0 & y2_star == 1 ~ param$theta0_1,
        y1_star == 1 & y2_star == 0 ~ param$theta1_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_1,
        y2_star == 0 & y3_star == 1 ~ param$theta0_1,
        y2_star == 1 & y3_star == 0 ~ param$theta1_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_1 = case_when(
        y1_star == 0 ~ -param$theta1_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1 ~ param$theta1_1/((param$theta1_1 + param$theta0_1)^2)
      ),
      d1_star_theta1_1 = case_when(
        y1_star == 0 ~ param$theta0_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1 ~ -param$theta0_1/((param$theta1_1 + param$theta0_1)^2)
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
      d3_star_theta0_1 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_1 = case_when(
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
        d2_star_theta0_1*joint_p/p2_star + 
        d3_star_theta0_1*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star + 
        d3_star_theta1_1*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_probs_temp2 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_2, 1 - mu_2),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_2,
        y1_star == 0 & y2_star == 1 ~ param$theta0_2,
        y1_star == 1 & y2_star == 0 ~ param$theta1_2,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_2
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_2,
        y2_star == 0 & y3_star == 1 ~ param$theta0_2,
        y2_star == 1 & y3_star == 0 ~ param$theta1_2,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_2
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_2 = case_when(
        y1_star == 0 ~ -param$theta1_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1 ~ param$theta1_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d1_star_theta1_2 = case_when(
        y1_star == 0 ~ param$theta0_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1 ~ -param$theta0_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d2_star_theta0_2 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_2 = case_when(
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
      joint_d_theta0_2 = 
        d1_star_theta0_2*joint_p/p1_star + 
        d2_star_theta0_2*joint_p/p2_star + 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_2 = 
        d1_star_theta1_2*joint_p/p1_star + 
        d2_star_theta1_2*joint_p/p2_star + 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    )  
  
  df_probs_temp3 <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu_3, 1 - mu_3),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta0_3,
        y1_star == 0 & y2_star == 1 ~ param$theta0_3,
        y1_star == 1 & y2_star == 0 ~ param$theta1_3,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta1_3
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta0_3,
        y2_star == 0 & y3_star == 1 ~ param$theta0_3,
        y2_star == 1 & y3_star == 0 ~ param$theta1_3,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta1_3
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_3 = case_when(
        y1_star == 0 ~ -param$theta1_3/((param$theta1_3 + param$theta0_3)^2),
        y1_star == 1 ~ param$theta1_3/((param$theta1_3 + param$theta0_3)^2)
      ),
      d1_star_theta1_3 = case_when(
        y1_star == 0 ~ param$theta0_3/((param$theta1_3 + param$theta0_3)^2),
        y1_star == 1 ~ -param$theta0_3/((param$theta1_3 + param$theta0_3)^2)
      ),
      d2_star_theta0_3 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta1_3 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta0_3 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta1_3 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_d_theta0_3 = 
        d1_star_theta0_3*joint_p/p1_star + 
        d2_star_theta0_3*joint_p/p2_star + 
        d3_star_theta0_3*joint_p/p3_star,
      joint_d_theta1_3 = 
        d1_star_theta1_3*joint_p/p1_star + 
        d2_star_theta1_3*joint_p/p2_star + 
        d3_star_theta1_3*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    )  
  
  
  df_grad1 <- df_probs_temp1 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1), 
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_pi1 = sum(joint_d_pi),
      joint_p1 = sum(joint_p),
      .groups = "drop") 
  
  df_grad2 <- df_probs_temp2 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta0_2 = sum(joint_d_theta0_2), 
      joint_d_theta1_2 = sum(joint_d_theta1_2),
      joint_d_pi2 = sum(joint_d_pi),
      joint_p2 = sum(joint_p),
      .groups = "drop") 
  
  df_grad3 <- df_probs_temp3 %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta0_3 = sum(joint_d_theta0_3), 
      joint_d_theta1_3 = sum(joint_d_theta1_3),
      joint_d_pi3 = sum(joint_d_pi),
      joint_p3 = sum(joint_p),
      .groups = "drop") 
  
  df_grad <- df_grad1 %>% 
    left_join(df_grad2, by = c('y1', 'y2', 'y3')) %>% 
    left_join(df_grad3, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      joint_p = param$p_1*joint_p1 + param$p_2*joint_p2 + (1 - param$p_1 - param$p_2)*joint_p3,
      joint_d_pi = param$p_1*joint_d_pi1 + param$p_2*joint_d_pi2 + (1 - param$p_1 - param$p_2)*joint_d_pi3
    ) %>% 
    mutate(
      joint_d_theta0_1 = param$p_1*joint_d_theta0_2*param$theta0_1*(1 - param$theta0_1), 
      joint_d_theta1_1 = param$p_1*joint_d_theta1_2*param$theta1_1*(1 - param$theta1_1),
      joint_d_theta0_2 = param$p_2*joint_d_theta0_2*param$theta0_2*(1 - param$theta0_2), 
      joint_d_theta1_2 = param$p_2*joint_d_theta1_2*param$theta1_2*(1 - param$theta1_2),
      joint_d_theta0_3 = (1 - param$p_1 - param$p_2)*joint_d_theta0_3*param$theta0_3*(1 - param$theta0_3), 
      joint_d_theta1_3 = (1 - param$p_1 - param$p_2)*joint_d_theta1_3*param$theta1_3*(1 - param$theta1_3),
      
      joint_d_p_1 = joint_p1*param$p_1*(1 - param$p_1) - joint_p3*param$p_1*(1 - param$p_1),
      joint_d_p_2 = joint_p2*param$p_2*(1 - param$p_2) - joint_p3*param$p_2*(1 - param$p_2),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta0_1 = weight*joint_d_theta0_1/joint_p,
      lgi_theta1_1 = weight*joint_d_theta1_1/joint_p,
      lgi_theta0_2 = weight*joint_d_theta0_2/joint_p,
      lgi_theta1_2 = weight*joint_d_theta1_2/joint_p,
      lgi_theta0_3 = weight*joint_d_theta0_3/joint_p,
      lgi_theta1_3 = weight*joint_d_theta1_3/joint_p,
      lgi_p_1 = weight*joint_d_p_1/joint_p,
      lgi_p_2 = weight*joint_d_p_2/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
    ) %>%
    select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_p_2, lgi_theta0_3, lgi_theta1_3, lgi_pi)
  
  if(pi0) df_gi <- df_gi %>% select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_p_2, lgi_theta0_3, lgi_theta1_3)
  return(df_gi)
  
}

calc_mle_3waves_ar1_fmm3 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_fmm3(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_fmm3 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_fmm3(param_transformed))
  return(lg)
}
