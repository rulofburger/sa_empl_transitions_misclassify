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
      p4_star = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ param$theta_0,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ param$theta_0 + param$theta_01,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ 1 - param$theta_1 - param$theta_10,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ 1 - param$theta_1,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ 1 - param$theta_0 - param$theta_01,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ param$theta_1 + param$theta_10,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ param$theta_1,
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      p4 = if_else(y4 == y4_star, 1 - param$pi, param$pi),
      joint_p = p12_star*p1*p2*p3_star*p3*p4_star*p4
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, y4) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'y4')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    # mutate(lli = log(joint_p)) %>% 
    pull(lli)
  
  # print(df_lli[1])
  # 
  # return(df_lli)
  
}


calc_lli_derivatives <- function(param_transformed, df_x) {
  
  param <- logit_inverse(param_transformed)
  Theta = param$theta_0*(-1 + param$theta_10) + param$theta_1*(-1 + param$theta_01)
  
  df_probs_temp <- df_template %>% 
    mutate(
      p12_star = case_when(
        y2_star == 0 & y1_star == 0 ~ param$theta_1*(param$theta_0 + param$theta_01 - 1)/Theta,
        y2_star == 0 & y1_star == 1 ~ -param$theta_1*param$theta_0/Theta,
        y2_star == 1 & y1_star == 0 ~ -param$theta_1*param$theta_0/Theta,
        y2_star == 1 & y1_star == 1 ~ param$theta_0*(param$theta_1 + param$theta_10 - 1)/Theta
      ),
      d12_star_theta_0 = case_when(
        y2_star == 0 & y1_star == 0 ~ param$theta_1*(1 + param$theta_1 - param$theta_10)*(-1 + param$theta_01)/(Theta^2),
        y2_star == 0 & y1_star == 1 ~ -(param$theta_1^2)*(-1 + param$theta_01)/(Theta^2),
        y2_star == 1 & y1_star == 0 ~ -(param$theta_1^2)*(-1 + param$theta_01)/(Theta^2),
        y2_star == 1 & y1_star == 1 ~ param$theta_1*(-1 + param$theta_1 + param$theta_10)*(-1 + param$theta_01)/(Theta^2)
      ),
      d12_star_theta_1 = case_when(
        y2_star == 0 & y1_star == 0 ~ param$theta_0*(-1 + param$theta_0 + param$theta_01)*(-1 + param$theta_10)/(Theta^2),
        y2_star == 0 & y1_star == 1 ~ -(param$theta_0^2)*(-1 + param$theta_10)/(Theta^2),
        y2_star == 1 & y1_star == 0 ~ -(param$theta_0^2)*(-1 + param$theta_10)/(Theta^2),
        y2_star == 1 & y1_star == 1 ~ param$theta_0*(1 + param$theta_0 - param$theta_01)*(-1 + param$theta_10)/(Theta^2)
      ),
      d12_star_theta_01 = case_when(
        y2_star == 0 & y1_star == 0 ~ param$theta_0*param$theta_1*(-1 - param$theta_1 + param$theta_10)/(Theta^2),
        y2_star == 0 & y1_star == 1 ~ param$theta_0*(param$theta_1^2)/(Theta^2),
        y2_star == 1 & y1_star == 0 ~ param$theta_0*(param$theta_1^2)/(Theta^2),
        y2_star == 1 & y1_star == 1 ~ -param$theta_0*param$theta_1*(-1 + param$theta_1 + param$theta_10)/(Theta^2)
      ),
      d12_star_theta_10 = case_when(
        y2_star == 0 & y1_star == 0 ~ -param$theta_0*param$theta_1*(-1 + param$theta_0 + param$theta_01)/(Theta^2),
        y2_star == 0 & y1_star == 1 ~ param$theta_1*(param$theta_0^2)/(Theta^2),
        y2_star == 1 & y1_star == 0 ~ param$theta_1*(param$theta_0^2)/(Theta^2),
        y2_star == 1 & y1_star == 1 ~ param$theta_0*param$theta_1*(-1 - param$theta_0 + param$theta_01)/(Theta^2)
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
      d3_star_theta_0 = case_when(
        y3_star == 1 & y2_star == 0 & y1_star == 0 ~ 1,
        y3_star == 1 & y2_star == 0 & y1_star == 1 ~ 1,
        y3_star == 1 & y2_star == 1 & y1_star == 0 ~ 0,
        y3_star == 1 & y2_star == 1 & y1_star == 1 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 0 ~ -1,
        y3_star == 0 & y2_star == 0 & y1_star == 1 ~ -1,
        y3_star == 0 & y2_star == 1 & y1_star == 0 ~ 0,
        y3_star == 0 & y2_star == 1 & y1_star == 1 ~ 0,
      ),
      d3_star_theta_1 = case_when(
        y3_star == 1 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 1 & y2_star == 0 & y1_star == 1 ~ 0,
        y3_star == 1 & y2_star == 1 & y1_star == 0 ~ -1,
        y3_star == 1 & y2_star == 1 & y1_star == 1 ~ -1,
        y3_star == 0 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 1 ~ 0,
        y3_star == 0 & y2_star == 1 & y1_star == 0 ~ 1,
        y3_star == 0 & y2_star == 1 & y1_star == 1 ~ 1,
      ),
      d3_star_theta_01 = case_when(
        y3_star == 1 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 1 & y2_star == 0 & y1_star == 1 ~ 1,
        y3_star == 1 & y2_star == 1 & y1_star == 0 ~ 0,
        y3_star == 1 & y2_star == 1 & y1_star == 1 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 1 ~ -1,
        y3_star == 0 & y2_star == 1 & y1_star == 0 ~ 0,
        y3_star == 0 & y2_star == 1 & y1_star == 1 ~ 0,
      ),
      d3_star_theta_10 = case_when(
        y3_star == 1 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 1 & y2_star == 0 & y1_star == 1 ~ 0,
        y3_star == 1 & y2_star == 1 & y1_star == 0 ~ -1,
        y3_star == 1 & y2_star == 1 & y1_star == 1 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 0 ~ 0,
        y3_star == 0 & y2_star == 0 & y1_star == 1 ~ 0,
        y3_star == 0 & y2_star == 1 & y1_star == 0 ~ 1,
        y3_star == 0 & y2_star == 1 & y1_star == 1 ~ 0,
      ),
      p4_star = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ param$theta_0,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ param$theta_0 + param$theta_01,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ 1 - param$theta_1 - param$theta_10,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ 1 - param$theta_1,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ 1 - param$theta_0 - param$theta_01,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ param$theta_1 + param$theta_10,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ param$theta_1,
      ),
      d4_star_theta_0 = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ 1,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ 1,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ 0,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ -1,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ -1,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ 0,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ 0,
      ),
      d4_star_theta_1 = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ 0,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ -1,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ -1,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ 0,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ 1,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ 1,
      ),
      d4_star_theta_01 = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ 1,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ 0,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ -1,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ 0,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ 0,
      ),
      d4_star_theta_10 = case_when(
        y4_star == 1 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 1 & y3_star == 0 & y2_star == 1 ~ 0,
        y4_star == 1 & y3_star == 1 & y2_star == 0 ~ -1,
        y4_star == 1 & y3_star == 1 & y2_star == 1 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 0 ~ 0,
        y4_star == 0 & y3_star == 0 & y2_star == 1 ~ 0,
        y4_star == 0 & y3_star == 1 & y2_star == 0 ~ 1,
        y4_star == 0 & y3_star == 1 & y2_star == 1 ~ 0,
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      p4 = if_else(y4 == y4_star, 1 - param$pi, param$pi),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      d4_pi = if_else(y4 == y4_star, -1, 1),
      joint_p = p12_star*p1*p2*p3_star*p3*p4_star*p4,
      joint_d_theta_0 = 
        d12_star_theta_0*p3_star*p4_star*p1*p2*p3*p4 +
        p12_star*d3_star_theta_0*p4_star*p1*p2*p3*p4 + 
        p12_star*p3_star*d4_star_theta_0*p1*p2*p3*p4,
      joint_d_theta_1 = 
        d12_star_theta_1*p3_star*p4_star*p1*p2*p3*p4 +
        p12_star*d3_star_theta_1*p4_star*p1*p2*p3*p4 + 
        p12_star*p3_star*d4_star_theta_1*p1*p2*p3*p4,
      joint_d_theta_01 = 
        d12_star_theta_01*p3_star*p4_star*p1*p2*p3*p4 +
        p12_star*d3_star_theta_01*p4_star*p1*p2*p3*p4 + 
        p12_star*p3_star*d4_star_theta_01*p1*p2*p3*p4,
      joint_d_theta_10 = 
        d12_star_theta_10*p3_star*p4_star*p1*p2*p3*p4 +
        p12_star*d3_star_theta_10*p4_star*p1*p2*p3*p4 + 
        p12_star*p3_star*d4_star_theta_10*p1*p2*p3*p4,
      joint_d_pi = 
        p12_star*p3_star*p4_star*d1_pi*p2*p3*p4 + 
        p12_star*p3_star*p4_star*p1*d2_pi*p3*p4 + 
        p12_star*p3_star*p4_star*p1*p2*d3_pi*p4 + 
        p12_star*p3_star*p4_star*p1*p2*p3*d4_pi
    ) 
  
  # df_probs <- df_probs_temp %>% 
  #   group_by(y1, y2, y3, y4) %>% 
  #   summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, y4) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0),
      joint_d_theta_1 = sum(joint_d_theta_1), 
      joint_d_theta_01 = sum(joint_d_theta_01), 
      joint_d_theta_10 = sum(joint_d_theta_10), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      .groups = "drop")
  
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'y4')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_01 = weight*joint_d_theta_01/joint_p,
      lgi_theta_10 = weight*joint_d_theta_10/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
      # lgi_theta_0 = joint_d_theta_0/joint_p,
      # lgi_theta_01 =joint_d_theta_01/joint_p,
      # lgi_theta_10 = joint_d_theta_10/joint_p,
      # lgi_theta_1 = joint_d_theta_1/joint_p,
      # lgi_pi = joint_d_pi/joint_p
      ) %>% 
    select(lgi_theta_0, lgi_theta_01, lgi_theta_10, lgi_theta_1, lgi_pi)

}
  



calc_ll <- function(param_transformed) {
  ll <- sum(calc_lli(param_transformed))
  return(ll)
}

calc_lg <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives(param_transformed))
  return(lg)
}