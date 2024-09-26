logit_transform <- function(param0) {
  param_input <- log(param0/(1 - param0))
  return(param_input)
}

logit_inverse <- function(param_input0) {
  param <- 1/(1 + exp(-param_input0))
  return(param)
}

calc_lli_3waves <- function(param_transformed, gradient = T) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_2, 1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_2, 1 - param$theta_2)),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  print(df_probs_temp)
  
  if(!gradient) {
    df_probs <- df_probs_temp %>% 
      group_by(y1, y2, y3) %>% 
      summarise(joint_p = sum(joint_p), .groups = "drop")
    
    df_lli <- df_estimate %>% 
      left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
      mutate(lli = weight*log(joint_p)) %>% 
      pull(lli)
    
  }
  
  if(gradient) {
  browser()
    
      df_grads <- df_probs_temp %>% 
        mutate(
          pi_help = (y1 == y1_star) + (y2 == y2_star) + (y3 == y3_star),
          d_probs_pi = case_when(
            pi_help == 0 ~ joint_p*(3/param$pi - 0/(1 - param$pi))*(1 - param$pi)*param$pi,
            pi_help == 1 ~ joint_p*(2/param$pi - 1/(1 - param$pi))*(1 - param$pi)*param$pi,
            pi_help == 2 ~ joint_p*(1/param$pi - 2/(1 - param$pi))*(1 - param$pi)*param$pi,
            pi_help == 3 ~ joint_p*(0/param$pi - 3/(1 - param$pi))*(1 - param$pi)*param$pi
          )
        ) %>% 
        mutate(
          theta_1_help0 = (y1_star == 1)*(y2_star == 1) + (y2_star == 1)*(y3_star == 1),
          theta_1_help1 = (y1_star == 1)*(y2_star == 0) + (y2_star == 1)*(y3_star == 0),
          theta_2_help0 = (y1_star == 0)*(y2_star == 0) + (y2_star == 0)*(y3_star == 0),
          theta_2_help1 = (y1_star == 0)*(y2_star == 1) + (y2_star == 0)*(y3_star == 1),
          d_probs_theta_1 = case_when(
            theta_1_help0 == 2                                ~ joint_p*(-2/(1 - param$theta_1)                   - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_1)*param$theta_1,
            theta_1_help0 == 1 & theta_1_help1 == 1           ~ joint_p*(-1/(1 - param$theta_1) + 1/param$theta_1 - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_1)*param$theta_1,
            theta_1_help0 == 1 & theta_1_help1 == 0           ~ joint_p*(-1/(1 - param$theta_1) + param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))*(1 - param$theta_1)*param$theta_1,
            theta_1_help1 == 1 & theta_1_help0 == 0 & y1_star == 1 ~ joint_p*(1/(param$theta_1) - 1/(param$theta_1 + param$theta_2))                                 *(1 - param$theta_1)*param$theta_1,
            theta_1_help1 == 1 & theta_1_help0 == 0 & y1_star == 0 ~ joint_p*(1/(param$theta_1) + param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))     *(1 - param$theta_1)*param$theta_1,
            theta_1_help1 == 0 & theta_1_help0 == 0           ~ joint_p*(param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))                         *(1 - param$theta_1)*param$theta_1
          ),
          d_probs_theta_2 = case_when(
            theta_2_help0 == 2                                ~ joint_p*(-2/(1 - param$theta_2)                   - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_2)*param$theta_2,
            theta_2_help0 == 1 & theta_2_help1 == 1           ~ joint_p*(-1/(1 - param$theta_2) + 1/param$theta_2 - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_2)*param$theta_2,
            theta_2_help0 == 1 & theta_2_help1 == 0           ~ joint_p*(-1/(1 - param$theta_2) + param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))*(1 - param$theta_2)*param$theta_2,
            theta_2_help1 == 1 & theta_2_help0 == 0 & y1_star == 0 ~ joint_p*(1/(param$theta_2) - 1/(param$theta_1 + param$theta_2))                                 *(1 - param$theta_2)*param$theta_2,
            theta_2_help1 == 1 & theta_2_help0 == 0 & y1_star == 1 ~ joint_p*(1/(param$theta_2) + param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))     *(1 - param$theta_2)*param$theta_2,
            theta_2_help1 == 0 & theta_2_help0 == 0           ~ joint_p*(param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))                         *(1 - param$theta_2)*param$theta_2
          )
        ) %>% 
        group_by(y1, y2, y3) %>% 
        summarise(
          joint_p = sum(joint_p),
          d_probs_theta_1 = sum(d_probs_theta_1),
          d_probs_theta_2 = sum(d_probs_theta_2),
          d_probs_pi = sum(d_probs_pi),
          .groups = "drop"
        )
      
      print(df_grads)
      
      df_lli <- df_estimate %>% 
        left_join(df_grads, by = c('y1', 'y2', 'y3')) %>% 
        mutate(
          lli = log(joint_p),
          lgi_theta_1 = d_probs_theta_1/joint_p,
          lgi_theta_2 = d_probs_theta_2/joint_p,
          lgi_pi = d_probs_pi/joint_p
        ) %>% 
        select(lli, lgi_theta_1, lgi_theta_2, lgi_pi)
      
  }
  
  return(df_lli)
  
}

calc_mle_3waves <- function(param_transformed) {
  ll <- sum(calc_lli_3waves(param_transformed, gradient = F))
  return(ll)
}

calc_mle_derivatives_3waves <- function(param_transformed) {
  lg <- calc_lli_3waves(param_transformed, gradient = T) %>% select(-lli) %>% colSums
  return(lg)
}


calc_lli_3waves_pi0 <- function(param_transformed, gradient = T) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  pi <- 0
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_2, 1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_2, 1 - param$theta_2)),
      p1 = if_else(y1 == y1_star, 1 - pi, pi),
      p2 = if_else(y2 == y2_star, 1 - pi, pi),
      p3 = if_else(y3 == y3_star, 1 - pi, pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  # print(df_probs_temp)
  
  if(!gradient) {
    df_probs <- df_probs_temp %>% 
      group_by(y1, y2, y3) %>% 
      summarise(joint_p = sum(joint_p), .groups = "drop")
    
    df_lli <- df_estimate %>% 
      left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
      mutate(lli = log(joint_p)) %>% 
      pull(lli)
    
  }
  
  
  # browser()
  
  if(gradient) {
    
    
    df_grads <- df_probs_temp %>% 
      mutate(
        pi_help = (y1 == y1_star) + (y2 == y2_star) + (y3 == y3_star),
        d_probs_pi = case_when(
          pi_help == 0 ~ joint_p*(3/pi - 0/(1 - pi))*(1 - pi)*pi,
          pi_help == 1 ~ joint_p*(2/pi - 1/(1 - pi))*(1 - pi)*pi,
          pi_help == 2 ~ joint_p*(1/pi - 2/(1 - pi))*(1 - pi)*pi,
          pi_help == 3 ~ joint_p*(0/pi - 3/(1 - pi))*(1 - pi)*pi
        )
      ) %>% 
      mutate(
        theta_1_help0 = (y1_star == 1)*(y2_star == 1) + (y2_star == 1)*(y3_star == 1),
        theta_1_help1 = (y1_star == 1)*(y2_star == 0) + (y2_star == 1)*(y3_star == 0),
        theta_2_help0 = (y1_star == 0)*(y2_star == 0) + (y2_star == 0)*(y3_star == 0),
        theta_2_help1 = (y1_star == 0)*(y2_star == 1) + (y2_star == 0)*(y3_star == 1),
        d_probs_theta_1 = case_when(
          theta_1_help0 == 2                                ~ joint_p*(-2/(1 - param$theta_1)                   - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_1)*param$theta_1,
          theta_1_help0 == 1 & theta_1_help1 == 1           ~ joint_p*(-1/(1 - param$theta_1) + 1/param$theta_1 - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_1)*param$theta_1,
          theta_1_help0 == 1 & theta_1_help1 == 0           ~ joint_p*(-1/(1 - param$theta_1) + param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))*(1 - param$theta_1)*param$theta_1,
          theta_1_help1 == 1 & theta_1_help0 == 0 & y1_star == 1 ~ joint_p*(1/(param$theta_1) - 1/(param$theta_1 + param$theta_2))                                 *(1 - param$theta_1)*param$theta_1,
          theta_1_help1 == 1 & theta_1_help0 == 0 & y1_star == 0 ~ joint_p*(1/(param$theta_1) + param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))     *(1 - param$theta_1)*param$theta_1,
          theta_1_help1 == 0 & theta_1_help0 == 0           ~ joint_p*(param$theta_2/(param$theta_1*(param$theta_1 + param$theta_2)))                         *(1 - param$theta_1)*param$theta_1
        ),
        d_probs_theta_2 = case_when(
          theta_2_help0 == 2                                ~ joint_p*(-2/(1 - param$theta_2)                   - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_2)*param$theta_2,
          theta_2_help0 == 1 & theta_2_help1 == 1           ~ joint_p*(-1/(1 - param$theta_2) + 1/param$theta_2 - 1/(param$theta_1 + param$theta_2))          *(1 - param$theta_2)*param$theta_2,
          theta_2_help0 == 1 & theta_2_help1 == 0           ~ joint_p*(-1/(1 - param$theta_2) + param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))*(1 - param$theta_2)*param$theta_2,
          theta_2_help1 == 1 & theta_2_help0 == 0 & y1_star == 0 ~ joint_p*(1/(param$theta_2) - 1/(param$theta_1 + param$theta_2))                                 *(1 - param$theta_2)*param$theta_2,
          theta_2_help1 == 1 & theta_2_help0 == 0 & y1_star == 1 ~ joint_p*(1/(param$theta_2) + param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))     *(1 - param$theta_2)*param$theta_2,
          theta_2_help1 == 0 & theta_2_help0 == 0           ~ joint_p*(param$theta_1/(param$theta_2*(param$theta_1 + param$theta_2)))                         *(1 - param$theta_2)*param$theta_2
        )
      ) %>% 
      group_by(y1, y2, y3) %>% 
      summarise(
        joint_p = sum(joint_p),
        d_probs_theta_1 = sum(d_probs_theta_1),
        d_probs_theta_2 = sum(d_probs_theta_2),
        .groups = "drop"
      )
    
    df_lli <- df_estimate %>% 
      left_join(df_grads, by = c('y1', 'y2', 'y3')) %>% 
      mutate(
        lli = log(joint_p),
        lgi_theta_1 = d_probs_theta_1/joint_p,
        lgi_theta_2 = d_probs_theta_2/joint_p
      ) %>% 
      select(lli, lgi_theta_1, lgi_theta_2)
    
  }
  
  return(df_lli)
  
}

calc_mle_3waves_pi0 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_pi0(param_transformed, gradient = F))
  return(ll)
}

calc_mle_derivatives_3waves_pi0 <- function(param_transformed) {
  lg <- calc_lli_3waves_pi0(param_transformed, gradient = T) %>% select(-lli) %>% colSums
  return(lg)
}

