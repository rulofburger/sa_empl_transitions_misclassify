logit_transform <- function(param0) {
  param_input <- log(param0/(1 - param0))
  return(param_input)
}

logit_inverse <- function(param_input0) {
  param <- 1/(1 + exp(-param_input0))
  return(param)
}


calc_lli_3waves_ar1 <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, 
                        mu, 
                        1 - mu),
      p2_star = if_else(y1_star == 1, 
                        if_else(y2_star == 1, 
                                1 - param$theta_1,
                                param$theta_1), # theta_1 ==> job exit rate
                        if_else(y2_star == 1, 
                                param$theta_2,  # theta_2 ==> job entry rate 
                                1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, 
                        if_else(y3_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y3_star == 1,
                                param$theta_2, 
                                1 - param$theta_2)),
      p1 = if_else(y1 == y1_star, 
                   1 - param$pi, 
                   param$pi),
      p2 = if_else(y2 == y2_star, 
                   1 - param$pi, 
                   param$pi),
      p3 = if_else(y3 == y3_star, 
                   1 - param$pi, 
                   param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
  df_lli

}


calc_lli_derivatives_3waves_ar1 <- function(param_transformed, df_x) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, 
                        mu, 
                        1 - mu),
      p2_star = if_else(y1_star == 1, 
                        if_else(y2_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y2_star == 1, 
                                param$theta_2, 
                                1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, 
                        if_else(y3_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y3_star == 1, 
                                param$theta_2, 
                                1 - param$theta_2)),
      p1 = if_else(y1 == y1_star, 
                   1 - param$pi, 
                   param$pi),
      p2 = if_else(y2 == y2_star, 
                   1 - param$pi, 
                   param$pi),
      p3 = if_else(y3 == y3_star, 
                   1 - param$pi, 
                   param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta_1 = case_when(
        y1_star == 1 ~ -param$theta_2/((param$theta_1 + param$theta_2)^2),
        y1_star == 0 ~ param$theta_2/((param$theta_1 + param$theta_2)^2)
      ),
      d1_star_theta_2 = case_when(
        y1_star == 1 ~ -param$theta_2/((param$theta_1 + param$theta_2)^2) + 1/(param$theta_1 + param$theta_2),
        y1_star == 0 ~ param$theta_2/((param$theta_1 + param$theta_2)^2) - 1/(param$theta_1 + param$theta_2)
      ),
      d2_star_theta_1 = case_when(
        y2_star == 0 & y1_star == 0 ~ 0,
        y2_star == 0 & y1_star == 1 ~ 1,
        y2_star == 1 & y1_star == 0 ~ 0,
        y2_star == 1 & y1_star == 1 ~ -1
      ),
      d2_star_theta_2 = case_when(
        y2_star == 0 & y1_star == 0 ~ -1,
        y2_star == 0 & y1_star == 1 ~ 0,
        y2_star == 1 & y1_star == 0 ~ 1,
        y2_star == 1 & y1_star == 1 ~ 0
      ),
      d3_star_theta_1 = case_when(
        y3_star == 0 & y2_star == 0 ~ 0,
        y3_star == 0 & y2_star == 1 ~ 1,
        y3_star == 1 & y2_star == 0 ~ 0,
        y3_star == 1 & y2_star == 1 ~ -1
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
        p1_star*p1*p2_star*p2*p3_star*d3_pi
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_theta_2 = sum(joint_d_theta_2), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      .groups = "drop")
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_theta_2 = weight*joint_d_theta_2/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p) %>% 
    select(lgi_theta_1, lgi_theta_2, lgi_pi)

  df_gi
  
}
  



calc_mle_3waves_ar1 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1(param_transformed))
  return(lg)
}


calc_lli_3waves_ar1_pi0 <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  pi <- 0
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, 
                        mu, 
                        1 - mu),
      p2_star = if_else(y1_star == 1, 
                        if_else(y2_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y2_star == 1, 
                                param$theta_2, 
                                1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, 
                        if_else(y3_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1),
                        if_else(y3_star == 1, 
                                param$theta_2, 
                                1 - param$theta_2)),
      p1 = if_else(y1 == y1_star, 
                   1 - pi, 
                   pi),
      p2 = if_else(y2 == y2_star, 
                   1 - pi, 
                   pi),
      p3 = if_else(y3 == y3_star, 
                   1 - pi, 
                   pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    # group_by(y1, y2, y3, y4) %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    # left_join(df_probs, by = c('y1', 'y2', 'y3', 'y4')) %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
  df_lli
  
}

calc_lli_derivatives_3waves_ar1_pi0 <- function(param_transformed, df_x) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_2/(param$theta_1 + param$theta_2)
  pi <- 0
  
  df_probs_temp <- df_template %>% 
    mutate(
      p1_star = if_else(y1_star == 1, 
                        mu, 
                        1 - mu),
      p2_star = if_else(y1_star == 1, 
                        if_else(y2_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y2_star == 1, 
                                param$theta_2, 
                                1 - param$theta_2)),
      p3_star = if_else(y2_star == 1, 
                        if_else(y3_star == 1, 
                                1 - param$theta_1, 
                                param$theta_1), 
                        if_else(y3_star == 1, 
                                param$theta_2, 1
                                - param$theta_2)),
      p1 = if_else(y1 == y1_star, 
                   1 - pi, 
                   pi),
      p2 = if_else(y2 == y2_star, 
                   1 - pi, 
                   pi),
      p3 = if_else(y3 == y3_star, 
                   1 - pi, 
                   pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta_1 = case_when(
        y1_star == 1 ~ -param$theta_2/((param$theta_1 + param$theta_2)^2),
        y1_star == 0 ~ param$theta_2/((param$theta_1 + param$theta_2)^2)
      ),
      d1_star_theta_2 = case_when(
        y1_star == 1 ~ -param$theta_2/((param$theta_1 + param$theta_2)^2) + 1/(param$theta_1 + param$theta_2),
        y1_star == 0 ~ param$theta_2/((param$theta_1 + param$theta_2)^2) - 1/(param$theta_1 + param$theta_2)
      ),
      d2_star_theta_1 = case_when(
        y2_star == 0 & y1_star == 0 ~ 0,
        y2_star == 0 & y1_star == 1 ~ 1,
        y2_star == 1 & y1_star == 0 ~ 0,
        y2_star == 1 & y1_star == 1 ~ -1
      ),
      d2_star_theta_2 = case_when(
        y2_star == 0 & y1_star == 0 ~ -1,
        y2_star == 0 & y1_star == 1 ~ 0,
        y2_star == 1 & y1_star == 0 ~ 1,
        y2_star == 1 & y1_star == 1 ~ 0
      ),
      d3_star_theta_1 = case_when(
        y3_star == 0 & y2_star == 0 ~ 0,
        y3_star == 0 & y2_star == 1 ~ 1,
        y3_star == 1 & y2_star == 0 ~ 0,
        y3_star == 1 & y2_star == 1 ~ -1
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
        p1_star*p1*p2_star*p2*d3_star_theta_2*p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_theta_2 = sum(joint_d_theta_2),
      joint_p = sum(joint_p),
      .groups = "drop")
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, 
              by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_theta_2 = weight*joint_d_theta_2/joint_p,
    ) %>% 
    select(lgi_theta_1, lgi_theta_2)
  
  df_gi
}


calc_mle_3waves_ar1_pi0 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_pi0(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_pi0 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_pi0(param_transformed))
  return(lg)
}
