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

# WITH MISCLASSIFICATION ERROR ====

calc_lli_3waves_ar1 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  if(pi0) param$pi <- 0
  
  df_probs_temp <- df_template %>% 
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
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)

}


calc_lli_derivatives_3waves_ar1 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  if(pi0) param$pi <- 0
  
  df_probs_temp <- df_template %>% 
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
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
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
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0), 
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0), 
      joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p
      ) %>%
    select(lgi_theta_0, lgi_theta_1, lgi_pi)
  
  if(pi0) df_gi <- df_gi %>% select(lgi_theta_0, lgi_theta_1)
  return(df_gi)
  
}

calc_hessian_3waves_ar1 <- function(param_transformed, pi0 = FALSE) {
 
   param <- logit_inverse(param_transformed)
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  if(pi0) param$pi <- 0
  
  df_probs_temp <- df_template %>% 
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
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star * p1 * p2_star * p2 * p3_star * p3,
      
      # Derivatives for Hessian terms
      d1_star_theta_0 = case_when(
        y1_star == 0 ~ -param$theta_1 / ((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ param$theta_1 / ((param$theta_1 + param$theta_0)^2)
      ),
      d1_star_theta_1 = case_when(
        y1_star == 0 ~ param$theta_0 / ((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ -param$theta_0 / ((param$theta_1 + param$theta_0)^2)
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
  
  # Second derivatives (Hessian terms)
  df_hessian_temp <- df_probs_temp %>% 
    mutate(
      d1_star_theta_0_theta_0 = case_when(
        y1_star == 0 ~ 2*param$theta_1/((param$theta_1 + param$theta_0)^3),
        y1_star == 1 ~ -2*param$theta_1/((param$theta_1 + param$theta_0)^3)
      ),
      d1_star_theta_1_theta_1 = case_when(
        y1_star == 0 ~ -2*param$theta_0/((param$theta_1 + param$theta_0)^3),
        y1_star == 1 ~ 2*param$theta_0/((param$theta_1 + param$theta_0)^3)
      ),
      d1_star_theta_0_theta_1 = case_when(
        y1_star == 0 ~ (param$theta_1 - param$theta_0)/((param$theta_1 + param$theta_0)^3),
        y1_star == 1 ~ (param$theta_0 - param$theta_1)/((param$theta_1 + param$theta_0)^3)
      ),
      d1_star_theta_0_scaled = d1_star_theta_0/p1_star,
      d1_star_theta_1_scaled = d1_star_theta_1/p1_star,
      d2_star_theta_0_scaled = d2_star_theta_0/p2_star,
      d2_star_theta_1_scaled = d2_star_theta_1/p2_star,
      d3_star_theta_0_scaled = d3_star_theta_0/p3_star,
      d3_star_theta_1_scaled = d3_star_theta_1/p3_star,
      d1_pi_scaled = d1_pi/p1,
      d2_pi_scaled = d2_pi/p2,
      d3_pi_scaled = d3_pi/p3,
      d1_star_theta_0_theta_0_scaled = d1_star_theta_0_theta_0/p1_star,
      d1_star_theta_1_theta_1_scaled = d1_star_theta_1_theta_1/p1_star,
      d1_star_theta_0_theta_1_scaled = d1_star_theta_0_theta_1/p1_star,
      joint_d2_theta_0_theta_0 = joint_p*(
        d1_star_theta_0_theta_0_scaled +
          2*d1_star_theta_0_scaled*d2_star_theta_0_scaled + 
          2*d1_star_theta_0_scaled*d3_star_theta_0_scaled + 
          2*d2_star_theta_0_scaled*d3_star_theta_0_scaled
      ),
      joint_d2_theta_1_theta_1 = joint_p*(
        d1_star_theta_1_theta_1_scaled +
          2*d1_star_theta_1_scaled*d2_star_theta_1_scaled + 
          2*d1_star_theta_1_scaled*d3_star_theta_1_scaled + 
          2*d2_star_theta_1_scaled*d3_star_theta_1_scaled
      ),
      joint_d2_pi_pi = joint_p*(
        2*d1_pi_scaled*d2_pi_scaled + 
          2*d1_pi_scaled*d3_pi_scaled + 
          2*d2_pi_scaled*d3_pi_scaled
      ),
      joint_d2_theta_0_theta_1 = joint_p*(
        d1_star_theta_0_theta_1_scaled +
          d1_star_theta_0_scaled*(d2_star_theta_1_scaled + d3_star_theta_1_scaled) + 
          d2_star_theta_0_scaled*(d1_star_theta_1_scaled + d3_star_theta_1_scaled) + 
          d3_star_theta_0_scaled*(d1_star_theta_1_scaled + d2_star_theta_1_scaled)
      ),
      joint_d2_theta_0_pi = joint_p*(
        d1_star_theta_0_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta_0_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta_0_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      ),
      joint_d2_theta_1_pi = joint_p*(
        d1_star_theta_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      )
    ) 

  df_hessian <- df_hessian_temp %>% 
    group_by(y1, y2, y3) %>%
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0),
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_pi = sum(joint_d_pi),
      joint_d2_theta_0_theta_0 = sum(joint_d2_theta_0_theta_0),
      joint_d2_theta_1_theta_1 = sum(joint_d2_theta_1_theta_1),
      joint_d2_theta_0_theta_1 = sum(joint_d2_theta_0_theta_1),
      joint_d2_theta_0_pi = sum(joint_d2_theta_0_pi),
      joint_d2_theta_1_pi = sum(joint_d2_theta_1_pi),
      joint_d2_pi_pi = sum(joint_d2_pi_pi),
      joint_p = sum(joint_p),
      .groups = "drop"
    ) %>% 
    mutate(
      joint_d2_theta_0_theta_0 = joint_d2_theta_0_theta_0*(param$theta_0^2)*((1 - param$theta_0)^2) + joint_d_theta_0*(param$theta_0)*(1 - param$theta_0)*(1 - 2*param$theta_0),
      joint_d2_theta_1_theta_1 = joint_d2_theta_1_theta_1*(param$theta_0^2)*((1 - param$theta_0)^2) + joint_d_theta_1*(param$theta_1)*(1 - param$theta_1)*(1 - 2*param$theta_1),
      joint_d2_pi_pi = joint_d2_pi_pi*(param$pi^2)*((1 - param$pi)^2) + joint_d_pi*(param$pi)*(1 - param$pi)*(1 - 2*param$pi),
      joint_d2_theta_0_theta_1 = joint_d2_theta_0_theta_1*(param$theta_1)*(1 - param$theta_1)*(param$theta_0)*(1 - param$theta_0),
      joint_d2_theta_0_pi = joint_d2_theta_0_pi*(param$theta_0)*(1 - param$theta_0)*(param$pi)*(1 - param$pi),
      joint_d2_theta_1_pi = joint_d2_theta_1_pi*(param$theta_1)*(1 - param$theta_1)*(param$pi)*(1 - param$pi),
      joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0), 
      joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi),
      joint_p = joint_p
    )
  
  df_hi <- df_estimate %>% 
    left_join(df_hessian, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lhi_theta_0_theta_0 = weight*(joint_d2_theta_0_theta_0/joint_p - joint_d_theta_0*joint_d_theta_0/(joint_p^2)),
      lhi_theta_1_theta_1 = weight*(joint_d2_theta_1_theta_1/joint_p - joint_d_theta_1*joint_d_theta_1/(joint_p^2)),
      lhi_pi_pi = weight*(joint_d2_pi_pi/joint_p - joint_d_pi*joint_d_pi/(joint_p^2)),
      lhi_theta_0_theta_1 = weight*(joint_d2_theta_0_theta_1/joint_p - joint_d_theta_0*joint_d_theta_1/(joint_p^2)),
      lhi_theta_0_pi = weight*(joint_d2_theta_0_pi/joint_p - joint_d_theta_0*joint_d_pi/(joint_p^2)),
      lhi_theta_1_pi = weight*(joint_d2_theta_1_pi/joint_p - joint_d_theta_1*joint_d_pi/(joint_p^2))
    ) %>%
    select(lhi_theta_0_theta_0, lhi_theta_1_theta_1, lhi_pi_pi, lhi_theta_0_theta_1, lhi_theta_0_pi, lhi_theta_1_pi)
  
  
  # Convert to matrix
  hessian_matrix <- matrix(c(
    sum(df_hi$lhi_theta_0_theta_0),
    sum(df_hi$lhi_theta_0_theta_1),
    sum(df_hi$lhi_theta_0_pi),
    sum(df_hi$lhi_theta_0_theta_1),
    sum(df_hi$lhi_theta_1_theta_1),
    sum(df_hi$lhi_theta_1_pi),
    sum(df_hi$lhi_theta_0_pi),
    sum(df_hi$lhi_theta_1_pi),
    sum(df_hi$lhi_pi_pi)
  ), nrow = 3, ncol = 3, byrow = TRUE)
  
  if(pi0) hessian_matrix <- hessian_matrix[1:2, 1:2] # Remove terms involving `pi` if `pi0 = TRUE`
  
  return(hessian_matrix)
}
  
calc_mle_3waves_ar1 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1(param_transformed))
  return(lg)
}

# WITH NO MISCLASSIFICATION ERROR ====

calc_lli_3waves_ar1_pi0 <- function(param_transformed) {
  calc_lli_3waves_ar1(param_transformed, pi0 = TRUE)
}

calc_lli_derivatives_3waves_ar1_pi0 <- function(param_transformed) {
  calc_lli_derivatives_3waves_ar1(param_transformed, pi0 = TRUE)
}

# 
# calc_lli_3waves_ar1_pi0 <- function(param_transformed) {
#   
#   param <- logit_inverse(param_transformed)
#   mu <- param$theta_0/(param$theta_1 + param$theta_0)
#   pi <- 0
#   
#   df_probs_temp <- df_template %>% 
#     mutate(
#       p1_star = if_else(y1_star == 1, mu, 1 - mu),
#       p2_star = case_when(
#         y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
#         y1_star == 0 & y2_star == 1 ~ param$theta_0,
#         y1_star == 1 & y2_star == 0 ~ param$theta_1,
#         y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
#       ),
#       p3_star = case_when(
#         y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
#         y2_star == 0 & y3_star == 1 ~ param$theta_0,
#         y2_star == 1 & y3_star == 0 ~ param$theta_1,
#         y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
#       ),      p1 = if_else(y1 == y1_star, 1 - pi, pi),
#       p2 = if_else(y2 == y2_star, 1 - pi, pi),
#       p3 = if_else(y3 == y3_star, 1 - pi, pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   df_probs <- df_probs_temp %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p = sum(joint_p), .groups = "drop")
#   
#   df_lli <- df_estimate %>% 
#     left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(lli = weight*log(joint_p)) %>%
#     pull(lli)
#   
# }
# 
# calc_lli_derivatives_3waves_ar1_pi0 <- function(param_transformed) {
#   
#   param <- logit_inverse(param_transformed)
#   mu <- param$theta_0/(param$theta_1 + param$theta_0)
#   pi <- 0
#   
#   df_probs_temp <- df_template %>% 
#     mutate(
#       p1_star = if_else(y1_star == 1, mu, 1 - mu),
#       p2_star = case_when(
#         y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
#         y1_star == 0 & y2_star == 1 ~ param$theta_0,
#         y1_star == 1 & y2_star == 0 ~ param$theta_1,
#         y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
#       ),
#       p3_star = case_when(
#         y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
#         y2_star == 0 & y3_star == 1 ~ param$theta_0,
#         y2_star == 1 & y3_star == 0 ~ param$theta_1,
#         y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
#       ),
#       p1 = if_else(y1 == y1_star, 1 - pi, pi),
#       p2 = if_else(y2 == y2_star, 1 - pi, pi),
#       p3 = if_else(y3 == y3_star, 1 - pi, pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta_0 = case_when(
#         y1_star == 0 ~ -param$theta_1/((param$theta_1 + param$theta_0)^2),
#         y1_star == 1 ~ param$theta_1/((param$theta_1 + param$theta_0)^2)
#       ),
#       d1_star_theta_1 = case_when(
#         y1_star == 0 ~ param$theta_0/((param$theta_1 + param$theta_0)^2),
#         y1_star == 1 ~ -param$theta_0/((param$theta_1 + param$theta_0)^2)
#       ),
#       d2_star_theta_0 = case_when(
#         y1_star == 0 & y2_star == 0 ~ -1,
#         y1_star == 0 & y2_star == 1 ~ 1,
#         y1_star == 1 & y2_star == 0 ~ 0,
#         y1_star == 1 & y2_star == 1 ~ 0
#       ),
#       d2_star_theta_1 = case_when(
#         y1_star == 0 & y2_star == 0 ~ 0,
#         y1_star == 0 & y2_star == 1 ~ 0,
#         y1_star == 1 & y2_star == 0 ~ 1,
#         y1_star == 1 & y2_star == 1 ~ -1
#       ),
#       d3_star_theta_0 = case_when(
#         y2_star == 0 & y3_star == 0 ~ -1,
#         y2_star == 0 & y3_star == 1 ~ 1,
#         y2_star == 1 & y3_star == 0 ~ 0,
#         y2_star == 1 & y3_star == 1 ~ 0
#       ),
#       d3_star_theta_1 = case_when(
#         y2_star == 0 & y3_star == 0 ~ 0,
#         y2_star == 0 & y3_star == 1 ~ 0,
#         y2_star == 1 & y3_star == 0 ~ 1,
#         y2_star == 1 & y3_star == 1 ~ -1
#       ),
#       d1_pi = if_else(y1 == y1_star, -1, 1),
#       d2_pi = if_else(y2 == y2_star, -1, 1),
#       d3_pi = if_else(y3 == y3_star, -1, 1),
#       joint_d_theta_0 = 
#         d1_star_theta_0*joint_p/p1_star + 
#         d2_star_theta_0*joint_p/p2_star + 
#         d3_star_theta_0*joint_p/p3_star,
#         joint_d_theta_1 = 
#         d1_star_theta_1*joint_p/p1_star + 
#         d2_star_theta_1*joint_p/p2_star + 
#         d3_star_theta_1*joint_p/p3_star
#     ) 
#   
#   df_grad <- df_probs_temp %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta_0 = sum(joint_d_theta_0),
#       joint_d_theta_1 = sum(joint_d_theta_1),
#       joint_p = sum(joint_p),
#       .groups = "drop") %>% 
#     mutate(
#       joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0),
#       joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1)
#     )
#   
#   df_gi <- df_estimate %>% 
#     left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(
#       lgi_theta_0 = weight*joint_d_theta_0/joint_p,
#       lgi_theta_1 = weight*joint_d_theta_1/joint_p
#     ) %>% 
#     select(lgi_theta_0, lgi_theta_1)
# }


calc_mle_3waves_ar1_pi0 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_pi0(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_pi0 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_pi0(param_transformed))
  return(lg)
}


# ASYMMETRIC MISCLASSIFICATION ERROR ====

calc_lli_3waves_ar1_asymmetric <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  
  df_probs_temp <- df_template %>% 
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
      p1 = case_when(
        y1 == 0 & y1_star == 0 ~ 1 - param$pi_0,
        y1 == 1 & y1_star == 0 ~ param$pi_0,
        y1 == 0 & y1_star == 1 ~ param$pi_1,
        y1 == 1 & y1_star == 1 ~ 1 - param$pi_1
      ),
      p2 = case_when(
        y2 == 0 & y2_star == 0 ~ 1 - param$pi_0,
        y2 == 1 & y2_star == 0 ~ param$pi_0,
        y2 == 0 & y2_star == 1 ~ param$pi_1,
        y2 == 1 & y2_star == 1 ~ 1 - param$pi_1
      ),
      p3 = case_when(
        y3 == 0 & y3_star == 0 ~ 1 - param$pi_0,
        y3 == 1 & y3_star == 0 ~ param$pi_0,
        y3 == 0 & y3_star == 1 ~ param$pi_1,
        y3 == 1 & y3_star == 1 ~ 1 - param$pi_1
      ),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}


calc_lli_derivatives_3waves_ar1_asymmetric <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  
  df_probs_temp <- df_template %>% 
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
      p1 = case_when(
        y1 == 0 & y1_star == 0 ~ 1 - param$pi_0,
        y1 == 1 & y1_star == 0 ~ param$pi_0,
        y1 == 0 & y1_star == 1 ~ param$pi_1,
        y1 == 1 & y1_star == 1 ~ 1 - param$pi_1
      ),
      p2 = case_when(
        y2 == 0 & y2_star == 0 ~ 1 - param$pi_0,
        y2 == 1 & y2_star == 0 ~ param$pi_0,
        y2 == 0 & y2_star == 1 ~ param$pi_1,
        y2 == 1 & y2_star == 1 ~ 1 - param$pi_1
      ),
      p3 = case_when(
        y3 == 0 & y3_star == 0 ~ 1 - param$pi_0,
        y3 == 1 & y3_star == 0 ~ param$pi_0,
        y3 == 0 & y3_star == 1 ~ param$pi_1,
        y3 == 1 & y3_star == 1 ~ 1 - param$pi_1
      ),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
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
      d1_pi_0 = case_when(
        y1 == 0 & y1_star == 0 ~ -1,
        y1 == 1 & y1_star == 0 ~ 1,
        y1 == 0 & y1_star == 1 ~ 0,
        y1 == 1 & y1_star == 1 ~ 0
      ),
      d2_pi_0 = case_when(
        y2 == 0 & y2_star == 0 ~ -1,
        y2 == 1 & y2_star == 0 ~ 1,
        y2 == 0 & y2_star == 1 ~ 0,
        y2 == 1 & y2_star == 1 ~ 0
      ),
      d3_pi_0 = case_when(
        y3 == 0 & y3_star == 0 ~ -1,
        y3 == 1 & y3_star == 0 ~ 1,
        y3 == 0 & y3_star == 1 ~ 0,
        y3 == 1 & y3_star == 1 ~ 0
      ),
      d1_pi_1 = case_when(
        y1 == 0 & y1_star == 0 ~ 0,
        y1 == 1 & y1_star == 0 ~ 0,
        y1 == 0 & y1_star == 1 ~ 1,
        y1 == 1 & y1_star == 1 ~ -1
      ),
      d2_pi_1 = case_when(
        y2 == 0 & y2_star == 0 ~ 0,
        y2 == 1 & y2_star == 0 ~ 0,
        y2 == 0 & y2_star == 1 ~ 1,
        y2 == 1 & y2_star == 1 ~ -1
      ),
      d3_pi_1 = case_when(
        y3 == 0 & y3_star == 0 ~ 0,
        y3 == 1 & y3_star == 0 ~ 0,
        y3 == 0 & y3_star == 1 ~ 1,
        y3 == 1 & y3_star == 1 ~ -1
      ),
      joint_d_theta_0 = 
        d1_star_theta_0*joint_p/p1_star + 
        d2_star_theta_0*joint_p/p2_star + 
        d3_star_theta_0*joint_p/p3_star,
      joint_d_theta_1 = 
        d1_star_theta_1*joint_p/p1_star + 
        d2_star_theta_1*joint_p/p2_star + 
        d3_star_theta_1*joint_p/p3_star,
      joint_d_pi_0 = 
        d1_pi_0*joint_p/p1 + 
        d2_pi_0*joint_p/p2 + 
        d3_pi_0*joint_p/p3,
      joint_d_pi_1 = 
        d1_pi_1*joint_p/p1 + 
        d2_pi_1*joint_p/p2 + 
        d3_pi_1*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0), 
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_pi_0 = sum(joint_d_pi_0),
      joint_d_pi_1 = sum(joint_d_pi_1),
      joint_p = sum(joint_p),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0), 
      joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1),
      joint_d_pi_0 = joint_d_pi_0*param$pi_0*(1 - param$pi_0),
      joint_d_pi_1 = joint_d_pi_1*param$pi_1*(1 - param$pi_1)
    )
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_pi_0 = weight*joint_d_pi_0/joint_p,
      lgi_pi_1 = weight*joint_d_pi_1/joint_p
    ) %>%
    select(lgi_theta_0, lgi_theta_1, lgi_pi_0, lgi_pi_1)
}




calc_mle_3waves_ar1_asymmetric <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar1_asymmetric(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar1_asymmetric <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar1_asymmetric(param_transformed))
  return(lg)
}


# AR(2) ESTIMATOR OVER 3 WAVES WITH NO MISCLASSIFICATION ERROR ====

calc_lli_3waves_ar2_pi0 <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  pi <- 0
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
      p1 = if_else(y1 == y1_star, 1 - pi, pi),
      p2 = if_else(y2 == y2_star, 1 - pi, pi),
      p3 = if_else(y3 == y3_star, 1 - pi, pi),
      joint_p = p12_star*p1*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop")
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
}



calc_lli_derivatives_3waves_ar2_pi0 <- function(param_transformed) {
  
  param <- logit_inverse(param_transformed)
  pi <- 0
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
      p1 = if_else(y1 == y1_star, 1 - pi, pi),
      p2 = if_else(y2 == y2_star, 1 - pi, pi),
      p3 = if_else(y3 == y3_star, 1 - pi, pi),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      joint_p = p12_star*p1*p2*p3_star*p3,
      joint_d_theta_0 = 
        d12_star_theta_0*p3_star*p1*p2*p3 +
        p12_star*d3_star_theta_0*p1*p2*p3 + 
        p12_star*p3_star*p1*p2*p3,
      joint_d_theta_1 = 
        d12_star_theta_1*p3_star*p1*p2*p3 +
        p12_star*d3_star_theta_1*p1*p2*p3 + 
        p12_star*p3_star*p1*p2*p3,
      joint_d_theta_01 = 
        d12_star_theta_01*p3_star*p1*p2*p3 +
        p12_star*d3_star_theta_01*p1*p2*p3 + 
        p12_star*p3_star*p1*p2*p3,
      joint_d_theta_10 = 
        d12_star_theta_10*p3_star*p1*p2*p3 +
        p12_star*d3_star_theta_10*p1*p2*p3 + 
        p12_star*p3_star*p1*p2*p3,
      joint_d_pi = 
        p12_star*p3_star*d1_pi*p2*p3 + 
        p12_star*p3_star*p1*d2_pi*p3 + 
        p12_star*p3_star*p1*p2*d3_pi + 
        p12_star*p3_star*p1*p2*p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0),
      joint_d_theta_1 = sum(joint_d_theta_1), 
      joint_d_theta_01 = sum(joint_d_theta_01), 
      joint_d_theta_10 = sum(joint_d_theta_10), 
      # joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      .groups = "drop")
  
  
  df_gi <- df_estimate %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_01 = weight*joint_d_theta_01/joint_p,
      lgi_theta_10 = weight*joint_d_theta_10/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      # lgi_pi = weight*joint_d_pi/joint_p
      # lgi_theta_0 = joint_d_theta_0/joint_p,
      # lgi_theta_01 =joint_d_theta_01/joint_p,
      # lgi_theta_10 = joint_d_theta_10/joint_p,
      # lgi_theta_1 = joint_d_theta_1/joint_p,
      # lgi_pi = joint_d_pi/joint_p
    ) %>% 
    select(lgi_theta_0, lgi_theta_01, lgi_theta_10, lgi_theta_1)
  
}

calc_mle_3waves_ar2_pi0 <- function(param_transformed) {
  ll <- sum(calc_lli_3waves_ar2_pi0(param_transformed))
  return(ll)
}

calc_mle_derivatives_3waves_ar2_pi0 <- function(param_transformed) {
  lg <- colSums(calc_lli_derivatives_3waves_ar2_pi0(param_transformed))
  return(lg)
}
