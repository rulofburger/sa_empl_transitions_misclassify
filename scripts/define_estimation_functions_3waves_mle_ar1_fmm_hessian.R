# DEFINE FUNCTIONS ====

#> General functions that restrict parameter values to unit interval

# logit_transform <- function(param0) {
#   param_input <- log(param0/(1 - param0))
#   return(param_input)
# }
# 
# logit_inverse <- function(param_input0) {
#   param <- 1/(1 + exp(-param_input0))
#   return(param)
# }

# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH 2 FMM groups ====
# 
# calc_lli_3waves_ar1_fmm2 <- function(param_transformed, pi0 = FALSE) {
#   
#   param <- logit_inverse(param_transformed)
#   mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
#   mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
#   if(pi0) param$pi <- 0
#   
#   df_probs_temp1 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_1,
#         y1_star == 0 & y2_star == 1, param$theta0_1,
#         y1_star == 1 & y2_star == 0, param$theta1_1,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_1
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_1,
#         y2_star == 0 & y3_star == 1, param$theta0_1,
#         y2_star == 1 & y3_star == 0, param$theta1_1,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_1
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   df_probs_temp2 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_2, 1 - mu_2),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_2,
#         y1_star == 0 & y2_star == 1, param$theta0_2,
#         y1_star == 1 & y2_star == 0, param$theta1_2,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_2
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_2,
#         y2_star == 0 & y3_star == 1, param$theta0_2,
#         y2_star == 1 & y3_star == 0, param$theta1_2,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_2
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   df_probs1 <- df_probs_temp1 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p1 = sum(joint_p), .groups = "drop")
#   
#   df_probs2 <- df_probs_temp2 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p2 = sum(joint_p), .groups = "drop")
#   
#   df_probs <- df_probs1 %>% 
#     left_join(df_probs2, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(joint_p = param$p_1*joint_p1 + (1 - param$p_1)*joint_p2) %>% 
#     mutate(joint_p = pmax(joint_p, 1e-10))
#   
#   df_lli <- df_estimate %>% 
#     left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(lli = weight*log(joint_p)) %>%
#     pull(lli)
#   
# }

# 
# calc_lli_derivatives_3waves_ar1_fmm2_bla <- function(param_transformed, pi0 = FALSE) {
#   
#   param <- logit_inverse(param_transformed)
#   mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
#   mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
#   if(pi0) param$pi <- 0
#   
#   df_probs_temp1 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_1,
#         y1_star == 0 & y2_star == 1, param$theta0_1,
#         y1_star == 1 & y2_star == 0, param$theta1_1,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_1
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_1,
#         y2_star == 0 & y3_star == 1, param$theta0_1,
#         y2_star == 1 & y3_star == 0, param$theta1_1,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_1
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta0_1 = fcase(
#         y1_star == 0, -param$theta1_1/((param$theta1_1 + param$theta0_1)^2),
#         y1_star == 1, param$theta1_1/((param$theta1_1 + param$theta0_1)^2)
#       ),
#       d1_star_theta1_1 = fcase(
#         y1_star == 0, param$theta0_1/((param$theta1_1 + param$theta0_1)^2),
#         y1_star == 1, -param$theta0_1/((param$theta1_1 + param$theta0_1)^2)
#       ),
#       d2_star_theta0_1 = fcase(
#         y1_star == 0 & y2_star == 0, -1,
#         y1_star == 0 & y2_star == 1, 1,
#         y1_star == 1 & y2_star == 0, 0,
#         y1_star == 1 & y2_star == 1, 0
#       ),
#       d2_star_theta1_1 = fcase(
#         y1_star == 0 & y2_star == 0, 0,
#         y1_star == 0 & y2_star == 1, 0,
#         y1_star == 1 & y2_star == 0, 1,
#         y1_star == 1 & y2_star == 1, -1
#       ),
#       d3_star_theta0_1 = fcase(
#         y2_star == 0 & y3_star == 0, -1,
#         y2_star == 0 & y3_star == 1, 1,
#         y2_star == 1 & y3_star == 0, 0,
#         y2_star == 1 & y3_star == 1, 0
#       ),
#       d3_star_theta1_1 = fcase(
#         y2_star == 0 & y3_star == 0, 0,
#         y2_star == 0 & y3_star == 1, 0,
#         y2_star == 1 & y3_star == 0, 1,
#         y2_star == 1 & y3_star == 1, -1
#       ),
#       d1_pi = fifelse(y1 == y1_star, -1, 1),
#       d2_pi = fifelse(y2 == y2_star, -1, 1),
#       d3_pi = fifelse(y3 == y3_star, -1, 1),
#       joint_d_theta0_1 = 
#         d1_star_theta0_1*joint_p/p1_star + 
#         d2_star_theta0_1*joint_p/p2_star + 
#         d3_star_theta0_1*joint_p/p3_star,
#       joint_d_theta1_1 = 
#         d1_star_theta1_1*joint_p/p1_star + 
#         d2_star_theta1_1*joint_p/p2_star + 
#         d3_star_theta1_1*joint_p/p3_star,
#       joint_d_pi = 
#         d1_pi*joint_p/p1 + 
#         d2_pi*joint_p/p2 + 
#         d3_pi*joint_p/p3
#     ) 
#   
#   df_probs_temp2 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_2, 1 - mu_2),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_2,
#         y1_star == 0 & y2_star == 1, param$theta0_2,
#         y1_star == 1 & y2_star == 0, param$theta1_2,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_2
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_2,
#         y2_star == 0 & y3_star == 1, param$theta0_2,
#         y2_star == 1 & y3_star == 0, param$theta1_2,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_2
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta0_2 = fcase(
#         y1_star == 0, -param$theta1_2/((param$theta1_2 + param$theta0_2)^2),
#         y1_star == 1, param$theta1_2/((param$theta1_2 + param$theta0_2)^2)
#       ),
#       d1_star_theta1_2 = fcase(
#         y1_star == 0, param$theta0_2/((param$theta1_2 + param$theta0_2)^2),
#         y1_star == 1, -param$theta0_2/((param$theta1_2 + param$theta0_2)^2)
#       ),
#       d2_star_theta0_2 = fcase(
#         y1_star == 0 & y2_star == 0, -1,
#         y1_star == 0 & y2_star == 1, 1,
#         y1_star == 1 & y2_star == 0, 0,
#         y1_star == 1 & y2_star == 1, 0
#       ),
#       d2_star_theta1_2 = fcase(
#         y1_star == 0 & y2_star == 0, 0,
#         y1_star == 0 & y2_star == 1, 0,
#         y1_star == 1 & y2_star == 0, 1,
#         y1_star == 1 & y2_star == 1, -1
#       ),
#       d3_star_theta0_2 = fcase(
#         y2_star == 0 & y3_star == 0, -1,
#         y2_star == 0 & y3_star == 1, 1,
#         y2_star == 1 & y3_star == 0, 0,
#         y2_star == 1 & y3_star == 1, 0
#       ),
#       d3_star_theta1_2 = fcase(
#         y2_star == 0 & y3_star == 0, 0,
#         y2_star == 0 & y3_star == 1, 0,
#         y2_star == 1 & y3_star == 0, 1,
#         y2_star == 1 & y3_star == 1, -1
#       ),
#       d1_pi = fifelse(y1 == y1_star, -1, 1),
#       d2_pi = fifelse(y2 == y2_star, -1, 1),
#       d3_pi = fifelse(y3 == y3_star, -1, 1),
#       joint_d_theta0_2 = 
#         d1_star_theta0_2*joint_p/p1_star + 
#         d2_star_theta0_2*joint_p/p2_star + 
#         d3_star_theta0_2*joint_p/p3_star,
#       joint_d_theta1_2 = 
#         d1_star_theta1_2*joint_p/p1_star + 
#         d2_star_theta1_2*joint_p/p2_star + 
#         d3_star_theta1_2*joint_p/p3_star,
#       joint_d_pi = 
#         d1_pi*joint_p/p1 + 
#         d2_pi*joint_p/p2 + 
#         d3_pi*joint_p/p3
#     )  
#   
#   df_grad1 <- df_probs_temp1 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta0_1 = sum(joint_d_theta0_1), 
#       joint_d_theta1_1 = sum(joint_d_theta1_1),
#       joint_d_pi1 = sum(joint_d_pi),
#       joint_p1 = sum(joint_p),
#       .groups = "drop") 
#   
#   df_grad2 <- df_probs_temp2 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta0_2 = sum(joint_d_theta0_2), 
#       joint_d_theta1_2 = sum(joint_d_theta1_2),
#       joint_d_pi2 = sum(joint_d_pi),
#       joint_p2 = sum(joint_p),
#       .groups = "drop") 
#   
#   df_grad <- df_grad1 %>% 
#     left_join(df_grad2, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(
#       joint_p = param$p_1*joint_p1 + (1 - param$p_1)*joint_p2,
#       joint_p = pmax(joint_p, 1e-10),
#       joint_d_pi = param$p_1*joint_d_pi1 + (1 - param$p_1)*joint_d_pi2
#       ) %>% 
#     mutate(
#       joint_d_theta0_1 = param$p_1*joint_d_theta0_1*param$theta0_1*(1 - param$theta0_1), 
#       joint_d_theta1_1 = param$p_1*joint_d_theta1_1*param$theta1_1*(1 - param$theta1_1),
#       joint_d_theta0_2 = (1 - param$p_1)*joint_d_theta0_2*param$theta0_2*(1 - param$theta0_2), 
#       joint_d_theta1_2 = (1 - param$p_1)*joint_d_theta1_2*param$theta1_2*(1 - param$theta1_2),
#       joint_d_p_1 = joint_p1*param$p_1*(1 - param$p_1) - joint_p2*param$p_1*(1 - param$p_1),
#       joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
#     )
#   
#   df_gi <- df_estimate %>% 
#     left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(
#       lgi_theta0_1 = weight*joint_d_theta0_1/joint_p,
#       lgi_theta1_1 = weight*joint_d_theta1_1/joint_p,
#       lgi_theta0_2 = weight*joint_d_theta0_2/joint_p,
#       lgi_theta1_2 = weight*joint_d_theta1_2/joint_p,
#       lgi_p_1 = weight*joint_d_p_1/joint_p,
#       lgi_pi = weight*joint_d_pi/joint_p
#     ) %>%
#     select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_pi)
#   
#   if(pi0) df_gi <- df_gi %>% select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2)
#   return(df_gi)
#   
# }

calc_hessian_3waves_ar1_fmm2 <- function(param_transformed, pi0 = FALSE) {
  
  param <- logit_inverse(param_transformed)
  mu_1  <- param$theta0_1/(param$theta0_1 + param$theta1_1)
  mu_2  <- param$theta0_2/(param$theta0_2 + param$theta1_2)
  if (pi0) param$pi <- 0
  
  df_probs_temp1 <- df_template %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - param$theta0_1,
        y1_star == 0 & y2_star == 1, param$theta0_1,
        y1_star == 1 & y2_star == 0, param$theta1_1,
        y1_star == 1 & y2_star == 1, 1 - param$theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - param$theta0_1,
        y2_star == 0 & y3_star == 1, param$theta0_1,
        y2_star == 1 & y3_star == 0, param$theta1_1,
        y2_star == 1 & y3_star == 1, 1 - param$theta1_1
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_1 = fcase(
        y1_star == 0,  -param$theta1_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1, param$theta1_1/((param$theta1_1 + param$theta0_1)^2)
      ),
      d1_star_theta1_1 = fcase(
        y1_star == 0, param$theta0_1/((param$theta1_1 + param$theta0_1)^2),
        y1_star == 1, -param$theta0_1/((param$theta1_1 + param$theta0_1)^2)
      ),
      d2_star_theta0_1 = fcase(
        y1_star == 0 & y2_star == 0, -1,
        y1_star == 0 & y2_star == 1, 1,
        y1_star == 1 & y2_star == 0, 0,
        y1_star == 1 & y2_star == 1, 0
      ),
      d2_star_theta1_1 = fcase(
        y1_star == 0 & y2_star == 0, 0,
        y1_star == 0 & y2_star == 1, 0,
        y1_star == 1 & y2_star == 0, 1,
        y1_star == 1 & y2_star == 1, -1
      ),
      d3_star_theta0_1 = fcase(
        y2_star == 0 & y3_star == 0, -1,
        y2_star == 0 & y3_star == 1, 1,
        y2_star == 1 & y3_star == 0, 0,
        y2_star == 1 & y3_star == 1, 0
      ),
      d3_star_theta1_1 = fcase(
        y2_star == 0 & y3_star == 0, 0,
        y2_star == 0 & y3_star == 1, 0,
        y2_star == 1 & y3_star == 0, 1,
        y2_star == 1 & y3_star == 1, -1
      ),
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
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
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_2, 1 - mu_2),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - param$theta0_2,
        y1_star == 0 & y2_star == 1, param$theta0_2,
        y1_star == 1 & y2_star == 0, param$theta1_2,
        y1_star == 1 & y2_star == 1, 1 - param$theta1_2
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - param$theta0_2,
        y2_star == 0 & y3_star == 1, param$theta0_2,
        y2_star == 1 & y3_star == 0, param$theta1_2,
        y2_star == 1 & y3_star == 1, 1 - param$theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
      d1_star_theta0_2 = fcase(
        y1_star == 0, -param$theta1_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1, param$theta1_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d1_star_theta1_2 = fcase(
        y1_star == 0, param$theta0_2/((param$theta1_2 + param$theta0_2)^2),
        y1_star == 1, -param$theta0_2/((param$theta1_2 + param$theta0_2)^2)
      ),
      d2_star_theta0_2 = fcase(
        y1_star == 0 & y2_star == 0, -1,
        y1_star == 0 & y2_star == 1, 1,
        y1_star == 1 & y2_star == 0, 0,
        y1_star == 1 & y2_star == 1, 0
      ),
      d2_star_theta1_2 = fcase(
        y1_star == 0 & y2_star == 0, 0,
        y1_star == 0 & y2_star == 1, 0,
        y1_star == 1 & y2_star == 0, 1,
        y1_star == 1 & y2_star == 1, -1
      ),
      d3_star_theta0_2 = fcase(
        y2_star == 0 & y3_star == 0, -1,
        y2_star == 0 & y3_star == 1, 1,
        y2_star == 1 & y3_star == 0, 0,
        y2_star == 1 & y3_star == 1, 0
      ),
      d3_star_theta1_2 = fcase(
        y2_star == 0 & y3_star == 0, 0,
        y2_star == 0 & y3_star == 1, 0,
        y2_star == 1 & y3_star == 0, 1,
        y2_star == 1 & y3_star == 1, -1
      ),
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
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
  
  # Second derivatives (Hessian terms)
  df_hessian_temp1 <- df_probs_temp1 %>% 
    fmutate(
      d1_star_theta0_1_theta0_1 = fcase(
        y1_star == 0, 2*param$theta1_1/((param$theta1_1 + param$theta0_1)^3),
        y1_star == 1, -2*param$theta1_1/((param$theta1_1 + param$theta0_1)^3)
      ),
      d1_star_theta1_1_theta1_1 = fcase(
        y1_star == 0, -2*param$theta0_1/((param$theta1_1 + param$theta0_1)^3),
        y1_star == 1, 2*param$theta0_1/((param$theta1_1 + param$theta0_1)^3)
      ),
      d1_star_theta0_1_theta1_1 = fcase(
        y1_star == 0, (param$theta1_1 - param$theta0_1)/((param$theta1_1 + param$theta0_1)^3),
        y1_star == 1, (param$theta0_1 - param$theta1_1)/((param$theta1_1 + param$theta0_1)^3)
      ),
      d1_star_theta0_1_scaled = d1_star_theta0_1/p1_star,
      d1_star_theta1_1_scaled = d1_star_theta1_1/p1_star,
      d2_star_theta0_1_scaled = d2_star_theta0_1/p2_star,
      d2_star_theta1_1_scaled = d2_star_theta1_1/p2_star,
      d3_star_theta0_1_scaled = d3_star_theta0_1/p3_star,
      d3_star_theta1_1_scaled = d3_star_theta1_1/p3_star,
      d1_pi_scaled = d1_pi/p1,
      d2_pi_scaled = d2_pi/p2,
      d3_pi_scaled = d3_pi/p3,
      d1_star_theta0_1_theta0_1_scaled = d1_star_theta0_1_theta0_1/p1_star,
      d1_star_theta1_1_theta1_1_scaled = d1_star_theta1_1_theta1_1/p1_star,
      d1_star_theta0_1_theta1_1_scaled = d1_star_theta0_1_theta1_1/p1_star,
      joint_d2_theta0_1_theta0_1 = joint_p*(
        d1_star_theta0_1_theta0_1_scaled +
          2*d1_star_theta0_1_scaled*d2_star_theta0_1_scaled + 
          2*d1_star_theta0_1_scaled*d3_star_theta0_1_scaled + 
          2*d2_star_theta0_1_scaled*d3_star_theta0_1_scaled
      ),
      joint_d2_theta1_1_theta1_1 = joint_p*(
        d1_star_theta1_1_theta1_1_scaled +
          2*d1_star_theta1_1_scaled*d2_star_theta1_1_scaled + 
          2*d1_star_theta1_1_scaled*d3_star_theta1_1_scaled + 
          2*d2_star_theta1_1_scaled*d3_star_theta1_1_scaled
      ),
      joint_d2_pi_pi = joint_p*(
        2*d1_pi_scaled*d2_pi_scaled + 
          2*d1_pi_scaled*d3_pi_scaled + 
          2*d2_pi_scaled*d3_pi_scaled
      ),
      joint_d2_theta0_1_theta1_1 = joint_p*(
        d1_star_theta0_1_theta1_1_scaled +
          d1_star_theta0_1_scaled*(d2_star_theta1_1_scaled + d3_star_theta1_1_scaled) + 
          d2_star_theta0_1_scaled*(d1_star_theta1_1_scaled + d3_star_theta1_1_scaled) + 
          d3_star_theta0_1_scaled*(d1_star_theta1_1_scaled + d2_star_theta1_1_scaled)
      ),
      joint_d2_theta0_1_pi = joint_p*(
        d1_star_theta0_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta0_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta0_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      ),
      joint_d2_theta1_1_pi = joint_p*(
        d1_star_theta1_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta1_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta1_1_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      )
    )  %>% 
    fgroup_by(y1, y2, y3) %>%
    fsummarise(
      joint_d_theta0_1           = fsum(joint_d_theta0_1),
      joint_d_theta1_1           = fsum(joint_d_theta1_1),
      joint_d_pi_1               = fsum(joint_d_pi),
      joint_d2_theta0_1_theta0_1 = fsum(joint_d2_theta0_1_theta0_1),
      joint_d2_theta1_1_theta1_1 = fsum(joint_d2_theta1_1_theta1_1),
      joint_d2_pi_pi_1           = fsum(joint_d2_pi_pi),
      joint_d2_theta0_1_theta1_1 = fsum(joint_d2_theta0_1_theta1_1),
      joint_d2_theta0_1_pi       = fsum(joint_d2_theta0_1_pi),
      joint_d2_theta1_1_pi       = fsum(joint_d2_theta1_1_pi),
      joint_p1                   = fsum(joint_p)) |> 
    fungroup()
  
  df_hessian_temp2 <- df_probs_temp2 %>% 
    fmutate(
      d1_star_theta0_2_theta0_2 = fcase(
        y1_star == 0, 2*param$theta1_2/((param$theta1_2 + param$theta0_2)^3),
        y1_star == 1, -2*param$theta1_2/((param$theta1_2 + param$theta0_2)^3)
      ),
      d1_star_theta1_2_theta1_2 = fcase(
        y1_star == 0, -2*param$theta0_2/((param$theta1_2 + param$theta0_2)^3),
        y1_star == 1, 2*param$theta0_2/((param$theta1_2 + param$theta0_2)^3)
      ),
      d1_star_theta0_2_theta1_2 = fcase(
        y1_star == 0, (param$theta1_2 - param$theta0_2)/((param$theta1_2 + param$theta0_2)^3),
        y1_star == 1, (param$theta0_2 - param$theta1_2)/((param$theta1_2 + param$theta0_2)^3)
      ),
      d1_star_theta0_2_scaled = d1_star_theta0_2/p1_star,
      d1_star_theta1_2_scaled = d1_star_theta1_2/p1_star,
      d2_star_theta0_2_scaled = d2_star_theta0_2/p2_star,
      d2_star_theta1_2_scaled = d2_star_theta1_2/p2_star,
      d3_star_theta0_2_scaled = d3_star_theta0_2/p3_star,
      d3_star_theta1_2_scaled = d3_star_theta1_2/p3_star,
      d1_pi_scaled = d1_pi/p1,
      d2_pi_scaled = d2_pi/p2,
      d3_pi_scaled = d3_pi/p3,
      d1_star_theta0_2_theta0_2_scaled = d1_star_theta0_2_theta0_2/p1_star,
      d1_star_theta1_2_theta1_2_scaled = d1_star_theta1_2_theta1_2/p1_star,
      d1_star_theta0_2_theta1_2_scaled = d1_star_theta0_2_theta1_2/p1_star,
      joint_d2_theta0_2_theta0_2 = joint_p*(
        d1_star_theta0_2_theta0_2_scaled +
          2*d1_star_theta0_2_scaled*d2_star_theta0_2_scaled + 
          2*d1_star_theta0_2_scaled*d3_star_theta0_2_scaled + 
          2*d2_star_theta0_2_scaled*d3_star_theta0_2_scaled
      ),
      joint_d2_theta1_2_theta1_2 = joint_p*(
        d1_star_theta1_2_theta1_2_scaled +
          2*d1_star_theta1_2_scaled*d2_star_theta1_2_scaled + 
          2*d1_star_theta1_2_scaled*d3_star_theta1_2_scaled + 
          2*d2_star_theta1_2_scaled*d3_star_theta1_2_scaled
      ),
      joint_d2_pi_pi = joint_p*(
        2*d1_pi_scaled*d2_pi_scaled + 
          2*d1_pi_scaled*d3_pi_scaled + 
          2*d2_pi_scaled*d3_pi_scaled
      ),
      joint_d2_theta0_2_theta1_2 = joint_p*(
        d1_star_theta0_2_theta1_2_scaled +
          d1_star_theta0_2_scaled*(d2_star_theta1_2_scaled + d3_star_theta1_2_scaled) + 
          d2_star_theta0_2_scaled*(d1_star_theta1_2_scaled + d3_star_theta1_2_scaled) + 
          d3_star_theta0_2_scaled*(d1_star_theta1_2_scaled + d2_star_theta1_2_scaled)
      ),
      joint_d2_theta0_2_pi = joint_p*(
        d1_star_theta0_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta0_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta0_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      ),
      joint_d2_theta1_2_pi = joint_p*(
        d1_star_theta1_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d2_star_theta1_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled) + 
          d3_star_theta1_2_scaled*(d1_pi_scaled + d2_pi_scaled + d3_pi_scaled)
      )
    ) %>% 
    fgroup_by(y1, y2, y3) %>%
    fsummarise(
      joint_d_theta0_2           = fsum(joint_d_theta0_2),
      joint_d_theta1_2           = fsum(joint_d_theta1_2),
      joint_d_pi_2               = fsum(joint_d_pi),
      joint_d2_theta0_2_theta0_2 = fsum(joint_d2_theta0_2_theta0_2),
      joint_d2_theta1_2_theta1_2 = fsum(joint_d2_theta1_2_theta1_2),
      joint_d2_pi_pi_2           = fsum(joint_d2_pi_pi),
      joint_d2_theta0_2_theta1_2 = fsum(joint_d2_theta0_2_theta1_2),
      joint_d2_theta0_2_pi       = fsum(joint_d2_theta0_2_pi),
      joint_d2_theta1_2_pi       = fsum(joint_d2_theta1_2_pi),
      joint_p2                   = fsum(joint_p)) |> 
    fungroup()
  
  df_hessian_temp <- df_hessian_temp1 %>% 
    join(df_hessian_temp2, 
         on       = c('y1', 'y2', 'y3'), 
         validate = F, 
         verbose  = F) %>% 
    fmutate(
      joint_d_pi                 = param$p_1*joint_d_pi_1 + (1 - param$p_1)*joint_d_pi_2,
      joint_d2_theta0_1_theta0_1 = param$p_1*joint_d2_theta0_1_theta0_1,
      joint_d2_theta1_1_theta1_1 = param$p_1*joint_d2_theta1_1_theta1_1,
      joint_d2_theta0_2_theta0_2 = (1 - param$p_1)*joint_d2_theta0_2_theta0_2,
      joint_d2_theta1_2_theta1_2 = (1 - param$p_1)*joint_d2_theta1_2_theta1_2,
      joint_d2_pi_pi             = param$p_1*joint_d2_pi_pi_1 + (1 - param$p_1)*joint_d2_pi_pi_2,
      joint_d2_theta0_1_theta1_1 = param$p_1*joint_d2_theta0_1_theta1_1,
      joint_d2_theta0_1_pi       = param$p_1*joint_d2_theta0_1_pi,
      joint_d2_theta1_1_pi       = param$p_1*joint_d2_theta1_1_pi,
      joint_d2_theta0_2_theta1_2 = (1 - param$p_1)*joint_d2_theta0_2_theta1_2,
      joint_d2_theta0_2_pi       = (1 - param$p_1)*joint_d2_theta0_2_pi,
      joint_d2_theta1_2_pi       = (1 - param$p_1)*joint_d2_theta1_2_pi,
      joint_p                    = param$p_1*joint_p1 + (1 - param$p_1)*joint_p2,
      joint_d_theta0_1           = param$p_1*joint_d_theta0_1,
      joint_d_theta1_1           = param$p_1*joint_d_theta1_1,
      joint_d_theta0_2           = (1 - param$p_1)*joint_d_theta0_2,
      joint_d_theta1_2           = (1 - param$p_1)*joint_d_theta1_2)
  
  
  df_hi <- df_estimate %>%
    join(df_hessian_temp, 
         on       = c('y1', 'y2', 'y3'), 
         validate = F, 
         verbose  = F) %>%
    fmutate(
      joint_d2_theta0_1_theta0_1 = (weight/joint_p)*(param$theta0_1*(1 - param$theta0_1))*((joint_d2_theta0_1_theta0_1 - joint_d_theta0_1*joint_d_theta0_1/joint_p)*param$theta0_1*(1 - param$theta0_1) + (1 - 2*param$theta0_1)*joint_d_theta0_1),
      joint_d2_theta1_1_theta1_1 = (weight/joint_p)*(param$theta1_1*(1 - param$theta1_1))*((joint_d2_theta1_1_theta1_1 - joint_d_theta1_1*joint_d_theta1_1/joint_p)*param$theta1_1*(1 - param$theta1_1) + (1 - 2*param$theta1_1)*joint_d_theta1_1),
      joint_d2_theta0_2_theta0_2 = (weight/joint_p)*(param$theta0_2*(1 - param$theta0_2))*((joint_d2_theta0_2_theta0_2 - joint_d_theta0_2*joint_d_theta0_2/joint_p)*param$theta0_2*(1 - param$theta0_2) + (1 - 2*param$theta0_2)*joint_d_theta0_2),
      joint_d2_theta1_2_theta1_2 = (weight/joint_p)*(param$theta1_2*(1 - param$theta1_2))*((joint_d2_theta1_2_theta1_2 - joint_d_theta1_2*joint_d_theta1_2/joint_p)*param$theta1_2*(1 - param$theta1_2) + (1 - 2*param$theta1_2)*joint_d_theta1_2),
      joint_d2_pi_pi             =  (weight/joint_p)*(param$pi*(1 - param$pi))*((joint_d2_pi_pi - joint_d_pi*joint_d_pi/joint_p)*param$pi*(1 - param$pi) + (1 - 2*param$pi)*joint_d_pi),
      joint_d2_theta0_1_theta1_1 = (weight/joint_p)*(joint_d2_theta0_1_theta1_1 - joint_d_theta0_1*joint_d_theta1_1/joint_p)*param$theta0_1*(1 - param$theta0_1)*param$theta1_1*(1 - param$theta1_1),
      joint_d2_theta0_1_pi       = (weight/joint_p)*(joint_d2_theta0_1_pi - joint_d_theta0_1*joint_d_pi/joint_p)*param$theta0_1*(1 - param$theta0_1)*param$pi*(1 - param$pi),
      joint_d2_theta1_1_pi       = (weight/joint_p)*(joint_d2_theta1_1_pi - joint_d_theta1_1*joint_d_pi/joint_p)*param$theta1_1*(1 - param$theta1_1)*param$pi*(1 - param$pi),
      joint_d2_theta0_2_theta1_2 = (weight/joint_p)*(joint_d2_theta0_2_theta1_2 - joint_d_theta0_2*joint_d_theta1_2/joint_p)*param$theta0_2*(1 - param$theta0_2)*param$theta1_2*(1 - param$theta1_2),
      joint_d2_theta0_2_pi       = (weight/joint_p)*(joint_d2_theta0_2_pi - joint_d_theta0_2*joint_d_pi/joint_p)*param$theta0_2*(1 - param$theta0_2)*param$pi*(1 - param$pi),
      joint_d2_theta1_2_pi       = (weight/joint_p)*(joint_d2_theta1_2_pi - joint_d_theta1_2*joint_d_pi/joint_p)*param$theta1_2*(1 - param$theta1_2)*param$pi*(1 - param$pi),
      joint_d2_theta0_1_theta0_2 = (weight/joint_p)*(-joint_d_theta0_1*joint_d_theta0_2/joint_p)*param$theta0_1*(1 - param$theta0_1)*param$theta0_2*(1 - param$theta0_2),
      joint_d2_theta0_1_theta1_2 = (weight/joint_p)*(-joint_d_theta0_1*joint_d_theta1_2/joint_p)*param$theta0_1*(1 - param$theta0_1)*param$theta1_2*(1 - param$theta1_2),
      joint_d2_theta1_1_theta0_2 = (weight/joint_p)*(-joint_d_theta1_1*joint_d_theta0_2/joint_p)*param$theta1_1*(1 - param$theta1_1)*param$theta0_2*(1 - param$theta0_2),
      joint_d2_theta1_1_theta1_2 = (weight/joint_p)*(-joint_d_theta1_1*joint_d_theta1_2/joint_p)*param$theta1_1*(1 - param$theta1_1)*param$theta1_2*(1 - param$theta1_2),
      joint_d2_theta0_1_p_1      = (weight/joint_p)*(joint_d_theta0_1/param$p_1)*(1 - param$p_1*(joint_p1 - joint_p2)/joint_p)*param$theta0_1*(1 - param$theta0_1)*param$p_1*(1 - param$p_1),
      joint_d2_theta1_1_p_1      = (weight/joint_p)*(joint_d_theta1_1/param$p_1)*(1 - param$p_1*(joint_p1 - joint_p2)/joint_p)*param$theta1_1*(1 - param$theta1_1)*param$p_1*(1 - param$p_1),
      joint_d2_theta0_2_p_1      = -(weight/joint_p)*(joint_d_theta0_2/(1 - param$p_1))*(1 + (1 - param$p_1)*(joint_p1 - joint_p2)/joint_p)*param$theta0_2*(1 - param$theta0_2)*param$p_1*(1 - param$p_1),
      joint_d2_theta1_2_p_1      = -(weight/joint_p)*(joint_d_theta1_2/(1 - param$p_1))*(1 + (1 - param$p_1)*(joint_p1 - joint_p2)/joint_p)*param$theta1_2*(1 - param$theta1_2)*param$p_1*(1 - param$p_1),
      joint_d2_pi_p_1            = (weight/joint_p)*(joint_d_pi_1 - joint_p1*joint_d_pi/joint_p - joint_d_pi_2 + joint_p2*joint_d_pi/joint_p)*param$pi*(1 - param$pi)*param$p_1*(1 - param$p_1),
      joint_d2_p_1_p_1           =  (weight/joint_p)*param$p_1*(1 - param$p_1)*(1 - 2*param$p_1)*(joint_p1 - joint_p2) - (weight/joint_p)*((param$p_1*(1 - param$p_1)*(joint_p1 - joint_p2))^2)/joint_p,
      joint_d_theta0_1           = param$theta0_1*(1 - param$theta0_1)*weight*joint_d_theta0_1/joint_p,
      joint_d_theta1_1           = param$theta1_1*(1 - param$theta1_1)*weight*joint_d_theta1_1/joint_p,
      joint_d_theta0_2           = param$theta0_2*(1 - param$theta0_2)*weight*joint_d_theta0_2/joint_p,
      joint_d_theta1_2           = param$theta1_2*(1 - param$theta1_2)*weight*joint_d_theta1_2/joint_p,
      joint_d_pi_1               = weight*joint_d_pi_1/joint_p,
      joint_d_pi_2               = weight*joint_d_pi_2/joint_p,
      joint_d_pi                 = weight*joint_d_pi/joint_p)
  
  # Convert to matrix
  hessian_matrix <- matrix(c(
    fsum(df_hi$joint_d2_theta0_1_theta0_1),
    fsum(df_hi$joint_d2_theta0_1_theta1_1),
    fsum(df_hi$joint_d2_theta0_1_p_1), 
    fsum(df_hi$joint_d2_theta0_1_theta0_2),
    fsum(df_hi$joint_d2_theta0_1_theta1_2),
    fsum(df_hi$joint_d2_theta0_1_pi),
    
    fsum(df_hi$joint_d2_theta0_1_theta1_1),  
    fsum(df_hi$joint_d2_theta1_1_theta1_1),
    fsum(df_hi$joint_d2_theta1_1_p_1), 
    fsum(df_hi$joint_d2_theta1_1_theta0_2),
    fsum(df_hi$joint_d2_theta1_1_theta1_2),
    fsum(df_hi$joint_d2_theta1_1_pi),  
    
    fsum(df_hi$joint_d2_theta0_1_p_1), 
    fsum(df_hi$joint_d2_theta1_1_p_1), 
    fsum(df_hi$joint_d2_p_1_p_1),
    fsum(df_hi$joint_d2_theta0_2_p_1), 
    fsum(df_hi$joint_d2_theta1_2_p_1),
    fsum(df_hi$joint_d2_pi_p_1),
    
    fsum(df_hi$joint_d2_theta0_1_theta0_2),
    fsum(df_hi$joint_d2_theta1_1_theta0_2),
    fsum(df_hi$joint_d2_theta0_2_p_1),
    fsum(df_hi$joint_d2_theta0_2_theta0_2),
    fsum(df_hi$joint_d2_theta0_2_theta1_2),
    fsum(df_hi$joint_d2_theta0_2_pi),
    
    fsum(df_hi$joint_d2_theta0_1_theta1_2),
    fsum(df_hi$joint_d2_theta1_1_theta1_2),
    fsum(df_hi$joint_d2_theta1_2_p_1),
    fsum(df_hi$joint_d2_theta0_2_theta1_2),
    fsum(df_hi$joint_d2_theta1_2_theta1_2),
    fsum(df_hi$joint_d2_theta1_2_pi),
    
    fsum(df_hi$joint_d2_theta0_1_pi),
    fsum(df_hi$joint_d2_theta1_1_pi),
    fsum(df_hi$joint_d2_pi_p_1),
    fsum(df_hi$joint_d2_theta0_2_pi),
    fsum(df_hi$joint_d2_theta1_2_pi),
    fsum(df_hi$joint_d2_pi_pi)), 
    nrow  = 6, 
    ncol  = 6, 
    byrow = TRUE)
  
  # if(pi0) hessian_matrix <- hessian_matrix[1:2, 1:2] # Remove terms involving `pi` if `pi0 = TRUE`
  # 
  return(hessian_matrix)
  
}
# 
# 
# calc_mle_3waves_ar1_fmm2 <- function(param_transformed) {
#   ll <- sum(calc_lli_3waves_ar1_fmm2(param_transformed))
#   return(ll)
# }
# 
# calc_mle_derivatives_3waves_ar1_fmm2 <- function(param_transformed) {
#   lg <- colSums(calc_lli_derivatives_3waves_ar1_fmm2(param_transformed))
#   return(lg)
# }
# 
# # FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH 3 FMM groups ====
# 
# calc_lli_3waves_ar1_fmm3 <- function(param_transformed, pi0 = FALSE) {
#   
#   param <- logit_inverse(param_transformed)
#   mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
#   mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
#   mu_3 <- param$theta0_3/(param$theta0_3 + param$theta1_3)
#   if(pi0) param$pi <- 0
#   
#   df_probs_temp1 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_1,
#         y1_star == 0 & y2_star == 1, param$theta0_1,
#         y1_star == 1 & y2_star == 0, param$theta1_1,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_1
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_1,
#         y2_star == 0 & y3_star == 1, param$theta0_1,
#         y2_star == 1 & y3_star == 0, param$theta1_1,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_1
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   df_probs_temp2 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_2, 1 - mu_2),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_2,
#         y1_star == 0 & y2_star == 1, param$theta0_2,
#         y1_star == 1 & y2_star == 0, param$theta1_2,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_2
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_2,
#         y2_star == 0 & y3_star == 1, param$theta0_2,
#         y2_star == 1 & y3_star == 0, param$theta1_2,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_2
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   df_probs_temp3 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_3, 1 - mu_3),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_3,
#         y1_star == 0 & y2_star == 1, param$theta0_3,
#         y1_star == 1 & y2_star == 0, param$theta1_3,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_3
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_3,
#         y2_star == 0 & y3_star == 1, param$theta0_3,
#         y2_star == 1 & y3_star == 0, param$theta1_3,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_3
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3
#     ) 
#   
#   
#   df_probs1 <- df_probs_temp1 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p1 = sum(joint_p), .groups = "drop")
#   
#   df_probs2 <- df_probs_temp2 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p2 = sum(joint_p), .groups = "drop")
#   
#   df_probs3 <- df_probs_temp3 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(joint_p3 = sum(joint_p), .groups = "drop")
#   
#   
#   df_probs <- df_probs1 %>% 
#     left_join(df_probs2, by = c('y1', 'y2', 'y3')) %>% 
#     left_join(df_probs3, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(joint_p = param$p_1*joint_p1 + param$p_2*joint_p2 + (1 - param$p_1 - param$p_2)*joint_p3)
#   
#   df_lli <- df_estimate %>% 
#     left_join(df_probs, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(lli = weight*log(joint_p)) %>%
#     pull(lli)
#   
# }
# 
# 
# calc_lli_derivatives_3waves_ar1_fmm3 <- function(param_transformed, pi0 = FALSE) {
#   
#   param <- logit_inverse(param_transformed)
#   mu_1 <- param$theta0_1/(param$theta0_1 + param$theta1_1)
#   mu_2 <- param$theta0_2/(param$theta0_2 + param$theta1_2)
#   mu_3 <- param$theta0_3/(param$theta0_3 + param$theta1_3)
#   if(pi0) param$pi <- 0
#   
#   df_probs_temp1 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_1,
#         y1_star == 0 & y2_star == 1, param$theta0_1,
#         y1_star == 1 & y2_star == 0, param$theta1_1,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_1
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_1,
#         y2_star == 0 & y3_star == 1, param$theta0_1,
#         y2_star == 1 & y3_star == 0, param$theta1_1,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_1
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta0_1 = fcase(
#         y1_star == 0, -param$theta1_1/((param$theta1_1 + param$theta0_1)^2),
#         y1_star == 1, param$theta1_1/((param$theta1_1 + param$theta0_1)^2)
#       ),
#       d1_star_theta1_1 = fcase(
#         y1_star == 0, param$theta0_1/((param$theta1_1 + param$theta0_1)^2),
#         y1_star == 1, -param$theta0_1/((param$theta1_1 + param$theta0_1)^2)
#       ),
#       d2_star_theta0_1 = fcase(
#         y1_star == 0 & y2_star == 0, -1,
#         y1_star == 0 & y2_star == 1, 1,
#         y1_star == 1 & y2_star == 0, 0,
#         y1_star == 1 & y2_star == 1, 0
#       ),
#       d2_star_theta1_1 = fcase(
#         y1_star == 0 & y2_star == 0, 0,
#         y1_star == 0 & y2_star == 1, 0,
#         y1_star == 1 & y2_star == 0, 1,
#         y1_star == 1 & y2_star == 1, -1
#       ),
#       d3_star_theta0_1 = fcase(
#         y2_star == 0 & y3_star == 0, -1,
#         y2_star == 0 & y3_star == 1, 1,
#         y2_star == 1 & y3_star == 0, 0,
#         y2_star == 1 & y3_star == 1, 0
#       ),
#       d3_star_theta1_1 = fcase(
#         y2_star == 0 & y3_star == 0, 0,
#         y2_star == 0 & y3_star == 1, 0,
#         y2_star == 1 & y3_star == 0, 1,
#         y2_star == 1 & y3_star == 1, -1
#       ),
#       d1_pi = fifelse(y1 == y1_star, -1, 1),
#       d2_pi = fifelse(y2 == y2_star, -1, 1),
#       d3_pi = fifelse(y3 == y3_star, -1, 1),
#       joint_d_theta0_1 = 
#         d1_star_theta0_1*joint_p/p1_star + 
#         d2_star_theta0_1*joint_p/p2_star + 
#         d3_star_theta0_1*joint_p/p3_star,
#       joint_d_theta1_1 = 
#         d1_star_theta1_1*joint_p/p1_star + 
#         d2_star_theta1_1*joint_p/p2_star + 
#         d3_star_theta1_1*joint_p/p3_star,
#       joint_d_pi = 
#         d1_pi*joint_p/p1 + 
#         d2_pi*joint_p/p2 + 
#         d3_pi*joint_p/p3
#     ) 
#   
#   df_probs_temp2 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_2, 1 - mu_2),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_2,
#         y1_star == 0 & y2_star == 1, param$theta0_2,
#         y1_star == 1 & y2_star == 0, param$theta1_2,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_2
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_2,
#         y2_star == 0 & y3_star == 1, param$theta0_2,
#         y2_star == 1 & y3_star == 0, param$theta1_2,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_2
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta0_2 = fcase(
#         y1_star == 0, -param$theta1_2/((param$theta1_2 + param$theta0_2)^2),
#         y1_star == 1, param$theta1_2/((param$theta1_2 + param$theta0_2)^2)
#       ),
#       d1_star_theta1_2 = fcase(
#         y1_star == 0, param$theta0_2/((param$theta1_2 + param$theta0_2)^2),
#         y1_star == 1, -param$theta0_2/((param$theta1_2 + param$theta0_2)^2)
#       ),
#       d2_star_theta0_2 = fcase(
#         y1_star == 0 & y2_star == 0, -1,
#         y1_star == 0 & y2_star == 1, 1,
#         y1_star == 1 & y2_star == 0, 0,
#         y1_star == 1 & y2_star == 1, 0
#       ),
#       d2_star_theta1_2 = fcase(
#         y1_star == 0 & y2_star == 0, 0,
#         y1_star == 0 & y2_star == 1, 0,
#         y1_star == 1 & y2_star == 0, 1,
#         y1_star == 1 & y2_star == 1, -1
#       ),
#       d3_star_theta0_2 = fcase(
#         y2_star == 0 & y3_star == 0, -1,
#         y2_star == 0 & y3_star == 1, 1,
#         y2_star == 1 & y3_star == 0, 0,
#         y2_star == 1 & y3_star == 1, 0
#       ),
#       d3_star_theta1_2 = fcase(
#         y2_star == 0 & y3_star == 0, 0,
#         y2_star == 0 & y3_star == 1, 0,
#         y2_star == 1 & y3_star == 0, 1,
#         y2_star == 1 & y3_star == 1, -1
#       ),
#       d1_pi = fifelse(y1 == y1_star, -1, 1),
#       d2_pi = fifelse(y2 == y2_star, -1, 1),
#       d3_pi = fifelse(y3 == y3_star, -1, 1),
#       joint_d_theta0_2 = 
#         d1_star_theta0_2*joint_p/p1_star + 
#         d2_star_theta0_2*joint_p/p2_star + 
#         d3_star_theta0_2*joint_p/p3_star,
#       joint_d_theta1_2 = 
#         d1_star_theta1_2*joint_p/p1_star + 
#         d2_star_theta1_2*joint_p/p2_star + 
#         d3_star_theta1_2*joint_p/p3_star,
#       joint_d_pi = 
#         d1_pi*joint_p/p1 + 
#         d2_pi*joint_p/p2 + 
#         d3_pi*joint_p/p3
#     )  
#   
#   df_probs_temp3 <- df_template %>% 
#     mutate(
#       p1_star = fifelse(y1_star == 1, mu_3, 1 - mu_3),
#       p2_star = fcase(
#         y1_star == 0 & y2_star == 0, 1 - param$theta0_3,
#         y1_star == 0 & y2_star == 1, param$theta0_3,
#         y1_star == 1 & y2_star == 0, param$theta1_3,
#         y1_star == 1 & y2_star == 1, 1 - param$theta1_3
#       ),
#       p3_star = fcase(
#         y2_star == 0 & y3_star == 0, 1 - param$theta0_3,
#         y2_star == 0 & y3_star == 1, param$theta0_3,
#         y2_star == 1 & y3_star == 0, param$theta1_3,
#         y2_star == 1 & y3_star == 1, 1 - param$theta1_3
#       ),
#       p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
#       p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
#       p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
#       joint_p = p1_star*p1*p2_star*p2*p3_star*p3,
#       d1_star_theta0_3 = fcase(
#         y1_star == 0, -param$theta1_3/((param$theta1_3 + param$theta0_3)^2),
#         y1_star == 1, param$theta1_3/((param$theta1_3 + param$theta0_3)^2)
#       ),
#       d1_star_theta1_3 = fcase(
#         y1_star == 0, param$theta0_3/((param$theta1_3 + param$theta0_3)^2),
#         y1_star == 1, -param$theta0_3/((param$theta1_3 + param$theta0_3)^2)
#       ),
#       d2_star_theta0_3 = fcase(
#         y1_star == 0 & y2_star == 0, -1,
#         y1_star == 0 & y2_star == 1, 1,
#         y1_star == 1 & y2_star == 0, 0,
#         y1_star == 1 & y2_star == 1, 0
#       ),
#       d2_star_theta1_3 = fcase(
#         y1_star == 0 & y2_star == 0, 0,
#         y1_star == 0 & y2_star == 1, 0,
#         y1_star == 1 & y2_star == 0, 1,
#         y1_star == 1 & y2_star == 1, -1
#       ),
#       d3_star_theta0_3 = fcase(
#         y2_star == 0 & y3_star == 0, -1,
#         y2_star == 0 & y3_star == 1, 1,
#         y2_star == 1 & y3_star == 0, 0,
#         y2_star == 1 & y3_star == 1, 0
#       ),
#       d3_star_theta1_3 = fcase(
#         y2_star == 0 & y3_star == 0, 0,
#         y2_star == 0 & y3_star == 1, 0,
#         y2_star == 1 & y3_star == 0, 1,
#         y2_star == 1 & y3_star == 1, -1
#       ),
#       d1_pi = fifelse(y1 == y1_star, -1, 1),
#       d2_pi = fifelse(y2 == y2_star, -1, 1),
#       d3_pi = fifelse(y3 == y3_star, -1, 1),
#       joint_d_theta0_3 = 
#         d1_star_theta0_3*joint_p/p1_star + 
#         d2_star_theta0_3*joint_p/p2_star + 
#         d3_star_theta0_3*joint_p/p3_star,
#       joint_d_theta1_3 = 
#         d1_star_theta1_3*joint_p/p1_star + 
#         d2_star_theta1_3*joint_p/p2_star + 
#         d3_star_theta1_3*joint_p/p3_star,
#       joint_d_pi = 
#         d1_pi*joint_p/p1 + 
#         d2_pi*joint_p/p2 + 
#         d3_pi*joint_p/p3
#     )  
#   
#   
#   df_grad1 <- df_probs_temp1 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta0_1 = sum(joint_d_theta0_1), 
#       joint_d_theta1_1 = sum(joint_d_theta1_1),
#       joint_d_pi1 = sum(joint_d_pi),
#       joint_p1 = sum(joint_p),
#       .groups = "drop") 
#   
#   df_grad2 <- df_probs_temp2 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta0_2 = sum(joint_d_theta0_2), 
#       joint_d_theta1_2 = sum(joint_d_theta1_2),
#       joint_d_pi2 = sum(joint_d_pi),
#       joint_p2 = sum(joint_p),
#       .groups = "drop") 
#   
#   df_grad3 <- df_probs_temp3 %>% 
#     group_by(y1, y2, y3) %>% 
#     summarise(
#       joint_d_theta0_3 = sum(joint_d_theta0_3), 
#       joint_d_theta1_3 = sum(joint_d_theta1_3),
#       joint_d_pi3 = sum(joint_d_pi),
#       joint_p3 = sum(joint_p),
#       .groups = "drop") 
#   
#   df_grad <- df_grad1 %>% 
#     left_join(df_grad2, by = c('y1', 'y2', 'y3')) %>% 
#     left_join(df_grad3, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(
#       joint_p = param$p_1*joint_p1 + param$p_2*joint_p2 + (1 - param$p_1 - param$p_2)*joint_p3,
#       joint_d_pi = param$p_1*joint_d_pi1 + param$p_2*joint_d_pi2 + (1 - param$p_1 - param$p_2)*joint_d_pi3
#     ) %>% 
#     mutate(
#       joint_d_theta0_1 = param$p_1*joint_d_theta0_1*param$theta0_1*(1 - param$theta0_1), 
#       joint_d_theta1_1 = param$p_1*joint_d_theta1_1*param$theta1_1*(1 - param$theta1_1),
#       joint_d_theta0_2 = param$p_2*joint_d_theta0_2*param$theta0_2*(1 - param$theta0_2), 
#       joint_d_theta1_2 = param$p_2*joint_d_theta1_2*param$theta1_2*(1 - param$theta1_2),
#       joint_d_theta0_3 = (1 - param$p_1 - param$p_2)*joint_d_theta0_3*param$theta0_3*(1 - param$theta0_3), 
#       joint_d_theta1_3 = (1 - param$p_1 - param$p_2)*joint_d_theta1_3*param$theta1_3*(1 - param$theta1_3),
#       
#       joint_d_p_1 = joint_p1*param$p_1*(1 - param$p_1) - joint_p3*param$p_1*(1 - param$p_1),
#       joint_d_p_2 = joint_p2*param$p_2*(1 - param$p_2) - joint_p3*param$p_2*(1 - param$p_2),
#       joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
#     )
#   
#   df_gi <- df_estimate %>% 
#     left_join(df_grad, by = c('y1', 'y2', 'y3')) %>% 
#     mutate(
#       lgi_theta0_1 = weight*joint_d_theta0_1/joint_p,
#       lgi_theta1_1 = weight*joint_d_theta1_1/joint_p,
#       lgi_theta0_2 = weight*joint_d_theta0_2/joint_p,
#       lgi_theta1_2 = weight*joint_d_theta1_2/joint_p,
#       lgi_theta0_3 = weight*joint_d_theta0_3/joint_p,
#       lgi_theta1_3 = weight*joint_d_theta1_3/joint_p,
#       lgi_p_1 = weight*joint_d_p_1/joint_p,
#       lgi_p_2 = weight*joint_d_p_2/joint_p,
#       lgi_pi = weight*joint_d_pi/joint_p
#     ) %>%
#     select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_p_2, lgi_theta0_3, lgi_theta1_3, lgi_pi)
#   
#   if(pi0) df_gi <- df_gi %>% select(lgi_theta0_1, lgi_theta1_1, lgi_p_1, lgi_theta0_2, lgi_theta1_2, lgi_p_2, lgi_theta0_3, lgi_theta1_3)
#   return(df_gi)
#   
# }
# 
# calc_mle_3waves_ar1_fmm3 <- function(param_transformed) {
#   ll <- sum(calc_lli_3waves_ar1_fmm3(param_transformed))
#   return(ll)
# }
# 
# calc_mle_derivatives_3waves_ar1_fmm3 <- function(param_transformed) {
#   lg <- colSums(calc_lli_derivatives_3waves_ar1_fmm3(param_transformed))
#   return(lg)
# }
