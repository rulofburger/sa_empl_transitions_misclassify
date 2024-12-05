# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

# Load functions
source("scripts/define_estimation_functions_4waves_gmm_ar2.R")
source("scripts/define_estimation_functions_4waves_mle_ar2.R")
source("scripts/define_estimation_functions_3waves_mle_ar1.R")
source("scripts/define_estimation_functions_4waves_mle_ar1.R")


df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 


# INGEST DATA ====

source("scripts/ingest_data_4waves.R")

df_estimate <- df_qlfs %>%
  # filter(period1 == 33)
  # filter(period1 == 31) # 0.01020117
filter(period1 >= 30 & period1 <= 32) #0.01153142
# filter(period1 >= 30 & period1 <= 33) #0.008116868 


df_estimate <- df_estimate %>% 
  mutate(weight_total = sum(weight)) %>% 
  mutate(weight = dim(df_estimate)[1]*weight/weight_total) %>% 
  select(y1, y2, y3, y4, weight, period1, period2, period3, period4)

df_x <- df_estimate %>% 
  mutate(
    y0000 = if_else(y1 == 0 & y2 == 0 & y3 == 0 & y4 == 0, 1L, 0L),
    y0001 = if_else(y1 == 0 & y2 == 0 & y3 == 0 & y4 == 1, 1L, 0L),
    y0010 = if_else(y1 == 0 & y2 == 0 & y3 == 1 & y4 == 0, 1L, 0L),
    y0011 = if_else(y1 == 0 & y2 == 0 & y3 == 1 & y4 == 1, 1L, 0L),
    y0100 = if_else(y1 == 0 & y2 == 1 & y3 == 0 & y4 == 0, 1L, 0L),
    y0101 = if_else(y1 == 0 & y2 == 1 & y3 == 0 & y4 == 1, 1L, 0L),
    y0110 = if_else(y1 == 0 & y2 == 1 & y3 == 1 & y4 == 0, 1L, 0L),
    y0111 = if_else(y1 == 0 & y2 == 1 & y3 == 1 & y4 == 1, 1L, 0L),
    y1000 = if_else(y1 == 1 & y2 == 0 & y3 == 0 & y4 == 0, 1L, 0L),
    y1001 = if_else(y1 == 1 & y2 == 0 & y3 == 0 & y4 == 1, 1L, 0L),
    y1010 = if_else(y1 == 1 & y2 == 0 & y3 == 1 & y4 == 0, 1L, 0L),
    y1011 = if_else(y1 == 1 & y2 == 0 & y3 == 1 & y4 == 1, 1L, 0L),
    y1100 = if_else(y1 == 1 & y2 == 1 & y3 == 0 & y4 == 0, 1L, 0L),
    y1101 = if_else(y1 == 1 & y2 == 1 & y3 == 0 & y4 == 1, 1L, 0L),
    y1110 = if_else(y1 == 1 & y2 == 1 & y3 == 1 & y4 == 0, 1L, 0L),
    y1111 = if_else(y1 == 1 & y2 == 1 & y3 == 1 & y4 == 1, 1L, 0L)
  ) %>% 
  select(y0000, y0001, y0010, y0011, y0100, y0101, y0110, y0111, y1000, y1001, y1010, y1011, y1100, y1101, y1110, y1111, weight) %>% 
  as.matrix()

# ESTIMATE MODELS ====

#> AR(1) model ====


df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 


#>> Three waves, no ME ----

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar1_pi0,
  grad = calc_mle_derivatives_3waves_ar1_pi0,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_pi0_se <- sqrt(((exp(model_mle_3w_ar1_pi0$estimate)/((1 + exp(model_mle_3w_ar1_pi0$estimate))^2))^2)*diag(vcov(model_mle_3w_ar1_pi0)))

print(logit_inverse(model_mle_3w_ar1_pi0$estimate))
print(model_mle_3w_ar1_pi0_se)
print(logit_inverse(model_mle_3w_ar1_pi0$estimate)/model_mle_3w_ar1_pi0_se)

#>> Three waves, ME ----

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, pi = 0.01)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar1 <- maxLik::maxLik(
  calc_mle_3waves_ar1,
  grad = calc_mle_derivatives_3waves_ar1,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar1_se <- sqrt(((exp(model_mle_3w_ar1$estimate)/((1 + exp(model_mle_3w_ar1$estimate))^2))^2)*diag(vcov(model_mle_3w_ar1)))

print(logit_inverse(model_mle_3w_ar1$estimate))
print(model_mle_3w_ar1_se)
print(logit_inverse(model_mle_3w_ar1$estimate)/model_mle_3w_ar1_se)

#>> Three waves, AR2, no ME ----

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615)
param_init_transformed <- logit_transform(param_init)

model_mle_3w_ar2_pi0 <- maxLik::maxLik(
  calc_mle_3waves_ar2_pi0,
  grad = calc_mle_derivatives_3waves_ar2_pi0,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_3w_ar2_pi0_se <- sqrt(((exp(model_mle_3w_ar2_pi0$estimate)/((1 + exp(model_mle_3w_ar2_pi0$estimate))^2))^2)*diag(vcov(model_mle_3w_ar2_pi0)))

print(logit_inverse(model_mle_3w_ar2_pi0$estimate))
print(model_mle_3w_ar2_pi0_se)
print(logit_inverse(model_mle_3w_ar2_pi0$estimate)/model_mle_3w_ar2_pi0_se)

#>> Four waves, ME ----

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085, pi = 0.01)
param_init_transformed <- logit_transform(param_init)

model_mle_4w_ar1 <- maxLik::maxLik(
  calc_mle_4waves_ar1,
  grad = calc_mle_derivatives_4waves_ar1,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_4w_ar1_se <- sqrt(((exp(model_mle_4w_ar1$estimate)/((1 + exp(model_mle_4w_ar1$estimate))^2))^2)*diag(vcov(model_mle_4w_ar1)))

print(logit_inverse(model_mle_4w_ar1$estimate))
print(model_mle_4w_ar1_se)
print(logit_inverse(model_mle_4w_ar1$estimate)/model_mle_4w_ar1_se)

#>> Four waves, no ME ----

param_init <- data.frame(theta_0 = 0.085, theta_1 = 0.085)
param_init_transformed <- logit_transform(param_init)

model_mle_4w_ar1_pi0 <- maxLik::maxLik(
  calc_mle_4waves_ar1_pi0,
  grad = calc_mle_derivatives_4waves_ar1_pi0,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_4w_ar1_pi0_se <- sqrt(((exp(model_mle_4w_ar1_pi0$estimate)/((1 + exp(model_mle_4w_ar1_pi0$estimate))^2))^2)*diag(vcov(model_mle_4w_ar1_pi0)))

print(logit_inverse(model_mle_4w_ar1_pi0$estimate))
print(model_mle_4w_ar1_pi0_se)
print(logit_inverse(model_mle_4w_ar1_pi0$estimate)/model_mle_4w_ar1_pi0_se)

#> AR(2) model ====

#>> GMM, ME ----

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615, pi = 0.01351413)
param_init_transformed <- logit_transform(param_init)
model_gmm_4w_ar2 <- gmm::gmm(calc_gmm_moments_4waves, df_x, unlist(param_init_transformed), gradv = calc_gmm_derivatives_4waves, traceIter = TRUE, prewhite = 0)
gmm_4w_ar2_se <- sqrt(((exp(model_gmm_4w_ar2$coefficients)/((1 + exp(model_gmm_4w_ar2$coefficients))^2))^2)*diag(vcov(model_gmm_4w_ar2)))

print(logit_inverse(model_gmm_4w_ar2$coefficients))
print(gmm_4w_ar2_se)
print(logit_inverse(model_gmm_4w_ar2$coefficients)/gmm_4w_ar2_se)

#>> GMM, no ME ----

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615)
param_init_transformed <- logit_transform(param_init)
model_gmm_4w_ar2_pi0 <- gmm::gmm(calc_gmm_moments_4waves_pi0, df_x, unlist(param_init_transformed), gradv = calc_gmm_derivatives_4waves_pi0, traceIter = TRUE, prewhite = 0)
model_gmm_4w_ar2_pi0_se <- sqrt(((exp(model_gmm_4w_ar2_pi0$coefficients)/((1 + exp(model_gmm_4w_ar2_pi0$coefficients))^2))^2)*diag(vcov(model_gmm_4w_ar2_pi0)))

print(logit_inverse(model_gmm_4w_ar2_pi0$coefficients))
print(model_gmm_4w_ar2_pi0_se)
print(logit_inverse(model_gmm_4w_ar2_pi0$coefficients)/model_gmm_4w_ar2_pi0_se)

#>> MLE, ME ----

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615, pi = 0.01351413)
param_init_transformed <- logit_transform(param_init)
model_mle_4w_ar2 <- maxLik::maxLik(
  calc_mle_4waves,
  grad = calc_mle_derivatives_4waves,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_4w_ar2_se <- sqrt(((exp(model_mle_4w_ar2$estimate)/((1 + exp(model_mle_4w_ar2$estimate))^2))^2)*diag(vcov(model_mle_4w_ar2)))

print(logit_inverse(model_mle_4w_ar2$estimate))
print(model_mle_4w_ar2_se)
print(logit_inverse(model_mle_4w_ar2$estimate)/model_mle_4w_ar2_se)

#>> MLE, no ME ----

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615)
param_init_transformed <- logit_transform(param_init)
model_mle_4w_ar2_pi0 <- maxLik::maxLik(
  calc_mle_4waves_pi0,
  grad = calc_mle_derivatives_4waves_pi0,
  start = param_init_transformed,
  method = "NM",
  reltol = 0,
  gradtol = 0
)
model_mle_4w_ar2_pi0_se <- sqrt(((exp(model_mle_4w_ar2_pi0$estimate)/((1 + exp(model_mle_4w_ar2_pi0$estimate))^2))^2)*diag(vcov(model_mle_4w_ar2_pi0)))

print(logit_inverse(model_mle_4w_ar2_pi0$estimate))
print(model_mle_4w_ar2_pi0_se)
print(logit_inverse(model_mle_4w_ar2_pi0$estimate)/model_mle_4w_ar2_pi0_se)


# MODEL INFERENCE ====

#> 3 WAVES ====
df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 

#> AR(1) model, no ME ====

param <- logit_inverse(model_mle_3w_ar1_pi0$estimate)
param$pi <- 0
mu <- param$theta_0/(param$theta_0 + param$theta_1)
df_probs_temp_3w_ar1_pi0 <- df_template %>% 
  mutate(
    p1_star = if_else(y1_star == 1, mu, 1 - mu),
    p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_0, 1 - param$theta_0)),
    p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_0, 1 - param$theta_0)),
    p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
    p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
    p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
    joint_p = p1*p1_star*p2*p2_star*p3*p3_star
  )  

#> AR(1) model, ME ====

param <- logit_inverse(model_mle_3w_ar1$estimate)
mu <- param$theta_0/(param$theta_0 + param$theta_1)
df_probs_temp_3w_ar1 <- df_template %>% 
  mutate(
    p1_star = if_else(y1_star == 1, mu, 1 - mu),
    p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_0, 1 - param$theta_0)),
    p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_0, 1 - param$theta_0)),
    p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
    p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
    p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
    joint_p = p1*p1_star*p2*p2_star*p3*p3_star
  )  

#> AR(2) model, no ME ====

param <- logit_inverse(model_mle_3w_ar2_pi0$estimate)
param$pi <- 0
Theta = param$theta_0*(param$theta_10 - 1) + param$theta_1*(param$theta_01 - 1)
df_probs_temp_3w_ar2_pi0 <- df_template %>% 
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
  )  %>% 
  group_by(y1, y2, y3) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop")



#> FOUR WAVES ----

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 


#> AR(1) model, no ME ====

param <- logit_inverse(model_mle_4w_ar1_pi0$estimate)
param$pi <- 0
mu <- param$theta_0/(param$theta_0 + param$theta_1)
df_probs_temp_ar1_pi0 <- df_template %>% 
  mutate(
    p1_star = if_else(y1_star == 1, mu, 1 - mu),
    p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_0, 1 - param$theta_0)),
    p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_0, 1 - param$theta_0)),
    p4_star = if_else(y3_star == 1, if_else(y4_star == 1, 1 - param$theta_1, param$theta_1), if_else(y4_star == 1, param$theta_0, 1 - param$theta_0)),
    p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
    p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
    p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
    p4 = if_else(y4 == y4_star, 1 - param$pi, param$pi),
    joint_p = p1*p1_star*p2*p2_star*p3*p3_star*p4_star*p4
  )  

df_probs_temp_ar1_pi0 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar1_pi0 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar1_pi0 %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")

#> AR(1) model, ME ====

param <- logit_inverse(model_mle_4w_ar1$estimate)
mu <- param$theta_0/(param$theta_1 + param$theta_0)
df_probs_temp_ar1 <- df_template %>% 
  mutate(
    p1_star = if_else(y1_star == 1, mu, 1 - mu),
    p2_star = if_else(y1_star == 1, if_else(y2_star == 1, 1 - param$theta_1, param$theta_1), if_else(y2_star == 1, param$theta_0, 1 - param$theta_0)),
    p3_star = if_else(y2_star == 1, if_else(y3_star == 1, 1 - param$theta_1, param$theta_1), if_else(y3_star == 1, param$theta_0, 1 - param$theta_0)),
    p4_star = if_else(y3_star == 1, if_else(y4_star == 1, 1 - param$theta_1, param$theta_1), if_else(y4_star == 1, param$theta_0, 1 - param$theta_0)),
    p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
    p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
    p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
    p4 = if_else(y4 == y4_star, 1 - param$pi, param$pi),
    joint_p = p1*p1_star*p2*p2_star*p3*p3_star*p4_star*p4
  )  

df_probs_temp_ar1 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar1 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

#> AR(2) model, no ME ====

param <- logit_inverse(model_mle_4w_ar2_pi0$estimate)
param$pi <- 0
Theta = param$theta_0*(param$theta_10 - 1) + param$theta_1*(param$theta_01 - 1)
df_probs_temp_ar2_pi0 <- df_template %>% 
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


df_probs_temp_ar2_pi0 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar2_pi0 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)


#> AR(2) model, ME ====

param <- logit_inverse(model_mle_4w_ar2$estimate)
Theta = param$theta_0*(param$theta_10 - 1) + param$theta_1*(param$theta_01 - 1)
df_probs_temp_ar2 <- df_template %>% 
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



df_probs_temp_ar2 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar2 %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_estimate %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(weight), .groups = "drop") %>%
  mutate(share = joint_p/sum(joint_p)) %>% 
  group_by(y1) %>% 
  summarise(share = sum(share), .groups = "drop")

df_estimate %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(weight), .groups = "drop") %>%
  mutate(share = joint_p/sum(joint_p)) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1 == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_estimate %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(weight), .groups = "drop") %>% 
  mutate(ever_found = pmax(y2, y3, y4)) %>% 
  group_by(y1, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1 == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)


df_estimate %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(weight = sum(weight)) %>% 
  ungroup() %>% 
  mutate(totsum = sum(weight)) %>% 
  mutate(joint_p = weight/totsum)

df_probs_temp_ar2 %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")
print(logit_inverse(model_mle_4w_ar1_pi0$estimate)/model_mle_4w_ar1_pi0_se)

# GOODNESS OF FIT ====

df_estimate_probs_3w <- df_estimate %>% 
  group_by(y1, y2, y3) %>% 
  summarise(weight = sum(weight)) %>% 
  ungroup() %>% 
  mutate(totsum = sum(weight)) %>% 
  mutate(empirical = weight/totsum) %>%
  select(-c('weight', 'totsum')) %>% 
  left_join(
    df_probs_temp_3w_ar1_pi0 %>% 
      group_by(y1, y2, y3) %>% 
      summarise(AR1_noME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3')
  ) %>% 
  left_join(
    df_probs_temp_3w_ar1 %>% 
      group_by(y1, y2, y3) %>% 
      summarise(AR1_ME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3')
  ) %>% 
  left_join(
    df_probs_temp_3w_ar2_pi0 %>% 
      group_by(y1, y2, y3) %>% 
      summarise(AR2_noME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3')
  )

df_estimate_probs_3w

df_estimate_probs_3w %>% 
  mutate(
    AR1_noME = AR1_noME - empirical,
    AR1_ME = AR1_ME - empirical,
    AR2_noME = AR2_noME - empirical
  ) %>% 
  summarise(
    AR1_noME = sum(AR1_noME^2),
    AR1_ME = sum(AR1_ME^2),
    AR2_noME = sum(AR2_noME^2)
  )

df_estimate_probs_4w <- df_estimate %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(weight = sum(weight)) %>% 
  ungroup() %>% 
  mutate(totsum = sum(weight)) %>% 
  mutate(empirical = weight/totsum) %>%
  select(-c('weight', 'totsum')) %>% 
  left_join(
    df_probs_temp_ar1_pi0 %>% 
      group_by(y1, y2, y3, y4) %>% 
      summarise(AR1_noME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3', 'y4')
  ) %>% 
  left_join(
    df_probs_temp_ar1 %>% 
      group_by(y1, y2, y3, y4) %>%
      summarise(AR1_ME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3', 'y4')
  ) %>% 
  left_join(
    df_probs_temp_ar2_pi0 %>% 
      group_by(y1, y2, y3, y4) %>%
      summarise(AR2_noME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3', 'y4')
  ) %>% 
  left_join(
    df_probs_temp_ar2 %>% 
      group_by(y1, y2, y3, y4) %>%
      summarise(AR2_ME = sum(joint_p), .groups = "drop"),
    by = c('y1', 'y2', 'y3', 'y4')
  )


df_estimate_probs_4w

df_estimate_probs_4w %>% 
  mutate(
    empirical = round(empirical, 4),
    AR1_noME = round(AR1_noME, 4),
    AR1_ME = round(AR1_ME, 4),
    AR2_noME = round(AR2_noME, 4),
    AR2_ME = round(AR2_ME, 4)
    ) 

df_estimate_probs_4w %>% 
  mutate(
    AR1_noME = AR1_noME - empirical,
    AR1_ME = AR1_ME - empirical,
    AR2_noME = AR2_noME - empirical,
    AR2_ME = AR2_ME - empirical
    ) %>% 
  summarise(
    AR1_noME = sum(AR1_noME^2),
    AR1_ME = sum(AR1_ME^2),
    AR2_noME = sum(AR2_noME^2),
    AR2_ME = sum(AR2_ME^2)
  )
