# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_0)

# DEFINE FUNCTIONS ====

# Load functions
source("scripts/define_estimation_functions_4waves_gmm_ar2.R")
source("scripts/define_estimation_functions_4waves_mle_ar2.R")
source("scripts/define_estimation_functions_3waves_mle_ar1.R")
source("scripts/define_estimation_functions_4waves_mle_ar1.R")


df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 


# INGEST DATA ====

dt_nigeria_panel_3waves <- readRDS("C:/Users/rulof/Gitlab/sa_empl_transitions_misclassify/data/raw/dt_nigeria_panel_3waves.rds") %>% 
  select(
    id = hhid,
    wave,
    weight = pweights_w3,
    # y = poor_ext
    y = poor_umic 
  ) %>% 
  mutate(y = if_else(y, 1L , 0L)) %>% 
  pivot_wider(id_cols = c("id", "weight"), values_from = "y", names_from = "wave", names_prefix = "y")

dt_nigeria_panel_4waves <- readRDS("C:/Users/rulof/Gitlab/sa_empl_transitions_misclassify/data/raw/dt_nigeria_panel_4waves.Rds") %>% 
  select(
    id = hhid,
    wave,
    weight = pweights_w4,
    # y = poor_ext
    y = poor_umic 
  ) %>% 
  mutate(y = if_else(y, 1L , 0L)) %>% 
  pivot_wider(id_cols = c("id", "weight"), values_from = "y", names_from = "wave", names_prefix = "y")

df_x <- dt_nigeria_panel_4waves %>% 
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


#>> Three waves, ME ----

df_estimate <- dt_nigeria_panel_3waves

param_init <- data.frame(theta_1 = 0.085, theta_0 = 0.085, pi = 0.01)
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

#>> Three waves, no ME ----

param_init <- data.frame(theta_1 = 0.085, theta_0 = 0.085)
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

#>> Four waves, ME ----

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 
df_estimate <- dt_nigeria_panel_4waves

param_init <- data.frame(theta_1 = 0.085, theta_0 = 0.085, pi = 0.01)
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

param_init <- data.frame(theta_1 = 0.085, theta_0 = 0.085)
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
model_gmm_4w_ar2_se <- sqrt(((exp(model_gmm_4w_ar2$coefficients)/((1 + exp(model_gmm_4w_ar2$coefficients))^2))^2)*diag(vcov(model_gmm_4w_ar2)))

print(logit_inverse(model_gmm_4w_ar2$coefficients))
print(model_gmm_4w_ar2_se)
print(logit_inverse(model_gmm_4w_ar2$coefficients)/model_gmm_4w_ar2_se)

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
  method = "BFGSR",
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
  method = "NR",
  reltol = 0,
  gradtol = 0
)
model_mle_4w_ar2_pi0_se <- sqrt(((exp(model_mle_4w_ar2_pi0$estimate)/((1 + exp(model_mle_4w_ar2_pi0$estimate))^2))^2)*diag(vcov(model_mle_4w_ar2_pi0)))

print(logit_inverse(model_mle_4w_ar2_pi0$estimate))
print(model_mle_4w_ar2_pi0_se)
print(logit_inverse(model_mle_4w_ar2_pi0$estimate)/model_mle_4w_ar2_pi0_se)


# MODEL INFERENCE ====

#> AR(1) model, no ME ====

param <- logit_inverse(model_mle_4w_ar1_pi0$estimate)
param$pi <- 0
mu <- param$theta_0/(param$theta_1 + param$theta_0)
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
  )  %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop")

df_probs_temp_ar1_pi0 %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar1_pi0 %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)


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
  )  %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop")

df_probs_temp_ar1 %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar1 %>% 
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
  )  %>% 
  group_by(y1_star, y2_star, y3_star, y4_star) %>%
  summarise(joint_p = sum(joint_p), .groups = "drop")

df_probs_temp_ar2_pi0 %>% 
  group_by(y1_star, y4_star) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)

df_probs_temp_ar2_pi0 %>% 
  mutate(ever_found = pmax(y2_star, y3_star, y4_star)) %>% 
  group_by(y1_star, ever_found) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  filter(y1_star == 0) %>% 
  mutate(tot_sum = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/tot_sum)


#> AR(1) model, ME ====

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
  summarise(weight = sum(weight)) %>% 
  ungroup() %>% 
  mutate(totsum = sum(weight)) %>% 
  mutate(joint_p = weight/totsum)

df_probs_temp_ar2 %>% 
  group_by(y1, y2, y3, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")



# CREATE LATEX TABLE ====

lm_model <- lm(data = df_estimate, y4 ~ y1 + y2)

#> AR(1) model, 3 waves, ME ====
coefficients <- unlist(logit_inverse(model_mle_3w_ar1$estimate))
std_errors_1 <- unlist(model_mle_3w_ar1_se)
obs <- 3*nrow(df_estimate)
model_1 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_1^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_1) <- "lm"
stargazer::stargazer(model_1, type = "text", coeff = coefficients, se = std_errors)


stargazer::stargazer(model_2,
          # coef = c(0,1),
          se = c(0,2),
          # omit = c(sequence),
          covariate.labels = c("a", "b"),
          dep.var.labels.include = FALSE,
          notes.append=FALSE, type = "text")


#> AR(1) model, 3 waves, no ME ====
coefficients <- unlist(logit_inverse(model_mle_3w_ar1_pi0$estimate))
std_errors_2 <- unlist(model_mle_3w_ar1_pi0_se)
obs <- 3*nrow(df_estimate)
model_2 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs),
  vcov = diag(std_errors_2^2)
)
class(model_2) <- "lm"
stargazer::stargazer(model_2, type = "text")

#> AR(1) model, 4 waves, ME ====
coefficients <- unlist(logit_inverse(model_mle_4w_ar1$estimate))
std_errors_3 <- unlist(model_mle_4w_ar1_se)

obs <- 4*nrow(df_estimate)
model_3 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_3^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_3) <- "lm"
stargazer::stargazer(model_3, type = "text")

#> AR(1) model, 4 waves, no ME ====
coefficients <- unlist(logit_inverse(model_mle_4w_ar1_pi0$estimate))
std_errors_4 <- unlist(model_mle_4w_ar1_pi0_se)
obs <- 4*nrow(df_estimate)
model_4 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_4^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_4) <- "lm"
stargazer::stargazer(model_4, type = "text")

stargazer::stargazer(model_2, model_1, model_4, model_3, 
                     type = "text", align = TRUE, digits = 3, digits.extra = 2, keep.stat = c("n", "ll"), model.numbers = T,
                     dep.var.labels.include = FALSE, dep.var.caption = "",
                     se = list(std_errors_2, std_errors_1, std_errors_4, std_errors_3),
                     add.lines=list(c("Waves", "3", "3", "4", "4"))
)

stargazer::stargazer(model_2, model_1, model_4, model_3, 
                     type = "latex", align = TRUE, digits = 3, digits.extra = 2, keep.stat = c("n", "ll"), model.numbers = T,
                     dep.var.labels.include = FALSE, dep.var.caption = "",
                     se = list(std_errors_2, std_errors_1, std_errors_4, std_errors_3),
                     add.lines=list(c("Waves", "3", "3", "4", "4")),
                     out = "output/tables/Nigeria/ar1_poor_umic.tex"
                     )

lm_model <- lm(data = df_estimate, y4 ~ y1 + y2 + y3 + y1:y2)

#> AR(2) model, MLE, ME ====
coefficients <- unlist(logit_inverse(model_mle_4w_ar2$estimate))
std_errors_5 <- unlist(model_mle_4w_ar2_se)
obs <- 4*nrow(df_estimate)
model_5 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + theta_10 + theta_01 + pi),  # Model terms object
  vcov = diag(std_errors_5^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_5) <- "lm"
stargazer::stargazer(model_5, type = "text")

#> AR(2) model, MLE, no ME ====
coefficients <- unlist(logit_inverse(model_mle_4w_ar2_pi0$estimate))
std_errors_6 <- unlist(model_mle_4w_ar2_pi0_se)
obs <- 4*nrow(df_estimate)
model_6 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_6^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_6) <- "lm"
stargazer::stargazer(model_6, type = "text")

#> AR(2) model, GMM, ME ====
coefficients <- unlist(logit_inverse(model_gmm_4w_ar2$coefficients))
std_errors_7 <- unlist(model_gmm_4w_ar2_se)
obs <- 4*nrow(df_estimate)
model_7 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_7^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_7) <- "lm"
stargazer::stargazer(model_7, type = "text")

#> AR(2) model, GMM, no ME ====
coefficients <- unlist(logit_inverse(model_gmm_4w_ar2_pi0$coefficients))
std_errors_8 <- unlist(model_gmm_4w_ar2_pi0_se)
obs <- 4*nrow(df_estimate)
model_8 <- list(
  coefficients = coefficients,  # Coefficient estimates
  fitted.values = rnorm(nrow(df_estimate)),  # Fitted values from the model
  df.residual = nrow(df_estimate) - length(coefficients),  # Residual degrees of freedom (N - p)
  call = as.call(list(quote(lm), y ~ theta_0 + theta_1 + pi)),  # Model call (use as.call() with a formula)
  terms = terms(Y ~ theta_0 + theta_1 + pi),  # Model terms object
  vcov = diag(std_errors_8^2),  # Variance-covariance matrix of coefficients
  effects = nrow(df_estimate),
  rank = length(coefficients),
  qr = lm_model$qr,
  residuals = rnorm(obs)
)
class(model_8) <- "lm"
stargazer::stargazer(model_8, type = "text")

stargazer::stargazer(model_6, model_5, model_8, model_7, 
                     type = "text", align = TRUE, digits = 4, keep.stat = c("n", "ll"), model.numbers = T,
                     dep.var.labels.include = FALSE, dep.var.caption = "",
                     se = list(std_errors_6, std_errors_5, std_errors_8, std_errors_7),
                     add.lines=list(c("Estimator", "ML", "ML", "GMM", "GMM"))
)

stargazer::stargazer(model_6, model_5, model_8, model_7, 
                     type = "latex", align = TRUE, digits = 4, keep.stat = c("n", "ll"), model.numbers = T,
                     dep.var.labels.include = FALSE, dep.var.caption = "",
                     se = list(std_errors_6, std_errors_5, std_errors_8, std_errors_7),
                     add.lines=list(c("Estimator", "ML", "ML", "GMM", "GMM")),
                     out = "output/tables/Nigeria/ar2_poor_umic.tex"
                     )

