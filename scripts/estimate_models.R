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

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615, pi = 0.01351413)
param_init_transformed <- logit_transform(param_init)
model_gmm_4w_ar2 <- gmm::gmm(calc_gmm_moments_4waves, df_x, unlist(param_init_transformed), gradv = calc_gmm_derivatives_4waves, traceIter = TRUE, prewhite = 0)
gmm_4w_ar2_se <- sqrt(((exp(model_gmm_4w_ar2$coefficients)/((1 + exp(model_gmm_4w_ar2$coefficients))^2))^2)*diag(vcov(model_gmm_4w_ar2)))

print(logit_inverse(model_gmm_4w_ar2$coefficients))
print(gmm_4w_ar2_se)
print(logit_inverse(model_gmm_4w_ar2$coefficients)/gmm_4w_ar2_se)


param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615)
param_init_transformed <- logit_transform(param_init)
model_gmm_4w_ar2_pi0 <- gmm::gmm(calc_gmm_moments_4waves_pi0, df_x, unlist(param_init_transformed), gradv = calc_gmm_derivatives_4waves_pi0, traceIter = TRUE, prewhite = 0)
model_gmm_4w_ar2_pi0_se <- sqrt(((exp(model_gmm_4w_ar2_pi0$coefficients)/((1 + exp(model_gmm_4w_ar2_pi0$coefficients))^2))^2)*diag(vcov(model_gmm_4w_ar2_pi0)))

print(logit_inverse(model_gmm_4w_ar2_pi0$coefficients))
print(model_gmm_4w_ar2_pi0_se)
print(logit_inverse(model_gmm_4w_ar2_pi0$coefficients)/model_gmm_4w_ar2_pi0_se)

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




df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 

param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.01)
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

param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085)
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

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 

param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.01)
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

param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085)
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