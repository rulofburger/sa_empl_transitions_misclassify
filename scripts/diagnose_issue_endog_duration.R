# METADATA =====================================================================
# DESCRIPTION: Investigate why point estimates for endogenous duration model are
# so weird
# CREATED: 2025-06-26 (rulofburger)

# INITIALISE ==================================================================

# Source environment/functions

# Load libraries
library(tidyverse)
library(haven)
library(data.table)
library(fastverse)

#> Load estimation functions defined in other scripts ----
source("scripts/define_estimation_functions_3waves_mle_ar1_duration_endog.R")

# DEFINE PARAMETERS ===========================================================

# Initial parameter values that seem sensible
param_initial <- data.frame(
  theta_0 = 0.02881745, 
  theta_1 = 0.02918875, 
  lambda_g = 6.439774, 
  lambda_h = 2.298665, 
  sigma_g = 0.2231302, 
  sigma_h = 0.2231302, 
  pi = 0.02990933
)
param_initial$mu <- param_initial$theta_0/(param_initial$theta_1 + param_initial$theta_0)

# MLE parameter values (with minimum constraints of 0.1 for sigmas and lambdas)
param_final <- data.frame(
  theta_0 = 0.06953751,
  theta_1 = 0.1148854,
  lambda_g = 6.276432,
  lambda_h = 2.458684,
  sigma_g = 0.01,
  sigma_h = 0.1882159,
  pi = 0.3313137
)
param_final$mu <- param_final$theta_0/(param_final$theta_1 + param_final$theta_0)


# DEFINE FUNCTIONS ============================================================

# Define variant of likelihood function which outputs additional moments

calc_lli_3waves_ar1_duration_endog_moments <- function(param, pi_fixed = NULL, only_pi = NULL) {
  
  
  df_probs_temp <- df_template_duration_endog %>% 
    mutate(
      p1_star = if_else(y1_star == 1, param$mu, 1 - param$mu),
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
      p_g1 = case_when(
        y1 == 0 ~ 1,
        y1 == 1 ~ gamlss.dist::dexGAUS(tenure1, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_g2 = case_when(
        y2 == 0 ~ 1,
        y2 == 1 & y2_star == 0 ~ gamlss.dist::dexGAUS(tenure2, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ dnorm(tenure2 - 0.125, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_g)),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ dnorm(tenure2 - tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)) / 
          (1 - pnorm(-tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure2 - 0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) /
          (1 - gamlss.dist::pexGAUS(-0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_g3 = case_when(
        y3 == 0 ~ 1,
        y3 == 1 & y3_star == 0 ~ gamlss.dist::dexGAUS(tenure3, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) /
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ dnorm(tenure3 - 0.125, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_g)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ dnorm(tenure3 - tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)) /
          (1 - pnorm(-tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ dnorm(tenure3 - 0.375, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.375, mean = 0, sd = param$sigma_g)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ dnorm(tenure3 - tenure1 - 0.5, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-tenure1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure3 - 0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)/
          (1 - gamlss.dist::pexGAUS(-0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_h1 = case_when(
        y1 == 1 ~ 1,
        y1 == 0 ~ gamlss.dist::dexGAUS(timegap1, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),
      p_h2 = case_when(
        y2 == 1 ~ 1,
        y2 == 0 & y2_star == 1 ~ gamlss.dist::dexGAUS(timegap2, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ dnorm(timegap2 - 0.125, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_h)),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ dnorm(timegap2 - timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap2 - 0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(-0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),
      p_h3 = case_when(
        y3 == 1 ~ 1,
        y3 == 0 & y3_star == 1 ~ gamlss.dist::dexGAUS(timegap3, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ dnorm(timegap3 - 0.125, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_h)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ dnorm(timegap3 - timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ dnorm(timegap3 - 0.375, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.375, mean = 0, sd = param$sigma_h)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ dnorm(timegap3 - timegap1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap3 - 0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(-0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),   
      p1_star = if_else(p1_star < minval, minval, p1_star),
      p2_star = if_else(p2_star < minval, minval, p2_star),
      p3_star = if_else(p3_star < minval, minval, p3_star),
      p1 = if_else(p1 < minval, minval, p1),
      p2 = if_else(p2 < minval, minval, p2),
      p3 = if_else(p3 < minval, minval, p3),
      p_g1 = if_else(p_g1 < minval, minval, p_g1),
      p_g2 = if_else(p_g2 < minval, minval, p_g2),
      p_g3 = if_else(p_g3 < minval, minval, p_g3),
      p_h1 = if_else(p_h1 < minval, minval, p_h1),
      p_h2 = if_else(p_h2 < minval, minval, p_h2),
      p_h3 = if_else(p_h3 < minval, minval, p_h3),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3*p_g1*p_g2*p_g3*p_h1*p_h2*p_h3,
      joint_trans = p1_star*p1*p2_star*p2*p3_star*p3,
      joint_g = p1_star*p2_star*p3_star*p_g1*p_g2*p_g3,
      joint_h = p1_star*p2_star*p3_star*p_h1*p_h2*p_h3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3) %>% 
    summarise(
      joint_p = sum(joint_p),
      joint_trans = sum(joint_trans),
      joint_g = sum(joint_g),
      joint_h = sum(joint_h),
      p1 = sum(p1_star*p2_star*p3_star*p1),
      p2 = sum(p1_star*p2_star*p3_star*p2),
      p3 = sum(p1_star*p2_star*p3_star*p3),
      p_g1 = sum(p1_star*p2_star*p3_star*p_g1),
      p_g2 = sum(p1_star*p2_star*p3_star*p_g2),
      p_g3 = sum(p1_star*p2_star*p3_star*p_g3),
      p_h1 = sum(p1_star*p2_star*p3_star*p_h1),
      p_h2 = sum(p1_star*p2_star*p3_star*p_h2),
      p_h3 = sum(p1_star*p2_star*p3_star*p_h3),
      joint_wave2_p_g = sum(p1_star*p1*p2_star*p2*p_g1*p_g2),
      joint_wave2_p_h = sum(p1_star*p1*p2_star*p2*p_h1*p_h2),
      .groups = "drop"
    ) 
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'timegap3', 'tenure1', 'tenure2', 'tenure3'))
}


# Define part of likelihood function that calculates probabilities analytically

calc_analytical_probs <- function(param) {
  
  df_template_duration_endog %>% 
    mutate(
      p1_star = if_else(y1_star == 1, param$mu, 1 - param$mu),
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
      p_g1 = case_when(
        y1 == 0 ~ 1,
        y1 == 1 ~ gamlss.dist::dexGAUS(tenure1, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_g2 = case_when(
        y2 == 0 ~ 1,
        y2 == 1 & y2_star == 0 ~ gamlss.dist::dexGAUS(tenure2, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ dnorm(tenure2 - 0.125, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_g)),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ dnorm(tenure2 - tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)) /
          (1 - pnorm(-tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure2 - 0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) /
          (1 - gamlss.dist::pexGAUS(-0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_g3 = case_when(
        y3 == 0 ~ 1,
        y3 == 1 & y3_star == 0 ~ gamlss.dist::dexGAUS(tenure3, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) /
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ dnorm(tenure3 - 0.125, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_g)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ dnorm(tenure3 - tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)) /
          (1 - pnorm(-tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ dnorm(tenure3 - 0.375, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.375, mean = 0, sd = param$sigma_g)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ dnorm(tenure3 - tenure1 - 0.5, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-tenure1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure3 - 0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)/
          (1 - gamlss.dist::pexGAUS(-0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_h1 = case_when(
        y1 == 1 ~ 1,
        y1 == 0 ~ gamlss.dist::dexGAUS(timegap1, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),
      p_h2 = case_when(
        y2 == 1 ~ 1,
        y2 == 0 & y2_star == 1 ~ gamlss.dist::dexGAUS(timegap2, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ dnorm(timegap2 - 0.125, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_h)),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ dnorm(timegap2 - timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap2 - 0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(-0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),
      p_h3 = case_when(
        y3 == 1 ~ 1,
        y3 == 0 & y3_star == 1 ~ gamlss.dist::dexGAUS(timegap3, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ dnorm(timegap3 - 0.125, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_h)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ dnorm(timegap3 - timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ dnorm(timegap3 - 0.375, mean = 0, sd = param$sigma_h) /
          (1 - pnorm(-0.375, mean = 0, sd = param$sigma_h)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ dnorm(timegap3 - timegap1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_h)^2)) /
          (1 - pnorm(-timegap1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_h)^2))),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap3 - 0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h) / 
          (1 - gamlss.dist::pexGAUS(-0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h))
      ),
      p1_star = if_else(p1_star < minval, minval, p1_star),
      p2_star = if_else(p2_star < minval, minval, p2_star),
      p3_star = if_else(p3_star < minval, minval, p3_star),
      p1 = if_else(p1 < minval, minval, p1),
      p2 = if_else(p2 < minval, minval, p2),
      p3 = if_else(p3 < minval, minval, p3),
      p_g1 = if_else(p_g1 < minval, minval, p_g1),
      p_g2 = if_else(p_g2 < minval, minval, p_g2),
      p_g3 = if_else(p_g3 < minval, minval, p_g3),
      p_h1 = if_else(p_h1 < minval, minval, p_h1),
      p_h2 = if_else(p_h2 < minval, minval, p_h2),
      p_h3 = if_else(p_h3 < minval, minval, p_h3),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3*p_g1*p_g2*p_g3*p_h1*p_h2*p_h3,
      joint_trans = p1_star*p1*p2_star*p2*p3_star*p3,
      joint_g = p1_star*p2_star*p3_star*p_g1*p_g2*p_g3,
      joint_h = p1_star*p2_star*p3_star*p_h1*p_h2*p_h3
    ) 
}

# Calculate wave 2 tenure with finer grid of tenure1 and tenure2 combinations

calc_tenure2_probs_grid <- function(param) {
  df_tenure_grid_wave2_initial <- expand.grid(
    tenure1 = sort(unique(c(seq(-10, 40, by = 1/12), 0.125, 0.25))),
    tenure2 = sort(unique(c(seq(-10, 40, by = 1/12), 0.125, 0.25))),
    y1 = c(0,1),
    y2 = 1,
    y1_star = c(0,1),
    y2_star = c(0,1)
  ) %>% 
    mutate(
      p1_star = if_else(y1_star == 1, param$mu, 1 - param$mu),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y1_star == 0 & y2_star == 1 ~ param$theta_0,
        y1_star == 1 & y2_star == 0 ~ param$theta_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p_g1 = case_when(
        y1 == 0 ~ 1,
        y1 == 1 ~ gamlss.dist::dexGAUS(tenure1, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p_g2 = case_when(
        y2 == 0 ~ 1,
        y2 == 1 & y2_star == 0 ~ gamlss.dist::dexGAUS(tenure2, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) / 
          (1 - gamlss.dist::pexGAUS(0, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ dnorm(tenure2 - 0.125, mean = 0, sd = param$sigma_g) /
          (1 - pnorm(-0.125, mean = 0, sd = param$sigma_g)),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ dnorm(tenure2 - tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)) /
          (1 - pnorm(-tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2))),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure2 - 0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g) /
          (1 - gamlss.dist::pexGAUS(-0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g))
      ),
      p1_star = if_else(p1_star < minval, minval, p1_star),
      p2_star = if_else(p2_star < minval, minval, p2_star),
      p1 = if_else(p1 < minval, minval, p1),
      p2 = if_else(p2 < minval, minval, p2),
      p_g1 = if_else(p_g1 < minval, minval, p_g1),
      p_g2 = if_else(p_g2 < minval, minval, p_g2),
      joint_p = p1_star*p1*p2_star*p2*p_g1*p_g2,
      joint_trans = p1_star*p1*p2_star*p2,
      joint_g = p_g1*p_g2
    )
}



# INGEST DATA ====

# Run script that loads 3 wave SA data as df_qlfs
source("scripts/ingest_data_3waves_SA.R")

# Limit survey rounds and calculate weights to be consistent within panel and to sum to 1
df_qlfs <- df_qlfs %>% 
  filter(period1 >= 30 & period1 <= 32)

df_qlfs_tenure_timegap <- df_qlfs %>% 
  filter(!is.na(timegap1) & !is.na(timegap2) & !is.na(timegap3)) %>% 
  filter(!is.na(tenure1) & !is.na(tenure2) & !is.na(tenure3))

df_qlfs_tenure_timegap <- df_qlfs_tenure_timegap %>% 
  mutate(weight_total = sum(weight))  %>% 
  mutate(weight = dim(df_qlfs_tenure_timegap)[1]*weight/weight_total) %>% 
  select(-weight_total)

df_estimate <- df_qlfs_tenure_timegap

df_covariate_combos <- df_estimate %>%
  select(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3) %>%
  unique() %>%
  na.omit %>%
  data.table()

df_template_employment <- data.table::CJ(
  y1_star = c(0, 1),
  y2_star = c(0, 1),
  y3_star = c(0, 1)
)

# Add a temporary key column without overwriting any existing data
df_covariate_combos[, k := 1]
df_template_employment[, k := 1]

# Perform the join while ensuring no columns are omitted
df_template_duration_endog <- df_covariate_combos[df_template_employment, on = "k", allow.cartesian = TRUE]

# Remove the temporary key column `k` from the result
df_template_duration_endog[, k := NULL]


# CREATE DATA OBJECTS =========================================================

# Calculate empirical likelihoods using likelihood function and survey data

df_empirical_individual_moments_final <- calc_lli_3waves_ar1_duration_endog_moments(param_final)

df_empirical_individual_moments_initial <- calc_lli_3waves_ar1_duration_endog_moments(param_initial)

df_empirical_moments_final <- df_empirical_individual_moments_final %>% 
  group_by(y1, y2, y3) %>% 
  summarise(
    lli = sum(weight*log(joint_p)),
    joint_trans = sum(weight*log(joint_trans)),
    joint_g = sum(weight*log(joint_g)),
    joint_h = sum(weight*log(joint_h)),
    p1 = sum(weight*log(p1)),
    p2 = sum(weight*log(p2)),
    p3 = sum(weight*log(p3)),
    p_g1 = sum(weight*log(p_g1)),
    p_g2 = sum(weight*log(p_g2)),
    p_g3 = sum(weight*log(p_g3)),
    p_h1 = sum(weight*log(p_h1)),
    p_h2 = sum(weight*log(p_h2)),
    p_h3 = sum(weight*log(p_h3))
  ) %>% 
  ungroup()

df_empirical_moments_initial <- df_empirical_individual_moments_initial %>% 
  group_by(y1, y2, y3) %>% 
  summarise(
    lli = sum(weight*log(joint_p)),
    joint_trans = sum(weight*log(joint_trans)),
    joint_g = sum(weight*log(joint_g)),
    joint_h = sum(weight*log(joint_h)),
    p1 = sum(weight*log(p1)),
    p2 = sum(weight*log(p2)),
    p3 = sum(weight*log(p3)),
    p_g1 = sum(weight*log(p_g1)),
    p_g2 = sum(weight*log(p_g2)),
    p_g3 = sum(weight*log(p_g3)),
    p_h1 = sum(weight*log(p_h1)),
    p_h2 = sum(weight*log(p_h2)),
    p_h3 = sum(weight*log(p_h3))
  ) %>% 
  ungroup()

df_empirical_individual_moments_pooled <- df_empirical_individual_moments_initial %>% 
  mutate(model = "Initial") %>% 
  bind_rows(
    df_empirical_individual_moments_final %>% 
      mutate(model = "Final")
  )


# Calculate analytical likelihoods using only likelihood functions

df_analytical_individual_moments_final <- calc_analytical_probs(param_final)

df_analytical_individual_moments_initial <- calc_analytical_probs(param_initial)

df_analytical_individual_moments_pooled <- df_analytical_individual_moments_initial %>% 
  mutate(model = "Initial") %>% 
  bind_rows(
    df_analytical_individual_moments_final %>% 
      mutate(model = "Final")
  )

df_analytical_moments_pooled <- df_analytical_individual_moments_pooled %>% 
  group_by(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3, model) %>% 
  summarise(
    joint_p = sum(joint_p),
    p1 = sum(p1_star*p2_star*p3_star*p1),
    p2 = sum(p1_star*p2_star*p3_star*p2),
    p3 = sum(p1_star*p2_star*p3_star*p3),
    p_g1 = sum(p1_star*p2_star*p3_star*p_g1),
    p_g2 = sum(p1_star*p2_star*p3_star*p_g2),
    p_g3 = sum(p1_star*p2_star*p3_star*p_g3),
    p_h1 = sum(p1_star*p2_star*p3_star*p_h1),
    p_h2 = sum(p1_star*p2_star*p3_star*p_h2),
    p_h3 = sum(p1_star*p2_star*p3_star*p_h3),
    .groups = "drop") 

# Calculate analytical likelihoods for wave 2 tenure with finer grid 

df_tenure2_probs_grid_initial <- calc_tenure2_probs_grid(param_initial) %>% 
  filter(tenure1 >= 0, tenure2 >= 0)

df_tenure2_probs_grid_final <- calc_tenure2_probs_grid(param_final) %>% 
  filter(tenure1 >= 0, tenure2 >= 0)

df_tenure2_probs_grid_pooled <- df_tenure2_probs_grid_initial %>% 
  mutate(model = "Initial") %>% 
  bind_rows(
    df_tenure2_probs_grid_final %>% 
      mutate(model = "Final")
  ) %>% 
  mutate(model = factor(model, levels = c("Initial", "Final")))


# ANALYSIS ====================================================================

#> Compare likelihood components across different outcome combinations ----

df_empirical_moments_final %>% 
  select(c('y1', 'y2', 'y3')) %>% 
  bind_cols(
    (df_empirical_moments_final %>% select(-c('y1', 'y2', 'y3')))/(df_empirical_moments_initial %>% select(-c('y1', 'y2', 'y3')))
  )


#> Density of g1 ----

# The distribution is the same, regardless of whether y1_star is 0 or 1

ggplot() +
  geom_line(
    data = df_tenure2_probs_grid_pooled %>% 
      filter(y1 == 1) %>% 
      group_by(tenure1, model) %>% 
      summarise(
        p_g1 = sum(p_g1*joint_trans)/sum(joint_trans),
        .groups = "drop"
      ),
    aes(x = tenure1, y = p_g1, color = model),
    linewidth = 1
  ) +
  geom_density(
    data = df_qlfs_tenure_timegap %>% filter(y1 == 1),
    aes(x = tenure1,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()


#> Density of h1 ----

ggplot() +
  # First: add your pre-computed density curve
  geom_line(
    data = df_analytical_moments_pooled %>% 
      filter(y1 == 0),
    aes(x = timegap1, y = p_h1, color = model),
    linewidth = 1
  ) +
  # Second: add the kernel density estimate from raw data
  geom_density(
    data = df_qlfs_tenure_timegap %>% filter(y1 == 0),
    aes(x = timegap1,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Unemployment duration",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()

#> Density of g2 ----

#>> Hypothetical distributions ---- 

# Case 1: y2_star == 0
df_tenure2_probs_grid_pooled %>% 
  filter(y2 == 1, y2_star == 0) %>% 
  group_by(model, tenure2) %>% 
  summarise(
    p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = tenure2, y = p_g2, color = model)) +
  geom_line() +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Comparison of Modelled Wave 2 Tenure Densities: y2_star == 0"
  ) +
  theme_minimal()

# Case 3: y2_star == 1 & y1_star == 1 & y1 == 0
df_tenure2_probs_grid_pooled %>% 
  filter(y2 == 1, y2_star == 1, y1_star == 1, y1 == 0) %>% 
  group_by(model, tenure2) %>% 
  summarise(
    p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = tenure2, y = p_g2, color = model)) +
  geom_line() +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Comparison of Modelled Wave 2 Tenure Densities: y2_star == 1, y1_star == 1, y1 == 0"
  ) +
  theme_minimal()

# Case 2: y2_star == 1 & y1_star == 0
df_tenure2_probs_grid_pooled %>% 
  filter(y2 == 1, y2_star == 1, y1_star == 0) %>% 
  # filter(tenure2 <= 5) %>% 
  group_by(model, tenure2) %>% 
  summarise(
    p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = tenure2, y = p_g2, color = model)) +
  geom_line() +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Comparison of Modelled Wave 2 Tenure Densities: y2_star == 1, y1_star == 0"
  ) +
  theme_minimal()

# Case4: y2_star == 1 & y1_star == 1 & y1 == 1
df_tenure2_probs_grid_pooled %>% 
  filter(y2 == 1, y2_star == 1 & y1_star == 1 & y1 == 1) %>% 
  group_by(model, tenure2) %>% 
  summarise(
    p_g2 = {
      t1_vals <- sort(tenure1)
      densities <- joint_g[order(tenure1)]
      sum(diff(t1_vals) * (head(densities, -1) + tail(densities, -1)) / 2)
    },
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = tenure2, y = p_g2, color = model)) +
  geom_line() +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Comparison of Modelled Wave 2 Tenure Densities: y2_star == 1 & y1_star == 1 & y1 == 1"
  ) +
  theme_minimal()

df_tenure2_probs_grid_pooled %>% 
  filter(y2 == 1, y2_star == 1 & y1_star == 1 & y1 == 1) %>% 
  filter(tenure1 != 0.125 & tenure2 != 0.125) %>% 
  ggplot(aes(x = tenure1, y = tenure2, fill = log(p_g2))) +
  geom_raster() +
  facet_grid(cols = vars(model)) +
  scale_fill_viridis_c(option = "plasma") +  # Or use other palettes like "magma", "viridis", etc.
  labs(
    title = "Joint (Log) Density of Waves 1 and 2 Tenure Heatmap",
    x = "Tenure 1",
    y = "Tenure 2",
    fill = "Log Density"
  ) +
  theme_minimal()


#>> Density of g2 if y1 == 0 ----

df_tenure2_probs_grid_pooled %>% 
  select(y1, y2, y1_star, y2_star, joint_trans, model) %>% 
  filter(y2 == 1, y1 == 0) %>% 
  group_by(y1_star, y2_star, model) %>% 
  summarise(joint_p = sum(joint_trans)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  mutate(sum_p = sum(joint_p)) %>% 
  mutate(conditional_p = joint_p/sum_p) %>% 
  mutate(
    case = case_when(
      y1_star == 0 & y2_star == 0 ~ "00",
      y1_star == 0 & y2_star == 1 ~ "01",
      y1_star == 1 & y2_star == 0 ~ "10",
      y1_star == 1 & y2_star == 1 ~ "11"
    )) %>%
  ggplot(aes(x = case, y = conditional_p, fill = case)) + 
  geom_col() + 
  facet_grid(cols = vars(model))  +
  labs(
    x = "Combination: y1_star and y2_star",
    y = "Conditional probability",
    fill = "Values of y1_star and y2_star",
    title = "Comparison of conditional probabilities"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

df_tenure2_probs_grid_pooled %>% 
  select(y1, y2, y1_star, y2_star, joint_trans, model) %>% 
  filter(y2 == 1, y1 == 0) %>% 
  group_by(y1_star, y2_star, model) %>% 
  summarise(joint_p = sum(joint_trans)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  mutate(sum_p = sum(joint_p)) %>% 
  mutate(conditional_p = joint_p/sum_p) %>% 
  mutate(
    case = case_when(
      y1_star == 0 & y2_star == 0 ~ "Case 1",
      y1_star == 0 & y2_star == 1 ~ "Case 2",
      y1_star == 1 & y2_star == 0 ~ "Case 1",
      y1_star == 1 & y2_star == 1 ~ "Case 3"
    )) %>%
  ggplot(aes(x = case, y = conditional_p, fill = case)) + 
  geom_col() + 
  facet_grid(cols = vars(model))  +
  labs(
    x = "Combination: y1_star and y2_star",
    y = "Conditional probability",
    fill = "Values of y1_star and y2_star",
    title = "Comparison of (grouped) conditional probabilities"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot() +
  geom_line(
    data = df_tenure2_probs_grid_pooled %>% 
      filter(y2 == 1, y1 == 0) %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
        .groups = "drop"
      ),
    aes(x = tenure2, y = p_g2, color = model),
    linewidth = 1
  ) +
  geom_density(
    data = df_qlfs_tenure_timegap %>% filter(y2 == 1, y1 == 0),
    aes(x = tenure2,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()

ggplot() +
  geom_line(
    data = df_tenure2_probs_grid_pooled %>% 
      filter(y2 == 1, y1 == 0) %>% 
      filter(tenure2 <= 10) %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
        .groups = "drop"
      ),
    aes(x = tenure2, y = p_g2, color = model),
    linewidth = 1
  ) +
  geom_density(
    data = df_qlfs_tenure_timegap %>% 
      filter(y2 == 1, y1 == 0) %>% 
      filter(tenure2 <= 10),
    aes(x = tenure2,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()

df_empirical_individual_moments_pooled %>% 
  select(model, y1, y2, tenure1, tenure2, joint_wave2_p_g) %>% 
  unique() %>% 
  pivot_wider(id_cols = c('y1', 'y2', 'tenure1', 'tenure2'), values_from = joint_wave2_p_g, names_from = model) %>% 
  mutate(diff = Final - Initial) %>% 
  filter(y1 == 0, y2 == 1) %>% 
  ggplot(aes(x = tenure2, y = diff)) +
  geom_line() +
  labs(
    title = "Final model fit improvement",
    x = "Tenure 2",
  ) +
  theme_minimal()

#>> Density of g2 if y1 == 1 ----

df_tenure2_probs_grid_pooled %>% 
  select(y1, y2, y1_star, y2_star, joint_trans, model) %>% 
  filter(y2 == 1, y1 == 1) %>% 
  group_by(y1_star, y2_star, model) %>% 
  summarise(joint_p = sum(joint_trans)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  mutate(sum_p = sum(joint_p)) %>% 
  mutate(conditional_p = joint_p/sum_p) %>% 
  mutate(
    case = case_when(
      y1_star == 0 & y2_star == 0 ~ "00",
      y1_star == 0 & y2_star == 1 ~ "01",
      y1_star == 1 & y2_star == 0 ~ "10",
      y1_star == 1 & y2_star == 1 ~ "11"
    )) %>%
  ggplot(aes(x = case, y = conditional_p, fill = case)) + 
  geom_col() + 
  facet_grid(cols = vars(model))  +
  labs(
    x = "Combination: y1_star and y2_star",
    y = "Conditional probability",
    fill = "Values of y1_star and y2_star",
    title = "Comparison of conditional probabilities"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


df_tenure2_probs_grid_pooled %>% 
  select(y1, y2, y1_star, y2_star, joint_trans, model) %>% 
  filter(y2 == 1, y1 == 1) %>% 
  group_by(y1_star, y2_star, model) %>% 
  summarise(joint_p = sum(joint_trans)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  mutate(sum_p = sum(joint_p)) %>% 
  mutate(conditional_p = joint_p/sum_p) %>% 
  mutate(
    case = case_when(
      y1_star == 0 & y2_star == 0 ~ "Case 1",
      y1_star == 0 & y2_star == 1 ~ "Case 2",
      y1_star == 1 & y2_star == 0 ~ "Case 1",
      y1_star == 1 & y2_star == 1 ~ "Case 4"
    )) %>%
  ggplot(aes(x = case, y = conditional_p, fill = case)) + 
  geom_col() + 
  facet_grid(cols = vars(model))  +
  labs(
    x = "Combination: y1_star and y2_star",
    y = "Conditional probability",
    fill = "Values of y1_star and y2_star",
    title = "Comparison of conditional probabilities"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot() +
  geom_line(
    data = df_tenure2_probs_grid_pooled %>% 
      filter(y2 == 1, y1 == 1 & !(y2_star == 1 & y1_star == 1)) %>% 
      group_by(model, tenure2, y1, y2, y1_star, y2_star) %>% 
      summarise(
        p_g2 = mean(p_g2),
        joint_trans = mean(joint_trans)
      ) %>% 
      ungroup() %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
        joint_trans = mean(joint_trans)
      ) %>% 
      bind_rows(
        df_tenure2_probs_grid_pooled %>% 
          filter(y2 == 1, y2_star == 1 & y1_star == 1 & y1 == 1) %>% 
          group_by(model, tenure2) %>% 
          summarise(
            p_g2 = {
              t1_vals <- sort(tenure1)
              densities <- joint_g[order(tenure1)]
              sum(diff(t1_vals) * (head(densities, -1) + tail(densities, -1)) / 2)
            },
            joint_trans = mean(joint_trans),
            .groups = "drop"
          )
      ) %>% 
      ungroup() %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
      ),
    aes(x = tenure2, y = p_g2, color = model),
    linewidth = 1
  ) + 
  geom_density(
    data = df_qlfs_tenure_timegap %>% filter(y2 == 1, y1 == 1),
    aes(x = tenure2,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()

ggplot() +
  geom_line(
    data = df_tenure2_probs_grid_pooled %>% 
      filter(y2 == 1, y1 == 1 & !(y2_star == 1 & y1_star == 1)) %>% 
      filter(tenure2 <= 10) %>% 
      group_by(model, tenure2, y1, y2, y1_star, y2_star) %>% 
      summarise(
        p_g2 = mean(p_g2),
        joint_trans = mean(joint_trans)
      ) %>% 
      ungroup() %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
        joint_trans = mean(joint_trans)
      ) %>% 
      bind_rows(
        df_tenure2_probs_grid_pooled %>% 
          filter(y2 == 1, y2_star == 1 & y1_star == 1 & y1 == 1) %>% 
          filter(tenure2 <= 10) %>% 
          group_by(model, tenure2) %>% 
          summarise(
            p_g2 = {
              t1_vals <- sort(tenure1)
              densities <- joint_g[order(tenure1)]
              sum(diff(t1_vals) * (head(densities, -1) + tail(densities, -1)) / 2)
            },
            joint_trans = mean(joint_trans),
            .groups = "drop"
          )
      ) %>% 
      ungroup() %>% 
      group_by(model, tenure2) %>% 
      summarise(
        p_g2 = sum(p_g2*joint_trans)/sum(joint_trans),
      ),
    aes(x = tenure2, y = p_g2, color = model),
    linewidth = 1
  ) + 
  geom_density(
    data = df_qlfs_tenure_timegap %>% 
      filter(y2 == 1, y1 == 1) %>% 
      filter(tenure2 <= 10),
    aes(x = tenure2,  weight = weight),
    color = "black",
    linetype = "dashed",
    adjust = 0.25
  ) +
  labs(
    x = "Tenure",
    y = "Density",
    color = "Model",
    title = "Overlay of Modelled and Empirical Tenure Densities"
  ) +
  theme_minimal()


df_empirical_individual_moments_pooled %>% 
  select(model, y1, y2, tenure1, tenure2, joint_wave2_p_g) %>% 
  unique() %>% 
  pivot_wider(id_cols = c('y1', 'y2', 'tenure1', 'tenure2'), values_from = joint_wave2_p_g, names_from = model) %>% 
  mutate(diff = Final - Initial) %>% 
  filter(y1 == 1, y2 == 1) %>% 
  filter(tenure1 < 8, tenure2 < 8) %>% 
  ggplot(aes(x = tenure1, y = tenure2, colour = diff)) +
  geom_point(size = 1) +
  scale_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Improvement"
  ) +
  labs(
    title = "Joint Density Heatmap",
    x = "Tenure 1",
    y = "Tenure 2",
    fill = "Density"
  ) +
  theme_minimal()
  

df_empirical_individual_moments_pooled %>% 
  select(model, y1, y2, tenure1, tenure2, joint_wave2_p_g) %>% 
  unique() %>% 
  pivot_wider(id_cols = c('y1', 'y2', 'tenure1', 'tenure2'), values_from = joint_wave2_p_g, names_from = model) %>% 
  mutate(diff = Final - Initial) %>% 
  filter(y1 == 1, y2 == 1) %>% 
  ggplot(aes(x = log(tenure1), y = log(tenure2), colour = diff)) +
  geom_point(size = 2) +
  scale_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "p_g"
  ) +
  labs(
    title = "Joint Density Heatmap",
    x = "Tenure 1",
    y = "Tenure 2",
    fill = "Density"
  ) +
  theme_minimal()

#> Density of h2 ----

df_qlfs_tenure_timegap %>% 
  filter(y2 == 0) %>% 
  mutate(
    case = case_when(
      y1 == 0 & y3 == 0 ~ "00",
      y1 == 0 & y3 == 1 ~ "01",
      y1 == 1 & y3 == 0 ~ "10",
      y1 == 1 & y3 == 1 ~ "11"
    )) %>%
  ggplot(aes(x = timegap2, group = case, color = case)) + 
  geom_density()