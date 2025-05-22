# METADATA =====================================================================
# DESCRIPTION: Estimate 3 wave models for SA
# CREATED: 2025-02-13 (rulofburger)

# SUMMARY: This script creates the estimates for the 3 wave SA employment
# transition models

#-------------------------------------------------------------------------------
# 1) INITIALISE-----------------------------------------------------------------
#-------------------------------------------------------------------------------
# Source environment/functions

# Load libraries
library(tidyverse)
library(haven)
library(data.table)
library(fastverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the
# parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

#> Load estimation functions defined in other scripts ----
source("scripts/define_estimation_functions_3waves_gmm_ar1_fmm.R")


# INGEST DATA ====

# Run script that loads 3 wave SA data as df_qlfs
source("scripts/ingest_data_3waves_SA.R")

# Limit survey rounds and calculate weights to be consistent within panel and to sum to 1
df_qlfs <- df_qlfs %>% 
  filter(period1 >= 30 & period1 <= 32) %>% 
  mutate(weight_total = sum(weight))  %>% 
  mutate(weight = dim(df_qlfs)[1]*weight/weight_total) 

# ESTIMATE MODELS ====

df_template <- data.table::CJ(
  y1      = c(0, 1), 
  y1_star = c(0, 1), 
  y2      = c(0, 1), 
  y2_star = c(0, 1), 
  y3      = c(0, 1), 
  y3_star = c(0, 1)
)

df_estimate <- df_qlfs

gmm_moments_wrapper <- function(param_unlist, data) {
  param_init <- data.frame(
    theta0_1 = param_unlist[1], 
    theta1_1 = param_unlist[2], 
    p_1 = param_unlist[3],
    theta0_2 = param_unlist[4], 
    theta1_2 = param_unlist[5], 
    pi = param_unlist[6]
  )
  gmm_moments(param_init)
}


# 1) Your named logit‐scale start
param_init_vec <- c(
  theta0_1 = 0.085,
  theta1_1 = 0.085,
  p_1      = 0.10,
  theta0_2 = 0.03,
  theta1_2 = 0.03,
  pi       = 0.03
)
param_init_transformed <- logit_transform(param_init_vec)

# 2) Wrapper returning a pure numeric matrix
gmm_moments_wrapper <- function(theta_tr, x) {
  names(theta_tr) <- names(param_init_vec)
  df_gi <- calc_lli_derivatives_3waves_ar1_fmm2(theta_tr, pi0 = FALSE)
  data.matrix(df_gi, rownames.force = FALSE)
}


# 3) Create a dummy numeric 'x' of the correct length
n <- nrow(df_estimate)
dummy_x <- matrix(1, n, 1)   # just a column of 1s

# 4) One‐step GMM with W=I, no pre‐whitening, i.i.d. vcov
res1 <- gmm::gmm(
  g        = gmm_moments_wrapper,
  x        = dummy_x,                # pure numeric matrix
  t0       = param_init_transformed,
  type     = "twoStep",              # identity in step 1, step 2 is skipped since S is singular
  wmatrix  = "ident",                # W = I
  vcov     = "iid",                  # use i.i.d. covariances
  prewhite = FALSE,                  # disable AR(1) pre‐whitening
  bw       = 0                       # no bandwidth smoothing
)

summary(res1)

# extract logit‐scale estimates and their SEs
logit_est  <- coef(res1)
logit_se   <- summary(res1)$coef[, "Std. Error"]

# back‐transform
prob_est   <- plogis(logit_est)
# delta‐method for SE: SE_prob ≈ se_logit * p*(1−p)
prob_se    <- logit_se * prob_est * (1 - prob_est)

# tabulate
results <- data.frame(
  param      = names(logit_est),
  logit_est  = logit_est,
  logit_se   = logit_se,
  prob_est   = prob_est,
  prob_se    = prob_se
)

print(results)
