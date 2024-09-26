# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

# Load functions
source("scripts/define_estimation_functions_4waves.R")

df_estimate <- df_qlfs %>% 
  filter(period1 == 33)
  
# INGEST DATA ====

source("scripts/ingest_data_4waves.R")


# DEFINE OBJECTS ====

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 


# ESTIMATE MODEL ====


# param_init <- data.frame(theta_1 = 0.03, theta_2 = 0.03, pi = 0.03)
# param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.001)
param_init <- data.frame(theta_0 = 0.04010916, theta_01 = 0.2442369, theta_10 = 0.227507, theta_1 = 0.0429738, pi = 0.01278907)

param_init_transformed <- logit_transform(param_init)

calc_ll(param_init_transformed)




estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  grad = calc_lg,
  start = param_init_transformed,
  control = list(tol = 0.001, reltol = 0.001, gradtol = 0.001)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1),
  method = "NR"
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1),
  method = "BFGS"
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1),
  method = "BFGSR"
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1),
  method = "CG"
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxLik(
  calc_ll,
  grad = calc_lg,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1),
  method = "NR"
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

b2 <- logit_inverse(estimate_ml$estimate)
grad <- exp()
vcov(estimate_ml)[2,2]


car::deltaMethod(~1/(1 + exp(-x1)), estimate_ml$estimate[2], vcov(estimate_ml)[2,2])

