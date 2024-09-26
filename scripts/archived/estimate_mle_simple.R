# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)

# EVALUATION PARAMETERS ====

include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

source("scripts/ingest_data.R")


# INGEST DATA ====

source("scripts/define_estimation_functions_3waves.R")



# DEFINE OBJECTS ====

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1)) 

df_estimate <- df_qlfs %>% 
  rename(
    y1 = Status_Q1,
    y2 = Status_Q2,
    y3 = Status_Q3
  ) %>% 
  select(y1, y2, y3) %>% 
  as.data.frame


# ESTIMATE MODEL ====


# param_init <- data.frame(theta_1 = 0.03, theta_2 = 0.03, pi = 0.03)
# param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.001)
param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.01)

param_init_transformed <- logit_transform(param_init)

calc_ll(param_init_transformed)
calc_lg(param_init_transformed)



estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  grad = calc_lg,
  start = param_init_transformed,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum

estimate_ml <- maxLik::maxBFGSR(
  calc_ll,
  start = estimate_ml$estimate,
  control = list(tol = 0.1, reltol = 0.1, gradtol = 0.1)
)

logit_inverse(estimate_ml$estimate)
estimate_ml$maximum
