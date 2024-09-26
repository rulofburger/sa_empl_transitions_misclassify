# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# INGEST DATA ====

source("scripts/ingest_data_4waves.R")

df_estimate <- df_qlfs %>%
  # filter(period1 >= 35 &  period1 <= 40)
  # filter(period1 == 32 | period1 == 33)
  filter(period1 >= 30 & period1 <= 32) #0.01153142


# DEFINE FUNCTIONS ====

# Load functions
source("scripts/define_estimation_functions_4waves_mle_ar2.R")

# Calculates probabilities for general AR2 process with misclassification error
create_probs_temp <- function(param) {
  
  df_probs_temp <- df_template %>% 
    mutate(
      Theta = param$theta_0*(param$theta_10 - 1) + param$theta_1*(param$theta_01 - 1),
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
  
}

#> Mean employment rates ====

# Calculates observed and actual steady-state employment state as function of theta_0 and pi (maintaining theta_1 = 0.03)
calc_s1_steady_state_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    filter(y1 == 1) %>% 
    pull(joint_p)
  
}

#> Two-wave transitions ====

# Calculates observed and actual two-wave transition rates as function of theta_0 and pi 
calc_p1p2_transition_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y2 == 1) %>% 
    filter(y1 == 0) %>% 
    pull(cond_p)
}

#> Three-wave transitions ====

# Calculates observed and actual three-wave transition rates as function of theta_0 and pi 
calc_p1p3_transition_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y3 == 1) %>% 
    filter(y1 == 0) %>% 
    pull(cond_p)
}

#> Two-wave dynamics ====

# Calculates observed and actual joint probability of NN (assuming symmetric transition rates)
calc_s1s2_00_probs_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    filter(y1 == 0) %>% 
    filter(y2 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    pull(joint_p)
}

# Calculates observed and actual joint probability of NE (assuming symmetric transition rates)
calc_s1s2_01_probs_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    filter(y1 == 0) %>% 
    filter(y2 == 1) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    pull(joint_p)
}

# Calculates observed probability of NN and NE
calc_s1s2_probs_theta0_theta1_pi <- function(theta_0, theta_1, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_1
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    filter(y1 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") 
}

#> Three-wave dynamics ====

# Calculates observed probability of NNN, NNE, NEN and NEE
calc_s1s2s3_probs_theta0_pi <- function(theta_0, theta_1, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_1
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2, y3) %>% 
    filter(y1 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") 
}

# DEFINE OBJECTS ====

df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 

# AR2 + ME
param_hat_ar2_me <-  data.frame(theta_0 = 0.036667917, theta_01 = 0.283934383, theta_10 = 0.327332178, theta_1 = 0.041882125, pi = 0.006396909)

# AR2
param_hat_ar2 <-  data.frame(theta_0 = 0.06304608, theta_01 = 0.1864498, theta_10 = 0.3093172, theta_1 = 0.05802624, pi = 0)


# EMPIRICAL ESTIMATES OBJECTS ====

#> Mean employment rates ====

# Tabulate mean employment rate across all 4 waves
table_empl_rates <- df_estimate %>%
  srvyr::as_survey(weights = c(weight)) %>%
  summarise(
    y1 = srvyr::survey_mean(y1, na.rm = T),
    y2 = srvyr::survey_mean(y2, na.rm = T),
    y3 = srvyr::survey_mean(y3, na.rm = T),
    y4 = srvyr::survey_mean(y4, na.rm = T)
  )

tibble(
  `Wave 1 (%)` = c(table_empl_rates[[1]], table_empl_rates[[2]])*100,
  `Wave 2 (%)` = c(table_empl_rates[[3]], table_empl_rates[[4]])*100,
  `Wave 3 (%)` = c(table_empl_rates[[5]], table_empl_rates[[6]])*100,
  `Wave 4 (%)` = c(table_empl_rates[[7]], table_empl_rates[[8]])*100
)

#> Two-wave transitions ====

# Tabulate 2-wave transition rates across all 3 successive wave-pairs
df_2w_transitions <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>%
  group_by(y1) %>% 
  summarise(y2 = srvyr::survey_mean(y2)) %>% 
  bind_cols(
    df_estimate %>% 
      srvyr::as_survey(weights = c(weight)) %>%
      group_by(y2) %>% 
      summarise(y3 = srvyr::survey_mean(y3)) %>% 
      select(y3, y3_se)
  ) %>% 
  bind_cols(
    df_estimate %>% 
      srvyr::as_survey(weights = c(weight)) %>%
      group_by(y3) %>% 
      summarise(y4 = srvyr::survey_mean(y4)) %>% 
      select(y4, y4_se)
  ) %>% 
  mutate(y2 = if_else(y1 == 1, 1 - y2, y2)) %>%
  mutate(y3 = if_else(y1 == 1, 1 - y3, y3)) %>% 
  mutate(y4 = if_else(y1 == 1, 1 - y4, y4)) %>% 
  select(-y1)

df_2w_transitions

tibble(
  `Wave 1 to 2 (%)` = c(df_2w_transitions[[1,1]], df_2w_transitions[[1,2]], df_2w_transitions[[2,1]], df_2w_transitions[[2,2]])*100,
  `Wave 2 to 3 (%)` = c(df_2w_transitions[[1,3]], df_2w_transitions[[1,4]], df_2w_transitions[[2,3]], df_2w_transitions[[2,4]])*100,
  `Wave 3 to 4 (%)` = c(df_2w_transitions[[1,5]], df_2w_transitions[[1,6]], df_2w_transitions[[2,5]], df_2w_transitions[[2,6]])*100
)

# AR1 model, estimated as average two-wave transition rates
param_hat_ar1 <-  data.frame(theta_0 = unlist(df_2w_transitions[1,4]), theta_01 = 0, theta_10 = 0, theta_1 = unlist(df_2w_transitions[2,4]), pi = 0)

#> Three-wave transitions ====

df_3w_transitions <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>%
  group_by(y1) %>% 
  summarise(y3 = srvyr::survey_mean(y3)) %>% 
  bind_cols(
    df_estimate %>% 
      srvyr::as_survey(weights = c(weight)) %>%
      group_by(y2) %>% 
      summarise(y4 = srvyr::survey_mean(y4)) %>% 
      select(y4)
  ) %>% 
  mutate(y3 = if_else(y1 == 1, 1 - y3, y3)) %>% 
  mutate(y4 = if_else(y1 == 1, 1 - y4, y4)) %>% 
  select(-c('y1', 'y3_se')) %>% 
  mutate(mean = rowMeans(., na.rm = TRUE))

df_3w_transitions

table_2w_3w_transitions <- data.frame(
  `Observed 2-wave` = scales::percent(as.numeric(unlist(df_2w_transitions[,1])), accuracy = 0.01),
  `Predicted 3-wave` = scales::percent(as.numeric(c(df_2w_transitions[1,1]*(1 - df_2w_transitions[2,3]) + (1 - df_2w_transitions[1,1])*df_2w_transitions[1,3], df_2w_transitions[2,1]*(1 - df_2w_transitions[1,3]) + (1 - df_2w_transitions[2,1])*df_2w_transitions[2,3])), accuracy = 0.01), 
  `Observed 3-wave` = scales::percent(as.numeric(unlist(df_3w_transitions[,1])), accuracy = 0.01)
) 
rownames(table_2w_3w_transitions) <- c("Job entry", "Job exit")

table_2w_3w_transitions
print(xtable::xtable(table_2w_3w_transitions, type = "latex"), file = 'output/tables/4 waves/table_2w_3w_transitions.tex')


#> Reduced-form estimates from two- and three-wave transitions ====

table_implied_structural_2w_3w <- data.frame(
  `Implied Actual Transition rate` = scales::percent(as.numeric(unlist(df_3w_transitions[,1] - df_2w_transitions[,1])), accuracy = 0.01),
  `Implied Miscl. rate` = scales::percent(as.numeric(unlist(df_2w_transitions[,1] - 0.5*df_3w_transitions[,1])), accuracy = 0.01)
) 
rownames(table_implied_structural_2w_3w) <- c("Job entry", "Job exit")

table_implied_structural_2w_3w
print(xtable::xtable(table_implied_structural_2w_3w, type = "latex"), file = 'output/tables/4 waves/table_implied_structural_2w_3w.tex')

# AR1 + ME
param_hat_ar1_me <-  data.frame(theta_0 = unlist(df_3w_transitions[1,3] - df_2w_transitions[1,3]), theta_01 = 0, theta_10 = 0, theta_1 = unlist(df_3w_transitions[2,3] - df_2w_transitions[2,3]), pi = mean(unlist(df_2w_transitions[,3] - 0.5*df_3w_transitions[,3])))



# EXPOSITIONAL GRAPHS & TABLES ====

#> Mean employment rates ====

figure_steady_state_theta0_pi <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_s1_steady_state_theta0_pi),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
  pi = rep(0:5, each = 11)/100
) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Two-wave job entry rate") + 
  ylab("Observed steady-state employment rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_steady_state_theta0_pi

ggsave(
  paste0('output/figures/4 waves/figure_steady_state_theta0_pi.jpeg'),
  width = 15,
  height = 7,
  plot = figure_steady_state_theta0_pi
)


#> Two-wave transitions ====

figure_theta_pi_w2 <- tibble(
    prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_p1p2_transition_theta0_pi),
    theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
    pi = rep(0:5, each = 11)/100
  ) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual two-wave transition rate") + 
  ylab("Observed two-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_w2

figure_theta_pi_w2 + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.25))

ggsave(
  paste0('output/figures/4 waves/figure_theta_pi_w2.jpeg'),
  width = 15,
  height = 7,
  plot = figure_theta_pi_w2
)


#> Three-wave transitions ====

figure_theta_pi_w3 <- tibble(
    prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_p1p3_transition_theta0_pi),
    theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
    pi = rep(0:5, each = 11)/100
  ) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual three-wave transition rate") + 
  ylab("Observed three-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_w3

ggsave(
  paste0('output/figures/4 waves/figure_theta_pi_w3.jpeg'),
  width = 15,
  height = 7,
  plot = figure_theta_pi_w3
)

#> Compare two- and three-wave transitions ====

figure_theta_pi_w2w3 <- ggpubr::ggarrange(
  figure_theta_pi_w2 + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.25)),
  figure_theta_pi_w3 + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.25)),
  legend = "right",
  common.legend = TRUE
)

figure_theta_pi_w2w3

ggsave(
  paste0('output/figures/4 waves/figure_theta_pi_w2w3.jpeg'),
  width = 15,
  height = 7,
  plot = figure_theta_pi_w2w3
)

#> Two-wave dynamics ====

# Calculates observed and actual joint probability of NN (assuming symmetric transition rates)
figure_theta_pi_s1s200 <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_s1s2_00_probs_theta0_pi),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
  pi = rep(0:5, each = 11)/100
) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual three-wave transition rate") + 
  ylab("Observed three-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_s1s200

# Calculates observed and actual joint probability of NE (assuming symmetric transition rates)
figure_theta_pi_s1s201 <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_s1s2_01_probs_theta0_pi),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
  pi = rep(0:5, each = 11)/100
) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual three-wave transition rate") + 
  ylab("Observed three-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_s1s201

# AR1 model with no misclassification
calc_s1s2_probs_theta0_theta1_pi(0.03, 0.03, 0)

# AR1 model with 3% misclassification
calc_s1s2_probs_theta0_theta1_pi(0.03, 0.03, 0.03)

calc_s1s2_probs_theta0_theta1_pi(0.03, 0.03, 0.03)$joint_p/calc_s1s2_probs_theta0_theta1_pi(0.03, 0.03, 0)$joint_p

#> Three-wave dynamics ====

calc_s1s2s3_probs_theta0_pi <- function(theta_0, theta_1, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_1
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2, y3) %>% 
    filter(y1 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") 
}

# AR1 model with no misclassification
calc_s1s2s3_probs_theta0_pi(0.03, 0.03, 0)

# AR1 model with 3% misclassification
calc_s1s2s3_probs_theta0_pi(0.03, 0.03, 0.03)

calc_s1s2s3_probs_theta0_pi(0.03, 0.03, 0.03)$joint_p/calc_s1s2s3_probs_theta0_pi(0.03, 0.03, 0)$joint_p

calc_s1s2s3_conNNE_probs_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2, y3) %>% 
    filter(y1 == 0) %>% 
    filter(y2 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1, y2) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y3 == 1) %>% 
    pull(cond_p)
}

figure_theta_pi_NNE <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_s1s2s3_conNNE_probs_theta0_pi),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
  pi = rep(0:5, each = 11)/100
) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual two-wave transition rate") + 
  ylab("Observed third-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_NNE

ggsave(
  paste0('output/figures/4 waves/figure_theta_pi_NNE.jpeg'),
  width = 15,
  height = 7,
  plot = figure_theta_pi_NNE
)

calc_s1s2s3_conENE_probs_theta0_pi <- function(theta_0, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2, y3) %>% 
    filter(y1 == 1) %>% 
    filter(y2 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1, y2) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y3 == 1) %>% 
    pull(cond_p)
}

figure_theta_pi_ENE <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_s1s2s3_conENE_probs_theta0_pi),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
  pi = rep(0:5, each = 11)/100
) %>%
  mutate(pi = as.character(paste0(pi*100, " %"))) %>% 
  mutate(prob = if_else(is.na(prob), 0, prob)) %>% 
  ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) + 
  geom_line() + 
  xlab("Actual two-wave transition rate") + 
  ylab("Observed third-wave transition rate") + 
  labs(colour = "Miscl. prob.") + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0,0.1,0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(text = element_text(size = 20))

figure_theta_pi_ENE

ggsave(
  paste0('output/figures/4 waves/figure_theta_pi_ENE.jpeg'),
  width = 15,
  height = 7,
  plot = figure_theta_pi_ENE
)

# AR1 model with no misclassification
calc_s1s2s3_probs_theta0_pi(unlist(df_2w_transitions[1,4]), unlist(df_2w_transitions[2,4]), 0)

# AR1 model with misclassification
calc_s1s2s3_probs_theta0_pi(
  as.numeric(unlist(df_3w_transitions[1, 3] - df_2w_transitions[1, 3])),
  as.numeric(unlist(df_3w_transitions[2, 3] - df_2w_transitions[2, 3])),
  mean(as.numeric(unlist(df_2w_transitions[, 3] - 0.5*df_3w_transitions[, 3])))
  )

# Empirical estimation
df_s1s2s3_probs <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>% 
  group_by(y1, y2, y3) %>% 
  summarise(n = srvyr::survey_total()) %>% 
  ungroup() %>% 
  mutate(tot = sum(n)) %>% 
  mutate(p = n/tot) %>%
  filter(y1 == 0) %>% 
  select(y1, y2, y3, p)
  
df_s1s2s3_probs

# AR2 model with no misclassification
create_probs_temp(param_hat_ar2) %>% 
  group_by(y1, y2, y3) %>% 
  filter(y1 == 0) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") 

# AR2 model with misclassification
create_probs_temp(param_hat_ar2_me) %>% 
  group_by(y1, y2, y3) %>% 
  filter(y1 == 0) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") 
   



#> Four-wave transitions ====

df_4w_transitions <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>%
  group_by(y1) %>% 
  summarise(y4 = srvyr::survey_mean(y4)) %>% 
  mutate(y4 = if_else(y1 == 1, 1 - y4, y4)) %>% 
  select(-c('y1', 'y4_se'))

df_4w_transitions

calc_s1s2s3s4_probs_theta0_theta1_pi <- function(theta_0, theta_1, pi) {
  
  param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_1
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2, y3, y4) %>% 
    filter(y1 == 0) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") 
}

# AR1 model with no misclassification
calc_s1s2s3s4_probs_theta0_theta1_pi(0.03, 0.03, 0)

# AR1 model with 3% misclassification
calc_s1s2s3s4_probs_theta0_theta1_pi(0.03, 0.03, 0.03)

calc_s1s2s3s4_probs_theta0_theta1_pi(0.03, 0.03, 0.03)$joint_p/calc_s1s2s3s4_probs_theta0_theta1_pi(0.03, 0.03, 0)$joint_p


# AR1 model with no misclassification
create_probs_temp(param_hat_ar1) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4) %>% 
  select(y1, y4, cond_p)

# AR1 model with misclassification
create_probs_temp(param_hat_ar1_me) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4) %>% 
  select(y1, y4, cond_p)

# AR2 model with no misclassification
create_probs_temp(param_hat_ar2) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4) %>% 
  select(y1, y4, cond_p)

# AR2 model with misclassification
create_probs_temp(param_hat_ar2_me) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4) %>% 
  select(y1, y4, cond_p)

table_2w_4w_transitions <- tibble(
  # `Observed 2-wave` = scales::percent(as.numeric(unlist(df_2w_transitions[,4])), accuracy = 0.01),
  `Predicted AR1` = scales::percent(as.numeric(c(
    (1 - df_2w_transitions[1,4])*(1 - df_2w_transitions[1,4])*df_2w_transitions[1,4] + 
      (1 - df_2w_transitions[1,4])*df_2w_transitions[1,4]*(1 - df_2w_transitions[2,4]) + 
      df_2w_transitions[1,4]*df_2w_transitions[2,4]*df_2w_transitions[1,4] + 
      df_2w_transitions[1,4]*(1 - df_2w_transitions[2,4])*(1 - df_2w_transitions[2,4]),
    (1 - df_2w_transitions[2,4])*(1 - df_2w_transitions[2,4])*df_2w_transitions[2,4] + 
      (1 - df_2w_transitions[2,4])*df_2w_transitions[2,4]*(1 - df_2w_transitions[1,4]) + 
      df_2w_transitions[2,4]*df_2w_transitions[1,4]*df_2w_transitions[2,4] + 
    df_2w_transitions[2,4]*(1 - df_2w_transitions[1,4])*(1 - df_2w_transitions[1,4]))), accuracy = 0.01), 
  `Observed 4-wave` = scales::percent(as.numeric(unlist(df_4w_transitions[,1])), accuracy = 0.01)
)
table_2w_4w_transitions

param <- data.frame(theta_0 = 0.0736, theta_01 = 0, theta_10 = 0, theta_1 = 0.0744, pi = 0)
create_probs_temp(param) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4)


param <- data.frame(theta_0 = 0.023, theta_01 = 0, theta_10 = 0, theta_1 = 0.023, pi = 0.0255)
create_probs_temp(param) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>% 
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 != y4)

figure_theta0_theta01_contour <- data.frame(
    theta_0 = rep(seq(0, 0.1, 0.01), 11),
    beta2 = rep(0:10, each = 11)/100
  ) %>% 
  mutate(theta_01 = 1 - theta_0/beta2) %>% 
  mutate(theta_01 = if_else(theta_01 < 0, NA_real_, theta_01)) %>% 
  mutate(beta2 = as.factor(beta2)) %>% 
  ggplot(aes(x = theta_0, y = theta_01, group = beta2, color = beta2)) + geom_line()
  
figure_theta0_theta01_contour


# theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
# pi = rep(0:5, each = 11)/100




#> AR(2)-dynamics ====

# Two-wave joint probs

calc_s1s2_probs_theta0_theta_01_pi <- function(theta_0, theta_01, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$theta_01 <- theta_01
  param$theta_10 <- theta_01
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    filter(y1 == 0 & y2 == 1) %>% 
    pull(joint_p)
}

figure_s1s2_probs_theta0_theta_01_pi0 <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)), 9), rep(seq(0,40,5), each = 11)/100, calc_s1s2_probs_theta0_theta_01_pi, pi = 0),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)), 9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  # mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = scales::percent(theta_01), color = scales::percent(theta_01))) +
  geom_line() + 
  xlab(expression(theta[0])) + 
  ylab("True probability of N*E*") + 
  labs(colour = expression(theta["01"])) + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.1, 0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.12)) + 
  theme(text = element_text(size = 20))

figure_s1s2_probs_theta0_theta_01_pi0

# ggsave(
#   paste0('output/figures/4 waves/figure_s1s2_probs_theta0_theta_01_pi0.jpeg'),
#   width = 15,
#   height = 7,
#   plot = figure_s1s2_probs_theta0_theta_01_pi0
# )

figure_s1s2_probs_theta0_theta_01_pi03 <- tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)), 9), rep(seq(0,40,5), each = 11)/100, calc_s1s2_probs_theta0_theta_01_pi, pi = 0.03),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)), 9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  # mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = scales::percent(theta_01), color = scales::percent(theta_01))) +
  geom_line() + 
  xlab(expression(theta[0])) + 
  ylab("Observed probability of NE") + 
  labs(colour = expression(theta["01"])) + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.1, 0.01)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.12)) + 
  theme(text = element_text(size = 20))

figure_s1s2_probs_theta0_theta_01_pi03

# ggsave(
#   paste0('output/figures/4 waves/figure_s1s2_probs_theta0_theta_01_pi03.jpeg'),
#   width = 15,
#   height = 7,
#   plot = figure_s1s2_probs_theta0_theta_01_pi03
# )

figure_s1s2_probs_theta0_theta_01_pi <- ggpubr::ggarrange(
  figure_s1s2_probs_theta0_theta_01_pi0,
  figure_s1s2_probs_theta0_theta_01_pi03,
  legend = "right",
  common.legend = TRUE
)

figure_s1s2_probs_theta0_theta_01_pi

ggsave(
  paste0('output/figures/4 waves/figure_s1s2_probs_theta0_theta_01_pi.jpeg'),
  width = 15,
  height = 7,
  plot = figure_s1s2_probs_theta0_theta_01_pi
)


param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0)
s1s2_ar2_0 <- create_probs_temp(param) %>% 
  group_by(y1, y2) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")
s1s2_ar2_0

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0.03)
s1s2_ar2_pi <- create_probs_temp(param) %>% 
  group_by(y1, y2) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")
s1s2_ar2_pi
(s1s2_ar2_pi %>% pull(joint_p))/(s1s2_ar2_0 %>% pull(joint_p))


#Two-wave transition rate

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0)
create_probs_temp(param) %>% 
  group_by(y1, y2) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y2 == 1) %>% 
  select(y1, y2, cond_p)

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.04, pi = 0.03)
create_probs_temp(param) %>% 
  group_by(y1, y2) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y2 == 1) %>% 
  select(y1, y2, cond_p)

#Three-wave transition rate

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0)
create_probs_temp(param) %>% 
  group_by(y1, y3) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y3 == 1) %>% 
  select(y1, y3, cond_p)

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0.03)
create_probs_temp(param) %>% 
  group_by(y1, y3) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y3 == 1) %>% 
  select(y1, y3, cond_p)


# Four-wave transition rate

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0)
create_probs_temp(param) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y4 == 1) %>% 
  select(y1, y4, cond_p)

param <- data.frame(theta_0 = 0.03, theta_01 = 0.3, theta_10 = 0.3, theta_1 = 0.03, pi = 0.03)
create_probs_temp(param) %>% 
  group_by(y1, y4) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
  group_by(y1) %>% 
  mutate(full_p = sum(joint_p)) %>%
  mutate(cond_p = joint_p/full_p) %>% 
  filter(y1 == 0, y4 == 1) %>% 
  select(y1, y4, cond_p)




calc_p1p2_transition_theta0_theta_01_pi <- function(theta_0, theta_01, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$theta_01 <- theta_01
  param$theta_10 <- theta_01
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y2) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y2 == 1) %>% 
    filter(y1 == 0) %>% 
    pull(cond_p)
  
}

calc_p1p3_transition_theta0_theta_01_pi <- function(theta_0, theta_01, pi) {
  
  param <- data.frame(theta_0 = 0.03, theta_01 = 0, theta_10 = 0, theta_1 = 0.03, pi = 0)
  param$theta_0 <- theta_0
  param$theta_1 <- theta_0
  param$theta_01 <- theta_01
  param$theta_10 <- theta_01
  param$pi <- pi
  
  create_probs_temp(param) %>% 
    group_by(y1, y3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") %>% 
    group_by(y1) %>% 
    mutate(full_p = sum(joint_p)) %>% 
    mutate(cond_p = joint_p/full_p) %>% 
    filter(y3 == 1) %>% 
    filter(y1 == 0) %>% 
    pull(cond_p)
  
}



# tibble(
#     prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),6), rep(0:5, each = 11)/100, calc_p1p3_transition),
#     theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),6),
#     pi = rep(0:5, each = 11)/100
#   ) %>%
#   mutate(pi = as.factor(pi)) %>%
#   mutate(prob = if_else(is.na(prob), 0, prob)) %>%
#   ggplot(aes(x = theta_0, y = prob, group = pi, color = pi)) +
#   geom_line() +
#   xlab("Actual two-wave transition rate") +
#   ylab("Observed two-wave transition rate") +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1))





tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)), 9), rep(seq(0,40,5), each = 11)/100, calc_p1p2_transition_theta0_theta_01_pi, pi = 0),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)), 9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = theta_01, color = theta_01)) +
  geom_line()

tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)), 9), rep(seq(0,40,5), each = 11)/100, calc_p1p2_transition_theta0_theta_01_pi, pi = 0.03),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)), 9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = theta_01, color = theta_01)) +
  geom_line()

tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),9), rep(seq(0,40,5), each = 11)/100, calc_p1p3_transition_theta0_theta_01_pi, pi = 0),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = theta_01, color = theta_01)) +
  geom_line()

tibble(
  prob = map2_dbl(rep(c(0.0001,seq(0.01,0.1,0.01)),9), rep(seq(0,40,5), each = 11)/100, calc_p1p3_transition_theta0_theta_01_pi, pi = 0.03),
  theta_0 = rep(c(0.0001,seq(0.01,0.1,0.01)),9),
  theta_01 = rep(seq(0,40,5), each = 11)/100
) %>%
  mutate(theta_01 = as.factor(theta_01)) %>%
  mutate(prob = if_else(is.na(prob), 0, prob)) %>%
  ggplot(aes(x = theta_0, y = prob, group = theta_01, color = theta_01)) +
  geom_line()

# Empirical estimation
df_s1s2p3_transitions <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>% 
  group_by(y1, y2) %>% 
  summarise(p = srvyr::survey_mean(y3)) %>% 
  mutate(p = if_else(y2 == 1, 1 - p, p)) %>% 
  arrange(y2, y1) %>% 
  select(y1, y2, p)

df_s1s2p3_transitions

param <- data.frame(theta_0 = df_s1s2p3_transitions$p[1], theta_01 = df_s1s2p3_transitions$p[2] - df_s1s2p3_transitions$p[1], theta_10 = df_s1s2p3_transitions$p[3] - df_s1s2p3_transitions$p[4], theta_1 = df_s1s2p3_transitions$p[4], pi = 0)
create_probs_temp(param) %>% 
  group_by(y1, y2) %>% 
  summarise(joint_p = sum(joint_p), .groups = "drop")


df_s1s2_probs <- df_estimate %>% 
  srvyr::as_survey(weights = c(weight)) %>% 
  group_by(y1, y2) %>% 
  summarise(n = srvyr::survey_total()) %>% 
  ungroup() %>% 
  mutate(tot = sum(n)) %>% 
  mutate(p = n/tot) %>%
  # filter(y1 == 0) %>% 
  select(y1, y2, p)

df_s1s2_probs