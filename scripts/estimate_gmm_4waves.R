# INITIALISE ====
# Source environment/functions

# Load libraries
library(tidyverse)


# EVALUATION PARAMETERS ====

# include_mu <- FALSE # Include mu as one of the independent parameters in the parameter vector (as opposed to having it be derived form theta_1 and theta_2)

# DEFINE FUNCTIONS ====

# Load functions
source("scripts/define_estimation_functions_4waves_gmm.R")


df_template <- data.table::CJ(y1 = c(0, 1), y1_star = c(0, 1), y2 = c(0, 1), y2_star = c(0, 1), y3 = c(0, 1), y3_star = c(0, 1), y4 = c(0, 1), y4_star = c(0, 1)) 


# INGEST DATA ====

source("scripts/ingest_data_4waves.R")

df_estimate <- df_qlfs %>%
  # filter(period1 == 33)
#   filter(period1 == 22)
filter(period1 >= 30 | period1 <= 31)

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
# param_init <- data.frame(theta_0 = 0.041898600, theta_01 = 0.308669212, theta_10 = 0.366688229, theta_1 = 0.044556080, pi = 0.009508421)

param_init_transformed <- logit_transform(param_init)
# gmm_init <- calc_gmm_moments(param_init_transformed, df_x)
# #
# sum(colMeans(gmm_init)^2)
# sum(colMeans(gmm_init))

model_gmm <- gmm::gmm(calc_gmm_moments, df_x, unlist(param_init_transformed), gradv = calc_gmm_derivatives, traceIter = TRUE, prewhite = 0)

print(logit_inverse(model_gmm$coefficients))

# gmm_hat <- calc_gmm_moments(res$coefficients, df_x)

# sum(colMeans(gmm_hat)^2)

summary(model_gmm)
print(gmm::specTest(model_gmm))

estimate_gmm_period <- function(period_begin, period_end) {
  
  print(paste0("Begin period:", period_begin))
  print(paste0("End period:", period_end))
  
  df_estimate <- df_qlfs %>% 
    filter(period1 >= period_begin & period1 <= period_end)
  
  # sd_unempl <- df_estimate %>%
  #   srvyr::as_survey(weights = c(weight)) %>%
  #   summarise(
  #     y1 = srvyr::survey_mean(y1, na.rm = T),
  #     y2 = srvyr::survey_mean(y2, na.rm = T),
  #     y3 = srvyr::survey_mean(y3, na.rm = T),
  #     y4 = srvyr::survey_mean(y4, na.rm = T)
  #   ) %>% 
  #   select(y1, y2, y3, y4) %>% 
  #   unlist() %>% as.vector() %>% 
  #   sd()
  # 
  # df_2w_transitions <- df_estimate %>% 
  #   srvyr::as_survey(weights = c(weight)) %>%
  #   group_by(y1) %>% 
  #   summarise(y2 = srvyr::survey_mean(y2)) %>% 
  #   bind_cols(
  #     df_estimate %>% 
  #       srvyr::as_survey(weights = c(weight)) %>%
  #       group_by(y2) %>% 
  #       summarise(y3 = srvyr::survey_mean(y3)) %>% 
  #       select(y3)
  #   ) %>% 
  #   bind_cols(
  #     df_estimate %>% 
  #       srvyr::as_survey(weights = c(weight)) %>%
  #       group_by(y3) %>% 
  #       summarise(y4 = srvyr::survey_mean(y4)) %>% 
  #       select(y4)
  #   ) %>% 
  #   mutate(y2 = if_else(y1 == 1, 1 - y2, y2)) %>%
  #   mutate(y3 = if_else(y1 == 1, 1 - y3, y3)) %>% 
  #   mutate(y4 = if_else(y1 == 1, 1 - y4, y4)) %>% 
  #   select(-c('y1', 'y2_se')) %>% 
  #   mutate(mean = rowMeans(., na.rm = TRUE))
  # 
  # sd_entry_w2 <- sd(df_2w_transitions[1,1:3])
  # sd_exit_w2 <- sd(df_2w_transitions[2,1:3])
  # 
  # df_3w_transitions <- df_estimate %>% 
  #   srvyr::as_survey(weights = c(weight)) %>%
  #   group_by(y1) %>% 
  #   summarise(y3 = srvyr::survey_mean(y3)) %>% 
  #   bind_cols(
  #     df_estimate %>% 
  #       srvyr::as_survey(weights = c(weight)) %>%
  #       group_by(y2) %>% 
  #       summarise(y4 = srvyr::survey_mean(y4)) %>% 
  #       select(y4)
  #   ) %>% 
  #   mutate(y3 = if_else(y1 == 1, 1 - y3, y3)) %>% 
  #   mutate(y4 = if_else(y1 == 1, 1 - y4, y4)) %>% 
  #   select(-c('y1', 'y3_se')) %>% 
  #   mutate(mean = rowMeans(., na.rm = TRUE))
  # 
  # sd_entry_w3 <- sd(df_3w_transitions[1,1:2])
  # sd_exit_w3 <- sd(df_3w_transitions[2,1:2])
  # 
  # df_ar2 <- data.frame(
  #   s1s2 =   df_estimate %>% 
  #     srvyr::as_survey(weights = c(weight)) %>%
  #     group_by(y1, y2) %>% 
  #     summarise(freq = srvyr::survey_total()) %>% 
  #     ungroup() %>% 
  #     mutate(total = sum(freq)) %>% 
  #     mutate(freq = freq/total) %>% 
  #     select(freq),
  #   s2s3 =   df_estimate %>% 
  #     srvyr::as_survey(weights = c(weight)) %>%
  #     group_by(y2, y3) %>% 
  #     summarise(freq = srvyr::survey_total()) %>% 
  #     ungroup() %>% 
  #     mutate(total = sum(freq)) %>% 
  #     mutate(freq = freq/total) %>% 
  #     select(freq),
  #   s3s4 =   df_estimate %>% 
  #     srvyr::as_survey(weights = c(weight)) %>%
  #     group_by(y3, y4) %>% 
  #     summarise(freq = srvyr::survey_total()) %>% 
  #     ungroup() %>% 
  #     mutate(total = sum(freq)) %>% 
  #     mutate(freq = freq/total) %>% 
  #     select(freq)
  # )
  # 
  # sd_entry_s00 <- sd(df_ar2[1,1:3])
  # sd_entry_s01 <- sd(df_ar2[2,1:3])
  # sd_entry_s10 <- sd(df_ar2[3,1:3])
  # sd_entry_s11 <- sd(df_ar2[4,1:3])

  
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
    select(y0000, y0001, y0010, y0011, y0100, y0101, y0110, y0111, y1000, y1001, y1010, y1011, y1100, y1101, y1110, y1111) %>% 
    as.matrix()

  # ESTIMATE MODEL ====
  
  
  # param_init <- data.frame(theta_1 = 0.03, theta_2 = 0.03, pi = 0.03)
  # param_init <- data.frame(theta_1 = 0.085, theta_2 = 0.085, pi = 0.001)
  param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615, pi = 0.01351413)
  # param_init <- data.frame(theta_0 = 0.041898600, theta_01 = 0.308669212, theta_10 = 0.366688229, theta_1 = 0.044556080, pi = 0.009508421)
  
  param_init_transformed <- logit_transform(param_init)
  # gmm_init <- calc_gmm_moments(param_init_transformed, df_x)
  # 
  # sum(colMeans(gmm_init)^2)
  
  estimate_gmm <- gmm::gmm(calc_gmm_moments, df_x, unlist(param_init_transformed), traceIter = TRUE, prewhite = 0)
  
  
  print(logit_inverse(estimate_gmm$coefficients))
  
  b2 <- logit_inverse(estimate_gmm$coefficients)
  grad <- exp()
  vcov(estimate_gmm)[2,2]
  
  coeff_list <- list(~1/(1 + exp(-x1)), ~1/(1 + exp(-x1)), ~1/(1 + exp(-x2)), ~1/(1 + exp(-x3)), ~1/(1 + exp(-x4)), ~1/(1 + exp(-x5)))
  

  
  se_list<- list(
    car::deltaMethod(estimate_gmm, "1/(1 + exp(-theta_0))")$SE,
    car::deltaMethod(estimate_gmm, "1/(1 + exp(-theta_01))")$SE,
    car::deltaMethod(estimate_gmm, "1/(1 + exp(-theta_10))")$SE,
    car::deltaMethod(estimate_gmm, "1/(1 + exp(-theta_1))")$SE,
    car::deltaMethod(estimate_gmm, "1/(1 + exp(-pi))")$SE
    )
  
  # gmm_hat <- calc_gmm_moments(res$coefficients, df_x)
  
  # sum(colMeans(gmm_hat)^2)
  
  summary(res)
  print(gmm::specTest(res))
  
  # object <- list(
  #   period_begin = period_begin,
  #   period_end = period_end,
  #   res = res, 
  #   sd_unempl = sd_unempl,
  #   sd_entry_w2 = sd_entry_w2,
  #   sd_exit_w2 = sd_exit_w2,
  #   sd_entry_w3 = sd_entry_w3,
  #   sd_exit_w3 = sd_exit_w3,
  #   sd_entry_s00 = sd_entry_s00,
  #   sd_entry_s01 = sd_entry_s01,
  #   sd_entry_s10 = sd_entry_s10,
  #   sd_entry_s11 = sd_entry_s11,
  #   obs = dim(df_x)[1]
  # )
  
  # return(object)
  }


# numlist <- names(table(df_qlfs$period1))
# counter <- 0
# for (i_begin in 1:43) {
#   for (i_end in i_begin:43) {
#     counter <- counter + 1
#     n_begin <- numlist[i_begin]
#     n_end <- numlist[i_end]
#     
#     # print(n_begin:n_end)
#     
#     result_object <- estimate_gmm_period(period_begin = n_begin, period_end = n_end)
# 
#     results_list[[counter]] <- result_object
#     names(results_list)[counter] <- paste0("result_", n_begin, "_", n_end)
#     
#   }
# }

# # result_35_46 <- estimate_gmm_period(period_begin = 35, period_end = 46)
# # result_35_42 <- estimate_gmm_period(period_begin = 35, period_end = 42)
# # result_35_35 <- estimate_gmm_period(period_begin = 35, period_end = 35)
# # result_36_36 <- estimate_gmm_period(period_begin = 36, period_end = 36)
# # result_37_37 <- estimate_gmm_period(period_begin = 37, period_end = 37)
# # result_38_38 <- estimate_gmm_period(period_begin = 38, period_end = 38)
# # result_39_39 <- estimate_gmm_period(period_begin = 39, period_end = 39)
# # result_40 <- estimate_gmm_period(period_begin = 40, period_end = 40)
# # result_41 <- estimate_gmm_period(period_begin = 41, period_end = 41)
# # result_42 <- estimate_gmm_period(period_begin = 42, period_end = 42)
# # result_43 <- estimate_gmm_period(period_begin = 43, period_end = 43)
# # result_44 <- estimate_gmm_period(period_begin = 44, period_end = 44)
# # result_45 <- estimate_gmm_period(period_begin = 45, period_end = 45)
# # result_9_21 <- estimate_gmm_period(period_begin = 9, period_end = 21)
# result_1 <- estimate_gmm_period(period_begin = 1, period_end = 1)
# result_2 <- estimate_gmm_period(period_begin = 2, period_end = 2)
# result_3 <- estimate_gmm_period(period_begin = 3, period_end = 3)
# result_4 <- estimate_gmm_period(period_begin = 4, period_end = 4)
# result_5 <- estimate_gmm_period(period_begin = 5, period_end = 5)
# result_6 <- estimate_gmm_period(period_begin = 6, period_end = 6)
# result_7 <- estimate_gmm_period(period_begin = 7, period_end = 7)
# result_8 <- estimate_gmm_period(period_begin = 8, period_end = 8)
# result_9 <- estimate_gmm_period(period_begin = 9, period_end = 9)
# result_10 <- estimate_gmm_period(period_begin = 10, period_end = 10)
# result_11 <- estimate_gmm_period(period_begin = 11, period_end = 11)
# result_12 <- estimate_gmm_period(period_begin = 12, period_end = 12)
# result_13 <- estimate_gmm_period(period_begin = 13, period_end = 13)
# result_14 <- estimate_gmm_period(period_begin = 14, period_end = 14)
# result_15 <- estimate_gmm_period(period_begin = 15, period_end = 15)
# result_16 <- estimate_gmm_period(period_begin = 16, period_end = 16)
# result_17 <- estimate_gmm_period(period_begin = 17, period_end = 17)
# result_18 <- estimate_gmm_period(period_begin = 18, period_end = 18)
# result_19 <- estimate_gmm_period(period_begin = 19, period_end = 19)
# result_20 <- estimate_gmm_period(period_begin = 20, period_end = 20)
# result_21 <- estimate_gmm_period(period_begin = 21, period_end = 21)
# result_22 <- estimate_gmm_period(period_begin = 22, period_end = 22)
# result_23 <- estimate_gmm_period(period_begin = 23, period_end = 23)
# result_24 <- estimate_gmm_period(period_begin = 24, period_end = 24)
# result_25 <- estimate_gmm_period(period_begin = 25, period_end = 25)
# 
# 
# result_29 <- estimate_gmm_period(period_begin = 29, period_end = 29)
# result_30 <- estimate_gmm_period(period_begin = 30, period_end = 30)
# result_31 <- estimate_gmm_period(period_begin = 31, period_end = 31)
# result_32 <- estimate_gmm_period(period_begin = 32, period_end = 32)
# result_33 <- estimate_gmm_period(period_begin = 33, period_end = 33)
# result_34 <- estimate_gmm_period(period_begin = 34, period_end = 34)
# result_35 <- estimate_gmm_period(period_begin = 35, period_end = 35)
# result_36 <- estimate_gmm_period(period_begin = 36, period_end = 36)
# result_37 <- estimate_gmm_period(period_begin = 37, period_end = 37)
# result_38 <- estimate_gmm_period(period_begin = 38, period_end = 38)
# result_39 <- estimate_gmm_period(period_begin = 39, period_end = 39)
# result_40 <- estimate_gmm_period(period_begin = 40, period_end = 40)
# result_41 <- estimate_gmm_period(period_begin = 41, period_end = 41)
# result_42 <- estimate_gmm_period(period_begin = 42, period_end = 42)
# result_43 <- estimate_gmm_period(period_begin = 43, period_end = 43)
# result_44 <- estimate_gmm_period(period_begin = 44, period_end = 44)
# result_45 <- estimate_gmm_period(period_begin = 45, period_end = 45)
# result_46 <- estimate_gmm_period(period_begin = 46, period_end = 46)
# 
# 
# 
# 
# result_13_14 <- estimate_gmm_period(period_begin = 13, period_end = 14)
# result_13_14 <- estimate_gmm_period(period_begin = 14, period_end = 15)
# result_12_14 <- estimate_gmm_period(period_begin = 12, period_end = 14)
# result_11_14 <- estimate_gmm_period(period_begin = 11, period_end = 14)
# result_10_14 <- estimate_gmm_period(period_begin = 10, period_end = 14)
# 
# result_14_16 <- estimate_gmm_period(period_begin = 14, period_end = 15)
# result_16_17 <- estimate_gmm_period(period_begin = 16, period_end = 17)
# result_17_18 <- estimate_gmm_period(period_begin = 17, period_end = 18)
# 
# result_31_33 <- estimate_gmm_period(period_begin = 31, period_end = 33)
# 
# result_43_44 <- estimate_gmm_period(period_begin = 43, period_end = 44)
# 
# result_45_46 <- estimate_gmm_period(period_begin = 45, period_end = 46)
# 
# result_31_32 <- estimate_gmm_period(period_begin = 31, period_end = 32)
# 
# # 
# # result_19_21 <- estimate_gmm_period(period_begin = 19, period_end = 21)
# # result_20_21 <- estimate_gmm_period(period_begin = 20, period_end = 21)
# # 
# # result_17_19 <- estimate_gmm_period(period_begin = 17, period_end = 19)
# # result_18_19 <- estimate_gmm_period(period_begin = 18, period_end = 19)
# # result_17_21 <- estimate_gmm_period(period_begin = 17, period_end = 21)
# # result_18_21 <- estimate_gmm_period(period_begin = 18, period_end = 21)
# # 
# # 
# # 
# # 
# # 
# # result_16_21 <- estimate_gmm_period(period_begin = 16, period_end = 21)
# # result_15_21 <- estimate_gmm_period(period_begin = 15, period_end = 21)
# # 
# # result_14_21 <- estimate_gmm_period(period_begin = 14, period_end = 21)
# # result_13_21 <- estimate_gmm_period(period_begin = 13, period_end = 21)
# # result_12_21 <- estimate_gmm_period(period_begin = 12, period_end = 21)
# # result_11_21 <- estimate_gmm_period(period_begin = 11, period_end = 21)
# # result_10_21 <- estimate_gmm_period(period_begin = 10, period_end = 21)
# result_9_21 <- estimate_gmm_period(period_begin = 9, period_end = 21)



# result_1_4 <- estimate_gmm_period(period_begin = 1, period_end = 4) 
# 
# results_list <- list(
#   result_1_1 = result_1_1,
#   result_1_2 = result_1_2
# )
# 
# results_list[[3]] <- result_1_3
# results_list[[4]] <- result_1_4
# names(results_list)[4] <- "result_1_4"
# 
# 
# names(table(df_qlfs$period1))

param_init <-  data.frame(theta_0 = 0.04440373, theta_01 = 0.2745643, theta_10 = 0.3716039, theta_1 = 0.04090615, pi = 0.01351413)
param_init_transformed <- logit_transform(param_init)
model_gmm <- gmm::gmm(calc_gmm_moments, df_x, unlist(param_init_transformed), weightsMatrix = diag(5))
print(logit_inverse(model_gmm$coefficients))

# gmm_hat <- calc_gmm_moments(res$coefficients, df_x)

# sum(colMeans(gmm_hat)^2)

summary(model_gmm)
print(gmm::specTest(model_gmm))

gb <- colMeans(model_gmm$gt)
j <- try(crossprod(gb,solve(model_gmm$v,gb))*model_gmm$n, silent=TRUE)
k <- noquote(cbind(j, ifelse(model_gmm$df>0,pchisq(j,model_gmm$df,
                                           lower.tail = FALSE),
                             "*******")))
J_test <- noquote(paste("J-Test: degrees of freedom is ",model_gmm$df,sep=""))
