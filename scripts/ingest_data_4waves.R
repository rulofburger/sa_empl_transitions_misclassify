
# Import
# df_qlfs <- haven::read_dta(
#   "data/raw/qlfs_raw_4waves.dta"
# )

df_qlfs <- readRDS(
    "data/raw/df_qlfs_A.rds"
  ) %>%
  filter(age1 > 17 & age1 < 56) %>% 
  filter(!is.na(employed1)) %>%
  filter(!is.na(employed2)) %>%
  filter(!is.na(employed3)) %>%
  filter(!is.na(employed4)) %>% 
  # filter(period1 >= 35 &  period1 <= 40) %>% 
  rename(
    y1 = employed1,
    y2 = employed2,
    y3 = employed3,
    y4 = employed4
  ) %>% 
  select(y1, y2, y3, y4, weight1, weight2, weight3, weight4, period1, period2, period3, period4) %>% 
  mutate(weight = (weight1*weight2*weight3*weight4)^0.25) %>% 
  as.data.frame

# Keep only age, educ, status for 3 waves
# df_qlfs <- df_qlfs %>%
#   select(contains(c("age", "educ", "status", "tenure", "timegap", "formal")))  # keep age, educ, empl status
#   # select(-contains(c("4"))) # remove columns on wage, and 4th wave
# #select(-contains(c("wage", "4"))) # remove columns on wage, and 4th wave

# Make Employment Binary



# 
# df_qlfs <- df_qlfs %>%
#   mutate(
#     Status_Q1 = ifelse(
#       status1 == 1, 1, 0
#     ),
#     Status_Q2 = ifelse(
#       status2 == 1, 1, 0
#     ),
#     Status_Q3 = ifelse(
#       status3 == 1, 1, 0
#     ),
#     Status_Q4= ifelse(
#       status4 == 1, 1, 0
#     ),
#     tenure1 = if_else(Status_Q1 == 0, 0, tenure1),
#     tenure2 = if_else(Status_Q2 == 0, 0, tenure2),
#     tenure3 = if_else(Status_Q3 == 0, 0, tenure3),
#     tenure4 = if_else(Status_Q4 == 0, 0, tenure4),
#     timegap1 = if_else(Status_Q1 == 1, 0, timegap1),
#     timegap2 = if_else(Status_Q2 == 1, 0, timegap2),
#     timegap3 = if_else(Status_Q3 == 1, 0, timegap3),
#     timegap4 = if_else(Status_Q4 == 1, 0, timegap4)
#     # tenure1 = ifelse(is.na(tenure1), 0, tenure1),
#     # tenure2 = ifelse(is.na(tenure2), 0, tenure2),
#     # tenure3 = ifelse(is.na(tenure3), 0, tenure3),
#     # timegap1 = ifelse(is.na(timegap1), 0,
#     #                   ifelse(timegap1>7, 0, timegap1)),
#     # timegap2 = ifelse(is.na(timegap2), 0,
#     #                   ifelse(timegap2>7, 0, timegap2)),
#     # timegap3 = ifelse(is.na(timegap3), 0,
#     #                   ifelse(timegap2>7, 0, timegap3))
#     
#   ) %>%
#   select(-c(status1, status2, status3, status4)) # remove other statuses


# Create timegap variables
# df_qlfs <- df_qlfs %>%
#   mutate(
#     TG_Q1 = case_when(
#       timegap1 == 0 ~ 0,
#       timegap1 == 1 ~ 1.5,
#       timegap1 == 2 ~ 4.5,
#       timegap1 == 3 ~ 7.5,
#       timegap1 == 4 ~ 10.5,
#       timegap1 == 5 ~ 24,
#       timegap1 == 6 ~ 48,
#       timegap1 == 7 ~ 90,
#       timegap1 == 8 ~ NA_real_,
#       timegap1 == 99 ~ NA_real_
#     )
#   )%>%
#   mutate(
#     TG_Q2 = case_when(
#       timegap2 == 0 ~ 0,
#       timegap2 == 1 ~ 1.5,
#       timegap2 == 2 ~ 4.5,
#       timegap2 == 3 ~ 7.5,
#       timegap2 == 4 ~ 10.5,
#       timegap2 == 5 ~ 24,
#       timegap2 == 6 ~ 48,
#       timegap2 == 7 ~ 90,
#       timegap2 == 8 ~ NA_real_,
#       timegap2 == 99 ~ NA_real_
#     )
#   )%>%
#   mutate(
#     TG_Q3 = case_when(
#       timegap3 == 0 ~ 0,
#       timegap3 == 1 ~ 1.5,
#       timegap3 == 2 ~ 4.5,
#       timegap3 == 3 ~ 7.5,
#       timegap3 == 4 ~ 10.5,
#       timegap3 == 5 ~ 24,
#       timegap3 == 6 ~ 48,
#       timegap3 == 7 ~ 90,
#       timegap3 == 8 ~ NA_real_,
#       timegap3 == 99 ~ NA_real_
#     )
#   )%>%
#   mutate(
#     TG_Q4 = case_when(
#       timegap4 == 0 ~ 0,
#       timegap4 == 1 ~ 1.5,
#       timegap4 == 2 ~ 4.5,
#       timegap4 == 3 ~ 7.5,
#       timegap4 == 4 ~ 10.5,
#       timegap4 == 5 ~ 24,
#       timegap4 == 6 ~ 48,
#       timegap4 == 7 ~ 90,
#       timegap4 == 8 ~ NA_real_,
#       timegap4 == 99 ~ NA_real_
#     )
#   )


# df_qlfs <- df_qlfs %>%
#   filter(!is.na(Status_Q1)) %>%
#   filter(!is.na(Status_Q2)) %>%
#   filter(!is.na(Status_Q3)) %>%
#   filter(!is.na(tenure1)) %>%
#   filter(!is.na(tenure2)) %>%
#   filter(!is.na(tenure3)) %>%
#   filter(!is.na(TG_Q1)) %>%
#   filter(!is.na(TG_Q2)) %>%
#   filter(!is.na(TG_Q3))


# 
# df_qlfs %>%
#   filter(employed1 == 1) %>%
#   ggplot(aes(y = employed2, x = log(tenure1 + 1.5))) +
#   geom_smooth(
#     method = "gam",
#     formula = y ~ s(x, bs = "cs", k = 5)
#   )

