
# Import
df_qlfs <- readRDS(
  "data/raw/df_qlfs_A.rds"
  ) %>%
  filter(age1 > 17 & age1 < 56) %>% 
  filter(!is.na(employed1)) %>%
  filter(!is.na(employed2)) %>%
  filter(!is.na(employed3)) %>%
  select(contains(c("employed", "age", "educ", "weight", "tenure", "timegap", "formal", "period", "contracttype", "neverworked"))) %>%  # keep age, educ, empl status
  select(-contains(c("4"))) %>% 
  rename(
    y1 = employed1,
    y2 = employed2,
    y3 = employed3
  ) %>% 
  select(-contains(c("4"))) %>% 
  mutate(weight = (weight1*weight2*weight3)^0.333) %>% 
  as.data.frame


# Make Employment Binary
df_qlfs <- df_qlfs %>%
  mutate(
    tenure1 = if_else(y1 == 0, 0, tenure1),
    tenure2 = if_else(y2 == 0, 0, tenure2),
    tenure3 = if_else(y3 == 0, 0, tenure3),
    timegap1 = if_else(y1 == 1, 0, timegap1),
    timegap2 = if_else(y2 == 1, 0, timegap2),
    timegap3 = if_else(y3 == 1, 0, timegap3),
    timegap1 = if_else(y1 == 0 & neverworked1 == 1, 0, timegap1),
    timegap2 = if_else(y2 == 0 & neverworked2 == 1, 0, timegap2),
    timegap3 = if_else(y3 == 0 & neverworked3 == 1, 0, timegap3),
  ) %>% 
  mutate(
    tenure1 = if_else(tenure1 < 0, NA_real_, tenure1),
    tenure2 = if_else(tenure2 < 0, NA_real_, tenure2),
    tenure3 = if_else(tenure3 < 0, NA_real_, tenure3)
    )


# Create timegap variables
df_qlfs <- df_qlfs %>%
  mutate(
    timegap1 = case_when(
      timegap1 == 0 ~ 0,
      timegap1 == 1 ~ 1.5,
      timegap1 == 2 ~ 4.5,
      timegap1 == 3 ~ 7.5,
      timegap1 == 4 ~ 10.5,
      timegap1 == 5 ~ 24,
      timegap1 == 6 ~ 48,
      timegap1 == 7 ~ 90,
      timegap1 == 8 ~ NA_real_,
      timegap1 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    timegap2 = case_when(
      timegap2 == 0 ~ 0,
      timegap2 == 1 ~ 1.5,
      timegap2 == 2 ~ 4.5,
      timegap2 == 3 ~ 7.5,
      timegap2 == 4 ~ 10.5,
      timegap2 == 5 ~ 24,
      timegap2 == 6 ~ 48,
      timegap2 == 7 ~ 90,
      timegap2 == 8 ~ NA_real_,
      timegap2 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    timegap3 = case_when(
      timegap3 == 0 ~ 0,
      timegap3 == 1 ~ 1.5,
      timegap3 == 2 ~ 4.5,
      timegap3 == 3 ~ 7.5,
      timegap3 == 4 ~ 10.5,
      timegap3 == 5 ~ 24,
      timegap3 == 6 ~ 48,
      timegap3 == 7 ~ 90,
      timegap3 == 8 ~ NA_real_,
      timegap3 == 99 ~ NA_real_
    )
  )


df_qlfs %>%
  filter(y1 == 1) %>%
  ggplot(aes(y = y2, x = log(tenure1 + 1.5))) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 5)
  )

