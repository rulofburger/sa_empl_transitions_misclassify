
# Import
# df_qlfs <- haven::read_dta(
#   "data/raw/qlfs_raw_4waves.dta"
# )


# Import
df_qlfs <- readRDS(
  "data/raw/df_qlfs_A.rds"
) %>%
  filter(age1 > 17 & age1 < 56) %>%  
  filter(!is.na(employed1)) %>%
  filter(!is.na(employed2)) %>%
  filter(!is.na(employed3)) %>%
  filter(!is.na(employed4)) %>%
  select(contains(c("employed", "age", "educ", "race", "female", "weight", "tenure", "timegap", "formal", "period", "contracttype", "neverworked"))) %>%  # keep age, educ, empl status
  # select(-contains(c("4"))) %>% 
  rename(
    y1 = employed1,
    y2 = employed2,
    y3 = employed3,
    y4 = employed4,
  ) %>% 
  mutate(weight = (weight1*weight2*weight3*weight4)^0.25) %>% 
  as.data.frame

